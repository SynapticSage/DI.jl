module video

    using ..Load
    using Glob, Printf
    using VideoIO
    using MATLAB
    using Infiltrator
    using ProtoStructs
    using DataFrames, DataFramesMeta
    using DimensionalData

    function __init__()
        mat"addpath('/home/ryoung/Code/pipeline/TrodesToMatlab')"
    end

    videoFolders=Dict(
                      "RY16"=>"/media/ryoung/GenuDrive/RY16_direct/videos/",
                      "RY22"=>"/media/ryoung/Ark/RY22_direct/videos/",
                     )

    function getVidCollection(animal::String, day::Int)::Vector
        # Get the list of files for the epoch, videoTS files
        # Load up the collection of video ts files
        globstr = "RY16video$day-*.mp4"
        glob(globstr, videoFolders[animal])

    end
    function getTsCollection(animal::String, day::Int)::Vector
        # Get the list of files for the epoch, vidoes and videoTS files
        # Load up the collection of video ts files
        globstr = "RY16timestamp$day-*.dat"
        glob(globstr, videoFolders[animal])
    end
    function ts2videots(animal::String, day::Int, timestamp::Real)
        tsCollection = getTsCollection(animal, day)
        ts2videots(timestamp, tsCollection)
    end

    """
    ts2videots

    get a video timestamp from timestamp
    """
    function ts2videots(timestamp::Real, tsCollection::Vector)
    end
    """
    ts2frame

    get a frame from a timestamp
    """
    function ts2frame()
    end

    function frameattime(vid, time; cropx=[], cropy=[], timecoord=nothing)
        if time == 0
            vid = seekstart(vid)
        else
            currtime = gettime(vid)
            if timecoord !== nothing
                currtime = currtime - minimum(timecoord)
            end
            vid = seek(vid, time - currtime)
        end
        seek(vid, time)
        img = read(vid)'
        if (length(cropx) & length(cropy)) > 0
            cropx = Int.(round.(pxtocm.(cropx)))
            cropy = Int.(round.(pxtocm.(cropy)))
            img=img[cropx[1] : cropx[2], 
                cropy[1] : cropy[2]]
        end
        #img = img[:, end:-1:begin]
        img
    end

    function load_videots(file::String; center_by_loadmintime::Bool=true)
        @info "opening $file" "readCameraModuleTimeStamps('$file')"
        ts = mat"readCameraModuleTimeStamps($file)"
        #ts = mat"readCameraModuleTimeStamps('/media/ryoung/GenuDrive/RY16_direct/videos/RY16timestamp36-01.dat')"
        center_by_loadmintime ? ts - Load.min_time_records[end] : ts
    end
    function load_videots(animal::String, day::Int, epoch::Int; 
            center_by_loadmintime::Bool=true)
        file = Load.video.getTsCollection("RY16",day)[epoch]
        @info  "file" file
        load_videots(file; center_by_loadmintime)
    end
    function load_video(file::String)
        stream    = VideoIO.open(file)
        VideoIO.openvideo(stream)
    end
    function load_video(animal::String, day::Int, epoch::Int)
        file = Load.video.getVidCollection(animal, day)
        file = file[epoch]
        @info  "file" file
        load_video(file)
    end
    
    mutable struct VideoArray
        # The actualk video reader object
        vid::VideoIO.VideoReader
        # Video statistics
        totaltime::Float64
        totalframe::Int32
        framespertime::Float64
        # A record of behavior time of the frames and distance of the behavior
        # time to the video time
        behaviortime::Vector{Float64}
        behaviorminusvidtime::Float64
        # Where we are in the video
        currvidtime::Float64
        # Pixel space cropping
        cropx::Tuple{Int16, Int16}
        cropy::Tuple{Int16, Int16}
        # Axis labels for the video axes
        xaxis::Vector{Float64}
        yaxis::Vector{Float64}
    end
    function getVideoArray(vidpath, vid, ts; cropx=nothing, cropy=nothing)
        # Image crop calculations
        # -----------------------
        sizeimg = size(read(vid))
        if cropx === nothing; cropx = (1, sizeimg[1])
        elseif cropx[1] == -Inf; cropx = (1, cropx[2])
        elseif cropx[2] == Inf; cropx = (cropx[1], sizeimg[1])
        end
        if cropy === nothing; cropy = (1, sizeimg[2])
        elseif cropy[1] == -Inf; cropy = (1, cropy[2])
        elseif cropy[2] == Inf; cropy = (cropy[1], sizeimg[2])
        end
        # Basic video statistics
        totalframe = length(ts)
        totaltime = VideoIO.get_duration(vidpath) 
        vidtime = gettime(vid)
        framespertime = totalframe/totaltime
        # Relatinoship to behavior times
        vididx = Int32(round(vidtime*framespertime)) + 1
        behaviorminusvidtime = ts[vididx] - vidtime
        # Xaxis and Yaxis
        xaxis = Load.pxtocm(collect(1:sizeimg[1]))
        yaxis = Load.pxtocm(collect(1:sizeimg[2]))
        # Create the video array object (a video object that can be indexed like an array)
        VideoArray(vid, totaltime, totalframe, totalframe/totaltime,
                   ts, behaviorminusvidtime,
                   vidtime,
                   cropx, cropy,
                   xaxis, yaxis
                  )
    end
    function cropfromtask(animal::String, day::Int, epoch::Int)
        task = Load.load_task(animal, day, epoch)
        task = @subset(task, :epoch .== epoch)
        @infiltrate
    end
    function indtotime(vid::VideoArray, i::Int)::Float64
        @fastmath i/vid.totalframe * vid.totaltime
    end
    function timetoind(vid::VideoArray, t::Float64)::Int32
        Int32(round(t/vid.totaltime * vid.totalframe))
    end
    function Base.getindex(vid::VideoArray, i::Int)
        t = indtotime(vid, i)
        seek(vid.vid, t)
        img = read(vid.vid)'
        updatetime!(vid)
        img
    end
    function Base.getindex(vid::VideoArray, t::Float64)
        seek(vid.vid, t)
        img = read(vid.vid)
        updatetime!(vid)
        DimArray(img, X(vid.xaxis), Y(vid.yaxis))
    end
    function Base.getindex(vid::VideoArray; time::Union{Float64,Float32,Int})
        vidtime = time - vid.behaviorminusvidtime
        DimArray(vid[vidtime], X(vid.xaxis), Y(vid.yaxis))
    end
    Base.close(vid::VideoArray) = close(vid.vid)
    function updatetime!(vid::VideoArray)
        @fastmath vid.currvidtime = gettime(vid.vid)
    end


    function load(animal::String, day::Int, epoch::Int)
        file = video.getVidCollection(animal, day)
        file = file[epoch]
        vid = load_video(file)
        ts = load_videots(animal, day, epoch)
        getVideoArray(file, vid, ts)
    end


    # ============================================================
    # ============================================================
    # ============================================================
    # ============================================================
    # ============================================================
    # ============================================================


    """
        get_path

    depcrecated
    """
    function get_path(animal, day, epoch; dayfactor=0, 
            guessdayfactor=true,
            source="deeplabcut")
        if guessdayfactor
            dayfactor = raw.animal_dayfactor[animal]
        end
        day += dayfactor
        if source == "deeplabcut"
            folder_path = "/Volumes/Colliculus/deeplabcut/" * 
                         "goalmaze_tape-Ryan-2020-05-28/videos/"
        end
        possible_files = 
        glob("$(animal)_$(@sprintf("%02d",day))_$(@sprintf("%02d",epoch))_*.mp4",
                 folder_path)
        videopath = possible_files[1]
        return videopath
    end

end
