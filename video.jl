module video

    using ..Load
    using  RecipesBase
    import  Utils
    using Glob, Printf
    using VideoIO
    using MATLAB
    using Infiltrator
    using ProtoStructs
    using DataFrames, DataFramesMeta
    using DimensionalData
    using FixedPointNumbers
    using Colors

    function __init__()
        mat"addpath('/home/ryoung/Code/pipeline/TrodesToMatlab')"
    end

    videoFolders=Dict(
                      "RY16"=>"/media/ryoung/GenuDrive/RY16_direct/videos/",
                      "RY22"=>"/media/ryoung/Ark/RY22_direct/videos/",
                     )

    function getVideoFiles(animal::String, day::Int)::Vector
        # Get the list of files for the epoch, videoTS files
        # Load up the collection of video ts files
        globstr = "RY16video$day-*.mp4"
        glob(globstr, videoFolders[animal])

    end
    function getTsFiles(animal::String, day::Int)::Vector
        # Get the list of files for the epoch, vidoes and videoTS files
        # Load up the collection of video ts files
        globstr = "RY16timestamp$day-*.dat"
        glob(globstr, videoFolders[animal])
    end
    function ts2videots(animal::String, day::Int, timestamp::Real)
        tsCollection = getTsFiles(animal, day)
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
        center_by_loadmintime ? ts .- Load.min_time_records[end] : ts
    end
    function load_videots(animal::String, day::Int, epoch::Int; 
            center_by_loadmintime::Bool=true)
        file = Load.video.getTsFiles("RY16",day)[epoch]
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

    @recipe function plotviddimarray(v::DimArray{RGB{FixedPointNumbers.N0f8}})
        seriestype --> :heatmap
        (collect.(v.dims)[2:-1:1]..., v.data)
    end
    @userplot PlotVid
    @recipe function plotvid(plt::PlotVid)
        v=plt.args[1]
        seriestype --> :heatmap
        (collect.(v.dims)[2:-1:1]..., v.data)
    end
    
    mutable struct VideoObj
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
        cropx::BitVector
        cropy::BitVector
        # Axis labels for the video axes
        xaxis::Vector{Float64}
        yaxis::Vector{Float64}
    end
    function getVideObj(vidpath, vid, ts; cropx=nothing, cropy=nothing,
        animal=nothing, day=nothing, epoch=nothing)
        cropx, cropy, xaxis, yaxis = computeimagecrop(vid;cropx,cropy,animal,day,epoch)
        # Basic video statistics
        totalframe = length(ts)
        totaltime = VideoIO.get_duration(vidpath) 
        vidtime = gettime(vid)
        framespertime = totalframe/totaltime
        # Relatinoship to behavior times
        vididx = Int32(round(vidtime*framespertime)) + 1
        behaviorminusvidtime = ts[vididx] - vidtime
        # Create the video array object (a video object that can be indexed like an array)
        VideoObj(vid, totaltime, totalframe, totalframe/totaltime,
                   ts, behaviorminusvidtime,
                   vidtime,
                   cropx, cropy,
                   xaxis, yaxis
                  )
    end
    function computeimagecrop(vid; cropx=nothing, cropy=nothing,
        animal=nothing, day=nothing, epoch=nothing)

        sizeimg = size(read(vid))
        # Image crop calculations
        # -----------------------

        vidxax = Load.pxtocm(1:sizeimg[1])
        vidyax = Load.pxtocm(1:sizeimg[2])

        if animal !== nothing && day !== nothing && epoch !== nothing &&
            cropx === nothing && cropy === nothing
            cropyN, cropxN = taskboundaries(animal, day, epoch)
            cropx = cropx === nothing ? cropxN : cropx
            cropy = cropy === nothing ? cropyN : cropy
            @info "crop" cropx  cropy
        else
           cropx = cropx === nothing ? (-Inf,Inf) : cropx
           cropy = cropy === nothing ? (-Inf,Inf) : cropx
        end
        # Xaxis and Yaxis
        sizeimg = size(read(vid))
        xaxis = Load.pxtocm(collect(1:sizeimg[1]))
        yaxis = Load.pxtocm(collect(1:sizeimg[2]))
        cropx = Utils.in_range(vidxax, cropx)
        cropy = Utils.in_range(vidyax, cropy)
        xaxis = xaxis[cropx]
        yaxis = yaxis[cropy]

        cropx, cropy, xaxis, yaxis
    end
    function taskboundaries(animal::String, day::Int, epoch::Int)
        task = Load.task.load_task(animal, day)
        task = @subset(task, :epoch .== epoch, :name .== "boundary")
        extrema(task.x), extrema(task.y)
    end
    function indtotime(vid::VideoObj, i::Int)::Float64
        @fastmath i/vid.totalframe * vid.totaltime
    end
    function timetoind(vid::VideoObj, t::Float64)::Int32
        Int32(round(t/vid.totaltime * vid.totalframe))
    end
    function Base.getindex(vid::VideoObj, i::Int)
        t = indtotime(vid, i)
        seek(vid.vid, t)
        img = read(vid.vid)
        updatetime!(vid)
        DimArray(img[vid.cropx, vid.cropy], (X(vid.xaxis), Y(vid.yaxis)))
    end
    function Base.getindex(vid::VideoObj, t::Float64)
        seek(vid.vid, t)
        img = read(vid.vid)
        updatetime!(vid)
        DimArray(Array(img[vid.cropx,vid.cropy]), (X(vid.xaxis), Y(vid.yaxis)))
    end
    function Base.getindex(vid::VideoObj; time::Union{Float64,Float32,Int})
        vidtime = time - vid.behaviorminusvidtime
        vid[vidtime]
    end
    Base.close(vid::VideoObj) = close(vid.vid)
    function updatetime!(vid::VideoObj)
        @fastmath vid.currvidtime = gettime(vid.vid)
    end

    struct VideoCollection
      vids::Vector{VideoObj}
      starts::Vector{Float64}
      stops::Vector{Float64}
    end
    function  getVideoCollection(vids::Vector{VideoObj}; 
            end_tolerance=0.2, begin_tolerance=0.2)
        starts = []
        stops = []
        for (v,vid) in enumerate(vids)
            v == 1 ?
                push!(starts,vid.behaviortime[begin]-begin_tolerance) :
                push!(starts,vid.behaviortime[begin])
            v == length(vids) ?
                push!(stops, vid.behaviortime[end]+end_tolerance) :
                push!(stops, vid.behaviortime[end])
        end
        VideoCollection(vids, starts, stops)
    end
    function Base.getindex(vids::VideoCollection, t::Union{Float32,Float64})::DimArray
        v = findfirst(t .>= vids.starts .&& t .<= vids.stops)
        v === nothing ? getNullVideoImage(vids.vids[1]) : vids.vids[v][time=t]
    end
    function Base.close(vids::VideoCollection)
        [close(vid) for vid in vids.vids]
        nothing
    end
    function getNullVideoImage(vid::VideoObj)
        DimArray(fill(NaN, length(vid.xaxis), length(vid.yaxis)),
                 (X(vid.xaxis), Y(vid.yaxis)))
    end

    function load(animal::String, day::Int, epoch::Int)::VideoObj
        file = video.getVideoFiles(animal, day)
        file = file[epoch]
        vid = load_video(file)
        ts = load_videots(animal, day, epoch)
        getVideObj(file, vid, ts; animal, day, epoch)
    end
    function load(animal::String, day::Int)::VideoCollection
        task = Load.task.load_task(animal, day)
        epochs = unique(task.epoch)
        vids = Vector{VideoObj}([])
        for epoch in epochs
            push!(vids, video.load(animal, Int(day), Int(epoch)))
        end
        getVideoCollection(vids)
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
