module video

    using Glob, Printf
    using VideoIO
    using MATLAB
    using ..Load
    using Infiltrator
    using ProtoStructs

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

    function load_videots(file::String)
        @info "opening $file" "readCameraModuleTimeStamps('$file')"
        ts = mat"readCameraModuleTimeStamps($file)"
        #ts = mat"readCameraModuleTimeStamps('/media/ryoung/GenuDrive/RY16_direct/videos/RY16timestamp36-01.dat')"
        ts
    end
    function load_videots(animal::String, day::Int, epoch::Int)
        file = Load.video.getTsCollection("RY16",day)[epoch]
        @info  "file" file
        load_videots(file)
    end
    function load_video(file::String)
        stream    = VideoIO.open(file)
        VideoIO.openvideo(stream)
    end
    function load_video(animal::String, day::Int, epoch::Int)
        file = Load.video.getVidCollection("RY16",day)
        file = file[epoch]
        @info  "file" file
        load_video(file)
    end
    
    struct VideoArray
        vid::VideoIO.VideoReader
        totaltime::Float64
        totalframe::Float64
        framespertime::Float64
        behaviortime::Vector{Float64}
        behaviorminusvidtime::Float64
        currvidtime::Float64
        currvidindex::Int32
    end
    function getVideoArray(vidpath, vid, ts)
        totalframe = length(ts)
        totaltime = VideoIO.get_duration(vidpath) 
        vidtime = gettime(vid)
        framespertime = totalframe/totaltime
        vididx = Int32(round(vidtime*framespertime))
        behaviorminusvidtime = ts[vididx] - vidtime

        VideoArray(vid, totaltime, totalframe, totalframe/totaltime,
                   behaviorminusvidtime, ts, vidtime, vididx)
    end
    function getindex(vid::VideoArray, i::Int)
        d = i - vid.currvidindex
        td = d / vid.framespertime
        seek(vid, td)
        vid.currvidindex += d
        vid.currvidtime  += td
        img = read(vid)'
        img
    end
    function getindex(vid::VideoArray; time::Float64)
        td = time - vid.currvidtime
        d  = td * vid.framespertime
        seek(vid, td)
        vid.currvidindex += d
        vid.currvidtime  += td
        img = read(vid)'
        img
    end

    function load(animal::String, day::Int, epoch::Int)
        file = video.getVidCollection("RY16",day)
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
