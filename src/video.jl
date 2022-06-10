module video
    using Glob, Printf
    using VideoIO
    include(".raw.jl")
    using .raw 
    using MATLAB

    function getVidCollection(animal::String, day::Int)::Vector
        # Get the list of files for the epoch, videoTS files
        # Load up the collection of video ts files
    end
    function getTsCollection(animal::String, day::Int)::Vector
        # Get the list of files for the epoch, vidoes and videoTS files
        # Load up the collection of video ts files
    end
    function ts2videots(animal::String, day::Int, timestamp::Real)
        tsCollection = getTsCollection(animal, day)
        ts2videots(timestamp, tsCollection)
    end
    function ts2videots(timestamp::Real, tsCollection::Vector)
    end
    function ts2frame()
    end

    function frameattime(vid, time; cropx=[], cropy=[])
        if time == 0
            vid = seekstart(vid)
        else
            currtime = gettime(vid)
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
    function load(pos...; kws...)
        videopath = get_path(pos...; kws...)
        stream    = VideoIO.open(videopath)
        vid       = VideoIO.openvideo(stream)
    end
end
export video
