using ..DI
using  RecipesBase
import  DIutils
using Glob, Printf
using VideoIO
using MATLAB
using Infiltrator
using ProtoStructs
using DataFrames, DataFramesMeta
using DimensionalData
using FixedPointNumbers
using Colors
using Statistics
import VideoIO

#getindexmethod = :const # searchsortednearest | const | constalignnearest

function __init__()
    mat"addpath('/home/ryoung/Code/pipeline/TrodesToMatlab')"
end

videoFolders=Dict(
                  "RY16"=>"/media/ryoung/GenuDrive/RY16_direct/videos/",
                  "RY22"=>"/media/ryoung/Ark/RY22_direct/videos/"
                 )

function getVideoFiles(animal::String, day::Int)::Vector
    # Get the list of files for the epoch, videoTS files
    # DI up the collection of video ts files
    globstr = "$(animal)video$day-*.mp4"
    glob(globstr, videoFolders[animal])

end
function getTsFiles(animal::String, day::Int)::Vector
    # Get the list of files for the epoch, vidoes and videoTS files
    # DI up the collection of video ts files
    globstr = "$(animal)timestamp$day-*.dat"
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

function load_videots(file::String; 
        initial_behavior_time::Union{Real,Nothing}=nothing,
        center_by_loadmintime::Bool=true)
    @info "opening $file" "readCameraModuleTimeStamps('$file')"
    ts = mat"readCameraModuleTimeStamps($file)"
    if initial_behavior_time !== nothing
        δ = initial_behavior_time - first(ts)
        ts = ts .+ δ
    end
    centering_time = (isempty(DI.min_time_records) ? 0 : DI.min_time_records[end])
    center_by_loadmintime ? (ts .- centering_time) : ts 
end
function load_videots(animal::String, day::Int, epoch::Int; 
        correct_behavior_ts_offset::Bool=true,
        center_by_loadmintime::Bool=true)
    file = DI.video.getTsFiles(animal,day)[epoch]
    if correct_behavior_ts_offset
        initial_behavior_time = firstbehaviortime(animal,day,epoch)
        initial_behavior_time = initial_behavior_time === nothing ? 
                                 0 : initial_behavior_time
    else
        initial_behavior_time = 0
    end
    @info  "file" file
    load_videots(file; center_by_loadmintime, initial_behavior_time)
end
function load_video(file::String)
    stream    = VideoIO.open(file)
    VideoIO.openvideo(stream)
end
function load_video(animal::String, day::Int, epoch::Int)
    file = DI.video.getVidCollection(animal, day)
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
    totaltimestamptime::Float64
    totalframe::Int32
    framespertime::Float64
    # A record of behavior time of the frames and distance of the behavior
    # time to the video time
    behaviortime::Vector{Float64}
    behaviorminusvidtime::Float64
    # Where we are in the video
    currvidtime::Float64
    # Search Table
    searchtable::DataFrame
    # Pixel space cropping
    cropx::BitVector
    cropy::BitVector
    # Axis labels for the video axes
    xaxis::Vector{Float64}
    yaxis::Vector{Float64}
end
function getVideObj(vidpath, vid, timestamps; cropx=nothing, cropy=nothing,
    animal=nothing, day=nothing, epoch=nothing)
    cropx, cropy, xaxis, yaxis = computeimagecrop(vid;cropx,cropy,animal,day,epoch)
    # Basic video statistics
    totalframe = length(timestamps)
    totaltime = VideoIO.get_duration(vidpath) 
    totaltimestamptime = diff(collect(extrema(timestamps)))[1]
    currentvideotime = gettime(vid) # relative time into video
    framespertime = totalframe/totaltime # how many frames to move per time
    # Setup a searchtable
    searchtable = DataFrame([collect(1:totalframe), collect(LinRange(0,totaltime,totalframe)), timestamps], 
                            [:frames, :linearts, :timestamp])
    # Relatinoship to behavior times
    vididx = max(Int32(round(currentvideotime*framespertime)), 1)
    behaviorminusvidtime = timestamps[vididx] - currentvideotime # the time stamp file minus the vidtime
    # Create the video array object (a video object that can be indexed like an array)
    VideoObj(vid, totaltime, totaltimestamptime, 
             totalframe, totalframe/totaltime,
               timestamps, behaviorminusvidtime,
               currentvideotime, searchtable,
               cropx, cropy,
               xaxis, yaxis
              )
end
function computeimagecrop(vid; cropx=nothing, cropy=nothing,
    animal=nothing, day=nothing, epoch=nothing)

    sizeimg = size(read(vid))
    # Image crop calculations
    # -----------------------

    vidxax = DI.pxtocm(1:sizeimg[1])
    vidyax = DI.pxtocm(1:sizeimg[2])

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
    xaxis = DI.pxtocm(collect(1:sizeimg[1]))
    yaxis = DI.pxtocm(collect(1:sizeimg[2]))
    cropx = DIutils.in_range(vidxax, cropx)
    cropy = DIutils.in_range(vidyax, cropy)
    xaxis = xaxis[cropx]
    yaxis = yaxis[cropy]

    cropx, cropy, xaxis, yaxis
end
"""
    taskboundaries(animal::String, day::Int, epoch::Int)
return the boundaries of the task
"""
function taskboundaries(animal::String, day::Int, epoch::Int)
    task = DI.task.load_task(animal, day)
    task = @subset(task, :epoch .== epoch, :name .== "boundary")
    extrema(task.x), extrema(task.y)
end
"""
    firstbehaviorepoch(animal::String, day::Int)
return the first behavioral epoch
"""
function firstbehaviorepoch(animal::String, day::Int)
    task = subset(DI.task.load_task(animal, day), :type => t -> t .== "run")
    minimum(task.epoch)
end
function firstbehaviortime(animal::String, day::Int, epoch::Union{Int,Nothing}=nothing)::Union{Nothing,Real}
    epoch = epoch === nothing ? firstbehaviorepoch(animal, day) : epoch
    beh = subset(DI.load_behavior(animal, day), :epoch=>e->e.==epoch)
    isempty(beh) ? nothing : first(beh.time)
end
function indtotime(vid::VideoObj, i::Int)::Float64
    @fastmath i/vid.totalframe * vid.totaltime
end
function timetoind(vid::VideoObj, t::Float64)::Int32
    Int32(round(t/vid.totaltime * vid.totalframe))
end
function VideoIO.gettime(vid::VideoObj)
    VideoIO.gettime(vid.vid) + vid.behaviorminusvidtime
end
function VideoIO.gettime(vid::VideoObj, i::Int)
    t = indtotime(vid, i)
    seek(vid.vid, t)
    updatetime!(vid)
    gettime(vid.vid) + vid.behaviorminusvidtime
end
function VideoIO.gettime(vid::VideoObj, t::T where T <: Union{Float32, Float64})
    seek(vid.vid, t)
    updatetime!(vid)
    gettime(vid.vid) + vid.behaviorminusvidtime
end
function Base.getindex(vid::VideoObj, i::Int)
    t = indtotime(vid, i) # translate index of video to time of video
    seek(vid.vid, t)  # seek to time
    img = read(vid.vid) # get image
    updatetime!(vid) # update time
    DimArray(img[vid.cropx, vid.cropy], (X(vid.xaxis), Y(vid.yaxis)))
end
function Base.getindex(vid::VideoObj, t::Float64)
    search = vid.searchtable[DIutils.searchsortednearest(vid.searchtable.timestamp, t),:]
    seek(vid.vid, search.linearts)
    img = read(vid.vid)
    updatetime!(vid)
    DimArray(Array(img[vid.cropx,vid.cropy]), (X(vid.xaxis), Y(vid.yaxis)))
end
function Base.getindex(vid::VideoObj; time::Union{Float64,Float32,Int})
    #vidtime = time - vid.behaviorminusvidtime
    vidtime = Float64(time)
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
function Base.getindex(vids::VideoCollection, t::Vector{<:Union{Float32,Float64}})::Vector{DimArray}
    Base.getindex.([vids], t)
end
function Base.getindex(vids::VideoCollection; epoch::Int)::DimArray
    t = mean([vids.starts[epoch] , vids.stops[epoch]])
    vids.vids[epoch][time=t]
end
function VideoIO.gettime(vids::VideoCollection, t::T where T <: Union{Float32,Float64})
    v = findfirst(t .>= vids.starts .&& t .<= vids.stops)
    v === nothing ? VideoIO.gettime(vids.vids[1]) : VideoIO.gettime(vids.vids[v], t)
end
function VideoIO.gettime(vids::VideoCollection, T::Vector{<:Union{Float32,Float64}})
    [VideoIO.gettime(vids, t) for t in T]
end
function VideoIO.gettime(vids::VideoCollection; epoch::Int)::DimArray
    t = mean([vids.starts[epoch] , vids.stops[epoch]])
    @error "not implemented"
    vids.vids[epoch][time=t]
end
function Base.close(vids::VideoCollection)
    [close(vid) for vid in vids.vids]
    nothing
end
function getNullVideoImage(vid::VideoObj)
    DimArray(fill(NaN, length(vid.xaxis), length(vid.yaxis)),
             (X(vid.xaxis), Y(vid.yaxis)))
end

function load_videoobj(animal::String, day::Int, epoch::Int; 
        vidkws=(;), tskws=(;))::VideoObj
    file = video.getVideoFiles(animal, day)
    isempty(file) ? @error("load_videoobj :: file list is empty") : nothing
    file = file[epoch]
    vid  = load_video(file; vidkws...)
    ts   = load_videots(animal, day, epoch; tskws...)
    getVideObj(file, vid, ts; animal, day, epoch)
end
function load_videocollection(animal::String, day::Int; 
        vidkws=(;), tskws=(;))::VideoCollection
    task = DI.task.load_task(animal, day)
    epochs = unique(task.epoch)
    vids = Vector{VideoObj}([])
    for epoch in epochs
        push!(vids, video.load_videoobj(animal, Int(day), Int(epoch); vidkws, tskws))
    end
    getVideoCollection(vids)
end


# ============================================================
# ============================================================
# ============================================================
# ============================================================
# ============================================================
# ============================================================
# """
#     get_path
#
# depcrecated
# """
# function get_path(animal, day, epoch; dayfactor=0, 
#         guessdayfactor=true,
#         source="deeplabcut")
#     if guessdayfactor
#         dayfactor = raw.animal_dayfactor[animal]
#     end
#     day += dayfactor
#     if source == "deeplabcut"
#         folder_path = "/Volumes/Colliculus/deeplabcut/" * 
#                      "goalmaze_tape-Ryan-2020-05-28/videos/"
#     end
#     possible_files = 
#     glob("$(animal)_$(@sprintf("%02d",day))_$(@sprintf("%02d",epoch))_*.mp4",
#              folder_path)
#     videopath = possible_files[1]
#     return videopath
# end

