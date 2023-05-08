import DIutils
using DrWatson, DataFrames, Infiltrator, StatsBase, ImageFiltering
using RollingFunctions, DataFramesMeta
export save_behavior, load_behavior, behaviorpath

move_thresh = 4 # cm/s

function behaviorpath(animal::String, day::Int, tag::String=""; type::String=DI.load_default)
    tag  = length(tag) == 0 ? tag : "_$tag"
    path = datadir("exp_raw", "visualize_raw_neural",
                     "$(animal)_$(day)_beh$tag.$type")
    if occursin("*",tag)
        path = glob(basename(path), dirname(path))
    end
    return path
end

"""
    load_behavior(animal::String, day::Int, tag::String=""; type::String=DI.load_default, kws...)

Load the behavior data for a given animal and day.
# Arguments
- `animal::String`: animal name
- `day::Int`: day number
- `tag::String="": tag to append to the filename
- `type::String=DI.load_default`: type of file to load. Default is to use the default
    file type for the given animal and day.
- `kws...`: keyword arguments to pass to `DI.load_table`
# Returns
- `beh::DataFrame`: behavior data
"""
function load_behavior(animal::String, day::Int, tag::String="";
                      type::String=DI.load_default, quick=false, kws...)

    function typeFunc(type, name)
        if occursin("Vec", string(name))
            type = ComplexF32;
        elseif name == "time" type = Float32;
        else
            type = nothing;
        end
        return type
    end
    if type == "csv"
        typemap = Dict(Int64=>Int16);

        load_kws = (;strict=false, missingstring=["NaNNaNi", "NaNNaNi,", ""], types=typeFunc, typemap=typemap, DI.csvkws...)
    else
        load_kws = (;)
    end
    beh = DI.load_table(animal, day, tag; tablepath=:behavior, type=type, 
                        load_kws=load_kws, kws...)
    if type == "csv"
        if beh.time isa Vector{String}
            beh.time = parse.(Float32, beh.time);
        end
        for col in names(beh)
            if (occursin("current", col) || occursin("egoVec", col)) &&
                eltype(typeof(beh[!,col])) <: Real
                try
                    replace!(beh[!,col], missing=>NaN)
                catch; @infiltrate; end
            end
        end
        @assert ("x" ∈ names(beh)) "Fuck"
    end
    if !quick && !startswith(animal, "super")
        postprocess!(beh)
    end
    return beh
end
function load_behavior(D::AbstractDict; kws...)
    pos = []
    @assert "animal" in keys(D)
    @assert "day" in keys(D)
    for k in ["animal", "day", "tag"]
        if k in keys(D)
            push!(pos, v)
        end
    end
    load_behavior(pos...; kws...)
end

"""
    save_behavior()
"""
function save_behavior(behavior::AbstractDataFrame, pos...; kws...)
    DI.save_table(behavior, pos...; tablepath=:behavior, kws...)
end


# ----------------------------------------------------------------------
# ADDITIONAL POST-PROCESSING FUNCTIONS
# 
# It's arguable these steps should happen before the data is imported to
# julia (my preprocessing pipeline is in matlab, because it uses our labs code base). 
# But for now the plan is to use these post-processing bandages.
# ----------------------------------------------------------------------

function postprocess!(beh::AbstractDataFrame)::Nothing
    register_epoch_homewell!(beh)
    annotate_hatraj!(beh)
    annotate_poke!(beh)
    munge_clean_velocity!(beh)
    immobility!(beh)
    nothing
end

function determine_epoch_homewell(beh::AbstractDataFrame)::AbstractDataFrame
    B = groupby(subset(dropmissing(beh,[:startWell,:stopWell]),
                       :startWell=>w->w .!= -1, :stopWell=>w->w .!=-1),
                :epoch)
    hws = combine(B, :stopWell=>mode, :startWell=>mode)
    @assert(all(hws.stopWell_mode .== hws.startWell_mode))
    hws.homewell = hws.stopWell_mode
    hws[!, Not([:startWell_mode, :stopWell_mode])]
end

function register_epoch_homewell!(beh::AbstractDataFrame)::DataFrame
    hws = determine_epoch_homewell(beh)
    DIutils.filtreg.register(hws, beh; on="epoch", transfer=["homewell"])
    beh
end

function annotate_hatraj!(beh::AbstractDataFrame)::Nothing
    beh[!,:hatraj]    = Vector{Union{String,Missing}}(missing, size(beh,1))
    beh[!,:hatrajnum] = Vector{Union{Int8,Missing}}(missing, size(beh,1))
    beh[!,:ha] = Vector{Union{Char,Missing}}(missing, size(beh,1))

    inds = transform(beh, 
                     :block => (b->
                                (!).(ismissing.(b)) .&& 
                                ((!).(isnan.(b))) .&& 
                                (b .!= -1))
                                 => :out).out

    B = groupby(beh[inds,:], [:epoch, :block])
    for block in B

        correct_trials = block.correct .== 1
        correct_block_traj = block[correct_trials, :]
        incorrect_block = block[(!).(correct_trials), :]

        correct_block = _assign_labels(correct_block_traj)

        # Assign incorrect block labels
        if !isempty(correct_block.time)
            closest_cor_samp = DIutils.searchsortedprevious.([correct_block.time], incorrect_block.time)
            incorrect_block.hatraj = "*" .* correct_block[closest_cor_samp,:hatraj]
        elseif isempty(correct_block) && !isempty(incorrect_block)
            incorrect_block = _assign_labels(incorrect_block)
            incorrect_block.hatraj = "*" .* incorrect_block.hatraj
        end
        
        # Assign values back to block
        if !isempty(correct_block) || !isempty(incorrect_block)
            block[:,:] .= sort(vcat(correct_block, incorrect_block), :time)
        end
        #@infiltrate length(unique(block.ha)) > 1

    end

    B = sort(combine(B,identity), [:epoch,:time])
    beh[inds,:hatraj] = B.hatraj
    beh[inds,:ha] = B.ha

    beh[inds,:hatrajnum] = DIutils.searchsortednearest.(
                    [sort(unique(beh.hatraj[inds]))], beh.hatraj[inds])

    nothing
end

"""
    _assign_labels

assigns `hatraj` and `hatrajnum` to a block of behavior
"""
function _assign_labels(block)
    home_trials  = block.stopWell .== block.homewell
    arena_trials = (!).(home_trials)

    if any(home_trials)
        uareantraj = sort(unique(block.traj[home_trials]))
        hometraj  = DIutils.searchsortednearest.([uareantraj], 
                                               block.traj[home_trials])
    else
        hometraj = []
    end

    if any(arena_trials)
        uareantraj = sort(unique(block.traj[arena_trials]))
        arenatraj  = DIutils.searchsortednearest.([uareantraj], 
                                               block.traj[arena_trials])
    else
        arenatraj = []
    end

    home_labels  = isempty(hometraj)   ? Vector{String}() : "H" .* string.(hometraj)
    arena_labels = isempty(arenatraj)  ? Vector{String}() : "A" .* string.(arenatraj)

    # Assign incorrect block labels
    block.hatraj[home_trials]  .= home_labels
    block.hatraj[arena_trials] .= arena_labels
    block.ha[home_trials]  .= 'H'
    block.ha[arena_trials] .= 'A'
    #@infiltrate length(home_trials) > 1

    return block
end


"""
    clean_velocity

Cleans velocity and adds some additional fields


# TODO
`NaN` cleaning ... by eye it appears much of the time NaN velocity
is 0 velocity. It may be the case that I should be using a more sophisticated
NaN fill procedure. e.g. ffill/bfill

"""
function munge_clean_velocity!(beh)
    R = real.(beh.velVec)
    I = imag.(beh.velVec)
    R[isnan.(R)] .= 0
    I[isnan.(I)] .= 0
    beh[!,:velVec] = R .+ (I)im

    ker = Kernel.gaussian((2,))
    #Plots.unicodeplots()
    #Plots.plot(ker, label="smoothing kernel", 
    #                title="length of time: $(round(length(ker)/30,sigdigits=2)) sec")
    beh[!,:smoothvel] = imfilter(beh.velVec, ker)
    @info "NOT moving defined as smoothvel < 4cm/s || dio poke"
    stillness = abs.(beh.smoothvel) .< move_thresh .|| (beh.poke .> 0)
    beh[!,:moving] = (!).(stillness)
    nothing
end

"""
    annotate_poke!

adds a `poke` field describing which DIO is poked using the individual
poke_X fields
"""
function annotate_poke!(beh::AbstractDataFrame; manualpokefix::Bool=false)::Nothing
    pn = sort([name for name in names(beh) if occursin("poke_", name)])
    poke_matrix = replace(Matrix(beh[!, pn]), NaN=>0)
    valid = DIutils.squeeze(any((!).(Matrix(ismissing.(beh[!, pn]))), dims=2))
    poke_matrix = BitMatrix(poke_matrix[valid,1:end])
    pn = replace([findfirst(row) for row in eachrow(poke_matrix)],nothing=>0)
    #locs = accumulate(|, [(!).(ismissing.(beh[!,col])) for col in pn])
    beh[!,:poke] = Vector{Union{Missing,Int8}}(missing, size(beh,1))
    if manualpokefix
        @warn("annotate_poke! ↘" *
              "`manualpokefix==true :: switching poke 1 and 2. "*
              "the poke metadata switched earlier in the pipeline..."*
              "this is a temporary bandage for that mistake"
             )
        tmp = copy(pn[:, 1])
        pn[:, 1] = pn[:, 2]
        pn[:, 2] = tmp
    end
    beh[valid, :poke] = pn
    nothing
end

# export debounce!
# """
#     debounce!(df, event; threshold=0.1)
# debounce the event in the dataframe `df` by `threshold` seconds.
# # Arguments
# - `df`: a dataframe with a `:time` column
# - `event`: the name of the event column
# - `threshold`: the minimum time between events
# # Returns
# - `df`: the dataframe with the debounced event
# """
# function debounce!(df::AbstractDataFrame, event; threshold=0.1)
#     dft = groupby(df, :time)
#     lastup, lastdown, isdown = [], [], false
#     (k, df) = first(zip(keys(dft), dft))
#     for (k, df) in zip(keys(dft), dft)
#         k    = k[1] # current time
#         isup = df[!,event][1] # current event
#         # store the current event
#         if isup == true
#             push!(lastup, k)
#         else
#             push!(lastdown, k)
#         end
#         if length(lastup) > 2 && length(lastdown) > 2
#             # if the event is down and the last two down events are
#             # closer than the threshold, then set all events up to
#             # the last down event to false
#             if !isup && (lastdown[end-1] - lastdown[end-2]) < threshold
#                 inds = lastup .>= lastdown[end-2]
#                 for k in lastup[inds]
#                     dft[(;time=k,)][!,event] .= false
#                 end
#                 # remove all the events that were affected
#                 lastup = lastup[.!inds]
#             # if the event is up and the last two up events are
#             # closer than the threshold, then set all events down to
#             # the last up event to true
#             elseif isup && (lastup[end-1] - lastup[end-2]) < threshold
#                 inds = lastdown .>= lastup[end-2]
#                 for k in lastdown[inds]
#                     dft[(;time=k,)][!,event] .= true
#                 end
#                 # remove all the events that were affected
#                 lastdown = lastdown[.!inds]
#             end
#         end
#     end
#     return combine(dft, identity)
# end
#
# export debounce_movement!
# function debounce_movement!(beh::DataFrame; threshold=0.1, 
#     replace_moving=true)
#     # Debounce the movement
#     println("Debouncing movement, fraction moving: ", 
#         mean(beh[:, :moving]))
#     beh[!,:movingdeb] = copy(beh[:, :moving])
#     println("About to debounce -- takes about a minute")
#     prog = Progress(length(groupby(beh, :epoch)); desc="debouncing")
#     epochs = groupby(beh, :epoch) |> collect
#     iters = enumerate(epochs) |> collect
#     Threads.@threads for (i,epoch) in iters
#         epochs[i] = debounce!(epoch, :movingdeb; threshold)
#         next!(prog)
#     end
#     beh = vcat(epochs...)
#     println("fraction before debounce: ", mean(beh[:, :moving]))
#     println("fraction moving after debounce: ", 
#         mean(beh[:, :movingdeb]))
#     @assert(!all(beh.moving .== beh.movingdeb),
#         "movingdeb should not be the same as moving")
#     beh.moving = beh[:, :movingdeb]
#     beh
# end

"""
    immobility!(beh::DataFrame; win=10, move_thresh=move_thresh)
smooth the speed of the animal using a rolling mean of `win` frames
and define the animal as moving if the speed is greater than
`move_thresh`. 
# Arguments
- `beh`: a dataframe with a `:speed` column
- `win`: the window size for the rolling mean
- `move_thresh`: the threshold for the speed to be considered moving
# Returns
- `beh`: the dataframe with the `:moving` column
"""
function immobility!(beh::AbstractDataFrame; win=10, 
    move_thresh=move_thresh, movebool_gaussian=1)::DataFrame

    rollit(x,win)  = [zeros(Int(win/2-1)); rollmean(x, win); 
                      zeros(Int(win/2))];
    beh[!,:speedsmooth] = rollit(beh.speed,win)
    # begin 
    #     tmp=@subset(beh, :traj .== 2) 
    #     @df tmp plot(:time, :speed)
    #     @df tmp plot!(:time, :speedsmooth)
    #     plot!(xlabel="Time (s)", ylabel="Speed (cm/s)")
    #     ylims!(0,10)
    # end
    stillness = (ismissing(beh.poke) .|| 
                (.!(ismissing.(beh.poke)) .&& (beh.poke .> 0))) .&&
                abs.(beh.speedsmooth) .< move_thresh
    replace!(stillness, missing=>true)
    beh[!,:moving_speedsmooth] = .!stillness
    beh[!,:movingsmooth] = beh[:, :moving_speedsmooth]
    try
        beh[!,:movingsmooth] = imfilter(beh.moving_speedsmooth,
                                Kernel.gaussian((movebool_gaussian*10,)))
    catch
        @infiltrate
    end
    B = groupby(
        @subset(beh, :traj .!= NaN, :startWell .!= -1, :stopWell .!= -1; 
                view=true), :traj)
    beh[!,:immobility] = beh.movingsmooth .< 0.5
    beh
end
