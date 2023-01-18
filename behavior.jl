module behavior

    import ..Load, Utils
    using DrWatson
    using DataFrames
    using Infiltrator
    using StatsBase
    export save_behavior, load_behavior, behaviorpath

    function behaviorpath(animal::String, day::Int, tag::String=""; type::String=Load.load_default)
        tag  = length(tag) == 0 ? tag : "_$tag"
        path = datadir("exp_raw", "visualize_raw_neural",
                         "$(animal)_$(day)_beh$tag.$type")
        if occursin("*",tag)
            path = glob(basename(path), dirname(path))
        end
        return path
    end

    function load_behavior(animal::String, day::Int, tag::String="";

        type::String=Load.load_default, kws...)

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

            load_kws = (;strict=false, missingstring=["NaNNaNi", "NaNNaNi,", ""], types=typeFunc, typemap=typemap, Load.csvkws...)
        else
            load_kws = (;)
        end
        beh = Load.load_table(animal, day, tag; tablepath=:behavior, type=type, 
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
        postprocess!(beh)
        return beh
    end

    function save_behavior(behavior::AbstractDataFrame, pos...; kws...)
        Load.save_table(behavior, pos...; tablepath=:behavior, kws...)
    end

    function postprocess!(beh::DataFrame)::Nothing
        register_epoch_homewell!(beh)
        annotate_hatraj!(beh)
        nothing
    end

    function determine_epoch_homewell(beh::DataFrame)::DataFrame
        B = groupby(subset(dropmissing(beh,[:startWell,:stopWell]),
                           :startWell=>w->w .!= -1, :stopWell=>w->w .!=-1),
                    :epoch)
        hws = combine(B, :stopWell=>mode, :startWell=>mode)
        @assert(all(hws.stopWell_mode .== hws.startWell_mode))
        hws.homewell = hws.stopWell_mode
        hws[!, Not([:startWell_mode, :stopWell_mode])]
    end

    function register_epoch_homewell!(beh::DataFrame)::DataFrame
        hws = determine_epoch_homewell(beh)
        Utils.filtreg.register(hws, beh; on="epoch", transfer=["homewell"])
        beh
    end

    function annotate_hatraj!(beh::DataFrame)::Nothing
        beh[!,:hatraj]    = Vector{Union{String,Missing}}(missing, size(beh,1))
        beh[!,:hatrajnum] = Vector{Union{Int8,Missing}}(missing, size(beh,1))

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
                closest_cor_samp = Utils.searchsortedprevious.([correct_block.time], incorrect_block.time)
                incorrect_block.hatraj = "*" .* correct_block[closest_cor_samp,:hatraj]
            elseif isempty(correct_block) && !isempty(incorrect_block)
                incorrect_block = _assign_labels(incorrect_block)
                incorrect_block.hatraj = "*" .* incorrect_block.hatraj
            end
            
            # Assign values back to block
            if !isempty(correct_block) || !isempty(incorrect_block)
                block[:,:] .= sort(vcat(correct_block, incorrect_block), :time)
            end

        end

        beh[inds,:hatraj] = sort(combine(B,identity), [:epoch,:time]).hatraj

        beh[inds,:hatrajnum] = Utils.searchsortednearest.(
                        [sort(unique(beh.hatraj[inds]))], beh.hatraj[inds])

        nothing
    end

    function _assign_labels(block)
        home_trials  = block.stopWell .== block.homewell
        arena_trials = (!).(home_trials)

        if any(home_trials)
            uareantraj = sort(unique(block.traj[home_trials]))
            hometraj  = Utils.searchsortednearest.([uareantraj], 
                                                   block.traj[home_trials])
        else
            hometraj = []
        end

        if any(arena_trials)
            uareantraj = sort(unique(block.traj[arena_trials]))
            arenatraj  = Utils.searchsortednearest.([uareantraj], 
                                                   block.traj[arena_trials])
        else
            arenatraj = []
        end

        home_labels  = isempty(hometraj)   ? Vector{String}() : "h" .* string.(hometraj)
        arena_labels = isempty(arenatraj)  ? Vector{String}() : "a" .* string.(arenatraj)

        # Assign incorrect block labels
        block.hatraj[home_trials]  .= home_labels
        block.hatraj[arena_trials] .= arena_labels

        return block
    end

end
