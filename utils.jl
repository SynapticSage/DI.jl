module utils

    export downsample
    export register, registerEventsToContinuous, filterAndRegister
    export filter, filterTables
    import Utils
    using Statistics, NaNStatistics
    findnearest = Utils.searchsortednearest
    using DataFrames
    using Infiltrator

    """
        downsample(animal, day)

    Downsamples the data

    Input
    =====
    raster : DataFrame
    beh : DataFrame


    Output
    ======
    (raster, behavior)

    """
    function downsample(x; dfactor=10)
        δ = dfactor;
        downsamp = x[begin:δ:end, :];
        downsamp
    end

    """
    function register(data::DataFrame...; transfer,
            on::String="time")::Vector{DataFrame} 
    """
    function register(data::DataFrame...; transfer,
            on::String="time", 
            tolerance::Union{Float64, Nothing}=0.9999,
            tolerance_violation=missing
        )::Vector{DataFrame} 
        if data isa Tuple
            data = [data...];
        end
        if tolerance === nothing
            @warn "No given tolerance"
        end
        # Get our columns into target
        @debug "→ → → → → → → → → → → → "
        @debug "Registration"
        @debug "→ → → → → → → → → → → → "
        for col ∈ transfer
            source = col[1].source
            target = col[1].target
            columns_to_transfer   = col[2]
            if columns_to_transfer == All()
                continue
            end
            @debug "columns=$columns_to_transfer from source->target on $on"

            match_on_source = data[source][:, on]
            match_on_target = data[target][:, on]

            match_on_target = convert.(Float64, match_on_target)
            match_on_source = convert.(Float64, match_on_source)

            match_on_source = (match_on_source,)
            indices_of_source_samples_in_target = findnearest.(match_on_source, match_on_target)
            if tolerance !== nothing
                δ = data[source][indices_of_source_samples_in_target, on] - data[target][:, on]
                out_of_tolerance = abs.(δ) .> tolerance
            else
                out_of_tolerance = zeros(Bool, size(data[target],1))
            end

            for item ∈ columns_to_transfer
                data[target][!, item] =
                data[source][indices_of_source_samples_in_target, item]
                if tolerance !== nothing && tolerance_violation === missing
                    data[target][!, item] = allowmissing(data[target][!,item])
                end
                if any(out_of_tolerance)
                    @info "mean out of tolerance => $(mean(out_of_tolerance))"
                    try
                        data[target][out_of_tolerance, item] .= tolerance_violation
                    catch
                        @infiltrate
                    end
                end
            end
        end
        @debug "← ← ← ← ← ← ← ← ← ← ← ← "
        return data
    end
    """
    register

    function register(source::DataFrame, target::DataFrame; 
            transfer, on::String="time")::Vector{DataFrame}

    register columns in a `source` to a `target` dataframe `on` a certain
    column
    """
    function register(source::DataFrame, target::DataFrame; kws...)::Union{Tuple, DataFrame}


        if :transfer ∈ keys(kws) && kws[:transfer] isa String
            kws=(;kws...,transfer = [kws[:transfer]])
        end

        if :transfer ∈  keys(kws) && kws[:transfer] isa Vector{String}
            addressing = (;source=1, target=2)
            kws = (;kws..., transfer=((addressing, kws[:transfer]),)) # create set of addressed transfer instructions
        end

        source, target, _ = register(source, target, DataFrame(); kws...)

        return source, target
    end

    function registerEventsToContinuous(events::DataFrame, target::DataFrame;
            transfer::Union{Vector{String},String}, on::String="time",
            targetEltype::Dict{String,<:Type}=Dict{String,Any}(),
            targetCast::Dict{String,<:Function}=Dict{String,Function}(),
            ifNonMissingAppend::Bool=false,
            eventStart::String="start",
            eventStop::String="stop")::DataFrame 

        # Get our columns into target
        @debug "→ → → → → → → → → → → → "
        @debug "Registration"
        @debug "→ → → → → → → → → → → → "
        
        for col ∈ transfer
            if col ∉ names(target)
                @debug "$col not in target"
                T = col ∈ keys(targetEltype) ? targetEltype[col] : eltype(events[!,col])
                @debug "col=$col => type=$T"
                @debug "creating new cool with type = $T"
                target[!,col] = Array{Union{Missing,T}}(missing, size(target,1))
                # Cast to user specified type?
                if col in keys(targetCast)
                    target[!,col] = targetCast[col](target[!,col])
                end
            end
        end

        match = target[!, on]
        p = Progress(size(events,1), desc="Registering")
        Threads.@threads for event in eachrow(events)
            start = event[eventStart]
            stop  = event[eventStop]
            indices_of_source_samples_in_target = (match .>= start) .&&
                                                  (match .< stop)
            #@debug "sum=$(sum(indices_of_source_samples_in_target))"
            for col ∈ transfer
                if col == All()
                    continue
                end
                specialStringBehavior = ifNonMissingAppend && 
                         eltype(target[!, col]) == Union{Missing,String} # has to be String, not InlineString
                #@debug "$col is eltype=$(eltype(target[!,col])) and specialStringBehavior=$specialStringBehavior"
                if specialStringBehavior
                    #@debug "special behavior"
                    missingVals    = indices_of_source_samples_in_target .&& ismissing.(target[!,col])
                    nonMissingVals = indices_of_source_samples_in_target .&& (!).(ismissing.(target[!,col]))
                    existing       = target[nonMissingVals, col]
                    target[nonMissingVals, col] .= existing .* event[col]
                    target[missingVals, col]    .= event[col]
                else
                    target[indices_of_source_samples_in_target, col] .= event[col]
                end
            end
            next!(p)
        end
        @debug "← ← ← ← ← ← ← ← ← ← ← ← "
        return target
    end
    """
    filter

    instructions to query/filter values in a set of dataframes
    """
    function filter(data::DataFrame...; filters::AbstractDict=Dict())::Vector{DataFrame}
        data = [data...]
        @debug "→ → → → → → → → → → → → "
        @debug "Filtration"
        @debug "→ → → → → → → → → → → → "
        for (cols, filt_for_cols) ∈ filters
            for i ∈ 1:length(data)
                @assert !(cols isa Bool)
                @assert !(filt_for_cols isa Bool)
                @assert !(filt_for_cols isa Vector{Bool})
                @assert !(cols isa Vector{Bool})
                if filt_for_cols isa Function
                    #println("Filter is a function")
                    inds = filt_for_cols(data[i][!, cols]);
                elseif typeof(filt_for_cols) <: Vector
                    #println("Filter is a set of functions")
                    inds = accumulate(.&, [ff(data[i][!, cols]) for ff
                                           in filt_for_cols])[end]
                else
                    @debug "cols = $cols"
                    @debug "Typeof(cols) = $(typeof(cols))"
                    @debug "filt_for_cols = $filt_for_cols"
                    throw(TypeError(:filt_for_cols, "",Vector,typeof(filt_for_cols)))
                end
                percent = mean(inds)*100
                @debug "data_$i filtration: $percent percent pass filter on $cols"
                x = data[i][findall(inds), :];
                data[i] = x;
            end
        end
        for i in 1:length(data)
            @debug "final_size(data($i)) = $(size(data[i]))"
        end
        @debug "← ← ← ← ← ← ← ← ← ← ← ← "
        return data
    end

    """
    `filterAndRegister`

    combination of  `raw.filter()` and `raw.register()` steps
    """
    function filterAndRegister(data::DataFrame...;
            filters::Union{Nothing,AbstractDict}=nothing, transfer=nothing,
            on="time")::Vector{DataFrame}
        if transfer != nothing
            data = register(data...; transfer=transfer, on=on)
            #utils.piso(data)
        end
        if filters != nothing
            data = utils.filter(data...; filters=filters)
        end
        return data    
    end

    """
    `filterTables`

    alias for `filterAndRegister`
    """
    filterTables(data::DataFrame...; filters::Union{Nothing,AbstractDict}=nothing,
                 lookupcols=nothing, lookupon="time") =
    filterAndRegister(data...;transfer=lookupcols, on=lookupon, filters=filters)

end
