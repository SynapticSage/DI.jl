module utils

    export downsample
    export register, registerEventsToContinuous, filterAndRegister
    export filter, filterTables
    import Utils
    import Filt
    using Statistics, NaNStatistics
    findnearest = Utils.searchsortednearest
    using DataFrames
    import DataFrames: ColumnIndex
    using Infiltrator

    CItype = Union{ColumnIndex, Vector{<:ColumnIndex}}

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
        register

    function register(source::DataFrame, target::DataFrame; 
            transfer, on::String="time")::Vector{DataFrame}

    register columns in a `source` to a `target` dataframe `on` a certain

    ### Input
    `source` -- the source dataframe to register
    `target` -- the recipient dataframe

    ### Output

    column
    """
    function register(source::DataFrame, target::DataFrame; 
            transfer::Vector{<:CItype},
            on::CItype="time",
            tolerance::Union{Float64, Nothing}=0.9999,
            tolerance_violation=missing,
            convert_type::Type=Float64
        )::Vector{DataFrame} 
        @info "Did IT!"
        if tolerance === nothing
            @warn "No given tolerance"
        end
        # Get our columns into target
        @debug "→ → → → → → → → → → → → "
        @debug "Registration"
        @debug "→ → → → → → → → → → → → "

        if transfer == All()
            return source, target
        end
        @debug "columns=$transfer from source->target on $on"

        match_on_source = source[:, on]
        match_on_target = target[:, on]

        match_on_target = typeof(match_on_target) == convert_type ?
                          match_on_target : convert(Vector{convert_type},
                                                    match_on_target)
        match_on_source = typeof(match_on_target) == convert_type ?
                          match_on_source : convert(Vector{convert_type},
                                                    match_on_source)

        # FindNearest
        match_on_source = (match_on_source,)
        indices_of_source_samples_in_target = findnearest.(match_on_source,
                                                           match_on_target)

        # Tolerance
        if tolerance !== nothing
            δ = source[indices_of_source_samples_in_target, on] -
                target[:, on]
            out_of_tolerance = abs.(δ) .> tolerance
        else
            out_of_tolerance = zeros(Bool, size(target,1))
        end

        @debug any(out_of_tolerance) ? 
        "mean out of tolerance=>$(mean(out_of_tolerance))" :
        "no out of tolerance"

        @infiltrate

        # Move columns
        target[!, transfer] =
        source[indices_of_source_samples_in_target, transfer]

        # If tolerance vio is nothing, drop out of tolerance entries
        if tolerance_violation === nothing
            target = target[Not(out_of_tolerance), :]
            indices_of_source_samples_in_target = 
                indices_of_source_samples_in_target[Not(out_of_tolerance)]
        elseif tolerance_violation === missing
            for item ∈ transfer
                target[!, item] = allowmissing(target[!,item])
            end
        end

        # Write in tolerance violations
        if  tolerance !== nothing &&
            tolerance_violation !== nothing && 
            any(out_of_tolerance)
            try
                target[out_of_tolerance, transfer] .= tolerance_violation
            catch
                @infiltrate
            end
        end

        @debug "← ← ← ← ← ← ← ← ← ← ← ← "
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
    function filter(data::DataFrame...; filters::AbstractDict=Dict(),
            filter_skipmissingcols::Bool=false,
        )::Vector{DataFrame}
        data = [data...]
        @debug "→ → → → → → → → → → → → "
        @debug "Filtration"
        @debug "→ → → → → → → → → → → → "
        @inbounds for (cols, filt_for_cols) ∈ filters
            for i ∈ 1:length(data)
                @assert !(cols isa Bool)
                @assert !(filt_for_cols isa Bool)
                @assert !(filt_for_cols isa Vector{Bool})
                @assert !(cols isa Vector{Bool})
                list_cols = cols isa Vector ? String.(cols) : String.([cols])
                if filter_skipmissingcols && any(c ∉ names(data[i]) for c ∈ list_cols)
                    @debug "data[$i] missing col ∈ cols=$cols, skip filter"
                    continue
                end
                if filt_for_cols isa Function
                    #println("Filter is a function")
                    inds = filt_for_cols(data[i][!, cols])
                elseif typeof(filt_for_cols) <: Vector
                    #println("Filter is a set of functions")
                    inds = @fastmath accumulate(.&, [ff(data[i][!, cols]) for ff
                                           in filt_for_cols])[end]
                else
                    @debug "cols = $cols"
                    @debug "Typeof(cols) = $(typeof(cols))"
                    @debug "filt_for_cols = $filt_for_cols"
                    throw(TypeError(:filt_for_cols, "",Vector,typeof(filt_for_cols)))
                end
                inds = replace(inds, missing=>false)
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
        filterAndRegister

    combination of  `raw.filter()` and `raw.register()` steps
    """
    function filterAndRegister(data::DataFrame...;
            filters::Union{Nothing,AbstractDict}=nothing, transfer=nothing,
            on="time", filter_skipmissingcols::Bool=false
            )::Vector{DataFrame}
        if filters !== nothing
            required_cols  = Filt.get_filter_req(filters)
            missing_fields = Vector{Vector{String}}(undef, length(data)-1)
            for i in 2:length(data)
                missing_fields[i-1] = setdiff(required_cols, names(data[i]))
            end
            if !(isempty(missing_fields))
                @debug "Adding missing_fields=$missing_fields to transfer=$transfer"
                transfer = unique(collect(Iterators.flatten([transfer, 
                                 missing_fields...])))
            end
        end
        if transfer !== nothing
            data = register(data...; transfer=transfer, on=on)
            #utils.piso(data)
        end
        if filters !== nothing
            data = utils.filter(data...; filters=filters, filter_skipmissingcols)
        end
        return collect(data)
    end

    """
    `filterTables`

    alias for `filterAndRegister`
    """
    filterTables(data::DataFrame...; filters::Union{Nothing,AbstractDict}=nothing,
                 lookupcols=nothing, lookupon="time") =
    filterAndRegister(data...;transfer=lookupcols, on=lookupon, filters=filters)

end
