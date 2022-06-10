module raw
    
    #Imports
    using Revise
    using DrWatson
    using DataFrames
    import CSV, Arrow
    using Statistics
    using Dates
    using Printf
    using ProgressMeter
    using Glob
    using Infiltrator
    using Colors

    # Fetch maze internal imports
    include("utils.jl")
    import .utils
    include("table.jl")
    import .table
    using .table: get_periods
    findnearest = utils.searchsortednearest
    __revise_mode__ = :eval

    load_default="arrow"

    components = [
        "./raw/behavior.jl",
         "./raw/lfp.jl",
         "./raw/task.jl",
         "./raw/decode.jl",
         "./raw/celltet.jl",
         "./raw/spikes.jl",
         "./raw/ripples.jl",
         "./raw/dlc.jl",
         "./raw/utils.jl"
    ]
    for comp in components
        comp = replace(comp,"./"=>"$(srcdir())/")
        println("$comp")
        include(comp)
    end

    load_functions = Dict(
        "cycles"   => load_cycles,
        "behavior" => load_behavior,
        "spikes"   => load_spikes,
        "lfp"      => load_lfp,
        "cells"    => load_cells,
        "tetrode"  => load_tetrode,
        "decode"   => load_decode,
        "task"     => load_task,
        "ripples"  => load_ripples
       )

    save_functions = Dict(
        "behavior" => save_behavior,
        "spikes"   => save_spikes,
        "lfp"      => save_lfp,
        "cells"    => save_cells,
        "task"     => save_task,
        "ripples"  => save_ripples
       )

    path_functions = Dict(
        "cycles"   => cyclepath,
        "behavior" => behaviorpath,
        "spikes"   => spikespath,
        "lfp"      => lfppath,
        "cells"    => cellpath,
        "task"     => taskpath,
        "ripples"  => ripplespath,
    )

    # Module-wide settings
    animal_dayfactor = Dict("RY16"=>33, "RY22"=>0)
    csvkws=(; silencewarnings=true, buffer_in_memory=true, ntasks=1, 
            strict=false, missingstring=["NaN", "", "NaNNaNi"])
    arrowkws = (;)
    pxtocm(x) = 0.1487 * x
    cmtopx(x) = x / 0.1487 
    default_set = ["spikes", "behavior", "ripples", "cells"]

    function normalize(data, time, data_source, mintime=nothing)
        # Determine a time normalizing function
        if mintime == nothing
            if "behavior" ∈ data_source
                mintime = minimum(data["behavior"].time)
            else
                mintime = 0
            end
        end

        (time .- mintime)./60, mintime
    end


    """
        load(animal, day)

    Loads up the datasets we plan to use to plot out the raw raster data

    =====
    Input
    =====
    animal : String
        Name of the animal
    day : Int
        Integer of the day


    ======
    Output
    ======
    (raster, behavior)

    """
    function load(args...; as="tuple",
            data_source=default_set)

        # Establish a load order
        if as == "dict"
            load_order = sort(data_source, by=source->source=="behavior", rev=true)
        elseif as == "tuple"
            if "behavior" in data_source
                load_order = ("behavior", setdiff(data_source, ["behavior"])...)
            end
        else
            throw(ArgumentError("as cannot be $as"))
        end
        

        # Load each data source
        data = Dict{String,Any}()
        m = nothing
        for source ∈ load_order
            data[source] = raw.load_functions[source](args...)
            @assert typeof(data[source]) != Nothing
            if "time" ∈ names(data[source])
                @info "Centering $source and **converting to minutes**"
                data[source].time, m = normalize(data, data[source].time, data_source, m);
            end
        end

        # How do we return, dict or tuple?
        if as == "tuple"
            return (data[name] for name in data_source)
        else
            return data
        end
        
    end

    function get_load_kws(type::String; load_kws...)
        if type == "csv"
            kws = (csvkws..., load_kws...)
        elseif type == "arrow"
            kws = (arrowkws..., load_kws...)
        end
    end
    function load_table_at_path(path::String, type::String; load_kws...)
        if type == "csv"
            data = CSV.read(path, DataFrame; load_kws...)
        elseif type == "arrow"
            data = DataFrame(Arrow.Table(path))
        end
    end
    function load_table(animal::String, day::Int, pos...; 
            tablepath=nothing, 
            load_kws::Union{Nothing,NamedTuple}=(;),
            kws...)
        if tablepath == nothing
            throw(ArgumentError("Must provide tablepath symbol... see path_functions dict in this module"))
        end

        tablepath = String(tablepath)
        path = path_functions[tablepath](animal, day, pos...; kws...)

        if :type in keys(kws)
            type = kws[:type]
        else
            type = load_default
        end

        data = load_table_at_path(path, type; load_kws)

        detect_complex_format_wrong = eltype.(eachcol(data)) .<: NamedTuple
        for i in findall(detect_complex_format_wrong)
            data[!,i] = getproperty.(data[!,i], :re) + 
                        getproperty.(data[!,i], :im)im
        end

        data
    end

    function save_table_at_path(data::AbstractDataFrame, path::String, type::String;
            save_kws...)
        if type == "csv"
            data |> CSV.write(path)
        elseif type == "arrow"
            Arrow.write(path, data)
        end
    end
    function save_table(data::AbstractDataFrame, pos...; tablepath=nothing, kws...)
        if tablepath == nothing
            throw(ArgumentError("Must provide tablepath symbol... see path_functions dict in this module"))
        end
        tablepath = String(tablepath)
        path = path_functions[tablepath](pos...; kws...)
        println("Saving $(tablepath) data at $path")
        if :type in keys(kws)
            type = kws[:type]
        else
            type = "csv"
        end
        save_table_at_path(data, path, type)
    end

    function tables_to_type(animal::String, day::Int; 
            inclusion=keys(save_functions),
            exclusion=[:lfp], from::String="csv", to::String="arrow",
            delete_prev::Bool=false)
        exclusion = String.(exclusion)
        @showprogress for key in inclusion
            key = String(key)
            if key ∈ exclusion
                continue
            end
            loadfunc, savefunc = load_functions[key], save_functions[key]
            data = loadfunc(animal, day)
            savefunc(data, animal, day; type=to)
            if delete_prev
                pathfunc = path_functions[key]
                path = pathfunc(animal, day; type=type)
                if isfile(path); rm(path); end
            end
        end
    end
    
    # ----------
    # OPERATIONS
    # ----------


    """
    Takes a set of dataframes and dicts and normalizes time variables of each
    via the normalize() function. 

    Timefields controls which fields to normalize (defaults to "time"), per
    datatype passed in. Right now each datatype gets a number 1 through N, and
    timefields is a dict from that number to the correct set of fields.


    """
    function normalize_time(data::Union{DataFrame, Dict}...;
            timefields::Dict=Dict(), factor=1, warnifnotime::Bool=false)
        data = [data...]
        if data[1] isa DataFrame
            tₘ = minimum(data[1].time)
        else
            tₘ = minimum(data[1]["time"])
        end
        for source ∈ 1:length(data)
            @debug "source = $source"
            if !(source in keys(timefields))
                tfs = ["time"]
                @debug "tfs = $tfs"
            else
                tfs = timefields[source]
                @debug "tfs = $tfs"
            end
            for tf in tfs
                if data[source] isa DataFrame && (tf ∈ names(data[source]))
                    data[source][!,tf] = (data[source][!,tf] .- tₘ)*factor
                elseif !(data[source] isa DataFrame) && (tf ∈ keys(data[source]))
                    data[source][tf] = (data[source][tf] .- tₘ)*factor
                else
                    if warnifnotime
                        @warn "No time_field=$tf for $source"
                    else
                        @error "No time_field=$tf for $source"
                    end
                end
            end
        end
        return data
    end

    function keep_overlapping_times(data::Union{DataFrame, Dict, Vector}...; 
            tf="time", returninds::Vector=[])

        data=[data...]
        
        # Get time extrema
        sz = []
        for d in data
            if d isa DataFrame && (tf ∈ names(d))
                push!(sz, extrema(d[!,tf]))
            elseif !(d isa DataFrame) && (tf ∈ keys(d))
                push!(sz, extrema(d[tf]))
            elseif d isa Vector
                push!(sz, extrema(d)) 
            end
        end

        # Compute overall extrema
        sz = cat([[s...] for s in sz]...; dims=2)'
        minmax = [maximum(sz[:,1]), minimum(sz[:,2])]


        # Contrain each object to live in the overall extrema
        ind_constrain(x) = (x .>= minmax[1]) .&& (x .< minmax[2])
        for (i, d) in zip(1:length(data), data)
            if d isa DataFrame && (tf ∈ names(d))
                inds = ind_constrain(d[!,tf])
            elseif !(d isa DataFrame) && (tf ∈ keys(d))
                inds = ind_constrain(d[tf]) 
            elseif d isa Vector
                inds = ind_constrain(d) 
            end

            
            if i in returninds
                data[i] = inds
            else
                if d isa DataFrame && (tf ∈ names(d))
                    data[i] = d[inds, :]
                elseif !(d isa DataFrame) && (tf ∈ keys(d))
                    time_length = length(d[tf])
                    dims = Dict(key=>findfirst(size(value).==time_length)
                                for (key,value) in d)
                    for (k,v) in d
                        I = Vector{Any}([Colon() for i in 1:ndims(v)])
                        if !(isempty(dims[k]))
                            I[dims[k]] = inds
                        end
                        data[i][k] = getindex(v, I...)
                    end
                elseif d isa Vector
                    data[i] = d[inds]
                end
            end
        end

        return data
    end

    function fix_complex(df::DataFrame)
        for (i, col) in enumerate(eachcol(df))
            if typeof(col) <: Arrow.Struct || 
                eltype(col) <: NamedTuple
                try
                    df[!,i] = fix_complex(df[!,i])
                catch
                    @infiltrate
                end
            end
        end
        return df
    end
    fix_complex(x::AbstractVector) = fix_complex.(x)
    fix_complex(x::NamedTuple)     = x.re + (x.im)im

    function fix_rgba(df::DataFrame)
        for (i, col) in enumerate(eachcol(df))
            if typeof(col) <: Arrow.Struct || 
                eltype(col) <: NamedTuple
                try
                    df[!,i] = fix_rgba(df[!,i])
                catch
                    @infiltrate
                end
            end
        end
        return df
    end
    fix_rgba(x::AbstractVector) = fix_rgba.(x)
    fix_rgba(x::NamedTuple)     = RGBA(x.r, x.g, x.b, x.alpha)



end
