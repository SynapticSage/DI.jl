__precompile__()
module DI
    
    using Revise
    using DrWatson

    using DataFrames
    import CSV, Arrow
    using Statistics, Dates, Printf, ProgressMeter, Glob, Colors,
          Infiltrator

    # Fetch maze internal imports
    using DIutils
    import DIutils.Table: get_periods
    __revise_mode__ = :eval

    load_default="arrow"

    include("behavior.jl")
    include("celltet.jl")
    include("decode.jl")
    include("dlc.jl")
    include("lfp.jl")
    include("path.jl")
    include("spikes.jl")
    include("ripples.jl")
    include("task.jl")
    include("utils.jl")
    include("video.jl")
    include("superanimal.jl")
    include("Filt.jl")
    include("Labels.jl")

    load_functions = Dict(
        "cycles"   => load_cycles,
        "behavior" => load_behavior,
        "beh"      => load_behavior,
        "spikes"   => load_spikes,
        "lfp"      => load_lfp,
        "cells"    => load_cells,
        "tetrode"  => load_tetrode,
        "decode"   => load_decode,
        "task"     => load_task,
        "ripples"  => load_ripples,
        "coh"      => load_coh,
        "avgcoh"   => load_avgcoh,
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
        "beh"      => behaviorpath,
        "spikes"   => spikespath,
        "lfp"      => lfppath,
        "cells"    => cellpath,
        "task"     => taskpath,
        "ripples"  => ripplespath,
        "coh"      => cohpath,
        "avgcoh"   => avgcohpath,
    )

    time_vars = Dict(
        "cycles"   => [:time, :start, :stop],
        "spikes"   => :time,
        "behavior" => :time,
        "beh"      => :time,
        "lfp"      => :time,
        "coh"      => :time,
        "avg"      => :time,
        "cells"    => [],
        "task"     => [],
        "ripples"  => ["time", "start", "stop"],
        "cycles"   => ["time", "start", "stop"]
    )

    # Module-wide settings
    default_args = ("super",0) # default animal to loadup
    animal_set = (("RY22", 21), ("RY16",36))
    animal_dayfactor = Dict("RY16"=>33, "RY22"=>0)
    csvkws=(; silencewarnings=true, buffer_in_memory=true, ntasks=1, 
            strict=false, missingstring=["NaN", "", "NaNNaNi"])
    arrowkws = (;)
    pxtocm(x) = 0.1487 * x
    cmtopx(x) = x / 0.1487 
    default_set = ["spikes", "behavior", "ripples", "cells"]
    total_set = ["spikes", "behavior", "ripples", "cells", "task"]
    min_time_records = [] # variable records minimum time of the most recent load
    function _set_mintime!(m::Real)
        push!(min_time_records, m)
    end

    """
        normalize(data, time, data_source, mintime=nothing)
    
    Normalizes the time variable to start at 0, or at mintime if passed
    
    # Arguments
    - data: Dict of data sources (dataframes of different types of data)
    - time: time variable to normalize
    - data_source: which data source is being normalized
    - mintime: optional, if passed, time is normalized to start at mintime
    """
    function normalize(data, time, data_source, mintime=nothing)
        # Determine a time normalizing function
        if mintime === nothing
            if "behavior" ∈ data_source
                mintime = minimum(data["behavior"].time)
            else
                mintime = 0
            end
        end

        time .- mintime, mintime
    end


    """
        load(animal, day; data_source=default_set, center=true)

    Loads up the datasets we plan to use to plot out the raw raster data

    # Input
    - `animal` - animal name
    - `day` - day number
    - `args...` - additional arguments to pass to the load functions
    - `as` - whether to return a tuple or a dict, default is tuple
    - `data_source` - which data sources to load, default is `default_set`
    - `center` - whether to center the time variable to start at 0, default is true

    # Output
    - `data` - a tuple or dict of data sources
    """
    function load(args...; as="tuple", data_source=default_set, center=true,
            kws...)
        if :data_sources in keys(kws)
            data_source = kws[:data_sources]
        end
        if "beh" in data_source
            data_source = setdiff(data_source, ["beh"])
            push!(data_source, "behavior")
        end
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

        if !isempty(data_source) && isempty(load_order)
            throw(ValueError("data_source must contain at least one element"))
        end
        

        # DI each data source
        data = Dict{String,Any}()
        m = center ? nothing : 0 # nothing asks normalize to find the mintime and remove it, a number passed means use that for min
        for source ∈ load_order
            data[source] = DI.load_functions[source](args...)
            @assert typeof(data[source]) != Nothing
            if "time" ∈ names(data[source])
                @info "Centering $source and **converting to minutes**"
                data[source].time, m = normalize(data, data[source].time, data_source, m);
            end
        end
        _set_mintime!(m)

        # How do we return, dict or tuple?
        if as == "tuple"
            return (data[name] for name in data_source)
        else
            return data
        end
    end
    load(;kws...) = DI.load(default_args[1], default_args[2]; kws...)

    function get_load_kws(type::String; load_kws...)
        if type == "csv"
            kws = (csvkws..., load_kws...)
        elseif type == "arrow"
            kws = (arrowkws..., load_kws...)
        end
    end
    function load_table_at_path(path::String, type::String; load_kws...)
        if !(isfile(path))
            @error "File not found" path
        end
        if type == "csv"
            data = CSV.read(path, DataFrame; load_kws...)
        elseif type == "arrow"
            data = fix_complex(copy(DataFrame(Arrow.Table(path))))
        end
    end

    function load_table(animal::String, day::Int, pos...; 
            tablepath=nothing, 
            append="",
            load_kws::Union{Nothing,NamedTuple}=(;),
            kws...)
        if tablepath === nothing
            throw(ArgumentError("Must provide tablepath symbol... see path_functions dict in this module"))
        end

        tablepath = String(tablepath)
        path = path_functions[tablepath](animal, day, pos...; kws...)

        type = if :type in keys(kws)
            kws[:type]
        else
            load_default
        end

        function gettable(path)
            data = load_table_at_path(path, type; load_kws...)
            detect_complex_format_wrong = eltype.(eachcol(data)) .<: NamedTuple
            for i in findall(detect_complex_format_wrong)
                data[!,i] = getproperty.(data[!,i], :re) + 
                            getproperty.(data[!,i], :im)im
            end
            data
        end

        if append != ""
            path = split(path, ".")
            path[1] = path[1] * "$append"
            path = join(path, ".")
        end
        data = if '*' in path
            data=[gettable(p) for p in glob(path)]
            vcat(data..., cols=:union)
        else
            gettable(path)
        end

        data
    end

    function save_table_at_path(data::AbstractDataFrame, path::String, type::String;
            save_kws...)
        if type == "csv"
            data |> CSV.write(path)
        elseif type == "arrow"
            Arrow.write(path, data; save_kws...)
        elseif type == "arrow-append"
            path = replace(path, "arrow-append"=>"arrow")
            # Append does not work; Try Arrow.Writer()
            #Arrow.append(path, data; save_kws...)

            # Janky hacky fix
            existing_data = fix_complex(copy(DataFrame(Arrow.Table(path))))
            data = vcat(existing_data, data, cols=:union)
            Arrow.write(path, data; save_kws...)

        end
    end

    """
        save_table(data, pos...; tablepath=nothing, append="", kws...)

    Saves a table to a path
    # Input
    - `data` - a DataFrame
    - `pos...` - positional arguments to pass to the path function
                see `path_functions` dict in this module
    - `tablepath` - the tablepath to use, see `path_functions` dict in this module
    - `append` - whether to append to the table, default is false
    - `kws...` - keyword arguments to pass to the path function
    """
    function save_table(data::AbstractDataFrame, pos...; tablepath=nothing,
            append::String="", kws...)
        if tablepath === nothing
            throw(ArgumentError("Must provide tablepath symbol... see path_functions dict in this module"))
        end
        tablepath = String(tablepath)
        path = path_functions[tablepath](pos...; kws...)
        println("Saving $(tablepath) data at $(replace(path,"arrow-append"=>"arrow"))")
        if append != ""
            path = split(path, ".")
            path[1] = path[1] * "$append"
            path = join(path, ".")
        end
        if :type in keys(kws)
            type = kws[:type]
        elseif occursin(".", path)
            type = String(split(path,".")[2])
        else
            type = "csv"
        end
        save_table_at_path(data, path, type)
    end

    """
        tables_to_type(animal::String, day::Int; 
            inclusion=keys(save_functions),
            exclusion=[:lfp], from::String="csv", to::String="arrow",
            delete_prev::Bool=false)

    Conversion program to turn one set of filetypes into another

    # Params

    - `animal` :: the animal's name
    - `day` :: the day of the data files
    - `inclusion` :: Datatypes to translate
    - `exclusion` :: Datatypes to not translate
    - `from` :: Type to convert from
    - `to` :: Type to convert to
    """
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
            loadfunc, savefunc, pathfunc = load_functions[key], 
                                           save_functions[key],
                                           path_functions[key]
            if !isfile(pathfunc(animal, day; type=from))
               @warn "source not exist" animal day key
               continue
            end
            data = loadfunc(animal, day; type=from)
            savefunc(data, animal, day; type=to)
            if delete_prev
                pathfunc = path_functions[key]
                path = pathfunc(animal, day; type=from)
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
            if nonmissingtype(typeof(col)) <: Arrow.Struct || 
                nonmissingtype(eltype(col)) <: NamedTuple
                    if nonmissingtype(eltype(df[!,i])) == eltype(df[!,i])
                        df[!,i] = fix_complex(df[!,i])
                    else
                        tmp = df[:,i]
                        inds = (!).(ismissing.(tmp))
                        df[!,i] = Vector{Union{Missing,Complex}}(missing, size(df,1))
                        df[inds,i] = fix_complex(tmp[inds])
                    end
                #catch
                #    #@infiltrate
                #    @warn "exception caught in fix_complex"
                #end
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
                    #@infiltrate
                    @warn "exception caught in fix_complex"
                end
            end
        end
        return df
    end
    fix_rgba(x::AbstractVector) = fix_rgba.(x)
    fix_rgba(x::NamedTuple)     = RGBA(x.r, x.g, x.b, x.alpha)

end
