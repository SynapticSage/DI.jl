module lfp
    using DataFrames
    using NetCDF
    using Infiltrator
    using ProgressMeter
    import ..Load
    using DrWatson
    export lfppath, load_lfp, save_lfp, load_cycles, save_cycles, cyclepath

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # LFP
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

    function lfppath(animal::String, day::Int; tet=nothing, type::String="nc")
        if tet == nothing
            netcdf = DrWatson.datadir("exp_raw", "visualize_raw_neural",
                                             "$(animal)_$(day)_rhythm.$type")
        else
            netcdf = DrWatson.datadir("exp_raw", "visualize_raw_neural",
                                             "$(animal)_$(day)_rhythm_$tet.$type")
        end
    end

    function load_lfp(pos...; tet=nothing, vars=nothing, kws...)
        if tet isa Vector
            lfp = [load_lfp(pos...; tet=t, vars=vars)
                   for t in tet]
            lfp = vcat(lfp...)
        else
            lfpPath = lfppath(pos...; tet=tet)
            @info lfpPath
            v = NetCDF.open(lfpPath)
            if "Var1" in keys(v.vars)
                v.vars["time"] = v.vars["Var1"]
                pop!(v.vars, "Var1")
            end
            keyset = keys(v.vars)
            if vars != nothing
                keyset = String.(vars)
            end
            lfp = Dict(var => Array(v.vars[var]) 
                       for var in keyset)
            lfp = DataFrame(Dict(var => vec(lfp[var]) 
                                 for var in keyset))
        end
        return lfp
    end

    function save_lfp(l::AbstractDataFrame, pos...; tet=nothing, kws...)
        function getkeys(lfpPath::String)
            ncFile = NetCDF.open(lfpPath)
            K = keys(ncFile.vars)
            NetCDF.close(ncFile) # id of the ncfile handle itself, may not be needed in new version
            K
        end
        lfpPath = lfppath(pos...; tet=l.tetrode[1])
        @debug "path=$lfpPath"
        if isfile(lfpPath)
            rm(lfpPath)
        end
        d=NcDim("sample", size(l,1))
        @infiltrate
        varlist = Vector{NcVar}([])
        for k in names(l)
            var = NetCDF.NcVar(k, d)
            var.nctype=NetCDF.getNCType(eltype(original_nc[k]))
            push!(varlist,var)
        end
        NetCDF.create(lfpPath, varlist)
        ncFile = NetCDF.open(lfpPath; mode=NC_WRITE)
        for (i,(key,value)) in enumerate(zip(names(l),eachcol(l)))
            @debug "file=$lfpPath ∃, but key=$key ∉ keys"
            NetCDF.putvar(ncFile[key], Array(value))
        end
    end

    function split_lfp_by_tet(pos...; lfp=nothing, vars=nothing, kws...)
        if lfp == nothing
            lfp = load_lfp(pos...; vars=vars)
        end
        if vars == nothing
            vars = names(lfp)
        end
        original_nc = NetCDF.open(lfppath(pos...))
        lfp = groupby(lfp, :tetrode)
        @showprogress for l in lfp
            save_lfp(l, pos...; tet=l.tetrode[begin], vars=vars, kws...)
        end
    end

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # Oscillation cycles
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

    function cyclepath(animal::String, day::Int, tetrode::Union{String,Int}; ext::String="csv")
        csv = DrWatson.datadir("exp_raw", "visualize_raw_neural",
                                  "$(animal)_$(day)_tet=$(tetrode)_cycles.$ext")
    end
    function save_cycles(cycles, pos...)
        cycles |> CSV.write(cyclepath(pos...))
    end
    function load_cycles(pos...)
        cycles = Load.load_table(pos...; tablepath=:cycles, type=type, 
                            load_kws=load_kws, kws...)
    end


end
