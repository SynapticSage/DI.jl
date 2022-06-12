using NetCDF

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
    cycles = load_table(pos...; tablepath=:cycles, type=type, 
                        load_kws=load_kws, kws...)
end

module lfp
    using DataFrames
    using Statistics
    using DirectionalStatistics
    using ImageFiltering
    include("../table.jl")
    function phase_to_radians(phase)
        phase = Float32.(phase)
        phase = 2*π*(phase .- minimum(phase))./diff([extrema(phase)...]) .- π
        #phase = convert.(Float16, phase)
    end
    function annotate_cycles(lfp; phase_col="phase", method="peak-to-peak")
        phase = lfp[!, phase_col]
        lfp.phase = phase_to_radians(lfp[:,"phase"])
        println("Method=$method")
        if method == "resets"
            Δₚ = [0; diff(phase)]
            reset_points = UInt32.(Δₚ .< (-1.5*π))
            cycle_labels = accumulate(+, reset_points)
        elseif method == "peak-to-peak"
            step_size = median(diff(phase))
            Δₚ = [0; diff(phase)]
            #falling_zero_point = [(phase[1:end-1] .>=0) .& (phase[2:end] .<0) ; false]
            rising_zero_point = [(phase[2:end] .>=0) .& (phase[1:end-1] .<0) ; false]
            cycle_labels = accumulate(+, rising_zero_point)
            lfp[!,"phase"] = mod2pi.(lfp[!,"phase"])
        elseif method == "trough-to-trough"
            step_size = median(diff(phase))
            Δₚ = [0; diff(phase)]
            falling_zero_point = [(phase[1:end-1] .>=0) .& (phase[2:end] .<0) ; false]
            #rising_zero_point = [(phase[2:end] .>=0) .& (phase[1:end-1] .<0) ; false]
            cycle_labels = accumulate(+, falling_zero_point)
            lfp[!,"phase′"] = mod.(lfp[!,"phase"] .- pi, 2*pi)
        else
            throw(ArgumentError("Unrecognized method=$method"))
        end
        lfp[!,"cycle"] = cycle_labels
        return lfp
    end
    function mean_lfp(lfp; mean_fields=["phase","amp","raw"], func=Circular.median)
        lfp = groupby(lfp, :time)
        non_mean_fields = setdiff(names(lfp), mean_fields)
        lfp =combine(lfp, mean_fields.=>func.=>mean_fields, 
                     non_mean_fields.=>first.=>non_mean_fields)
        return lfp[!, Not(:tetrode)]
    end
    function gauss_lfp(lfp; fields=["phase"], gaussian=3)
        kernel = Kernel.gaussian((gaussian,))
        for field in fields
            lfp[!,field] = imfilter(lfp[!,field], kernel)
        end
        return lfp
    end
    """
    weighted_lfp

    creates an amplitude weighted average of fields across tetrodes
    """
    function weighted_lfp(lfp; mean_fields=["phase","amp","raw"],
            weighting="amp")
        lfp = groupby(lfp, :time)
        non_mean_fields = setdiff(names(lfp), mean_fields)
        new = DataFrame()
        @time for lf in lfp
            item = DataFrame(lf[1,:])
            for field in mean_fields
                item[!, field] .= sum(lf[!, field] .* lf[!, weighting])/sum(lf.amp)
            end
            append!(new, item)
        end
        return new[!, Not(:tetrode)]
    end

    function unstack_tetrode(df; measure::Symbol=:phase)
        unstack(df[!, [:time, :tetrode, measure]], :tetrode, measure)
    end

    function get_cycle_table(lfp, pos...; kws...)
        @assert "cycle" in names(lfp)
        tab = table.get_periods(lfp, "cycle", :amp=>mean, pos...; kws...)
        return tab
    end
    getTet(L::DataFrame, T::Int) = filter(:tetrode=> t->t==T, L)

end
export lfp

