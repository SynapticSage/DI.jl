#=
    lfp

Module for loading up and saving lfp related datasets 
=#

using DataFrames, NetCDF, Infiltrator, ProgressMeter, DrWatson
export lfppath, load_lfp, save_lfp, load_cycles, save_cycles, cyclepath

default_tetrodes = Dict(
    #"RY22" => 16, # lot of cells, theta = ass
    #"RY22" => 7, # good theta
    "RY22" => 18, # good theta
    "RY16" => 5,  # good theta 
    "super" => :default
    )

ca1ref_tetrodes = Dict(
    #"RY22" => 16, # lot of cells, theta = ass
    #"RY22" => 7, # good theta
    "RY22" => 20, # good theta
    "RY16" => 17,  # good theta 
    "super" => :ca1ref
    )

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# LFP
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

export lfppath
"""
    lfppath(animal::String, day::Int; tet=nothing, type::String="nc", 
                 ref=nothing)

Obtains path to lfp file
"""
function lfppath(animal::String, day::Int; tet=nothing, type::String="nc", 
                 ref=nothing, write=false, kws...)

    pathstring(ref,::Nothing) = DrWatson.datadir("exp_raw",
                                 "visualize_raw_neural", 
                                 "$(animal)_$(day)_rhythm$ref.$type")
    pathstring(ref,tet::T where T<:Union{Int,UInt8}) = DrWatson.datadir("exp_raw",
                          "visualize_raw_neural", 
                          "$(animal)_$(day)_rhythm$(ref)_$(Int(tet)).$type")
    pathstring(ref,tet::T where T<:Union{String,Symbol}) = DrWatson.datadir("exp_raw",
                          "visualize_raw_neural", 
                          "$(animal)_$(day)_rhythm$(ref)_$(string(tet)).$type")

    if ref === nothing
        ref = if isfile(pathstring("ref",tet))
            ref = true;
        elseif isfile(pathstring("",tet))
            ref = false;
        else
            if write == false
                error("No rhythm file $(pathstring("ref",tet))")
            else
                ref = false;
            end
        end
    end
    if ref !== nothing && ref !== false
        pathstring("ref",tet)
    else
        pathstring("",tet)
    end
end

function load_lfp(pos...; tet=nothing, vars=nothing, 
        subtract_earlytime=false, kws...)
    if tet == :default
        animal = pos[1]
        tet = default_tetrodes[animal]
    # elseif tet == :ca1ref
    #     animal = pos[1]
    #     tet = ca1ref_tetrodes[animal]
    elseif tet isa String
        tets  = DI.load_tetrode(pos...)
        tets  = groupby(tets,:area)[(;area=tet)]
        tet  = collect(tets.tetrode)
    end
    if tet isa Vector
        lfp = [(println(t); load_lfp(pos...; tet=t, vars=vars))
               for t in tet]
        lfp = vcat(lfp...)
    else
        lfpPath = lfppath(pos...; tet=tet, kws...)
        @info lfpPath
        v = NetCDF.open(lfpPath)
        if "Var1" in keys(v.vars)
            v.vars["time"] = v.vars["Var1"]
            pop!(v.vars, "Var1")
        end
        keyset = keys(v.vars)
        if vars !== nothing
            keyset = String.(vars)
        end
        lfp = Dict()
        @showprogress "lfp" for var in keyset
            q = try
                var => Array(v.vars[var]) 
            catch
                continue
                # @info "issue, tryin another method" var
                # @time var => Array(v.vars[var][1:length(v.vars[var]))
                # @info "success" var
            end
            push!(lfp, q)
        end
        lfp = DataFrame(Dict(var => vec(lfp[var]) 
            for var in keys(lfp)))
    end
    if subtract_earlytime
        lfp[!,:time] .-= DI.min_time_records[end]
    end
    return lfp
end

function save_lfp(l::AbstractDataFrame, pos...; tet=nothing, kws...)
    if tet == :default
        animal = pos[1]
        tet = default_tetrodes[animal]
    # elseif tet == :ca1ref
    #     animal = pos[1]
    #     tet = ca1ref_tetrodes[animal]
    end
    function getkeys(lfpPath::String)
        ncFile = NetCDF.open(lfpPath)
        K = keys(ncFile.vars)
        NetCDF.close(ncFile) # id of the ncfile handle itself, may not be needed in new version
        K
    end
    if length(unique(l.tetrode)) == 1
        tet = l.tetrode[1];
    end
    lfpPath = lfppath(pos...; kws..., tet, write=true)
    @info "saving" lfpPath
    if isfile(lfpPath)
        rm(lfpPath)
    end
    d=NcDim("sample", size(l,1))
    varlist = Vector{NcVar}([])

    file = lfppath(pos...;kws..., write=true)
    dir, base = dirname(file), basename(file)
    dir = islink(dir) ? readlink(dir) : dir
    file = joinpath(dir, base)
    for k in names(l)
        var = NetCDF.NcVar(k, d)
        #var.nctype=NetCDF.getNCType(eltype(original_nc[k]))
        var.nctype=NetCDF.getNCType(eltype(l[!,k]))
        push!(varlist,var)
    end
    NetCDF.create(lfpPath, varlist)
    ncFile = NetCDF.open(lfpPath; mode=NC_WRITE)
    for (i,(key,value)) in enumerate(zip(names(l),eachcol(l)))
        @debug "file=$lfpPath ∃, but key=$key ∉ keys"
        NetCDF.putvar(ncFile[key], Array(value))
    end
end

"""
    split_lfp_by_tet(pos...; lfp=nothing, vars=nothing, kws...)

split up a dataframe of tetrodes into different files
"""
function split_lfp_by_tet(pos...; lfp::Union{DataFrame,Nothing}=nothing, 
                          vars=nothing, kws...)
    if lfp === nothing
        lfp = load_lfp(pos...; vars=vars)
    end
    if vars === nothing
        vars = names(lfp)
    end
    #lfp = NetCDF.open(lfppath(pos...;kws...))
    lfp = groupby(lfp, :tetrode)
    @showprogress for l in lfp
        save_lfp(l, pos...; tet=l.tetrode[begin], vars=vars, kws...)
    end
end

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Oscillation cycles
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

function cyclepath(animal::String, day::Int, tetrode::Union{String,Int}; type::String="csv")
    DrWatson.datadir("exp_raw", "visualize_raw_neural",
                              "$(animal)_$(day)_tet=$(tetrode)_cycles.$type")
end
function save_cycles(cycles, pos...)
    cycles |> CSV.write(cyclepath(pos...))
end
function load_cycles(pos...; type="csv", kws...)
    DI.load_table(pos...; tablepath=:cycles, type=type, kws...)
end

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Coherence
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
export cohpath
function cohpath(animal::String, day::Int; type::String="arrow")
    DrWatson.datadir("exp_raw", "visualize_raw_neural",
                              "$(animal)_$(day)_coh.$type")
end
export load_coh
function load_coh(pos...; type="arrow", kws...)
    DI.load_table(pos...; tablepath=:coh, type, kws...)
end

export avgcohpath
function avgcohpath(animal::String, day::Int; type::String="arrow")
    DrWatson.datadir("exp_raw", "visualize_raw_neural",
                              "$(animal)_$(day)_avgcoh.$type")
end
export load_avgcoh
function load_avgcoh(pos...; type="arrow", kws...)
    DI.load_table(pos...; tablepath=:avgcoh, type, kws...)
end
export save_avgcoh
function save_avgcoh(pos...; type="arrow", kws...)
    DI.save_table(pos...; tablepath=:avgcoh, type, kws...)
end

