#=
    lfp
Module for loading up and saving lfp related datasets 
=#
using DataFrames, NetCDF, Infiltrator, ProgressMeter, DrWatson
export lfppath, load_lfp, save_lfp, load_cycles, save_cycles, cyclepath
# TODO: EVENTUALLY, transfer all of this metadata into JSON or YAML files
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
"""
    get_tetrode_set
function for selecting named sets of tetrodes per animal
"""
function get_tetrode_set(tetrode_set; cells::Union{DataFrame,Nothing}=nothing, 
                         animal="super_clean", day=0)
    if tetrode_set     == :all # all tetrodes
        cells = cells === nothing ? DI.load_cells(animal, day) : cells
        tets = Dict(animal=>unique(subset(cells, :animal.==animal).tetrode) 
                    for animal in unique(cells.animal))
    elseif tetrode_set == :best # tetrodes with best pop phase locking across animals
        tets = Dict(
            "RY16" => [56],
            "RY22" => [29]
        )
    elseif tetrode_set == :synced 
        tets = Dict(
            "RY16" => [1],
            "RY22" => [29]
        )
    elseif tetrode_set == :ca1ref
        tets = Dict(
            "RY16" => [:ca1ref],
            "RY22" => [:ca1ref]
        )
    elseif tetrode_set == :pyr # pyramidal layer tetrodes
        DI.annotate_pyrlayer!(cells)
        tets = Dict(animal=>subset(cells, :animal=>a->a.==animal,
            :pyrlayer=>p->p.!=true, view=true).tetrode |> unique for animal in
                unique(cells.animal))
    end
end
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# LFP

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

export lfppath
"""
    lfppath(animal::String, day::Int; tet=nothing, type::String="nc", 
                 ref=nothing)

Obtains path to lfp file

# Arguments
 ---------
- `animal` -- 
"""
function lfppath(animal::String, day::Int; tet=nothing, type::String="nc", 
                 ref=nothing, write=false, tag=nothing, kws...)

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
    file = if ref !== nothing && ref !== false
        pathstring("ref",tet)
    else
        pathstring("",tet)
    end
    file = tag === nothing ? file : replace(file, ".nc" => "_$tag.nc")
end

"""
new doc
"""
function load_lfp(pos...; tet=nothing, vars=nothing, 
    subtract_earlytime=false, kws...)
    if :tetrode in propertynames(kws)
        tet = kws[:tetrode]
    end
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
        lfp = Vector{DataFrame}(undef,length(tet))
        prog = ProgressMeter.Progress(length(tet), 
            desc="Loading LFP, tetrode list")
        for (i,t) in enumerate(tet)|>collect
            println(t);
            lfp[i] = load_lfp(pos...; tet=t, vars=vars)
            next!(prog)
        end 
        lfp = vcat(lfp...)
    else
        lfpPath = lfppath(pos...; tet=tet, kws...)
        println(kws)
        if :append in keys(kws)
            append = kws[:append]
            append_str = !startswith(append,"_") ? "_$append" : append
            lfpPath = replace(lfpPath, ".nc" => "$append_str.nc")
        elseif :tag in keys(kws)
            tag = kws[:tag]
            tag_str = !startswith(tag,"_") ? "_$tag" : tag
            lfpPath = replace(lfpPath, ".nc" => "$tag_str.nc")
        end
        @info lfpPath
        v = NetCDF.open(lfpPath)
        if "Var1" in keys(v.vars)
            v.vars["time"] = v.vars["Var1"]
            pop!(v.vars, "Var1")
        end
        keyset = keys(v.vars)
        if vars !== nothing
            keyset = intersect(keyset, String.(vars))
        end
        println("Loading $keyset")
        lfp = Dict()
        prog = ProgressMeter.Progress(length(keyset), desc="Loading tet=$tet")
        for var in keyset
            q = try
                t = eltype(v.vars[var])
                var => Array{t}(v.vars[var]) 
            catch
                continue
                # @info "issue, tryin another method" var
                # @time var => Array(v.vars[var][1:length(v.vars[var]))
                # @info "success" var
            end
            push!(lfp, q)
            ProgressMeter.next!(prog)
        end
        lfp = DataFrame(Dict(var => vec(lfp[var]) 
            for var in keys(lfp)))
    end
    if subtract_earlytime
        lfp[!,:time] .-= DI.min_time_records[end]
    end
    if :broadraw in propertynames(lfp)
        lfp.broadraw = convert(Vector{Float16}, lfp.broadraw)
    end
    return lfp
end

"""
    save_lfp(l::AbstractDataFrame, pos...; tet=nothing, kws...)
Saves lfp to netcdf file
# Arguments
- `l::AbstractDataFrame` : lfp data
- `pos...` : position arguments to `lfppath`
- `tet=nothing` : tetrode to save; if nothing, uses default tetrode
- `kws...` : keyword arguments to `lfppath`
             - `type::String="nc"` : file type
             - `ref=nothing` : whether to use reference file
             - tag=nothing : tag to add to file name
            see `lfppath` for more details
"""
function save_lfp(l::AbstractDataFrame, pos...; tet=nothing, 
    handlemissing=:error, kws...)
    if tet == :default
        animal = pos[1]
        tet = default_tetrodes[animal]
    # elseif tet == :ca1ref
    #     animal = pos[1]
    #     tet = ca1ref_tetrodes[animal]
    end
    if :broadraw in propertynames(l)
        l.broadraw = convert(Vector{Float32}, l.broadraw)
    end
    function getkeys(lfpPath::String)
        ncFile = NetCDF.open(lfpPath)
        K = keys(ncFile.vars)
        NetCDF.close(ncFile) # id of the ncfile handle itself, may not be needed in new version
        K
    end
    if length(unique(l.tetrode)) == 1 && :append ∉ keys(kws)
        tet = l.tetrode[1];
    end
    lfpPath = lfppath(pos...; kws..., tet, write=true)
    lfpPath = if :tag in keys(kws)
        replace(lfpPath, ".nc" => "_$(kws[:tag]).nc")
    elseif :append in keys(kws)
        replace(lfpPath, ".nc" => "_$(kws[:append]).nc")
    else
        lfpPath
    end
    if any(Missing .<: eltype.(eachcol(first(l,10))))
        if handlemissing == :error
            error("missing values in lfp, use handlemissing=:drop " *
                  "if you want to drop them")
        elseif handlemissing == :drop
            l = dropmissing(l)
        end
    end
    @info "saving" lfpPath
    if isfile(lfpPath)
        @info "removing old file"
        rm(lfpPath)
    end
    d=NcDim("sample", size(l,1))
    varlist = Vector{NcVar}([])

    file = lfppath(pos...;kws..., write=true)
    dir, base = dirname(file), basename(file)
    dir = islink(dir) ? readlink(dir) : dir
    file = joinpath(dir, base)
    k = (names(l) |> collect)[1]
    for k in names(l)
        try
            var = NetCDF.NcVar(k, d)
            #var.nctype=NetCDF.getNCType(eltype(original_nc[k]))
            var.nctype=NetCDF.getNCType(eltype(l[!,k]))
            push!(varlist,var)
        catch
            @error "Check type of $k"
        end
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

"""
    load_lfp_tetrodeset(animal, day, tetrode_set; convert=true, cells=nothing)

Load LFP for a tetrode set
# Arguments
- `animal` : animal name
- `day` : day
- `tetrode_set` : tetrode set
- `convert=true` : whether to 0 time to first behavioral event
- `cells=nothing` : cells to load
# Returns
- `lfp` : lfp dataframe
"""
function load_lfp_tetrodeset(animal, day, tetrode_set; convert::Bool=true, 
    cells::Union{AbstractDataFrame, Nothing}=nothing)
    lfp_convert  = DI.get_0time_pre_superanimal()
    all_tetrodes = DI.get_tetrode_set(tetrode_set; cells, animal, day)
    println("Loading LFP tetrodes = $tetrode_set") 
    printstyled("tet = ", all_tetrodes[animal], "\n", blink=true)
    lfp = DI.load_lfp(animal, day, tet=all_tetrodes[animal])
    if convert
        lfp.time .-= lfp_convert[animal];
    end
    lfp
end

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Oscillation cycles
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

function cyclepath(animal::String, day::Int, 
    tetrode::Union{String,Int,Nothing}=nothing,
    tag=nothing; type::String="arrow", kws...)
    if tag !== nothing
        tag = "_$(tag)"
    elseif :append in keys(kws)
        tag = "_$(kws[:append])"
    else
        tag = ""
    end
    if tetrode !== nothing
        tetstr = "_tet=$(tetrode)"
    else
        tetstr = ""
    end
    DrWatson.datadir("exp_raw", "visualize_raw_neural",
        "$(animal)_$(day)$(tetstr)_cycles$tag.$type")
end
function save_cycles(cycles, pos...; kws...)
    path=cyclepath(pos...; kws...)
    Arrow.write(path, cycles; compress=:lz4)
end
function load_cycles(pos...; type="arrow", kws...)
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
