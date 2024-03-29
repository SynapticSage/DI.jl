using DataFrames
using DrWatson
using ProgressMeter
using CSV
import DIutils: Table
using Glob
using Infiltrator

export cellpath
export load_cells
export save_cells
export save_cell_taginfo
export load_cells, load_tetrode, save_cells, save_cell_table, save_cell_taginfo
export cell_resort

"""
how to construct the path for a single cell/unit table
"""
function cellpath(animal::String, day::Int, tag::String=""; type="csv", kws...)
    if tag != "" && tag != "*"
        if !(startswith(tag,"_"))
            tag = "_$tag"
        end
    end
    csvFile = DrWatson.datadir("exp_raw", "visualize_raw_neural",
                               "$(animal)_$(day)_cell$tag.$type")
end

"""
how to construct the path for a single cell/unit table

# Arguments
- `animal`: animal name
- `day`: day number
- `tag`: tag for the cell/unit table; default is empty string; if tag is "*",
then the path is constructed by globbing the tag
# Returns
- `paths`: a vector of paths
"""
function cellpaths(animal::String, day::Int, tag::String=""; kws...)
    path = cellpath(animal, day, tag; kws...)
    if occursin("*", path)
        base, dir = basename(path), dirname(path)
        @debug "base=$base, dir=$dir"
        paths = glob(base, dir)
    else
        paths = [path]
    end
    return paths
end
"""
    cellpath(animal::String, dat::Int, tag::String; type="csv", kws...)

how to construct the path for a single cell/unit table

# Arguments
- `animal`: animal name
- `day`: day number
- `tag`: tag for the cell/unit table; default is empty string; if tag is "*",
then the path is constructed by globbing the tag
# Returns
- `paths`: a vector of paths
"""
function cellpaths(animal::String, dat::Int, tag::Vector{String}; kws...)
    paths = Vector{String}(undef, length(tag))
    for (i,t) in enumerate(tag)
        paths[i] = cellpath(animal, dat, t; kws...)
    end
    return Iterators.flatten(paths) |> collect
end


"""
    load_cells(pos...; type="arrow", kws...)

load cell/unit table(s) from path(s)

# Arguments
- `pos`: positional arguments for `cellpath`
- `type`: type of the table; default is "arrow"
- `kws`: keyword arguments for `cellpath`; which includes tag field 
         see `cellpath` for more details
"""
function load_cells(pos...; type="arrow", kws...)
    paths  = cellpaths(pos...; type, kws...)
    cells = DataFrame()
    @showprogress 0.1 "loading cell files" for path in paths
        cell = DI.load_table_at_path(path, type)
        cells = isempty(cells) ? cell : outerjoin(cells, cell, on=:unit, makeunique=true)
        Table.clean.clean_duplicate_cols(cells)
    end
    annotate_interneuron!(cells)
    return cells
end

function save_cells(cells::AbstractDataFrame, pos...; kws...)
    DI.save_table(cells, pos...; tablepath=:cells, kws...)
end

"""
    annotate_interneuron

for now my simplistic method is just to thresh around 5-7 hz
"""
function annotate_interneuron!(cells::DataFrame; thresh=6)
    cells[!,:interneuron] = cells.meanrate .> 6
    cells[!,:celltype] = Vector{Symbol}(undef, size(cells,1))
    cells.celltype[findall(cells.interneuron)] .= :int
    cells.celltype[findall((!).(cells.interneuron))] .= :pyr
    @info "cell types"  interneurons=sum(cells.interneuron) pyr=size(cells,1)-sum(cells.interneuron)
    nothing
end


"""
    save_cell_taginfo(cells::AbstractDataFrame, animal::String, day::Int, 
        tag::String; kws...)

convenience wrapper to save_cells, ensuring you don't forget to tag the data
if you meant to

in general, if you have some bonus cell information, it's a good idea to save
it using a tag. optionally, above, a cell table can be loaded from multiple
tags.

# Arguments
- `cells`: cell table
- `animal`: animal name
- `day`: day number
- `tag`: tag for the cell/unit table; default is empty string; if tag is "*",
then the path is constructed by globbing the tag
# Returns
- `paths`: a vector of paths
"""
function save_cell_taginfo(cells::AbstractDataFrame, animal::String, day::Int, 
    tag::String; kws...)
    if :type ∉ propertynames(kws)
        kws = (kws..., :type=>"arrow")
    end
    DI.save_table(cells, animal, day, tag; tablepath=:cells, kws...)
end

function cells_to_type(animal::String, day::Int, tag::String="*", 
        from::String="csv", to::String="arrow")
    paths = cellpaths(animal, day, tag)
    for path in paths
        data = load_table_at_path(path, from; load_kws...)
        path = replace(path, "."*from=>"."*to)
        savekws=(;)
        save_table_at_path(data, path, type; save_kws...)
    end
end

function fill_missing_cellinds(cells::DataFrame, missing_val=missing)
    cells = sort(cells, :unit)
    @error "Not implemented"
end

function load_tetrode(animal, day)
    cells = load_cells(animal,day)[!, Not([:csi, :meanrate, :propbursts, :tag])]
    groups = groupby(cells,"tetrode")
    tetrodes = DataFrame()
    for group = groups
        n_cells = size(group,1)
        numspikes = sum(group.numspikes)
        row = DataFrame(group[1,:])
        row[!, :n_cells] .= n_cells;
        row[!, :numspikes] .= numspikes;
        append!(tetrodes, row);
    end
    out = if "cell" in names(tetrodes)
            tetrodes[!, Not(:cell)]
        else
            tetrodes
        end
    return out
end

function cell_resort(cells::DataFrame, pos...; kws...)
    cells = sort(cells, pos...; kws...)
    if :origunit ∉ propertynames(cells)
        cells.origunit = cells.unit
    end
    cells.tmpunit = cells.unit
    cells.unit = 1:size(cells,1)
    cells
end
function cell_resort(cells::DataFrame, spikes::DataFrame, pos...; kws...)
    cells = cell_resort(cells, pos...; kws...)
    trades = Dict(x=>y for (x,y) in zip(cells.tmpunit, cells.unit))
    spikes.unit = replace(spikes.unit, trades...)
    cells, spikes
end


"""
    annotate_pyrlayer!

    pyramidal layer based on area == "CA1" and number of cells >= 4
    with meanrate < 6
"""
function annotate_pyrlayer!(cells)
    cells[!,:pyrlayer] = fill(false, size(cells,1))
    groups = :animal in propertynames(cells) ? [:animal,:tetrode] :
                                                :tetrode
    tetrodes = groupby(cells, groups)
    for tetrode in tetrodes
        if size(subset(tetrode, :meanrate => m->m .< 6),1) >= 4 && 
                                 tetrode[1,:area] == "CA1"
            tetrode[!,:pyrlayer] .= true
        else
            tetrode[!,:pyrlayer] .= false
        end
    end
end
