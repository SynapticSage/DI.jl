module celltet

    using DataFrames
    import ..Load
    using DrWatson
    export load_cells, load_tetrode, save_cells, save_cell_table, save_cell_taginfo
    using ProgressMeter
    using CSV
    import Table

    export cellpath
    export load_cells
    export save_cells
    export save_cell_taginfo

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

    function cellpaths(animal::String, day::Int, tag::String=""; type="csv", kws...)
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


    function load_cells(pos...; kws...)
        paths  = cellpaths(pos...; kws...)
        cells = DataFrame()
        @showprogress 0.1 "loading cell files" for path in paths
            cell = CSV.read(path, DataFrame; Load.csvkws...)
            cells = isempty(cells) ? cell : outerjoin(cells, cell, on=:unit, makeunique=true)
            Table.clean_duplicate_cols(cells)
        end
        return cells
    end

    function save_cells(cells::AbstractDataFrame, pos...; kws...)
        save_table(cells, pos...; tablepath=:cells, kws...)
    end
    """
    convenience wrapper to save_cells, ensuring you don't forget to tag the data
    if you meant to
    """
    function save_cell_taginfo(cells::AbstractDataFrame, animal::String, day::Int, tag::String; kws...)
        save_table(cells, animal, day, tag; tablepath=:cells, kws...)
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

    function fill_missing_cellinds(cells, missing_val=missing)
        cells = sort(cells, :unit)
        @error "Not implemented"
    end

    function load_tetrode(animal, day)
        cells = load_cells(animal,day)
        groups = groupby(cells,"tetrode")
        tetrodes = DataFrame()
        for group = groups
            n_cells = size(group,1)
            row = DataFrame(group[1,:])
            row[!, :n_cells] .= n_cells;
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

end
