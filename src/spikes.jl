module spikes

    export load_spikes, save_spikes, spikespath, column_save_spikes,
    column_load_spikes
    import ..DI
    using DrWatson
    using DataFrames
    using DataFrames: ColumnIndex
    
    CItype = Union{ColumnIndex, Vector{<:ColumnIndex}}

    index_vars = ["day","epoch","time"]

    function spikespath(animal::String, day::Int; type::String=DI.load_default)
        rawSpikingCSV = DrWatson.datadir("exp_raw",
                                         "visualize_raw_neural",
                                         "$(animal)_$(day)_labeled_spiking.$type"
                                        )
    end

    function column_load_spikes(column::CItype, animal::String, day::Int;
            type::String=DI.load_default, 
            kws...)
        if type == "csv"
            typemap = Dict(Int64=>Int16);
            load_kws = (;strict=false, typemap=typemap, missingstring=["NaN", "NaNNaNi", "NaNNaNi,", ""], DI.csvkws...)
        else
            load_kws = (;)
        end
        if column isa Vector
            append = "_" * join(String.(column),"-") 
        else
            append = "_" * String(column)
        end
        spikes = DI.load_table(animal, day; tablepath=:spikes, type=type,
                            append,
                            load_kws=load_kws, kws...)

    end
    function load_spikes(animal::String, day::Int;
            type::String=DI.load_default, beh=nothing, additional_columns=[],
            kws...)
        if type == "csv"
            typemap = Dict(Int64=>Int16);
            load_kws = (;strict=false, typemap=typemap, missingstring=["NaN", "NaNNaNi", "NaNNaNi,", ""], DI.csvkws...)
        else
            load_kws = (;)
        end
        spikes = DI.load_table(animal, day; tablepath=:spikes, type=type, 
                            load_kws=load_kws, kws...)

        # And let's add some brain area specific numbering
        if type == "csv"
            groups = groupby(spikes, "area");
            for g = 1:length(groups)
                unit = unique(groups[g].unit)
                areaunit = 1:length(unit);
                groups[g].areaunit = map(x -> Dict(unit .=> areaunit)[x],
                                         groups[g].unit)
            end
            spikes = combine(groups, x->x)
        end
        return spikes
    end

    function save_spikes(spikes::AbstractDataFrame, pos...; kws...)
        DI.save_table(spikes, pos...; tablepath=:spikes, kws...)
    end
    
    """
    Check point a spiking property computed, that
    can reload later
    """
    function column_save_spikes(column::CItype, spikes::AbstractDataFrame,
    pos...; column_transform=nothing, kws...)
        if column_transform !== nothing
            spikes = transform(spikes, column_transform...)
        end
        if column isa Vector
            append = "_" * join(String.(column),"-") 
            column = union(index_vars, column)
        else
            append = "_" * String(column)
            column = union(index_vars, [column])
        end
        DI.save_table(spikes[!,column], pos...; 
                        tablepath=:spikes,
                        append, kws...)
    end

end
