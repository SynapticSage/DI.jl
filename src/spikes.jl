
function spikespath(animal::String, day::Int; type::String=load_default)
    rawSpikingCSV = DrWatson.datadir("exp_raw",
                                     "visualize_raw_neural",
                                     "$(animal)_$(day)_labeled_spiking.$type"
                                    )
end

function load_spikes(animal::String, day::Int; type::String=load_default, beh=Nothing, kws...)
    if type == "csv"
        typemap = Dict(Int64=>Int16);
        load_kws = (;strict=false, typemap=typemap, missingstring=["NaN", "NaNNaNi", "NaNNaNi,", ""], csvkws...)
    else
        load_kws = (;)
    end
    spikes = load_table(animal, day; tablepath=:spikes, type=type, 
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
    save_table(spikes, pos...; tablepath=:spikes, kws...)
end
