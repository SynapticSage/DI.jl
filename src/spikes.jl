
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
# spikespath(animal, day, tag) = replace(spikespath(animal, day), 
#     r"\.(csv|arrow)$"=>"_$(tag).\\1")

function load_spikes(animal::String, day::Int, tag=nothing;
        type::String=DI.load_default, beh=nothing, additional_columns=[],
        kws...)
    if type == "csv"
        typemap = Dict(Int64=>Int16);
        load_kws = (;strict=false, typemap=typemap, 
            missingstring=["NaN", "NaNNaNi", "NaNNaNi,", ""], DI.csvkws...)
    else
        load_kws = (;)
    end
    if tag !== nothing
        kws=(;append= tag == "*" ? tag : "_" * tag)
    end
    spikes = DI.load_table(animal, day; tablepath=:spikes, type=type, 
                        load_kws=load_kws, kws...)

    # And let's add some brain area specific numbering
    if type == "csv"
        groups = groupby(spikes, "area");
        for g in eachindex(groups)
            unit = unique(groups[g].unit)
            areaunit = 1:length(unit);
            groups[g].areaunit = map(x -> Dict(unit .=> areaunit)[x],
                                     groups[g].unit)
        end
        spikes = combine(groups, x->x)
    end
    return spikes
end

function save_spikes(spikes::AbstractDataFrame, animal::String, day::Int, 
                     tag=nothing; kws...)
    if tag !== nothing
        kws=(;append="_" * tag)
    end
    DI.save_table(spikes, animal, day; tablepath=:spikes, kws...)
end

"""
    save_spikes_taginfo(spikes::AbstractDataFrame, pos...; tag=nothing, 
                        append=nothing, kws...)
# Arguments
- `spikes`: a dataframe of spikes
- `pos...`: the position of the spikes table
- `tag`: the tag for the spikes table
- `append`: the append for the spikes table
- `kws...`: other keyword arguments
# Returns
- `spikes`: a dataframe of spikes
"""
function save_spikes_taginfo(spikes::AbstractDataFrame, animal::String, 
    day::Int, tag=nothing; 
    append=nothing, kws...)
    if tag !== nothing
        append = "_" * tag
    elseif append !== nothing
        append = startswith(append,"_") ? "_" * append : append
    else
        @error "Either tag or append must be specified"
    end
    DI.save_table(spikes, animal, day; tablepath=:spikes, append=append, kws...)
end


# # Alternative ways to save spiking data with tags ... based on the columns
# # of the dataframe that one would like to save
#
# """
# Checkpoint a spiking property computed based on the column of interest
# """
# function column_save_spikes(column::CItype, spikes::AbstractDataFrame,
# pos...; column_transform=nothing, kws...)
#     if column_transform !== nothing
#         spikes = transform(spikes, column_transform...)
#     end
#     if column isa Vector
#         append = "_" * join(String.(column),"-") 
#         column = union(index_vars, column)
#     else
#         append = "_" * String(column)
#         column = union(index_vars, [column])
#     end
#     DI.save_table(spikes[!,column], pos...; 
#                     tablepath=:spikes,
#                     append, kws...)
# end
#
# """
# Load a spiking property computed based on the column of interest, if a special
# save version of that column exists
# """
# function column_load_spikes(column::CItype, animal::String, day::Int;
#         type::String=DI.load_default, 
#         kws...)
#     if type == "csv"
#         typemap = Dict(Int64=>Int16);
#         load_kws = (;strict=false, typemap=typemap, missingstring=["NaN", "NaNNaNi", "NaNNaNi,", ""], DI.csvkws...)
#     else
#         load_kws = (;)
#     end
#     if column isa Vector
#         append = "_" * join(String.(column),"-") 
#     else
#         append = "_" * String(column)
#     end
#     spikes = DI.load_table(animal, day; tablepath=:spikes, type=type,
#                         append,
#                         load_kws=load_kws, kws...)
# end
