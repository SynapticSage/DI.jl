export load_ripples, save_ripples, ripplespath
using DataFrames
import ..DI
import DrWatson
import CSV, Arrow

function ripplespath(animal, day; type::String=DI.load_default)
    rippleFile = DrWatson.datadir("exp_raw",
                               "visualize_raw_neural",
                               "$(animal)_$(day)_ripple.$type") 
end

"""
    load_ripples(animal::String, day::Int; type::String=DI.load_default,
                      kws...)
"""
function load_ripples(animal::String, day::Int; type::String=DI.load_default,
                      kws...)
    if type == "csv"
        typemap = Dict(Int64=>Int16);
        load_kws = (;strict=false, typemap=typemap, missingstring=["NaN",
            "NaNNaNi", "NaNNaNi,", ""], DI.csvkws...)
    else
        load_kws = (;)
    end
    ripples = DI.load_table(animal, day; tablepath=:ripples, type=type,
        load_kws=load_kws, kws...)
    postprocess_ripples(ripples)
end
load_ripples(;kws...) = load_ripples(DI.default_args...;kws...)

function save_ripples(ripples::AbstractDataFrame, pos...; kws...)
    DI.save_table(ripples, pos...; tablepath=:ripples, kws...)
end


function postprocess_ripples(ripples::AbstractDataFrame)
    goodinds = (ripples.stop - ripples.start) .< 1
    ripples[goodinds, :]
end
