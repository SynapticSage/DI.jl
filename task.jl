
function taskpath(animal, day; type::String=load_default)
    csvFile = DrWatson.datadir("exp_raw",
                               "visualize_raw_neural",
                               "$(animal)_$(day)_task")
    csvFile = "$csvFile.$type"
end

function load_task(animal::String, day::Int; type::String=load_default, kws...)
    if type == "csv"
        typemap = Dict(Int64=>Int16);
        load_kws = (;strict=false, typemap=typemap, missingstring=["NaN", "NaNNaNi", "NaNNaNi,", ""], csvkws...)
    else
        load_kws = (;)
    end
    task = load_table(animal, day; tablepath=:task, type=type, load_kws=load_kws, kws...)
    if type == "csv"
        transform!(task, :x => pxtocm => :x, :y => pxtocm => :y)
    end
    return task
end

function save_task(behavior::AbstractDataFrame, pos...; kws...)
    save_table(behavior, pos...; tablepath=:task, kws...)
end

function well_locations(task; day=nothing, epoch=nothing)
end
