module task

    export load_task, save_task, taskpath
    import ..Load
    using DrWatson
    using DataFrames

    function taskpath(animal, day; type::String=Load.load_default)
        csvFile = DrWatson.datadir("exp_raw",
                                   "visualize_raw_neural",
                                   "$(animal)_$(day)_task")
        csvFile = "$csvFile.$type"
    end

    function load_task(animal::String, day::Int; type::String=Load.load_default, kws...)
        if type == "csv"
            typemap = Dict(Int64=>Int16);
            load_kws = (;strict=false, typemap=typemap, missingstring=["NaN", "NaNNaNi", "NaNNaNi,", ""], Load.csvkws...)
        else
            load_kws = (;)
        end
        task = Load.load_table(animal, day; tablepath=:task, type=type, load_kws=load_kws, kws...)
        if type == "csv"
            transform!(task, :x => pxtocm => :x, :y => pxtocm => :y)
        end
        return task
    end

    function save_task(behavior::AbstractDataFrame, pos...; kws...)
        Load.save_table(behavior, pos...; tablepath=:task, kws...)
    end

end
