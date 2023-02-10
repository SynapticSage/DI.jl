module ripples

    export load_ripples, save_ripples, ripplespath
    using DataFrames
    import ..Load
    import DrWatson
    import CSV, Arrow

    function ripplespath(animal, day; type::String=Load.load_default)
        rippleFile = DrWatson.datadir("exp_raw",
                                   "visualize_raw_neural",
                                   "$(animal)_$(day)_ripple.$type") 
    end

    function load_ripples(animal::String, day::Int; type::String=Load.load_default, kws...)
        if type == "csv"
            typemap = Dict(Int64=>Int16);
            load_kws = (;strict=false, typemap=typemap, missingstring=["NaN", "NaNNaNi", "NaNNaNi,", ""], Load.csvkws...)
        else
            load_kws = (;)
        end
        ripples = Load.load_table(animal, day; tablepath=:ripples, type=type, load_kws=load_kws, kws...)
    end
    load_ripples(;kws...) = load_ripples(Load.default_args...;kws...)

    function save_ripples(ripples::AbstractDataFrame, pos...; kws...)
        Load.save_table(ripples, pos...; tablepath=:ripples, kws...)
    end

end
