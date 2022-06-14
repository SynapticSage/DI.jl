module ripples

    export load_ripples, save_ripples, ripplespath
    import ..Load: load_default, load_kws, csvkws, save_table, load_table

    function ripplespath(animal, day; type::String=load_default)
        rippleFile = DrWatson.datadir("exp_raw",
                                   "visualize_raw_neural",
                                   "$(animal)_$(day)_ripple.$type") 
    end

    function load_ripples(animal::String, day::Int; type::String=load_default, kws...)
        if type == "csv"
            typemap = Dict(Int64=>Int16);
            load_kws = (;strict=false, typemap=typemap, missingstring=["NaN", "NaNNaNi", "NaNNaNi,", ""], csvkws...)
        else
            load_kws = (;)
        end
        ripples = load_table(animal, day; tablepath=:ripples, type=type, load_kws=load_kws, kws...)
    end

    function save_ripples(ripples::AbstractDataFrame, pos...; kws...)
        save_table(ripples, pos...; tablepath=:ripples, kws...)
    end

end
