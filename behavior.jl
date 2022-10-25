module behavior

    import ..Load
    using DrWatson
    using DataFrames
    using Infiltrator
    export save_behavior, load_behavior, behaviorpath

    function behaviorpath(animal::String, day::Int, tag::String=""; type::String=Load.load_default)
        tag  = length(tag) == 0 ? tag : "_$tag"
        path = datadir("exp_raw", "visualize_raw_neural",
                         "$(animal)_$(day)_beh$tag.$type")
        if occursin("*",tag)
            path = glob(basename(path), dirname(path))
        end
        return path
    end

    function load_behavior(animal::String, day::Int, tag::String="";
        type::String=Load.load_default, kws...)
        function typeFunc(type, name)
            if occursin("Vec", string(name))
                type = ComplexF32;
            elseif name == "time" type = Float32;
            else
                type = nothing;
            end
            return type
        end
        if type == "csv"
            typemap = Dict(Int64=>Int16);

            load_kws = (;strict=false, missingstring=["NaNNaNi", "NaNNaNi,", ""], types=typeFunc, typemap=typemap, Load.csvkws...)
        else
            load_kws = (;)
        end
        beh = Load.load_table(animal, day, tag; tablepath=:behavior, type=type, 
                            load_kws=load_kws, kws...)
        if type == "csv"
            if beh.time isa Vector{String}
                beh.time = parse.(Float32, beh.time);
            end
            for col in names(beh)
                if (occursin("current", col) || occursin("egoVec", col)) &&
                    eltype(typeof(beh[!,col])) <: Real
                    try
                        replace!(beh[!,col], missing=>NaN)
                    catch; @infiltrate; end
                end
            end
            @assert ("x" ∈ names(beh)) "Fuck"
        end
        return beh
    end

    function save_behavior(behavior::AbstractDataFrame, pos...; kws...)
        Load.save_table(behavior, pos...; tablepath=:behavior, kws...)
    end

end
