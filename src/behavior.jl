using DataFrames

function behaviorpath(animal::String, day::Int, tag::String=""; type::String=load_default)
    tag  = length(tag) == 0 ? tag : "_$tag"
    path = datadir("exp_raw", "visualize_raw_neural",
                     "$(animal)_$(day)_beh$tag.$type")
    if occursin("*",tag)
        path = glob(basename(path), dirname(path))
    end
    return path
end

function load_behavior(animal::String, day::Int, tag::String="";
    type::String=load_default, kws...)
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

        load_kws = (;strict=false, missingstring=["NaNNaNi", "NaNNaNi,", ""], types=typeFunc, typemap=typemap, csvkws...)
    else
        load_kws = (;)
    end
    beh = load_table(animal, day, tag; tablepath=:behavior, type=type, 
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
    save_table(behavior, pos...; tablepath=:behavior, kws...)
end

module behavior

    export annotate_pastFutureGoals
    export annotate_relative_xtime!
    using StatsBase
    using DataFrames
    using ProgressMeter

    function annotate_pastFutureGoals(beh::DataFrame; doPrevPast::Bool=false)

        beh = groupby(beh, [:epoch, :traj])
        for (g,group) in enumerate(beh)
            if g!=length(beh)
                group.futureStopWell .= beh[g+1].stopWell[1]
            end
            if g!=1
                group.pastStopWell .= beh[g-1].stopWell[1]
            end
        end
        beh = sort(combine(beh, identity), :time)
        replace!(beh.pastStopWell, missing=>-1)
        replace!(beh.futureStopWell, missing=>-1)
        if doPrevPast
            beh.hw .= argmax(StatsBase.fit(StatsBase.Histogram,
                                          filter(b->b!=-1,beh.stopWell), 1:6).weights)
            beh = groupby(beh, [:epoch, :block])
            for (g,group) in enumerate(beh)
                # TODO
            end
            replace!(beh.prevPastStopWell, missing=>-1)
            beh = sort(combine(beh, identity),:time)
        end
        beh
    end

    function annotate_poke(beh::DataFrame)
        pn = sort([name for name in names(beh) if occursin("poke_", name)])
        poke_matrix = replace(Matrix(beh[!, pn]), NaN=>0)
        poke_matrix = BitMatrix(poke_matrix)
        pn = replace([findfirst(row) for row in eachrow(poke_matrix)],nothing=>0)
        beh.poke = pn
        beh
    end

    function annotate_relative_xtime!(beh::DataFrame, x=:traj, on=:time)
        find_relative(t,m,M) = (t.-m)./(M-m)
        beh = groupby(beh, x)
        @showprogress 0.1 "adding rel $x on $on" for group in beh
            m, M = extrema(group[!,on])
            group[!, String(x)*"rel"*String(on)] = find_relative(group[!,on], m, M)
        end
        beh = combine(beh, identity)
    end
    
end
export behavior

