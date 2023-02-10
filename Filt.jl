module Filt

    using DataFrames
    using DataStructures
    export get_filters, get_filters_precache, get_filter_req, 
             required_precache_functions, get_filter_prereq, precache
    using Base: merge
    import DIutils
    using Infiltrator

    function SPEED(x)
        abs.(x) .> 4
    end
    function SPEED_LIB(x)
        abs.(x) .> 2
    end
    function STILL(x)
        abs.(x) .> 2
    end
    speed     = OrderedDict("velVec"=>SPEED)
    speed_lib = OrderedDict("velVec"=>SPEED_LIB)
    still     = OrderedDict("velVec"=>STILL)

    function CORRECT(x)
        x.==1
    end
    function INCORRECT(x)
        x.==0
    end
    function NONTASK(x)
        (x.!=0) .&& (x.!=1)
    end
    function TASK(x)
        (x.==0) .|| (x.==1)
    end
    function HOME(X)
        X .== 'H'
    end
    function ARENA(X)
        X .== 'A'
    end
    correct   = OrderedDict("correct" => CORRECT)
    incorrect = OrderedDict("correct" => INCORRECT)
    nontask   = OrderedDict("correct" => NONTASK, "cuemem" => NONTASK)
    task      = OrderedDict("correct" => TASK)
    home      = OrderedDict("ha" => HOME)
    arena     = OrderedDict("ha" => ARENA)
    # Alias
    error     = incorrect

    function MEMORY(x)
        x.==1
    end
    function CUE(x)
        x.==0
    end
    cue       = OrderedDict("cuemem" => CUE)
    memory       = OrderedDict("cuemem" => MEMORY)

    notnan(x)       = OrderedDict(x  => x->((!).(isnan).(x)))
    minmax(x, m, M) = OrderedDict(x  => x-> x .>= m .&& x .<= M)
    max(x, M)       = OrderedDict(x  => x-> x .<= M)
    min(x, m)       = OrderedDict(x  => x-> x .>= m)

    """
    This is kind of more of a spike count filter than a cell count
    """
    req_spikes = 50
    cellcount       = OrderedDict([:time,:unit] =>
                          x -> groupby_summary_cond(x, :unit,
                                                  x->x.spikecount .> req_spikes, # > 50 spikes
                                                  nrow=>:spikecount))
    spikecount = cellcount

    req_traj = 4
    trajectory_diversity  = OrderedDict([:unit, :period] => 
            x -> groupby_summary_cond(x, :unit,
                                    x->x.trajdiversity .>= req_traj, # >= 3 trajectories
                                    :period =>
                                    (x->length(unique(x))) => :trajdiversity)
           )


    #function test_filt(spikes)
    #    println(all(combine(groupby(spikes[cellcount(spikes),:],:unit),nrow=>:count)[:,:count] .> 50))
    #    x = groupby_summary_filt(spikes, :unit, x->x.count.>50, nrow=>:count)
    #    println(all(combine(groupby(x,:unit),nrow=>:count)[:,:count] .> 50))
    #end

    """
    CROSS-VALIDATION BASED
    """
    function FIRSTHALF(x)
        1:size(x,1) .< round(size(x,1)/2)
    end
    function LASTHALF(x)
        1:size(x,1) .> round(size(x,1)/2)
    end
    function MOD(x, n::Int, v::Int)
        mod.(1:size(x,1), n) .== v
    end
    firsthalf      = OrderedDict("time" => FIRSTHALF)
    lasthalf       =  OrderedDict("time" => LASTHALF)
    evens          =  OrderedDict("time" => x->MOD(x, 2, 0))
    odds           =  OrderedDict("time" => x->MOD(x, 2, 1))

    function apply(filt::OrderedDict, x::DataFrame)
        x = copy(x)
        for (k,v) in filt
            x = x[v(x[!,k]), :]
        end
        x
    end

    """
    Currently matches N filters with matching keys

    ... this could do a lot more, like function as a swapin for the
    actual merge method for Dicts, and search for matching keys
    """
    function filtmerge(D::AbstractDict...)
        newD = Dict{keytype(D[1])}{Any}()
        K = Tuple(keys(D[1]))[1]
        newD[K] = [d[K] for d in D]
        newD
    end

    function groupby_summary_filt(df, splitby, summary_condition, combine_args...)
        groups = groupby(df, splitby, sort=true)
        summaries = combine(groups, combine_args...)
        summaries[!,"condition"] = summary_condition(summaries)
        groups = groups[ summaries.condition ]
        return combine(groups, identity)
    end

    """
             groupby_summary_cond(df::DataFrame, splitby, summary_condition,
                                  combine_args...)

    splits(splitby) a df and filters by a summary condition.

    ### Algorithm
    after it splits, we transform the splits with combine_args and apply a
    condition to the resultant combined table to filter the original dataframe.

    The original and combined dataframe are kept in register by an index
    variable.
    """
    function groupby_summary_cond(df::DataFrame, splitby, summary_condition,
                                  combine_args...)::BitVector
        columns = names(df)
        if splitby isa Vector{Symbol} || splitby isa Symbol
            columns = Symbol.(columns)
        end
        if all(in.(splitby, [columns]))
            df.condition = falses(size(df,1))
            df[!,:index] = 1:size(df,1)
            groups = groupby(df, splitby, sort=true)
            summaries = combine(groups, combine_args...)
            summaries[!,"condition"] = summary_condition(summaries)
            summaries = groupby(summaries, splitby)
            for (summary,group) in zip(summaries,groups)
                @assert summary[1,splitby] == group[1,splitby]
                if summary.condition[1]
                    group[!,:condition] .= true
                else
                    group[!,:condition] .= false
                end
            end
            return sort(combine(groups, identity), :index)[!,:condition]
        else
            return trues(size(df,1))
        end
    end
    function groupby_summary_condition_column!(df::DataFrame, splitby,
            summary_condition::Function, combine_args...; name=:condition,
            store_summary_stats::Bool=false)::DataFrame
        columns = names(df)
        if splitby isa Vector{Symbol} || splitby isa Symbol
            columns = Symbol.(columns)
        end
        if all(in.(splitby, [columns]))
            df[!, name] = falses(size(df,1))
            df[!,:index] = 1:size(df,1)
            groups = groupby(df, splitby, sort=true)
            summaries = combine(groups, combine_args...)
            summaries[!,name] = summary_condition(summaries)
            summaries = groupby(summaries, splitby)
            splitby_str = typeof(splitby) <: Vector ?
                    String.(splitby) : [String(splitby)]
            for (summary,group) in zip(summaries,groups)
                @assert summary[1,splitby] == group[1,splitby]
                if store_summary_stats
                    for col in setdiff(names(summary), splitby_str)
                        val = summary[!,col][1]
                        if !(ismissing(val))
                            group[!,col] .= val
                        end
                    end
                end
                if summary[!, name][1]
                    group[!,name] .= true
                else
                    group[!,name] .= false
                end
            end
            sort!(df, :index)
            select!(df, Not(:index))
        else
            df[!, name] =  trues(size(df,1))
            df
        end
    end

    """
        filtergroup

    filters groups/bundles of splits by a functoin computed per group. each
    pair is a column selector => function that evaluates of the whole group.
    """
    function filtergroup!(df::DataFrame, split, pairs::Pair...)
        function filtrations(x)
            answer = trues(size(df,1), 1)
            for (col, lambda) in pairs
                answer .&= lambda(df[!, col])
            end
            answer
        end
        cols = [col for (col, lambda) in pairs]

        Filt.groupby_summary_condition_column!(df, split, filtrations, cols...)
        deleteat!(df, df.condition .!= true)
    end
    filtercells(df::DataFrame, split, pairs::Pair...) = filtercells!(copy(df), split, pairs...)
    cellfilter!(df::DataFrame, pairs::Pair...) = filtergroup!(df, :unit, pairs...)
    cellfilter(df::DataFrame, pairs::Pair...) = filtergroup(df, :unit, pairs...)


    # Basic conditions
    conditions = (:all, :task, :correct, :error, :nontask, :memory, :cue,
                  :firsthalf, :lasthalf, :evens, :odds, :arena, :home)
    # Combination conditions
    cuecorr_conditions = (:cue_correct, :cue_error, :mem_correct, :mem_error)
    all_conditions = (conditions..., cuecorr_conditions...)

    # Create a set of predefined filter combinations
    """
        get_filters(initial::Union{Vector, Tuple})::OrderedDict

    get a dict with possible filtration options for the data based
    on behavior. each points to a function that will fitler your
    data frame for the condition indicated by the key name.

    If a column that one of these filters looks for is not in the dataframe,
    the function filtreg.filter() can either ignore or require it, e.g.
    spikecount. Behavior dataframes do not have a spikecount, unless, they've
    been munged to count some type of spiking event in the camera frame. The
    filter function that uses these will default to skipping it.
    
    # Params
    initial
        these are initial conditions that setup the :all times filter.
        each condition (e.g. :task, :correct etc) are also filtered
        by these initial conditions.
        default = (
        speed_lib => liberal speed filter (>2cm/s),
        spikecount => spike count > 50
                  )
    """
    function get_filters(initial=(speed_lib, spikecount);
                        keyfilter=nothing, keyfilterstr=nothing)
        initial = merge(initial...)
        filters = OrderedDict{Symbol,Union{OrderedDict,Nothing}}()

        filters[:all]         = initial
        filters[:odds]        = merge(filters[:all], Filt.odds)
        for key in setdiff(conditions, [:all])
            filters[key] = merge(filters[:all], getproperty(Filt, key))
        end
        for key in cuecorr_conditions
            key1, key2 = Symbol.(split(string(key), "_"))
            key1 = key1 == :mem ? :memory : key1
            filters[key] = merge(filters[:all], 
                                 getproperty(Filt, key1),
                                 getproperty(Filt, key2)
                                )
        end

        # HOME
        for key in union(conditions, cuecorr_conditions)
            newkey = Symbol("home_" * String(key))
            filters[newkey] = merge(filters[key], filters[:arena])
        end
        # ARENA
        for key in union(conditions, cuecorr_conditions)
            newkey = Symbol("arena_" * String(key))
            filters[newkey] = merge(filters[key], filters[:home])
        end

        filters[:none]        = nothing

        if keyfilter !== nothing
            K = filter(keyfilter, keys(filters))
            OrderedDict(k=>filters[k] for k in K)
        elseif keyfilterstr !==nothing
            K = Symbol.(filter(keyfilterstr, String.(collect(keys(filters)))))
            OrderedDict(k=>filters[k] for k in K)
        else
            filters
        end
    end

    home_conditions = [key for key in keys(get_filters())
                       if occursin("home_",string(key))]
    arena_conditions = [key for key in keys(get_filters())
                       if occursin("home_",string(key))]



                                                              
    #            ,---.                         |    o          
    #            |---',---.,---.,---.,---.,---.|---..,---.,---.
    #            |    |    |---'|    ,---||    |   |||   ||   |
    #            `    `    `---'`---'`---^`---'`   '``   '`---|
    #                                                     `---'

    """
        get_filters_precache

    precached filters

    See `get_filters`
    """
    get_filters_precache(;kws...) = get_filters(
                                         (speed_lib, spikecount,
                                          trajdiversitycached);
                                         kws...
                                        )

    function get_filters_desc()::OrderedDict
        filt_desc = OrderedDict(:all => "2cm/s")
        filters = collect(keys(get_filters()))
        for filt in filters
            #@info "get_filters_desc" filt
            filt_desc[filt] = join(("$(string(filt))", filt_desc[:all]), " :: ")
        end
        filt_desc
    end


    """
        get_filter_req

    Gets the fields required by the set of filters to operate
    """
    function get_filter_req(filts::AbstractDict)
        sets = [String.(f) for f ∈ keys(filts)]
        sets = [s isa Vector ? s : [s] for s in sets]
        unique(collect(Iterators.flatten(sets)))
    end

    function get_filter_prereq(filts::AbstractDict)
        [item for item in get_filter_req(filts)
         if Symbol(item) ∈ keys(precache_funcs)]
    end

    function required_precache_functions(filts::AbstractDict)
        [precache_funcs[Symbol(item)] for item in get_filter_prereq(filts)]
    end

    """

    returns the fields that are *required* to execute precaching faithfully
    """
    function required_precache_fields(filts::AbstractDict)
        Iterators.flatten([precache_field_reqs[Symbol(item)] for item in get_filter_prereq(filts)
                           if Symbol(item) ∈ keys(precache_field_reqs)])
    end

    function filters_use_precache(filts::AbstractDict)::Bool
        any(x -> Symbol(x) ∈ keys(precache_funcs), keys(filts))
    end

    """
    
    are precache inputs missing?
    """
    function missing_precache_input_cols(spikes::DataFrame,
                                  filts::AbstractDict)::Bool
        any(key -> key ∉ propertynames(spikes), 
            required_precache_fields(filts))
    end

    """
    
    are precache outputs missing?
    """
    function missing_precache_output_cols(spikes::DataFrame,
                                  filts::AbstractDict)::Bool
        any(key -> key ∉ propertynames(spikes), 
            intersect(keys(filts), keys(precache_funcs)))
    end

    """
        precache

    Precaches any filters that support precaching
    """
    function precache(spikes::DataFrame, beh::DataFrame, filts::AbstractDict; 
            kws...)
        reqfields = String.(required_precache_fields(filts))
        _, spikes = DIutils.filtreg.register(beh,spikes; 
                             transfer=reqfields)
        precache(spikes, filts; kws...)
    end
    function precache(spikes::DataFrame, filts::AbstractDict; kws...)::DataFrame
        funcs = required_precache_functions(filts)
        sp = spikes
        @assert spikes === sp
        for func in funcs
            spikes = func(spikes; kws...)
            @assert spikes === sp
        end
        @assert spikes === sp
        spikes
    end
    precache! = precache

    # Cache version
    function cachespikecount(spikes::DataFrame; kws...)
          groupby_summary_condition_column!(spikes, :unit,
                                  x -> x.spikecount .> req_spikes, # > 50 spikes
                                  nrow => :spikecount; name=:spikecount_cache, 
                                  kws...)
    end
    spikecountcached = OrderedDict(:spikecount_cache => x -> x) # stored true or false

    # Cache version
    function cachetrajdiversity(spikes::DataFrame; kws...)
        groupby_summary_condition_column!(spikes, :unit,
                                x -> x.trajdiversity .>= req_traj, 
                                :period => (x->length(unique(x))) => :trajdiversity
                                ; name=:trajdiversity_cache, kws...)
    end
    trajdiversitycached       = OrderedDict(:trajdiversity_cache =>
                                         x -> x) # stored true or false

    """
    stores functions used to precache each field
    """
    precache_funcs = Dict(:trajdiversity_cache => cachetrajdiversity, 
                         :spikecount_cache     => cachespikecount)
    """
        precache_field_reqs

    stores the fields required to precache
    """
    precache_field_reqs = Dict(:trajdiversity_cache => [:period])

    # Comparisons
    comparisons = (
       (:nontask, :task), (:cue, :memory), (:mem_error, :mem_correct),
       (:cue_error, :cue_correct), (:error, :correct), (:home, :arena)
      )
    """
    Obtain a list of comparison
    """
    function get_comparisons(possible_filts::Union{Nothing,Vector{Symbol},Vector{String}}=nothing; strfilter=nothing)
        if possible_filts === nothing
            possible_filts = keys(get_filters())
        else
            if eltype(possible_filts) == Symbol
                possible_filts = String.(possible_filts)
            end
        end
        C = collect(comparisons)
        custom_append!(c1, c2) = if c1 ∈ possible_filts && c2 ∈ possible_filts
            push!(C, (Symbol(c1), Symbol(c2)))
        else
            @info "skipping $c1 $c2" c1∈possible_filts c2∈possible_filts
        end
        for (c1,c2) in comparisons
            c1, c2 = String.((c1, c2))
            if occursin(c1, "home") || occursin(c2, "arena")
                continue
            end
            if c1 != "nontask"
                custom_append!("arena_"*c1, "arena_"*c2)
                custom_append!("home_"*c1, "home_"*c2)
                custom_append!("arena_"*c1, "home_"*c1)
            end
            if c2 != "nontask"
                custom_append!("arena_"*c2, "home_"*c2)
            end
        end
        C
    end

end
