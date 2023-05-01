import ..DI
using Infiltrator
using DataFrames
import DataStructures: OrderedDict

export make_superanimal
"""
make_superanimal

#params

"""
function make_superanimal(animal_day_pairs::Tuple{String, Int}...;
        data_sources=DI.total_set, tag=nothing, numsuperanim::Int=0)
    for source in data_sources
        #data=[]
        loadmethod = DI.load_functions[source]
        savemethod = DI.save_functions[source]
        for (i,(animal, day)) in enumerate(animal_day_pairs)
            datum = loadmethod(animal, day; type="arrow")
            datum[!,:animal] .= animal
            datum[!,:day]    .= day
            type = i == 1 ? "arrow" : "arrow-append"
            savemethod(datum, "super", numsuperanim; type)
        end
    end
    GC.gc()
end

export make_superanimal_lfp
"""
make_superanimal

#params

"""
function make_superanimal_lfp(animal_day_pairs::Tuple{String, Int}...;
        tag=nothing, numsuperanim::Int=0)
    for tet in [:ca1ref, :pfcref, :default]
        #data=[]
        loadmethod = DI.load_lfp
        savemethod = DI.save_lfp
        data = []
        for (i,(animal, day)) in enumerate(animal_day_pairs)
            datum = loadmethod(animal, day; tet)
            datum[!,:animal] .= animal
            datum[!,:day]    .= day
            push!(data, datum)
        end
        data = vcat(data...)
        @assert(!isempty(data), "No data for $tet")
        DI.save_lfp(data, "super", numsuperanim; tet)
    end
    GC.gc()
end

"""
center_times

centers the times in the stored data per animal/day

Arguments
---------
- numsuperanim::Int=0
- data_sources=DI.total_set
- append="_cleaned"
- stagger_units::Bool=true -- stagger the units, avoid unit collisions
- stagger_time::Bool=true -- stagger the times, avoid time collisions
"""
function superanimal_clean_times_and_neurons(numsuperanim::Int=0; 
        data_sources=DI.total_set, append="_clean", 
        stagger_units::Bool=true,
        stagger_time::Bool=true,
        skip_list=["lfp"], 
        beh=nothing, cells=nothing,
        if_has_old_unit=:error)

    println("PROCESSING SETS: ", data_sources)
    println("SKIP LIST: ", skip_list)

    # Load up beh and cells
    if beh === nothing || cells === nothing
        beh, cells = DI.load_behavior("super",numsuperanim) ,
                     DI.load_cells("super",   numsuperanim)
    end
    # Figure out the time min and max per animal/day
    time_stats = combine(groupby(beh, [:animal,:day], sort=false), 
                       :time=>minimum, :time=>maximum)
    # shift the time end points down by one
    time_stats[:,:time_maximum_prev] .= [0; time_stats[1:end-1,:time_maximum]]
    # take the cumulative sum over the end points to figure out how
    # to stagger time across days/animals
    time_stats  = transform(time_stats, 
                          :time_maximum_prev=>cumsum=>:time_prev_maximum)
    time_stats.correct_factor = fill(Float64(0.0), size(time_stats,1))
    # return it to the grouped format for below
    time_stats_orig = time_stats

    # Figure out the unit min and max per animal/day
    unit_stats = combine(groupby(cells, [:animal,:day]), 
                       :unit=>minimum, :unit=>maximum)
    # shift the unit end points down by one
    unit_stats[:,:unit_maximum] .= [0; unit_stats[1:end-1,:unit_maximum]]
    # take the cumulative sum over the end points to figure out how
    # to stagger time across days/animals
    unit_stats  = transform(unit_stats, 
                          :unit_maximum=>cumsum=>:unit_prev_maximum)
    unit_stats.unit_minimum .= 0
    unit_stats_orig = unit_stats

    time_extrema_pre  = OrderedDict()
    time_extrema_post = OrderedDict()

    # Load each datasource and correct
    for source in data_sources
        if source in skip_list
            @info "skipping lfp"
            continue
        end
        loadmethod = DI.load_functions[source]
        savemethod = DI.save_functions[source]
        time_vars  = DI.time_vars[source]
        time_vars  = time_vars isa Vector ? time_vars : [time_vars]
        datum      = loadmethod("super", numsuperanim)
        time_stats = TS = copy(time_stats_orig)
        unit_stats = US = copy(unit_stats_orig)
        time_stats  = groupby(time_stats, [:animal, :day])
        unit_stats  = groupby(unit_stats, [:animal,:day])
        TS_vars = [:time_minimum, :time_maximum, :time_prev_maximum,
                   :correct_factor]
        println("datum: ", source, 
                "\nkey_order: ", keys(time_stats)|>collect.|>NamedTuple)
        if source == "cells"
            if :old_unit in propertynames(datum)
                if has_old_unit == :error
                    @error ":old_unit already in datum -- please check"
                elseif has_old_unit  == :continue
                    @warn ":old_unit already in datum"
                    continue
                else
                    throw(ValueError("Not a valid option"))
                end
            end
            datum[!,:old_unit] = datum[:,:unit]
        end
        groups = groupby(datum, [:animal,:day])
        key = ((keys(time_stats)) |> collect)[2]
        # PRINT pre-stagger states
        println(">>>>>>>>>>>>>>>>>>>>>>>")
        if :unit in propertynames(datum)
            println("Pre-stagger unit stats")
            println(combine(groupby(datum, [:animal,:day]), 
            :unit=>minimum, :unit=>maximum))
        end
        if :time in propertynames(datum)
            println("Pre-stagger time stats")
            println(combine(groupby(datum, [:animal,:day]), 
            :time=>minimum, :time=>maximum))
        end
        println(">>>>>>>>>>>>>>>>>>>>>>>")
        prev_end = -Inf
        for (k,key) in enumerate(keys(time_stats)) # iterate over animal, day
            mt, nt, group = time_stats[key], 
                        unit_stats[NamedTuple(key)], 
                        groups[    NamedTuple(key)]
            # 
            for timefield in time_vars
                if stagger_time
                    if key ∉ keys(time_extrema_pre)
                        time_extrema_pre = push!(time_extrema_pre, 
                                                 key=>extrema(group[!,timefield]))
                    elseif time_extrema_pre[key] == extrema(group[!,timefield])
                        @warn(
                    "time_extrema_pre[$key] != extrema(group[!,$timefield])")   
                        time_extrema_pre[key] = extrema(group[!,timefield])
                    end
                    Anext, Bprev = mt.time_minimum[1],
                                   mt.time_prev_maximum[1]
                    C = Anext - Bprev
                    group[!,timefield] .-= C # add previous max time of prev dataset
                    TS[1:k, TS_vars]   .-= C
                    time_extrema_post = push!(time_extrema_pre, 
                                             key=>extrema(group[!,timefield]))
                    if time_extrema_post[key][1] < prev_end
                        @infiltrate
                    end
                    prev_end = time_extrema_post[key][2]
                else
                    group[!,timefield] .-= mt.time_minimum      # center by behavior 0
                    mt.correct_factor   -= mt.time_minimum
                end
            end
            if stagger_units && :unit in propertynames(group)
                println("\tunit found in ", source, " ", key)
                # println(nt)
                delta = (nt.unit_prev_maximum)
                println("Applying delta=", delta)
                if !((nt.unit_minimum[1], nt.unit_maximum[1]) == extrema(group.unit))
                    # "you may need to redo the base super dataset")
                end
                group.unit .+= delta
                println("min of group[2] =", minimum(groups[2].unit))
            end
            # if nt.unit_prev_maximum > 0
            #     @assert minimum(group.unit) == (nt.unit_prev_maximum.+1)[1]
            # end
        end
        # Print post-stagger states
        println("<<<<<<<<<<<<<<<<<<<<<<<<")
        if :unit in propertynames(datum)
            println("Post-stagger unit stats")
            println(combine(groupby(datum, [:animal,:day]), 
            :unit=>minimum, :unit=>maximum))
        end
        if :time in propertynames(datum)
            println("Post-stagger time stats")
            println(combine(groupby(datum, [:animal,:day]), 
            :time=>minimum, :time=>maximum))
        end
        println("<<<<<<<<<<<<<<<<<<<<<<<<")
        # Sort and save
        if stagger_time && :time in propertynames(datum)
            sort!(datum, [:time])
        else
            sort!(datum, [x for x in [:animal, :day, :time] 
                if x in propertynames(datum)])
        end
        # If cells, save both csv and arrow, otherwise default
        if source == "cells"
            @assert datum.old_unit != datum.unit
            savemethod(datum, "super$append", numsuperanim)
            savemethod(datum, "super$append", numsuperanim; 
                              type="arrow")
        else
            savemethod(datum, "super$append", numsuperanim)
        end
        # Save the time correction factor
        if :unit in propertynames(datum)
            @assert(all(datum.unit .< 500), "unit is too large")
        end
    end
    # Save the time correction factor
    # try
    #     DI.save_table(combine(time_stats,identity), "super$append", numsuperanim; 
    #                       name="time_correction_factor")
    # catch e
    #     @infiltrate
    # end
    # GC.gc()
    @info "be sure to run conversion to arrow if you change anything : tables_to_type()"
end

"""
    superanimal_timeconversion
if both an uncleaned and cleaned versions of superanimal data, this
can be used to find the time conversion factor between the two
"""
function superanimal_timeconversion(superanimal, day; append="_clean")
    beh2=DI.load_behavior(superanimal*append, day)
    andays=vcat((unique(eachrow(beh2[!,[:animal, :day]])).|>DataFrame)...)
    animals, days = andays.animal, andays.day
    beh1 = vcat([
(x=DI.load_behavior(animal, day)[!,[:time]];x.animal.=animal;x) for (animal,day) in
            zip(animals, days)]...; cols=:intersect)
    sort!(beh1, [:animal, :time])
    sort!(beh2, [:animal, :time])
    g1=groupby(beh1, :animal)
    g2=groupby(beh2, :animal)
    animals = intersect(keys(g1)|>collect.|>NamedTuple,
                        keys(g2)|>collect.|>NamedTuple)
    Dict(animal[1]=>g2[animal].time[1]-g1[animal].time[1] for animal in animals)
end

