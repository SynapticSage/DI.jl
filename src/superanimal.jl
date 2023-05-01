module superanimal
    import ..DI
    using Infiltrator
    using DataFrames

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
            if_has_old_unit=:error)

        beh, cells = DI.load_behavior("super",numsuperanim) ,
                     DI.load_cells("super",   numsuperanim)
        time_stats = combine(groupby(beh, [:animal,:day], sort=false), 
                           :time=>minimum, :time=>maximum)

        time_stats[:,:time_maximum] .= [0; time_stats[1:end-1,:time_maximum]]
        time_stats  = transform(time_stats, 
                              :time_maximum=>cumsum=>:time_prev_maximum)
        time_stats.correct_factor = fill(Float64(0.0), size(time_stats,1))
        time_stats  = groupby(time_stats, [:animal, :day])

        unit_stats = combine(groupby(cells, [:animal,:day]), 
                           :unit=>minimum, :unit=>maximum)
        unit_stats[:,:unit_maximum] .= [0; unit_stats[1:end-1,:unit_maximum]]
        unit_stats  = transform(unit_stats, 
                              :unit_maximum=>cumsum=>:unit_prev_maximum)
        unit_stats.unit_minimum .= 0
        unit_stats  = groupby(unit_stats, [:animal,:day])

        for source in data_sources
            if source in skip_list
                @info "skipping lfp"
                continue
            end
            @info source
            loadmethod = DI.load_functions[source]
            savemethod = DI.save_functions[source]
            time_vars  = DI.time_vars[source]
            time_vars  = time_vars isa Vector ? time_vars : [time_vars]
            datum      = loadmethod("super", numsuperanim)
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
            if :unit in propertynames(datum)
                println(combine(groupby(datum, [:animal,:day]), 
                :unit=>minimum, :unit=>maximum))
            end
            for key in keys(time_stats) # iterate over animal, day
                mt, nt, group = time_stats[key], 
                            unit_stats[NamedTuple(key)], 
                            groups[NamedTuple(key)]
                # 
                for timefield in time_vars
                    group[!,timefield] .-= mt.time_minimum      # center by behavior 0
                    mt.correct_factor -= mt.time_minimum
                    if stagger_time
                        group[!,timefield] .+= mt.time_prev_maximum - mt.time_minimum # add previous max time of prev dataset
                        mt.correct_factor += mt.time_prev_maximum - mt.time_minimum
                    end
                end
                if stagger_units && :unit in propertynames(group)
                    @info "unit found" source key 
                    println(nt)
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
            if stagger_time && :time in propertynames(datum)
                sort!(datum, [:time])
            else
                sort!(datum, [x for x in [:animal, :day, :time] 
                    if x in propertynames(datum)])
            end
            @infiltrate
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
        try
            DI.save_table(combine(time_stats,identity), "super$append", numsuperanim; 
                              name="time_correction_factor")
        catch e
            @infiltrate
        end
        GC.gc()
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

end
