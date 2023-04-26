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
    function center_and_stagger_times_and_neurons(numsuperanim::Int=0; 
            data_sources=DI.total_set, append="_cleaned", 
            stagger_units::Bool=true,
            stagger_time::Bool=true,
            if_has_old_unit=:error)

        beh, cells = DI.load_behavior("super",numsuperanim) ,
                     DI.load_cells("super",   numsuperanim)
        time_stats = combine(groupby(beh, [:animal,:day]), 
                           :time=>minimum, :time=>maximum)

        time_stats[:,:time_maximum] .= [0; time_stats[1:end-1,:time_maximum]]
        time_stats  = transform(time_stats, 
                              :time_maximum=>cumsum=>:time_prev_maximum)
        time_stats  = groupby(time_stats, [:animal, :day])

        unit_stats = combine(groupby(cells, [:animal,:day]), 
                           :unit=>minimum, :unit=>maximum)
        unit_stats[:,:unit_maximum] .= [0; unit_stats[1:end-1,:unit_maximum]]
        unit_stats  = transform(unit_stats, 
                              :unit_maximum=>cumsum=>:unit_prev_maximum)
        unit_stats.unit_minimum .= 0
        unit_stats  = groupby(unit_stats, [:animal,:day])


        for source in data_sources
            if source == "lfp"
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
            for key in keys(time_stats) # iterate over animal, day
                mt, nt, group = time_stats[key], 
                            unit_stats[NamedTuple(key)], 
                            groups[NamedTuple(key)]
                # 
                for timefield in time_vars
                    group[!,timefield] .-= mt.time_minimum      # center by behavior 0
                    if stagger_time
                        group[!,timefield] .+= mt.time_prev_maximum - mt.time_minimum # add previous max time of prev dataset
                    end
                end
                if stagger_units && :unit in propertynames(group)
                    @info "unit found" source key
                    group.unit .-= nt.unit_minimum      # center by behavior 0
                    group.unit .+= nt.unit_prev_maximum # add previous max time of prev dataset
                end
            end
            if stagger_time && :time in propertynames(datum)
                sort!(datum, [:time])
            else
                sort!(datum, [x for x in [:animal, :day, :time] 
                    if x in propertynames(datum)])
            end
            if source == "cells"
                @assert datum.old_unit != datum.unit
            end
            if source == "cells"
                savemethod(datum, "super$append", numsuperanim)
                savemethod(datum, "super$append", numsuperanim; 
                           type="arrow")
            else
                savemethod(datum, "super$append", numsuperanim)
            end
        end
        GC.gc()
        @info "be sure to run conversion to arrow if you change anything : tables_to_type()"
    end

    function fix_complex()
    end

end
