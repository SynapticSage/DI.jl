module superanimal
    import ..Load
    using Infiltrator
    using DataFrames

    export make_superanimal
    """
    make_superanimal

    #params

    """
    function make_superanimal(animal_day_pairs::Tuple{String, Int}...;
            data_sources=Load.total_set, tag=nothing, numsuperanim::Int=0)
        for source in data_sources
            #data=[]
            loadmethod = Load.load_functions[source]
            savemethod = Load.save_functions[source]
            for (i,(animal, day)) in enumerate(animal_day_pairs)
                datum = loadmethod(animal, day; type="arrow")
                datum[!,:animal] .= animal
                datum[!,:day]    .= day
                type = i == 1 ? "arrow" : "arrow-append"
                savemethod(datum, "super", numsuperanim; type)
            end
        end
    end

    """
    center_times

    centers the times in the stored data per animal/day
    """
    function center_and_stagger_times_and_neurons(numsuperanim::Int=0; 
            data_sources=Load.total_set, append="", 
            stagger_units::Bool=true,stagger_time::Bool=true)

        beh, cells = Load.load_behavior("super",numsuperanim) ,
                     Load.load_cells("super", numsuperanim)
        time_stats = combine(groupby(beh, [:animal,:day]), 
                           :time=>minimum, :time=>maximum)

        time_stats[:,:time_maximum] .= [0; time_stats[1:end-1,:time_maximum]]
        time_stats  = transform(time_stats, 
                              :time_maximum=>cumsum=>:time_prev_maximum)
        time_stats  = groupby(time_stats, [:animal,:day])

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
            loadmethod = Load.load_functions[source]
            savemethod = Load.save_functions[source]
            time_vars  = Load.time_vars[source]
            time_vars  = time_vars isa Vector ? time_vars : [time_vars]
            datum = loadmethod("super", numsuperanim)
            groups = groupby(datum, [:animal,:day])
            for key in keys(time_stats)
                @infiltrate
                mt, nt, group = time_stats[key], 
                            unit_stats[NamedTuple(key)], 
                            groups[NamedTuple(key)]
                for timefield in time_vars
                    group[!,timefield] .-= mt.time_minimum      # center by behavior 0
                    if stagger_time
                        group[!,timefield] .+= mt.time_prev_maximum # add previous max time of prev dataset
                    end
                end
                if stagger_units && :unit in propertynames(group)
                    group.unit .-= nt.unit_minimum      # center by behavior 0
                    group.unit .+= nt.unit_prev_maximum # add previous max time of prev dataset
                end
                @infiltrate
            end
            @infiltrate
            savemethod(datum, "super$append", numsuperanim)
        end
    end
end
