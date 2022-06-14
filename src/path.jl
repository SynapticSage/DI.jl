module path

    using DrWatson

    import CSV
    using DataFrames

    using Plots, Measures
    using Statistics
    using ImageFiltering
    using ColorSchemes
    
    import ..Load

    function load_pathtable(animal, day)
        f=CSV.read(datadir("paths.csv"), DataFrame; csvkws...)
        init = Dates.Time("00:00:00")
        function s(x)
            x = Second.(x .- init)
            x = [e.value for e in x]
        end
        transform!(f, 
                   :start=>(x->s(x))=>:start, 
                   :end=>(x->s(x))=>:end, 
                   [:end,:start]=>((e,s)->Dates.Second.(e.-s))=>:duration)
        return f
    end

    export load_pathtable

    function get_frame(video,sub)
        frame = Load.video.frameattime(video, median(sub[sub.delta .> 1,:].time))
        frame = frame'
        diff1, diff2 = ImageFiltering.Kernel.sobel()
        frameX=ImageFiltering.imfilter(frame, diff1)
        frameY=ImageFiltering.imfilter(frame, diff2)
        frame = (1.2.*frame) .+ (3 .* (frameX .+ frameY))
    end


    function get_path(position, path)
        sub = position[(position.time .>= path.start[1]) .& 
                       (position.time .<= path.end[1]) .&
                       (position.likelihood .> 0.5) .&
                       (position.delta .< 6) ,:]
    end

    function subplot_block(position, paths; video=nothing)
        views = []
        for path in eachrow(paths)
            color = get(ColorSchemes.coolwarm, (path.btraj[1]-1)/8)
            sub = get_path(position, path)
            if isempty(sub)
                print("sub is empty")
                continue
            end
            if path.to == "a"
                m = :rtriangle
            else
                m = :ltriangle
            end
            if path.btraj[1] <= 4
                cm = "cue"
            else
                cm = "mem"
            end
            title = "$(path.btraj[1]) $(path.to[1]) $cm"
            if video != nothing
                frame = get_frame(video, sub)
                h=Plots.heatmap(frame, framestyle=:none)
                Plots.plot!(sub.X, sub.Y, c=color, m=m, markersize=1,
                            markerstrokewidth=0,
                            markerstrokealpha=0,
                            linestyle=:dot, title=title,
                            markercolor=color, label=nothing,
                            alpha=(sub.likelihood))
            else
                h=Plots.plot(sub.X, sub.Y, c=color, m=m, markersize=1,
                             markerstrokewidth=0,
                             markerstrokealpha=0,
                             linestyle=:dot, title=title, titlefontsize=1,
                             markercolor=color, label=nothing,
                             alpha=(sub.likelihood))
            end
            h.attr[:dpi] = 300
            push!(views, h)
        end
        #push!(views, Plots.plot([NaN],title="Dynamic RY22 (DLC needs more training)"))
        P = Plots.plot(views..., margin=-4mm, link=:all, titlefontsize=3)
        P.attr[:dpi] = 300
        return P, views
    end

    function plot_block(position, paths; video=nothing)
        P = Plots.plot(framestyle=:none)
        for (i,path) in enumerate(eachrow(paths))
            color = get(ColorSchemes.coolwarm, (path.btraj[1]-1)/8)
            sub = get_path(position, path)
            if isempty(sub)
                print("sub is empty")
                continue
            end
            if path.to == "a"
                m = :rtriangle
            else
                m = :ltriangle
            end
            if (i == 1) && (video != nothing)
                P = Plots.heatmap!(get_frame(video, sub), framestyle=:none)
            end
            Plots.plot!(P, sub.X, sub.Y, label=path.btraj[1], c=color,
                        markercolor=color, alpha=(0.8.*sub.likelihood),
                        markerstrokewidth=0, markerstrokealpha=0,
                        linestyle=:dot, linewidth=0.5, markersize=2, m=m)
            #for (x,y) in zip(sub.X, sub.Y)
            #    Plots.annotate!(P, x, y, text(path.btraj[1], 3))
            #end
        end
        P.attr[:dpi] = 300
        P
    end

    function newvar()
    end

end
