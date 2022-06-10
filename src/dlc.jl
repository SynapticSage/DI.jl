module dlc
    using Glob, Printf
    using CSV, DataFrames
    using ..raw 
    function get_path(animal, day, epoch; dayfactor=0, 
            guessdayfactor=true, filtered=false, source="deeplabcut")
        if guessdayfactor
            dayfactor = raw.animal_dayfactor[animal]
        end
        day += dayfactor
        if source == "deeplabcut"
            folder_path = "/Volumes/Colliculus/deeplabcut/" * 
                         "goalmaze_tape-Ryan-2020-05-28/videos/"
        end
        possible_files = 
            glob("$(animal)_$(@sprintf("%02d",day))_$(@sprintf("%02d",epoch))_*.csv",
                 folder_path)
        #print(possible_files)
        if filtered
            possible_files = [file for file in possible_files if
                              occursin(file,"filtered")]
        else
            possible_files = [file for file in possible_files if
                              !(occursin(file,"filtered"))]
        end
        videopath = possible_files[1]
        return videopath
    end
    function load(pos...; kws...)
        df = CSV.read(get_path(pos...;kws...), DataFrame; header=3, skipto=4,
                     csvkws...)
        transform!(df, :coords=>(x->x*(1/30))=>:time, 
                      [:x,:x_1]=>((a,b)->(a.+b)./2)=>:X,
                      [:y,:y_1]=>((a,b)->(a.+b)./2)=>:Y,
               )

        transform!(df, [:X, :Y]=>((X,Y)->sqrt.(X.^2 .+ Y.^2))=>:vec)
        transform!(df, :vec => (x->[0;diff(x)]) => :delta)
    end
end
export dlc
