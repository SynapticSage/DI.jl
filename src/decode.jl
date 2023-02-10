
using Printf
using NetCDF
export decodedir, decodepath, load_decode

function decodedir(;method::String="sortedspike", transition="empirical",
        binstate="notbinned", n_split=4, downsamp=1, speedup=1.0)
    base = "/Volumes/FastData/decode/"
    paramfolder = "$method.$transition.$binstate.n_split=$n_split.downsamp=$downsamp.speedup=$(@sprintf("%1.1f", speedup))"
    return joinpath(base, paramfolder)
end
function decodepath(animal="RY16", day=36, epoch=7; type="test", split=1, ext::String="nc", kws...)
    dir = decodedir(;kws...)
    file = "$(animal)decode$(@sprintf("%02d",day))-$(@sprintf("%02d", epoch))split=$split.$type.$ext"
    fullpath = joinpath(dir, file)
    if !(isfile(fullpath))
        @warn "file=$fullpath does not exist"
    end
    return fullpath
end

function load_decode(filename)
    @info "Loading $filename"
    v = NetCDF.open(filename)
    if "Var1" in keys(v.vars)
        v.vars["time"] = v.vars["Var1"]
        pop!(v.vars, "Var1")
    end
    decode = Dict(var => Array(v.vars[var]) 
               for var in keys(v.vars))
    return decode
end
