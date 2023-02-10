# Labels
module Labels
import DataStructures: OrderedDict
corerr = correct = Dict(0=>"correct", 1=>"error")
cuemem = tsk = Dict(-1=>"nontask", 0=>"cue",     1=>"mem")
cortsk = OrderedDict([0,1]=>"CUE correct", [1,1]=>"MEM correct",
                  [0,0]=>"CUE error",   [1,0]=>"MEM error")
end
