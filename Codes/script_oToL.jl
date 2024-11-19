topoFile = ""
nLevelList = 1:10
nInits = 100000
nIters = 1000

include(ENV["BMODEL"]*"/bmodel.jl")
using Base.Threads

# topoFile = "EMT_RACIPE2.topo"
Threads.@threads for nLevels in nLevelList
    if (nLevels == 0)
        y1 = @elapsed x = bmodel_reps(topoFile; nInit = nInits, 
              nIter = nIters, shubham = false)
    else
        y1 = @elapsed x = bmodel_reps(topoFile; nInit = nInits, 
            nIter = nIters, shubham = true, nLevels = nLevels)
    end
    println(topoFile, " - nLevel ", nLevels, " : ", y1," seconds.")
end
