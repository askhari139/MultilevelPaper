topoFile = ""
turnOffNodes = Int64[]
nLevelList = [0,1,2]
nInits = 100000
nIters = 1000

include(ENV["BMODEL"]*"/bmodel.jl")
using Base.Threads

# topoFile = "EMT_RACIPE2.topo"
for nLevels in nLevelList
    if (nLevels == 0)
        y1 = @elapsed x = bmodel_reps(topoFile; nInit = nInits, 
              nIter = nIters, shubham = false, vaibhav = true,
              turnOffNodes = turnOffNodes)
    else
        y1 = @elapsed x = bmodel_reps(topoFile; nInit = nInits, 
            nIter = nIters, shubham = true, nLevels = nLevels, 
            vaibhav = true, turnOffNodes = turnOffNodes)
    end
    println(topoFile, " - nLevel ", nLevels, " : ", y1," seconds.")
end
