topoFile = ""
turnOffNodes = Int64[]
nLevelList = [0,1,2]
nInits = 100000
nIters = 1000

include(ENV["BMODEL"]*"/bmodel.jl")
using Base.Threads

# topoFile = "EMT_RACIPE2.topo"
for nLevels in nLevelList
        y1 = @elapsed x = scanNodeTurnOff(topoFile; nInit = nInits, 
            nIter = nIters, nLevels = nLevels, tSet = turnOffNodes)
    println(topoFile, " - nLevel ", nLevels, " : ", y1," seconds.")
end
