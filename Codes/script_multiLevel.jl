nLevelList = 1:10
nInits = 100000
nIters = 1000

include("/Users/kishorehari/Desktop/Boolean.jl/bmodel.jl")
using Base.Threads

fileList = readdir()
topoFiles = String[]
for i in fileList
	if endswith(i, "topo")
		push!(topoFiles, i)
	end
end


println(Threads.nthreads())
# topoFile = "EMT_RACIPE2.topo"
Threads.@threads for topoFile in topoFiles
	for nLevels in nLevelList
		y1 = @elapsed x = bmodel_reps(topoFile; nInit = nInits, 
		  nIter = nIters, shubham = true, nLevels = nLevels)
		println(topoFile, " - nLevel ", nLevels, " : ", y1," seconds.")
	end
end
