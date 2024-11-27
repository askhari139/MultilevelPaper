library(funcsKishore)
library(compiler)
library(tictoc)
library(cowplot)

BmodelFolder <- Sys.getenv("BMODEL")
# mlDir <- "~/MultiLevelPaper" - this should come from the pipline.R script

simFolder <- paste0(mlDir, "/SimData")
multiLevelSim <- paste0(simFolder, "/MultiLevel")
turnOffSim <- paste0(simFolder, "/TurnOff")
signalFold <- paste0(turnOffSim,"/Signals")
sequentialFold <- paste0(turnOffSim,"/Sequential")
singleNodeFold <- paste0(turnOffSim,"/SingleNode")
allNodeFold <- paste0(turnOffSim,"/AllNodes")


codeFolder <- paste0(mlDir, "/Codes")
topoFolder <- paste0(mlDir, "/TopoFiles")

figures <- paste0(mlDir, "/Figures")


### julia scripts

oneTopoMultiLevel <- paste0(codeFolder, "/script_oTmL.jl")
multiTopoOneLevel <- paste0(codeFolder, "/script_mToL.jl")
oneTopoOneLevel <- paste0(codeFolder, "/script_oToL.jl")
turnOffScript <- paste0(codeFolder, "/script_turnOff.jl")
turnOffScanScript <- paste0(codeFolder, "/script_turnOffScan.jl")

### create folders

folderList <- c(simFolder, multiLevelSim, turnOffSim, signalFold, 
    sequentialFold, singleNodeFold, allNodeFold, codeFolder, topoFolder, figures)

sapply(folderList, function(f) {
    if (!dir.exists(f)) {
        dir.create(f)
    }
})
