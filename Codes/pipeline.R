source("~/MultiLevelPaper/codes/setupScript.R") # Change it to the appropriate folder
source(paste0(codeFolder, "/functions.R"))
source(paste0(codeFolder, "/FigurePlots.r"))
source(paste0(codeFolder, "/FigureFunctions.r"))
source(paste0(codeFolder, "/turnOffFuncs.r"))
numThreads <- 4 # Change to the appropriate number of threads. Required for simulations with large initial conditions

generateFigures <- function(net, levelVec = c(1:5, 100)) { 
  # dir.create("Figures")
  ### Simulate and Format ------
  # Formatting the simulation output. For RACIPE, copy the simulation output from previous data.
  pwd <- getwd()
  setwd(paste0(multiLevelSim, "/", net))
  fileList <- sapply(levelVec, function(l) {
    if (l == 0) {
      topoFile <- paste0(net, ".topo")
      formatter(topoFile, l)
      fl <- paste0(net, "_finFlagFreq_format.csv")
      return(fl)
    }
    fl <- paste0(net, "_shubham_", l, "_finFlagFreq.csv")
    if (!file.exists(fl)) {
      print(paste0("File does not exist for level ", l, ". Simulate!"))
      return(NA)
    }
    topoFile <- paste0(net, ".topo")
    formatter(topoFile, l)
    fl <- paste0(net, "_shubham_", l, "_finFlagFreq_format.csv")
    return(fl)
  })
  if (any(is.na(fileList))) {
    print("Some files are missing. Check the messages above and Simulate!")
    return(NA)
  }
  
  ### Plotting ------
  
  # segregatedPlots(net, levelVec)
  Figure2(net)
  Figure2S1(net)
  Figure3(net)
  Figure4(net)
  Figure5(net)
  setwd(pwd)
}

topoFiles <- c("EMT_RACIPE2_noPeri.topo", "EMT_RACIPE2.topo"#, 
          # "EMT_RACIPE.topo", "EMT_RACIPE2.topo", 
          #"Gonadal_noPeri.topo", 
          # "melanoma_noPeri.topo", 
          #"Pluripotency_noPeri.topo", 
          # "EMT_MET_frust_noPeri.topo", 
          # "SIL_noPeri.topo", "SIL2_noPeri.topo",
          #"ThreeTeams_0.3.topo", "ThreeTeams_0.9.topo"
          )
nets <- topoFiles %>% str_remove(".topo")
simulateNetworks(nets, levelVec = c(1:10), largeIC = F, juliaThreading = F, 
   numThreads = numThreads)
turnOffSim(nets, nLevelList = c(1:4,10), numThreads = numThreads)
sapply(nets, generateFigures, levelVec = c(1:10), new = T)
# sapply(nets, generateFigures, levelVec = c(0:10))
# sapply(nets, generateFigures, levelVec = c(0:10, 100), new = F)
# generateFigures("EMT_RACIPE", levelVec = c(0:5, 100), new = F)
