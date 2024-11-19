source("/Users/kishorehari/Desktop/Postdoc/MultiLevel/CoreData/setupScript.R")
source(paste0(mlCore, "/functions.R"))
source(paste0(mlCore, "/FigurePlots.R"))

setwd(signalTurnOff)

freqFiles <- list.files(".", "_finFlagFreq.csv")
topoFiles <- list.files(".", ".topo")
topoFins <- lapply(topoFiles, function(t) {
    net <- str_remove(t, ".topo")
    fl <- freqFiles[freqFiles %>% str_detect(paste0(net, "_sh"))]
    sapply(fl, function(f) {
        print(f)
        formatterGeneral(f, t, compute = T)
    })
})
