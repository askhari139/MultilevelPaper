topoFiles <- list.files(".", ".topo")
nets <- str_remove(topoFiles, ".topo")
BmodelSetup(".")
sapply(nets, function(net) {
    DirectoryNav(net)
    setwd("..")
    topoFile <- paste0(net, ".topo")
    file.copy(topoFile, paste0(net, "/", topoFile))
    file.copy("script_multiLevel.jl", paste0(net, "/script_multiLevel.jl"))
    setwd(net)
    system(paste0("julia -t 4 script_multiLevel.jl "))
    setwd("..")
})