signals <- function(topoDf) {
    nodes <- topoDf %>% select(Source, Target) %>% unlist %>% unique
    signals <- data.frame(nodes = nodes, inSource = nodes %in% topoDf$Source, 
                       inTarget = nodes %in% topoDf$Target) %>%
            filter(!(inTarget)) %>% pull(nodes)
    list(topoDf %>% filter(!(Source %in% signals | Target %in% signals)), signals)
}

getSignalNodes <- function(topoFile) {
    topoDf <- read_delim(topoFile, delim = " ")
    s <- signals(topoDf)
    sAll <- c(s[[2]])
    while(length(s[[2]]) > 0) {
        s <- signals(s[[1]])
        sAll <- c(sAll, s)
    }
    return(sAll)
}

folderSetup <- function(topoFile, dr) {
    net <- topoFile %>% str_remove(".topo")
    drSignal <- paste0(dr, "/", net)
    DirectoryNav(drSignal)
    file.copy(paste0(topoFolder, "/", topoFile), ".", overwrite = T)
    file.copy(paste0(topoFolder, "/", net, "_nodes.txt"), ".", overwrite = T)
    t <- file.copy(paste0(topoFolder, "/", net, ".teams"), ".", overwrite = T)
    if (!t) {
        getGsVec()
    }
}

checkLevels <- function(nLevelList) {
    fl <- list.files(".", "format.csv")
    fT <- str_detect(fl, "shubham")
    nLevelsPresent <- c()
    if (!any(fT) & length(fT)) {
        nLevelsPresent <- c(0)
    }
    fl <- fl[fT]
    nLevels <- fl %>% str_extract("shubham_\\d+") %>% str_remove("shubham_") %>% as.numeric
    nLevelsPresent <- nLevelsPresent %>% c(nLevels)
    setdiff(nLevelList, nLevelsPresent)
}

signalNodeTurnOff <- function(topoFile, nLevelList=c(0:10, 20), numThreads = 100) {
    folderSetup(topoFile, signalFold)
    nLevelList <- checkLevels(nLevelList)
    if (length(nLevelList) == 0) {
        return()
    }
    net <- topoFile %>% str_remove(".topo")
    file.copy(turnOffScript, "script.jl", overwrite = T)
    signalNodes <- getSignalNodes(topoFile)
    if (length(signalNodes) == 0) {
        return()
    }
    nodes <- readLines(paste0(net, "_nodes.txt"))
    signalIDs <- which(nodes %in% signalNodes)
    tOffNodes <- paste0("[", paste0(signalIDs, collapse = ","), "]")
    nLevs <- paste0("[", paste0(nLevelList, collapse = ","), "]")
    sc <- readLines("script.jl")
    sc[1] <- paste0("topoFile = \"", topoFile, "\"")
    sc[2] <- paste0("turnOffNodes = ", tOffNodes)
    sc[3] <- paste0("nLevelList = ", nLevs)
    writeLines(sc, "script.jl")
    system(paste0("julia -t ", numThreads, " script.jl"))
    fList <- list.files(".", "finFlagFreq.csv")
    sapply(fList, function(f){
        formatterGeneral(f, topoFile)
    })
}

singleNodeTurnOff <- function(topoFile, nLevelList=c(0:10, 20), numThreads = 100) {
    folderSetup(topoFile, singleNodeFold)
    nLevelList <- checkLevels(nLevelList)
    if (length(nLevelList) == 0) {
        return()
    }
    net <- topoFile %>% str_remove(".topo")
    file.copy(turnOffScanScript, "script.jl", overwrite = T)
    nodes <- readLines(paste0(net, "_nodes.txt"))
    nNodes <- length(nodes)
    tOffNodes <- paste0("[]")
    nLevs <- paste0("[", paste0(nLevelList, collapse = ","), "]")
    sc <- readLines("script.jl")
    sc[1] <- paste0("topoFile = \"", topoFile, "\"")
    # sc[2] <- paste0("turnOffNodes = ", tOffNodes)
    sc[3] <- paste0("nLevelList = ", nLevs)
    writeLines(sc, "script.jl")
    system(paste0("julia -t ", numThreads, " script.jl"))
    fList <- list.files(".", "finFlagFreq.csv")
    sapply(fList, function(f){
        file.rename(f, f %>% str_replace("_turnOff.*", "_turnOff_signal_finFlagFreq.csv"))
        formatterGeneral(f, topoFile)
    })
}

allNodeTurnOff <- function(topoFile, nLevelList=c(0:10, 20), numThreads = 100) {
    folderSetup(topoFile, allNodeFold)
    nLevelList <- checkLevels(nLevelList)
    if (length(nLevelList) == 0) {
        return()
    }
    net <- topoFile %>% str_remove(".topo")
    file.copy(turnOffScript, "script.jl", overwrite = T)
    nodes <- readLines(paste0(net, "_nodes.txt"))
    nNodes <- length(nodes)
    # tOffNodes <- paste0("[", paste0(1:nNodes, collapse = ","), "]")
    nLevs <- paste0("[", paste0(nLevelList, collapse = ","), "]")
    sc <- readLines("script.jl")
    sc[1] <- paste0("topoFile = \"", topoFile, "\"")
    # sc[2] <- paste0("turnOffNodes = ", tOffNodes)
    sc[3] <- paste0("nLevelList = ", nLevs)
    writeLines(sc, "script.jl")
    system(paste0("julia -t ", numThreads, " script.jl"))
    fList <- list.files(".", "finFlagFreq.csv")
    sapply(fList, function(f){
        formatterGeneral(f, topoFile)
    })
    setwd(mlDir)
}

getCostMetric <- function(df, costMetric) {
    if (costMetric == "hybridness") {
        d <- df %>% filter(Phenotype == "Hybrid") %>% 
            pull(Avg0) %>% sum
    }
    if (costMetric == "meanFrust") {
        d <- df %>% mutate(frust = Avg0*frust0) %>%
            pull(frust) %>% sum
    }
    if(costMetric == "terminalness") {
        d <- df %>% filter(Phenotype != "Hybrid") %>% 
            pull(Avg0) %>% sum
    }
    return(d)
}


sequentialTurnOff <- function(topoFile, nLevelList = c(0:10, 20), numThreads = 100) {
    costMetrics <- c("hybridness", "meanFrust", "terminalness")
    folderSetup(topoFile, sequentialFold)
    nLevelList <- checkLevels(nLevelList)
    if (length(nLevelList) == 0) {
        return()
    }
    net <- topoFile %>% str_remove(".topo")
    file.copy(turnOffScanScript, "script.jl", overwrite = T)
    nodes <- readLines(paste0(net, "_nodes.txt"))
    teams <- readLines(paste0(net, ".teams")) %>% str_split(",")
    nNodes <- length(nodes)
    nLevs <- paste0("[", paste0(nLevelList, collapse = ","), "]")
    sc <- readLines("script.jl")
    sc[4] <- paste0("nInits = 1000")
    sapply(costMetrics, function(costMetric) {
        sapply(nLevelList, function(nLev) {
            mV <- 10
            turnOffs <- c()
            turnOffIDs <- c()
            dList <- list()
            while(mV != 0 & length(turnOffs) < nNodes) {
                sc[1] <- paste0("topoFile = \"", topoFile, "\"")
                sc[2] <- paste0("turnOffNodes = ", ifelse(length(turnOffIDs), paste0("[", paste0(turnOffIDs, collapse = ","), "]"), "Int64[]"))
                sc[3] <- paste0("nLevelList = ", "[", nLev, "]")
                writeLines(sc, "script.jl")
                system(paste0("julia -t ", numThreads, " script.jl"))
                fl <- list.files(".", "scanNode")
                formatterGeneral(fl, topoFile, nLev)
                file.remove(fl)
                fl <- fl %>% str_replace(".csv", "_format.csv")
                df <- read_csv(fl, show_col_types = F)
                tNodes <- df %>% pull(turnOffNode) %>% unique
                metricVals <- sapply(tNodes, function(t) {
                    getCostMetric(df %>% filter(turnOffNode == t), costMetric)
                })
                # browser()
                nodeDf <- data.frame(turnOffNode = tNodes, metricVals = metricVals)
                df <- df %>% left_join(nodeDf, by = "turnOffNode")
                mV <- max(metricVals)
                tNode <- tNodes[which(metricVals == mV)]
                tNode <- sample(tNode, 1)
                tNode <- tNode %>% str_split("_") %>% unlist
                tNode <- tNode[length(tNode)]
                tID <- which(nodes == tNode)
                turnOffs <- c(turnOffs, tNode)
                turnOffIDs <- c(turnOffIDs, tID)
                df$Node <- tNode
                df$NodeID <- tID
                df[costMetric] <- mV
                dList <- c(dList, list(df))
                file.remove(fl)
            }
            df <- bind_rows(dList)
            DirectoryNav(costMetric)
            write_csv(df, paste0(net, ifelse(nLev == 0, "", paste0("_shubham_", nLev)), 
                "_seqTurnOff_format.csv"))
            setwd("..")
        })
    })
    
    setwd(mlDir)
}

turnOffSim <- function(nets, nLevelList = c(0:10, 20), numThreads = 100) {
    sapply(nets, function(net) {
        topoFile <- paste0(net, ".topo")
        signalNodeTurnOff(topoFile, nLevelList, numThreads)
        # singleNodeTurnOff(topoFile, nLevelList, numThreads)
        allNodeTurnOff(topoFile, nLevelList, numThreads)
        sequentialTurnOff(topoFile, nLevelList, numThreads)
    })
}