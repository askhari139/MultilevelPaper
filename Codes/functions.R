peripherals <- function(topoDf) {
  nodes <- topoDf %>% select(Source, Target) %>% unlist %>% unique
  peripherals <- data.frame(nodes = nodes, inSource = nodes %in% topoDf$Source, 
                       inTarget = nodes %in% topoDf$Target) %>%
            filter(!(inSource & inTarget)) %>% pull(nodes)
  list(topoDf %>% filter(!(Source %in% peripherals | Target %in% peripherals)), peripherals)
}
removePeripherals <- function(topoFile) {
  topoDf <- read.delim(topoFile, sep = "")
  ls <- peripherals(topoDf)
  while(length(ls[[2]]) > 0) {
    topoDf <- ls[[1]]
    ls <- peripherals(topoDf)
  }
  write_delim(topoDf, topoFile %>% str_replace(".topo", "_noPeri.topo"), delim = " ")
}


simulateNetworks <- function(topoFiles, levelVec = c(1:10, 20, 50, 100), 
                             largeIC = T, juliaThreading = F, numThreads = 100) {
  nets <- topoFiles %>% str_remove(".topo")
  sapply(nets, function(net) {
    # browser()
    topoFile <- paste0(net, ".topo")
    ### Simulations ---------
    ### Multi level 
    DirectoryNav(paste0(multiLevelSim, "/", net))
    if (!file.exists(paste0(net, ".topo"))) {
      file.copy(paste0(topoFolder, "/", net, ".topo"), 
                paste0(net, ".topo"))
      
    }
    if (!file.exists(paste0(net, ".teams"))) {
      t <- file.copy(paste0(topoFolder, "/", net, ".teams"), 
                paste0(net, ".teams"))
      if (!t) {
        getGsVec()
      }
    }
    levelsAll <- levelVec
    filz <- paste0(net, "_shubham_", levelsAll, "_finFlagFreq.csv")
    levelsLeft <- levelsAll[!(filz %>% file.exists)]
    
    if (largeIC) {
      if (!dir.exists("simulation")) {
        dir.create("simulation")
      }
      sapply(1:100, function(i) {
        file.copy(paste0(topoFolder, "/", net, ".topo"), 
                  paste0("simulation/t_", str_pad(i, 3, pad = "0"), ".topo"))
      })
      DirectoryNav("simulation")
      if (juliaThreading) {
        for (nLevels in levelsLeft) {
          file.copy(multiTopoOneLevel, "script_mToL.jl")
          sc <- readLines("script_mToL.jl")
          sc[1] <- paste0("nLevels = ", nLevels)
          sc[3] <- "nIters = 1000"
          writeLines(sc, "script_mToL.jl")
          system(paste0("julia -t ", numThreads, " script_mToL.jl"))
        }
      }
      else {
        file.copy(oneTopoMultiLevel, "script.jl")
        sc <- readLines("script.jl")
        sc[2] <- paste0("nLevelList = ", "[", paste0(levelsLeft, collapse = ","), "]")
        sc[3] <- "nInits = 100000"
        topoFiles <- list.files(".", pattern = ".topo")
        plan(multisession, workers = numThreads)
        future_sapply(topoFiles, function(topoFile) {
          sc[1] <- paste0("topoFile = \"", topoFile, "\"")
          writeLines(sc, str_replace(topoFile, ".topo", "_script.jl"))
          system(paste0("julia ", str_replace(topoFile, ".topo", "_script.jl")))

        }, future.seed = T)
      }
      for(nLevels in levelVec) {
        numSimFiles <- list.files(".", paste0(nLevels, "_finFlagFreq"))
        # while(T) {# to ensure all files are present
          
        #   if (length(numSimFiles) == 100) {
        #     break
        #   }
        #   Sys.sleep(1)
        #   iter <- iter + 1
        #   if (iter > 100) break
        # }
        if (length(numSimFiles) == 0) {
          print("Simulations did not happen for" %>% paste0(net, "_", nLevels))
          setwd("..")
          next
        }
        df <- lapply(numSimFiles, function(f) {
          read_csv(f, show_col_types = F)
        }) %>% bind_rows %>%
          group_by(states, flag) %>%
          summarise(Avg0 = sum(Avg0), 
                    frust0 = mean(frust0),
                    time = mean(time), 
                    .groups = "drop") %>%
          mutate(Avg0 = Avg0/sum(Avg0))

        setwd("..")
        write_csv(df %>% arrange(desc(Avg0)), paste0(net, "_shubham_", nLevels, "_finFlagFreq.csv"))
      }
        
      
      file.copy(paste0("../", net, "_nodes.txt"), paste0(net, "_nodes.txt"))
    }
    else {
      file.copy(oneTopoMultiLevel, "script.jl")
      sc <- readLines("script.jl")
      sc[1] <- paste0("topoFile = \"", topoFile, "\"")
      sc[2] <- paste0("nLevelList = ", "[", paste0(levelsLeft, collapse = ","), "]")
      sc[3] <- "nInits = 100000"
      writeLines(sc, "script.jl")
      system(paste0("julia -t ", numThreads, " script.jl"))
    }
  })
}

frustCalc <- function(state, nodeOrder, topoDf, nLevels = 1, levelKey = c(-1,1))
{#browser()
    state <- str_split(state, "_") %>% unlist %>% as.integer
    state <- levelKey[state + 1]
    names(state) <- nodeOrder
    nEdges <- nrow(topoDf)
    frust <- topoDf %>%
        mutate(Source = state[Source], Target = state[Target], 
            Type = ifelse(Type == 2, -1, 1)) %>%
        mutate(Sum = Source*Target*Type) %>%
        mutate(frust = ifelse(Sum<=0, Sum, 0)) %>%
        select(frust) %>%
        unlist %>% sum
    frust/nEdges
}
frustCalc <- cmpfun(frustCalc)

numConverter <- function(n, numLevels) {
    nodeLevels <- 1:numLevels
    nodeCompare <- c(-1*rev(nodeLevels), 0, nodeLevels)/numLevels
    nodeLevels <- c(-1*rev(nodeLevels), nodeLevels)/numLevels
    # discritize the number line into levels
    n <- as.numeric(n)
    n1 <- 1
    s <- sum(n > nodeCompare)
    if (s > 2*numLevels) n1 <- 1
    else {
        if (s <= 1) n1 <- -1
        else n1 <- nodeLevels[s]
    }
    return(n1)
}
numConverter <- cmpfun(numConverter)

isingConvert <- function(state, numLevels = 2) {
    if (numLevels == 0 || numLevels == 1)
        return(state)
    nLevs <- 1:(2*numLevels)
    isingLevs <- c(rep("0", numLevels), rep("1", numLevels))
    names(isingLevs) <- nLevs
    state <- state %>% str_split("_") %>% unlist %>% as.integer
    state <- isingLevs[state + 1] %>% paste0(collapse = "_")
    return(state)
}
isingConvert <- cmpfun(isingConvert)

formatterGeneral <- function(f, topoFile, nLevels = NULL, compute = F) {
    turnOff <- ifelse(str_detect(f, "_turnOff"), T, F)
    if (is.null(nLevels)) {
        nLevels <- ifelse(str_detect(f, "_shubham_"), 
            as.numeric(str_extract(f, "_shubham_(\\d+)_") %>% 
                str_remove("_shubham_") %>% str_remove("_")), 1)
    }
    levelKey <- c(-1*rev(1:nLevels), 1:nLevels)/nLevels
    # browser()
    if(turnOff) levelKey <- c(levelKey, 0)
    formattedFile <- str_replace(f, ".csv", "_format.csv")
    if (file.exists(formattedFile) && !compute) {
        print(paste0("Formatter file exists ", formattedFile,"."))
        return(0)
    }
    net <- topoFile %>% str_remove(".topo")
    nodes <- net %>%
        paste0("_nodes.txt") %>%
        readLines
    topoDf <- read.delim(topoFile, sep = "")
    freqFile <- read_csv(f, show_col_types = F) %>%
        separate(states, c(nodes), sep = "_", remove = F)
    teams <- paste0(net, ".teams")
    if (!file.exists(teams)) {
        getGsVec()
    }
    teams <- readLines(paste0(net, ".teams")) %>%
        str_split(",")
    eTeam <- teams[[1]]
    mTeam <- teams[[2]]
    freqFile <- freqFile %>%
        mutate(across(all_of(nodes), .fns = function(x) {
            x <- x %>% as.integer
            levelKey[x + 1]
        }))
    freqFile <- freqFile %>%
        mutate(eScore = freqFile %>% select(all_of(eTeam)) %>% rowMeans,
               mScore = freqFile %>% select(all_of(mTeam)) %>% rowMeans) %>%
        mutate(emScore = eScore-mScore) %>%
        mutate(Phenotype = case_when(
          abs(emScore) == 2 ~ "Terminal",
          eScore == 1 | mScore == 1 ~ "Incomplete Terminal",
          .default = "Hybrid"
        ))
    write_csv(freqFile, formattedFile, quote = "none")
}

formatter <- function(topoFile, nLevels = 2, compute = F) {
    net <- topoFile %>% str_remove(".topo")
    levelKey <- c(-1*rev(1:nLevels), 1:nLevels)/nLevels
    freqFile <- paste0(net, ifelse(nLevels == 0, "", paste0("_shubham_", nLevels)), 
      "_finFlagFreq.csv")
    formattedFile <- str_replace(freqFile, ".csv", "_format.csv")
    if (file.exists(formattedFile) && !compute) {
        print(paste0("Formatter file exists for nLevels ", nLevels,"."))
        return(0)
    }
    if(!file.exists(freqFile)) {
        print(paste0("File not found for nLevels ", nLevels,". Simulate!"))
        return(NA)
    }
    nodes <- net %>%
        paste0("_nodes.txt") %>%
        readLines
    topoDf <- read.delim(paste0(net, ".topo"), sep = "")
    teams <- readLines(paste0(net, ".teams")) %>%
      str_split(",")
    nTeams <- length(teams)
    if (nTeams == 2) {
        eTeam <- teams[[1]]
        mTeam <- teams[[2]]
    }

    stateLevels <- c(-1*rev(1:nLevels), 1:nLevels)/nLevels
    if (nLevels == 0) stateLevels <- c(-1, 1)
    teamNodes <- unlist(teams)
    nodesNew <- c(teamNodes, nodes[!(nodes %in% teamNodes)])
    
    freqFile <- freqFile %>%
        read_csv(show_col_types = F) %>%
        separate(states, c(nodes), sep = "_", remove = F) %>%
        mutate(across(all_of(nodesNew), .fns = function(x) {
          x <- x %>% as.integer
          stateLevels[x + 1]
        })) %>%
        mutate(State = paste0("S", rownames(.)))
    if (nTeams == 2)
    freqFile <- freqFile %>%
        mutate(eScore = freqFile %>% select(all_of(eTeam)) %>% rowMeans,
               mScore = freqFile %>% select(all_of(mTeam)) %>% rowMeans) %>%
        mutate(emScore = eScore-mScore) %>%
        mutate(Phenotype = case_when(
          abs(emScore) == 2 ~ "Terminal",
          eScore == 1 | mScore == 1 ~ "Incomplete Terminal",
          .default = "Hybrid"
        ))
      else {
        for (i in 1:nTeams) {
          freqFile[[paste0("Team_", i, "_Score")]] <- freqFile %>% 
            select(all_of(teams[[i]])) %>% rowMeans
        }
        freqFile <- freqFile %>%
          mutate(Sum = freqFile %>% select(ends_with("_Score")) %>% rowSums, 
          Prod = freqFile %>% select(ends_with("_Score")) %>% rowSums) %>%
          mutate(Phenotype = "Hybrids")
        for (k in 1:nTeams) {
          freqFile$Phenotype[freqFile$Prod == 0 & freqFile$Sum == k] <- paste0(k, " Teams")
        }
      }
    write_csv(freqFile,
              formattedFile,
              quote = "none")
    return(0)
}
formatter <- cmpfun(formatter)

headNtail <- function(df, n = 5) {
    if (nrow(df) <= 2*n) return(df)
    rbind.data.frame(head(df, n), tail(df, n))
}
