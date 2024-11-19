getBreaks <- function(vec, nBreaks = 5) {
    # get nBreaks equidistant items in a vector
    n <- length(vec)
    if (n < nBreaks) return(vec)
    breaks <- seq(1, n, length.out = nBreaks) %>% round
    # breaks <- breaks[1:(nBreaks-1)]
    return(vec[breaks])
}

segmentPlots <- function(df, relative = F) {
    # browser()
    df <- df %>% 
        # mutate(nLevels = ifelse(nLevels == "RACIPE", 11, nLevels)) %>%
        mutate(nLevels = as.integer(as.character(nLevels))) %>%
        filter(!is.na(Avg0), !is.na(nLevels)) %>%
        select(emScore, nLevels, Avg0, Phenotype) %>%
        mutate(PhenThreshSum = Phenotype, emScoreSum = emScore) %>%
        mutate(yBeg = round(emScoreSum, 2)) %>%
        group_by(yBeg, nLevels, emScoreSum) %>%
        summarise(Avg0 = sum(Avg0), 
            Phenotype = unique(Phenotype), 
            emScore = mean(emScore), .groups = "drop") %>%
        # filter(Phenotype == "Hybrid") %>%
        mutate(xBeg = nLevels - 2, xEnd = nLevels, yEnd = yBeg)
    if (relative) {
        df <- df %>%
            group_by(nLevels) %>%
            mutate(yBeg = 2*yBeg/max(abs(yBeg)), yEnd = yBeg)
    }
    turnOffLevels <- df$nLevels %>% unique %>% as.integer %>% sort
    df <- df %>% mutate(nLevels = factor(nLevels, levels = turnOffLevels))
    write_csv(df, "seqTurnOff_segment.csv")
    p <- ggplot()
    sapply(turnOffLevels, function(t) {
        p <<- p + geom_segment(data = df %>% filter(nLevels == t), 
            aes(x = xBeg, xend = xEnd, y = yBeg,
            yend = yEnd, color = log10(Avg0)), size = 2)
    })
    p <- p + theme_Publication() + 
    theme(legend.position = "right", 
        legend.direction = "vertical", legend.key.height = unit(0.8, "cm"),
        axis.text = element_text(size = rel(1.2)), 
        axis.title = element_text(size = rel(1.5))) +
    scale_color_viridis_c() +
    labs(x = "Number of levels", 
        y = ifelse(relative, "Relative Score", "Score"), 
        color = "Log\nFrequency")
    return(p)
}

stateExpressionHeatmap <- function(dat, nodeOrder, key, maxStates = 100,
    freqWidth = F, returnFixed = T) {
    # browser()
    dat <- dat %>%
            arrange(emScore, desc(Avg0)) %>%
            mutate(stateNum = 1:nrow(.))
    d <- dat %>% #filter(stateNum <= maxStates) %>%
            select(all_of(nodeOrder), stateNum, Avg0, emScore) %>%
            gather(key = "Node", value = "State", -Avg0, -emScore, -stateNum) %>%
            mutate(Node = factor(Node, levels = nodeOrder)) %>%
            arrange(desc(stateNum)) %>%
            mutate(height = log10(Avg0)) %>%
            mutate(height = (height-min(height))/(max(height)-min(height))) %>%
            mutate(height = 0.1 + 0.8*height)
    d1 <- d %>% filter(stateNum <= maxStates)
    p <- ggplot(d1, aes(x = Node, y = stateNum, fill = State))
    if (freqWidth) {
        p <- p + geom_tile(height = d1$height)
    }
    else {
        p <- p + geom_tile()
    }
    p1 <- p +
        theme_Publication() +
        labs(x = "", y = "", fill = "State") +
        scale_fill_gradient2(low = "blue", mid = "white", 
            high = "red", midpoint = 0) +
        scale_x_discrete(expand = c(0,0)) +
        scale_y_continuous(expand = c(0,0)) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
                legend.position = "right", legend.direction = "vertical",
                legend.key.height = unit(0.8, "cm"), 
                axis.text.y = element_blank())
    ggsave(key, 
        width = 8, height = 2.5 + 0.2*max(d1$stateNum), limitsize = F)

    p <- ggplot(d, aes(x = Node, y = stateNum, fill = State))
    if (freqWidth) {
        p <- p + geom_tile(height = d$height)
    }
    else {
        p <- p + geom_tile()
    }
    p <- p +
        theme_Publication() +
        labs(x = "", y = "", fill = "State") +
        scale_fill_gradient2(low = "blue", mid = "white", 
            high = "red", midpoint = 0) +
        scale_x_discrete(expand = c(0,0)) +
        scale_y_continuous(expand = c(0,0), 
            n.breaks = min(length(dat$stateNum), 10)) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
                legend.position = "right", legend.direction = "vertical",
                legend.key.height = unit(0.8, "cm"))
    ggsave(key %>% str_replace(".png", "_fixed.png"), 
        width = 9, height = 5.2)
    return(p)
}

seqTurnPlot <- function(topoFile, costMetric, nLevels = c(1:2)) {
    net <- topoFile %>% str_remove(".topo")
    dr <- paste0(sequentialFold, "/", net, "/", costMetric)
    if (!dir.exists(dr)) 
        print("Invalid directory. Either no simulation or wrong cost metric")
    pList <- sapply(nLevels, function(nLev) {
        fl <- paste0(dr, "/", net, "_shubham_", nLev, "_seqTurnOff_format.csv")
        if (!file.exists(fl)) {
            print(paste0("File does not exist for level ", nLev, ". Simulate!"))
            return(NA)
        }
        df <- read_csv(paste0(dr, "/", net, "_shubham_", nLev, 
            "_seqTurnOff_format.csv"), show_col_types = F)
        nodes <- df$Node %>% unique
        d <- lapply(1:length(nodes), function(i) {
            node <- nodes[i]
            d <- df %>% filter(Node == node) %>%
                filter(str_detect(turnOffNode, paste0(node, "$"))) %>%
                mutate(ID = i) %>%
                filter(flag == 1 | Avg0 > 0.1*max(Avg0))
            d
        }) %>% bind_rows %>%
            mutate(ID = factor(ID, levels = 1:length(nodes)))
        p <- ggplot(d, aes(x = ID, y = Avg0, fill = Phenotype)) +
            geom_point(aes(color = Phenotype), size = 1.5, 
                position = position_dodge(width = 0.9)) +
            geom_boxplot() +
            
            theme_Publication() +
            labs(x = "# Nodes turned off", y = "SSF") +
            scale_y_log10() +
            guides(color = F) +
            theme(legend.position = "top", legend.text = element_text(size = rel(1)), 
                legend.title = element_text(size = rel(1)), 
                axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))    
        return(list(p))    
    })
    return(pList)
}

Figure2 <- function(net) {
    figFolder <- paste0(figures, "/Figure2/", net)
    if (!dir.exists(figFolder)) dir.create(figFolder, recursive = T)
    ### Expression heatmaps ------
    fls <- c(paste0(net, "_shubham_2_finFlagFreq_format.csv"), 
        paste0(net, "_shubham_0_finFlagFreq_format.csv"))
    nms <- c("4", "2")
    names(fls) <- nms
    teams <- readLines(paste0(multiLevelSim, "/", net, "/", net, ".teams")) %>% 
        str_split(",") %>% unlist
    tAll <- unlist(teams)
    topoDf <- read.delim(paste0(topoFolder, "/", net, ".topo"), sep = "")
    nodes <- c(topoDf$Source, topoDf$Target) %>% unique
    nodeOrder <- c(tAll, nodes[!nodes %in% tAll])
    pList <- list()
    df <- lapply(nms, function(nm) {
        fl <- fls[nm]
        d <- read_csv(paste0(multiLevelSim, "/", net, "/", fl), 
            show_col_types = F) %>%
            filter(flag == 1 | Avg0 > 0.05*max(Avg0))
        
        p1 <- stateExpressionHeatmap(d, nodeOrder, 
            paste0(figFolder, "/l_", nm, "_st_exp.png"))
        p2 <- ggplot(d, aes(x = eScore, y = mScore, color = Phenotype)) +
            geom_point(size = 2.3) +
            # geom_density2d() +
            labs(x = "E score", y = "M score") +
            theme_Publication() +
            theme(legend.position = c(0.8, 0.86), legend.direction = "vertical")
        ggsave(paste0(figFolder, "/l_", nm, 
            "_eScore_mScore.png"), width = 5.5, height = 5.2)
        pList <<- c(pList, list(p1))
        pList <<- c(pList, list(p2))
        d <- d %>% mutate(nLevels = paste0(nm, "-Levels"))
        d
    }) %>% bind_rows()
    dfSumm <- df %>% 
        group_by(nLevels, Phenotype) %>%
        summarise(nStates = n(), Freq = sum(Avg0), .groups = "drop") %>%
        mutate(Phenotype = str_replace(Phenotype, " ", "\n")) %>%
        mutate(Phenotype = factor(Phenotype, 
            levels = c("Terminal", "Incomplete\nTerminal", "Hybrid")))
    p3 <- ggplot(dfSumm, aes(x = Phenotype, y = nStates, fill = nLevels)) +
        geom_bar(stat = "identity", 
            position = position_dodge(preserve = "single")) +
        scale_y_log10() +
        theme_Publication() + 
        theme(legend.position = "top", legend.text = element_text(size = rel(1.1)), 
            legend.title = element_text(size = rel(1.1))) +
        labs(x = "Phenotype", y = "Number of states", fill = "Model")
    ggsave(paste0(figFolder, 
            "/nStates_levels.png"), width = 5.5, height = 5.4)
    p4 <- ggplot(dfSumm, aes(x = Phenotype, y = Freq, fill = nLevels)) +
        geom_bar(stat = "identity", 
            position = position_dodge(preserve = "single")) +
        theme_Publication() + 
        theme(legend.position = "top", legend.text = element_text(size = rel(1.1)), 
            legend.title = element_text(size = rel(1.1))) +
        labs(x = "Phenotype", y = "Total SSF", fill = "Model")
    ggsave(paste0(figFolder, 
            "/Freq_levels.png"), width = 5.5, height = 5.4)
    pList <- c(pList, list(p3))
    pList <- c(pList, list(p4))
    figFolder <- paste0(figures, "/Figure2")
    # browser()
    print("getting the compiled figure")
    png(paste0(figFolder, "/", net, ".png"), 
        width = 21, height = 10, units = "in", res = 300)
    print({plot_grid(pList[[3]], pList[[4]], pList[[5]], pList[[1]], pList[[2]], pList[[6]], 
        rel_widths = c(8/18, 5/18, 5/18, 8/18, 5/18, 5/18), 
        ncol = 3, 
        labels = c("A (i)", "B (i)", "C", "(ii)", "(ii)", "D"), 
        label_size = 20)})
    dev.off()
}

Figure2S1 <- function(net) {

}

Figure3 <- function(net) {
    # browser()
    figFold <- paste0(figures, "/Figure3/", net)
    if(!dir.exists(figFold)) dir.create(figFold, recursive = T)
    ### Expression heatmaps ------
    fls <- c(paste0(net, "_shubham_2_finFlagFreq_format.csv"), 
        paste0(net, "_shubham_0_finFlagFreq_format.csv"))
    nms <- c("4", "2")
    names(fls) <- nms
    teams <- readLines(paste0(multiLevelSim, "/", net, "/", net, ".teams")) %>% 
        str_split(",") %>% unlist
    tAll <- unlist(teams)
    topoDf <- read.delim(paste0(topoFolder, "/", net, ".topo"), sep = "")
    nodes <- c(topoDf$Source, topoDf$Target) %>% unique
    nodeOrder <- c(tAll, nodes[!nodes %in% tAll])
    pList <- list()
    df <- lapply(nms, function(nm) {
        fl <- fls[nm]
        d <- read_csv(paste0(multiLevelSim, "/", net, "/", fl), 
            show_col_types = F) %>%
            filter(flag == 1 | Avg0 > 0.05*max(Avg0))
        p1 <- ggplot(d, aes(x = emScore, y = frust0, color = Phenotype)) +
            geom_point(size = 2) +
            # geom_density2d() +
            labs(x = "Phenotype score", y = "Frustration") +
            theme_Publication() +
            theme(legend.position = "none")
        ggsave(paste0(figFold, "/l_", nm,
            "_emScore_frust.png"), width = 5.5, height = 5.2)
        p2 <- ggplot(d, aes(x = emScore, y = log10(Avg0), color = Phenotype)) +
            geom_point(size = 2) +
            # geom_density2d() +
            labs(x = "Phenotype score", y = "Log10(SSF)") +
            theme_Publication() +
            theme(legend.position = "none")
        ggsave(paste0(figFold, "/l_", nm,
            "_emScore_ssf.png"), width = 5.5, height = 5.2)
        pList <<- c(pList, list(p1))
        pList <<- c(pList, list(p2))
        d <- d %>% mutate(nLevels = paste0(nm, "-Levels")) %>% 
            filter(Phenotype == "Hybrid") %>%
            arrange(desc(Avg0))
        p3 <- stateExpressionHeatmap(d %>% head(5), nodeOrder, 
            paste0(figFold, "/l_", nm, "_hybrid_exp.png"))
        pList <<- c(pList, list(p3))
        d <- d %>% headNtail(5)
        d %>%
            mutate(Category = c(rep("High SSF", min(5, nrow(d))), rep("Low SSF", max(0, nrow(d)-5))))
    }) %>% bind_rows()
    p4 <- ggplot(df, aes(x = log10(Avg0), y = frust0, color = nLevels, shape = Category)) +
        geom_point(size = 3) +
        labs(x = "Log10(SSF)", y = "Frustration", color = "Model", shape = "") +
        theme_Publication() +
        theme(legend.position = "top", legend.text = element_text(size = rel(1.1)), 
            legend.title = element_text(size = rel(1.1))) + 
        guides(shape = F)
    ggsave(paste0(figFold, "/SSF_frust.png"), width = 5.5, height = 5.2)
    pList <- c(pList, list(p4))
    figFold <- paste0(figures, "/Figure3")
    png(paste0(figFold, "/", net, ".png"), 
        width = 17, height = 10.5, units = "in", res = 300)
    D <- plot_grid(pList[[6]] + 
        theme(axis.text.x = element_blank(), 
        legend.position = "top", legend.direction = "horizontal", 
        legend.key.width = unit(1, "cm"), legend.key.height = unit(0.2, "cm")), 
        pList[[3]] + 
        theme(legend.position = "none"), ncol = 1, rel_heights = c(2.5, 3))
    print({plot_grid(
            pList[[4]], pList[[5]], pList[[7]], pList[[1]], pList[[2]], D,
        rel_widths = rep(1/3, 6), 
        ncol = 3, 
        labels = c("A (i)", "B (i)", "C", "(ii)", "(ii)", "D"), 
        label_size = 20)})
    dev.off()

}

Figure5 <- function(net) {
    figureFolder <- paste0(figures, "/Figure5/", net)
    if (!dir.exists(figureFolder)) dir.create(figureFolder, recursive = T)
    levelVec <- c(0, 2:10)
    df <- lapply(levelVec, function(l) {
        fl <- paste0(multiLevelSim, "/", net, "/", 
            net, "_shubham_", l, "_finFlagFreq_format.csv")
        if (!file.exists(fl)) {
            print(paste0("File does not exist for level ", l, ". Simulate!"))
            return(NA)
        }
        d <- read_csv(fl, show_col_types = F) %>%
            filter(flag == 1 | Avg0 > 0.05*max(Avg0))
        d <- d %>% mutate(nLevels = ifelse(l == 0, 2, 2*l))
        hc <- hclust(d %>% select(eScore, mScore) %>% dist, method = "ward.D2")
        d <- d %>% mutate(cluster = cutree(hc, 3))
        dClustPhen <- d %>% group_by(cluster) %>%
            summarise(emScore = mean(emScore), mScore = mean(mScore), 
            eScore = mean(eScore), .groups = "drop") %>%
            mutate(PhenotypeClust = case_when(
                eScore == max(eScore) | mScore == max(mScore) ~ "Terminal",
                .default = "Hybrid"
            ))
        phen <- dClustPhen %>% pull(PhenotypeClust)
        names(phen) <- dClustPhen$cluster
        d <- d %>% mutate(PhenotypeClust = phen[cluster]) %>% select(-cluster)
        d
    }) %>% bind_rows() %>%
        mutate(nLevels = factor(nLevels, levels = 2*(1:10)))
    p4 <- segmentPlots(df)
    ggsave(paste0(figureFolder, "/segment.png"), 
        width = 5.5, height = 5.4)
    p3 <- segmentPlots(df, relative = T)
    ggsave(paste0(figureFolder, "/segment_relative.png"), 
        width = 5.5, height = 5.4)
    nSS <- df %>% group_by(nLevels) %>% summarise(nStates = n(), .groups = "drop")
    p1 <- ggplot(nSS, aes(x = nLevels, y = nStates)) +
        geom_point() + geom_line(group = 1) +
        theme_Publication() +
        scale_y_log10() +
        labs(x = "Number of levels", y = "Number of states")
    ggsave(paste0(figureFolder, "/nStates.png"),
        width = 5.5, height = 5.2)
    
    df10 <- df %>% filter(nLevels == 20) %>%
        arrange(desc(Avg0)) %>%
        mutate(CumSSF = cumsum(Avg0)) %>%
        mutate(stateID = 1:nrow(.))
    p2 <- ggplot(df10, aes(x = emScore, y = log10(Avg0))) +
        geom_point() +
        # geom_density2d() +
        labs(x = "Phenotype score", y = "Log10(SSF)") +
        theme_Publication()
    ggsave(paste0(figureFolder, "/l_20_emScore_SSF.png"),
        width = 5.5, height = 5.2)
    
    p5 <- ggplot(df10, aes(x = stateID, y = CumSSF)) +
        geom_point() +
        geom_line() +
        scale_x_log10() +
        labs(x = "Number of Steady States", y = "Cumulative SSF") +
        theme_Publication()
    ggsave(paste0(figureFolder, "/l_20_emScore_cumSSF.png"),
        width = 5.5, height = 5.1)

    ## cluster by scores
    p6 <- ggplot(df10, aes(x = PhenotypeClust, y = log10(Avg0))) +
        geom_boxplot() +
        theme_Publication() +
        labs(x = "Cluster", y = "Log10(SSF)")
    ggsave(paste0(figureFolder, "/l_20_cluster_SSF.png"),
        width = 5.5, height = 5.2)
    
    dfAll <- df %>%
        group_by(nLevels, PhenotypeClust) %>%
        summarise(nStates = n(), Freq = sum(Avg0), .groups = "drop") %>%
        mutate(PhenotypeClust = factor(PhenotypeClust, 
            levels = c("Terminal", "Hybrid")))
    ggplot(dfAll, aes(x = nLevels, y = Freq, fill = PhenotypeClust)) +
        geom_bar(stat = "identity", 
            position = position_dodge(preserve = "single")) +
        theme_Publication() + 
        theme(legend.position = "top", legend.text = element_text(size = rel(1.1)), 
            legend.title = element_text(size = rel(1.1))) +
        labs(x = "Number of Levels", y = "Total SSF", fill = "Phenotype")
    
    ggsave(paste0(figureFolder, "/Freq_levels.png"), 
        width = 5.5, height = 5.4)
    ### cowplot compilation
    figureFolder <- paste0(figures, "/Figure5")
    png(paste0(figureFolder, "/", net,".png"), width = 17, height = 10, 
        units = "in", res = 300)
    print({plot_grid(p1, p2, p3, p4, p5, p6, 
        rel_widths = rep(1/3, 6), ncol = 3, 
        labels = c("A", "C (i)", "D", "B", "(ii)", "E"), 
        label_size = 20)})
    dev.off()
}

### sequential turn off plot



Figure4 <- function(net, nLevels = c(1:2)) {
    figFold <- paste0(figures, "/Figure4/", net)
    if (!dir.exists(figFold)) dir.create(figFold, recursive = T)
    costMetrics <- c("hybridness", "meanFrust", "terminalness")
    pList <- list()
    dfList <- list()
    # browser()
    d <- paste0(multiLevelSim, "/", net, "/")
    # no turnoff
    sapply(nLevels, function(nLev) {
        l <- paste0(d, net, "_shubham_", nLev, "_finFlagFreq_format.csv")
        if (!file.exists(l)) {
            simulateNetworks(net, nLev, largeIC = F, juliaThreading = T, numThreads = numThreads)
            formatter(paste0(net, ".topo"), nLev)
        }
        df <- read_csv(l, show_col_types = F) %>%
            filter(flag == 1 | Avg0 > 0.05*max(Avg0))
        teams <- readLines(paste0(d, net, ".teams")) %>% str_split(",") %>% unlist
        nodes <- readLines(paste0(d, net, "_nodes.txt"))
        nodeOrder <- c(teams, nodes[!nodes %in% teams])
        p <- stateExpressionHeatmap(df, nodeOrder, 
            paste0(figFold, "/l_", nLev, "_no_turnOff.png"), freqWidth = T)
        if (nLev == 1) {
            pList <<- c(pList, list(p + theme(axis.text.x = element_blank(), axis.text.y = element_blank())+ labs(x = "", y = "")))
        }
        if (nLev == 2) {
            dfList <<- c(dfList, list(df %>% mutate(TurnOff = "None")))
        }
    })

    # signal turn off
    # browser()
    d <- paste0(signalFold, "/", net, "/")
    pwd <- getwd()
    setwd(d)
    fl <- list.files(".", "finFlagFreq_format.csv")
    sapply(nLevels, function(nLev) {
        l <- fl[str_detect(fl, paste0("_shubham_", nLev, "_"))]
        if (length(l) == 0) {
            print(paste0("Signal file does not exist for level ", nLev, ". Simulate!"))
            return(NA)
        }
        df <- read_csv(l, show_col_types = F) %>%
            filter(flag == 1 | Avg0 > 0.05*max(Avg0)) %>%
            filter(Avg0 > 0.01*max(Avg0))
        teams <- readLines(paste0(d, net, ".teams")) %>% str_split(",") %>% unlist
        nodes <- readLines(paste0(d, net, "_nodes.txt"))
        nodeOrder <- c(teams, nodes[!nodes %in% teams])
        p <- stateExpressionHeatmap(df, nodeOrder, 
            paste0(figFold, "/l_", nLev, "_signal_turnOff.png"), freqWidth = T)
        if (nLev == 1) {
            pList <<- c(pList, list(p + theme(axis.text.x = element_blank(), axis.text.y = element_blank(), legend.position = "none", plot.margin = unit(c(1,0,0,0), "cm")) + 
                    labs(x = "", y = "", title = "Signal node turn-off 2-level")))
        }
        if (nLev == 2) {
            dfList <<- c(dfList, list(df %>% mutate(TurnOff = "Signals")))
        }
    })

    ## all turn off
    d <- paste0(allNodeFold, "/", net, "/")
    sapply(nLevels, function(nLev) {
        df <- read_csv(paste0(d, net, "_shubham_", nLev, 
            "_turnOff_All_finFlagFreq_format.csv")) %>%
            filter(flag == 1 | Avg0 > 0.05*max(Avg0))
        teams <- readLines(paste0(d, net, ".teams")) %>% str_split(",") %>% unlist
        nodes <- readLines(paste0(d, net, "_nodes.txt"))
        nodeOrder <- c(teams, nodes[!nodes %in% teams])
        p <- stateExpressionHeatmap(df, nodeOrder, 
            paste0(figFold, "/l_", nLev, "_all_turnOff.png"), freqWidth = T)
        if (nLev == 1) 
            pList <<- c(pList, list(p + theme(axis.text.x = element_blank(), axis.text.y = element_blank(), legend.position = "none", plot.margin = unit(c(0,0,0,0), "cm")) + 
            labs(x = "", y = "", title = "All node turn-off 2-level")))
        if (nLev == 2) {
            pList <<- c(pList, list(p + theme(axis.text.y = element_blank(), axis.text.x = element_text(size = rel(0.7)), legend.position = "none", plot.margin = unit(c(0,0,0,0), "cm")) + 
            labs(x = "", y = "", title = "All node turn-off 4-level")))
            dfList <<- c(dfList, list(df %>% mutate(TurnOff = "All")))
        }
    })
    df <- bind_rows(dfList) %>%
        mutate(TurnOff = factor(TurnOff, levels = c("None", "Signals", "All"))) %>%
        mutate(Phenotype = str_replace(Phenotype, " ", "\n"))
    p <- ggplot(df, aes(x = Phenotype, y = log10(Avg0), fill = TurnOff, color = TurnOff)) +
        geom_boxplot(alpha = 0.5) +
        geom_point(position = position_jitterdodge(), alpha = 0.65) + 
        theme_Publication() +
        theme(legend.position = c(0.6, 0.1), axis.text = element_text(size = rel(0.8)), axis.title.y = element_text(size = rel(0.7)), 
            legend.text = element_text(size = rel(0.8)), legend.title = element_text(size = rel(0.8)), plot.margin = unit(c(0,0,0,0), "cm")) +
        labs(x = "", y = "Log10(SSF)", color = "Turn off") +
        guides(fill = F)
    ggsave(paste0(figFold, "/turnOff_SSF.png"),
        width = 5.5, height = 5.2)
    pList <- c(pList, list(p))

    sapply(costMetrics, function(costMetric) {
        pL <- seqTurnPlot(net, costMetric, nLevels = nLevels)
        ggsave(paste0(figFold, "/", costMetric, "_", nLevels[1], ".png"), plot = pL[[1]], 
            width = 12, height = 5)
        ggsave(paste0(figFold, "/", costMetric, "_", nLevels[2], ".png"), plot = pL[[2]],
            width = 12, height = 5)
        if (costMetric == "meanFrust") {
            pList <<- c(pList, list(pL[[1]] + theme(axis.title = element_text(size = rel(1.1)), axis.text = element_text(size = rel(0.85)))))
            pList <<- c(pList, list(pL[[2]] + theme(axis.title = element_text(size = rel(1.1)), axis.text = element_text(size = rel(0.85)))))
        }
    })
    
    png(paste0(figures, "/Figure4/", net, ".png"), width = 17, height = 9, units = "in", res = 300)
    print({plot_grid(pList[[1]], 
        plot_grid(pList[[2]], pList[[3]], pList[[4]], pList[[5]], 
            rel_heights = c(1, 0.6, 1, 1), ncol = 1, labels = c("B", "C", "D", "E"), 
            label_size = 20),
            plot_grid(pList[[6]], pList[[7]], rel_heights = c(1, 1), ncol = 1, labels = c("F", "G"), label_size = 20), 
            rel_widths = c(1,1,1.5), ncol = 3, labels = c("A", "", ""), label_size = 20)})
    dev.off()
    # pList <- seqTurnPlot(net, "meanFrust", nLevels = c(1:2))

    
}