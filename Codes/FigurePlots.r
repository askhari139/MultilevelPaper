library(compiler)



compareHybridSS <- function(net, levelVec = c(1,2), figureDir = ".") {
    wd <- getwd()
    key <- paste0(net, "_shubham")
    DirectoryNav(paste0(multiLevelSim, "/", net))
    files <- paste0(key, "_", levelVec, "_finFlagFreq_format.csv")
    # levelVec <- ifelse(levelVec == 100, "RACIPE", levelVec)
    plotKey <- paste0(levelVec, collapse = "_")
    countKey <- c("Hybrid", "Incomplete Terminal", "Terminal", "Total")
    names(countKey) <- c("nHybrid", "nIncomplete", "nTerminal", "nTotal")
    freqKey <- countKey
    names(freqKey) <- c("freqHybrid", "freqIncomplete", "freqTerminal")
    dat <- sapply(files, function(x) {
        d <- read_csv(x, show_col_types = F)
        nTotal <- nrow(d)
        nHybrid <- nrow(d %>% filter(Phenotype == "Hybrid"))
        nIncomplete <- nrow(d %>% filter(str_detect(Phenotype, "Incomplete")))
        nTerminal <- nrow(d %>% filter(Phenotype == "Terminal"))
        nSteady <- nrow(d %>% filter(flag == 1))
        freqHybrid <- d %>% filter(Phenotype == "Hybrid") %>%
            select(Avg0) %>% unlist %>% sum
        freqIncomplete <- d %>% filter(str_detect(Phenotype, "Incomplete")) %>%
            select(Avg0) %>% unlist %>% sum
        freqTerminal <- d %>% filter(Phenotype == "Terminal") %>%
            select(Avg0) %>% unlist %>% sum
        freqSteady <- d %>% filter(flag == 1) %>%
            select(Avg0) %>% unlist %>% sum
        frustHybrid <- d %>% filter(Phenotype == "Hybrid") %>% 
            select(frust0) %>% unlist %>% mean
        frustIncomplete <- d %>% filter(str_detect(Phenotype, "Incomplete")) %>%
            select(frust0) %>% unlist %>% mean
        scoreHybrid <- d %>% filter(Phenotype == "Hybrid") %>% 
            select(emScore) %>% unlist %>% abs %>% mean
        scoreIncomplete <- d %>% filter(str_detect(Phenotype, "Incomplete")) %>%
            select(emScore) %>% unlist %>% abs %>% mean
        scoreHybridSd <- d %>% filter(Phenotype == "Hybrid") %>% 
            select(emScore) %>% unlist %>% abs %>% sd
        scoreIncompleteSd <- d %>% filter(str_detect(Phenotype, "Incomplete")) %>%
            select(emScore) %>% unlist %>% abs %>% sd
        # isingFrustHybrid <- d %>% filter(Phenotype == "Hybrid") %>% 
        #     select(isingFrust) %>% unlist %>% mean
        # isingFrustIncomplete <- d %>% filter(str_detect(Phenotype, "Incomplete")) %>%
        #     select(isingFrust) %>% unlist %>% mean
        c(nTotal,nHybrid, nIncomplete, nTerminal, nSteady,freqHybrid, 
          freqIncomplete, freqTerminal, freqSteady, frustHybrid, frustIncomplete,
          scoreHybrid, scoreIncomplete, scoreHybridSd, scoreIncompleteSd)
    }) %>% t %>% as.data.frame %>%
        set_names(c("nTotal", "nHybrid", "nIncomplete", 
                    "nTerminal", "nSteady", "freqHybrid", "freqIncomplete", 
                    "freqTerminal", "freqSteady",
            "frustHybrid", "frustIncomplete", "scoreHybrid", "scoreIncomplete",
                     "scoreHybridSd", "scoreIncompleteSd")) %>%
        mutate(levels = factor(levelVec*2, levels = levelVec*2))
    datPhenCount1 <- dat %>% select(levels, nHybrid, nTotal, 
                                    nIncomplete, nTerminal) %>%
        gather(key = "Phenotype", value = "Count", -levels) %>%
        mutate(Phenotype = factor(countKey[Phenotype], levels = countKey))
    datPhenFreq1 <- dat %>% select(levels, freqHybrid, freqIncomplete, 
                                   freqTerminal) %>%
        gather(key = "Phenotype", value = "Frequency", -levels) %>%
        mutate(Phenotype = factor(freqKey[Phenotype], 
            levels = freqKey))
    
    datViol <- lapply(files, function(x) {
        nLevels <- x %>% str_remove(paste0(key, "_")) %>%
            str_remove("_finFlagFreq_format.csv") %>%
            as.character
        # nLevels <- ifelse(nLevels == "100", "RACIPE", nLevels)
        d <- read_csv(x, show_col_types = F)
        d %>% select(Phenotype, Avg0, frust0, emScore, 
                     eScore, mScore) %>%
            mutate(nLevels = nLevels)
    }) %>% bind_rows %>%
        mutate(nLevels = factor(nLevels, levels = levelVec*2))
    # DirectoryNav("Figures")
    setwd(figureDir)
    DirectoryNav("multiCompare")
    ggplot(datPhenCount1, aes(x = levels, y = Count, fill = Phenotype)) +
        geom_bar(stat = "identity", position = position_dodge(preserve = "single")) +
        labs(x = "Number of levels", y = "Number of states") +
        theme_Publication() +
        theme(legend.position = c(0.2, 0.8),
            legend.direction = "vertical")
    ggsave(paste0(net,"_", plotKey, "_PhenCount.png"), width = 5.5, height = 5.2)

    ggplot(datPhenCount1, aes(x = levels, y = Count, fill = Phenotype)) +
        geom_bar(stat = "identity", position = position_dodge(preserve = "single")) +
        labs(x = "Number of levels", y = "Number of states") +
        theme_Publication() +
        scale_y_log10() +
        theme(legend.position = c(0.2, 0.8),
            legend.direction = "vertical")
    ggsave(paste0(net,"_", plotKey, "_PhenCountLog.png"), width = 5.5, height = 5.2)

    ggplot(datPhenFreq1, aes(x = levels, y = Frequency, fill = Phenotype)) +
        geom_bar(stat = "identity", position = position_dodge(preserve = "single")) +
        labs(x = "Number of levels", y = "Log10 (SSF)") +
        theme_Publication() +
        theme(legend.position = "top")
    ggsave(paste0(net,"_", plotKey, "_PhenFreqAbs.png"), width = 5.5, height = 5.5)
    # setwd(wd)
### voilins Data ----

    
    ggplot(datViol, aes(x = nLevels, y = Avg0, fill = Phenotype)) +
        geom_violin(width = 1) +
        scale_y_log10() +
        labs(x = "Number of levels", y = "Log10 (SSF)", fill = "Phenotype") +
        theme_Publication() +
        theme(legend.position = c(0.8, 0.8),
            legend.direction = "vertical")
    ggsave(paste0(net,"_", plotKey, "_FreqDist.png"), width = 5.5, height = 5.2)
    ggplot(datViol, aes(x = nLevels, y = Avg0, fill = Phenotype)) +
        geom_boxplot() +
        scale_y_log10() +
        labs(x = "Number of levels", y = "Log10 (SSF)", fill = "Phenotype") +
        theme_Publication() +
        theme(legend.position = c(0.8, 0.8),
            legend.direction = "vertical")
    ggsave(paste0(net,"_", plotKey, "_FreqDistBox.png"), width = 5.5, height = 5.2)

    ggplot(datViol, aes(x = nLevels, y = frust0, fill = Phenotype)) +
        geom_violin(width = 1) +
        labs(x = "Number of levels", y = "Frustration", fill = "Phenotype") +
        theme_Publication() +
        theme(legend.position = c(0.8, 0.8),
            legend.direction = "vertical")
    ggsave(paste0(net,"_", plotKey, "_FrustDistAbs.png"), width = 5.5, height = 5.2)
    ggplot(datViol, aes(x = nLevels, y = frust0, fill = Phenotype)) +
        geom_boxplot() +
        labs(x = "Number of levels", y = "Frustration", fill = "Phenotype") +
        theme_Publication() +
        theme(legend.position = c(0.8, 0.8),
            legend.direction = "vertical")
    ggsave(paste0(net,"_", plotKey, "_FrustDistBox.png"), width = 5.5, height = 5.2)

    ggplot(datViol, aes(x = nLevels, y = abs(emScore), fill = Phenotype)) +
        geom_violin() +
        labs(x = "Number of levels", y = "Absolute value of emScore", fill = "Phenotype") +
        theme_Publication() +
        theme(legend.position = "top")
    ggsave(paste0(net,"_", plotKey, "_emScoreDist.png"), width = 5.5, height = 5.2)

    ggplot(datViol, aes(x = nLevels, y = abs(emScore), fill = Phenotype)) +
        geom_boxplot() +
        labs(x = "Number of levels", y = "Absolute value of emScore", fill = "Phenotype") +
        theme_Publication() +
        theme(legend.position = "top")
    ggsave(paste0(net,"_", plotKey, "_emScoreDistBox.png"), width = 5.5, height = 5.2)
    setwd("..")
    sapply(levelVec, function(l) {
        DirectoryNav(paste0("nLevels_", l))
        d <- datViol %>% filter(nLevels == l)
        ggplot(d, aes(x = eScore, y = mScore, color = frust0)) +
            geom_point() +
            scale_color_viridis_c() +
            labs(x = "eScore", y = "mScore", color = "Frustration") +
            theme_Publication() +
            theme(legend.position = "top",
            legend.key.width = unit(0.6, "cm"))
        ggsave(paste0(net,"_", plotKey, "_emScoreScatter_", l, ".png"), 
            width = 6, height = 5.5)
        setwd("..")
    })
### -------
    setwd(wd)
}

scoreSegregation <- function(net) {
    wd <- getwd()
    files <- list.files(pattern = paste0(net, "_shubham_", ".*", "_finFlagFreq_format.csv"))
    files <- c(files, paste0(net, "_finFlagFreq_format.csv"))
    files <- paste0(wd, "/", files)
    nodes <- readLines(paste0(net, "_nodes.txt"))
    teams <- readLines(paste0(net, ".teams")) %>%
        str_split(",")
    nodeOrder <- c(unlist(teams), nodes[!nodes %in% unlist(teams)] %>% sort %>% unique)
    scoringMethods <- c("stateSum")
    DirectoryNav(paste0("Figures/", net))
    # sapply(scoringMethods, function(m) {
        # DirectoryNav(m)
    sapply(files, function(x) {
        fName <- x %>% str_remove(".*/")
        nLevels <- fName %>% str_remove(paste0(net, "_shubham_")) %>%
            str_remove("_finFlagFreq_format.csv") %>%
            as.integer
        if (is.na(nLevels)) {
            nLevels <- 1
        }

        d <- read_csv(x, show_col_types = F)
        d1 <- d %>% select(all_of(nodeOrder))
        d <- d %>%
                mutate(eScore = d1 %>% select(all_of(teams[[1]])) %>%
                            sapply(as.numeric) %>% rowSums,
                        mScore = d1 %>% select(all_of(teams[[2]])) %>%
                            sapply(as.numeric) %>% rowSums) %>%
                mutate(eScore = eScore/length(teams[[1]]),
                        mScore = mScore/length(teams[[2]]),
                        emScore = eScore - mScore)
        d <- d %>%
            mutate(PhenAbs = ifelse(abs(emScore) == 2, "Terminal", "Hybrid"),
                PhenAbs = ifelse((eScore == -1 | mScore == -1) & abs(emScore) < 2, 
                        paste0("Incomplete ", PhenAbs), PhenAbs),
                PhenAbs = ifelse(PhenAbs == "Incomplete Hybrid", 
                    "Incomplete Terminal", PhenAbs)) %>%
            mutate(PhenThresh = ifelse(abs(emScore) <= 1, "Hybrid", "Terminal"),
                PhenThresh = ifelse((eScore == -1 | mScore == -1) & abs(emScore) < 1, 
                    paste0("Incomplete ", PhenThresh), PhenThresh),
                PhenThresh = ifelse(PhenThresh == "Incomplete Hybrid",
                    "Incomplete Terminal", PhenThresh))
        
        write_csv(d, fName, quote = "none")
    })
        # setwd("..")
    # })
    setwd(wd)
}

stateExpressionHeatmap <- function(dat, nodeOrder, key, maxStates = 100) {
    # browser()
    d <- dat %>%
            arrange(emScore, desc(Avg0)) %>%
            mutate(stateNum = 1:nrow(.))
    d <- d %>% #filter(stateNum <= maxStates) %>%
            select(all_of(nodeOrder), stateNum, Avg0, emScore) %>%
            gather(key = "Node", value = "State", -Avg0, -emScore, -stateNum) %>%
            mutate(Node = factor(Node, levels = nodeOrder)) %>%
            arrange(desc(stateNum))
    d1 <- d %>% filter(stateNum <= maxStates)
    ggplot(d1, aes(x = Node, y = stateNum, fill = State)) +
        geom_tile() +
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

    p <- ggplot(d, aes(x = Node, y = stateNum, fill = State)) +
        geom_tile() +
        theme_Publication() +
        labs(x = "", y = "", fill = "State") +
        scale_fill_gradient2(low = "blue", mid = "white", 
            high = "red", midpoint = 0) +
        scale_x_discrete(expand = c(0,0)) +
        scale_y_continuous(expand = c(0,0)) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
                legend.position = "right", legend.direction = "vertical",
                legend.key.height = unit(0.8, "cm"))
    ggsave(key %>% str_replace(".png", "_fixed.png"), 
        width = 10, height = 6.5)
}

frustScorePlots <- function(dat, key) {
    if (nrow(dat) < 5) {
        print("Not enough states to plot")
        return(0)
    }
    #######Phen score vs frustration-------
    ggplot(dat, aes(x = emScore, y = frust0)) +
        geom_point() +
        geom_density2d() +
        labs(x = "Phenotype score", y = "Frustration") +
        theme_Publication()
    ggsave(paste0(key, "_frustScore.png"), width = 5.5, height = 5)
    ggplot(dat, aes(x = emScore, y = log10(Avg0))) +
        geom_point() +
        geom_density2d() +
        labs(x = "Phenotype score", y = "Log10 (SSF)") +
        theme_Publication()
    ggsave(paste0(key, "_freqScore.png"), width = 5.5, height = 5)

    ggplot(dat, aes(x = emScore, y = frust0, color = PhenAbs)) +
        geom_point() +
        # geom_density2d() +
        labs(x = "Phenotype score", y = "Frustration", color = "Phenotype") +
        theme_Publication() +
        theme(legend.position = "top", legend.direction = "horizontal")
    ggsave(paste0(key, "_frustScoreAbsPhen.png"), width = 5.5, height = 5)
    ggplot(dat, aes(x = emScore, y = log10(Avg0), color = PhenAbs)) +
        geom_point() +
        # geom_density2d() +
        labs(x = "Phenotype score", y = "Log10 (SSF)", color = "Phenotype") +
        theme_Publication() +
        theme(legend.position = "top", legend.direction = "horizontal")
    ggsave(paste0(key, "_freqScoreAbsPhen.png"), width = 5.5, height = 5)

    ggplot(dat, aes(x = emScore, y = frust0, color = PhenThresh)) +
        geom_point() +
        # geom_density2d() +
        labs(x = "Phenotype score", y = "Frustration", color = "Phenotype") +
        theme_Publication() +
        theme(legend.position = "top", legend.direction = "horizontal")
    ggsave(paste0(key, "_frustScoreThreshPhen.png"), width = 5.5, height = 5)
    ggplot(dat, aes(x = emScore, y = log10(Avg0), color = PhenThresh)) +
        geom_point() +
        # geom_density2d() +
        labs(x = "Phenotype score", y = "Log10 (SSF)", color = "Phenotype") +
        theme_Publication() +
        theme(legend.position = "top", legend.direction = "horizontal")
    ggsave(paste0(key, "_freqScoreThreshPhen.png"), width = 5.5, height = 5)

    # frequency vs frustration for hybrid states, threshold and abs
    ggplot(dat %>% filter(PhenAbs == "Hybrid"), 
        aes(x = log10(Avg0), y = frust0)) +
        geom_point() +
        geom_density2d() +
        labs(x = "Log10 (SSF)", y = "Frustration") +
        theme_Publication()
    ggsave(paste0(key, "_frustScoreAbsHybrid.png"), width = 5.5, height = 5)

    ggplot(dat %>% filter(PhenThresh == "Hybrid"), 
        aes(x = log10(Avg0), y = frust0)) +
        geom_point() +
        geom_density2d() +
        labs(x = "Log10 (SSF)", y = "Frustration") +
        theme_Publication()
    ggsave(paste0(key, "_frustScoreThreshHybrid.png"), width = 5.5, height = 5)



    #### eScore vs mScore ----
    ggplot(dat, aes(x = eScore, y = mScore)) +
        geom_point() +
        geom_density2d() +
        labs(x = "E score", y = "M score") +
        theme_Publication()
    ggsave(paste0(key, "_eScore_mScore.png"), width = 5.5, height = 5)

    ggplot(dat, aes(x = eScore, y = mScore, color = PhenAbs)) +
        geom_point() +
        # geom_density2d() +
        labs(x = "E score", y = "M score", color = "Phenotype") +
        theme_Publication() +
        theme(legend.position = "top", legend.direction = "horizontal")
    ggsave(paste0(key, "_eScore_mScoreAbsPhen.png"), width = 5.5, height = 5)

    ggplot(dat, aes(x = eScore, y = mScore, color = PhenThresh)) +
        geom_point() +
        # geom_density2d() +
        labs(x = "E score", y = "M score", color = "Phenotype") +
        theme_Publication() +
        theme(legend.position = "top", legend.direction = "horizontal")
    ggsave(paste0(key, "_eScore_mScoreThreshPhen.png"), width = 5.5, height = 5)

    ggplot(dat, aes(x = eScore, y = mScore, color = frust0)) +
        geom_point() +
        # geom_density2d() +
        labs(x = "E score", y = "M score", color = "Frustration") +
        theme_Publication() +
        theme(legend.position = "top", legend.direction = "horizontal",
            legend.key.width = unit(0.6, "cm"))
    ggsave(paste0(key, "_eScore_mScoreThreshFrust.png"), width = 5.5, height = 5)

    #### eScore vs mScore Hybrid phenotype----
    ggplot(dat %>% filter(PhenAbs == "Hybrid"), 
        aes(x = eScore, y = mScore)) +
        geom_point() +
        geom_density2d() +
        labs(x = "E score", y = "M score") +
        theme_Publication()
    ggsave(paste0(key, "_eScore_mScoreAbsHybrid.png"), width = 5.5, height = 5)

    ggplot(dat %>% filter(PhenAbs == "Hybrid"), 
        aes(x = eScore, y = mScore, color = frust0)) +
        geom_point() +
        # geom_density2d() +
        labs(x = "E score", y = "M score", color = "Frustration") +
        theme_Publication() +
        theme(legend.position = "top", legend.direction = "horizontal",
            legend.key.width = unit(0.6, "cm"))
    ggsave(paste0(key, "_eScore_mScoreAbsHybridFrust.png"), width = 5.5, height = 5)

    ggplot(dat %>% filter(PhenThresh == "Hybrid"), 
        aes(x = eScore, y = mScore)) +
        geom_point() +
        geom_density2d() +
        labs(x = "E score", y = "M score") +
        theme_Publication()
    ggsave(paste0(key, "_eScore_mScoreThreshHybrid.png"), width = 5.5, height = 5)

    ggplot(dat %>% filter(PhenThresh == "Hybrid"), 
        aes(x = eScore, y = mScore, color = frust0)) +
        geom_point() +
        # geom_density2d() +
        labs(x = "E score", y = "M score", color = "Frustration") +
        theme_Publication() +
        theme(legend.position = "top", legend.direction = "horizontal",
            legend.key.width = unit(0.6, "cm"))
    ggsave(paste0(key, "_eScore_mScoreThreshHybridFrust.png"), width = 5.5, height = 5)
}

isingMultiCompare <- function(net, nLevels = 2, compute = F) {
    fl <- paste0(net, "_compare_ising_", nLevels, ".csv")
    if (file.exists(fl) && !compute) return(0)
    ising <- paste0(net, "_finFlagFreq_format.csv")
    shubham <- paste0(net, "_shubham_",nLevels,"_finFlagFreq_format.csv")
    ising <- ising %>% read_csv(show_col_types = F)
    shubham <- shubham %>% read_csv(show_col_types = F)
    isingStates <- ising %>% select(states, Avg0, PhenThresh, PhenAbs, frust0) %>%
        mutate(isingFrust = frust0) %>% select(-frust0)
    shubhamStates <- shubham %>% select(states, Avg0, PhenThresh, PhenAbs, frust0, 
                                        eScore, mScore) %>%
        mutate(statesMult = states, shubhamFreq = Avg0) %>%
        mutate(states = states %>% sapply(isingConvert, numLevels = nLevels)) %>%
        mutate(shubhamPhenThresh = PhenThresh, 
            shubhamPhenAbs = PhenAbs) %>% select(-PhenThresh, -PhenAbs, -Avg0)
    df <- merge(isingStates, shubhamStates, by = "states", all = T) %>%
        mutate(Label = ifelse(!is.na(Avg0), "isingState", "newState"),
               PhenThresh = ifelse(is.na(Avg0), "Absent", PhenThresh),
               PhenAbs = ifelse(is.na(Avg0), "Absent", PhenAbs)) %>%
               mutate(ConversionAbs = paste0(PhenAbs, "_To_", shubhamPhenAbs),
               ConversionThresh = paste0(PhenThresh, "_To_", shubhamPhenThresh))
    write_csv(df, paste0(net, "_compare_ising_", nLevels, ".csv"), quote = "none")
    return(0)
}

isingMultiCompare <- cmpfun(isingMultiCompare)

isingComparePlots <- function(net, levelVec = c(1,2,100), hybridTresh = 0.01) {
    dat <- lapply(levelVec, function(l) {
        isingMultiCompare(net, l, compute = T)
        
        d <- read_csv(paste0(net, "_compare_ising_", l, ".csv"), show_col_types = F)
        DirectoryNav(paste0("nLevels_", l))
    #### Individual comparison plots ----------------------------------------------
        # top 5 and bottom 5 hybrid states by frequency both threshold and abs
        if (nrow(d %>% filter(PhenAbs == "Hybrid")) && 
                nrow(d %>% filter(shubhamPhenAbs == "Hybrid"))) {
            d1Ising <- d %>% filter(PhenAbs == "Hybrid") %>%
                arrange(Avg0) %>%
                select(Avg0, isingFrust) %>%
                set_names(c("Freq", "Frust")) %>%
                headNtail(5) %>%
                mutate(Formalism = "Ising")
            d1Multi <- d %>% filter(shubhamPhenAbs == "Hybrid") %>%
                arrange(shubhamFreq) %>%
                select(shubhamFreq, frust0) %>%
                set_names(c("Freq", "Frust")) %>%
                headNtail(5) %>%
                mutate(Formalism = "MultiLevel")
            d1 <- rbind.data.frame(d1Ising, d1Multi)
            ggplot(d1, aes(x = log10(Freq), y = Frust, color = Formalism)) +
                geom_jitter() +
                labs(x = "Log10 (SSF)", y = "Frustration", 
                    color = "Formalism") +
                theme_Publication() +
                theme(legend.position = "top", legend.direction = "horizontal",
                    legend.key.width = unit(0.6, "cm"))
            ggsave(paste0(net, "_frustFreqAbsHybridFreqExtreme_", l, ".png"), 
                width = 5.5, height = 5)

            # total frequency of states across formalisms

            d1Ising <- d %>%
                select(states, PhenAbs, Avg0) %>%
                filter(PhenAbs %in% c("Hybrid", "Terminal", "Incomplete Terminal")) %>%
                unique %>% select(-states) %>%
                set_names(c("Phenotype", "Freq")) %>%
                group_by(Phenotype) %>%
                summarise(Freq = sum(Freq), nStates = n(), .groups = "drop") %>%
                mutate(Formalism = "Ising")
            d1Multi <- d %>%
                select(shubhamPhenAbs, shubhamFreq, statesMult) %>%
                filter(shubhamPhenAbs %in% c("Hybrid", "Terminal", "Incomplete Terminal")) %>%
                unique %>% select(-statesMult) %>%
                set_names(c("Phenotype", "Freq")) %>%
                group_by(Phenotype) %>%
                summarise(Freq = sum(Freq), nStates = n(), .groups = "drop") %>%
                mutate(Formalism = "MultiLevel")
            d1 <- rbind.data.frame(d1Ising, d1Multi)
            ggplot(d1, aes(x = Phenotype, y = Freq, fill = Formalism)) +
                geom_bar(stat = "identity", position = position_dodge(preserve = "single")) +
                labs(x = "Phenotype", y = "Frequency", fill = "Formalism") +
                theme_Publication() +
                theme(legend.position = "top", legend.direction = "horizontal",
                    legend.key.width = unit(0.6, "cm"))
            ggsave(paste0(net, "_totalFreqAbs_", l, ".png"), width = 5.5, height = 5)
            ggplot(d1, aes(x = Phenotype, y = nStates, fill = Formalism)) +
                geom_bar(stat = "identity", position = position_dodge(preserve = "single")) +
                scale_y_log10() +
                labs(x = "Phenotype", y = "# Steady states", fill = "Formalism") +
                theme_Publication() +
                theme(legend.position = "top", legend.direction = "horizontal",
                    legend.key.width = unit(0.6, "cm"))
            ggsave(paste0(net, "_totalCountAbs_", l, ".png"), width = 5.5, height = 5)
        }

        if (nrow(d %>% filter(PhenThresh == "Hybrid")) && 
                nrow(d %>% filter(shubhamPhenThresh == "Hybrid"))) {
            d1Ising <- d %>% filter(PhenThresh == "Hybrid") %>%
                arrange(Avg0) %>%
                select(Avg0, isingFrust) %>%
                set_names(c("Freq", "Frust")) %>%
                headNtail(5) %>%
                mutate(Formalism = "Ising")
            d1Multi <- d %>% filter(shubhamPhenThresh == "Hybrid") %>%
                arrange(shubhamFreq) %>%
                select(shubhamFreq, frust0) %>%
                set_names(c("Freq", "Frust")) %>%
                headNtail(5) %>%
                mutate(Formalism = "MultiLevel")
            d1 <- rbind.data.frame(d1Ising, d1Multi)
            ggplot(d1, aes(x = log10(Freq), y = Frust, color = Formalism)) +
                geom_jitter() +
                labs(x = "Log10 (SSF)", y = "Frustration", 
                    color = "Formalism") +
                theme_Publication() +
                theme(legend.position = "top", legend.direction = "horizontal",
                    legend.key.width = unit(0.6, "cm"))
            ggsave(paste0(net, "_frustFreqThreshHybridFreqExtreme_", l, ".png"), 
                width = 5.5, height = 5)
            d1Ising <- d %>%
                select(states, PhenAbs, Avg0) %>%
                filter(PhenAbs %in% c("Hybrid", "Terminal", "Incomplete Terminal")) %>%
                unique %>% select(-states) %>%
                set_names(c("Phenotype", "Freq")) %>%
                group_by(Phenotype) %>%
                summarise(Freq = sum(Freq), nStates = n(), .groups = "drop") %>%
                mutate(Formalism = "Ising")
            d1Multi <- d %>%
                select(shubhamPhenAbs, shubhamFreq, statesMult) %>%
                filter(shubhamPhenAbs %in% c("Hybrid", "Terminal", "Incomplete Terminal")) %>%
                unique %>% select(-statesMult) %>%
                set_names(c("Phenotype", "Freq")) %>%
                group_by(Phenotype) %>%
                summarise(Freq = sum(Freq), nStates = n(), .groups = "drop") %>%
                mutate(Formalism = "MultiLevel")
            d1 <- rbind.data.frame(d1Ising, d1Multi)
            ggplot(d1, aes(x = Phenotype, y = Freq, fill = Formalism)) +
                geom_bar(stat = "identity", position = position_dodge(preserve = "single")) +
                labs(x = "Phenotype", y = "Frequency", fill = "Formalism") +
                theme_Publication() +
                theme(legend.position = "top", legend.direction = "horizontal",
                    legend.key.width = unit(0.6, "cm"))
            ggsave(paste0(net, "_totalFreqThresh_", l, ".png"), width = 5.5, height = 5)
            ggplot(d1, aes(x = Phenotype, y = nStates, fill = Formalism)) +
                geom_bar(stat = "identity", position = position_dodge(preserve = "single")) +
                scale_y_log10() +
                labs(x = "Phenotype", y = "# Steady states", fill = "Formalism") +
                theme_Publication() +
                theme(legend.position = "top", legend.direction = "horizontal",
                    legend.key.width = unit(0.6, "cm"))
            ggsave(paste0(net, "_totalCountThresh_", l, ".png"), width = 5.5, height = 5)
        }

        # conversion heatmap
        
        d1 <- d %>%
            group_by(PhenThresh, shubhamPhenThresh) %>%
            summarise(numStates = n(), 
                shubhamFreq = sum(shubhamFreq),
                .groups = "drop")
        ggplot(d1, aes(x = PhenThresh, y = shubhamPhenThresh, fill = numStates)) +
            geom_tile() +
            scale_fill_viridis_c() +
            theme_Publication() +
            labs(x = "Ising Phenotype", y = "MultiLevel Phenotype", fill = "# States") +
            theme(legend.position = "right", legend.direction = "vertical",
                legend.key.height = unit(0.8, "cm"),
                axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
        ggsave(paste0(net, "_conversionCountThresh_", l, ".png"), width = 6, height = 5)
        ggplot(d1, aes(x = PhenThresh, y = shubhamPhenThresh, fill = shubhamFreq)) +
            geom_tile() +
            scale_fill_viridis_c() +
            theme_Publication() +
            labs(x = "Ising Phenotype", y = "MultiLevel Phenotype",
                fill = "Frequency") +
            theme(legend.position = "right", legend.direction = "vertical",
                legend.key.height = unit(0.8, "cm"),
                axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
        ggsave(paste0(net, "_conversionFreqThresh_", l, ".png"), width = 6, height = 5)
        
        d1 <- d %>%
            group_by(PhenAbs, shubhamPhenAbs) %>%
            summarise(numStates = n(), 
                shubhamFreq = sum(shubhamFreq),
                .groups = "drop")
        ggplot(d1, aes(x = PhenAbs, y = shubhamPhenAbs, fill = numStates)) +
            geom_tile() +
            scale_fill_viridis_c() +
            theme_Publication() +
            labs(x = "Ising Phenotype", y = "MultiLevel Phenotype", fill = "# States") +
            theme(legend.position = "right", legend.direction = "vertical",
                legend.key.height = unit(0.8, "cm"),
                axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
        ggsave(paste0(net, "_conversionCountAbs_", l, ".png"), width = 6, height = 5)
        ggplot(d1, aes(x = PhenAbs, y = shubhamPhenAbs, fill = shubhamFreq)) +
            geom_tile() +
            scale_fill_viridis_c() +
            theme_Publication() +
            labs(x = "Ising Phenotype", y = "MultiLevel Phenotype", 
                fill = "Frequency") +
            theme(legend.position = "right", legend.direction = "vertical",
                legend.key.height = unit(0.8, "cm"),
                axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
        ggsave(paste0(net, "_conversionFreqAbs_", l, ".png"), width = 6, height = 5)
        d <- d %>% filter(!is.na(PhenAbs))
        ggplot(d %>% mutate(ConversionAbs = str_replace_all(ConversionAbs, "_", "\n")),
            aes(x = ConversionAbs, y = shubhamFreq)) +
            geom_violin() +
            scale_y_log10() +
            theme_Publication() +
            labs(x = "Conversion", y = "Log10 (SSF)") +
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                size = rel(0.8)))
        ggsave(paste0(net, "_conversionFreqAbsViolin_", l, ".png"), width = 6, height = 5)
        ggplot(d %>% mutate(ConversionAbs = str_replace_all(ConversionAbs, "_", "\n")),
            aes(x = ConversionAbs, y = shubhamFreq)) +
            geom_boxplot() +
            scale_y_log10() +
            theme_Publication() +
            labs(x = "Conversion", y = "Log10 (SSF)") +
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                size = rel(0.8)))
        ggsave(paste0(net, "_conversionFreqAbsBox_", l, ".png"), width = 6, height = 5)

        ggplot(d %>% mutate(ConversionThresh = str_replace_all(ConversionThresh, "_", "\n")),
            aes(x = ConversionThresh, y = shubhamFreq)) +
            geom_violin() +
            theme_Publication() +
            labs(x = "Conversion", y = "Log10 (SSF)") +
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
        ggsave(paste0(net, "_conversionFreqThreshViolin_", l, ".png"), width = 6, height = 5)

        ggplot(d %>% mutate(ConversionAbs = str_replace_all(ConversionAbs, "_", "\n")),
            aes(x = ConversionAbs, y = frust0)) +
            geom_violin() +
            theme_Publication() +
            labs(x = "Conversion", y = "Frustration") +
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
        ggsave(paste0(net, "_conversionFrustAbsViolin_", l, ".png"), width = 6, height = 5)

        ggplot(d %>% mutate(ConversionThresh = str_replace_all(ConversionThresh, "_", "\n")),
            aes(x = ConversionThresh, y = frust0)) +
            geom_violin() +
            theme_Publication() +
            labs(x = "Conversion", y = "Frustration") +
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
        ggsave(paste0(net, "_conversionFrustThreshViolin_", l, ".png"), width = 6, height = 5)
    setwd("..")
    ### Data return --------------------------------------------------------------
    return(d %>% mutate(nLevels = as.character(l)))
    }) %>% bind_rows
    ### Combined comparison plots ------------------------------------------------
    levelVec <- levelVec %>% as.character
    dat <- dat %>% mutate(nLevels = factor(nLevels, levels = levelVec*2))
    DirectoryNav("multiCompare")
    ggplot(dat, aes(x = ConversionAbs, y = shubhamFreq, fill = nLevels)) +
        geom_violin() +
        theme_Publication() +
        labs(x = "Conversion", y = "Log10 (SSF)", fill = "Number of levels") +
        theme(legend.position = "top")
    ggsave(paste0(net, "_conversionFreqAbsViolinAll.png"), width = 10, height = 5)

    ggplot(dat, aes(x = ConversionThresh, y = shubhamFreq, fill = nLevels)) +
        geom_violin() +
        theme_Publication() +
        labs(x = "Conversion", y = "Log10 (SSF)", fill = "Number of levels") +
        theme(legend.position = "top")
    ggsave(paste0(net, "_conversionFreqThreshViolinAll.png"), width = 10, height = 5)

    ggplot(dat, aes(x = ConversionAbs, y = frust0, fill = nLevels)) +
        geom_violin() +
        theme_Publication() +
        labs(x = "Conversion", y = "Frustration", fill = "Number of levels") +
        theme(legend.position = "top",
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    ggsave(paste0(net, "_conversionFrustAbsViolinAll.png"), width = 10, height = 5)

    ggplot(dat, aes(x = ConversionThresh, y = frust0, fill = nLevels)) +
        geom_violin() +
        theme_Publication() +
        labs(x = "Conversion", y = "Frustration", fill = "Number of levels") +
        theme(legend.position = "top",
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    ggsave(paste0(net, "_conversionFrustThreshViolinAll.png"), width = 10, height = 5)

    ### Hybrid stability comparison ------------------------------------------------
    datHybrid <- dat %>%
        filter(!is.na(shubhamFreq)) %>%
        mutate(stabClass = ifelse(shubhamFreq >= hybridTresh, "High", "Low")) %>%
        mutate(nLevels = factor(nLevels, levels = levelVec*2))
    datAbs <- datHybrid %>% filter(PhenAbs == "Hybrid") %>%
        group_by(nLevels, stabClass) %>%
        summarise(nStates = n(), Frequency = sum(shubhamFreq), meanFrust = mean(frust0),
            .groups = "drop")
    totals <- datAbs %>%
        group_by(nLevels) %>%
        summarise(tStates = sum(nStates), tFrequency = sum(Frequency), tFrust = mean(meanFrust),
            .groups = "drop")
    datAbs <- merge(datAbs, totals, by = "nLevels")

    datThresh <- datHybrid %>% filter(PhenThresh == "Hybrid") %>%
        group_by(nLevels, stabClass) %>%
        summarise(nStates = n(), Frequency = sum(shubhamFreq), meanFrust = mean(frust0),
            .groups = "drop")
    totals <- datThresh %>%
        group_by(nLevels) %>%
        summarise(tStates = sum(nStates), tFrequency = sum(Frequency), tFrust = mean(meanFrust),
            .groups = "drop")
    datThresh <- merge(datThresh, totals, by = "nLevels")
    ggplot(datAbs, aes(x = nLevels, y = Frequency, fill = stabClass)) +
        geom_bar(stat = "identity", position = position_dodge(preserve = "single")) +
        labs(x = "Number of levels", y = "Frequency", fill = "Stability") +
        theme_Publication() +
        theme(legend.position = "top", legend.direction = "horizontal",
            legend.key.width = unit(0.6, "cm"))
    ggsave(paste0(net, "_hybridFreqAbsAll.png"), width = 5.5, height = 5)

    ggplot(datThresh, aes(x = nLevels, y = Frequency, fill = stabClass)) +
        geom_bar(stat = "identity", position = position_dodge(preserve = "single")) +
        labs(x = "Number of levels", y = "Frequency", fill = "Stability") +
        theme_Publication() +
        theme(legend.position = "top", legend.direction = "horizontal",
            legend.key.width = unit(0.6, "cm"))
    ggsave(paste0(net, "_hybridFreqThreshAll.png"), width = 5.5, height = 5)

    ggplot(datAbs, aes(x = nLevels, y = nStates, fill = stabClass)) +
        geom_bar(stat = "identity", position = position_dodge(preserve = "single")) +
        scale_y_log10() +
        labs(x = "Number of levels", y = "Number of States", fill = "Stability") +
        theme_Publication() +
        theme(legend.position = "top", legend.direction = "horizontal",
            legend.key.width = unit(0.6, "cm"))
    ggsave(paste0(net, "_hybridStatesAbsAll.png"), width = 5.5, height = 5)

    ggplot(datThresh, aes(x = nLevels, y = nStates, fill = stabClass)) +
        geom_bar(stat = "identity", position = position_dodge(preserve = "single")) +
        scale_y_log10() +
        labs(x = "Number of levels", y = "Number of States", fill = "Stability") +
        theme_Publication() +
        theme(legend.position = "top", legend.direction = "horizontal",
            legend.key.width = unit(0.6, "cm"))
    ggsave(paste0(net, "_hybridStatesThreshAll.png"), width = 5.5, height = 5)

    setwd("..")

    ### Segment plots --------------------------------------------------------------
    setwd("multiCompare")
    p <- segmentPlots(dat)
    ggsave(paste0(net, "_segmentPlots.png"), width = 6.5, height = 5.7)
    p1 <- segmentPlots(dat, relative = T)
    ggsave(paste0(net, "_segmentPlotsRelative.png"), width = 6.5, height = 5.7)
    setwd("..")
}


singleFilePlots <- function(fl, nodeOrder, hybridTresh = 0.01) {
    print(fl)
    d <- read_csv(fl, show_col_types = F)
    nLevels <- str_extract(fl, "_shubham_[0-9]+") %>%
        str_remove("_shubham_") %>% as.numeric
    net <- str_extract(fl, ".*(?=_shubham)")
    MasterKey <- paste0(net, "_", nLevels)
    DirectoryNav(paste0("nLevels_", nLevels))
    #### heatmaps ----
    stateExpressionHeatmap(d, nodeOrder, paste0(MasterKey, "_expressionHeatmap.png"))
    # high frequency hybrids
    d1 <- d %>% filter(PhenAbs == "Hybrid") %>%
        filter(Avg0 >= hybridTresh)
    if (nrow(d1) > 2)
    stateExpressionHeatmap(d1, nodeOrder, paste0(MasterKey, "_topHybridAbs.png"))
    d1 <- d %>% filter(PhenThresh == "Hybrid") %>%
        filter(Avg0 >= hybridTresh)
    if (nrow(d1) > 2)
    stateExpressionHeatmap(d1, nodeOrder, paste0(MasterKey, "_topHybridThresh.png"))
    print("state expression done.")
    #### frust vs score plots ----
    frustScorePlots(d, MasterKey)
    print("frust score done.")
    setwd("..")
}

segregatedPlots <- function(net, levelVec = c(1,2, 100)) {
    wd <- getwd()
    nodes <- readLines(paste0(net, "_nodes.txt"))
    teams <- readLines(paste0(net, ".teams")) %>%
        str_split(",")
    nodeOrder <- c(unlist(teams), nodes[!nodes %in% unlist(teams)] %>% sort %>% unique)
    scoreSegregation(net)
    scoringMethods <- c("stateSum")
    if (!dir.exists(paste0("Figures/", net))) {
        dir.create(paste0("Figures/", net), recursive = T)
    }
    setwd(paste0("Figures/", net))
    sapply(scoringMethods, function(x) {
        DirectoryNav(x)
        filez <- paste0(net, "_shubham_", levelVec, "_finFlagFreq_format.csv")
        # sapply(filez, function(f) {
        #     singleFilePlots(f, nodeOrder)
        # })
        # compareHybridSS(net, levelVec, ".")
        isingComparePlots(net, levelVec)
        setwd("..")
    })
    setwd(wd)
}

