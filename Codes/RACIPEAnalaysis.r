library(funcsKishore)

frustCalcIsing <- function(state, nodeOrder, topoDf)
{#browser()
    names(state) <- nodeOrder
    nEdges <- nrow(topoDf)
    frust <- topoDf %>%
        mutate(Source = state[Source], Target = state[Target], 
            Type = ifelse(Type == 2, -1, 1)) %>%
        mutate(Sum = Source*Target*Type) %>%
        mutate(frust = ifelse(Sum<=0, 1, 0)) %>%
        select(frust) %>%
        unlist %>% sum
    frust/nEdges
}
frustCalcIsing <- cmpfun(frustCalcIsing)

frustCalcMulti <- function(state, nodeOrder, topoDf)
{#browser()
    names(state) <- nodeOrder
    nEdges <- nrow(topoDf)
    frust <- topoDf %>%
        mutate(Source = state[Source], Target = state[Target], 
            Type = ifelse(Type == 2, -1, 1)) %>%
        mutate(Sum = Source*Target*Type) %>%
        mutate(frust = ifelse(Sum<=0, abs(Sum), 0)) %>%
        select(frust) %>%
        unlist %>% sum
    frust/nEdges
}
frustCalcMulti <- cmpfun(frustCalcMulti)

formatRACIPE <- function(topoFile) {
    setwd("RACIPE")
    net <- topoFile %>% str_remove(".topo")
    topoDf <- read.delim(topoFile, sep = "", stringsAsFactors = F)
    deletables <- paste0(net, c("_solution_", "_T_test_", ".cfg"))
    fileList <- sapply(deletables, function(x) {
        fAll <- list.files(pattern = x)
        sapply(fAll, file.remove)
    })
    if(!file.exists(paste0(net, "_solution.dat"))) {
        print(paste0("Simulate!"))
        setwd("..")
        return(NA)
    }
    paramNames <- read.delim(paste0(net, ".prs"), sep = "") %>%
        select(Parameter) %>%
        unlist
    nodes <- paramNames[str_detect(paramNames, "Prod_of")] %>%
        str_remove("Prod_of_")
    dat <- read_delim(paste0(net, "_solution.dat"), 
        delim = "\t", col_names = F, col_types = "d")
    colnames(dat) <- c("paramID", "nStates", "Basin", nodes)
    parameters <- read_delim(paste0(net, "_parameters.dat"),
        delim = "\t", col_names = F, col_types = "d") %>%
        set_names(c("paramID", "nStates", paramNames))
    normFactor <- (parameters %>%
        select(all_of(paste0("Prod_of_", nodes))) %>% set_names(nodes))/(
            parameters %>% select(all_of(paste0("Deg_of_", nodes))) %>% set_names(nodes)) 
    normFactor <- normFactor %>%
        mutate(paramID = parameters$paramID) %>%
        merge(dat %>% select(paramID), by = "paramID")
    
    sapply(nodes, function(x) {
        dat[[x]] <<- round((2^dat[[x]])/normFactor[[x]], 2)
        return(NULL)
    })
    if (!file.exists(paste0(net, ".teams"))) {
        getGsVec(topoFiles = topoFile)
    }
    teams <- readLines(paste0(net, ".teams")) %>%
        str_split(",")
    dat <- dat %>%
        mutate(eScore = dat %>% select(all_of(teams[[1]])) %>%
                   sapply(as.numeric) %>% rowSums,
               mScore = dat %>% select(all_of(teams[[2]])) %>%
                   sapply(as.numeric) %>% rowSums) %>%
        mutate(eScore = eScore/length(teams[[1]]),
               mScore = mScore/length(teams[[2]]),
            emScore = eScore - mScore) %>%
        mutate(PhenAbs = ifelse(abs(emScore) == 1, "Terminal", "Hybrid"),
                PhenAbs = ifelse(eScore*mScore == 0 & abs(emScore) < 1, 
                    paste0("Partial ", PhenAbs), PhenAbs),
                PhenAbs = ifelse(PhenAbs == "Partial Hybrid", "Partial Terminal", PhenAbs)) %>%
        mutate(PhenThresh = ifelse(abs(emScore) <= 0.5, "Hybrid", "Terminal"),
                PhenThresh = ifelse(eScore*mScore == 0 & abs(emScore) < 0.5, 
                    paste0("Partial ", PhenThresh), PhenThresh),
                PhenThresh = ifelse(PhenThresh == "Partial Hybrid",
                    "Partial Terminal", PhenThresh))
    nMax <- dat %>% filter(paramID == 1) %>%
        select(Basin) %>% unlist %>% sum
    dat <- dat %>% mutate(Basin = Basin/nMax)
    write_csv(dat, paste0(net, "_compiled.csv"))
    nodeOrder <- c(unlist(teams), nodes[!nodes %in% unlist(teams)] %>% unique)
    dat <- dat %>%
        mutate(states = NA, Avg0 = Basin) %>%
        group_by(across(all_of(nodeOrder))) %>%
        summarise(states = NA, Avg0 = sum(Avg0), eScore = mean(eScore),
            mScore = mean(mScore), emScore = mean(emScore),
            PhenAbs = first(PhenAbs), PhenThresh = first(PhenThresh),
            .groups = "drop")
    dat <- dat %>%
            mutate(
            frust0 = apply(dat %>% select(all_of(nodeOrder)), 1, function(x) {
                frustCalcMulti(x, nodeOrder, topoDf)
            }),
            isingFrust = apply(dat %>% select(all_of(nodeOrder)), 1, function(x) {
                frustCalcIsing(x, nodeOrder, topoDf)
            }),
            Avg0 = Avg0/nrow(parameters)) %>%
        select(states, all_of(nodeOrder), Avg0, frust0,
            isingFrust, eScore, mScore, emScore, PhenAbs, PhenThresh)
    write_csv(dat, paste0("../", net, "_shubham_100_finFlagFreq_format.csv"))
    setwd("..")
    return(paste0(net, " Done!"))
}

heatmapPlot <- function(topoFile) {
    net <- topoFile %>% str_remove(".topo")
    if(!file.exists(paste0(net, "_compiled.csv"))) {
        d <- formatRACPIE(topoFile)
        if(is.na(d)) {
            return(NA)
        }
    }
    dat <- read_csv(paste0(net, "_compiled.csv"))
    paramNames <- read.delim(paste0(net, ".prs"), sep = "") %>%
        select(Parameter) %>%
        unlist
    nodes <- paramNames[str_detect(paramNames, "Prod_of")] %>%
        str_remove("Prod_of_")
    teams <- readLines(paste0(net, ".teams")) %>%
        str_split(",") %>% unlist
    DirectoryNav("Plots")
    nodes <- c(teams, nodes[!nodes %in% teams] %>% unique)
    heatMapDat <- dat %>% arrange(emScore) %>%
        mutate(ID = 1:nrow(.)) %>%
        mutate(ID = factor(ID %>% as.character, levels = ID %>% as.character)) %>%
        select(ID, all_of(nodes)) %>%
        gather(key = "Node", value = "Expression", -ID)
    ggplot(heatMapDat, aes(x = Node, y = ID, fill = Expression)) +
        geom_tile() +
        scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0.5) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1),
            axis.text.y = element_blank(),
            legend.position = "right", legend.direction = "vertical",
            legend.key.height = unit(0.8, "cm")) +
        labs(x = "Node", y = "Parameter ID", fill = "Expression")
    ggsave(paste0(net, "_heatmap.png"), width = 6, height = 10)
     
    setwd("..")
}