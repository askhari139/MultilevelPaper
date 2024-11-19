
library(funcsKishore)

nTeams <- 3

nNodes <- 10
nodes <- list(Team1 = paste0("T1_", 1:nNodes), 
    Team2 = paste0("T2_", 1:nNodes), 
    Team3 = paste0("T3_", 1:nNodes))
densities <- c(0.3, 0.9)

sapply(densities, function(density) {
    df <- lapply(1:nTeams, function(t1) {
        d <- lapply(1:nTeams, function(t2) {
            n1 <- nodes[[t1]]
            n2 <- nodes[[t2]]
            edge <- 2
            if (t1 == t2) {
                edge <- 1
            }
            n <- round(nNodes*nNodes*density)
            Source <- sample(n1, n, replace = TRUE)
            Target <- sample(n2, n, replace = TRUE)
            data.frame(Source = Source, Target = Target, Type = edge) %>% unique
        }) %>% bind_rows
    }) %>% bind_rows
    write_delim(df, paste0("ThreeTeams_", density, ".topo"), delim = " ")
    nodesAll <- unique(c(df$Source, df$Target))
    teams <- sapply(1:nTeams, function(t) {
        nd <- nodes[[t]]
        nd <- nd[nd %in% nodesAll]
        nd %>% paste(collapse = ",")
    })
    writeLines(teams, paste0("ThreeTeams_", density, ".teams"))
})
