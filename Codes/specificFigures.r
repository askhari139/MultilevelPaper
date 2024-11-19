# source("/Users/kishorehari/Desktop/Postdoc/MultiLevel/CoreData/setupScript.R")

setwd(mlNewIneq)
setwd("Figures/EMT_RACIPE2_noPeri/stateSum")
df <- read_csv("EMT_RACIPE2_noPeri_shubham_2finFlagFreq_format.csv")

# hybrid only
d1 <- df %>% filter(PhenAbs == "Hybrid")
d1 <- d1 %>% select(Avg0) %>% mutate(Avg0 = as.numeric(Avg0)) %>%
    mutate(Avg0 = log10(Avg0))
ggplot(d1, aes(x = Avg0)) + geom_histogram(binwidth = 0.5, aes(y = ..count../sum(..count..))) + 
    labs(x = "log10(Hybrid SSF)", y = "Probability") + 
    theme_Publication()
ggsave("nLevels_2/HybridHist.png", width = 5.5, height = 5)

# all
d1 <- df
d1 <- d1 %>% select(Avg0) %>% mutate(Avg0 = as.numeric(Avg0)) %>%
    mutate(Avg0 = log10(Avg0))
ggplot(d1, aes(x = Avg0)) + geom_histogram(binwidth = 0.5, aes(y = ..count../sum(..count..))) + 
    labs(x = "log10(SSF)", y = "Probability") + 
    theme_Publication()
ggsave("nLevels_2/AllHist.png", width = 5.5, height = 5)

# terminal only
d1 <- df %>% filter(str_detect(PhenAbs, "Terminal"))
d1 <- d1 %>% select(Avg0) %>% mutate(Avg0 = as.numeric(Avg0)) %>%
    mutate(Avg0 = log10(Avg0))
ggplot(d1, aes(x = Avg0)) + geom_histogram(binwidth = 0.5, aes(y = ..count../sum(..count..))) +
    labs(x = "log10(Terminal SSF)", y = "Probability") + 
    theme_Publication()
ggsave("nLevels_2/TerminalHist.png", width = 5.5, height = 5)


### number of levels vs em score max

setwd(mlNewIneq)

setwd("Figures/EMT_RACIPE2_noPeri/stateSum")

nLevels <- c(1:10, 20, 50, 100)
flz <- paste0("EMT_RACIPE2_noPeri_shubham_", nLevels, "_finFlagFreq_format.csv")
scores <- sapply(flz, function(f) {
    read_csv(f) %>% select(emScore) %>% unlist %>% max
})

scores <- data.frame(nLevels = 1/nLevels, scores = scores)
ggplot(scores, aes(x = nLevels, y = scores)) + geom_point() + 
    geom_smooth(se = T) +
    labs(x = "Precision of the formalism", y = "Max emScore") + 
    theme_Publication()
ggsave("nLevelsVsEmScore.png", width = 5.5, height = 5)
