##################################################################
# Titel:      How Corpus-Based Methods Can Support Language Learning 
#             and Teaching: An Analysis of Amplifier Use by L1-Speakers 
#             and Learners of English - Part 3
# R version:  3.5.1 (2018-07-02) -- "Feather Spray"
# Autor:      Martin Schweinberger
# Date:       2019-02-02
# Contact:    martin.schweinberger.hh@gmail.com
# Disclaimer: If you have questions,suggestions or you found errors
#             or in case you would to provide feedback, questions
#             write an email to martin.schweinberger.hh@gmail.com.
# Citation:   If you use this script or results thereof, please cite it as:
#             Schweinberger, Martin. 2018. "How Corpus-Based Methods 
#             Can Support Language Learning and Teaching: An Analysis 
#             of Amplifier Use by L1-Speakers and Learners of English - Part 3",
#             unpublished R script, The University of Queensland.
###############################################################
rm(list=ls(all=T))                                      # clean current workspace
setwd("D:\\Uni\\Projekte\\02-Intensification\\Amp2L")   # set wd
library(plyr)                                           # load packages
library(Rling)                                          # load packages
options(stringsAsFactors = F)                           # set options
options(scipen = 999)                                   # set options
options(max.print=10000)                                # set options
imageDirectory<-"images"                                # define image directory
###############################################################
# read in data
ampicle <- read.table("ampicle_clean04.txt", sep = "\t", header=TRUE)
# remove superfluous columns
ampicle$Subfile <- NULL
ampicle$Speaker <- NULL
head(ampicle)  # inspect data

amplocness <- read.table("amplocness_clean04.txt", sep = "\t", header=TRUE)
head(amplocness)  # inspect data

colnames(ampicle); colnames(amplocness)

# combine data sets
ampicle <- rbind(ampicle, amplocness)
head(ampicle)  # inspect data
str(ampicle)

###############################################################
# factorize variables
clfct <- c("Function", "Priming", "Gradabilty", "SemanticCategory", "Emotionality")
ampicle[clfct] <- lapply(ampicle[clfct], factor)
###############################################################
#              SEMANTIC VECTOR SPACE MODEL 1
library(tm)
library(dplyr)
svsmd <- ampicle %>%
  select(Adjective, Variant)
t2 <- table(svsmd$Adjective, svsmd$Variant) # tabulate data
t3 <- t(t2)
t3icle <- t3
t3icle <- t3icle[, 3: ncol(t3icle)]
# remove Adjectives that were not Amplified
t3icle <- t3icle[rowSums(t3icle) > 0, ]
# save row and column names
colnamesicle <- colnames(t3icle)
rownamesicle <- rownames(t3icle)
# turn dataframe Amplifiedo matrix
svsmicle <- as.matrix(t3icle)
# convert token frequency to type frequency
svsmicle <- apply(svsmicle, 1, function(x) { 
  x <- ifelse(x > 1, 1, x) } )
svsmicle <- t(svsmicle)
#svsmicle <- svsmicle[, colSums(svsmicle) >= 2]
#svsmicle <- svsmicle[rowSums(svsmicle) >= 2, ]
svsmicle

# compute expected values
svsmicle.exp <- chisq.test(svsmicle)$expected
# calculate PMI and PPMI
svsmicle.PMI <- log2(svsmicle/svsmicle.exp)
svsmicle.PPMI <- ifelse(svsmicle.PMI < 0, 0, svsmicle.PMI)
# calculate cosine similarity
svsmicle.tmp1 <- svsmicle.PPMI
svsmicle.cos <- cossim(svsmicle.tmp1)
#round(svsmicle.cos, 2)
###############################################################
#               CLUSTER SEMANTIC VECTORS
# load library
library(cluster)
# find max value that is not 1
svsmicle.cos.test <- apply(svsmicle.cos, 1, function(x){
  x <- ifelse(x == 1, 0, x) } )
maxval <- max(svsmicle.cos.test)
# create distance matrix
svsmicle.dist <- 1 - (svsmicle.cos/maxval)
clustd <- as.dist(svsmicle.dist)
# create distance matrix
clustd <- dist(svsmicle.cos, method = "manhattan") 
# alternative methods
# eucledian - not good when dealing with many dimensions
# manhattan - most popular choice
# method - here the difference between points dominates
# canberra - for count data
# binary - for binary data only!
# minkowski - is not a true distance measure

# find optimal number of clusters
asw <- as.vector(unlist(sapply(2:nrow(svsmicle)-1, function(x) pam(clustd, k = x)$silinfo$avg.width)))
# determine the optimal number of clusters (max width is optimal)
optclust <- which(asw == max(asw))+1 # optimal number of clusters

# inspect clustering with optimal number of clusters
svsmicle.clust <- pam(clustd, 2)
svsmicle.clust$clustering

# create cluster object
# alternative methods: "single", "ward.D2", "averAge", "mcquitty", "median", "centroid"
ampiclehclust <- hclust(clustd, method="ward.D")    
# plot cluster solution
png("images/Clusticle.png",  width = 680, height = 480) # save plot
plot(ampiclehclust, main = "", xlab = "", ylab = "", cex = .8)
rect.hclust(ampiclehclust, k = 2, border= "orange")
dev.off()
# load libraries for nicer dendrograms
library(factoextra)
library(dendextend)
# plot with colored clusters
opar <- par(mar = c(5, 4, 4, 2) + 0.1)      # make a copy of current settings
npar <- par(mar = c(3, 1, 1, 20))
png("images/Clusticle2.png",  width = 400, height = 800) # save plot
npar <- par(mar = c(3, 1, 1, 20))
# plot with colored clusters
fviz_dend(ampiclehclust, k = 2, cex = 1, horiz = T,  type = "rectangle",
          k_colors = c("grey60", "grey20"), 
          rect_border = "white",#c("orange", "orange", "orange", "grey30"), 
          rect_fill = F, main = "", labels_track_height=10, rect = T, margin = c(5,20))
dev.off()
par(opar)          # restore original settings 
# plot as unrooted tree
png("images/PhyClustAmpicle.png",  width = 680, height = 480) 
fviz_dend(ampiclehclust, k = 2, color_labels_by_k = T, type = "phylogenic", repel = TRUE, cex = .9,
          k_colors = c("grey70",  "grey30"))
dev.off()
###############################################################
# Unrooted clustering
# library ape
library(ape)
# convert 'hclust' to 'phylo' object
phylo_tree = as.phylo(ampiclehclust)
# get edges
graph_edges = phylo_tree$edge
# library igraph
library(igraph)
# get graph from edge list
graph_net = graph.edgelist(graph_edges)
# extract layout (x-y coords)
graph_layout = layout.auto(graph_net)
# number of observations
nobs = nrow(svsmicle.cos)
# save plot
png("images/UClustAmpicle.png",  width = 680, height = 480) 
# start plot
plot(graph_layout[,1], graph_layout[,2], type = "n", axes = FALSE,
     xlab = "", ylab = "")
# draw tree branches
segments(
  x0 = graph_layout[graph_edges[,1],1], 
  y0 = graph_layout[graph_edges[,1],2],
  x1 = graph_layout[graph_edges[,2],1],
  y1 = graph_layout[graph_edges[,2],2],
  col = "gray90", lwd = 2
)
# add labels
text(graph_layout[1:nobs,1], graph_layout[1:nobs,2],
     phylo_tree$tip.label, cex = .9, xpd = TRUE, font = 1)
dev.off()
###############################################################
#                 WARNING
#             DATA REDUCTION
# exclude amplifiers that are very dissimilar to main group of amplifiers
rmvamp <- c("considerably", "significantly", "profoundly", "strickingly", "substantially", "positively")
nrow(ampicle)

ampicle <- ampicle[!ampicle$Variant %in% rmvamp, ]
ampicle <- ampicle[!ampicle$Variant == "", ]
ampicle <- ampicle[complete.cases(ampicle),]
nrow(ampicle)

###############################################################
# prepare data for plotting
# create data frame with relevant variables
pd <- data.frame(ampicle$Language, ampicle$Function, ampicle$Adjective, ampicle$Amplified, ampicle$Variant)
# clean col names
colnames(pd) <- gsub("ampicle.", "", colnames(pd))
# convert Age column
LanguageLbs <- names(table(pd$Language))
# multiply Amplified * 100 to get percent for Variant
pd$Amplified <- ifelse(pd$Amplified == 1, 100, 0)
# convert Amplified Amplifiedo a numeric variables
pd$Amplified <- as.numeric(pd$Amplified)
famps <- names(table(pd$Variant))[which(table(pd$Variant) > 250)]
# reclassify Adjectives - infreq. Adjectives are collapsed Amplifiedo category other
pd$Variant <- ifelse(pd$Variant  %in% famps, pd$Variant , "other")
# create variables 
pd$other <- ifelse(pd$Variant == "other", 100, 0)
pd$extremely <- ifelse(pd$Variant == "extremely", 100, 0)
pd$really <- ifelse(pd$Variant == "really", 100, 0)
pd$so <- ifelse(pd$Variant == "so", 100, 0)
pd$very <- ifelse(pd$Variant == "very", 100, 0)
pd$zero <- ifelse(pd$Variant == "0", 100, 0)
pd <- pd[complete.cases(pd),]
pd$Language <- factor(pd$Language, levels=c("English", "Bulgarian", "Czech",
                                            "Dutch", "Finnish", "Flemish", "French", 
                                            "German", "Italian", "Polish", "Russian", 
                                            "Spanish", "Swedish"))
###############################################################
# p0
p0d <- pd
# start plot: Amplified
p0 <- ggplot(p0d, aes(Language, Amplified)) +
  geom_point(aes(reorder(Language, Amplified, function(Amplified) -mean(Amplified)), y=Amplified), size = NA) +
  stat_summary(fun.y = mean, geom = "point", size = 2) +
  stat_summary(fun.y = mean, geom = "line", size = .5) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.2, size = .5) +
  theme_set(theme_bw(base_size = 10)) +
  coord_cartesian(ylim = c(0, 20)) +
  theme(legend.position="none", axis.text.x = element_text(size=8, angle=90)) +
  labs(x = "Language", y = "Percent (Amplified Adjectives)") +
  scale_color_manual(values = c("grey30", "grey30"))
imageFile <- paste(imageDirectory,"AmplifiedLanguage.png",sep="/")
ggsave(file = imageFile)
p0

###############################################################
# p1
p1d <- pd
# start plot: Amplified
p1 <- ggplot(p1d, aes(Language, Amplified)) +
  facet_grid(vars(Function)) +
  geom_point(aes(reorder(Language, Amplified, function(Amplified) -mean(Amplified)), y=Amplified), size = NA) +
  stat_summary(fun.y = mean, geom = "point", size = 2) +
  stat_summary(fun.y = mean, geom = "line", size = .5) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.2, size = .5) +
  theme_set(theme_bw(base_size = 10)) +
  coord_cartesian(ylim = c(0, 20)) +
  theme(legend.position="none", axis.text.x = element_text(size=8, angle=90)) +
  labs(x = "Language", y = "Percent (Amplified Adjectives)") +
  scale_color_manual(values = c("grey30", "grey30"))
# activate (remove #) to save
imageFile <- paste(imageDirectory,"AmplifiedLanguageFunction.png",sep="/")
ggsave(file = imageFile)
# activate (remove #) to show
p1

###############################################################
library(dplyr)
# p2
p2d <- pd
# remove non-amplified instances
p2d <- p2d[p2d$Amplified != 0,]
p2d <- p2d %>%
  group_by(Function, Language) %>%
  summarise(mean_very = mean(very),
            mean_really = mean(really),
            mean_so = mean(so),
            mean_extremely = mean(extremely),
            mean_other = mean(other))
# start plot: all
p2 <- ggplot(p2d, aes(Language, mean_very)) +
  facet_grid(vars(Function)) +
  geom_point(aes(y = mean_very, color = "very"), size=1) +
  geom_point(aes(y = mean_really, color = "really"), size=1) +
  geom_point(aes(y = mean_so, color = "so"), size=1) +
  geom_point(aes(y = mean_extremely, color = "extremely"), size=1) +
  geom_point(aes(y = mean_other, color = "other"), size=1) +
  guides(color=guide_legend(override.aes=list(fill=NA))) +
  scale_colour_manual(values=c("gray50", "goldenrod2", "gray70", "indianred4", "grey30"),
                      name="Variant", 
                      breaks=c("other", "extremely", "really", "so", "very"), 
                      labels = c("other", "extremely", "really", "so", "very")) +
  theme_set(theme_light(base_size = 10)) +
  theme(legend.position="top") +
  coord_cartesian(ylim = c(0, 100)) +
  labs(x = "Language", y = "Percent of Amplification") +
  guides(size = FALSE)+
  guides(alpha = FALSE)
ggsave(file = paste(imageDirectory,"VariantLanguageFunction.png",sep="/"))
p2

###############################################################
# p3
p3d <- pd
# remove non-amplified instances
p3d <- p3d[p3d$Amplified != 0,]
# start plot: all with zero
p3 <- ggplot(p3d, aes(x = Language, y = very)) +
  facet_grid(vars(Function)) +
  geom_smooth(aes(y = very, color = "very", linetype = "very"), size=.25, se = F) +
  geom_smooth(aes(y = really, color = "really", linetype = "really"), size=.25, se = F) +
  geom_smooth(aes(y = so, color = "so", linetype = "so"), size=.25, se = F) +
  geom_smooth(aes(y = extremely, color = "extremely", linetype = "extremely"), size=.25, se = F) +
  geom_smooth(aes(y = other, color = "other", linetype = "other"), size=.25, se = F) +
  guides(color=guide_legend(override.aes=list(fill=NA))) +
  scale_linetype_manual(values=c("longdash","twodash", "dashed", "dotdash", "solid"),
                        name="Variant",
                        breaks = c("other", "extremely", "really", "so", "very"), 
                        labels = c("other", "extremely", "really", "so", "very")) +
  scale_colour_manual(values=c("gray50", "goldenrod2", "gray70", "indianred4", "grey30"),
                      name="Variant", 
                      breaks=c("other", "extremely", "really", "so", "very"), 
                      labels = c("other", "extremely", "really", "so", "very")) +
  theme_set(theme_light(base_size = 10)) +
  theme(legend.position="top") +
  coord_cartesian(ylim = c(0, 100)) +
  labs(x = "Age", y = "Percent of Amplification") +
  guides(size = FALSE)+
  guides(alpha = FALSE)+
ggsave(file = paste(imageDirectory,"VariantPdFunction.png",sep="/"), width = 10, height = 10, units = c("cm"),  dpi = 320)
p3

###############################################################
#            WARNING: DATA REDUCTION
###############################################################
tbd <- ampicle
# recode adjectives
ntfrqadj <- names(table(tbd$Adjective))[which(table(tbd$Adjective) <= 2000)]
tbd$Adjective <- ifelse(tbd$Adjective %in% ntfrqadj, "other", tbd$Adjective)
###############################################################
###             TABULARIZATION
###############################################################
# tb1
tb1 <- tbd %>%
  group_by(Variant) %>%
  summarize(
    Tokens = n()
    ) %>%
  mutate(Percent = round(Tokens/sum(Tokens)*100, 2)) %>%
  mutate(PercentAmplifiers = c("", round(Tokens[2:length(Tokens)]/sum(Tokens[2:length(Tokens)])*100, 2)))
tb1 <- tb1[order(tb1$Percent, decreasing = T),]
tb1

# save data to disc
write.table(tb1, "Table1.txt", sep = "\t", row.names = F)
###############################################################
# tb2
tb2 <- select(tbd, Variant, Language)
freqamp <- c("0", "very", "so", "really", "completely")
tb2$Variant <- ifelse(tb2$Variant %in% freqamp, tb2$Variant, "other")
tb2 <- ftable(tb2$Variant, tb2$Language)
tb2

# save data to disc
write.table(tb2, "Table2.txt", sep = "\t", row.names = F)
###############################################################
# tb3
tb3 <- tbd[tbd$Variant != "0",]
tb3 <- tb3 %>%
  group_by(Adjective, Language) %>%
  mutate(
    very = Variant == "very",
    so = Variant == "so",
    really = Variant == "really",
    extremely = Variant == "extremely",
    other = Variant == "other") %>%
  summarise(very = sum(very),
            really = sum(really),
            so = sum(so),
            extremely = sum(extremely),
            other = sum(other))
tb3

# save data to disc
write.table(tb3, "Table3.txt", sep = "\t", row.names = T)
###############################################################
# tb4
# reorder data
tb4 <- tbd[tbd$Variant != "0",]
tb4 <- tb4 %>%
  group_by(Adjective, Language) %>%
  mutate(
    very = Variant == "very",
    so = Variant == "so",
    really = Variant == "really",
    extremely = Variant == "extremely",
    other = Variant == "other") %>%
  summarise(very = sum(very),
            really = sum(really),
            so = sum(so),
            extremely = sum(extremely),
            other = sum(other))
tb4$RowTotal <- rowSums(tb4[,c(3:ncol(tb4))])
tb4 <- tb4 %>%
  mutate(
    verypercent = round(very/RowTotal*100, 2),
    sopercent = round(so/RowTotal*100, 2),
    reallypercent = round(really/RowTotal*100, 2),
    extremelypercent = round(extremely/RowTotal*100, 2),
    otherpercent = round(other/RowTotal*100, 2))
tb4 <- tb4 %>%
  select(Adjective, Language, verypercent, sopercent, reallypercent, 
         extremelypercent, otherpercent)
tb4

# save data to disc
write.table(tb4, "Table4.txt",  sep = "\t", row.names = T)
###############################################################
tb5 <- tbd[tbd$Variant != "0",]
tb5 <- select(tb5, Variant, Language)
freqamp <- c("0", "very", "so", "really", "completely")
tb5$Variant <- ifelse(tb5$Variant %in% freqamp, tb5$Variant, "other")
tb5 <- ftable(tb5)

png("images/mosaic.png",  width = 400, height = 800) 
mosaicplot(tb5)
dev.off()
png("images/assocplot.png",  width = 500, height = 2000) 
assocplot(as.matrix(tb5))
dev.off()

###############################################################
tb6 <- tbd[tbd$Variant != "0",]
tb6 <- select(tb6, Variant, Language)
freqamp <- c("very", "so", "really", "completely")
tb6$Variant <- ifelse(tb6$Variant %in% freqamp, tb6$Variant, "other")
library(cfa)             # activate library
cfatb <- cfa(tb6)        # perform cfa
write.table(cfatb$table, "CFA.txt",  sep = "\t", row.names = T)
###############################################################
tb7 <- tbd[tbd$Variant != "0",]
tb7 <- select(tb7, Variant, Language)
freqamp <- c("very", "so", "really", "completely")
tb7$Variant <- ifelse(tb7$Variant %in% freqamp, tb7$Variant, "other")
tb7 <- tb7 %>%
  count(Language, Variant) %>%
  group_by(Language)
tball <- tb7 %>%
  group_by(Language) %>%
  mutate(
    Amplifiers = sum(n)
  )
tb7$AmplifiersAll <- tball$Amplifiers
tb7$Amplifiers <- tb7$AmplifiersAll-tb7$n
tb7$AmplifiersAll <- NULL
# Chi-Square for 2*k Tables
source("D:\\R/x2.2k.r")  # load function
# very
tbvery <- tb7 %>%
  filter(Variant == "very")
tbvery$Language <- NULL
tbvery$Variant <- NULL
tbvery <- as.data.frame(tbvery)
tbvery$n <- as.numeric(tbvery$n)
tbvery$Amplifiers <- as.numeric(tbvery$Amplifiers)
# perform tests for very
BulgarianGBvery <- x2.2k(tbvery, 1, 4)
CzechGBvery <- x2.2k(tbvery, 2, 4)
DutchGBvery <- x2.2k(tbvery, 3, 4)
FinnishGBvery <- x2.2k(tbvery, 5, 4)
FlemishGBvery <- x2.2k(tbvery, 6, 4)
FrenchGBvery <- x2.2k(tbvery, 7, 4)
GermanGBvery <- x2.2k(tbvery, 8, 4)
ItalianGBvery <- x2.2k(tbvery, 9, 4)
PolishGBvery <- x2.2k(tbvery, 10, 4)
RussianGBvery <- x2.2k(tbvery, 11, 4)
SpanishGBvery <- x2.2k(tbvery, 12, 4)
SwedishGBvery <- x2.2k(tbvery, 13, 4)
# so
tbso <- tb7 %>%
  filter(Variant == "so")
tbso$Language <- NULL
tbso$Variant <- NULL
tbso <- as.data.frame(tbso)
tbso$n <- as.numeric(tbso$n)
tbso$Amplifiers <- as.numeric(tbso$Amplifiers)
# perform tests for so
BulgarianGBso <- x2.2k(tbso, 1, 4)
CzechGBso <- x2.2k(tbso, 2, 4)
DutchGBso <- x2.2k(tbso, 3, 4)
FinnishGBso <- x2.2k(tbso, 5, 4)
FlemishGBso <- x2.2k(tbso, 6, 4)
FrenchGBso <- x2.2k(tbso, 7, 4)
GermanGBso <- x2.2k(tbso, 8, 4)
ItalianGBso <- x2.2k(tbso, 9, 4)
PolishGBso <- x2.2k(tbso, 10, 4)
RussianGBso <- x2.2k(tbso, 11, 4)
SpanishGBso <- x2.2k(tbso, 12, 4)
SwedishGBso <- x2.2k(tbso, 13, 4)
# really
tbreally <- tb7 %>%
  filter(Variant == "really")
tbreally$Language <- NULL
tbreally$Variant <- NULL
tbreally <- as.data.frame(tbreally)
tbreally$n <- as.numeric(tbreally$n)
tbreally$Amplifiers <- as.numeric(tbreally$Amplifiers)
# perform tests for really
BulgarianGBreally <- x2.2k(tbreally, 1, 4)
CzechGBreally <- x2.2k(tbreally, 2, 4)
DutchGBreally <- x2.2k(tbreally, 3, 4)
FinnishGBreally <- x2.2k(tbreally, 5, 4)
FlemishGBreally <- x2.2k(tbreally, 6, 4)
FrenchGBreally <- x2.2k(tbreally, 7, 4)
GermanGBreally <- x2.2k(tbreally, 8, 4)
ItalianGBreally <- x2.2k(tbreally, 9, 4)
PolishGBreally <- x2.2k(tbreally, 10, 4)
RussianGBreally <- x2.2k(tbreally, 11, 4)
SpanishGBreally <- x2.2k(tbreally, 12, 4)
SwedishGBreally <- x2.2k(tbreally, 13, 4)
# completely
tbcompletely <- tb7 %>%
  filter(Variant == "completely")
tbcompletely$Language <- NULL
tbcompletely$Variant <- NULL
tbcompletely <- as.data.frame(tbcompletely)
tbcompletely$n <- as.numeric(tbcompletely$n)
tbcompletely$Amplifiers <- as.numeric(tbcompletely$Amplifiers)
# perform tests for completely
BulgarianGBcompletely <- x2.2k(tbcompletely, 1, 4)
CzechGBcompletely <- x2.2k(tbcompletely, 2, 4)
DutchGBcompletely <- x2.2k(tbcompletely, 3, 4)
FinnishGBcompletely <- x2.2k(tbcompletely, 5, 4)
FlemishGBcompletely <- x2.2k(tbcompletely, 6, 4)
FrenchGBcompletely <- x2.2k(tbcompletely, 7, 4)
GermanGBcompletely <- x2.2k(tbcompletely, 8, 4)
ItalianGBcompletely <- x2.2k(tbcompletely, 9, 4)
PolishGBcompletely <- x2.2k(tbcompletely, 10, 4)
RussianGBcompletely <- x2.2k(tbcompletely, 11, 4)
SpanishGBcompletely <- x2.2k(tbcompletely, 12, 4)
SwedishGBcompletely <- x2.2k(tbcompletely, 13, 4)
# other
tbother <- tb7 %>%
  filter(Variant == "other")
tbother$Language <- NULL
tbother$Variant <- NULL
tbother <- as.data.frame(tbother)
tbother$n <- as.numeric(tbother$n)
tbother$Amplifiers <- as.numeric(tbother$Amplifiers)
# perform tests for other
BulgarianGBother <- x2.2k(tbother, 1, 4)
CzechGBother <- x2.2k(tbother, 2, 4)
DutchGBother <- x2.2k(tbother, 3, 4)
FinnishGBother <- x2.2k(tbother, 5, 4)
FlemishGBother <- x2.2k(tbother, 6, 4)
FrenchGBother <- x2.2k(tbother, 7, 4)
GermanGBother <- x2.2k(tbother, 8, 4)
ItalianGBother <- x2.2k(tbother, 9, 4)
PolishGBother <- x2.2k(tbother, 10, 4)
RussianGBother <- x2.2k(tbother, 11, 4)
SpanishGBother <- x2.2k(tbother, 12, 4)
SwedishGBother <- x2.2k(tbother, 13, 4)
###############################################################
# function to extract values from list of tables
x2ktb <- function(x){
  x2kdf <- sapply(x, function(y){
    Description <- y$`Description`
    X2 <- y$`Chi-Squared`
    DF <- y$df
    P <- y$`p-value`
    Phi <- y$Phi
    Report <- y$Report
    x2kdf <- matrix(c(Description, X2, DF, P, Phi, Report), nrow= 1, byrow = T)
  })
  x2kdf <- t(x2kdf)
  colnames(x2kdf) <- c("Description", "X2", "DF", "P", "Phi", "Report")
  x2kdf <- as.data.frame(x2kdf)
  x2kdf$X2 <- as.numeric(x2kdf$X2)
  x2kdf$DF <- as.numeric(x2kdf$DF)
  x2kdf$P <- as.numeric(x2kdf$P)
  x2kdf$Phi <- as.numeric(x2kdf$Phi)
  return(x2kdf)
}
# create list of x2 tables
x2list <- list(BulgarianGBvery, CzechGBvery, DutchGBvery, FinnishGBvery, 
            FlemishGBvery, FrenchGBvery, GermanGBvery, ItalianGBvery, 
            PolishGBvery, RussianGBvery, SpanishGBvery, SwedishGBvery, 
            BulgarianGBso, CzechGBso, DutchGBso, FinnishGBso, FlemishGBso, 
            FrenchGBso, GermanGBso, ItalianGBso, PolishGBso, RussianGBso, 
            SpanishGBso, SwedishGBso, BulgarianGBreally, CzechGBreally, 
            DutchGBreally, FinnishGBreally, FlemishGBreally, FrenchGBreally, 
            GermanGBreally, ItalianGBreally, PolishGBreally, RussianGBreally, 
            SpanishGBreally, SwedishGBreally, BulgarianGBcompletely, CzechGBcompletely, 
            DutchGBcompletely, FinnishGBcompletely, FlemishGBcompletely, 
            FrenchGBcompletely, GermanGBcompletely, ItalianGBcompletely, 
            PolishGBcompletely, RussianGBcompletely, SpanishGBcompletely, 
            SwedishGBcompletely, BulgarianGBother, CzechGBother, DutchGBother, 
            FinnishGBother, FlemishGBother, FrenchGBother, GermanGBother, 
            ItalianGBother, PolishGBother, RussianGBother, SpanishGBother, SwedishGBother 
)
x2results <- x2ktb(x2list)
# add names
x2results$Language <- rep(c("Bulgarian", "Czech", "Dutch", "Finnish", "Flemish", "French", 
                        "German", "Italian", "Polish", "Russian", "Spanish", "Swedish"), 5)
x2results$Variant <- rep(c("very", "really", "so", "completely", "other"), each = 12)
x2results$Bonferroni <- rep(.05/nrow(x2results))
# extract sig. results
x2resultssig <- x2results %>% 
  filter(P < .05)
# inspect results
x2resultssig

# extract bonferroni corrected results
x2resultscorr <- x2resultssig %>% 
  filter(P < Bonferroni)
# inspect results
x2resultscorr

###############################################################
# tabulate x2 results for viz
tb8 <- x2results %>%
  select(X2, DF, P, Phi, Language, Variant, Bonferroni) %>%
  arrange(Language, Variant)
tb9 <- left_join(tb7, tb8, by = c("Language", "Variant"))
tb9 <- tb9  %>%
  mutate(
    Ratio = n/Amplifiers
  )
English_completely <- tb9 %>%
  filter(Language == "English")%>%
  filter(Variant == "completely")
English_other <- tb9 %>%
  filter(Language == "English")%>%
  filter(Variant == "other")
English_really <- tb9 %>%
  filter(Language == "English")%>%
  filter(Variant == "really")
English_so <- tb9 %>%
  filter(Language == "English")%>%
  filter(Variant == "so")
English_very <- tb9 %>%
  filter(Language == "English")%>%
  filter(Variant == "very")
completely <- tb9 %>%
  filter(Language != "English")%>%
  filter(Variant == "completely")
other <- tb9 %>%
  filter(Language != "English")%>%
  filter(Variant == "other")
really <- tb9 %>%
  filter(Language != "English")%>%
  filter(Variant == "really")
so <- tb9 %>%
  filter(Language != "English")%>%
  filter(Variant == "so")
very <- tb9 %>%
  filter(Language != "English")%>%
  filter(Variant == "very")
completely_comparative <- completely$Ratio-English_completely$Ratio
other_comparative <- other$Ratio-English_other$Ratio
really_comparative <- really$Ratio-English_really$Ratio
so_comparative <- so$Ratio-English_so$Ratio
very_comparative <- very$Ratio-English_very$Ratio
Language <- rep(very$Language, 5)
Difference <- c(completely_comparative, other_comparative,
                  really_comparative, so_comparative, very_comparative)
Variant <- c(rep(c("completely", "other", "really", "so", "very"), each = 12))
tb10 <- data.frame(Language, Difference, Variant)
tb10$Colour <- ifelse(tb10$Difference > 0, "grey80", "indianred4")

p4 <- ggplot(tb10, aes(Language, Difference, fill = Colour)) +                            
  facet_grid(vars(Variant), scales = "free_y") +
  geom_bar(stat = "identity")+#, aes(fill = Colour)) +                        
  theme_bw() +                                                                
  guides(fill=FALSE) +                                                        
  labs(x = "Language", y = "Difference to L1 English") + 
  scale_fill_manual(values=c("grey80", "indianred4"))                                                   
imageFile <- paste(imageDirectory,"DifferenceGB.png",sep="/")
ggsave(file = imageFile)
p4

###############################################################
#              SEMANTIC VECTOR SPACE MODEL 2
library(tidyr)
# evaluation how strongly really and very correlate
svsmd <- ampicle %>%
  select(Language, Adjective, Variant) %>%
  filter(Variant != "0") %>%
  group_by(Language) %>%
  count(Adjective, Variant) %>%
  spread(Language, n, fill = 0, convert = FALSE) %>%
  mutate(VariantAdjective = paste(Variant, Adjective, sep = ":")) %>%
  select(VariantAdjective, everything())
svsmd <- as.data.frame(svsmd)
rownames(svsmd) <- svsmd$VariantAdjective
svsmd <- svsmd %>%
  select(-VariantAdjective, -Variant, -Adjective) 
svsmd <- as.data.frame(svsmd)
svsmd <- svsmd[which(rowSums(svsmd) > 0),]
svsmd <- as.matrix(t(svsmd))
# compute expected values
svsmd.exp <- chisq.test(svsmd)
svsmd.exp <- svsmd.exp$expected
# calculate PMI and PPMI
svsmd.PMI <- log2(svsmd/svsmd.exp)
svsmd.PPMI <- ifelse(svsmd.PMI < 0, 0, svsmd.PMI)
# calculate cosine similarity
svsmd.tmp1 <- svsmd.PPMI
svsmd.cos <- cossim(svsmd.tmp1)
#round(svsmd.cos, 2)
###############################################################
#               CLUSTER SEMANTIC VECTORS
# load library
library(cluster)
# find max value that is not 1
svsmd.cos.test <- apply(svsmd.cos, 1, function(x){
  x <- ifelse(x == 1, 0, x) } )
maxval <- max(svsmd.cos.test)
# create distance matrix
svsmd.dist <- 1 - (svsmd.cos/maxval)
clustd <- as.dist(svsmd.dist)
# create distance matrix
clustd <- dist(svsmd.cos, method = "canberra") 
# alternative methods
# eucledian - not good when dealing with many dimensions
# manhattan - most popular choice
# method - here the difference between points dominates
# canberra - for count data
# binary - for binary data only!
# minkowski - is not a true distance measure

# find optimal number of clusters
asw <- as.vector(unlist(sapply(2:12, function(x) pam(clustd, k = x)$silinfo$avg.width)))
# determine the optimal number of clusters (max width is optimal)
optclust <- which(asw == max(asw))+1 # optimal number of clusters

# inspect clustering with optimal number of clusters
svsmd.clust <- pam(clustd, optclust)
svsmd.clust$clustering

# create cluster object
# alternative methods: "single", "ward.D2", "averAge", "mcquitty", "median", "centroid"
svsmdclust <- hclust(clustd, method="ward.D")    
# plot cluster solution
png("images/ClusteringSVSM.png",  width = 680, height = 480) # save plot
plot(svsmdclust, main = "", xlab = "", ylab = "", cex = .8)
rect.hclust(svsmdclust, k = optclust, border= "orange")
dev.off()
# Unrooted clustering
# load libraries for nicer dendrograms
library(factoextra)
library(dendextend)
library(ape)
# convert 'hclust' to 'phylo' object
phylo_tree = as.phylo(svsmdclust)
# get edges
graph_edges = phylo_tree$edge
# library igraph
library(igraph)
# get graph from edge list
graph_net = graph.edgelist(graph_edges)
# extract layout (x-y coords)
graph_layout = layout.auto(graph_net)
# number of observations
nobs = nrow(svsmd.cos)
# save plot
png("images/UnrootedClusteringSVSM.png",  width = 680, height = 480) 
# start plot
plot(graph_layout[,1], graph_layout[,2], type = "n", axes = FALSE,
     xlab = "", ylab = "")
# draw tree branches
segments(
  x0 = graph_layout[graph_edges[,1],1], 
  y0 = graph_layout[graph_edges[,1],2],
  x1 = graph_layout[graph_edges[,2],1],
  y1 = graph_layout[graph_edges[,2],2],
  col = "gray90", lwd = 2
)
# add labels
text(graph_layout[1:nobs,1], graph_layout[1:nobs,2],
     phylo_tree$tip.label, cex = .9, xpd = TRUE, font = 1)
dev.off()
###############################################################
# tabulate data
t1 <- tapply(ampicle$Amplified, list(ampicle$Adjective, ampicle$Variant), table)
t2 <- apply(t1, 1, function(x) ifelse(is.na(x) == T, 0, x))
t3 <- t(t2)
t3icle <- t3
#t3icle <- t3icle[, 2: ncol(t3icle)]
# remove Adjectives that were not Amplified
t3icle <- t3icle[rowSums(t3icle) > 0, ]
# save row and column names
colnamesicle <- colnames(t3icle)
rownamesicle <- rownames(t3icle)
# turn dataframe Amplifiedo matrix
svsmicle <- as.matrix(t3icle)
# convert token frequency to type frequency
#svsmicle <- apply(svsmicle, 1, Function(x) { x <- ifelse(x > 1, 1, x) } )
svsmicle <- t(svsmicle)
#svsmicle <- svsmicle[, colSums(svsmicle) >= 2]
#svsmicle <- svsmicle[rowSums(svsmicle) >= 2, ]
svsmicle

# determine overall n in data
n_icle <- sum(svsmicle)
n_icle

# correlate amplifiers based on collocation
r_icle <- cor(t(svsmicle))
r_icle

# extract correlation coefficient r for really and very
r_reallyveryicle <- r_icle[which(attr(r_icle, "dimnames")[[1]] == "very"), which(attr(r_icle, "dimnames")[[2]] == "really")]
r_reallyveryicle

# load required library
library(psych)
z_reallyveryicle <- fisherz(r_reallyveryicle)
z_reallyveryicle

# the z value can be tested for significance using the r-test from the psych library
#r.test(n=100,r12=.5,r34=.4, n2=80) 
###############################################################
#               PLOTTING LEXICAL DIVESRITY
# Function for extracting lexdiv values
lexdiv <- function(x){
  Varianttokentbicle <- table(x$Variant, x$Adjective)
  Varianttokentbicle <- Varianttokentbicle[2:nrow(Varianttokentbicle), ]
  #Varianttokentbicle <- Varianttokentbicle[rowSums(Varianttokentbicle) > 1, ]
  # extract typefrequency of tokenectives
  Varianttokentbicletyp <- t(apply(Varianttokentbicle, 1, function(x) ifelse(x > 1, 1, x)  ))
  # claculate lexical diversity measure
  lexdivicle <- rowSums(Varianttokentbicletyp)/rowSums(Varianttokentbicle)
  lexdivicle <- lexdivicle[order(lexdivicle)]
  return(lexdivicle)
}
# apply Function to data
lexdivicle <- lexdiv(tbd)

# ggplot2 p5
lexdivdf <- data.frame(1:length(lexdivicle), names(lexdivicle), round(lexdivicle, 2))
colnames(lexdivdf) <- c("id", "amp", "lexdiv")

# example extraction
x <- tbd
Varianttokentbicle <- table(x$Variant, x$Adjective)
Varianttokentbicle <- Varianttokentbicle[2:nrow(Varianttokentbicle), ]
Varianttokentbicletyp <- t(apply(Varianttokentbicle, 1, function(x) ifelse(x > 1, 1, x)  ))
lexdivicle <- rowSums(Varianttokentbicletyp)/rowSums(Varianttokentbicle)
###############################################################
#            PLOT LEXICAL DIVERSITY Across Language
p6d <- tbd %>%
  select(Language, Adjective, Variant)
vld <- c("completely", "really", "so", "very", "0")
p6d$Variant <- ifelse(p6d$Variant %in% vld, p6d$Variant, "other")
# extract tokenfrequency of Amplifiedensifiers
names(table(p6d$Language))

Bulgarian <- subset(p6d, Language == "Bulgarian")
Czech <- subset(p6d, Language == "Czech")
Dutch <- subset(p6d, Language == "Dutch")
English <- subset(p6d, Language == "English")
Finnish <- subset(p6d, Language == "Finnish")
Flemish <- subset(p6d, Language == "Flemish")
French <- subset(p6d, Language == "French")
German <- subset(p6d, Language == "German")
Italian <- subset(p6d, Language == "Italian")
Polish <- subset(p6d, Language == "Polish")
Russian <- subset(p6d, Language == "Russian")
Spanish <- subset(p6d, Language == "Spanish")
Swedish <- subset(p6d, Language == "Swedish")
# apply Function to data sets
lexdivBulgarian <- lexdiv(Bulgarian)
lexdivCzech <- lexdiv(Czech)
lexdivDutch <- lexdiv(Dutch)
lexdivEnglish <- lexdiv(English)
lexdivFinnish <- lexdiv(Finnish)
lexdivFlemish <- lexdiv(Flemish)
lexdivFrench <- lexdiv(French)
lexdivGerman <- lexdiv(German)
lexdivItalian <- lexdiv(Italian)
lexdivPolish <- lexdiv(Polish)
lexdivRussian <- lexdiv(Russian)
lexdivSpanish <- lexdiv(Spanish)
lexdivSwedish <- lexdiv(Swedish)
# find common items
cmnamps <- Reduce(intersect, list(names(lexdivBulgarian), names(lexdivCzech), 
                                  names(lexdivDutch), names(lexdivEnglish), 
                                  names(lexdivFinnish), names(lexdivFlemish), 
                                  names(lexdivFrench), names(lexdivGerman), 
                                  names(lexdivItalian), names(lexdivPolish), 
                                  names(lexdivRussian), names(lexdivSpanish), 
                                  names(lexdivSwedish)))
# extract lex div values for amps which occur in all Age groups
lexdivvls <- data.frame(lexdivBulgarian[which(names(lexdivBulgarian) %in% cmnamps)][order(names(lexdivBulgarian[which(names(lexdivBulgarian) %in% cmnamps)]))], 
                        lexdivCzech[which(names(lexdivCzech) %in% cmnamps)][order(names(lexdivCzech[which(names(lexdivCzech) %in% cmnamps)]))],
                        lexdivDutch[which(names(lexdivDutch) %in% cmnamps)][order(names(lexdivDutch[which(names(lexdivDutch) %in% cmnamps)]))],
                        lexdivEnglish[which(names(lexdivEnglish) %in% cmnamps)][order(names(lexdivEnglish[which(names(lexdivEnglish) %in% cmnamps)]))],
                        lexdivFinnish[which(names(lexdivFinnish) %in% cmnamps)][order(names(lexdivFinnish[which(names(lexdivFinnish) %in% cmnamps)]))],
                        lexdivFlemish[which(names(lexdivFlemish) %in% cmnamps)][order(names(lexdivFlemish[which(names(lexdivFlemish) %in% cmnamps)]))],
                        lexdivFrench[which(names(lexdivFrench) %in% cmnamps)][order(names(lexdivFrench[which(names(lexdivFrench) %in% cmnamps)]))],
                        lexdivGerman[which(names(lexdivGerman) %in% cmnamps)][order(names(lexdivGerman[which(names(lexdivGerman) %in% cmnamps)]))],
                        lexdivItalian[which(names(lexdivItalian) %in% cmnamps)][order(names(lexdivItalian[which(names(lexdivItalian) %in% cmnamps)]))],
                        lexdivPolish[which(names(lexdivPolish) %in% cmnamps)][order(names(lexdivPolish[which(names(lexdivPolish) %in% cmnamps)]))],
                        lexdivRussian[which(names(lexdivRussian) %in% cmnamps)][order(names(lexdivRussian[which(names(lexdivRussian) %in% cmnamps)]))],
                        lexdivSpanish[which(names(lexdivSpanish) %in% cmnamps)][order(names(lexdivSpanish[which(names(lexdivSpanish) %in% cmnamps)]))],
                        lexdivSwedish[which(names(lexdivSwedish) %in% cmnamps)][order(names(lexdivSwedish[which(names(lexdivSwedish) %in% cmnamps)]))])
# transpose data
lexdivvlst <- t(lexdivvls)
# combine lexdiv tables
p6d <- data.frame(1:nrow(lexdivvlst), names(table(ampicle$Language)), lexdivvlst)
colnames(p6d)[1:2] <- c("id", "Language")
rownames(p6d) <- 1:nrow(p6d)
p6d$Language <- factor(p6d$Language, levels=c("English", "Bulgarian", "Czech",
                                            "Dutch", "Finnish", "Flemish", "French", 
                                            "German", "Italian", "Polish", "Russian", 
                                            "Spanish", "Swedish"))
# start plot: Amplified
p6 <- ggplot(p6d, aes(x = Language, y = other, label = Language), size = 8) +
  geom_point(aes(y = completely, color = "completely"), size=1) +  
  geom_point(aes(y = other, color = "other"), size=1) +
  geom_point(aes(y = really, color = "really"), size=1) +
  geom_point(aes(y = very, color = "very"), size=1) +
  geom_point(aes(y = so, color = "so", lty = "so"), size=1) +
  guides(color=guide_legend(override.aes=list(fill=NA))) +
  scale_colour_manual(values=c("gray50", "goldenrod2", "gray70", "indianred4", "grey30"),
                      name="", 
                      breaks=c("completely", "other", "really",  "so",  "very"), 
                      labels = c("completely", "other", "really",  "so",  "very")) +
  theme(legend.position="top", axis.text.x = element_text(size=8, angle=90)) +
  theme_light(base_size = 10) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Language", y = "Lexical Diversity") +
  ggsave(file = paste(imageDirectory,"LexDivLanguage.png", sep="/"), width = 15,  height = 7.5, units = c("cm"),  dpi = 320)
p6

# recode variant 
tb11d <- tbd %>%
  select(Language, Adjective, Variant)
vld <- c("completely", "really", "so", "very", "0")
tb11d$Variant <- ifelse(tb11d$Variant %in% vld, tb11d$Variant, "other")
# tabulate
tb11 <- tb11d %>%
  group_by(Language, Variant) %>%
  summarise(
    VariantFrequency = n(),
    AdjectiveFrequency = length(names(table(Adjective))),
    LexicalDiversity = round(AdjectiveFrequency/VariantFrequency, 2)
  )
tb11

###############################################################
#       CORRESPONDENCE ANALYSIS
library("factoextra")        # activate packages for CA
#library("gplots")           # activate packages for CA
cad <- ampicle %>%
  select(Language, Amplified, Adjective, Variant) %>%
  filter(Variant != "0")
frqvntr <- names(table(cad$Variant))[which(table(cad$Variant) > 50)]
cad$Variant <- ifelse(cad$Variant %in% frqvntr, cad$Variant, "other")
frqadj <- names(table(cad$Adjective))[which(table(cad$Adjective) > 50)]
cad$Adjective <- ifelse(cad$Adjective %in% frqadj, cad$Adjective, "other")
# cerate language subsets
Bulgarian <- subset(cad, Language == "Bulgarian")
Czech <- subset(cad, Language == "Czech")
Dutch <- subset(cad, Language == "Dutch")
English <- subset(cad, Language == "English")
Finnish <- subset(cad, Language == "Finnish")
Flemish <- subset(cad, Language == "Flemish")
French <- subset(cad, Language == "French")
German <- subset(cad, Language == "German")
Italian <- subset(cad, Language == "Italian")
Polish <- subset(cad, Language == "Polish")
Russian <- subset(cad, Language == "Russian")
Spanish <- subset(cad, Language == "Spanish")
Swedish <- subset(cad, Language == "Swedish")
# English
CAEnglish <- table(English$Adjective, English$Variant)   
CAEnglish2 <- t(apply(CAEnglish, 1, function(x){ x <- ifelse(x > 1, 1, x)})) 
CAEnglish <- CAEnglish[, colSums(CAEnglish2) >= 2] # remove adj. amplified by only 1 variant
CAEnglish <- as.table(as.matrix(CAEnglish))        # convert data into table
library("FactoMineR")                              # activate packages for CA
res.ca <- CA(CAEnglish, graph = FALSE)             # create correspondence object
eig.val <- get_eigenvalue(res.ca)                  # extract eigencvalues
fviz_screeplot(res.ca, addlabels = TRUE, ylim = c(0, 60)) # % explained variances per dim.
fviz_contrib(res.ca, choice = "row", axes = 1, top = 10)  # Contributions to PC1
fviz_contrib(res.ca, choice = "col", axes = 2, top = 10)  # Contributions to PC2
# plot correspondence analysis results
png("images/CAFactorMapFvisAmpEnglish.png",  width = 680, height = 480) 
fviz_ca_biplot(res.ca, repel = T, col.ind = "cos2", col.row = "gray", col.var = "darkgray", title = "English")
dev.off()
# alternative plot
png("images/CAFactorMapAmpEnglish.png",  width = 680, height = 480) 
plot(res.ca, shadow = T, cex = 1, selectRow = "cos2 0.1", selectCol = "cos2 0.1", col.row = "gray50", title = "English")
dev.off()
###############################################################
# Bulgarian
CABulgarian <- table(Bulgarian$Adjective, Bulgarian$Variant)   
CABulgarian2 <- t(apply(CABulgarian, 1, function(x){ x <- ifelse(x > 1, 1, x)})) 
CABulgarian <- CABulgarian[, colSums(CABulgarian2) >= 2] 
CABulgarian <- as.table(as.matrix(CABulgarian))
library("FactoMineR")
res.ca <- CA(CABulgarian, graph = FALSE)
png("images/CAFactorMapFvisAmpBulgarian.png",  width = 680, height = 480) 
fviz_ca_biplot(res.ca, repel = T, col.ind = "cos2", col.row = "gray", col.var = "darkgray", title = "Bulgarian")
dev.off()
###############################################################
# Czech
CACzech <- table(Czech$Adjective, Czech$Variant)   
CACzech2 <- t(apply(CACzech, 1, function(x){ x <- ifelse(x > 1, 1, x)})) 
CACzech <- CACzech[, colSums(CACzech2) >= 2] 
CACzech <- as.table(as.matrix(CACzech))
library("FactoMineR")
res.ca <- CA(CACzech, graph = FALSE)
png("images/CAFactorMapFvisAmpCzech.png",  width = 680, height = 480) 
fviz_ca_biplot(res.ca, repel = T, col.ind = "cos2", col.row = "gray", col.var = "darkgray", title = "Czech")
dev.off()
###############################################################
# Dutch
CADutch <- table(Dutch$Adjective, Dutch$Variant)   
CADutch2 <- t(apply(CADutch, 1, function(x){ x <- ifelse(x > 1, 1, x)})) 
CADutch <- CADutch[, colSums(CADutch2) >= 2] 
CADutch <- as.table(as.matrix(CADutch))
library("FactoMineR")
res.ca <- CA(CADutch, graph = FALSE)
png("images/CAFactorMapFvisAmpDutch.png",  width = 680, height = 480) 
fviz_ca_biplot(res.ca, repel = T, col.ind = "cos2", col.row = "gray", col.var = "darkgray", title = "Dutch")
dev.off()
###############################################################
# Finnish
CAFinnish <- table(Finnish$Adjective, Finnish$Variant)   
CAFinnish2 <- t(apply(CAFinnish, 1, function(x){ x <- ifelse(x > 1, 1, x)})) 
CAFinnish <- CAFinnish[, colSums(CAFinnish2) >= 2] 
CAFinnish <- as.table(as.matrix(CAFinnish))
library("FactoMineR")
res.ca <- CA(CAFinnish, graph = FALSE)
png("images/CAFactorMapFvisAmpFinnish.png",  width = 680, height = 480) 
fviz_ca_biplot(res.ca, repel = T, col.ind = "cos2", col.row = "gray", col.var = "darkgray", title = "Finnish")
dev.off()
###############################################################
# Flemish
CAFlemish <- table(Flemish$Adjective, Flemish$Variant)   
CAFlemish2 <- t(apply(CAFlemish, 1, function(x){ x <- ifelse(x > 1, 1, x)})) 
CAFlemish <- CAFlemish[, colSums(CAFlemish2) >= 2] 
CAFlemish <- as.table(as.matrix(CAFlemish))
library("FactoMineR")
res.ca <- CA(CAFlemish, graph = FALSE)
png("images/CAFactorMapFvisAmpFlemish.png",  width = 680, height = 480) 
fviz_ca_biplot(res.ca, repel = T, col.ind = "cos2", col.row = "gray", col.var = "darkgray", title = "Flemish")
dev.off()
###############################################################
# French
CAFrench <- table(French$Adjective, French$Variant)   
CAFrench2 <- t(apply(CAFrench, 1, function(x){ x <- ifelse(x > 1, 1, x)})) 
CAFrench <- CAFrench[, colSums(CAFrench2) >= 2] 
CAFrench <- as.table(as.matrix(CAFrench))
library("FactoMineR")
res.ca <- CA(CAFrench, graph = FALSE)
png("images/CAFactorMapFvisAmpFrench.png",  width = 680, height = 480) 
fviz_ca_biplot(res.ca, repel = T, col.ind = "cos2", col.row = "gray", col.var = "darkgray", title = "French")
dev.off()
###############################################################
# German
CAGerman <- table(German$Adjective, German$Variant)   
CAGerman2 <- t(apply(CAGerman, 1, function(x){ x <- ifelse(x > 1, 1, x)})) 
CAGerman <- CAGerman[, colSums(CAGerman2) >= 2] 
CAGerman <- as.table(as.matrix(CAGerman))
library("FactoMineR")
res.ca <- CA(CAGerman, graph = FALSE)
png("images/CAFactorMapFvisAmpGerman.png",  width = 680, height = 480) 
fviz_ca_biplot(res.ca, repel = T, col.ind = "cos2", col.row = "gray", col.var = "darkgray", title = "German")
dev.off()
###############################################################
# Italian
CAItalian <- table(Italian$Adjective, Italian$Variant)   
CAItalian2 <- t(apply(CAItalian, 1, function(x){ x <- ifelse(x > 1, 1, x)})) 
CAItalian <- CAItalian[, colSums(CAItalian2) >= 2] 
CAItalian <- as.table(as.matrix(CAItalian))
library("FactoMineR")
res.ca <- CA(CAItalian, graph = FALSE)
png("images/CAFactorMapFvisAmpItalian.png",  width = 680, height = 480) 
fviz_ca_biplot(res.ca, repel = T, col.ind = "cos2", col.row = "gray", col.var = "darkgray", title = "Italian")
dev.off()
###############################################################
# Polish
CAPolish <- table(Polish$Adjective, Polish$Variant)   
CAPolish2 <- t(apply(CAPolish, 1, function(x){ x <- ifelse(x > 1, 1, x)})) 
CAPolish <- CAPolish[, colSums(CAPolish2) >= 2] 
CAPolish <- as.table(as.matrix(CAPolish))
library("FactoMineR")
res.ca <- CA(CAPolish, graph = FALSE)
png("images/CAFactorMapFvisAmpPolish.png",  width = 680, height = 480) 
fviz_ca_biplot(res.ca, repel = T, col.ind = "cos2", col.row = "gray", col.var = "darkgray", title = "Polish")
dev.off()
###############################################################
# Russian
CARussian <- table(Russian$Adjective, Russian$Variant)   
CARussian2 <- t(apply(CARussian, 1, function(x){ x <- ifelse(x > 1, 1, x)})) 
CARussian <- CARussian[, colSums(CARussian2) >= 2] 
CARussian <- as.table(as.matrix(CARussian))
library("FactoMineR")
res.ca <- CA(CARussian, graph = FALSE)
png("images/CAFactorMapFvisAmpRussian.png",  width = 680, height = 480) 
fviz_ca_biplot(res.ca, repel = T, col.ind = "cos2", col.row = "gray", col.var = "darkgray", title = "Russian")
dev.off()
###############################################################
# Spanish
CASpanish <- table(Spanish$Adjective, Spanish$Variant)   
CASpanish2 <- t(apply(CASpanish, 1, function(x){ x <- ifelse(x > 1, 1, x)})) 
CASpanish <- CASpanish[, colSums(CASpanish2) >= 2] 
CASpanish <- as.table(as.matrix(CASpanish))
library("FactoMineR")
res.ca <- CA(CASpanish, graph = FALSE)
png("images/CAFactorMapFvisAmpSpanish.png",  width = 680, height = 480) 
fviz_ca_biplot(res.ca, repel = T, col.ind = "cos2", col.row = "gray", col.var = "darkgray", title = "Spanish")
dev.off()
###############################################################
# Swedish
CASwedish <- table(Swedish$Adjective, Swedish$Variant)   
CASwedish2 <- t(apply(CASwedish, 1, function(x){ x <- ifelse(x > 1, 1, x)})) 
CASwedish <- CASwedish[, colSums(CASwedish2) >= 2] 
CASwedish <- as.table(as.matrix(CASwedish))
library("FactoMineR")
res.ca <- CA(CASwedish, graph = FALSE)
png("images/CAFactorMapFvisAmpSwedish.png",  width = 680, height = 480) 
fviz_ca_biplot(res.ca, repel = T, col.ind = "cos2", col.row = "gray", col.var = "darkgray", title = "Swedish")
dev.off()
###############################################################
#           COVARYING COLLEXEME ANALYSIS
# collex Function
collex <- function(data = data, cv1 = cv1){
  # set up rslttb
  rslttb <- matrix(c("amp", "Adjective", "namp", "nAdjective", "obs", "exp", "prob", "cs", "or", "p"), ncol = 10)
  colnames(rslttb) <- c("Amp", "Adjective", "N(Amp)", "N(Adjective)", "OBS", "EXP", 
                        "Probability", "CollStrength", "OddsRatio", "p")
  rvs <- 1:nrow(data)
  # define column values
  cv0 <- 1
  # set up table
  sapply(rvs, function(x){
    # extract values
    obs <- data[x,cv1] # freq Adjective with amp
    fAdjective <- sum(data[x,]) # freq Adjective
    n <- sum(data[,cv1]) # freq amp
    fall <- sum(data) # freq amps and Adjectives
    # calculate exp
    exp <- (fAdjective*n)/fall
    prob <- exp/n
    # create table to extract odds ratios
    m <- matrix(c(obs, (n-obs), (fAdjective-obs), (fall-fAdjective)), nrow = 2, byrow = T)
    o <- fisher.test(m)
    # perform binomial test
    rslt <- binom.test(obs, n, prob)
    # set up table with results
    rslttb <- list(c(colnames(data)[cv1], 
                     rownames(data)[x], 
                     n, 
                     fAdjective,
                     obs,
                     round(exp, 1),
                     round(prob, 2),
                     round(abs(log(as.vector(unlist(rslt$p.value, 10)), 10)), 2), 
                     round(as.vector(unlist(o$estimate)), 2),
                     round(as.vector(unlist(rslt$p.value)), 6)
    ))
    # return results
    return(rslttb)
  } )
}
###############################################################
#                 CVCLA
# determine which amplifiers are used in all languages
tst <- ampicle[ampicle$Amplified == 1,]
tst2 <- table(tst$Language, tst$Variant)
frqvrnt <- colnames(tst2)[which(colSums(tst2) >50)]
tst$Variant <- ifelse(tst$Variant %in% frqvrnt, tst$Variant, "other")
tst2 <- table(tst$Language, tst$Variant)
tst3 <- apply(tst2, 1, function(x){x <- ifelse(x > 1, 1, x)})
Amplifiers <- rownames(tst3)[which(rowSums(tst3) == 13)]
###############################################################
#                 CVCLA English
Englishtb1 <- tapply(English$Amplified, list(English$Adjective, English$Variant), table)
Englishtb2 <- t(apply(Englishtb1, 1, function(x) ifelse(is.na(x) == T, 0, x)))
Englishtb3 <- Englishtb2[rowSums(Englishtb2) > 0, ]
# apply collex function
absolutely  <- collex(data = Englishtb3, cv1 = which(colnames(Englishtb3) == Amplifiers[1]))
completely  <- collex(data = Englishtb3, cv1 = which(colnames(Englishtb3) == Amplifiers[2]))
extremely  <- collex(data = Englishtb3, cv1 = which(colnames(Englishtb3) == Amplifiers[3]))
highly  <- collex(data = Englishtb3, cv1 = which(colnames(Englishtb3) == Amplifiers[4]))
other  <- collex(data = Englishtb3, cv1 = which(colnames(Englishtb3) == Amplifiers[5]))
particularly  <- collex(data = Englishtb3, cv1 = which(colnames(Englishtb3) == Amplifiers[6]))
perfectly  <- collex(data = Englishtb3, cv1 = which(colnames(Englishtb3) == Amplifiers[7]))
really  <- collex(data = Englishtb3, cv1 = which(colnames(Englishtb3) == Amplifiers[8]))
so  <- collex(data = Englishtb3, cv1 = which(colnames(Englishtb3) == Amplifiers[9]))
totally  <- collex(data = Englishtb3, cv1 = which(colnames(Englishtb3) == Amplifiers[10]))
very  <- collex(data = Englishtb3, cv1 = which(colnames(Englishtb3) == Amplifiers[11]))
# extract informaltion
absolutely <- matrix(unlist(absolutely),ncol=10,byrow=TRUE)
completely <- matrix(unlist(completely),ncol=10,byrow=TRUE)
extremely <- matrix(unlist(extremely),ncol=10,byrow=TRUE)
highly <- matrix(unlist(highly),ncol=10,byrow=TRUE)
other <- matrix(unlist(other),ncol=10,byrow=TRUE)
particularly <- matrix(unlist(particularly),ncol=10,byrow=TRUE)
perfectly <- matrix(unlist(perfectly),ncol=10,byrow=TRUE)
really <- matrix(unlist(really),ncol=10,byrow=TRUE)
so <- matrix(unlist(so),ncol=10,byrow=TRUE)
totally <- matrix(unlist(totally),ncol=10,byrow=TRUE)
very <- matrix(unlist(very),ncol=10,byrow=TRUE)
# set up table with results
EnglishCollextb <- rbind(absolutely, completely, extremely, highly, other, 
                         particularly, perfectly, really, so, totally, very)
# write Function to process collextab (input = collextab)
collextbedit <- function(collextab){
  # convert Amplifiedo data frame
  collexdf <- as.data.frame(collextab)
  # add colnames
  colnames(collexdf) <- c("Amp", "Adjective", "N(Amp)", "N(Adjective)", "OBS", "EXP", 
                          "Probability", "CollStrength", "OddsRatio", "p")
  # add attraction column
  collexdf$attr <- ifelse(as.numeric(collexdf$OBS) > as.numeric(collexdf$EXP), "attr", "repel")
  # modify CollStrength column 
  collexdf$CollStrength <- ifelse(collexdf$attr == "repel", 
                                  paste("-", collexdf$CollStrength, sep =""), collexdf$CollStrength)
  # perform bonferroni correction
  corr05 <- 0.05/nrow(collexdf)
  collexdf$corr05 <- rep(corr05, nrow(collexdf))
  corr01 <- 0.01/nrow(collexdf)
  collexdf$corr01 <- rep(corr01, nrow(collexdf))
  corr001 <- 0.001/nrow(collexdf)
  collexdf$corr001 <- rep(corr001, nrow(collexdf))
  # calculate corrected significance status
  collexdf$sig <- as.vector(unlist(sapply(collexdf$p, function(x){
    x <- ifelse(x <= corr001, "p<.001",
                ifelse(x <= corr01, "p<.01",
                       ifelse(x <= corr001, "p<.001", "n.s."))) } )))
  return(collexdf)
}
# apply collextbedit Function
EnglishCollex <- collextbedit(EnglishCollextb)
###############################################################
#                 CVCLA Bulgarian
Bulgariantb1 <- tapply(Bulgarian$Amplified, list(Bulgarian$Adjective, Bulgarian$Variant), table)
Bulgariantb2 <- t(apply(Bulgariantb1, 1, function(x) ifelse(is.na(x) == T, 0, x)))
Bulgariantb3 <- Bulgariantb2[rowSums(Bulgariantb2) > 0, ]
# apply collex function
absolutely  <- collex(data = Bulgariantb3, cv1 = which(colnames(Bulgariantb3) == Amplifiers[1]))
completely  <- collex(data = Bulgariantb3, cv1 = which(colnames(Bulgariantb3) == Amplifiers[2]))
extremely  <- collex(data = Bulgariantb3, cv1 = which(colnames(Bulgariantb3) == Amplifiers[3]))
highly  <- collex(data = Bulgariantb3, cv1 = which(colnames(Bulgariantb3) == Amplifiers[4]))
other  <- collex(data = Bulgariantb3, cv1 = which(colnames(Bulgariantb3) == Amplifiers[5]))
particularly  <- collex(data = Bulgariantb3, cv1 = which(colnames(Bulgariantb3) == Amplifiers[6]))
perfectly  <- collex(data = Bulgariantb3, cv1 = which(colnames(Bulgariantb3) == Amplifiers[7]))
really  <- collex(data = Bulgariantb3, cv1 = which(colnames(Bulgariantb3) == Amplifiers[8]))
so  <- collex(data = Bulgariantb3, cv1 = which(colnames(Bulgariantb3) == Amplifiers[9]))
totally  <- collex(data = Bulgariantb3, cv1 = which(colnames(Bulgariantb3) == Amplifiers[10]))
very  <- collex(data = Bulgariantb3, cv1 = which(colnames(Bulgariantb3) == Amplifiers[11]))
# extract informaltion
absolutely <- matrix(unlist(absolutely),ncol=10,byrow=TRUE)
completely <- matrix(unlist(completely),ncol=10,byrow=TRUE)
extremely <- matrix(unlist(extremely),ncol=10,byrow=TRUE)
highly <- matrix(unlist(highly),ncol=10,byrow=TRUE)
other <- matrix(unlist(other),ncol=10,byrow=TRUE)
particularly <- matrix(unlist(particularly),ncol=10,byrow=TRUE)
perfectly <- matrix(unlist(perfectly),ncol=10,byrow=TRUE)
really <- matrix(unlist(really),ncol=10,byrow=TRUE)
so <- matrix(unlist(so),ncol=10,byrow=TRUE)
totally <- matrix(unlist(totally),ncol=10,byrow=TRUE)
very <- matrix(unlist(very),ncol=10,byrow=TRUE)
# set up table with results
BulgarianCollextb <- rbind(absolutely, completely, extremely, highly, other, 
                           particularly, perfectly, really, so, totally, very)
# apply collextbedit Function
BulgarianCollex <- collextbedit(BulgarianCollextb)
###############################################################
#                 CVCLA Czech
Czechtb1 <- tapply(Czech$Amplified, list(Czech$Adjective, Czech$Variant), table)
Czechtb2 <- t(apply(Czechtb1, 1, function(x) ifelse(is.na(x) == T, 0, x)))
Czechtb3 <- Czechtb2[rowSums(Czechtb2) > 0, ]
# apply collex function
absolutely  <- collex(data = Czechtb3, cv1 = which(colnames(Czechtb3) == Amplifiers[1]))
completely  <- collex(data = Czechtb3, cv1 = which(colnames(Czechtb3) == Amplifiers[2]))
extremely  <- collex(data = Czechtb3, cv1 = which(colnames(Czechtb3) == Amplifiers[3]))
highly  <- collex(data = Czechtb3, cv1 = which(colnames(Czechtb3) == Amplifiers[4]))
other  <- collex(data = Czechtb3, cv1 = which(colnames(Czechtb3) == Amplifiers[5]))
particularly  <- collex(data = Czechtb3, cv1 = which(colnames(Czechtb3) == Amplifiers[6]))
perfectly  <- collex(data = Czechtb3, cv1 = which(colnames(Czechtb3) == Amplifiers[7]))
really  <- collex(data = Czechtb3, cv1 = which(colnames(Czechtb3) == Amplifiers[8]))
so  <- collex(data = Czechtb3, cv1 = which(colnames(Czechtb3) == Amplifiers[9]))
totally  <- collex(data = Czechtb3, cv1 = which(colnames(Czechtb3) == Amplifiers[10]))
very  <- collex(data = Czechtb3, cv1 = which(colnames(Czechtb3) == Amplifiers[11]))
# extract informaltion
absolutely <- matrix(unlist(absolutely),ncol=10,byrow=TRUE)
completely <- matrix(unlist(completely),ncol=10,byrow=TRUE)
extremely <- matrix(unlist(extremely),ncol=10,byrow=TRUE)
highly <- matrix(unlist(highly),ncol=10,byrow=TRUE)
other <- matrix(unlist(other),ncol=10,byrow=TRUE)
particularly <- matrix(unlist(particularly),ncol=10,byrow=TRUE)
perfectly <- matrix(unlist(perfectly),ncol=10,byrow=TRUE)
really <- matrix(unlist(really),ncol=10,byrow=TRUE)
so <- matrix(unlist(so),ncol=10,byrow=TRUE)
totally <- matrix(unlist(totally),ncol=10,byrow=TRUE)
very <- matrix(unlist(very),ncol=10,byrow=TRUE)
# set up table with results
CzechCollextb <- rbind(absolutely, completely, extremely, highly, other, 
                       particularly, perfectly, really, so, totally, very)
# apply collextbedit Function
CzechCollex <- collextbedit(CzechCollextb)
###############################################################
#                 CVCLA Dutch
Dutchtb1 <- tapply(Dutch$Amplified, list(Dutch$Adjective, Dutch$Variant), table)
Dutchtb2 <- t(apply(Dutchtb1, 1, function(x) ifelse(is.na(x) == T, 0, x)))
Dutchtb3 <- Dutchtb2[rowSums(Dutchtb2) > 0, ]
# apply collex function
absolutely  <- collex(data = Dutchtb3, cv1 = which(colnames(Dutchtb3) == Amplifiers[1]))
completely  <- collex(data = Dutchtb3, cv1 = which(colnames(Dutchtb3) == Amplifiers[2]))
extremely  <- collex(data = Dutchtb3, cv1 = which(colnames(Dutchtb3) == Amplifiers[3]))
highly  <- collex(data = Dutchtb3, cv1 = which(colnames(Dutchtb3) == Amplifiers[4]))
other  <- collex(data = Dutchtb3, cv1 = which(colnames(Dutchtb3) == Amplifiers[5]))
particularly  <- collex(data = Dutchtb3, cv1 = which(colnames(Dutchtb3) == Amplifiers[6]))
perfectly  <- collex(data = Dutchtb3, cv1 = which(colnames(Dutchtb3) == Amplifiers[7]))
really  <- collex(data = Dutchtb3, cv1 = which(colnames(Dutchtb3) == Amplifiers[8]))
so  <- collex(data = Dutchtb3, cv1 = which(colnames(Dutchtb3) == Amplifiers[9]))
totally  <- collex(data = Dutchtb3, cv1 = which(colnames(Dutchtb3) == Amplifiers[10]))
very  <- collex(data = Dutchtb3, cv1 = which(colnames(Dutchtb3) == Amplifiers[11]))
# extract informaltion
absolutely <- matrix(unlist(absolutely),ncol=10,byrow=TRUE)
completely <- matrix(unlist(completely),ncol=10,byrow=TRUE)
extremely <- matrix(unlist(extremely),ncol=10,byrow=TRUE)
highly <- matrix(unlist(highly),ncol=10,byrow=TRUE)
other <- matrix(unlist(other),ncol=10,byrow=TRUE)
particularly <- matrix(unlist(particularly),ncol=10,byrow=TRUE)
perfectly <- matrix(unlist(perfectly),ncol=10,byrow=TRUE)
really <- matrix(unlist(really),ncol=10,byrow=TRUE)
so <- matrix(unlist(so),ncol=10,byrow=TRUE)
totally <- matrix(unlist(totally),ncol=10,byrow=TRUE)
very <- matrix(unlist(very),ncol=10,byrow=TRUE)
# set up table with results
DutchCollextb <- rbind(absolutely, completely, extremely, highly, other, 
                       particularly, perfectly, really, so, totally, very)
# apply collextbedit Function
DutchCollex <- collextbedit(DutchCollextb)
###############################################################
#                 CVCLA Finnish
Finnishtb1 <- tapply(Finnish$Amplified, list(Finnish$Adjective, Finnish$Variant), table)
Finnishtb2 <- t(apply(Finnishtb1, 1, function(x) ifelse(is.na(x) == T, 0, x)))
Finnishtb3 <- Finnishtb2[rowSums(Finnishtb2) > 0, ]
# apply collex function
absolutely  <- collex(data = Finnishtb3, cv1 = which(colnames(Finnishtb3) == Amplifiers[1]))
completely  <- collex(data = Finnishtb3, cv1 = which(colnames(Finnishtb3) == Amplifiers[2]))
extremely  <- collex(data = Finnishtb3, cv1 = which(colnames(Finnishtb3) == Amplifiers[3]))
highly  <- collex(data = Finnishtb3, cv1 = which(colnames(Finnishtb3) == Amplifiers[4]))
other  <- collex(data = Finnishtb3, cv1 = which(colnames(Finnishtb3) == Amplifiers[5]))
particularly  <- collex(data = Finnishtb3, cv1 = which(colnames(Finnishtb3) == Amplifiers[6]))
perfectly  <- collex(data = Finnishtb3, cv1 = which(colnames(Finnishtb3) == Amplifiers[7]))
really  <- collex(data = Finnishtb3, cv1 = which(colnames(Finnishtb3) == Amplifiers[8]))
so  <- collex(data = Finnishtb3, cv1 = which(colnames(Finnishtb3) == Amplifiers[9]))
totally  <- collex(data = Finnishtb3, cv1 = which(colnames(Finnishtb3) == Amplifiers[10]))
very  <- collex(data = Finnishtb3, cv1 = which(colnames(Finnishtb3) == Amplifiers[11]))
# extract informaltion
absolutely <- matrix(unlist(absolutely),ncol=10,byrow=TRUE)
completely <- matrix(unlist(completely),ncol=10,byrow=TRUE)
extremely <- matrix(unlist(extremely),ncol=10,byrow=TRUE)
highly <- matrix(unlist(highly),ncol=10,byrow=TRUE)
other <- matrix(unlist(other),ncol=10,byrow=TRUE)
particularly <- matrix(unlist(particularly),ncol=10,byrow=TRUE)
perfectly <- matrix(unlist(perfectly),ncol=10,byrow=TRUE)
really <- matrix(unlist(really),ncol=10,byrow=TRUE)
so <- matrix(unlist(so),ncol=10,byrow=TRUE)
totally <- matrix(unlist(totally),ncol=10,byrow=TRUE)
very <- matrix(unlist(very),ncol=10,byrow=TRUE)
# set up table with results
FinnishCollextb <- rbind(absolutely, completely, extremely, highly, other, 
                         particularly, perfectly, really, so, totally, very)
# apply collextbedit Function
FinnishCollex <- collextbedit(FinnishCollextb)
###############################################################
#                 CVCLA Flemish
Flemishtb1 <- tapply(Flemish$Amplified, list(Flemish$Adjective, Flemish$Variant), table)
Flemishtb2 <- t(apply(Flemishtb1, 1, function(x) ifelse(is.na(x) == T, 0, x)))
Flemishtb3 <- Flemishtb2[rowSums(Flemishtb2) > 0, ]
# apply collex function
absolutely  <- collex(data = Flemishtb3, cv1 = which(colnames(Flemishtb3) == Amplifiers[1]))
completely  <- collex(data = Flemishtb3, cv1 = which(colnames(Flemishtb3) == Amplifiers[2]))
extremely  <- collex(data = Flemishtb3, cv1 = which(colnames(Flemishtb3) == Amplifiers[3]))
highly  <- collex(data = Flemishtb3, cv1 = which(colnames(Flemishtb3) == Amplifiers[4]))
other  <- collex(data = Flemishtb3, cv1 = which(colnames(Flemishtb3) == Amplifiers[5]))
particularly  <- collex(data = Flemishtb3, cv1 = which(colnames(Flemishtb3) == Amplifiers[6]))
perfectly  <- collex(data = Flemishtb3, cv1 = which(colnames(Flemishtb3) == Amplifiers[7]))
really  <- collex(data = Flemishtb3, cv1 = which(colnames(Flemishtb3) == Amplifiers[8]))
so  <- collex(data = Flemishtb3, cv1 = which(colnames(Flemishtb3) == Amplifiers[9]))
totally  <- collex(data = Flemishtb3, cv1 = which(colnames(Flemishtb3) == Amplifiers[10]))
very  <- collex(data = Flemishtb3, cv1 = which(colnames(Flemishtb3) == Amplifiers[11]))
# extract informaltion
absolutely <- matrix(unlist(absolutely),ncol=10,byrow=TRUE)
completely <- matrix(unlist(completely),ncol=10,byrow=TRUE)
extremely <- matrix(unlist(extremely),ncol=10,byrow=TRUE)
highly <- matrix(unlist(highly),ncol=10,byrow=TRUE)
other <- matrix(unlist(other),ncol=10,byrow=TRUE)
particularly <- matrix(unlist(particularly),ncol=10,byrow=TRUE)
perfectly <- matrix(unlist(perfectly),ncol=10,byrow=TRUE)
really <- matrix(unlist(really),ncol=10,byrow=TRUE)
so <- matrix(unlist(so),ncol=10,byrow=TRUE)
totally <- matrix(unlist(totally),ncol=10,byrow=TRUE)
very <- matrix(unlist(very),ncol=10,byrow=TRUE)
# set up table with results
FlemishCollextb <- rbind(absolutely, completely, extremely, highly, other, 
                         particularly, perfectly, really, so, totally, very)
# apply collextbedit Function
FlemishCollex <- collextbedit(FlemishCollextb)
###############################################################
#                 CVCLA French
Frenchtb1 <- tapply(French$Amplified, list(French$Adjective, French$Variant), table)
Frenchtb2 <- t(apply(Frenchtb1, 1, function(x) ifelse(is.na(x) == T, 0, x)))
Frenchtb3 <- Frenchtb2[rowSums(Frenchtb2) > 0, ]
# apply collex function
absolutely  <- collex(data = Frenchtb3, cv1 = which(colnames(Frenchtb3) == Amplifiers[1]))
completely  <- collex(data = Frenchtb3, cv1 = which(colnames(Frenchtb3) == Amplifiers[2]))
extremely  <- collex(data = Frenchtb3, cv1 = which(colnames(Frenchtb3) == Amplifiers[3]))
highly  <- collex(data = Frenchtb3, cv1 = which(colnames(Frenchtb3) == Amplifiers[4]))
other  <- collex(data = Frenchtb3, cv1 = which(colnames(Frenchtb3) == Amplifiers[5]))
particularly  <- collex(data = Frenchtb3, cv1 = which(colnames(Frenchtb3) == Amplifiers[6]))
perfectly  <- collex(data = Frenchtb3, cv1 = which(colnames(Frenchtb3) == Amplifiers[7]))
really  <- collex(data = Frenchtb3, cv1 = which(colnames(Frenchtb3) == Amplifiers[8]))
so  <- collex(data = Frenchtb3, cv1 = which(colnames(Frenchtb3) == Amplifiers[9]))
totally  <- collex(data = Frenchtb3, cv1 = which(colnames(Frenchtb3) == Amplifiers[10]))
very  <- collex(data = Frenchtb3, cv1 = which(colnames(Frenchtb3) == Amplifiers[11]))
# extract informaltion
absolutely <- matrix(unlist(absolutely),ncol=10,byrow=TRUE)
completely <- matrix(unlist(completely),ncol=10,byrow=TRUE)
extremely <- matrix(unlist(extremely),ncol=10,byrow=TRUE)
highly <- matrix(unlist(highly),ncol=10,byrow=TRUE)
other <- matrix(unlist(other),ncol=10,byrow=TRUE)
particularly <- matrix(unlist(particularly),ncol=10,byrow=TRUE)
perfectly <- matrix(unlist(perfectly),ncol=10,byrow=TRUE)
really <- matrix(unlist(really),ncol=10,byrow=TRUE)
so <- matrix(unlist(so),ncol=10,byrow=TRUE)
totally <- matrix(unlist(totally),ncol=10,byrow=TRUE)
very <- matrix(unlist(very),ncol=10,byrow=TRUE)
# set up table with results
FrenchCollextb <- rbind(absolutely, completely, extremely, highly, other, 
                        particularly, perfectly, really, so, totally, very)
# apply collextbedit Function
FrenchCollex <- collextbedit(FrenchCollextb)
###############################################################
#                 CVCLA German
Germantb1 <- tapply(German$Amplified, list(German$Adjective, German$Variant), table)
Germantb2 <- t(apply(Germantb1, 1, function(x) ifelse(is.na(x) == T, 0, x)))
Germantb3 <- Germantb2[rowSums(Germantb2) > 0, ]
# apply collex function
absolutely  <- collex(data = Germantb3, cv1 = which(colnames(Germantb3) == Amplifiers[1]))
completely  <- collex(data = Germantb3, cv1 = which(colnames(Germantb3) == Amplifiers[2]))
extremely  <- collex(data = Germantb3, cv1 = which(colnames(Germantb3) == Amplifiers[3]))
highly  <- collex(data = Germantb3, cv1 = which(colnames(Germantb3) == Amplifiers[4]))
other  <- collex(data = Germantb3, cv1 = which(colnames(Germantb3) == Amplifiers[5]))
particularly  <- collex(data = Germantb3, cv1 = which(colnames(Germantb3) == Amplifiers[6]))
perfectly  <- collex(data = Germantb3, cv1 = which(colnames(Germantb3) == Amplifiers[7]))
really  <- collex(data = Germantb3, cv1 = which(colnames(Germantb3) == Amplifiers[8]))
so  <- collex(data = Germantb3, cv1 = which(colnames(Germantb3) == Amplifiers[9]))
totally  <- collex(data = Germantb3, cv1 = which(colnames(Germantb3) == Amplifiers[10]))
very  <- collex(data = Germantb3, cv1 = which(colnames(Germantb3) == Amplifiers[11]))
# extract informaltion
absolutely <- matrix(unlist(absolutely),ncol=10,byrow=TRUE)
completely <- matrix(unlist(completely),ncol=10,byrow=TRUE)
extremely <- matrix(unlist(extremely),ncol=10,byrow=TRUE)
highly <- matrix(unlist(highly),ncol=10,byrow=TRUE)
other <- matrix(unlist(other),ncol=10,byrow=TRUE)
particularly <- matrix(unlist(particularly),ncol=10,byrow=TRUE)
perfectly <- matrix(unlist(perfectly),ncol=10,byrow=TRUE)
really <- matrix(unlist(really),ncol=10,byrow=TRUE)
so <- matrix(unlist(so),ncol=10,byrow=TRUE)
totally <- matrix(unlist(totally),ncol=10,byrow=TRUE)
very <- matrix(unlist(very),ncol=10,byrow=TRUE)
# set up table with results
GermanCollextb <- rbind(absolutely, completely, extremely, highly, other, 
                        particularly, perfectly, really, so, totally, very)
# apply collextbedit Function
GermanCollex <- collextbedit(GermanCollextb)
###############################################################
#                 CVCLA Italian
Italiantb1 <- tapply(Italian$Amplified, list(Italian$Adjective, Italian$Variant), table)
Italiantb2 <- t(apply(Italiantb1, 1, function(x) ifelse(is.na(x) == T, 0, x)))
Italiantb3 <- Italiantb2[rowSums(Italiantb2) > 0, ]
# apply collex function
absolutely  <- collex(data = Italiantb3, cv1 = which(colnames(Italiantb3) == Amplifiers[1]))
completely  <- collex(data = Italiantb3, cv1 = which(colnames(Italiantb3) == Amplifiers[2]))
extremely  <- collex(data = Italiantb3, cv1 = which(colnames(Italiantb3) == Amplifiers[3]))
highly  <- collex(data = Italiantb3, cv1 = which(colnames(Italiantb3) == Amplifiers[4]))
other  <- collex(data = Italiantb3, cv1 = which(colnames(Italiantb3) == Amplifiers[5]))
particularly  <- collex(data = Italiantb3, cv1 = which(colnames(Italiantb3) == Amplifiers[6]))
perfectly  <- collex(data = Italiantb3, cv1 = which(colnames(Italiantb3) == Amplifiers[7]))
really  <- collex(data = Italiantb3, cv1 = which(colnames(Italiantb3) == Amplifiers[8]))
so  <- collex(data = Italiantb3, cv1 = which(colnames(Italiantb3) == Amplifiers[9]))
totally  <- collex(data = Italiantb3, cv1 = which(colnames(Italiantb3) == Amplifiers[10]))
very  <- collex(data = Italiantb3, cv1 = which(colnames(Italiantb3) == Amplifiers[11]))
# extract informaltion
absolutely <- matrix(unlist(absolutely),ncol=10,byrow=TRUE)
completely <- matrix(unlist(completely),ncol=10,byrow=TRUE)
extremely <- matrix(unlist(extremely),ncol=10,byrow=TRUE)
highly <- matrix(unlist(highly),ncol=10,byrow=TRUE)
other <- matrix(unlist(other),ncol=10,byrow=TRUE)
particularly <- matrix(unlist(particularly),ncol=10,byrow=TRUE)
perfectly <- matrix(unlist(perfectly),ncol=10,byrow=TRUE)
really <- matrix(unlist(really),ncol=10,byrow=TRUE)
so <- matrix(unlist(so),ncol=10,byrow=TRUE)
totally <- matrix(unlist(totally),ncol=10,byrow=TRUE)
very <- matrix(unlist(very),ncol=10,byrow=TRUE)
# set up table with results
ItalianCollextb <- rbind(absolutely, completely, extremely, highly, other, 
                         particularly, perfectly, really, so, totally, very)
# apply collextbedit Function
ItalianCollex <- collextbedit(ItalianCollextb)
###############################################################
#                 CVCLA Polish
Polishtb1 <- tapply(Polish$Amplified, list(Polish$Adjective, Polish$Variant), table)
Polishtb2 <- t(apply(Polishtb1, 1, function(x) ifelse(is.na(x) == T, 0, x)))
Polishtb3 <- Polishtb2[rowSums(Polishtb2) > 0, ]
# apply collex function
absolutely  <- collex(data = Polishtb3, cv1 = which(colnames(Polishtb3) == Amplifiers[1]))
completely  <- collex(data = Polishtb3, cv1 = which(colnames(Polishtb3) == Amplifiers[2]))
extremely  <- collex(data = Polishtb3, cv1 = which(colnames(Polishtb3) == Amplifiers[3]))
highly  <- collex(data = Polishtb3, cv1 = which(colnames(Polishtb3) == Amplifiers[4]))
other  <- collex(data = Polishtb3, cv1 = which(colnames(Polishtb3) == Amplifiers[5]))
particularly  <- collex(data = Polishtb3, cv1 = which(colnames(Polishtb3) == Amplifiers[6]))
perfectly  <- collex(data = Polishtb3, cv1 = which(colnames(Polishtb3) == Amplifiers[7]))
really  <- collex(data = Polishtb3, cv1 = which(colnames(Polishtb3) == Amplifiers[8]))
so  <- collex(data = Polishtb3, cv1 = which(colnames(Polishtb3) == Amplifiers[9]))
totally  <- collex(data = Polishtb3, cv1 = which(colnames(Polishtb3) == Amplifiers[10]))
very  <- collex(data = Polishtb3, cv1 = which(colnames(Polishtb3) == Amplifiers[11]))
# extract informaltion
absolutely <- matrix(unlist(absolutely),ncol=10,byrow=TRUE)
completely <- matrix(unlist(completely),ncol=10,byrow=TRUE)
extremely <- matrix(unlist(extremely),ncol=10,byrow=TRUE)
highly <- matrix(unlist(highly),ncol=10,byrow=TRUE)
other <- matrix(unlist(other),ncol=10,byrow=TRUE)
particularly <- matrix(unlist(particularly),ncol=10,byrow=TRUE)
perfectly <- matrix(unlist(perfectly),ncol=10,byrow=TRUE)
really <- matrix(unlist(really),ncol=10,byrow=TRUE)
so <- matrix(unlist(so),ncol=10,byrow=TRUE)
totally <- matrix(unlist(totally),ncol=10,byrow=TRUE)
very <- matrix(unlist(very),ncol=10,byrow=TRUE)
# set up table with results
PolishCollextb <- rbind(absolutely, completely, extremely, highly, other, 
                        particularly, perfectly, really, so, totally, very)
# apply collextbedit Function
PolishCollex <- collextbedit(PolishCollextb)
###############################################################
#                 CVCLA Russian
Russiantb1 <- tapply(Russian$Amplified, list(Russian$Adjective, Russian$Variant), table)
Russiantb2 <- t(apply(Russiantb1, 1, function(x) ifelse(is.na(x) == T, 0, x)))
Russiantb3 <- Russiantb2[rowSums(Russiantb2) > 0, ]
# apply collex function
absolutely  <- collex(data = Russiantb3, cv1 = which(colnames(Russiantb3) == Amplifiers[1]))
completely  <- collex(data = Russiantb3, cv1 = which(colnames(Russiantb3) == Amplifiers[2]))
extremely  <- collex(data = Russiantb3, cv1 = which(colnames(Russiantb3) == Amplifiers[3]))
highly  <- collex(data = Russiantb3, cv1 = which(colnames(Russiantb3) == Amplifiers[4]))
other  <- collex(data = Russiantb3, cv1 = which(colnames(Russiantb3) == Amplifiers[5]))
particularly  <- collex(data = Russiantb3, cv1 = which(colnames(Russiantb3) == Amplifiers[6]))
perfectly  <- collex(data = Russiantb3, cv1 = which(colnames(Russiantb3) == Amplifiers[7]))
really  <- collex(data = Russiantb3, cv1 = which(colnames(Russiantb3) == Amplifiers[8]))
so  <- collex(data = Russiantb3, cv1 = which(colnames(Russiantb3) == Amplifiers[9]))
totally  <- collex(data = Russiantb3, cv1 = which(colnames(Russiantb3) == Amplifiers[10]))
very  <- collex(data = Russiantb3, cv1 = which(colnames(Russiantb3) == Amplifiers[11]))
# extract informaltion
absolutely <- matrix(unlist(absolutely),ncol=10,byrow=TRUE)
completely <- matrix(unlist(completely),ncol=10,byrow=TRUE)
extremely <- matrix(unlist(extremely),ncol=10,byrow=TRUE)
highly <- matrix(unlist(highly),ncol=10,byrow=TRUE)
other <- matrix(unlist(other),ncol=10,byrow=TRUE)
particularly <- matrix(unlist(particularly),ncol=10,byrow=TRUE)
perfectly <- matrix(unlist(perfectly),ncol=10,byrow=TRUE)
really <- matrix(unlist(really),ncol=10,byrow=TRUE)
so <- matrix(unlist(so),ncol=10,byrow=TRUE)
totally <- matrix(unlist(totally),ncol=10,byrow=TRUE)
very <- matrix(unlist(very),ncol=10,byrow=TRUE)
# set up table with results
RussianCollextb <- rbind(absolutely, completely, extremely, highly, other, 
                         particularly, perfectly, really, so, totally, very)
# apply collextbedit Function
RussianCollex <- collextbedit(RussianCollextb)
###############################################################
#                 CVCLA Spanish
Spanishtb1 <- tapply(Spanish$Amplified, list(Spanish$Adjective, Spanish$Variant), table)
Spanishtb2 <- t(apply(Spanishtb1, 1, function(x) ifelse(is.na(x) == T, 0, x)))
Spanishtb3 <- Spanishtb2[rowSums(Spanishtb2) > 0, ]
# apply collex function
absolutely  <- collex(data = Spanishtb3, cv1 = which(colnames(Spanishtb3) == Amplifiers[1]))
completely  <- collex(data = Spanishtb3, cv1 = which(colnames(Spanishtb3) == Amplifiers[2]))
extremely  <- collex(data = Spanishtb3, cv1 = which(colnames(Spanishtb3) == Amplifiers[3]))
highly  <- collex(data = Spanishtb3, cv1 = which(colnames(Spanishtb3) == Amplifiers[4]))
other  <- collex(data = Spanishtb3, cv1 = which(colnames(Spanishtb3) == Amplifiers[5]))
particularly  <- collex(data = Spanishtb3, cv1 = which(colnames(Spanishtb3) == Amplifiers[6]))
perfectly  <- collex(data = Spanishtb3, cv1 = which(colnames(Spanishtb3) == Amplifiers[7]))
really  <- collex(data = Spanishtb3, cv1 = which(colnames(Spanishtb3) == Amplifiers[8]))
so  <- collex(data = Spanishtb3, cv1 = which(colnames(Spanishtb3) == Amplifiers[9]))
totally  <- collex(data = Spanishtb3, cv1 = which(colnames(Spanishtb3) == Amplifiers[10]))
very  <- collex(data = Spanishtb3, cv1 = which(colnames(Spanishtb3) == Amplifiers[11]))
# extract informaltion
absolutely <- matrix(unlist(absolutely),ncol=10,byrow=TRUE)
completely <- matrix(unlist(completely),ncol=10,byrow=TRUE)
extremely <- matrix(unlist(extremely),ncol=10,byrow=TRUE)
highly <- matrix(unlist(highly),ncol=10,byrow=TRUE)
other <- matrix(unlist(other),ncol=10,byrow=TRUE)
particularly <- matrix(unlist(particularly),ncol=10,byrow=TRUE)
perfectly <- matrix(unlist(perfectly),ncol=10,byrow=TRUE)
really <- matrix(unlist(really),ncol=10,byrow=TRUE)
so <- matrix(unlist(so),ncol=10,byrow=TRUE)
totally <- matrix(unlist(totally),ncol=10,byrow=TRUE)
very <- matrix(unlist(very),ncol=10,byrow=TRUE)
# set up table with results
SpanishCollextb <- rbind(absolutely, completely, extremely, highly, other, 
                         particularly, perfectly, really, so, totally, very)
# apply collextbedit Function
SpanishCollex <- collextbedit(SpanishCollextb)
###############################################################
#                 CVCLA Swedish
Swedishtb1 <- tapply(Swedish$Amplified, list(Swedish$Adjective, Swedish$Variant), table)
Swedishtb2 <- t(apply(Swedishtb1, 1, function(x) ifelse(is.na(x) == T, 0, x)))
Swedishtb3 <- Swedishtb2[rowSums(Swedishtb2) > 0, ]
# apply collex function
absolutely  <- collex(data = Swedishtb3, cv1 = which(colnames(Swedishtb3) == Amplifiers[1]))
completely  <- collex(data = Swedishtb3, cv1 = which(colnames(Swedishtb3) == Amplifiers[2]))
extremely  <- collex(data = Swedishtb3, cv1 = which(colnames(Swedishtb3) == Amplifiers[3]))
highly  <- collex(data = Swedishtb3, cv1 = which(colnames(Swedishtb3) == Amplifiers[4]))
other  <- collex(data = Swedishtb3, cv1 = which(colnames(Swedishtb3) == Amplifiers[5]))
particularly  <- collex(data = Swedishtb3, cv1 = which(colnames(Swedishtb3) == Amplifiers[6]))
perfectly  <- collex(data = Swedishtb3, cv1 = which(colnames(Swedishtb3) == Amplifiers[7]))
really  <- collex(data = Swedishtb3, cv1 = which(colnames(Swedishtb3) == Amplifiers[8]))
so  <- collex(data = Swedishtb3, cv1 = which(colnames(Swedishtb3) == Amplifiers[9]))
totally  <- collex(data = Swedishtb3, cv1 = which(colnames(Swedishtb3) == Amplifiers[10]))
very  <- collex(data = Swedishtb3, cv1 = which(colnames(Swedishtb3) == Amplifiers[11]))
# extract informaltion
absolutely <- matrix(unlist(absolutely),ncol=10,byrow=TRUE)
completely <- matrix(unlist(completely),ncol=10,byrow=TRUE)
extremely <- matrix(unlist(extremely),ncol=10,byrow=TRUE)
highly <- matrix(unlist(highly),ncol=10,byrow=TRUE)
other <- matrix(unlist(other),ncol=10,byrow=TRUE)
particularly <- matrix(unlist(particularly),ncol=10,byrow=TRUE)
perfectly <- matrix(unlist(perfectly),ncol=10,byrow=TRUE)
really <- matrix(unlist(really),ncol=10,byrow=TRUE)
so <- matrix(unlist(so),ncol=10,byrow=TRUE)
totally <- matrix(unlist(totally),ncol=10,byrow=TRUE)
very <- matrix(unlist(very),ncol=10,byrow=TRUE)
# set up table with results
SwedishCollextb <- rbind(absolutely, completely, extremely, highly, other, 
                         particularly, perfectly, really, so, totally, very)
# apply collextbedit Function
SwedishCollex <- collextbedit(SwedishCollextb)
###########################################################################
# combine covar collex data frames
collexdfLanguages <- rbind(EnglishCollex[,c(1:11,15)], BulgarianCollex[,c(1:11,15)],
                           CzechCollex[,c(1:11,15)], DutchCollex[,c(1:11,15)],
                           FinnishCollex[,c(1:11,15)], FlemishCollex[,c(1:11,15)],
                           FrenchCollex[,c(1:11,15)], GermanCollex[,c(1:11,15)],
                           ItalianCollex[,c(1:11,15)], PolishCollex[,c(1:11,15)],
                           RussianCollex[,c(1:11,15)], SpanishCollex[,c(1:11,15)], 
                           SwedishCollex[,c(1:11,15)])
Language <- c(rep("English", nrow(EnglishCollex)), rep("Bulgarian", nrow(BulgarianCollex)), 
         rep("Czech", nrow(CzechCollex)), rep("Dutch", nrow(DutchCollex)),
         rep("Finnish", nrow(FinnishCollex)), rep("Flemish", nrow(FlemishCollex)),
         rep("French", nrow(FrenchCollex)), rep("German", nrow(GermanCollex)),
         rep("Italian", nrow(ItalianCollex)), rep("Polish", nrow(PolishCollex)),
         rep("Russian", nrow(RussianCollex)), rep("Spanish", nrow(SpanishCollex)),
         rep("Swedish", nrow(SwedishCollex)))
# create data frame
covarcoldf <- data.frame(Language, collexdfLanguages)
#convert Amplifiedo numeric
covarcoldf$Probability <- as.numeric(covarcoldf$Probability)
covarcoldf$CollStrength <- as.numeric(covarcoldf$CollStrength)
covarcoldf$OddsRatio <- as.numeric(covarcoldf$OddsRatio)
# inspect data
str(covarcoldf); head(covarcoldf); summary(covarcoldf$CollStrength)

###########################################################################
# rename data
p10d <- covarcoldf[, c(1:3, 9)] # Language, Amp, Adjective, CollStrength
#extract Adjectives present in all Languages
AmpAdjectiveLanguageftb <- ftable(p10d$Amp, p10d$Adjective, p10d$Language)
Amp <- unlist(attr(AmpAdjectiveLanguageftb, "row.vars")[1])
Adjective <- unlist(attr(AmpAdjectiveLanguageftb, "row.vars")[2])
Language <- unlist(attr(AmpAdjectiveLanguageftb, "col.vars")[1])
Adjectiver <- rep(Adjective, length(Amp))
ampr <- rep(Amp, each = length(Adjective))
# create new id variable
p10d$LanguageAdjective <- paste(p10d$Language, "_", p10d$Adjective, sep = "")
# reorder data frame
p10tb <- reshape(p10d, idvar = "LanguageAdjective", timevar = "Amp",direction = "wide")
colnames(p10tb)

# Language, Adjective, collstrength:quite, collstrength:so, collstrength:really, collstrength:very
p10tb <- p10tb[, c(2:4, 7, 10, 13, 16, 19, 22, 25, 28, 31, 34)] 
colnames(p10tb) <- c("Language", "Adjective",  "absolutely", "completely", "extremely",
                     "highly", "other", "particularly", "perfectly", "really", "so",
                     "totally", "very")
p10tb$Adjective <- as.factor(p10tb$Adjective)
p10tb$Language <- factor(p10tb$Language, levels=c("English", "Bulgarian", "Czech",
                                            "Dutch", "Finnish", "Flemish", "French", 
                                            "German", "Italian", "Polish", "Russian", 
                                            "Spanish", "Swedish"))
p10tb <- p10tb %>%
  select(Language, Adjective, very, so, really, extremely, completely, other) # select columns
frqadjp10 <- table(p10tb$Adjective)[which(table(p10tb$Adjective) == 13)]
frqadjp10 <- names(frqadjp10)[1:6]
p10tba <- p10tb %>%
  filter(Adjective %in% frqadjp10) # select rows
head(p10tba)
# start plot: all
p10a <- ggplot(p10tba, aes(x = Language, y = very)) +
  facet_wrap(vars(Adjective)) +
  geom_point(aes(y = very, color = "very"), size=.5) +
  geom_point(aes(y = so, color = "so"), size=.5) +
  geom_point(aes(y = really, color = "really"), size=.5) +
  geom_point(aes(y = extremely, color = "extremely"), size=.5) +
  geom_point(aes(y = completely, color = "completely"), size=.5) +
  geom_point(aes(y = other, color = "other"), size=.5) +
  guides(color=guide_legend(override.aes=list(fill=NA))) +
  scale_colour_manual(values=c("gray50", "goldenrod2", "gray70", "indianred4", "cornflowerblue", "black"),
                      name="Variant", 
                      breaks=c("very", "so", "really", "extremely", "completely", "other"), 
                      labels = c("very", "so", "really", "extremely", "completely", "other")) +
  theme_set(theme_light(base_size = 8)) +
  theme(legend.position="top", axis.text.x = element_text(size=8, angle=90)) +
  coord_cartesian(ylim = c(-2.5, 11)) +
  labs(x = "Language", y = "Collocation Strength (LOG(p), 10)") +
  guides(size = FALSE)+
  guides(alpha = FALSE)
ggsave(file = paste(imageDirectory,"CovarcollAmpIcleFreqAdjective1.png",sep="/"), width = 15, height = 10, units = c("cm"),  dpi = 320)
p10a

###
frqadjp10 <- table(p10tb$Adjective)[which(table(p10tb$Adjective) == 13)]
frqadjp10 <- names(frqadjp10)[7:12]
p10tbb <- p10tb %>%
  filter(Adjective %in% frqadjp10) # select rows
head(p10tbb)

# start plot: all
p10b <- ggplot(p10tbb, aes(x = Language, y = very)) +
  facet_wrap(vars(Adjective)) +
  geom_point(aes(y = very, color = "very"), size=1) +
  geom_point(aes(y = so, color = "so"), size=1) +
  geom_point(aes(y = really, color = "really"), size=1) +
  geom_point(aes(y = extremely, color = "extremely"), size=1) +
  geom_point(aes(y = completely, color = "completely"), size=1) +
  geom_point(aes(y = other, color = "other"), size=1) +
  guides(color=guide_legend(override.aes=list(fill=NA))) +
  scale_colour_manual(values=c("gray50", "goldenrod2", "gray70", "indianred4", "cornflowerblue", "black"),
                      name="Variant", 
                      breaks=c("very", "so", "really", "extremely", "completely", "other"), 
                      labels = c("very", "so", "really", "extremely", "completely", "other")) +
  theme_set(theme_light(base_size = 8)) +
  theme(legend.position="top", axis.text.x = element_text(size=8, angle=90)) +
  coord_cartesian(ylim = c(-2.5, 4)) +
  labs(x = "Language", y = "Collocation Strength (LOG(p), 10)") +
  guides(size = FALSE)+
  guides(alpha = FALSE)
ggsave(file = paste(imageDirectory,"CovarcollAmpIcleFreqAdjective2.png",sep="/"), width = 15, height = 10, units = c("cm"),  dpi = 320)
p10b














###########################################################################
#                  CHANGES IN Adjective FREQ
tb12 <- ampicle[ampicle$Variant != "0",] 
frqadj <- names(table(tb12$Adjective))[which(table(tb12$Adjective) > 100)]
tb12$Adjective <- ifelse(tb12$Adjective %in% frqadj, tb12$Adjective, "other")
tb12$Language <- factor(tb12$Language, 
                        levels=c("English", "Bulgarian", "Czech", "Dutch", "Finnish", 
                                 "Flemish", "French", "German", "Italian", "Polish", 
                                 "Russian", "Spanish", "Swedish"))
tb12 <- tb12%>%
  select(Language, Adjective, Variant) %>%
  filter(Variant != "0") %>%
  group_by(Language) %>%
  count(Adjective) %>%
  spread(Language, n, fill = 0, convert = FALSE) %>%
  select(Adjective, everything())
tb12 <- as.data.frame(tb12)
rownames(tb12) <- tb12$Adjective
tb12$Adjective <- NULL
tb12 <- apply(tb12, 2, function(x) {x <- round(x/sum(x)*100, 1)})
tb12
  
# save data to disc
write.table(tb12, "Table5.txt", sep = "\t", row.names = F)
###########################################################################
# start plot: Adjective
p11d <- as.data.frame(t(tb12))
p11 <- ggplot(p11d, aes(x = Language, y = other), size = 8) +
  geom_point(aes(y = other, color = "other"), size=1) +
  geom_point(aes(y = good, color = "good"), size=1) +
  geom_point(aes(y = different, color = "different"), size=1) +
  geom_point(aes(y = difficult, color = "difficult"), size=1) +
  geom_point(aes(y = hard, color = "hard"), size=1) +
  geom_point(aes(y = important, color = "important"), size=1) +
  geom_point(aes(y = little, color = "little"), size=1) +
  geom_point(aes(y = strong, color = "strong"), size=1) +
  guides(color=guide_legend(override.aes=list(fill=NA))) +
scale_colour_manual(values=c("grey30", "grey60", "goldenrod2",  "indianred4", "grey30", "goldenrod2", "grey60", "blue"),
                      name="", 
                      breaks=c("other", "good", "different", "difficult", "hard", "important", "little", "strong"), 
                      labels = c("other", "good", "different", "difficult", "hard", "important", "little", "strong")) +
  theme(legend.position="top") +
  theme_light(base_size = 8) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(x = "Language", y = "Percent of Adjectives") +
ggsave(file = paste(imageDirectory,"AdjectivefreqLanguage.png", sep="/"), width = 15,  height = 7.5, units = c("cm"),  dpi = 320)
p11


###########################################################################
#                  REGRESSION DATA SET
###########################################################################
# only amplified instances
reallyicle <- ampicle[ampicle$Amplified == 1,]
# inspect data
str(reallyicle); colnames(reallyicle)

# remove superfluous columns
reallyicle$ID <- NULL
reallyicle$File <- NULL
reallyicle$PreContext <- NULL
reallyicle$PostContext <- NULL
reallyicle$PreContextLong <- NULL
reallyicle$Amplified <- NULL
reallyicle$so <- NULL
reallyicle$extremely <- NULL
reallyicle$Variant <- NULL 
# define vector for data inspection
cltb <- c("Language", "Adjective", "Emotionality", "Function", "Priming", 
          "very", "really", "Gradabilty", "SemanticCategory")
# tabulate data
lapply(reallyicle[cltb], table)

# check data
reallyicle[which(reallyicle$SemanticCategory == "NoSemType"),]

# recategorzie Adjective
fradj <- names(table(reallyicle$Adjective))[which(table(reallyicle$Adjective) > 10)]
reallyicle$Adjective <- ifelse(reallyicle$Adjective %in% fradj, reallyicle$Adjective, "other")
# recategorzie Priming
reallyicle$Priming <- ifelse(reallyicle$Priming == "noprime", "NoPrime",
                             ifelse(reallyicle$Priming == "prime", "Prime", reallyicle$Priming))
# recategorzie SemanticCategory
reallyicle$SemanticCategory <- as.character(reallyicle$SemanticCategory)
reallyicle$SemanticCategory <- ifelse(reallyicle$SemanticCategory == "Age", "NoSemType",
                                      ifelse(reallyicle$SemanticCategory == "Color", "NoSemType", 
                                             ifelse(reallyicle$SemanticCategory == "Difficulty", "NoSemType",reallyicle$SemanticCategory)))
# factorize variables
clfct <- c("Language", "Adjective", "Emotionality",  "Function", "Priming",  
           "Gradabilty", "SemanticCategory")
reallyicle[clfct] <- lapply(reallyicle[clfct], factor)
# define vector for data inspection
cltb <- c("Language", "Adjective", "Emotionality",  "Function", "Priming",  
          "very", "really", "Gradabilty", "SemanticCategory")
# tabulate data
lapply(reallyicle[cltb], table)

# inspect data
nrow(reallyicle); str(reallyicle); colnames(reallyicle)

###############################################################
write.table(reallyicle, "ampicle_statz05.txt", row.names= F, sep = "\t")
###############################################################
#                   END PART 3
###############################################################
