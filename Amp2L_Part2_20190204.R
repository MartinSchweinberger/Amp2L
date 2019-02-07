##################################################################
# Titel:      How Corpus-Based Methods Can Support Language Learning 
#             and Teaching: An Analysis of Amplifier Use by L1-Speakers 
#             and Learners of English - Part 2
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
#             of Amplifier Use by L1-Speakers and Learners of English - Part 2",
#             unpublished R script, The University of Queensland.
###############################################################
rm(list=ls(all=T))                                      # clean current workspace
setwd("D:\\Uni\\Projekte\\02-Intensification\\Amp2L")   # set wd
library(stringr)                                        # load packages
options(stringsAsFactors = F)                           # set options
options(scipen = 999)                                   # set options
imageDirectory<-"images"                                # define image directory
corpus.locness <- "D:\\Uni\\Korpora\\Original\\LOCNESS" # specify path to corpus
###############################################################
# define files to load
corpus.files = list.files(path = corpus.locness, pattern = ".txt", all.files = T,
                           full.names = T, recursive = T, ignore.case = T, include.dirs = T)
Content <- sapply(corpus.files, function(x) {
  x <- scan(x, what = "char", sep = "", quote = "\"", quiet = T, skipNul = T)
  x <- gsub(" {2,}", " ", x)
  x <- str_trim(x, side = "both")
  x <- paste(x, sep = " ", collapse = " ")
})
Language <- rep("English", length(Content))
File <- gsub(".*/", "", corpus.files)
File <- gsub(".txt", "", File)
CleanContent <- gsub(" {2,}", " ", Content)
# combine vectors to data frame
locness <- as.data.frame(cbind(Language, File, Content, CleanContent))
# check data
head(locness[,1:2]); nrow(locness)

###############################################################
# save data to disc
write.table(locness, "locness_raw01.txt", sep = "\t", row.names = F, col.names = T, quote = T)
# WARNING: DO NOT ACTIVATE!
#locness <- read.table("locness_raw01.txt", sep = "\t", header = T)#, skipNul = T, quote = "\"", fill = T)
###############################################################
# split data into smaller chunks
pos01 <- locness$CleanContent[1]
pos02 <- locness$CleanContent[2]
pos03 <- locness$CleanContent[3]
pos04 <- locness$CleanContent[4]
pos05 <- locness$CleanContent[5]
pos06 <- locness$CleanContent[6]
pos07 <- locness$CleanContent[7]
pos08 <- locness$CleanContent[8]
pos09 <- locness$CleanContent[9]
pos10 <- locness$CleanContent[10]
pos11 <- locness$CleanContent[11]
pos12 <- locness$CleanContent[12]
pos13 <- substring(locness$CleanContent[13], 1, 200000)
pos14 <- substring(locness$CleanContent[13], 200001, 400000)
pos15 <- substring(locness$CleanContent[13], 400001, 600000)
pos16 <- substring(locness$CleanContent[13], 600001, 800000)
pos17 <- substring(locness$CleanContent[13], 800001, nchar(locness$CleanContent[13]))
pos18 <- locness$CleanContent[14]
# reload libraries
source("D:\\R/POStagObject.R") # for pos-tagging objects in R
library(NLP)
library(openNLP)
library(openNLPmodels.en)
# pos tagging data
locnesspos01 <- POStag(object = pos01)
locnesspos01 <- as.vector(unlist(locnesspos01))
writeLines(locnesspos01, con = "locnesspos01.txt", sep = "\n", useBytes = FALSE)
# chunk 2
locnesspos02 <- POStag(object = pos02)
locnesspos02 <- as.vector(unlist(locnesspos02))
writeLines(locnesspos02, con = "locnesspos02.txt", sep = "\n", useBytes = FALSE)
# chunk 03
locnesspos03 <- POStag(object = pos03)
locnesspos03 <- as.vector(unlist(locnesspos03))
writeLines(locnesspos03, con = "locnesspos03.txt", sep = "\n", useBytes = FALSE)
# chunk 04
locnesspos04 <- POStag(object = pos04)
locnesspos04 <- as.vector(unlist(locnesspos04))
writeLines(locnesspos04, con = "locnesspos04.txt", sep = "\n", useBytes = FALSE)
# chunk 05
locnesspos05 <- POStag(object = pos05)
locnesspos05 <- as.vector(unlist(locnesspos05))
writeLines(locnesspos05, con = "locnesspos05.txt", sep = "\n", useBytes = FALSE)
# chunk 06
locnesspos06 <- POStag(object = pos06)
locnesspos06 <- as.vector(unlist(locnesspos06))
writeLines(locnesspos06, con = "locnesspos06.txt", sep = "\n", useBytes = FALSE)
# chunk 07
locnesspos07 <- POStag(object = pos07)
locnesspos07 <- as.vector(unlist(locnesspos07))
writeLines(locnesspos07, con = "locnesspos07.txt", sep = "\n", useBytes = FALSE)
# chunk 08
locnesspos08 <- POStag(object = pos08)
locnesspos08 <- as.vector(unlist(locnesspos08))
writeLines(locnesspos08, con = "locnesspos08.txt", sep = "\n", useBytes = FALSE)
# chunk 09
locnesspos09 <- POStag(object = pos09)
locnesspos09 <- as.vector(unlist(locnesspos09))
writeLines(locnesspos09, con = "locnesspos09.txt", sep = "\n", useBytes = FALSE)
# chunk 010
locnesspos10 <- POStag(object = pos10)
locnesspos10 <- as.vector(unlist(locnesspos10))
writeLines(locnesspos10, con = "locnesspos10.txt", sep = "\n", useBytes = FALSE)
# chunk 011
locnesspos11 <- POStag(object = pos11)
locnesspos11 <- as.vector(unlist(locnesspos11))
writeLines(locnesspos11, con = "locnesspos11.txt", sep = "\n", useBytes = FALSE)
# chunk 12
locnesspos12 <- POStag(object = pos12)
locnesspos12 <- as.vector(unlist(locnesspos12))
writeLines(locnesspos12, con = "locnesspos12.txt", sep = "\n", useBytes = FALSE)
# chunk 13
locnesspos13 <- POStag(object = pos13)
locnesspos13 <- as.vector(unlist(locnesspos13))
writeLines(locnesspos13, con = "locnesspos13.txt", sep = "\n", useBytes = FALSE)
# chunk 14
locnesspos14 <- POStag(object = pos14)
locnesspos14 <- as.vector(unlist(locnesspos14))
writeLines(locnesspos14, con = "locnesspos14.txt", sep = "\n", useBytes = FALSE)
# chunk 15
locnesspos15 <- POStag(object = pos15)
locnesspos15 <- as.vector(unlist(locnesspos15))
writeLines(locnesspos15, con = "locnesspos15.txt", sep = "\n", useBytes = FALSE)
# chunk 16
locnesspos16 <- POStag(object = pos16)
locnesspos16 <- as.vector(unlist(locnesspos16))
writeLines(locnesspos16, con = "locnesspos16.txt", sep = "\n", useBytes = FALSE)
# chunk 17
locnesspos17 <- POStag(object = pos17)
locnesspos17 <- as.vector(unlist(locnesspos17))
writeLines(locnesspos17, con = "locnesspos17.txt", sep = "\n", useBytes = FALSE)
# chunk 18
locnesspos18 <- POStag(object = pos18)
locnesspos18 <- as.vector(unlist(locnesspos18))
writeLines(locnesspos18, con = "locnesspos18.txt", sep = "\n", useBytes = FALSE)
# list pos tagged elements
postag.files = c("locnesspos01.txt", "locnesspos02.txt", "locnesspos03.txt",  
                 "locnesspos04.txt","locnesspos05.txt", "locnesspos06.txt",  
                 "locnesspos07.txt", "locnesspos08.txt", "locnesspos09.txt", 
                 "locnesspos10.txt", "locnesspos11.txt", "locnesspos12.txt",
                 "locnesspos13.txt", "locnesspos14.txt", "locnesspos15.txt", 
                 "locnesspos16.txt", "locnesspos17.txt", "locnesspos18.txt")
# load pos tagged elements
locnesspos <- sapply(postag.files, function(x) {
  x <- scan(x, what = "char", sep = "\n", quote = "", quiet = T, skipNul = T)
  x <- gsub(" {2,}", " ", x)
  x <- str_trim(x, side = "both")
  x <- str_replace_all(x, fixed("\n"), " ")
})
locnesspos13 <- paste(locnesspos[13:17], sep = " ", collapse = " ")
locnesspos <- c(locnesspos[1:12], locnesspos13, locnesspos[18])                           
# unlist pos tagged elements
locness$locnesspos <- unlist(locnesspos)
###############################################################
# extract number of adjs per line
pstggd <- locness$locnesspos
lpstggd <- strsplit(pstggd, " ")
nlpstggd <- as.vector(unlist(sapply(lpstggd, function(x){
  x <- x[grep("[A-Z]{0,1}[a-z]{1,}\\/JJ[A-Z]{0,1}", x)]
  x <- length(x) } )))
rp <- nlpstggd
rp <- ifelse(rp == 0, 1, rp)
# load function for concordancing
source("D:\\R/ConcR_2.3_loadedfiles.R")
# set parameters for concordancing
pattern <- "[A-Z]{0,1}[a-z]{1,}\\/JJ[A-Z]{0,1}"
context <- 50
# extract all adjectives (concordance)
concjjlocness <- ConcR(locness$locnesspos, pattern, context, all.pre = FALSE)
# repeat rows in data frame as often as there are adjectives in it (if 0 adj, repeat once)
corpuslocnessadjdf <- locness[rep(seq(nrow(locness)), rp),]
# combine data sets
corpuslocnessadj <- data.frame(1:nrow(corpuslocnessadjdf), corpuslocnessadjdf, concjjlocness)
# remove rows without Tokens
amplocness <- corpuslocnessadj[is.na(corpuslocnessadj$Token) == F,]
# add clean column names
colnames(amplocness) <- c("ID", "Language", "File", "Content", "ContentClean",
                          "PosTaggedContent", "OriginalString", "PreContext", "Adjective", "PostContext")
# remove superfluous columns
amplocness$Content <- NULL
amplocness$ContentClean <- NULL
amplocness$PosTaggedContent <- NULL
amplocness$OriginalString <- NULL
# clean adjectives
amplocness$Adjective <- gsub(".*/JJR", "remove", amplocness$Adjective)
amplocness$Adjective <- gsub(".*/JJS", "remove", amplocness$Adjective)
amlocness <- amplocness[!amplocness$Adjective == "remove",]
amplocness$Adjective <- str_replace_all(amplocness$Adjective, fixed("/JJ"), "")
# add Vraiant column
amplocness$Variant <- gsub(".* ", "", str_trim(amplocness$PreContext, side = "both")) 
# inspect data
nrow(amplocness); head(amplocness)

###############################################################
# define amplifiers
amplifiers <- c("absolutely", "actually", "aggressively", 
                "amazingly", "appallingly", "awful", "awfully", 
                "badly", "bloody", "certainly", "clearly",
                "complete", "dead", "completely", "considerably", 
                "crazy", "decidedly", "definitely",  "distinctly", 
                "dreadfully", "enormously", "entirely", "especially", 
                "exactly", "exceedingly", "exceptionally", 
                "excruciatingly", "extraordinarily", "extremely",
                "fiercely", "firmly", "frightfully", "fucking", 
                "fully", "genuinely", "greatly",
                "grossly", "heavily", "highly", "hopelessly", 
                "horrendously", "hugely",
                "immediately", "immensely", "incredibly", 
                "infinitely", "intensely", "irrevocably",
                "mad", "mega", "mighty", #"most", "much", 
                "obviously", "openly", "overwhelmingly", "particularly", 
                "perfectly", "plenty", "positively", "precisely", 
                "pretty", "profoundly", "purely", 
                #"quite", 
                "real", "really", "remarkably", "seriously", 
                "shocking",   "significant", "significantly", "so", 
                "specially", "specifically", "strikingly",
                "strongly", "substantially", "super", "surely", 
                "terribly", "terrifically", 
                #"too",
                "total", "totally", "traditionally", #"true", 
                "truly", "ultra", "utterly", "very",
                "viciously", 
                #"well", 
                "wholly", "wicked", "wildly")
# clean ice locness data
amplocness$Function <- str_trim(amplocness$PostContext, side = "both")
amplocness$Function <- tolower(amplocness$Function)
amplocness$Function <- gsub(" {2,}", " ", amplocness$Function)
amplocness$Function <- sub(" ", "", amplocness$Function)
amplocness$Function <- gsub(" .*", "", amplocness$Function)
amplocness$Function <- gsub(".*/nn.*", "Attributive", amplocness$Function)
amplocness$Function <- ifelse(amplocness$Function == "Attributive", amplocness$Function, "Predicative")
# shorten post Context
amplocness$PostContext <- substr(amplocness$PostContext, 1, ifelse((nchar(amplocness$PostContext)+25) <25, maamplocness(nchar(amplocness$PostContext)), 25))
# pre Context
amplocness$PreContext <- str_trim(amplocness$PreContext, side = "both")
amplocness$PreContextLong <- amplocness$PreContext
amplocness$PreContextLong <- substr(amplocness$PreContextLong, ifelse(nchar(amplocness$PreContextLong)-25 <=0, 1, 
                                                                      nchar(amplocness$PreContextLong)-25), nchar(amplocness$PreContextLong))
amplocness$PreContext <- gsub(".* ", "", amplocness$PreContext)
# amplifier variant
amplocness$PreContext <- gsub("\\/.*", "", amplocness$PreContext)
amplocness$Variant <- ifelse(amplocness$PreContext %in% amplifiers, amplocness$PreContext, "0")
# amplified y/n
amplocness$Amplified <- ifelse(amplocness$Variant == "0", 0, 1) 
# adjective
amplocness$Adjective <- tolower(amplocness$Adjective)
# inspect data
nrow(amplocness); head(amplocness); table(amplocness$Variant)

# define forms that require removal
sups <- c(".*most.*", ".*more.*") 
negs <- c(".*not.*", ".*never.*", ".*n't.*")
downtoners <- c(".*sort/.*", ".*kind/.*", ".* bit/.*", ".*somewhat.*", ".*fairly.*", 
                ".*rather.*", ".*reasonably.*", ".*slightly.*", ".*comparatively.*", ".*semi.*", 
                ".*relatively.*", ".*little.*", ".*somehow.*", ".*almost.*", ".*partly.*", 
                ".*hardly.*", ".* less.*", ".*barely.*", ".* just/.*")
specialforms <- c(".* too.*", ".*quite.*")
PostContextdowntoners <- c(".*enough.*")
nonpropadj <- c("only", "much", "many", "cheaper", "cheaperr", "bests", "larger", "such", "other")
# check length of dataset
str(amplocness); head(amplocness); nrow(amplocness)#; table(amplocness$pint); head(amplocness$PreContextLong); head(amplocness$PreContextLong)

# code priming
prim1 <- c(rep(0, 1), amplocness$Variant[1:length(amplocness$Variant)-1])
prim2 <- c(rep(0, 2), amplocness$Variant[1:(length(amplocness$Variant)-2)])
prim3 <- c(rep(0, 3), amplocness$Variant[1:(length(amplocness$Variant)-3)])
primtb <- cbind(amplocness$Variant, prim1, prim2, prim3)

amplocness$Priming <- as.vector(unlist(apply(primtb, 1, function(x){
  x <- ifelse(x[1]== "0" , "noprime",
              ifelse(x[1] == x[2] | x[1] == x[3] | x[1] == x[4], "prime", "noprime"))
})))
# remove items that were not intensified by a minimum of 2 intensifier variants
nrow(amplocness)

# find items to be removed
supsidx <- unique(grep(paste(sups,collapse="|"), amplocness$PreContextLong, value=F))
negsidx <- unique(grep(paste(negs,collapse="|"), amplocness$PreContextLong, value=F))
downtonersidx <- unique(grep(paste(downtoners,collapse="|"), amplocness$PreContextLong, value=F))
specialformsidx <- unique(grep(paste(specialforms,collapse="|"), amplocness$PreContextLong, value=F))
PostContextdowntonersidx <- unique(grep(paste(PostContextdowntoners,collapse="|"), amplocness$PostContext, value=F))
nonpropadjidx <- unique(grep(paste(nonpropadj,collapse="|"), amplocness$Adjective, value=F))
# combine indices
idxs <- unique(c(supsidx, negsidx, downtonersidx, specialformsidx, PostContextdowntonersidx, nonpropadjidx))
# remove forms that require removal
amplocness <- amplocness[-idxs,]
# remove empty values
amplocness <- amplocness[!amplocness$Variant == "", ]

###############################################################
# save raw data to disc
write.table(amplocness, "amplocness_woneg02.txt", sep = "\t", row.names = F)
###############################################################

pintadjtb <- table(amplocness$Adjective, amplocness$Variant)
#pintadjtb <- pintadjtb[2:nrow(pintadjtb),]
pintadjtb <- pintadjtb[,2:ncol(pintadjtb)]
pintadjtb2 <- apply(pintadjtb, 1, function(x){
  x <- ifelse(x > 1, 1, x)})
pintadjtb3 <- colSums(pintadjtb2)
pintadjschildes <- names(pintadjtb3)[which(pintadjtb3 >=2 )]
amplocness <- amplocness[amplocness$Adjective %in% pintadjschildes, ]
nrow(amplocness)

# inspect adjectives
names(table(amplocness$Adjective))

# create vector with false adjectives
rmvadj <- c("uhr")
amplocness$remove <- ifelse(amplocness$Adjective %in% rmvadj, "remove", amplocness$Adjective)
amplocness <- amplocness[amplocness$remove != "remove",]
# remove superfluous columns
amplocness$remove <- NULL
colnames(amplocness)

# inspecta data
nrow(amplocness); length(table(amplocness$Adjective)); head(amplocness)

###############################################################
# save raw data to disc
write.table(amplocness, "amplocness_semlocnessan03.txt", sep = "\t", row.names = F)
###############################################################
# tabluate amplifiers
table(amplocness$Variant)[order(table(amplocness$Variant), decreasing = T)]

# create columns with intensifier freqs
amplocness$very <- ifelse(amplocness$Variant == "very", 1, 0) 
amplocness$so <- ifelse(amplocness$Variant == "so", 1, 0) 
amplocness$really <- ifelse(amplocness$Variant == "really", 1, 0) 
amplocness$extremely <- ifelse(amplocness$Variant == "extremely", 1, 0)
amplocness$other <- ifelse(amplocness$Variant == "extremely" |
                             amplocness$Variant == "really" | 
                             amplocness$Variant == "so" |
                             amplocness$Variant == "very", 1, 0)
# inspect results
head(amplocness)

###############################################################
# code freq of adj type by age group
frqadjtb <- table(amplocness$Language, amplocness$Adjective)
relfreqadjtb <- round(prop.table(frqadjtb, margin = 1)*100, 5)
relfreqadjdf <- as.data.frame(relfreqadjtb)
colnames(relfreqadjdf)[1:2] <- c("Language", "Adjective")
# add freq by date to data
amplocness <- merge(amplocness, relfreqadjdf, by=c("Language", "Adjective"))
# reorder data
amplocness <- amplocness[order(amplocness$ID),]
# inspect data
head(amplocness)

###############################################################
# code gradability
# gradability - manual classification
# done by trained grad students
# add gradability
ngrd_manual <- c("abject", "able", "abrasive", "abstract", "absurd", "abundant", "abusive",
                 "accurate", "acrimonious", "active", "advanced", "adverse", "affectionate", "afraid",
                 "aged", "aggressive", "agile", "agitated", "aimless", "airy", "alert", "alleged",
                 "allusive", "amazing", "ambitious", "amused", "amusing", "ancient", "angry",
                 "annoying", "anxious", "appalling", "apparent", "appealing", "applicable", "applied",
                 "appreciative", "apprehensive", "approachable", "appropriate", "approving",
                 "arduous", "arrogant", "ashamed", "associated", "astute", "athletic", "atrocious",
                 "attitudinal", "attractive", "authentic", "authoritarian   ", "authoritative",
                 "available", "aware", "awesome", "awful", "awkward", "awry", "bad", "bare", "base",
                 "battered", "beautiful", "beloved", "benevolent", "benign", "besetting", "bitter",
                 "bizarre", "bleak", "bleary", "bloody", "blotchy", "bold", "boppy", "bored",
                 "boring", "bossy", "brave", "brief", "bright", "brilliant", "broad", "browsing",
                 "brutal", "bubbly", "burly", "buzzy", "callous", "calm", "campy", "candid",
                 "capable", "careful", "careless", "casual", "cautious", "ceremonial", "challenging",
                 "changed", "charismatic", "charming", "cheap", "circumspect", "civic", "civil",
                 "civilised", "classy", "clever", "cocky", "cold", "collective", "colossal", "colourful",
                 "comfortable", "commandeered", "committed", "compatible", "compelling", "competent",
                 "competitive", "complex", "complicated", "conceivable", "concentrated", "concerned",
                 "confident", "confidential", "confused", "confusing", "considerable", "constructive",
                 "consultative", "contrived", "controversial", "convenient", "conventional", "converted",
                 "convinced", "cool", "corrupt", "cosy", "coy", "cramped", "crass", "crazy", "creative",
                 "criminal", "crippling", "critical", "cross", "crowded", "crucial", "cruel", "cumbersome",
                 "curious", "cushy", "cute", "cynical", "damaged", "damaging", "damp", "dangerous",
                 "daring", "darkened", "darn", "daunting", "dear", "debatable", "decent", "dedicated",
                 "deep", "defective", "defensive", "delicate", "delicious", "delighted", "delightful",
                 "dense", "dependent", "depressed", "desirable", "despairing", "desperate", "despicable",
                 "despondent", "destructive", "detailed", "detrimental", "devilish", "difficult",
                 "dirty", "disabled", "disadvantaged", "disappointed", "disappointing", "disastrous",
                 "disenchanted", "disgraceful", "disgusting", "dishonest", "disparaging", "distant",
                 "distinguished   ", "distorted", "distressed", "disturbed", "disturbing", "dizzy",
                 "dodgy", "dominant", "dotty", "double", "doubtful", "downhill", "dramatic", "dreadful",
                 "driving", "drunk", "drunken", "ductile", "dull", "dumb", "dusty", "dylan", "dynamic",
                 "dynamical", "eager", "early", "earnest", "earthy", "easterly", "eastern", "easy",
                 "eccentric", "economic", "edible", "effective", "efficient", "elderly", "elegant",
                 "eligible", "elitist", "elusive", "embarrassed", "embarrassing", "emergent", "eminent",
                 "emotional", "emotive", "encouraging", "energetic", "enlightening", "enormous",
                 "entertaining", "enthusiastic", "epic", "erudite", "estimated", "estranged", "everyday",
                 "evil", "exact", "exceptional", "excessive", "excited", "exciting", "expensive",
                 "experienced", "expert", "explicit", "express", "expressive", "extended", "extensive",
                 "extraordinary", "extravagant", "extroverted", "fabulous", "facile", "factual", "faint",
                 "familiar", "famous", "fanatic", "fancy", "fantastic", "fascinating", "fast", "fastidious",
                 "fat", "favourable", "favoured", "fearful", "feisty", "fergal", "ferocious", "fertile",
                 "fierce", "fiery", "filthy", "fine", "finished", "finite", "firm", "fitting", "fizzy",
                 "flexible", "fluffy", "fluttering", "foggy", "foolish", "forceful", "formalised",
                 "formidable", "fortunate", "frank", "frantic", "fraudulent", "fraught", "frenzied",
                 "frequent", "friendly", "frightening", "frightful", "frustrated", "frustrating",
                 "fulsome", "fun", "funny", "furious", "generous", "gentle", "giant", "gifted",
                 "gigantic", "glad", "glib", "glorious", "glossy", "good", "goodhearted", "gorgeous",
                 "gracious", "gradual", "grand", "grandiose", "grateful", "grave", "greasy", "great",
                 "grim", "groggy", "groovy", "gross", "grubby", "guilty", "gutless", "habitual",
                 "handsome", "handy", "hapless", "happy", "hard", "hardy", "harmful", "harmless",
                 "harmonic", "harsh", "hazardous", "hazy", "heavy", "hectic", "helpful", "hideous",
                 "high", "hilarious", "holy", "honest", "honorable", "honorary", "honourable",
                 "hooked", "hopeful", "hopeless", "horrendous", "horrible", "horrific", "hostile",
                 "hot", "huge", "humble", "humorous", "hungry", "hurt", "hysterical", "idealistic",
                 "igneous", "ignorant", "imaginative", "immature", "immediate", "immense", "imperative",
                 "important", "impotent", "impressive", "inane", "incompetent", "inconsistent",
                 "incorporate", "incorporated", "increased", "incredible", "incredulous", "indecent",
                 "independent", "individual", "individualistic ", "ineffective", "ineffectual",
                 "inept", "inevitable", "inexorable", "inexpensive", "inexperienced", "infamous",
                 "infertile", "informal", "infuriating", "injured", "innovative", "insatiable",
                 "insecure", "insidious", "inspirational   ", "inspired", "instructive", "insuperable",
                 "integrated", "intellectual", "intelligent", "intense", "intensive", "intimate",
                 "intolerant", "invaluable", "inventive", "ironic", "irresponsible", "irritable",
                 "irritating", "itchy", "jealous", "joyful", "justified", "justifying", "keen",
                 "labour", "ladylike", "lame", "large", "late", "layered", "lazy", "lean", "legitimate",
                 "leisurely", "less", "liberal", "liberating", "light", "likely", "limp", "little",
                 "loath", "locating", "lone", "lonely", "long", "loony", "loud", "lousy", "lovely",
                 "low", "loyal", "lucky", "lumbering", "luminous", "lumpy", "lunatic", "lush",
                 "mad", "magic", "magnificent", "major", "mandatory", "manipulated", "marginal",
                 "marvellous", "massive", "matrimonial", "mean", "meaningful", "measurable", "medical",
                 "medicinal", "mediocre", "mere", "mighty", "mild", "minatory", "minded", "minor",
                 "minted", "miraculous", "miscellaneous", "misleading", "mixed", "mock", "modal",
                 "modern", "modest", "modesty", "momentous", "monetary", "monstrous", "moral", "motivating",
                 "muddy", "muggy", "multiple", "mutual", "mystical", "mythical", "naive", "narrow", "nasty",
                 "naughty", "near", "nearby", "neat", "necessary", "neglected", "negligent", "nervous",
                 "net", "new", "nice", "noble", "noisy", "normal", "northern", "nostalgic", "notable",
                 "noted", "noteworthy", "noxious", "numerous", "objective", "obnoxious", "obscure",
                 "observant", "odd", "off", "oily", "okay", "old", "oldfashioned", "operatic", "optimistic",
                 "orderly", "ordinary", "orientated", "oriented", "other", "outdated", "outrageous",
                 "outstanding", "over", "overhanging", "overwhelming", "painful", "parky", "parlous",
                 "passionate", "pathetic", "patronising", "patterned", "peaked", "peculiar", "perforated",
                 "perishable", "pernicious", "perplexed", "perplexing", "persistent", "personal",
                 "persuasive", "perverted", "pessimistic", "petite", "petty", "phenomenal", "picturesque",
                 "pinkish", "plain", "pleasant", "pleased", "pleasing", "pleasurable", "plenty",
                 "poetic", "polite", "poor", "popular", "possessive", "potent", "potential", "powerful",
                 "practical", "pragmatic", "preachy", "precarious", "precious", "precise", "predatory",
                 "predictable", "prepared", "prescriptive", "pressing", "prestigious", "presumptuous",
                 "pretentious", "pretty", "prevalent", "primitive", "privileged", "prodigious", "productive",
                 "professional", "profitable", "profligate", "progressive", "prominent", "promotional",
                 "prone", "proper", "proportionate", "prospective", "prosperous", "protective", "proud",
                 "provocative", "prudential", "psycho", "psychotic", "public", "puerile", "purposeful",
                 "quaint", "qualitative", "queer", "quick", "quiet", "racist", "radical", "rainy",
                 "rampant", "rank", "rapid", "rapt", "rare", "rational", "rattled", "raw", "reactionary",
                 "reactive", "ready", "realistic", "reasonable", "recognisable", "recognised",
                 "recreational", "reddish", "reduced", "refreshing", "regretful", "regular", "relaxed",
                 "relaxing", "relentless", "relevant", "reliable", "reluctant", "remote", "required",
                 "resourceful", "respected", "responsible", "restless", "revealing", "rich", "ridiculous",
                 "risky", "robust", "rocky", "romantic", "rotten", "rough", "rowdy", "rude", "rumbling",
                 "rusty", "sacred", "sad", "safe", "sandy", "sarcastic", "satisfied", "satisfying",
                 "savage", "scarce", "scared", "sceptical", "scientific", "scrappy", "scratchy",
                 "scruffy", "scurrilous", "secret", "secular", "secure", "sedate", "seduced", "seedy",
                 "seismic", "selfconfessed", "selfish", "selfreliant", "senile", "sensational",
                 "sensible", "sensitive", "sentimental", "serious", "severe", "sexist", "sexual",
                 "sexy", "shadowy", "shaky", "shaped", "sharp", "shiny", "shitty", "shocking",
                 "short", "sick", "sickly", "silly", "simple", "sizeable", "skilful", "skilled",
                 "sleepy", "slight", "slim", "slippery", "sloppy", "slow", "small", "smart", "snoopy",
                 "snotty", "sociable", "social", "soft", "soggy", "solid", "sophisticated", "sore",
                 "sorry", "sour", "south", "spare", "sparkling", "spectacular", "spectral", "spiritual",
                 "spiteful", "splendid", "sporting", "starkly", "startling", "staunch", "steady",
                 "steamy", "steep", "stellar", "sticky", "stiff", "stimulating", "stoical", "stormy",
                 "strange", "strategic", "stressful", "stretched", "strict", "striking", "strong",
                 "structured", "stubborn", "stunning", "stupid", "subject", "subtle", "successful",
                 "suffering", "suitable", "sunny", "super", "superficial", "superior", "supernatural",
                 "supportive", "suppressed", "sure", "surplus", "surprised", "surprising", "susceptible",
                 "suspicious", "sustainable", "sweaty", "sweet", "swift", "sympathetic", "tacky",
                 "tactic", "talented", "tall", "tantalising", "tasteful", "tedious", "teensy", "temperate",
                 "tempting", "tended", "tense", "tentative", "terrible", "terrific", "theatrical",
                 "theoretical", "thermal", "thick", "thickened", "thin", "thirsty", "thoughtful",
                 "threatening", "thriving", "tight", "tiny", "tired", "titanic", "tony", "top", "topical",
                 "torrential", "tortious", "tortured", "torturous", "tough", "touring", "tragic",
                 "transcendental", "transferable", "traumatic", "treacherous", "tremendous", "trendy",
                 "tricky", "trim", "triumphal", "trivial", "troubled", "twee", "twisted", "typical",
                 "ugly", "ulterior", "unable", "unattractive", "unaware", "unbeknown", "unbelievable",
                 "uncaring", "uncertain", "unclear", "unctuous", "undecided", "undeniable", "undifferentiated",
                 "undignified", "uneven", "unexpected", "unfair", "unfamiliar", "unfavourable",
                 "unfit", "unflattering", "unforced", "unfortunate", "ungrateful", "unhappy", "unholy",
                 "unified", "unknown", "unlikely", "unlucky", "unorthodox", "unpleasant", "unreal",
                 "unseemly", "unsmiling", "unsocial", "unsound", "unstable", "unusual", "upset",
                 "uptight", "urban", "urbanised", "urgent", "useful", "vague", "vain", "valiant",
                 "valuable", "variable", "varied", "vast", "venerated", "vengeful", "versatile",
                 "vested", "veteran", "viable", "vigorous", "vile", "violent", "virtual", "visionary",
                 "visual", "vital", "vivid", "vocal", "volatile", "vulnerable", "wakeful", "warm",
                 "wayward", "weak", "weakly", "wealthy", "weary", "wee", "weird", "wet", "wicked",
                 "wide", "widespread", "wild", "willing", "wise", "wishful", "witty", "wobbly",
                 "wonderful", "wondrous", "worried", "worthwhile", "worthy", "wounded", "young",
                 "yukky", "yummy")
grd_manual <- c("delusive", "abdominal", "aboriginal", "absent", "absolute", "academic",
                "accented", "acceptable", "accessible", "accomplished", "accountable", "acoustical",
                "acrylic", "actual", "additional", "adequate", "adjacent", "administrative",
                "adolescent", "advantageous", "aerial", "affected", "affirmative", "affordable",
                "african", "aggregate", "agricultural", "albanian", "alive", "allergic", "alternative",
                "ambiguous", "american", "analogous", "analytical", "ancestral", "anecdotal",
                "angled", "anglican", "announced", "annual", "anonymous", "antarctic", "apocryphal",
                "aqueous", "arbitrary", "archaeological", "archaic", "arctic", "armed", "armoured",
                "artificial", "artistic", "asian", "asthmatic", "atmospheric", "atomic", "locnesssie",
                "locnesstralian", "locnesstrian", "authorised", "automatic", "autonomous", "average", "awake",
                "back", "backward", "baked", "balanced", "bald", "bankrupt", "basic", "bearded",
                "beneficial", "best", "biblical", "bibliographic", "binding", "biodegradeable",
                "biographical", "biological", "black", "blank", "blatant", "blind", "blonde", "blue",
                "bodily", "booed", "botanical", "bottom", "british", "broke", "broken", "brown",
                "bucketful", "budgetary", "bureaucratic", "burnt", "businesslike", "busy", "californian",
                "canonical", "capitalistic", "captive", "cardiac", "catholic", "cellular", "central",
                "centralised", "centred", "certain", "characteristic  ", "chartered", "cheated",
                "chemical", "chilean", "chinese", "chivalrous", "christian", "chromatic", "chronological",
                "churchy", "classic", "classical", "clean", "clear", "close", "closed", "coarse",
                "coated", "coherent", "cohesive", "coincidental", "colloquial", "coloured", "coming",
                "commercial", "common", "compact", "comparable", "complete", "compound", "comprehensive",
                "compulsory", "computerised", "conceptual", "concrete", "confessional", "confirmed",
                "conscious", "conservative", "consistent", "constant", "constituent", "contemporary",
                "contestable", "continual", "contraceptive", "contrary", "cooked", "cooking",
                "corporate", "correct", "cracked", "crushed", "cubic", "cultural", "curly", "current",
                "customary", "cut", "daily", "dark", "dead", "deadly", "deaf", "decisive", "definite",
                "definitive", "deliberate", "democratic", "demographic", "determined", "diagnostic",
                "diagonal", "dietetic", "different", "digestive", "digital", "diplomatic", "direct",
                "discursive", "displaced", "disqualified", "distinct", "distinctive", "diverse", "divine",
                "domestic", "down", "downward", "dry", "dual", "dubious", "dummy", "dutch", "east",
                "educational", "effluent", "egalitarian", "electable", "electric", "electrical", "electronic",
                "elemental", "empty", "endemic", "endless", "english", "enough", "enrolled", "entailed",
                "entire", "equal", "equatorial", "equestrian", "equitable", "equivalent", "eritrean",
                "essential", "estonian", "ethic", "ethiopian", "ethnic", "european", "ewen", "exalted",
                "excellent", "executive", "exiguous", "existent", "existing", "exotic", "expected",
                "experimental", "explosive", "exponential", "external", "extinct", "extra", "extreme",
                "fair", "fake", "false", "far", "fatal", "favourite", "federal", "federated",
                "fellow", "female", "feminist", "feudal", "few", "fictional", "final", "financial",
                "first", "fixed", "flagged", "flannelled", "flat", "fleet", "flowing", "fluent",
                "fluid", "focused", "folded", "folding", "following", "foreign", "foremost", "formal",
                "forthcoming", "forward", "fossil", "foster", "founding", "fragile", "free", "french",
                "fresh", "frisian", "front", "frontal", "frosted", "fucking", "full", "fundamental",
                "funded", "further", "future", "gaelic", "gay", "general", "generational", "generic",
                "genuine", "geographical", "geological", "geotechnical", "german", "germanic", "glandular",
                "global", "gold", "golden", "governmental", "granulitic", "graphical", "gray", "greek",
                "green", "grey", "guaranteed", "half", "halved", "halving", "healthy", "hereditary",
                "heterogeneous   ", "heterogenious", "hidden", "historic", "historical", "holistic",
                "homosexual", "hooped", "horizontal", "hourly", "human", "humanitarian", "humiliating",
                "hungary", "hydroplaning", "hypocritical", "hypothetical", "iambic", "ideal", "identical",
                "ideological", "idle", "ill", "illegal", "imaginable", "immune", "imperial", "implicit",
                "implied", "impossible", "improved", "inaccessible", "inaccurate", "inadequate", "inclusive",
                "incoming", "incorrect", "incumbent", "indian", "indifferent", "indigenous", "indispensable",
                "indisputable", "industrial", "inefficient", "inescapable", "inexplicable", "infallible",
                "inflatable", "inflated", "informed", "infrequent", "inherent", "initial", "innate",
                "inner", "innocent", "innumerable", "inorganic", "inside", "insignificant", "instant",
                "instrumental", "insufficient", "intact", "integral", "intentional", "interactive",
                "intercultural", "interested", "interesting", "internal", "international", "interrupted",
                "intervening", "intriguing", "intrinsic", "inverted", "iraq", "irish", "irrelevant",
                "irrespective", "islamic", "italian", "japanese", "jewish", "joint", "journalistic",
                "judicial", "judicious", "junior", "just", "last", "latter", "leading", "learned",
                "learnt", "left", "lefthand", "legal", "legged", "legislative", "lesbian", "liable",
                "lime", "limited", "linear", "linguistic", "liquid", "literary", "liturgical", "live",
                "loaded", "local", "logarithmic", "logical", "logistic", "lost", "macrocyclic",
                "magnetic", "main", "male", "marine", "marked", "married", "masqueraded", "masterly",
                "materialistic   ", "maternal", "mathematical", "mature", "maximum", "mechanistic",
                "medieval", "mega", "melodic", "mental", "messy", "metamorphic", "meterological",
                "metrical", "metropolitan", "mexican", "micro", "microeconomic", "mid", "middle",
                "militaristic", "military", "milky", "minimal", "minimalist", "minimum", "ministerial",
                "missionary", "mobile", "moderate", "molecular", "molten", "monotonous", "mundane",
                "muscovite", "musical", "mutant", "naked", "narrative", "nasal", "natal", "national",
                "nationwide", "native", "natural", "nautical", "naval", "nazi", "needy", "negative",
                "neurotic", "next", "nitric", "north", "noticeable", "now", "nuclear", "obligatory",
                "obvious", "occasional", "octave", "official", "olympic", "ongoing", "only", "onward",
                "open", "operational", "opposed", "opposite", "optical", "optimum", "optional", "oral",
                "orange", "orchestral", "orchestrated", "organic", "original", "outside", "overlapping",
                "pacific", "painless", "pakistani", "parallel", "paramount", "parental", "parliamentary",
                "partial", "particular", "partisan", "passive", "past", "pastoral", "paternal",
                "paternalistic", "patriarchal", "patriotic", "perfect", "peripheral", "permanent",
                "permissive", "pertinent", "peruvian", "philosophical", "phonetic", "physical", "pink",
                "plastic", "pluralistic", "polar", "political", "politicised", "polynesian", "pornographic",
                "portable", "positive", "possible", "practicable", "preconceived", "preferential",
                "preferred", "pregnant", "preliminary", "presbyterian", "present", "presidential",
                "previous", "prewarned", "priceless", "primary", "prime", "principal", "prior", "pristine",
                "private", "privatised", "probable", "procedural", "programmed", "prolonged", "pronged",
                "proportional", "provincial", "psychiatric", "psychic", "pure", "purple", "quantifiable",
                "quantitative", "racial", "radioactive", "random", "readable", "real", "rear", "recent",
                "recycled", "red", "redundant", "reformed", "regional", "registered", "regulated",
                "regulatory", "reissued", "related", "relational", "relative", "remarkable", "remedial",
                "reportable", "reported", "residential", "respective", "resulting", "retrospective",
                "reusable", "reverse", "revolutionary", "ridged", "right", "rightful", "righthand",
                "rigid", "rigorous", "romanian", "rotary", "round", "royal", "ruined", "rural",
                "russian", "same", "samoan", "sane", "saturated", "scandinavian", "scholastic",
                "scottish", "scriptural", "seasonal", "secondary", "securing", "selected", "selective",
                "selfstyled", "senior", "senseless", "separate", "separated", "serial", "sheer",
                "siberian", "significant", "silver", "similar", "simultaneous", "sincere", "singaporean",
                "single", "skinned", "sleeveless", "sliced", "smokefree", "smooth", "sober", "socialist",
                "sociodemographic", "socioeconomic", "sole", "solitary", "soluble", "southern",
                "southwest", "sovereign", "soviet", "spanish", "special", "specialised", "specific",
                "spinal", "spontaneous", "spurious", "square", "stable", "stagnant", "standard",
                "stated", "stationary", "statistical", "statutory", "steely", "stereo", "stolen",
                "straight", "stratospheric", "striped", "structural", "subconscious", "subordinate",
                "subset", "substantial", "substantive", "suburban", "sudden", "sufficient", "suggestive",
                "sundry", "superheated", "supplementary", "supreme", "surgical", "sustained", "swedish",
                "swiss", "swollen", "symbolic", "synthetic", "technical", "technological", "temporary",
                "terminal", "territorial", "textual", "textural", "thematic", "thorough", "thoroughgoing",
                "timely", "total", "totalitarian", "toxic", "traditional", "transmitted", "traversable",
                "true", "twin", "ultimate", "unacceptable", "unaffected", "unallocated", "unannounced",
                "unanswered", "unbeaten", "unbiased", "unblemished", "unchanged", "uncoordinated",
                "under", "undisclosed", "undone", "unemployed", "unequal", "unexpired", "unfilled",
                "unfurnished", "unique", "universal", "unlimited", "unnatural", "unnecessary", "unoccupied",
                "unofficial", "unplayable", "unpopular", "unprecedented", "unprejudiced", "unpretentious",
                "unpromising", "unreceptive", "unregulated", "unrelated", "unresolved", "unrhymed",
                "unseeded", "unseen", "unselective", "unselfish", "unspecified", "unspoilt", "unstressed",
                "unsubsidised", "untold", "untrue", "unvarnished", "unwanted", "unwarranted", "unwilling",
                "unwrinkled", "upward", "usable", "useless", "usual", "utter", "vacant", "valid",
                "various", "veiled", "venetian", "verbal", "verbatim", "verifiable", "vertical",
                "volcanic", "voluntary", "weekly", "west", "western", "white", "whole", "wilful",
                "wooden", "woollen", "written", "wrong", "yellow", "youthful", "religious")
# gradability - data driven classification
# determine which adjectives are gradable
alladj <- names(table(amplocness$Adjective)) 
grd1 <- names(table(amplocness$Adjective[amplocness$Variant == "very"])) # adjs with very are gradable
grd2 <- names(table(amplocness$Adjective[amplocness$Variant == "extremely"])) # adjs with extremely are gradable
ngrd1 <- names(table(amplocness$Adjective[amplocness$Variant == "completely"])) # adjs with completely are gradable
ngrd2 <- names(table(amplocness$Adjective[amplocness$Variant == "total"])) # adjs with total are gradable
ngrd3 <- names(table(amplocness$Adjective[amplocness$Variant == "totally"])) # adjs with totally are gradable
ngrd4 <- names(table(amplocness$Adjective[amplocness$Variant == "utterly"])) # adjs with utterly are gradable
ngrd5 <- names(table(amplocness$Adjective[amplocness$Variant == "absolutely"])) # adjs with absolutely are gradable
# create vector with gradable adjectives
grdadj <- intersect(grd1, grd2)
# create vector with non-gradable adjectives
ngrdadj <- c(ngrd1, ngrd2, ngrd3, ngrd4, ngrd5)
ngrdadj <- names(table(ngrdadj))
# find elements that occur in both groups
bthgrdngrd1 <- grdadj[grdadj %in% ngrdadj]
bthgrdngrd2 <- ngrdadj[ngrdadj %in% grdadj]
bthgrdngrd <- names(table(c(bthgrdngrd1, bthgrdngrd2)))
# extract adjs that are clearly gradable 
grd_datadriven <- grdadj[!grdadj %in% bthgrdngrd]
# extract adjs that are clearly not gradable 
ngrd_datadriven <- ngrdadj[!ngrdadj %in% bthgrdngrd]
# find elements that are neither in gradable nor in nongradable
gradnongrad <- names(table(c(grd_datadriven, ngrd_datadriven)))
nagrad <- alladj[!alladj %in% gradnongrad]
naadjs <- names(table(c(nagrad, bthgrdngrd)))
###############################################################
# combine data driven and manual classification
# adj unclassified by datadriven now assigned value based on manual coding
grd_add <- naadjs[naadjs %in% grd_manual]
ngrd_add <- naadjs[naadjs %in% ngrd_manual]
# combine data driven and manual coding
grdlocness <- names(table(c(grd_datadriven, grd_add)))
ngrdlocness <- names(table(c(ngrd_datadriven, ngrd_add)))
# check which adjs are still unclassified
nagradlocness1 <- alladj[!alladj %in% grdlocness]
nagradlocness <- nagradlocness1[!nagradlocness1 %in% ngrdlocness]
# inspect unclassified adj
nagradlocness

# inspect length of vectors
length(alladj); length(grdlocness); length(ngrdlocness); length(nagradlocness)

# add gradability coding
amplocness$Gradabilty <- ifelse(amplocness$Adjective %in% grdlocness, "Gradable", amplocness$Adjective)
amplocness$Gradabilty <- ifelse(amplocness$Gradabilty %in% ngrdlocness, "NotGradable", amplocness$Gradabilty)
amplocness$Gradabilty <- ifelse(amplocness$Gradabilty  == "Gradable" | 
                                  amplocness$Gradabilty  == "NotGradable", 
                                amplocness$Gradabilty, "GradabilityUndetermined")
# inspect data
head(amplocness); table(amplocness$Gradabilty)

# save table of gradable and non-gradable adj
gradtb <- table(amplocness$Adjective, amplocness$Gradabilty)
gradtb <- data.frame(gradtb)
colnames(gradtb) <- c("adj", "grad", "freq")
head(gradtb)

# dave table to file
write.table(gradtb, "gradtb.txt", sep = "\t", row.names = T)
###############################################################
# add semantic types (tagliamonte 2008, based on dixon 1977)
# dimension = semdim (e.g. big, large, little, small, long, short, wide, narrow, thick)
# difficulty = semdif (e.g. difficult, simple)
# physical property = (e.g. hard, soft, heavy, light, rough, smooth, hot, sweet)
# color = semcol (e.g. black, white, red)
# human propensity: semhup (e.g. jealous, happy, kind, clever, generous, gay, rude)
# age = semage (e.g. new, young, old) 
# value (e.g. good, bad, proper, perfect, excellent, delicious, poor), 
# speed Speed (fast, quick, slow)
# position (e.g. right, left, near, far)
# other

# age
semage <- c("actual", "adolescent", "aged", "ancestral", "ancient", "annual", "archaeological", "archaic", 
            "biographical", "contemporary", "elderly", "foster", "generational", "historic", "historical", 
            "immature", "junior", "late", "mature", "medieval", "modern", "old", "oldfashioned", "outdated", 
            "past", "preliminary", "present", "primary", "prime", "prior", "puerile", "recent", "seasonal", 
            "senile", "senior", "temporary", "topical", "veteran", "young", "youthful")
# color
semcol <- c("colourful", "darkened", "pinkish", "reddish", "black", "blue", "brown",
            "coloured", "dark", "gold", "golden", "gray", "green", "grey", "lime", "marine",
            "orange", "pink", "purple", "red", "silver", "white", "yellow")
semdif <- c("basic", "complicated", "difficult", "easy", "elusive", "facile", "precarious",
            "risky", "simple", "stressful", "tricky", "twisted", "unpromising")
# dimension
semdim <- c("adjacent", "angled", "arctic", "back", "backward", "big", "bottom", "brief", "bright", 
            "broad", "central", "centralised", "centred", "close", "compact", "deep", "diagonal", 
            "direct", "distant", "distorted", "down", "downward", "early", "east", "easterly", "eastern", 
            "endemic", "endless", "equatorial", "european", "ewen", "far", "few", "first", "flat", 
            "foreign", "foremost", "forthcoming", "forward", "free", "front", "frontal", "further", 
            "geographical", "giant", "gigantic", "global", "grand", "half", "halved", "halving", 
            "high", "horizontal", "huge", "inner", "inside", "internal", "international", "large", 
            "last", "latter", "left", "linear", "little", "local", "locating", "long", "low", "massive", 
            "micro", "mid", "middle", "minimal", "minimalist", "minimum", "minor", "misleading", 
            "narrow", "national", "nationwide", "native", "near", "nearby", "next", "north", "northern", 
            "off", "onward", "orientated", "outside", "over", "overhanging", "overlapping", "pacific", 
            "parallel", "paramount", "peripheral", "petite", "polar", "proportional", "provincial", 
            "public", "rear", "regional", "remote", "reverse", "round", "rural", "separate", 
            "separated", "short", "sizeable", "slight", "small", "south", "southern", "southwest", 
            "spinal", "square", "steep", "stratospheric", "suburban", "super", "tall", "teensy", 
            "terminal", "territorial", "thick", "thickened", "thin", "tight", "tiny", "titanic", 
            "top", "torrential", "touring", "tremendous", "under", "universal", "unseeded", "upward", 
            "urban", "urbanised", "vast", "vertical", "warped", "wee", "west", "western", "wide", "widespread")
semhup <- c("able", "abrasive", "abusive", "academic", "accomplished", "advanced", "adverse", "afraid", 
            "aggressive", "aimless", "amused", "amusing", "analytical", "angry", "anxious", "appreciative", 
            "apprehensive", "ashamed", "astute", "aware", "benevolent", "besetting", "blind", "bold", "bossy", 
            "brave", "brutal", "busy", "callous", "candid", "capable", "careful", "challenging", 
            "charismatic", "cheated", "clever", "cocky", "compelling", "competent", "competitive", 
            "concerned", "confident", "consultative", "convinced", "creative", "cross", "cruel", "cute", 
            "cynical", "delighted", "depressed", "despairing", "desperate", "despondent", "disappointed", 
            "dodgy", "dotty", "dubious", "dull", "dumb", "eager", "elitist", "embarrassed", "emotional", "encouraging", 
            "entertaining", "enthusiastic", "erudite", "evil", "excited", "fanatic", "fearful", "ferocious", 
            "fierce", "foolish", "forceful", "fortunate", "fraudulent", "friendly", "frustrated", "fun", 
            "funny", "furious", "generous", "gifted", "glad", "goodhearted", "gracious", "grateful", 
            "grim", "gross", "gutless", "hairy", "hapless", "happy", "hopeful", "hopeless", "horrible", "hostile", 
            "hysterical", "ignorant", "ill", "imperative", "incompetent", "inexorable", "inexperienced", 
            "infallible", "informed", "insatiable", "insidious", "intellectual", "intelligent", "intriguing", 
            "inventive", "jealous", "joyful", "keen", "lazy", "learned", "learnt", "loath", "lone", "lonely", 
            "lucky", "lunatic", "mad", "mean", "minded", "motivating", "nasty", "nervous", "nice", 
            "optimistic", "passionate", "patronising", "pessimistic", "pleased", "polite", "poor", 
            "preachy", "prepared", "presumptuous", "primitive", "procedural", "professional", "promotional", 
            "proud", "prudential", "psycho", "puzzled", "rapt", "rational", "regretful", "relentless", 
            "resourceful", "respected", "rich", "romantic", "rowdy", "rude", "sad", "sane", "sarcastic", 
            "satisfied", "satisfying", "scared", "sceptical", "selective", "selfish", "sensitive", 
            "sentimental", "sick", "silly", "skilful", "skilled", "smart", "snotty", "sociable", "sophisticated", 
            "sorry", "sovereign", "spiteful", "staunch", "strategic", "strict", "stubborn", "stupid", 
            "suffering", "superior", "supportive", "suspicious", "tactic", "talented", "technical", "treacherous", 
            "troubled", "unable", "unanswered", "unaware", "uncaring", "ungrateful", "unhappy", "unsmiling", 
            "unsocial", "upset", "valiant", "valid", "vengeful", "vile", "wicked", "willing", "wise", "witty", 
            "worried")
# physical property
semphy <- c("cheap", "clear", "cold", "comfortable", "cool", "dark", "different", "dry", "flexible", "hard", 
            "heavy", "hot", "light", "neat", "obvious", "quick", "quiet", "real", "same", "scarce", "similar", 
            "slow", "strong", "sweet", "tidy", "warm")
# Value
semval <- c("amazing", "appropriate", "awful", "bad", "beautiful", "bizarre", "boring", "brilliant", 
            "competitive", "counterproductive", "dangerous", "easy", "effective", "efficient", "essential", "excellent", 
            "exciting", "expensive", "fantastic", "fat", "good", "great", "important", "inadequate", "interesting", 
            "new", "original", "painful", "pathetic", "popular", "positive", "relevant", "ridiculous", "right", 
            "scary", "serious", "simple", "special", "strange", "sure", "terrible", "tough", "trendy", "true", "ugly",  
            "uncomfortable", "unrealistic", "unusual", "useful", "useless", "weird", "wealthy", "worser", "worthwhile", "wrong", "yummy")
# add semantic type classification
amplocness$SemanticCategory <- ifelse(amplocness$Adjective %in% semage, "Age", amplocness$Adjective)
amplocness$SemanticCategory <- ifelse(amplocness$SemanticCategory %in% semcol, "Color", amplocness$SemanticCategory)
amplocness$SemanticCategory <- ifelse(amplocness$SemanticCategory %in% semdif, "Difficulty", amplocness$SemanticCategory)
amplocness$SemanticCategory <- ifelse(amplocness$SemanticCategory %in% semdim, "Dimension", amplocness$SemanticCategory)
amplocness$SemanticCategory <- ifelse(amplocness$SemanticCategory %in% semhup, "HumanPropensity", amplocness$SemanticCategory)
amplocness$SemanticCategory <- ifelse(amplocness$SemanticCategory %in% semphy, "PhysicalProperty", amplocness$SemanticCategory)
amplocness$SemanticCategory <- ifelse(amplocness$SemanticCategory %in% semval, "Value", amplocness$SemanticCategory)
amplocness$SemanticCategory <- ifelse(amplocness$SemanticCategory == "Age" | amplocness$SemanticCategory == "Color" | amplocness$SemanticCategory == "Difficulty" | 
                                        amplocness$SemanticCategory == "Dimension" | amplocness$SemanticCategory == "HumanPropensity" |
                                        amplocness$SemanticCategory == "PhysicalProperty" | amplocness$SemanticCategory == "Value",  amplocness$SemanticCategory, "NoSemType")
# table sem class of tokens
table(amplocness$SemanticCategory)

# check classification
names(table(amplocness$Adjective[which(amplocness$SemanticCategory == "NoSemType")]))

# inspect data
head(amplocness); nrow(amplocness); length(table(amplocness$Adjective))

###############################################################
# load library
library(syuzhet)
# code emotion
class_emo <- get_nrc_sentiment(amplocness$Adjective)
# process sentiment
amplocness$Emotionality <- as.vector(unlist(apply(class_emo, 1, function(x){
  x <- ifelse(x[9] == 1, "NegativeEmotional",
              ifelse(x[10] == 1, "PositiveEmotional", "NonEmotional")) } )))
# revert order of factor Emotionality
amplocness$Emotionality <- factor(amplocness$Emotionality, levels = c("NonEmotional", "NegativeEmotional", "PositiveEmotional"))
# inspect data
head(amplocness); str(amplocness); nrow(amplocness)

###############################################################
# save raw data to disc
write.table(amplocness, "amplocness_clean04.txt", sep = "\t", row.names = F)
###############################################################
#                        END PART 2
###############################################################
