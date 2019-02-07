##################################################################
# Titel:      How Corpus-Based Methods Can Support Language Learning 
#             and Teaching: An Analysis of Amplifier Use by L1-Speakers 
#             and Learners of English - Part 1
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
#             of Amplifier Use by L1-Speakers and Learners of English - Part 1",
#             unpublished R script, The University of Queensland.
###############################################################
rm(list=ls(all=T))                                      # clean current workspace
setwd("D:\\Uni\\Projekte\\02-Intensification\\Amp2L")   # set wd
library(stringr)                                        # load packages
options(stringsAsFactors = F)                           # set options
options(scipen = 999)                                   # set options
imageDirectory<-"images"                                # define image directory
corpus.icle <- "D:\\Uni\\Korpora\\Original\\ICLE\\DATA" # specify path to corpus
###############################################################
# define files to load
corpus.files = list.files(path = corpus.icle, pattern = NULL, all.files = T,
                          full.names = T, recursive = T, ignore.case = T, include.dirs = T)
icle <- lapply(corpus.files, function(x) {
  x <- scan(x, what = "char", sep = "", quote = "", quiet = T, skipNul = T)
  x <- gsub(" {2,}", " ", x)
  x <- str_trim(x, side = "both")
  Speaker<- x[1] 
  Content <- x[2:length(x)]
  Content <- paste(Content, sep = " ", collapse = " ")
  Content <- unlist(Content)
  icle <- cbind(Speaker, Content)
  return(icle)
})
ID <- sapply(icle,function(x) x[1])
LanguageSpeaker <- gsub("-[0-9]{2,}.*", "", ID)
LanguageSpeaker <- gsub("<ICLE-", "", LanguageSpeaker)
Language <- gsub("-.*", "", LanguageSpeaker)
File <- gsub(".*-", "", ID)
File <- gsub("\\..*", "", File)
Speaker <- gsub(".*-", "", LanguageSpeaker)
Content <- sapply(icle,function(x) x[2])
CleanContent <- gsub("\"", " ", Content, fixed = T)
CleanContent <- gsub(" {2,}", " ", CleanContent)
# combine vectors to data frame
icle <- as.data.frame(cbind(ID, Language, File, Speaker, Content, CleanContent))
# check data
head(icle[,1:4]); nrow(icle); icle[1320:1330,1:4]
###############################################################
# save data to disc
write.table(icle, "icle_raw01.txt", sep = "\t", row.names = F, col.names = T, quote = T)
# WARNING: DO NOT ACTIVATE!
#icle <- read.table("icle_raw01.txt", sep = "\t", header = T)#, skipNul = T, quote = "\"", fill = T)
###############################################################
# split data into smaller chunks
pos01 <- icle$CleanContent[1:500]
pos02 <- icle$CleanContent[501:1000]
pos03 <- icle$CleanContent[1001:1500]
pos04 <- icle$CleanContent[1501:2000]
pos05 <- icle$CleanContent[2001:2500]
pos06 <- icle$CleanContent[2501:3000]
pos07 <- icle$CleanContent[3001:3500]
pos08 <- icle$CleanContent[3501:nrow(icle)]
# reload libraries
source("D:\\R/POStagObject.R") # for pos-tagging objects in R
library(NLP)
library(openNLP)
library(openNLPmodels.en)
# pos tagging data
iclepos01 <- POStag(object = pos01)
iclepos01 <- as.vector(unlist(iclepos01))
writeLines(iclepos01, con = "iclepos01.txt", sep = "\n", useBytes = FALSE)
# chunk 2
iclepos02 <- POStag(object = pos02)
iclepos02 <- as.vector(unlist(iclepos02))
writeLines(iclepos02, con = "iclepos02.txt", sep = "\n", useBytes = FALSE)
# chunk 03
iclepos03 <- POStag(object = pos03)
iclepos03 <- as.vector(unlist(iclepos03))
writeLines(iclepos03, con = "iclepos03.txt", sep = "\n", useBytes = FALSE)
# chunk 04
iclepos04 <- POStag(object = pos04)
iclepos04 <- as.vector(unlist(iclepos04))
writeLines(iclepos04, con = "iclepos04.txt", sep = "\n", useBytes = FALSE)
# chunk 05
iclepos05 <- POStag(object = pos05)
iclepos05 <- as.vector(unlist(iclepos05))
writeLines(iclepos05, con = "iclepos05.txt", sep = "\n", useBytes = FALSE)
# chunk 06
iclepos06 <- POStag(object = pos06)
iclepos06 <- as.vector(unlist(iclepos06))
writeLines(iclepos06, con = "iclepos06.txt", sep = "\n", useBytes = FALSE)
# chunk 07
iclepos07 <- POStag(object = pos07)
iclepos07 <- as.vector(unlist(iclepos07))
writeLines(iclepos07, con = "iclepos07.txt", sep = "\n", useBytes = FALSE)
# chunk 08
iclepos08 <- POStag(object = pos08)
iclepos08 <- as.vector(unlist(iclepos08))
writeLines(iclepos08, con = "iclepos08.txt", sep = "\n", useBytes = FALSE)
# list pos tagged elements
postag.files = c("iclepos01.txt", "iclepos02.txt", "iclepos03.txt",  "iclepos04.txt", 
                 "iclepos05.txt", "iclepos06.txt",  "iclepos07.txt", "iclepos08.txt")
# load pos tagged elements
iclepos <- sapply(postag.files, function(x) {
  x <- scan(x, what = "char", sep = "\n", quote = "", quiet = T, skipNul = T)
  x <- gsub(" {2,}", " ", x)
  x <- str_trim(x, side = "both")
  x <- str_replace_all(x, fixed("\n"), " ")
})
# unlist pos tagged elements
icle$iclepos <- unlist(iclepos)
###############################################################
# extract number of adjs per line
pstggd <- icle$iclepos
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
concjjicle <- ConcR(icle$iclepos, pattern, context, all.pre = FALSE)
# repeat rows in data frame as often as there are adjectives in it (if 0 adj, repeat once)
corpusicleadjdf <- icle[rep(seq(nrow(icle)), rp),]
# combine data sets
corpusicleadj <- data.frame(1:nrow(corpusicleadjdf), corpusicleadjdf, concjjicle)
# remove rows without Tokens
ampicle <- corpusicleadj[is.na(corpusicleadj$Token) == F,]
# add clean column names
colnames(ampicle) <- c("ID", "File", "Language", "Subfile", "Speaker", "Content", "ContentClean",
                       "PosTaggedContent", "OriginalString", "PreContext", "Adjective", "PostContext")
# clean adjectives
ampicle$Adjective <- gsub(".*/JJR", "remove", ampicle$Adjective)
ampicle$Adjective <- gsub(".*/JJS", "remove", ampicle$Adjective)
amicle <- ampicle[!ampicle$Adjective == "remove",]
ampicle$Adjective <- str_replace_all(ampicle$Adjective, fixed("/JJ"), "")
# add Vraiant column
ampicle$Variant <- gsub(".* ", "", str_trim(ampicle$PreContext, side = "both")) 
# inspect data
nrow(ampicle); head(ampicle)

###############################################################
# define amplifiers
amplifiers <- c("absolutely", #"actually", "aggressively", 
                #"amazingly", "appallingly", 
                "awful", "awfully", 
                #"badly", 
                "bloody", "certainly", "clearly",
                #"complete", 
                "dead", "completely", 
                #"considerably", 
                "crazy", "decidedly", "definitely",  "distinctly", 
                "dreadfully", "enormously", "entirely", #"especially", 
                #"exactly", 
                "exceedingly", "exceptionally", 
                "excruciatingly", "extraordinarily", "extremely",
                "fiercely", "firmly", "frightfully", "fucking", 
                "fully", "genuinely", "greatly",
                "grossly", "heavily", "highly", "hopelessly", 
                #"horrendously", "hugely",
                #"immediately", "immensely", 
                "incredibly", 
                #"infinitely", "intensely", "irrevocably",
                "mad", "mega", "mighty", #"most", "much", 
                "obviously", "openly", "overwhelmingly", "particularly", 
                "perfectly", "plenty", "positively", "precisely", 
                "pretty", "profoundly", "purely", 
                #"quite", 
                "real", "really", "remarkably", "seriously", 
                #"shocking",   
                "significant", "significantly", "so", 
                #"specially", 
                "specifically", "strikingly",
                "strongly", "substantially", "super", "surely", 
                "terribly", "terrifically", 
                #"too",
                "total", "totally", "traditionally", #"true", 
                "truly", "ultra", "utterly", "very",
                "viciously", 
                #"well", 
                "wholly", "wicked", "wildly")
# clean ice icle data
ampicle$Function <- str_trim(ampicle$PostContext, side = "both")
ampicle$Function <- tolower(ampicle$Function)
ampicle$Function <- gsub(" {2,}", " ", ampicle$Function)
ampicle$Function <- sub(" ", "", ampicle$Function)
ampicle$Function <- gsub(" .*", "", ampicle$Function)
ampicle$Function <- gsub(".*/nn.*", "Attributive", ampicle$Function)
ampicle$Function <- ifelse(ampicle$Function == "Attributive", ampicle$Function, "Predicative")
# language
ampicle$Language <- ifelse(ampicle$Language == "GE", "German", ampicle$Language)
ampicle$Language <- ifelse(ampicle$Language == "BG", "Bulgarian", ampicle$Language)
ampicle$Language <- ifelse(ampicle$Language == "CZ", "Czech", ampicle$Language)
ampicle$Language <- ifelse(ampicle$Language == "DB", "Flemish", ampicle$Language)
ampicle$Language <- ifelse(ampicle$Language == "DN", "Dutch", ampicle$Language)
ampicle$Language <- ifelse(ampicle$Language == "FIN", "Finnish", ampicle$Language)
ampicle$Language <- ifelse(ampicle$Language == "FR", "French", ampicle$Language)
ampicle$Language <- ifelse(ampicle$Language == "IT", "Italian", ampicle$Language)
ampicle$Language <- ifelse(ampicle$Language == "PO", "Polish", ampicle$Language)
ampicle$Language <- ifelse(ampicle$Language == "SP", "Spanish", ampicle$Language)
ampicle$Language <- ifelse(ampicle$Language == "SW", "Swedish", ampicle$Language)
ampicle$Language <- ifelse(ampicle$Language == "FSW", "Swedish", ampicle$Language)
ampicle$Language <- ifelse(ampicle$Language == "RU", "Russian", ampicle$Language)
ampicle$Language <- ifelse(ampicle$Language == "MOSC", "Russian", ampicle$Language)
languages <- c("German", "Bulgarian", "Czech", "Flemish", "Dutch", "Finnish", 
               "French", "Italian", "Polish", "Spanish", "Swedish", "Russian")
ampicle$Language <- ifelse(ampicle$Language %in% languages, ampicle$Language, "remove")
ampicle <- ampicle[ampicle$Language != "remove",]
# shorten post Context
ampicle$PostContext <- substr(ampicle$PostContext, 1, ifelse((nchar(ampicle$PostContext)+25) <25, maampicle(nchar(ampicle$PostContext)), 25))
# pre Context
ampicle$PreContext <- str_trim(ampicle$PreContext, side = "both")
ampicle$PreContextLong <- ampicle$PreContext
ampicle$PreContextLong <- substr(ampicle$PreContextLong, ifelse(nchar(ampicle$PreContextLong)-25 <=0, 1, 
                                                                nchar(ampicle$PreContextLong)-25), nchar(ampicle$PreContextLong))
ampicle$PreContext <- gsub(".* ", "", ampicle$PreContext)
# amplifier variant
ampicle$PreContext <- gsub("\\/.*", "", ampicle$PreContext)
ampicle$Variant <- ifelse(ampicle$PreContext %in% amplifiers, ampicle$PreContext, "0")
# amplified y/n
ampicle$Amplified <- ifelse(ampicle$Variant == "0", 0, 1) 
# adjective
ampicle$Adjective <- tolower(ampicle$Adjective)
# inspect data
nrow(ampicle); head(ampicle); table(ampicle$Variant)

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
str(ampicle); head(ampicle); nrow(ampicle)#; table(ampicle$pint); head(ampicle$PreContextLong); head(ampicle$PreContextLong)

# code priming
prim1 <- c(rep(0, 1), ampicle$Variant[1:length(ampicle$Variant)-1])
prim2 <- c(rep(0, 2), ampicle$Variant[1:(length(ampicle$Variant)-2)])
prim3 <- c(rep(0, 3), ampicle$Variant[1:(length(ampicle$Variant)-3)])
primtb <- cbind(ampicle$Variant, prim1, prim2, prim3)

ampicle$Priming <- as.vector(unlist(apply(primtb, 1, function(x){
  x <- ifelse(x[1]== "0" , "noprime",
              ifelse(x[1] == x[2] | x[1] == x[3] | x[1] == x[4], "prime", "noprime"))
})))
# remove items that were not intensified by a minimum of 2 intensifier variants
nrow(ampicle)

# find items to be removed
supsidx <- unique(grep(paste(sups,collapse="|"), ampicle$PreContextLong, value=F))
negsidx <- unique(grep(paste(negs,collapse="|"), ampicle$PreContextLong, value=F))
downtonersidx <- unique(grep(paste(downtoners,collapse="|"), ampicle$PreContextLong, value=F))
specialformsidx <- unique(grep(paste(specialforms,collapse="|"), ampicle$PreContextLong, value=F))
PostContextdowntonersidx <- unique(grep(paste(PostContextdowntoners,collapse="|"), ampicle$PostContext, value=F))
nonpropadjidx <- unique(grep(paste(nonpropadj,collapse="|"), ampicle$Adjective, value=F))
# combine indices
idxs <- unique(c(supsidx, negsidx, downtonersidx, specialformsidx, PostContextdowntonersidx, nonpropadjidx))
# remove forms that require removal
ampicle <- ampicle[-idxs,]
# remove empty values
ampicle <- ampicle[!ampicle$Variant == "", ]
# remove superfluous columns
colnames(ampicle)
ampicle$Content <- NULL
ampicle$ContentClean <- NULL
ampicle$PosTaggedContent <- NULL
ampicle$OriginalString <- NULL
###############################################################
# save raw data to disc
write.table(ampicle, "ampicle_woneg02.txt", sep = "\t", row.names = F)
###############################################################

pintadjtb <- table(ampicle$Adjective, ampicle$Variant)
#pintadjtb <- pintadjtb[2:nrow(pintadjtb),]
pintadjtb <- pintadjtb[,2:ncol(pintadjtb)]
pintadjtb2 <- apply(pintadjtb, 1, function(x){
  x <- ifelse(x > 1, 1, x)})
pintadjtb3 <- colSums(pintadjtb2)
pintadjschildes <- names(pintadjtb3)[which(pintadjtb3 >=2 )]
ampicle <- ampicle[ampicle$Adjective %in% pintadjschildes, ]
nrow(ampicle)

# inspect adjectives
names(table(ampicle$Adjective))

# create vector with false adjectives
rmvadj <- c("uhr")
ampicle$remove <- ifelse(ampicle$Adjective %in% rmvadj, "remove", ampicle$Adjective)
ampicle <- ampicle[ampicle$remove != "remove",]
# remove superfluous columns
ampicle$remove <- NULL
colnames(ampicle)

# inspecta data
nrow(ampicle); length(table(ampicle$Adjective)); head(ampicle)

###############################################################
# save raw data to disc
write.table(ampicle, "ampicle_semiclean03.txt", sep = "\t", row.names = F)
###############################################################
# tabluate amplifiers
table(ampicle$Variant)[order(table(ampicle$Variant), decreasing = T)]

# create columns with intensifier freqs
ampicle$very <- ifelse(ampicle$Variant == "very", 1, 0) 
ampicle$so <- ifelse(ampicle$Variant == "so", 1, 0) 
ampicle$really <- ifelse(ampicle$Variant == "really", 1, 0) 
ampicle$extremely <- ifelse(ampicle$Variant == "extremely", 1, 0)
ampicle$other <- ifelse(ampicle$Variant == "extremely" |
                          ampicle$Variant == "really" | 
                          ampicle$Variant == "so" |
                          ampicle$Variant == "very", 1, 0)
# inspect results
head(ampicle)

###############################################################
# code freq of adj type by age group
frqadjtb <- table(ampicle$Language, ampicle$Adjective)
relfreqadjtb <- round(prop.table(frqadjtb, margin = 1)*100, 5)
relfreqadjdf <- as.data.frame(relfreqadjtb)
colnames(relfreqadjdf)[1:2] <- c("Language", "Adjective")
# add freq by date to data
ampicle <- merge(ampicle, relfreqadjdf, by=c("Language", "Adjective"))
# reorder data
ampicle <- ampicle[order(ampicle$ID),]
# inspect data
head(ampicle)

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
                "artificial", "artistic", "asian", "asthmatic", "atmospheric", "atomic", "iclesie",
                "icletralian", "icletrian", "authorised", "automatic", "autonomous", "average", "awake",
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
alladj <- names(table(ampicle$Adjective)) 
grd1 <- names(table(ampicle$Adjective[ampicle$Variant == "very"])) # adjs with very are gradable
grd2 <- names(table(ampicle$Adjective[ampicle$Variant == "extremely"])) # adjs with extremely are gradable
ngrd1 <- names(table(ampicle$Adjective[ampicle$Variant == "completely"])) # adjs with completely are gradable
ngrd2 <- names(table(ampicle$Adjective[ampicle$Variant == "total"])) # adjs with total are gradable
ngrd3 <- names(table(ampicle$Adjective[ampicle$Variant == "totally"])) # adjs with totally are gradable
ngrd4 <- names(table(ampicle$Adjective[ampicle$Variant == "utterly"])) # adjs with utterly are gradable
ngrd5 <- names(table(ampicle$Adjective[ampicle$Variant == "absolutely"])) # adjs with absolutely are gradable
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
grdicle <- names(table(c(grd_datadriven, grd_add)))
ngrdicle <- names(table(c(ngrd_datadriven, ngrd_add)))
# check which adjs are still unclassified
nagradicle1 <- alladj[!alladj %in% grdicle]
nagradicle <- nagradicle1[!nagradicle1 %in% ngrdicle]
# inspect unclassified adj
nagradicle

# inspect length of vectors
length(alladj); length(grdicle); length(ngrdicle); length(nagradicle)

# add gradability coding
ampicle$Gradabilty <- ifelse(ampicle$Adjective %in% grdicle, "Gradable", ampicle$Adjective)
ampicle$Gradabilty <- ifelse(ampicle$Gradabilty %in% ngrdicle, "NotGradable", ampicle$Gradabilty)
ampicle$Gradabilty <- ifelse(ampicle$Gradabilty  == "Gradable" | 
                               ampicle$Gradabilty  == "NotGradable", 
                             ampicle$Gradabilty, "GradabilityUndetermined")
# inspect data
head(ampicle); table(ampicle$Gradabilty)

# save table of gradable and non-gradable adj
gradtb <- table(ampicle$Adjective, ampicle$Gradabilty)
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
ampicle$SemanticCategory <- ifelse(ampicle$Adjective %in% semage, "Age", ampicle$Adjective)
ampicle$SemanticCategory <- ifelse(ampicle$SemanticCategory %in% semcol, "Color", ampicle$SemanticCategory)
ampicle$SemanticCategory <- ifelse(ampicle$SemanticCategory %in% semdif, "Difficulty", ampicle$SemanticCategory)
ampicle$SemanticCategory <- ifelse(ampicle$SemanticCategory %in% semdim, "Dimension", ampicle$SemanticCategory)
ampicle$SemanticCategory <- ifelse(ampicle$SemanticCategory %in% semhup, "HumanPropensity", ampicle$SemanticCategory)
ampicle$SemanticCategory <- ifelse(ampicle$SemanticCategory %in% semphy, "PhysicalProperty", ampicle$SemanticCategory)
ampicle$SemanticCategory <- ifelse(ampicle$SemanticCategory %in% semval, "Value", ampicle$SemanticCategory)
ampicle$SemanticCategory <- ifelse(ampicle$SemanticCategory == "Age" | ampicle$SemanticCategory == "Color" | ampicle$SemanticCategory == "Difficulty" | 
                                     ampicle$SemanticCategory == "Dimension" | ampicle$SemanticCategory == "HumanPropensity" |
                                     ampicle$SemanticCategory == "PhysicalProperty" | ampicle$SemanticCategory == "Value",  ampicle$SemanticCategory, "NoSemType")
# table sem class of tokens
table(ampicle$SemanticCategory)

# check classification
names(table(ampicle$Adjective[which(ampicle$SemanticCategory == "NoSemType")]))

# inspect data
head(ampicle); nrow(ampicle); length(table(ampicle$Adjective))

###############################################################
# load library
library(syuzhet)
# code emotion
class_emo <- get_nrc_sentiment(ampicle$Adjective)
# process sentiment
ampicle$Emotionality <- as.vector(unlist(apply(class_emo, 1, function(x){
  x <- ifelse(x[9] == 1, "NegativeEmotional",
              ifelse(x[10] == 1, "PositiveEmotional", "NonEmotional")) } )))
# revert order of factor Emotionality
ampicle$Emotionality <- factor(ampicle$Emotionality, levels = c("NonEmotional", "NegativeEmotional", "PositiveEmotional"))
# inspect data
head(ampicle); str(ampicle); nrow(ampicle)

###############################################################
# save raw data to disc
write.table(ampicle, "ampicle_clean04.txt", sep = "\t", row.names = F)
###############################################################
#                        END PART 1
###############################################################
