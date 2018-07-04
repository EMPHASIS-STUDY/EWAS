###############################################################
##### SARAS Kids: Phenotype Data Processing Cognitive #########
###############################################################

# 03/08/2017

# load packages
library(haven)
library(dplyr)
library(ggplot2)

# Data were sent by Patsy on 05/06/2017. 
# File name is sarascog01.dta (STATA v14 file). 
# Only children in the per-protocol analysis are processed.
# sernos 981, 1536, 1939, 2587, 1548, 3156, 4250, 4500, 4513, 4580, 4701, 5218 
# Patsy and Vaness retrieved the data and sent a new version of the file sarascog02.dta
# serno 2587 has Down syndrom and does not have cognitive test taken

dat <- read_dta("sarascog02.dta")
# check for duplicates
length(unique(dat$serno))

# add information on allocation group, sex and age. PP analysis only
datList <- read_dta("SARASKidsBloodCountPP.dta")
datList <- datList %>% dplyr::select(serno, AnalysisType, allocation, age, sex)
# check if all sernos in datList are in dat
datList$serno[!(datList$serno %in% dat$serno)] # 12 sernos are not in dat
datList$allocation <- as_factor(datList$allocation)
datList$sex <- as_factor(datList$sex)


# merge the two dataset
dat <- inner_join(dat, datList)
rm(datList)
dim(dat) # 708 children has it is supposed to be

####################################################################################################################
#### Univariate Analysis ###########################################################################################
####################################################################################################################

# 1. Check each variable singularly
# 2. Check if same measures are consistent
# 3. Compute the mean of each observation
# 4. Check for outliers using the average computed in 3. and the criteria established in the analysis 
#    plan (median +/- 5MAD)
# 5. Retain the outliers but create a variable with 0/1 (named flagvariablename) where 1 indicated
#    the outlier based on the analysis plan.

###################################################################################################################
#### ANTLANTIS SCORE ##############################################################################################
###################################################################################################################

with(dat, summary(matlantis))
with(dat, sd(matlantis, na.rm = T))
dat[is.na(dat$matlantis),] 
# serno 1762 does not have any cognitive measure recorded except for Koh's block (flagged to Patsy).

hist(dat$matlantis)
# check for outliers
cutoffU <- median(dat$matlantis, na.rm = T) + 5*mad(dat$matlantis, na.rm = T)
cutoffL <- median(dat$matlantis, na.rm = T) - 5*mad(dat$matlantis, na.rm = T)
sum(dat$matlantis > cutoffU | dat$matlantis < cutoffL, na.rm = T)
# no outliers were found
rm(cutoffU, cutoffL)

ggplot(dat, aes(x = matlantis)) + geom_histogram(colour="black", fill="white", binwidth = 4) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 120)) + ggtitle("SARAS Kids: Atlantis Score") +
  xlab("Atlantis Score")

###################################################################################################################
#### CODING #######################################################################################################
###################################################################################################################

with(dat, summary(mwisc))
with(dat, sd(mwisc, na.rm = T))
dat[is.na(dat$mwisc),] 
# serno 1762 does not have any cognitive measure recorded except for Koh's block (flagged to Patsy).

hist(dat$mwisc)
# check for outliers
cutoffU <- median(dat$mwisc, na.rm = T) + 5*mad(dat$mwisc, na.rm = T)
cutoffL <- median(dat$mwisc, na.rm = T) - 5*mad(dat$mwisc, na.rm = T)
sum(dat$mwisc > cutoffU | dat$mwisc < cutoffL, na.rm = T)
# no outliers were found
rm(cutoffU, cutoffL)

ggplot(dat, aes(x = mwisc)) + geom_histogram(colour="black", fill="white", binwidth = 4) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 120)) + ggtitle("SARAS Kids: Coding") +
  xlab("Coding")

###################################################################################################################
#### KOH'S BLOCKS #################################################################################################
###################################################################################################################

with(dat, summary(mfinaliq))
with(dat, sd(mfinaliq, na.rm = T))

hist(dat$mfinaliq)
# check for outliers
cutoffU <- median(dat$mfinaliq, na.rm = T) + 5*mad(dat$mfinaliq, na.rm = T)
cutoffL <- median(dat$mfinaliq, na.rm = T) - 5*mad(dat$mfinaliq, na.rm = T)
sum(dat$mfinaliq > cutoffU | dat$mfinaliq< cutoffL, na.rm = T)
# 26 were identified as outliers based on the criteria estabished in the analysis plan
dat[dat$mfinaliq > cutoffU, ]
# create a flag variable for Koh's IQ
dat$flagmfinaliq <- factor(ifelse(dat$mfinaliq > cutoffU | dat$mfinaliq < cutoffL, "Outlier", "No Outlier"))
table(dat$flagmfinaliq)

# summary statistics (no outliers removed)
summary(dat$mfinaliq)
ggplot(dat, aes(x = mfinaliq)) + geom_histogram(colour="black", fill="white", binwidth = 6) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 220)) + ggtitle("SARAS Kids: Koh's Block Design") +
  xlab("Koh's Block Design")
# after removing the outliers
dat %>% filter(flagmfinaliq == "No Outlier") %>% 
  summarise(n = n(), Min = min(mfinaliq), Mean = mean(mfinaliq), SD = 
              sd(mfinaliq), Median = median(mfinaliq), Q1 = 
              quantile(mfinaliq, probs = 0.25), Q3 = 
              quantile(mfinaliq, probs = 0.75), Max = max(mfinaliq))
ggplot(subset(dat, flagmfinaliq == "No Outlier"), aes(x = mfinaliq)) + 
  geom_histogram(colour="black", fill="white", binwidth = 6) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 220)) + 
  ggtitle("SARAS Kids (No Outliers): Koh's Block Design") +
  xlab("Koh's Block Design")


###################################################################################################################
#### PATTERN REASONING ############################################################################################
###################################################################################################################

with(dat, summary(mpattern))
with(dat, sd(mpattern, na.rm = T))
dat[is.na(dat$mpattern),] 
# serno 1762 does not have any cognitive measure recorded except for Koh's block (flagged to Patsy).

hist(dat$mpattern)
# check for outliers
cutoffU <- median(dat$mpattern, na.rm = T) + 5*mad(dat$mpattern, na.rm = T)
cutoffL <- median(dat$mpattern, na.rm = T) - 5*mad(dat$mpattern, na.rm = T)
sum(dat$mpattern > cutoffU | dat$mpattern < cutoffL, na.rm = T)
# 0 outliers identified
rm(cutoffU, cutoffL)

ggplot(dat, aes(x = mpattern)) + geom_histogram(colour="black", fill="white", binwidth = 1) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 150)) + ggtitle("SARAS Kids: Pattern Score") +
  xlab("Pattern Score")


###################################################################################################################
#### VERBAL FLUENCY (ANIMAL) ######################################################################################
###################################################################################################################

with(dat, summary(mverbala))
with(dat, sd(mverbala, na.rm = T))
dat[is.na(dat$mverbala),] 
# serno 1762 does not have any cognitive measure recorded except for Koh's block (flagged to Patsy).

hist(dat$mverbala)
# check for outliers
cutoffU <- median(dat$mverbala, na.rm = T) + 5*mad(dat$mverbala, na.rm = T)
cutoffL <- median(dat$mverbala, na.rm = T) - 5*mad(dat$mverbala, na.rm = T)
sum(dat$mverbala > cutoffU | dat$mverbala < cutoffL, na.rm = T)
# 0 outliers identified
rm(cutoffU, cutoffL)

ggplot(dat, aes(x = mverbala)) + geom_histogram(colour="black", fill="white", binwidth = 1) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 120)) + ggtitle("SARAS Kids: Verbal Fluency - Animals") +
  xlab("Verbal Fluency - Animals")


###################################################################################################################
#### VERBAL FLUENCY (NAMES) #######################################################################################
###################################################################################################################

with(dat, summary(mverbalf))
with(dat, sd(mverbalf, na.rm = T))
dat[is.na(dat$mverbalf),] 
# serno 1762 does not have any cognitive measure recorded except for Koh's block (flagged to Patsy).

hist(dat$mverbalf)
# check for outliers
cutoffU <- median(dat$mverbalf, na.rm = T) + 5*mad(dat$mverbalf, na.rm = T)
cutoffL <- median(dat$mverbalf, na.rm = T) - 5*mad(dat$mverbalf, na.rm = T)
sum(dat$mverbalf > cutoffU | dat$mverbalf < cutoffL, na.rm = T)
# 0 outliers identified
rm(cutoffU, cutoffL)

ggplot(dat, aes(x = mverbalf)) + geom_histogram(colour="black", fill="white", binwidth = 1) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 120)) + ggtitle("SARAS Kids: Verbal Fluency - Names") +
  xlab("Verbal Fluency - Names")

###################################################################################################################
#### WORD ORDER (NAMES) ##########################################################################################
###################################################################################################################

with(dat, summary(mwordord))
with(dat, sd(mwordord, na.rm = T))
dat[is.na(dat$mwordord),] 
# serno 1762 does not have any cognitive measure recorded except for Koh's block (flagged to Patsy).

hist(dat$mwordord)
# check for outliers
cutoffU <- median(dat$mwordord, na.rm = T) + 5*mad(dat$mwordord, na.rm = T)
cutoffL <- median(dat$mwordord, na.rm = T) - 5*mad(dat$mwordord, na.rm = T)
sum(dat$mwordord > cutoffU | dat$mwordord < cutoffL, na.rm = T)
# 0 outliers identified
rm(cutoffU, cutoffL)

ggplot(dat, aes(x = mwordord)) + geom_histogram(colour="black", fill="white", binwidth = 1.2) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 160)) + ggtitle("SARAS Kids: Word Order") +
  xlab("Word Order")

###################################################################################################################
#### DERIVED MEASURES #############################################################################################
###################################################################################################################

#### Compute z-scores ####

# average verbal fluency (animal and people)
cor(dat$mverbala, dat$mverbalf, method = "spearman", use = "complete.obs")
with(dat, plot(mverbala, mverbalf))
dat$cVerbalFluency <- with(dat, (mverbala + mverbalf)/2)
summary(dat$cVerbalFluency)

# z-scores need to be age and sex adjusted
zscore <- function(x){rstandard(lm(x ~ as.factor(sex) + age, data = dat))}

# create variables
dat <- data.frame(dat, zmatlantis = NA, zmwisc = NA, zmpattern = NA, zcVerbalFluency = NA, zmwordord = NA)
zmatlantis <- data.frame(zscore(dat$matlantis)) 
dat[rownames(dat) %in% rownames(zmatlantis), "zmatlantis"] <- zscore(dat$matlantis) 
zmwisc <- data.frame(zscore(dat$mwisc)) 
dat[rownames(dat) %in% rownames(zmwisc), "zmwisc"] <- zscore(dat$mwisc) 
dat$zcKoh <- zscore(dat$mfinaliq)
zmpattern <- data.frame(zscore(dat$mpattern)) 
dat[rownames(dat) %in% rownames(zmpattern), "zmpattern"] <- zscore(dat$mpattern) 
zcVerbalFluency <- data.frame(zscore(dat$cVerbalFluency)) 
dat[rownames(dat) %in% rownames(zcVerbalFluency), "zcVerbalFluency"] <- zscore(dat$cVerbalFluency) 
zmwordord <- data.frame(zscore(dat$mwordord)) 
dat[rownames(dat) %in% rownames(zmwordord), "zmwordord"] <- zscore(dat$mwordord) 

rm(zmatlantis, zmwisc, zmpattern, zcVerbalFluency, zmwordord)

# check summary and histogram
summary(dat$zmatlantis)
summary(dat$zmwisc)
summary(dat$zcKoh)
summary(dat$zmpattern)
summary(dat$zcVerbalFluency)
summary(dat$zmwordord)

par(mfrow = c(2,3))
hist(dat$zmatlantis, main = "Atlantis", xlab = "Atlantis")
hist(dat$zmwisc, main = "Coding", xlab = "Coding")
hist(dat$zcKoh, main = "Koh's Block", xlab = "Koh's Block")
hist(dat$zmpattern, main = "Pattern Reasoning", xlab = "Pattern Reasoning")
hist(dat$zcVerbalFluency, main = "Verbal Fluency", xlab = "Verbal Fluency")
hist(dat$zmwordord, main = "Word Order")


#### Compute mental processing index ####
# mental processing Index is the mean z-scores of the 6 cognitive domains.
zscores <- c("zmatlantis", "zmwisc", "zcKoh", "zmpattern", "zcVerbalFluency", "zmwordord")
dat$cMPI <- rowMeans(dat[,names(dat) %in% zscores], na.rm = T)

with(dat, summary(cMPI))

ggplot(dat, aes(x = cMPI))  + geom_histogram(colour="black", fill="white", binwidth = 0.5) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 250)) + ggtitle("SARAS Kids: Mental Process Index") +
  xlab("Mental Process Index")

# check for outliers 
cutoffU <- median(dat$cMPI, na.rm = T) + 5*mad(dat$cMPI, na.rm = T)
cutoffL <- median(dat$cMPI, na.rm = T) - 5*mad(dat$cMPI, na.rm = T)
sum(dat$cMPI > cutoffU | dat$cMPI < cutoffL, na.rm = T)
# 0 were identified as outliers
rm(cutoffU, cutoffL)


###################################################################################################################
#### SAVE DATASET #################################################################################################
###################################################################################################################

names(dat)

# select only relevant variable (based on the analysis plan)
dat <- dat %>% dplyr::select(serno, sex, age, allocation, matlantis, mwisc, mfinaliq, mpattern, mverbala,
                             mverbalf, mwordord, cVerbalFluency, zmatlantis, zmwisc, zcKoh, zmpattern, 
                             zcVerbalFluency, zmwordord, cMPI, flagmfinaliq)


names(dat) <- c("serno", "csex", "cage", "allocation",  
                "cAtlantis", "cCoding", "cKoh", "cPattern", "cVFAnimals", "cVFNames", "cWordOrder", 
                "cVerbalFluency", "zcAtlantis", "zcCoding", "zcKoh", "zcPattern", "zcVerbalFluency", 
                "zcWordOrder", "cMPI", "flagcKoh")

attributes(dat$cKoh)$label <- "saras - Koh's score (child IQ)"
attributes(dat$zcAtlantis)$label <- "saras - Atlantis z-score"
attributes(dat$zcCoding)$label <- "saras - coding z-score"
attributes(dat$zcKoh)$label <- "saras - Koh's z-score (child IQ)"
attributes(dat$zcPattern)$label <- "saras - pattern z-score"
attributes(dat$zcWordOrder)$label <- "saras - words order z-score"
attributes(dat$zcVerbalFluency)$label <- "saras - verbal fluency z-score"
attributes(dat$cMPI)$label <- "saras - mental processing index"
attributes(dat$flagcKoh)$label <- "saras - Koh's score outlier?"


# save as csv, stata and spss file
write.table(dat, "./CleanedDataset/SARASCognitivev1_20170804.csv", row.names = FALSE)
write_dta(dat, "./CleanedDataset/SARASCognitivev1_20170804.dta")
write_sav(dat, "./CleanedDataset/SARASCognitivev1_20170804.sav")
