###################################################################################################################
##### PMMST: Phenotype Data Processing Clinic Visit 2 Cognitive  ##################################################
###################################################################################################################

# 13/07/2017

# load packages
library(readxl)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(GenABEL)
library(haven)

# Data were sent by Ayden on 31/05/2017. 
# File is in a zipped folder. Folder name is Emphasis_GMB_CV_data.zip. File name is
# Emphasis_CV2_Cognitive_170531 (excel file)
# On June 7, 2017 queries were sent to Ayden. Ayden came back with answers to most of those queries on
# 30/60/2017. A new version of the data was sent (Emphasis_CV2_Cognitive_170630.xlsx)

# Each cognitive component was reported in different excel spreadsheet.

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

# load dataset
datAtlantis <- read_excel("./PMMST Outcome Data/Emphasis_CV2_Cognitive_170630.xlsx", 
                          sheet = "Atlantis", na = "null", skip = 1, col_names = TRUE)

# change every ID number in upper case
datAtlantis$cISubjectID <- toupper(datAtlantis$cISubjectID)

# sanity check on imported dataset
dim(datAtlantis)
# last row contains NA only (excel imported an extra row). Remove it
datAtlantis <- datAtlantis[-nrow(datAtlantis),]

str(datAtlantis)
# use data dictionary sent by Ayden on 31/05/2017 to rename variables and levels 
# cognitive function were taken at visit 2. Check all the subject have this recorded
# as visit 2.
table(datAtlantis$cVisNo)
# Check all subject have a visit date
addmargins(table(datAtlantis$dVisdate))
# visit date format: YYYY-MM-DD

# transform school attended in a factor. Levels are as following:
# 1. English
# 2. Arabic
# 3. Both
# 4. None
addmargins(table(datAtlantis$cSchoolAttended))
datAtlantis$cSchoolAttended <- factor(datAtlantis$cSchoolAttended,
                                      labels=c("English","Arabic","Both", "None"))
addmargins(table(datAtlantis$cSchoolAttended))

# tTime is the time of the visit. Once imported in R the date 1899-12-31 is added.
# As this is clearly a wrong date eliminate and keep only time, and merge the time with
# the correct date
TimeOnly <- as.character(datAtlantis$tTime)
TimeOnly <- unlist(strsplit(TimeOnly, " "))[seq(2, length(TimeOnly), by = 2)]
datAtlantis$tTime <- as.POSIXct(paste(datAtlantis$dVisdate, TimeOnly), tz = "GMT")
# remove time only variable
rm(TimeOnly)

# check score imported correctly using the score sum provided in the Excel sheet
sum(datAtlantis$nCumulative_score)
sum(datAtlantis$nRaw_Score, na.rm = T)
sum(datAtlantis$nAtlantis_score, na.rm = T)

# check missing Atlantis score
sum(is.na(datAtlantis$nAtlantis_score))
# 1 subject is missing Atlantis score
datAtlantis[is.na(datAtlantis$nAtlantis_score),"cISubjectID", ]
# subject number EMP2851

# take only the relevant variables
datAtlantis <- datAtlantis %>% dplyr::select(cISubjectID, cVisNo, 
                                             dVisdate, cSchoolAttended, 
                                             nAtlantis_score)

hist(datAtlantis$nAtlantis_score)
summary(datAtlantis$nAtlantis_score)

# check for outliers
cutoffU <- median(datAtlantis$nAtlantis_score, na.rm = T) + 5*mad(datAtlantis$nAtlantis_score, na.rm = T)
cutoffL <- median(datAtlantis$nAtlantis_score, na.rm = T) - 5*mad(datAtlantis$nAtlantis_score, na.rm = T)
sum(datAtlantis$nAtlantis_score > cutoffU | datAtlantis$nAtlantis_score < cutoffL, na.rm = T)
# no outliers were found
rm(cutoffU, cutoffL)

ggplot(datAtlantis, aes(x = nAtlantis_score)) + geom_histogram(colour="black", fill="white", binwidth = 4) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 35)) + ggtitle("PMMST: Atlantis Score") +
  xlab("Atlantis Score")

###################################################################################################################
#### CODING #######################################################################################################
###################################################################################################################

# load dataset
datCoding <- read_excel("./PMMST Outcome Data/Emphasis_CV2_Cognitive_170630.xlsx", sheet = "Coding", na = "NULL",
                        col_names = TRUE)
# sanity check on imported dataset
dim(datCoding)

str(datCoding)
# Check all the subject have this recorded as visit 2.
table(datCoding$dVisNo)
# Check all subject have a visit date
addmargins(table(datCoding$dvisDate))
# visit date format: YYYY-MM-DD

# change every ID number in upper case
datCoding$cISubjectID <- toupper(datCoding$cISubjectID)

# transform all complessby variable in uppercase for consistency among levels
datCoding$cCheck_complenessby <- toupper(datCoding$cCheck_complenessby)
table(datCoding$cCheck_complenessby)

# rename cComment so that is clear that comment refer to the coding part only
names(datCoding)[names(datCoding) == "cComment"] <- "CodingComment"

# Time taken to complete the task
summary(datCoding$tTimetaken)
# check the missing Time for possible reasons
datCoding %>% filter(is.na(tTimetaken)) %>% select(cISubjectID, CodingComment)
# children missing cannot write

# take only the relevant variables
datCoding <- datCoding %>% 
  select(cISubjectID, nNumber_of_correctitem, CodingComment)
head(datCoding)

summary(datCoding$nNumber_of_correctitem)
hist(datCoding$nNumber_of_correctitem)

# check for outliers 
cutoffU <- median(datCoding$nNumber_of_correctitem, na.rm = T) + 5*mad(datCoding$nNumber_of_correctitem, na.rm = T)
cutoffL <- median(datCoding$nNumber_of_correctitem, na.rm = T) - 5*mad(datCoding$nNumber_of_correctitem, na.rm = T)
sum(datCoding$nNumber_of_correctitem > cutoffU | datCoding$nNumber_of_correctitem < cutoffL, na.rm = T)
# no outliers were found
rm(cutoffU, cutoffL)

ggplot(datCoding, aes(x = nNumber_of_correctitem)) + geom_histogram(colour="black", fill="white", binwidth = 4) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 65)) + ggtitle("PMMST: Coding") +
  xlab("Coding")


# Join Atlantis and Coding
# check if coding and atalntis agree on subject ID and visit date
datAC <- inner_join(datAtlantis, datCoding)
# remove single dataset
rm(datAtlantis, datCoding)
# look at the join dataset
dim(datAC)


###################################################################################################################
#### KOH'S BLOCKS #################################################################################################
###################################################################################################################

# load dataset
datKoh <- read_excel("./PMMST Outcome Data/Emphasis_CV2_Cognitive_170630.xlsx", sheet = "TKOH_Block", na = "NULL",
                        col_names = TRUE)
# sanity check on imported dataset
dim(datKoh)

str(datKoh)
# change every ID number in upper case
datKoh$cISubjectID <- toupper(datKoh$cISubjectID)
# Check all the subject have this recorded as visit 2.
table(datKoh$cVisNo)
# Check all the subject have a visit date recorded
addmargins(table(datKoh$dVisdate))
# visit date format: YYYY-MM-DD

# Whenever a time variable is in the table, R adds a wrong date (1899-12-31). Need to delete the date
# check all the variables that include time
# transform time variables in character.
datKoh <- lapply(datKoh, function(x) if(inherits(x, "POSIXct")) as.character(x) else x)
datKoh <- do.call(cbind.data.frame, datKoh)
# take the column we need to strip the date
cols <- names(datKoh)[grep("minute", names(datKoh))]
for(col in cols){
  # for each colum strip the incorrect date
  s1 <- sapply(strsplit(as.character(datKoh[,col]), split=" "), function(x) (x[2]))
  # take the seconds
  s1 <- sapply(strsplit(s1, split=":"), function(x) (x[2]))
  # if seconds = "00" then the children took a minute, otherwise second/60 (round to 3 digits)
  s1 <- round(ifelse(s1 == "00", 1, as.numeric(s1)/60), 3)
  # substitute date/time variable with time only
  datKoh[,col] <- s1
}
# remove unses vectors
rm(col); rm(cols); rm(s1)
# re-transform visit date in POSIXct.
datKoh$dVisdate <- as.POSIXct(datKoh$dVisdate) 
str(datKoh)
# look at the field workers
table(datKoh$cFW)

# check missing Koh score
sum(is.na(datKoh$nIQ))
# no subject are missing koh score

# get only the variables of interest
datKoh <- datKoh %>% dplyr::select(cISubjectID, dVisdate, nIQ)
head(datKoh)

summary(datKoh$nIQ)
hist(datKoh$nIQ)

# check for outliers 
cutoffU <- median(datKoh$nIQ, na.rm = T) + 5*mad(datKoh$nIQ, na.rm = T)
cutoffL <- median(datKoh$nIQ, na.rm = T) - 5*mad(datKoh$nIQ, na.rm = T)
sum(datKoh$nIQ > cutoffU | datKoh$nIQ < cutoffL, na.rm = T)
# 16 were identified as outliers based on the criteria estabished in the analysis plan
datKoh[datKoh$nIQ > cutoffU, ]

# create a flag variable for Koh's IQ
datKoh$flagKoh <- factor(ifelse(datKoh$nIQ > cutoffU, "Outlier", "No Outlier"))
table(datKoh$flagKoh)
str(datKoh)
rm(cutoffU, cutoffL)

# summary statistics (no outliers removed)
summary(datKoh$nIQ)
ggplot(datKoh, aes(x = nIQ)) + geom_histogram(colour="black", fill="white", binwidth = 6) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 150)) + ggtitle("PMMST: Koh's Block Design") +
  xlab("Koh's Block Design")
# after removing the outliers
datKoh %>% filter(flagKoh == "No Outlier") %>% 
  summarise(n = n(), Min = min(nIQ), Mean = mean(nIQ), SD = 
              sd(nIQ), Median = median(nIQ), Q1 = 
              quantile(nIQ, probs = 0.25), Q3 = 
              quantile(nIQ, probs = 0.75), Max = max(nIQ))
ggplot(subset(datKoh, flagKoh == "No Outlier"), aes(x = nIQ)) + 
  geom_histogram(colour="black", fill="white", binwidth = 1) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 40)) + 
  ggtitle("PMMST (No Outliers): Koh's Block Design") +
  xlab("Koh's Block Design")

# merge Koh dataset with Atlantis and Coding
datKoh$cISubjectID <- as.character(datKoh$cISubjectID)
names(datKoh)[names(datKoh) == "dVisdate"] <- "Visit"
datACK <- inner_join(datAC, datKoh)
# check visit dates are the same
sum(as.character(datACK$dVisdate) == as.character(datACK$Visit))
# drop one visit date
datACK <- subset(datACK, select = - Visit)
head(datACK)
names(datACK)
# remove not useful data
rm(datAC, datKoh)


###################################################################################################################
#### PATTERN REASONING ############################################################################################
###################################################################################################################

# load dataset
datPattern <- read_excel("./PMMST Outcome Data/Emphasis_CV2_Cognitive_170630.xlsx", 
                         sheet = "TPattern_reasoning", na = "NULL", col_names = TRUE)

# sanity check on imported dataset
dim(datPattern)
str(datPattern)
# Check all the subject have this recordedas visit 2.
table(datPattern$cVisNo)
# Check all subject have a visit date
addmargins(table(datPattern$dVisdate))
# visit date format: YYYY-MM-DD

# check missing Pattern score
sum(is.na(datPattern$nCommulative_score))
# no subject has missing pattern score

# change every ID number in upper case
datPattern$cISubjectID <- toupper(datPattern$cISubjectID)
names(datPattern)

# get only the variables of interest
datPattern <- datPattern %>% dplyr::select(cISubjectID, dVisdate, nCommulative_score)
head(datPattern)

summary(datPattern$nCommulative_score)
hist(datPattern$nCommulative_score)

# check for outliers 
cutoffU <- median(datPattern$nCommulative_score, na.rm = T) + 5*mad(datPattern$nCommulative_score, na.rm = T)
cutoffL <- median(datPattern$nCommulative_score, na.rm = T) - 5*mad(datPattern$nCommulative_score, na.rm = T)
sum(datPattern$nCommulative_score > cutoffU | datPattern$nCommulative_score < cutoffL, na.rm = T)
# 2 were identified as outliers based on the criteria estabished in the analysis plan
datPattern[datPattern$nCommulative_score > cutoffU, ]

# create a flag variable for pattern score
datPattern$flagPatternScore <- factor(ifelse(datPattern$nCommulative_score > cutoffU, "Outlier", "No Outlier"))
table(datPattern$flagPatternScore)
str(datPattern)
rm(cutoffU, cutoffL)

# summary statistics (no outliers removed)
summary(datPattern$nCommulative_score)
ggplot(datPattern, aes(x = nCommulative_score)) + geom_histogram(colour="black", fill="white", binwidth = 1) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 65)) + ggtitle("PMMST: Pattern Score") +
  xlab("Pattern Score")
# after removing the outliers
datPattern %>% filter(flagPatternScore == "No Outlier") %>% 
  summarise(n = n(), Min = min(nCommulative_score), Mean = mean(nCommulative_score), SD = 
              sd(nCommulative_score), Median = median(nCommulative_score), Q1 = 
              quantile(nCommulative_score, probs = 0.25), Q3 = 
              quantile(nCommulative_score, probs = 0.75), Max = max(nCommulative_score))
ggplot(subset(datPattern, flagPatternScore == "No Outlier"), aes(x = nCommulative_score)) + 
  geom_histogram(colour="black", fill="white", binwidth = 1) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 65)) + 
  ggtitle("PMMST (No Outliers): Pattern Score") +
  xlab("Pattern Score")

# merge Pattern dataset with Atlantis, Coding and Koh
names(datPattern)[names(datPattern) == "dVisdate"] <- "Visit"
datACKP <- inner_join(datACK, datPattern)
# check visit dates are the same
sum(as.character(datACKP$dVisdate) == as.character(datACKP$Visit))
# check the one that is different
datACKP[as.character(datACKP$dVisdate) != as.character(datACKP$Visit), c("cISubjectID", "dVisdate", "Visit")]
# one subject has to different dates. Correct date is 2016-10-07 (dVisdate)
# drop one visit date
datACKP <- subset(datACKP, select = - Visit)
head(datACKP)
names(datACKP)
# remove not useful data
rm(datACK, datPattern)


###################################################################################################################
#### VERBAL FLUENCY (ANIMAL) ######################################################################################
###################################################################################################################

# load dataset
datVFAnimal <- read_excel("./PMMST Outcome Data/Emphasis_CV2_Cognitive_170630.xlsx", 
                         sheet = "TVerbalFluency_animals", na = "NULL", col_names = TRUE)

# change every ID number in upper case
datVFAnimal$cISubjectID <- toupper(datVFAnimal$cISubjectID)

# sanity check on imported dataset
dim(datVFAnimal)
str(datVFAnimal)
# Check all the subject have this recorded as visit 2.
table(datVFAnimal$cVisNo)
# Check all subject have a visit date
addmargins(table(datVFAnimal$dVisitDate))
# visit date format: YYYY-MM-DD

# check missing Pattern score
sum(is.na(datVFAnimal$cTotalAnimal))
# no subject has missing verbal animal fluency

summary(datVFAnimal$cTotalAnimal)
hist(datVFAnimal$cTotalAnimal)

# check for outliers 
cutoffU <- median(datVFAnimal$cTotalAnimal, na.rm = T) + 5*mad(datVFAnimal$cTotalAnimal, na.rm = T)
cutoffL <- median(datVFAnimal$cTotalAnimal, na.rm = T) - 5*mad(datVFAnimal$cTotalAnimal, na.rm = T)
sum(datVFAnimal$cTotalAnimal > cutoffU | datVFAnimal$cTotalAnimal < cutoffL, na.rm = T)
# 0 were identified as outliers
rm(cutoffU, cutoffL)

# take only the variables needed
datVFAnimal <- datVFAnimal %>% dplyr::select(cISubjectID, dVisitDate, cTotalAnimal)

summary(datVFAnimal$cTotalAnimal)
ggplot(datVFAnimal, aes(x = cTotalAnimal)) + geom_histogram(colour="black", fill="white", binwidth = 1) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 60)) + ggtitle("PMMST: Verbal Fluency - Animals") +
  xlab("Verbal Fluency - Animals")


# merge verbal fluency (animal) dataset with Atlantis, Coding, Koh and Pattern
names(datVFAnimal)[names(datVFAnimal) == "dVisitDate"] <- "Visit"
datACKPA <- inner_join(datACKP, datVFAnimal)
# check visit dates are the same
sum(as.character(datACKPA$dVisdate) == as.character(datACKPA$Visit))
# drop one visit date
datACKPA <- subset(datACKPA, select = - Visit)
head(datACKPA)
names(datACKPA)
# remove not useful data
rm(datACKP, datVFAnimal)



###################################################################################################################
#### VERBAL FLUENCY (NAMES) #######################################################################################
###################################################################################################################

# load dataset
datVFNames <- read_excel("./PMMST Outcome Data/Emphasis_CV2_Cognitive_170630.xlsx", 
                          sheet = "TVerbalFluency_names", na = "NULL", col_names = TRUE)

# change every ID number in upper case
datVFNames$cISubjectID <- toupper(datVFNames$cISubjectID)

# sanity check on imported dataset
dim(datVFNames)
str(datVFNames)
# Check all the subject have this recorded as visit 2.
table(datVFNames$cVisNo)
# Check all subject have a visit date
addmargins(table(datVFNames$dVisitDate))
# visit date format: YYYY-MM-DD

# check missing Pattern score
sum(is.na(datVFNames$cTotalName))
# no subject has missing verbal names fluency

summary(datVFNames$cTotalName)
hist(datVFNames$cTotalName)

# check for outliers 
cutoffU <- median(datVFNames$cTotalName, na.rm = T) + 5*mad(datVFNames$cTotalName, na.rm = T)
cutoffL <- median(datVFNames$cTotalName, na.rm = T) - 5*mad(datVFNames$cTotalName, na.rm = T)
sum(datVFNames$cTotalName > cutoffU | datVFNames$cTotalName < cutoffL, na.rm = T)
# 0 were identified as outliers
rm(cutoffU, cutoffL)

# take only the variables needed
datVFNames <- datVFNames %>% dplyr::select(cISubjectID, dVisitDate, cTotalName)

ggplot(datVFNames, aes(x = cTotalName)) + geom_histogram(colour="black", fill="white", binwidth = 1) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 50)) + ggtitle("PMMST: Verbal Fluency - Names") +
  xlab("Verbal Fluency - Names")

# merge verbal fluency (names) dataset with Atlantis, Coding, Koh and Pattern and verabl fluency animal
names(datVFNames)[names(datVFNames) == "dVisitDate"] <- "Visit"
datACKPAN <- inner_join(datACKPA, datVFNames)
# check visit dates are the same
sum(as.character(datACKPAN$dVisdate) == as.character(datACKPAN$Visit))
# drop one visit date
datACKPAN <- subset(datACKPAN, select = - Visit)
head(datACKPAN)
names(datACKPAN)
# remove not useful data
rm(datACKPA, datVFNames)



###################################################################################################################
#### WORD ORDER (NAMES) ##########################################################################################
###################################################################################################################

# load dataset
datWord <- read_excel("./PMMST Outcome Data/Emphasis_CV2_Cognitive_170630.xlsx", 
                         sheet = "TWordOrder", na = "NULL", col_names = TRUE)

# sanity check on imported dataset
dim(datWord)
str(datWord)
# Check all the subject have this recordedas visit 2.
table(datWord$cVisNo)
# Check all subject have a visit date
addmargins(table(datWord$dVisdate))
# visit date format: YYYY-MM-DD

# check missing word order
sum(is.na(datWord$cScore))
# no subject has missing word order

# change every ID number in upper case
datWord$cISubjectID <- toupper(datWord$cISubjectID)

summary(datWord$cScore)
hist(datWord$cScore)

# check for outliers 
cutoffU <- median(datWord$cScore, na.rm = T) + 5*mad(datWord$cScore, na.rm = T)
cutoffL <- median(datWord$cScore, na.rm = T) - 5*mad(datWord$cScore, na.rm = T)
sum(datWord$cScore > cutoffU | datWord$cScore < cutoffL, na.rm = T)
# 0 were identified as outliers
rm(cutoffU, cutoffL)

# take only the variables needed
datWord <- datWord %>% dplyr::select(cISubjectID, dVisdate, cScore)

ggplot(datWord, aes(x = cScore)) + geom_histogram(colour="black", fill="white", binwidth = 1.5) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 80)) + ggtitle("PMMST: Word Order") +
  xlab("Word Order")

# merge word order with all the others cognitive domain
names(datWord)[names(datWord) == "dVisdate"] <- "Visit"
dat <- inner_join(datACKPAN, datWord)
# check visit dates are the same
sum(as.character(dat$dVisdate) == as.character(dat$Visit))
# check the one that is different
dat[as.character(dat$dVisdate) != as.character(dat$Visit), c("cISubjectID", "dVisdate", "Visit")]
# 2016-10-04 is the correct visit date as confirmed by Ayden on June 30, 2017
# drop one visit date
dat <- subset(dat, select = - Visit)
head(dat)
names(dat)
# remove not useful data
rm(datACKPAN, datWord)



####################################################################################################################
#### Add Information on sex, age, dob and allocation group #########################################################
####################################################################################################################

# information on allocation group and updated sex information are in the file PMMST_pdata_redux.xlsx
datSexGroup <- read_excel("./PMMST Outcome Data/PMMST_pdata_redux_170718.xlsx", col_names = TRUE)
# select variables of interest 
datSexGroup <- datSexGroup %>% dplyr::select(cSubjectID, DateOfBirth, Sex, MasterGroupNo, Age)
# merge with cognitive function dataset
# change name of subject ID to agree with dat
names(datSexGroup)[names(datSexGroup) == "cSubjectID"] <- "cISubjectID" 
# transform everything to uppercase
datSexGroup$cISubjectID <- toupper(datSexGroup$cISubjectID)
sum(dat$cISubjectID %in% datSexGroup$cISubjectID)
# EMP180Q has only CV2 data add a comment
datSexGroup$GeneralComments <- NA
datSexGroup[datSexGroup$cISubjectID == "EMP180Q", "GeneralComments"] <- "Child has only visit 2 data"
dat <- inner_join(dat, datSexGroup)
rm(datSexGroup)
# drop not useful variable
dat <- subset(dat, select = - c(cVisNo, dvisDate, dVisNo, tTimetaken, cCheck_complenessby))
# compute age at visit 2 from dob and date of visit
dat$cageV2 <- with(dat, (dVisdate - DateOfBirth)/365.25)
dat$cageV2 <- as.numeric(dat$cageV2)
summary(dat$cageV2)



###################################################################################################################
#### DERIVED MEASURES #############################################################################################
###################################################################################################################

#### Compute z-scores ####

# average verbal fluency (animal and people)
cor(dat$cTotalAnimal, dat$cTotalName, method = "spearman")
with(dat, plot(cTotalAnimal, cTotalName))
dat$cVerbalFluency <- with(dat, (cTotalAnimal + cTotalName)/2)
summary(dat$cVerbalFluency)

# z-scores need to be age and sex adjusted
zscore <- function(x){rstandard(lm(x ~ as.factor(Sex) + cageV2, data = dat))}

# create variables
dat <- data.frame(dat, zcAtlantis = NA, zcCoding = NA)
zcAtlantis <- data.frame(zscore(dat$nAtlantis_score)) 
dat[rownames(dat) %in% rownames(zcAtlantis), "zcAtlantis"] <- zscore(dat$nAtlantis_score) 
zcCoding <- data.frame(zscore(dat$nNumber_of_correctitem)) 
dat[rownames(dat) %in% rownames(zcCoding), "zcCoding"] <- zscore(dat$nNumber_of_correctitem) 
dat$zcKoh <- zscore(dat$nIQ)
dat$zcPattern <- zscore(dat$nCommulative_score)
dat$zcVerbalFluency <- zscore(dat$cVerbalFluency)
dat$zcWordOrder <- zscore(dat$cScore)


# check summary and histogram
summary(dat$zcAtlantis)
summary(dat$zcCoding)
summary(dat$zcKoh)
summary(dat$zcPattern)
summary(dat$zcVerbalFluency)
summary(dat$zcWordOrder)

par(mfrow = c(2,3))
hist(dat$zcAtlantis, main = "Atlantis", xlab = "Atlantis")
hist(dat$zcCoding, main = "Coding", xlab = "Coding")
hist(dat$zcKoh, main = "Koh's Block", xlab = "Koh's Block")
hist(dat$zcPattern, main = "Pattern Reasoning", xlab = "Pattern Reasoning")
hist(dat$zcVerbalFluency, main = "Verbal Fluency", xlab = "Verbal Fluency")
hist(dat$zcWordOrder, main = "Word Order")

rm(zscore); rm(zcAtlantis); rm(zcCoding)

#### Compute mental processing index ####
# mental processing Index is the mean z-scores of the 6 cognitive domains.
zscores <- c("zcAtlantis", "zcCoding", "zcKoh", "zcPattern", "zcVerbalFluency", "zcWordOrder")
dat$cMPI <- rowMeans(dat[,names(dat) %in% zscores], na.rm = T)

with(dat, summary(cMPI))

ggplot(dat, aes(x = cMPI))  + geom_histogram(colour="black", fill="white", binwidth = 0.5) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) + ggtitle("PMMST: Mental Process Index") +
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
dat <- dat %>% dplyr::select(cISubjectID, Sex, Age, cageV2, MasterGroupNo, cSchoolAttended, nAtlantis_score, 
                      nNumber_of_correctitem, nIQ, nCommulative_score, 
                      cTotalAnimal, cTotalName,  cScore,  cVerbalFluency, zcAtlantis, zcCoding,
                      zcKoh, zcPattern, zcVerbalFluency, zcWordOrder, cMPI, flagKoh, flagPatternScore, 
                      CodingComment, GeneralComments)


names(dat) <- c("serno", "csex", "cageV1", "cageV2", "allocation", "cSchoolAttended", 
                "cAtlantis", "cCoding", "cKoh", "cPattern", "cVFAnimals", "cVFNames", "cWordOrder", 
                "cVerbalFluency", "zcAtlantis", "zcCoding",
                "zcKoh", "zcPattern", "zcVerbalFluency", "zcWordOrder", "cMPI", "flagcKoh", "flagcPattern", 
                "CodingComment", "GeneralComment")

dat$cageV1 <- as.numeric(dat$cageV1)

# add attributes
attributes(dat$serno)$label <- "pmmst - subject ID"
attributes(dat$csex)$label <- "pmmst - child sex"
attributes(dat$cageV1)$label <- "pmmst - child age at CV1 (years)"
attributes(dat$cageV2)$label <- "pmmst - child age at CV2 (years)"
attributes(dat$allocation)$label <- "pmmst - allocation group"
attributes(dat$cSchoolAttended)$label <- "pmmst - school attended by child"
attributes(dat$cAtlantis)$label <- "pmmst - Atlantis score"
attributes(dat$cCoding)$label <- "pmmst - coding score"
attributes(dat$cKoh)$label <- "pmmst - Koh's score (child IQ)"
attributes(dat$cPattern)$label <- "pmmst - pattern score"
attributes(dat$cVFAnimals)$label <- "pmmst - verbal fluency (animals)"
attributes(dat$cVFNames)$label <- "pmmst - verbal fluency (names)"
attributes(dat$cWordOrder)$label <- "pmmst - words order"
attributes(dat$cVerbalFluency)$label <- "pmmst - verbal fluency (average animals and names)"
attributes(dat$zcAtlantis)$label <- "pmmst - Atlantis z-score"
attributes(dat$zcCoding)$label <- "pmmst - coding z-score"
attributes(dat$zcKoh)$label <- "pmmst - Koh's z-score (child IQ)"
attributes(dat$zcPattern)$label <- "pmmst - pattern z-score"
attributes(dat$zcWordOrder)$label <- "pmmst - words order z-score"
attributes(dat$zcVerbalFluency)$label <- "pmmst - verbal fluency z-score"
attributes(dat$cMPI)$label <- "pmmst - mental processing index"
attributes(dat$flagcKoh)$label <- "pmmst - Koh's score outlier?"
attributes(dat$flagcPattern)$label <- "pmmst - pattern score outlier?"
attributes(dat$CodingComment)$label <- "pmmst - comments on coding score only"
attributes(dat$GeneralComment)$label <- "pmmst - comments"

# save as csv, stata and spss file
write.table(dat, "./Cleaned Dataset/PMMSTCognitivev1_20170719.csv", row.names = FALSE)
write_dta(dat, "./Cleaned Dataset/PMMSTCognitivev1_20170719.dta")
write_sav(dat, "./Cleaned Dataset/PMMSTCognitivev1_20170719.sav")
