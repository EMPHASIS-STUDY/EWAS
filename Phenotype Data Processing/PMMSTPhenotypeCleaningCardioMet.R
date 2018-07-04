##########################################################################
##### PMMST: Phenotype Data Processing Cardio-metabolic risk markers #####
##########################################################################

# 04/07/2017

# load packages
library(haven)
library(dplyr)
library(readxl)
library(data.table)
library(ggplot2)
library(broom)

#### BLOOD PRESSURE ####

# Data were sent by Ayden on 31/05/2017. 
# File wasin a zipped folder. Folder name is Emphasis_GMB_CV_data.zip. File name is
# Emphasis_CV1_AnthropsBP_170531 (excel file)
# On June 7, 2017 queries were sent to Ayden. Ayden came back with answers to most of those queries on
# 30/60/2017. A new version of the data was sent (Emphasis_CV1_AnthropsBP_170630.xlsx)


# load datasets with anthropometry
datBP <- read_excel("./PMMST Outcome Data/Emphasis_CV1_AnthropsBP_170630.xlsx", 
                    sheet = "AnthropBP", na = "NULL", col_names = TRUE)

# sanity check on imported dataset
dim(datBP)
str(datBP)
# keep only variable needed:
datBP <- datBP %>% dplyr::select(cISubjectID, dVisitDate, nSystolicBP1, nSystolicBP2, 
                                 nSystolicBP3, nDiastolicBP1, nDiastolicBP2, nDiastolicBP3)
# Transform every ID in upper case
datBP$cISubjectID <- toupper(datBP$cISubjectID)

# add information on sex, age, allocation group
# add sex and age into the dataset
datsex_age <- read_excel("./PMMST Outcome Data/PMMST_pdata_redux_170718.xlsx")
# a new dataset containing an updated version of sex and age was sent by Ayden on July 18, 2017.
# The new dataset has 2 sex mismatches that were corrected after methylation analysis.
# form datage select only the variables needed for the analysis
datsex_age <- datsex_age %>% dplyr::select(cSubjectID, Age, Sex, MasterGroupNo)
# rename cSubjectID to match the name in dat
names(datsex_age)[names(datsex_age) == "cSubjectID"] <- "cISubjectID"
names(datsex_age)[names(datsex_age) == "Age"] <- "age"
names(datsex_age)[names(datsex_age) == "Sex"] <- "sex"
datBP$cISubjectID <- as.character(datBP$cISubjectID)
datsex_age$age <- as.numeric(datsex_age$age)
# merge
datBP <- inner_join(datBP, datsex_age)
rm(datsex_age)

#########################################################################################################
#### Univariate Analysis ################################################################################
#########################################################################################################

# 1. Check each variable singularly
# 2. Check if same measures are consistent
# 3. Compute the mean of each observation
# 4. Check for outliers using the average computed in 3. and the criteria established in the analysis 
#    plan (median +/- 5MAD)
# 5. Retain the outliers but create a variable with 0/1 (named flagvariablename) where 1 indicated
#    the outlier based on the analysis plan.

#### Systolic Blood Pressure ####
with(datBP, summary(nSystolicBP1))
with(datBP, sd(nSystolicBP1, na.rm = T))
# some measures seem too high
datBP[which(datBP$nSystolicBP1 >= 200),]
# set those to missing
datBP[datBP$cISubjectID == "EMP057T", c("nSystolicBP1", "nSystolicBP2", "nSystolicBP3")] <- NA
datBP[datBP$cISubjectID == "EMP291B", c("nSystolicBP1", "nSystolicBP2", "nSystolicBP3")] <- NA

with(datBP, summary(nSystolicBP2))
with(datBP, sd(nSystolicBP2, na.rm = T))
with(datBP, summary(nSystolicBP3))
with(datBP, sd(nSystolicBP3, na.rm = T))
# compare the distribution of SBP 1, 2 and 3
par(mfrow = c(3,1))
with(datBP, hist(nSystolicBP1))
with(datBP, hist(nSystolicBP2))
with(datBP, hist(nSystolicBP3))

# compute average SBP
datBP$cavgsbp <- with(datBP, (nSystolicBP1 + nSystolicBP2 + nSystolicBP3)/3)
summary(datBP$cavgsbp)
par(mfrow = c(1,1))
hist(datBP$cavgsbp)
attributes(datBP$cavgsbp)$label <- "pmmst - average SBP (mmHg)"

# check for outliers
cutoffU <- median(datBP$cavgsbp, na.rm = T) + 5*mad(datBP$cavgsbp, na.rm = T)
cutoffL <- median(datBP$cavgsbp, na.rm = T) - 5*mad(datBP$cavgsbp, na.rm = T)
sum(datBP$cavgsbp > cutoffU, na.rm = T)
sum(datBP$cavgsbp < cutoffL, na.rm = T)
# no outliers

# plot
ggplot(datBP, aes(x = cavgsbp)) + geom_histogram(colour="black", fill="white", binwidth = 2.5) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 50)) + 
  ggtitle("PMMST: Systolic Blood Pressure (mmHg)") +
  xlab("SBP (mmHg)")


#### Diastolic Blood Pressure ####
with(datBP, summary(nDiastolicBP1))
with(datBP, sd(nDiastolicBP1, na.rm = T))
with(datBP, summary(nDiastolicBP2))
with(datBP, sd(nDiastolicBP2, na.rm = T))
with(datBP, summary(nDiastolicBP3))
with(datBP, sd(nDiastolicBP3, na.rm = T))
# compare the distribution of SBP 1, 2 and 3
par(mfrow = c(3,1))
with(datBP, hist(nDiastolicBP1))
with(datBP, hist(nDiastolicBP2))
with(datBP, hist(nDiastolicBP3))

# differences in the 3 measures of BP are possible. Only check unplausible values
datBP[datBP$nDiastolicBP1 > 100,]
datBP[datBP$nDiastolicBP1 < 30,]
datBP[datBP$nDiastolicBP2 > 100,]
datBP[datBP$nDiastolicBP2 < 30,]
datBP[datBP$nDiastolicBP3 > 100,]
datBP[datBP$nDiastolicBP3 < 30,]

# compute average DBP
datBP$cavgdbp <- with(datBP, (nDiastolicBP1 + nDiastolicBP2 + nDiastolicBP3)/3)
summary(datBP$cavgdbp)
par(mfrow = c(1,1))
hist(datBP$cavgdbp)
attributes(datBP$cavgdbp)$label <- "pmmst - average DBP (mmHg)"

# check for outliers
cutoffU <- median(datBP$cavgdbp, na.rm = T) + 5*mad(datBP$cavgdbp, na.rm = T)
cutoffL <- median(datBP$cavgdbp, na.rm = T) - 5*mad(datBP$cavgdbp, na.rm = T)
sum(datBP$cavgdbp > cutoffU, na.rm = T)
sum(datBP$cavgdbp < cutoffL, na.rm = T)

# plot
ggplot(datBP, aes(x = cavgdbp)) + geom_histogram(colour="black", fill="white", binwidth = 2.5) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 50)) + 
  ggtitle("PMMST: Diastolic Blood Pressure (mmHg)") +
  xlab("DBP (mmHg)")


#########################################################################################################
#### Derive Measures ####################################################################################
#########################################################################################################

# 1. Compute derived measures listed in the analysis plan
# 2. Check if same measures are consistent
# 3. Check for outlier using the criteria established in the analysis plan (median +/- 5MAD)
# 4. Retain the outliers but create a variable with 0/1 (named flagvariablename) where 1 indicated
#    the outlier based on the analysis plan.

#### Compute >99th percentile blood pressure ####

# based on sex, height and BP percentile 
# document with instruction on how to compute percentiles was downloaded on 28/06/2017
# (hbp_ped.pdf)

# Refer to the most recent CDC growth charts, which are available online, and convert the
# height to a height Z-score relative to children of the same age; this is denoted by Zht.

# zhfa was already computed import the data and take the variable we need
datHFA <- read.csv("PMMSTZscore_z.csv")
datHFA <- datHFA %>% dplyr::select(cISubjectID, zhfa)
datBP <- inner_join(datBP, datHFA)
rm(datHFA)

# Compute the expected BP (μ) for boys/girls of age y years and height zhfa inches
# create the function
meanBP <- function(age, height, alpha, b1, b2, b3, b4, g1, g2, g3, g4){
  alpha + b1*(age - 10) + b2*(age - 10)^2 + b3*(age - 10)^3 + b4*(age - 10)^5 + 
    g1*height + (g2*height)^2 + (g3*height)^3 + (g4*height)^4
}

# Systolic Blood Pressure 
meanSBP <- rep(NA, nrow(datBP))
for(i in 1:nrow(datBP)){
  if(datBP$sex[i] == "M" & !is.na(datBP$zhfa[i])){
    meanSBP[i] <- meanBP(age = datBP$age[i], height = datBP$zhfa[i], 
                         alpha = 102.19768, b1 = 1.82416, b2 = 0.12776, b3 = 0.12776, b4 = -0.00135, 
                         g1 = 2.73157, g2 = -0.19618, g3 = -0.04659, g4 = 0.00947)}
  if(datBP$sex[i] == "F"  & !is.na(datBP$zhfa[i])){
    meanSBP[i] <- meanBP(age = datBP$age[i], height = datBP$zhfa[i], 
                         alpha = 102.01027, b1 = 1.94397, b2 = 0.00598, b3 = -0.000789, b4 = -0.00059, 
                         g1 = 2.03526, g2 = 0.02534, g3 = -0.01884, g4 = 0.00121)}
}
# add to the dataset
datBP$meanSBP <- meanSBP

# Diastolic Blood Pressure 
meanDBP <- rep(NA, nrow(datBP))
for(i in 1:nrow(datBP)){
  if(datBP$sex[i] == "M" & !is.na(datBP$zhfa[i])){
    meanDBP[i] <- meanBP(age = datBP$age[i], height = datBP$zhfa[i], 
                         alpha = 61.01217, b1 = 0.68314, b2 = -0.09835, b3 = 0.01711, b4 = 0.00045, 
                         g1 = 1.46993, g2 = -0.07849, g3 = -0.03144, g4 = 0.00967)}
  if(datBP$sex[i] == "F"  & !is.na(datBP$zhfa[i])){
    meanDBP[i] <- meanBP(age = datBP$age[i], height = datBP$zhfa[i], 
                         alpha = 60.50510, b1 = 1.01301, b2 = 0.01157, b3 = 0.00424, b4 = -0.00137, 
                         g1 = 1.16641, g2 = 0.12795, g3 = -0.03869, g4 = -0.00079)}
}
# add to the dataset
datBP$meaDBP <- meanDBP


# Convert the observed BP to a Z-score (Zbp)
# add a sigma in the dataset based on the child sex
datBP$sigmaSBP <- NA
datBP$sigmaDBP <- NA
for(i in 1:nrow(datBP)){
  if(datBP$sex[i] == "M"){
    datBP$sigmaSBP[i] <- 10.7128
    datBP$sigmaDBP[i] <- 11.6032}
  if(datBP$sex[i] == "F"){
    datBP$sigmaSBP[i] <- 10.4855
    datBP$sigmaDBP[i] <- 10.9573
  }
}
# create z-scores
datBP$zsys <- with(datBP, (cavgsbp - meanSBP)/sigmaSBP)
datBP$zdia <- with(datBP, (cavgdbp - meanDBP)/sigmaDBP)

# convert z-scores to percentile
datBP$psys <- pnorm(datBP$zsys)*100
datBP$pdia <- pnorm(datBP$zdia)*100

# check those whose percentile is greater than 99
datBP$cSBPper99 <- ifelse(datBP$psys >= 99, 1, 0)
table(datBP$cSBPper99)
datBP$cSBPper99 <- factor(datBP$cSBPper99, labels = c("No", "Yes"))

datBP$cDBPper99 <- ifelse(datBP$pdia >= 99, 1, 0)
table(datBP$cDBPper99)
datBP$cDBPper99 <- factor(datBP$cDBPper99, labels = c("No", "Yes"))


#### BASELINE AND 30 MINUTES INSULIN ####

# Data were sent by Ebrima Comma 24/08/2017. 
# File name is BiochemistryResults-2017-08-22 (excel file)

# load dataset
datIns <- read_excel("./PMMST Outcome Data/BiochemistryResults-2017-08-22.xlsx", 
                     na = "NULL", col_names = TRUE, sheet = "InsulinNoTime")%>% 
  dplyr::select(cISubjectID, Description, cResult)
table(datIns$Description)
# take baseline values
baseline <- datIns %>% dplyr::filter(Description == "BaseLine") %>% select(cISubjectID, cResult)
# missing subjects
datBP[(datBP$cISubjectID %in% baseline$cISubjectID) == FALSE,]
names(baseline) <- c("cISubjectID", "cInsBaseline") 
# merge with BP dat
dat <- left_join(datBP, baseline)
rm(datBP, baseline)
# do similar thing with 30 minutes insulin
min30ins <- datIns %>% dplyr::filter(Description == "30 Minutes") %>% 
  select(cISubjectID, cResult)
# missing subjects
dat[(dat$cISubjectID %in% min30ins$cISubjectID) == FALSE,]
names(min30ins) <- c("cISubjectID", "cIns30min") 
# merge with BP dat
dat <- left_join(dat, min30ins)
rm(min30ins, datIns)

# insulin had to be in pmol/l (need to multiply by 6).
dat$cInsBaseline <- 6*dat$cInsBaseline
dat$cIns30min <- 6*dat$cIns30min

# check for outliers in fasting insulin
summary(dat$cInsBaseline)
cutoffU <- median(dat$cInsBaseline, na.rm = T) + 5*mad(dat$cInsBaseline, na.rm = T)
cutoffL <- median(dat$cInsBaseline, na.rm = T) - 5*mad(dat$cInsBaseline, na.rm = T)
sum(dat$cInsBaseline > cutoffU, na.rm = T)
sum(dat$cInsBaseline < cutoffL, na.rm = T)
hist(dat$cInsBaseline)
# four outliers fating insulin create flag variable
dat$flagcInsBaseline <- factor(ifelse(dat$cInsBaseline > cutoffU, "Outlier", "No Outlier"))
table(dat$flagcInsBaseline)
rm(cutoffU, cutoffL)

# check for outliers in 30 minutes insulin
summary(dat$cIns30min)
cutoffU <- median(dat$cIns30min, na.rm = T) + 5*mad(dat$cIns30min, na.rm = T)
cutoffL <- median(dat$cIns30min, na.rm = T) - 5*mad(dat$cIns30min, na.rm = T)
sum(dat$cIns30min > cutoffU, na.rm = T)
sum(dat$cIns30min < cutoffL, na.rm = T)
hist(dat$cIns30min)
# four outliers fating insulin create flag variable
dat$flagcIns30min <- factor(ifelse(dat$cIns30min > cutoffU, "Outlier", "No Outlier"))
table(dat$flagcIns30min)
rm(cutoffU, cutoffL)



#### OTHER CARDIO-METABOLIC MARKERS ####

# Data were sent by Modupeth on July 3, 2018
# A new version of the data was sent (EMPHASIS_BioChemFinal_Musa_July2018.xlsx)

# load datasets with cardio-metabolic markers
datCardio <- read_excel("./PMMST Outcome Data/EMPHASIS_BioChemFinal_Musa_July2018.xlsx",
                        col_names = TRUE)
# keep only relevant variables
datCardio <- datCardio %>% 
  dplyr::select(SubjectID, TimpePoint, Test, Value, Unit, Flag, Comment)
# subject EMP291P should be EMP291B
# datCardio[datCardio$SubjectID == "EMP291P","SubjectID"] <- "EMP291B"
# sort by subject ID
datCardio <- datCardio[order(datCardio$SubjectID),] 

# look at available measures
table(datCardio$Test) 
# The following are available:
# total cholesterol
# glucose
# HDL cholesterol
# LDL cholesterol
# Triglycerides
# look at available times
table(datCardio$TimpePoint) 


#### Total cholesterol ####

datChol <- datCardio %>% filter(Test == "CHOL2")
# there are 277 measure of cholesterol
table(datChol$SubjectID)
length(unique(datChol$SubjectID))
# missing total cholesterol
datCardio$SubjectID[!(datCardio$SubjectID %in% datChol$SubjectID)]
# tranform value into numeric
datChol$Value <- as.numeric(datChol$Value)
summary(datChol$Value)
hist(datChol$Value)

# check for outliers
cutoffU <- median(datChol$Value, na.rm = T) + 5*mad(datChol$Value, na.rm = T)
cutoffL <- median(datChol$Value, na.rm = T) - 5*mad(datChol$Value, na.rm = T)
sum(datChol$Value > cutoffU, na.rm = T)
sum(datChol$Value < cutoffL, na.rm = T)
# 1 outlier based on the criterium established in the analysis plan.
# Flag the outlier and create summary statististics with and without the outlier
datChol$flagChol <- ifelse(datChol$Value < cutoffL | datChol$Value > cutoffU, 1, 0)
datChol$flagChol <- factor(datChol$flagChol, labels = c("No Outlier", "Outlier"))
table(datChol$flagChol)
rm(cutoffU); rm(cutoffL)
attributes(datChol$flagChol)$label <- "pmmst - total cholesterol outlier?"
# summary and plot (no outlier removed)
summary(datChol$Value)
ggplot(datChol, aes(x = Value)) + geom_histogram(colour="black", fill="white", binwidth = 0.5) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 95)) + ggtitle("PMMST: Total Cholesterol (mmol/L)") +
  xlab("Total Cholesterol (mmol/L)")
# excluding the outliers
datChol %>% filter(flagChol == "No Outlier") %>% 
  summarise(n = n(), Min = min(Value), Mean = mean(Value), SD = 
              sd(Value), Median = median(Value), Q1 = 
              quantile(Value, probs = 0.25), Q3 = 
              quantile(Value, probs = 0.75), Max = max(Value))
ggplot(subset(datChol, flagChol == "No Outlier"), aes(x = Value)) + 
  geom_histogram(colour="black", fill="white", binwidth = 0.5) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 95)) + 
  ggtitle("PMMST (No Outliers): Total Cholesterol (mmol/L)") +
  xlab("Total Cholesterol (mmol/L)")
# select only relevant variables
datChol <- datChol %>% dplyr::select(SubjectID, Value, flagChol)
names(datChol) <- c("cISubjectID", "cTotChol", "flagcTotChol")
attributes(datChol$cTotChol)$label <- "pmmst - total cholesterol (mmol/L)"
# merge with the full data 
dat <- full_join(dat, datChol)
# check those without measures of total cholesterol
temp <- dat[is.na(dat$cTotChol),]
rm(temp, datChol)


#### HDL Cholesterol ####
datChol <- datCardio %>% filter(Test == "HDLC3")
# there are 277 measure of HDL cholesterol
length(unique(datChol$SubjectID))
# missing HDL cholesterol
dat$cISubjectID[!(dat$cISubjectID %in% datChol$SubjectID)]
# tranform value into numeric
datChol$Value <- as.numeric(datChol$Value)
summary(datChol$Value)
hist(datChol$Value)
# check for outliers
cutoffU <- median(datChol$Value, na.rm = T) + 5*mad(datChol$Value, na.rm = T)
cutoffL <- median(datChol$Value, na.rm = T) - 5*mad(datChol$Value, na.rm = T)
sum(datChol$Value > cutoffU, na.rm = T)
sum(datChol$Value < cutoffL, na.rm = T)
# no outliers for HDL cholesterol
ggplot(datChol, aes(x = Value)) + geom_histogram(colour="black", fill="white", binwidth = 0.2) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 95)) + ggtitle("PMMST: HDL Cholesterol (mmol/L)") +
  xlab("HDL Cholesterol (mmol/L)")
# select only relevant variables
datChol <- datChol %>% dplyr::select(SubjectID, Value)
names(datChol) <- c("cISubjectID", "cHDLChol")
attributes(datChol$cHDLChol)$label <- "pmmst - hdl cholesterol (mmol/L)"
# merge with the full data 
dat <- left_join(dat, datChol)
rm(datChol)

# look into LDL cholesterol (need to take baseline only)
datChol <- datCardio %>% filter(Test == "LDL_C" & TimpePoint == "B-L")
# there are 268 measure of LDL cholesterol
length(unique(datChol$SubjectID))
# missing LDL cholesterol
dat$cISubjectID[!(dat$cISubjectID %in% datChol$SubjectID)]
# tranform value into numeric
datChol$Value <- as.numeric(datChol$Value)
summary(datChol$Value)
hist(datChol$Value)
# check for outliers
cutoffU <- median(datChol$Value, na.rm = T) + 5*mad(datChol$Value, na.rm = T)
cutoffL <- median(datChol$Value, na.rm = T) - 5*mad(datChol$Value, na.rm = T)
sum(datChol$Value > cutoffU, na.rm = T)
sum(datChol$Value < cutoffL, na.rm = T)
# one outlier identified from the analysis plan
# Flag the outlier and create summary statististics with and without the outlier
datChol$flagChol <- ifelse(datChol$Value < cutoffL | datChol$Value > cutoffU, 1, 0)
datChol$flagChol <- factor(datChol$flagChol, labels = c("No Outlier", "Outlier"))
table(datChol$flagChol)
rm(cutoffU); rm(cutoffL)
attributes(datChol$flagChol)$label <- "pmmst - LDL cholesterol outlier?"
# summary and plot (no outlier removed)
summary(datChol$Value)
ggplot(datChol, aes(x = Value)) + geom_histogram(colour="black", fill="white", binwidth = 0.25) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 65)) + 
  ggtitle("PMMST: LDL Cholesterol (mmol/L)") +
  xlab("LDL Cholesterol (mmol/L)")
# excluding the outliers
datChol %>% filter(flagChol == "No Outlier") %>% 
  summarise(n = n(), Min = min(Value), Mean = mean(Value), SD = 
              sd(Value), Median = median(Value), Q1 = 
              quantile(Value, probs = 0.25), Q3 = 
              quantile(Value, probs = 0.75), Max = max(Value))
ggplot(subset(datChol, flagChol == "No Outlier"), aes(x = Value)) + 
  geom_histogram(colour="black", fill="white", binwidth = 0.25) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,65)) + 
  ggtitle("PMMST (No Outliers): LDL Cholesterol (mmol/L)") +
  xlab("LDL Cholesterol (mmol/L)")
# select only relevant variables
datChol <- datChol %>% dplyr::select(SubjectID, Value, flagChol)
names(datChol) <- c("cISubjectID", "cLDLChol", "flagcLDLChol")
attributes(datChol$cLDLChol)$label <- "pmmst - LDL cholesterol (mmol/L)"
# merge with the full data 
dat <- left_join(dat, datChol)



#### Triglycerides ####

datTri <- datCardio %>% filter(Test == "TRIGL")
# there are 277 measure of triglycerides cholesterol
table(datTri$SubjectID)
length(unique(datTri$SubjectID))
# tranform value into numeric
datTri$Value <- as.numeric(datTri$Value)
summary(datTri$Value)
hist(datTri$Value)
# check for outliers
cutoffU <- median(datTri$Value, na.rm = T) + 5*mad(datTri$Value, na.rm = T)
cutoffL <- median(datTri$Value, na.rm = T) - 5*mad(datTri$Value, na.rm = T)
sum(datTri$Value > cutoffU, na.rm = T)
sum(datTri$Value < cutoffL, na.rm = T)
# no outlier observed
summary(datTri$Value)
ggplot(datTri, aes(x = Value)) + geom_histogram(colour="black", fill="white", binwidth = 0.1) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 65)) + ggtitle("PMMST: Triglycerides (mmol/L)") +
  xlab("Triglycerides (mmol/L)")
# select only relevant variables
datTri <- datTri %>% dplyr::select(SubjectID, Value)
names(datTri) <- c("cISubjectID", "cTrigl")
attributes(datTri$cTrigl)$label <- "pmmst - triglycerides (mmol/L)"
# merge with the full data 
dat <- left_join(dat, datTri)
# check those without measures of total cholesterol
temp <- dat[is.na(dat$cTrigl),]
rm(temp, datTri)


#### Fasting Glucose ####

datGlu <- datCardio %>% filter(Test == "GLU2")
# look at baseline values
datBaseline <- datGlu %>% filter(TimpePoint == "B-L")
# missing subjects
dat[(dat$cISubjectID %in% datBaseline$SubjectID) == FALSE,]
datBaseline <- datBaseline %>% select(SubjectID, Value) 
names(datBaseline) <- c("cISubjectID", "cGluBaseline") 
# merge with the full data 
dat <- left_join(dat, datBaseline)
rm(datBaseline)
dat$cGluBaseline <- as.numeric(dat$cGluBaseline)
summary(dat$cGluBaseline)
hist(dat$cGluBaseline)

# check for outliers
cutoffU <- median(dat$cGluBaseline, na.rm = T) + 5*mad(dat$cGluBaseline, na.rm = T)
cutoffL <- median(dat$cGluBaseline, na.rm = T) - 5*mad(dat$cGluBaseline, na.rm = T)
sum(dat$cGluBaseline > cutoffU, na.rm = T)
sum(dat$cGluBaseline < cutoffL, na.rm = T)
# no outliers for fasting glucose
ggplot(dat, aes(x = cGluBaseline)) + 
  geom_histogram(colour="black", fill="white", binwidth = 0.2) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 65)) + 
  ggtitle("PMMST: Fasting Glucose (mmol/L)") +
  xlab("Fasting Glucose (mmol/L)")


#### 30min Glucose ####

dat30 <- datGlu %>% filter(TimpePoint == "30")
# missing subjects
dat[(dat$cISubjectID %in% dat30$SubjectID) == FALSE,]
dat30 <- dat30 %>% select(SubjectID, Value) 
names(dat30) <- c("cISubjectID", "cGlu30min") 
# merge with the full data 
dat <- left_join(dat, dat30)
rm(dat30)
dat$cGlu30min <- as.numeric(dat$cGlu30min)
summary(dat$cGlu30min)
hist(dat$cGlu30min)
dat[dat$cGlu30min < 2 & !is.na(dat$cGlu30min), ]
temp <- dat[is.na(dat$cGlu30min), ]
# 30 minutes glucose of EMP210C is too low (0.34). Set as missing
dat[dat$cGlu30min < 2 & !is.na(dat$cGlu30min), "cGlu30min"] <- NA
# check for outliers
cutoffU <- median(dat$cGlu30min, na.rm = T) + 5*mad(dat$cGlu30min, na.rm = T)
cutoffL <- median(dat$cGlu30min, na.rm = T) - 5*mad(dat$cGlu30min, na.rm = T)
sum(dat$cGlu30min > cutoffU, na.rm = T)
sum(dat$cGlu30min < cutoffL, na.rm = T)
# no outliers for 30 minutes glucose
ggplot(dat, aes(x = cGlu30min)) + 
  geom_histogram(colour="black", fill="white", binwidth = 0.5) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 65)) + 
  ggtitle("PMMST: 30min Glucose (mmol/L)") +
  xlab("30min Glucose (mmol/L)")


#### 120min Glucose ####
dat120 <- datGlu %>% filter(TimpePoint == "120")
# missing subjects
dat[(dat$cISubjectID %in% dat120$SubjectID) == FALSE,]
dat120 <- dat120 %>% select(SubjectID, Value) 
names(dat120) <- c("cISubjectID", "cGlu120min") 
# merge with the full data 
dat <- left_join(dat, dat120)
rm(dat120)
dat$cGlu120min <- as.numeric(dat$cGlu120min)
summary(dat$cGlu120min)

dat[dat$cGlu120min < 2 & !is.na(dat$cGlu120min), ]
# check for outliers
cutoffU <- median(dat$cGlu120min, na.rm = T) + 5*mad(dat$cGlu120min, na.rm = T)
cutoffL <- median(dat$cGlu120min, na.rm = T) - 5*mad(dat$cGlu120min, na.rm = T)
sum(dat$cGlu120min > cutoffU, na.rm = T)
sum(dat$cGlu120min < cutoffL, na.rm = T)
# 1 outlier for 120 minutes glucose
# create a flag variable for 120 minutes glucose
dat$flagcGlu120min <- ifelse(dat$cGlu120min < cutoffL, 1, 0)
dat$flagcGlu120min <- factor(dat$flagcGlu120min, labels = c("No Outlier", "Outlier"))
table(dat$flagcGlu120min)
rm(cutoffU); rm(cutoffL)
attributes(dat$flagcGlu120min)$label <- "pmmst - 120min glucose outlier?"
# summary and plot (no outlier removed)
summary(dat$cGlu120min)
ggplot(dat, aes(x = cGlu120min)) + 
  geom_histogram(colour="black", fill="white", binwidth = 0.5) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 65)) + 
  ggtitle("PMMST: 120min Glucose (mmol/L)") +
  xlab("120min Glucose (mmol/L)")
# excluding the outliers
dat %>% filter(flagcGlu120min == "No Outlier") %>% 
  summarise(n = n(), Min = min(cGlu120min), Mean = mean(cGlu120min), SD = 
              sd(cGlu120min), Median = median(cGlu120min), Q1 = 
              quantile(cGlu120min, probs = 0.25), Q3 = 
              quantile(cGlu120min, probs = 0.75), Max = max(cGlu120min))
ggplot(subset(dat, flagcGlu120min == "No Outlier"), aes(x = cGlu120min)) + 
  geom_histogram(colour="black", fill="white", binwidth = 0.5) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 65)) + 
  ggtitle("PMMST (No Outliers): 120min Glucose (mmol/L)") +
  xlab("120min Glucose (mmol/L)")


# look at the relation between 0, 30 minutes and 120 minutes glucose.
temp <- dat %>% select(cISubjectID, cGluBaseline, cGlu30min, cGlu120min)
# glucose should go up between baseline and 30 minutes
temp$diff <- temp$cGlu30min - temp$cGluBaseline
# and go down again between 30 minutes and 120 minutes
temp$diff2 <- temp$cGlu120min - temp$cGlu30min
rm(temp)

rm(datCardio, datChol, datGlu)

#########################################################################################################
#### Derive Measures ####################################################################################
#########################################################################################################

# 1. Compute derived measures listed in the analysis plan
# 2. Check if same measures are consistent
# 3. Check for outlier using the criteria established in the analysis plan (median +/- 5MAD)
# 4. Retain the outliers but create a variable with 0/1 (named flagvariablename) where 1 indicated
#    the outlier based on the analysis plan.


#### Compute HOMA-IR ####

# save serno, fasting insulin and glucose in a separate file to be imported in Excel for HOMA computation. 

# to run only the first time
# fileForHOMA <- dat %>% dplyr::select(cISubjectID, cGluBaseline, cInsBaseline)
# fasting insulin needs to be at least 20 for the calculator to work. Transform all the values <20 to 20
# fileForHOMA$cInsBaseline[fileForHOMA$cInsBaseline < 20] <- 20
# head(fileForHOMA)
# write.table(fileForHOMA, "PMMSTfileForHOMA.csv", sep = ",", row.names = F)
# rm(fileForHOMA)
# an excel spreadsheet named childHOMA.xlxs was created with the HOMA values
# load the HOMA values and add them to the dataset

# start again from here
datHOMA <-  read_excel("ChildHOMA.xlsx", sheet = "Gambia", na = "NA")
datHOMA <- datHOMA %>% dplyr::select(cISubjectID, chomas, chomair)
# merge datHOMA with datCardio
dat <- inner_join(dat, datHOMA)
rm(datHOMA)
dat$chomair <- as.numeric(dat$chomair)
dat$chomas <- as.numeric(dat$chomas)

# summary of HOMA-IR
summary(dat$chomair)
sd(dat$chomair, na.rm = T)
hist(dat$chomair)

# check for outliers
cutoffU <- median(dat$chomair, na.rm = T) + 5*mad(dat$chomair, na.rm = T)
cutoffL <- median(dat$chomair, na.rm = T) - 5*mad(dat$chomair, na.rm = T)
sum(dat$chomair > cutoffU, na.rm = T)
sum(dat$chomair < cutoffL, na.rm = T)
# based on the analysis plan there are 53 outliers. All with highe HOMA-IR values

# create a flag variable for triglycerides
dat$flagchomair <- ifelse(dat$chomair > cutoffU, 1, 0)
dat$flagchomair <- factor(dat$flagchomair, labels = c("No Outlier", "Outlier"))
table(dat$flagchomair)
rm(cutoffU); rm(cutoffL)
attributes(dat$flagchomair)$label <- "pmmst - HOMA-IR outlier?"

# summary and plot (no outlier removed)
summary(dat$chomair)
ggplot(dat, aes(x = chomair)) + geom_histogram(colour="black", fill="white", binwidth = 0.08) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 200)) + ggtitle("PMMST: HOMA-IR") +
  xlab("HOMA-IR")
# excluding the outliers
dat %>% filter(flagchomair == "No Outlier") %>% 
  summarise(n = n(), Min = min(chomair), Mean = mean(chomair), SD = 
              sd(chomair), Median = median(chomair), Q1 = 
              quantile(chomair, probs = 0.25), Q3 = 
              quantile(chomair, probs = 0.75), Max = max(chomair))
ggplot(subset(dat, flagchomair == "No Outlier"), aes(x = chomair)) + 
  geom_histogram(colour="black", fill="white", binwidth = 0.01) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) + 
  ggtitle("PMMST (No Outliers): HOMA-IR") +
  xlab("HOMA-IR")

# HOMA-IR is the reciprocal of HOMA-S (100/%S). Check:
plot(dat$chomair, 1/dat$chomas)
cor(dat$chomair, 1/dat$chomas, use = "complete.obs")


#### Compute Insulinogenic Index ####
# computed as (Insulin30 - Insulin0)/(Glucose30-Glucose0)
# In the Mumbai data there were quite a lot of children whose 30 minute glucose 
# is lower than their fasting value. Clive suggests using the following to calculate the insulogenic index, to deal with this problem:
# (log [insulin30] – log [insulin0]) / (log [glucose30] – log [glucose0])
# Here we use the same formula for consistency

dat$cinsindex <- with(dat, (log(cIns30min) - log(cInsBaseline))/(log(cGlu30min) - log(cGluBaseline)))
summary(dat$cinsindex)
# 2 children (EMP172L and EMP262V) have infinite values for Insulinogenic Index. 
# This happened because their fasting
# insulin is 0. Set those two children as missing (might need to change once their fasting
# insulin in checked with the data manager.
dat[dat$cISubjectID %in% c("EMP172L", "EMP262V"), "cinsindex"] <- NA
sum(!is.na(dat$cinsindex)) 
# 262 have insulinogenic index
# how many have insulinogenic index negative?
dat[dat$cinsindex <0 & !is.na(dat$cinsindex), c("cISubjectID", "cinsindex")]
# 15 out of 262
(15/262)*100
hist(dat$cinsindex)

# check for outliers
cutoffU <- median(dat$cinsindex, na.rm = T) + 5*mad(dat$cinsindex, na.rm = T)
cutoffL <- median(dat$cinsindex, na.rm = T) - 5*mad(dat$cinsindex, na.rm = T)
sum(dat$cinsindex > cutoffU, na.rm = T)
sum(dat$cinsindex < cutoffL, na.rm = T)
# based on the analysis plan there are 28 outliers

# create a flag variable for insulinogenic index
dat$flagcinsindex <- ifelse(dat$cinsindex > cutoffU | dat$cinsindex < cutoffL, 1, 0)
dat$flagcinsindex <- factor(dat$flagcinsindex, labels = c("No Outlier", "Outlier"))
table(dat$flagcinsindex)
rm(cutoffU); rm(cutoffL)
attributes(dat$flagcinsindex)$label <- "pmmst - Insulinogenic index outlier?"

# summary and plot (no outlier removed)
summary(dat$cinsindex)
ggplot(dat, aes(x = cinsindex)) + 
  geom_histogram(colour="black", fill="white", binwidth = 30) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 300)) + 
  ggtitle("PMMST: Insulinogenic Index") +
  xlab("Insulinogenic Index")

# excluding the outliers
dat %>% filter(flagcinsindex == "No Outlier") %>% 
  summarise(n = n(), Min = min(cinsindex), Mean = mean(cinsindex), SD = 
              sd(cinsindex), Median = median(cinsindex), Q1 = 
              quantile(cinsindex, probs = 0.25), Q3 = 
              quantile(cinsindex, probs = 0.75), Max = max(cinsindex))
ggplot(subset(dat, flagcinsindex == "No Outlier"), aes(x = cinsindex)) + 
  geom_histogram(colour="black", fill="white", binwidth = 1) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 75)) + 
  ggtitle("PMMST (No Outliers): Insulinogenic Index") +
  xlab("Insulinogenic Index")


#### Compute Disposition Index ####

# disposition index is calculated as Insulogenic index/HOMA-IR
dat$cdispindex <- with(dat, cinsindex/chomair)
summary(dat$cdispindex)
sd(dat$cdispindex, na.rm = T)
sum(!is.na(dat$cdispindex)) # 262

# check for outliers
cutoffU <- median(dat$cdispindex, na.rm = T) + 5*mad(dat$cdispindex, na.rm = T)
cutoffL <- median(dat$cdispindex, na.rm = T) - 5*mad(dat$cdispindex, na.rm = T)
sum(dat$cdispindex > cutoffU, na.rm = T)
sum(dat$cdispindex < cutoffL, na.rm = T)
# based on the analysis plan there are 23 outliers

# create a flag variable for triglycerides
dat$flagcdispindex <- ifelse(dat$cdispindex > cutoffU | dat$cdispindex < cutoffL, 1, 0)
dat$flagcdispindex <- factor(dat$flagcdispindex, labels = c("No Outlier", "Outlier"))
table(dat$flagcdispindex)
rm(cutoffU); rm(cutoffL)
attributes(dat$flagcdispindex)$label <- "pmmst - Disposition index outlier?"

# summary and plot (no outlier removed)
summary(dat$cdispindex)
ggplot(dat, aes(x = cdispindex)) + 
  geom_histogram(colour="black", fill="white", binwidth = 60) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 300)) + 
  ggtitle("PMMST: Disposition Index") +
  xlab("Disposition Index")

# excluding the outliers
dat %>% filter(flagcdispindex == "No Outlier") %>% 
  summarise(n = n(), Min = min(cdispindex), Mean = mean(cdispindex), SD = 
              sd(cdispindex), Median = median(cdispindex), Q1 = 
              quantile(cdispindex, probs = 0.25), Q3 = 
              quantile(cdispindex, probs = 0.75), Max = max(cdispindex))
ggplot(subset(dat, flagcdispindex == "No Outlier"), aes(x = cdispindex)) + 
  geom_histogram(colour="black", fill="white", binwidth = 5) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) + 
  ggtitle("PMMST (No Outliers): Disposition Index") +
  xlab("Disposition Index")

# check insulinogenic index X HOMA-%S
plot(dat$cdispindex, dat$cinsindex*dat$chomas)
# same!

#### Compute Metabolic Syndrome ####

# Create a 0/1 variable where 1 is for children who are above the sex-specific within-population upper quartile 
# for android fat on DXA and SBP and triglycerides and HOMA-IR, and below the lower quartile for HDL-cholesterol.

# add android fat from DXA
datDXA <- read_dta("./Cleaned Dataset/PMMSTGrowthBodyCompv1_20180131.dta")
datDXA <- datDXA %>% filter(cvisit == 2) %>% dplyr::select(serno, candroidfat)
names(datDXA) <- c("cISubjectID", "candroidfat")
# add android fat to dataset
dat <- left_join(dat, datDXA)
rm(datDXA)


# sex-specific within-population higher than the upper quartile for android fat on DXA
DXAquart <- dat %>% group_by(sex) %>% summarise(Q3 = quantile(candroidfat, prob = 0.75, na.rm = T))
qDXA <- c()
for(i in 1:nrow(dat)){
  if(dat$sex[i] == "M" & !is.na(dat$sex[i]) & !is.na(dat$candroidfat[i]) & dat$candroidfat[i] > DXAquart[1,"Q3"]){qDXA[i] <- 1
  } else if(dat$sex[i] == "M" & !is.na(dat$sex[i]) & !is.na(dat$candroidfat[i]) & dat$candroidfat[i] <= DXAquart[1,"Q3"]){
    qDXA[i] <- 0
  } else if(dat$sex[i] == "F" & !is.na(dat$sex[i]) & !is.na(dat$candroidfat[i]) & dat$candroidfat[i] > DXAquart[2,"Q3"]){
    qDXA[i] <- 1
  } else if(dat$sex[i] == "F" & !is.na(dat$sex[i]) & !is.na(dat$candroidfat[i]) & dat$candroidfat[i] <= DXAquart[2,"Q3"]){
    qDXA[i] <- 0
  } else{ qDXA[i] <- NA}
}
dat$qDXA <- qDXA
rm(qDXA)
rm(DXAquart)

# sex-specific within-population higher than the upper quartile for SBP
SBPquart <- dat %>% group_by(sex) %>% summarise(Q3 = quantile(cavgsbp, prob = 0.75, na.rm = T))
qSBP <- c()
for(i in 1:nrow(dat)){
  if(dat$sex[i] == "M" & !is.na(dat$sex[i]) & !is.na(dat$cavgsbp[i]) & dat$cavgsbp[i] > SBPquart[1,"Q3"]){qSBP[i] <- 1
  } else if(dat$sex[i] == "M" & !is.na(dat$sex[i]) & !is.na(dat$cavgsbp[i]) & dat$cavgsbp[i] <= SBPquart[1,"Q3"]){
    qSBP[i] <- 0
  } else if(dat$sex[i] == "F" & !is.na(dat$sex[i]) & !is.na(dat$cavgsbp[i]) & dat$cavgsbp[i] > SBPquart[2,"Q3"]){
    qSBP[i] <- 1
  } else if(dat$sex[i] == "F" & !is.na(dat$sex[i]) & !is.na(dat$cavgsbp[i]) & dat$cavgsbp[i] <= SBPquart[2,"Q3"]){
    qSBP[i] <- 0
  } else{ qSBP[i] <- NA}
}
dat$qSBP <- qSBP
rm(qSBP)
rm(SBPquart)

# sex-specific within-population upper quartile for triglycerides
tgquart <- dat %>% group_by(sex) %>% summarise(Q3 = quantile(cTrigl, prob = 0.75, na.rm = T))
qtg <- c()
for(i in 1:nrow(dat)){
  if(dat$sex[i] == "M" & !is.na(dat$sex[i])  & !is.na(dat$cTrigl[i]) & dat$cTrigl[i] > tgquart[1,"Q3"]){qtg[i] <- 1
  } else if(dat$sex[i] == "M" & !is.na(dat$sex[i])  & !is.na(dat$cTrigl[i]) & dat$cTrigl[i] <= tgquart[1,"Q3"]){
    qtg[i] <- 0
  } else if(dat$sex[i] == "F" & !is.na(dat$sex[i])  & !is.na(dat$cTrigl[i]) & dat$cTrigl[i] > tgquart[2,"Q3"]){
    qtg[i] <- 1
  } else if(dat$sex[i] == "F" & !is.na(dat$sex[i])  & !is.na(dat$cTrigl[i]) & dat$cTrigl[i] <= tgquart[2,"Q3"]){
    qtg[i] <- 0
  } else{ qtg[i] <- NA}
}
dat$qtg <- qtg
rm(qtg)
rm(tgquart)

# sex-specific within-population higher than the upper quartile for HOMA-IR
HOMAquart <- dat %>% group_by(sex) %>% summarise(Q3 = quantile(chomair, prob = 0.75, na.rm = T))
qHOMA <- c()
for(i in 1:nrow(dat)){
  if(dat$sex[i] == "M" & !is.na(dat$sex[i])  & !is.na(dat$chomair[i]) & dat$chomair[i] > HOMAquart[1,"Q3"]){
    qHOMA[i] <- 1
  } else if(dat$sex[i] == "M" & !is.na(dat$sex[i]) & !is.na(dat$chomair[i]) & dat$chomair[i] <= HOMAquart[1,"Q3"]){
    qHOMA[i] <- 0
  } else if(dat$sex[i] == "F" & !is.na(dat$sex[i]) & !is.na(dat$sex[i])  & !is.na(dat$chomair[i]) & dat$chomair[i] > HOMAquart[2,"Q3"]){
    qHOMA[i] <- 1
  } else if(dat$sex[i] == "F" & !is.na(dat$sex[i]) & !is.na(dat$chomair[i]) & dat$chomair[i] <= HOMAquart[2,"Q3"]){
    qHOMA[i] <- 0
  } else{qHOMA[i] <- NA}
}
dat$qHOMA <- qHOMA
rm(qHOMA)
rm(HOMAquart)

# sex-specific within-population below the lower quartile for HDL
HDLquart <- dat %>% group_by(sex) %>% summarise(Q1 = quantile(cHDLChol, prob = 0.25, na.rm = T))
qHDL <- c()
for(i in 1:nrow(dat)){
  if(dat$sex[i] == "M" & !is.na(dat$sex[i]) & !is.na(dat$cHDLChol[i]) & dat$cHDLChol[i] < HDLquart[1,"Q1"]){qHDL[i] <- 1
  } else if(dat$sex[i] == "M" & !is.na(dat$sex[i]) & !is.na(dat$cHDLChol[i]) & dat$cHDLChol[i] >= HDLquart[1,"Q1"]){
    qHDL[i] <- 0
  } else if(dat$sex[i] == "F" & !is.na(dat$sex[i]) & !is.na(dat$cHDLChol[i]) & dat$cHDLChol[i] < HDLquart[2,"Q1"]){
    qHDL[i] <- 1
  } else if(dat$sex[i] == "F" & !is.na(dat$sex[i]) & !is.na(dat$cHDLChol[i]) & dat$cHDLChol[i] >= HDLquart[2,"Q1"]){
    qHDL[i] <- 0
  } else{qHDL[i] <- NA}
}
dat$qHDL <- qHDL
rm(qHDL)
rm(HDLquart)

# build metabolic syndrome
dat$sumQ <- with(dat, qDXA + qSBP + qtg + qHOMA + qHDL)
summary(dat$sumQ)
dat$cmetsynd <- ifelse(dat$sumQ == 5, 1, 0) 
table(dat$cmetsynd)
# 2 children have metabolic syndrome

###################################################################################################################
#### SAVE DATASET #################################################################################################
###################################################################################################################

names(dat)

# select only relevant variable (based on the analysis plan)
dat <- dat %>% dplyr::select(cISubjectID, sex, age, MasterGroupNo, cavgsbp, cavgdbp, cSBPper99,
                             cDBPper99, cInsBaseline, cIns30min, cGluBaseline, cGlu30min, 
                             cGlu120min, cTotChol, cHDLChol, cLDLChol, cTrigl, chomair,
                             cinsindex, cdispindex, cmetsynd,
                             flagcInsBaseline, flagcIns30min, flagcGlu120min, 
                             flagcTotChol, flagcLDLChol, flagchomair, flagcinsindex, 
                             flagcdispindex)          


names(dat) <- c("serno", "csex", "cage", "allocation", "cavgSBP", "cavgDBP", 
                "cSBPper99", "cDBPper99", "cins0", "cins30", "cglu0",
                "cglu30", "cglu120", "cchol", "chdl", "cldl", "ctg", 
                "choma", "cinsindex", "cdispindex",
                "cmetsynd", "flagcins0", "flagcins30", "flagcglu120", "flagcchol", "flagcldl",
                "flagchoma", "flagcinsindex", "flagcdispindex")

table(dat$csex)
dat$csex <- factor(dat$csex, labels = c("Girl", "Boy"))

# add attributes
attributes(dat$serno)$label <- "pmmst - subject ID"
attributes(dat$csex)$label <- "pmmst - sex"
attributes(dat$cage)$label <- "pmmst - age (years)"
attributes(dat$allocation)$label <- "pmmst - allocation group"
attributes(dat$cSBPper99)$label <- "pmmst - SBP greater than 99th percentile"
attributes(dat$cDBPper99)$label <- "pmmst - DBP greater than 99th percentile"
attributes(dat$choma)$label <- "pmmst - HOMA-IR"
attributes(dat$cmetsynd)$label <- "pmmst - metabolic syndrome"
attributes(dat$cglu0)$label <- "pmmst - fasting glucose (mmol/l)"
attributes(dat$cglu30)$label <- "pmmst - 30min glucose (mmol/l)"
attributes(dat$cglu120)$label <- "pmmst - 120 glucose (mmol/l)"
attributes(dat$ctg)$label <- "pmmst - triglycerides (mmol/l)"
attributes(dat$chdl)$label <- "pmmst - HDL cholesterol (mmol/l)"
attributes(dat$cldl)$label <- "pmmst - LDL cholesterol (mmol/l)"
attributes(dat$cins0)$label <- "pmmst - fasting insulin (pmol/l)"
attributes(dat$cins30)$label <- "pmmst - 30min insulin (pmol/l)"
attributes(dat$cinsindex)$label <- "pmmst - insulinogenic index"
attributes(dat$cdispindex)$label <- "pmmst - disposition index"
attributes(dat$flagcins0)$label <- "pmmst - fasting insulin outlier?"
attributes(dat$flagcins30)$label <- "pmmst - 30min insulin outlier?"


# save as csv as stata ans as spss file
write.table(dat, "./Cleaned Dataset/PMMSTCardioMetv1_20180704.csv", row.names = FALSE)
write_dta(dat, "./Cleaned Dataset/PMMSTCardioMetv1_20180704.dta")
write_sav(dat, "./Cleaned Dataset/PMMSTCardioMetv1_20180704.sav")
