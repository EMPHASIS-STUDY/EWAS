##############################################################################
##### PMMST: Phenotype Data Processing Clinic Visit 1 Anthropometry  #########
##############################################################################

# 30/06/2017

# load packages
library(readxl)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(haven)

# Data were sent by Ayden on 31/05/2017. 
# File wasin a zipped folder. Folder name is Emphasis_GMB_CV_data.zip. File name is
# Emphasis_CV1_AnthropsBP_170531 (excel file)
# On June 7, 2017 queries were sent to Ayden. Ayden came back with answers to most of those queries on
# 30/60/2017. A new version of the data was sent (Emphasis_CV1_AnthropsBP_170630.xlsx)

# load datasets with antrhopometry and close all connections
dat <- read_excel("./PMMST Outcome Data/Emphasis_CV1_AnthropsBP_170630.xlsx", 
                  sheet = "AnthropBP", na = "NULL", col_names = TRUE)

# sanity check on imported dataset
dim(dat)
str(dat)

# use data dictionary sent by Ayden on 31/05/2017.
# Check all the subject have anthropometry recorded
# as visit 1.
table(dat$cVisNo)
# Check all subject have a visit date
addmargins(table(dat$dVisitDate))
# visit date format: YYYY-MM-DD

# Transform every ID in upper case
dat$cISubjectID <- toupper(dat$cISubjectID)

########################################################################################################################
#### Univariate Analysis ###############################################################################################
########################################################################################################################

# 1. Check each variable singularly
# 2. Check if same measures are consistent
# 3. Compute the mean of each observation
# 4. Check for outliers using the average computed in 3. and the criteria established in the analysis 
#    plan (median +/- 5MAD)
# 5. Retain the outliers but create a variable with 0/1 (named flagvariablename) where 1 indicated
#    the outlier based on the analysis plan.

#### Weight ####

with(dat, summary(nWeight1))
with(dat, sd(nWeight1))
with(dat, summary(nWeight2))
with(dat, sd(nWeight2))
# compare the distribution of weight 1 and 2
par(mfrow = c(2,1))
with(dat, hist(nWeight1))
with(dat, hist(nWeight2))
# generate a variable looking at the difference between the two measures
diff <- dat$nWeight1 - dat$nWeight2
summary(diff)
rm(diff)

# average the weight measures and create an average weight variable
dat$avgWeight <- (dat$nWeight1 +  dat$nWeight2)/2
summary(dat$avgWeight)
sd(dat$avgWeight)
par(mfrow = c(1,1))
hist(dat$avgWeight)

# check for outliers
cutoffU <- median(dat$avgWeight) + 5*mad(dat$avgWeight)
cutoffL <- median(dat$avgWeight) - 5*mad(dat$avgWeight)
sum(dat$avgWeight > cutoffU)
sum(dat$avgWeight < cutoffL)
# no outliers based on the criteria established in the Analysis Plan
rm(cutoffU); rm(cutoffL)

ggplot(dat, aes(x = avgWeight)) + geom_histogram(colour="black", fill="white", binwidth = 1.5) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 60)) + ggtitle("PMMST: Weight (kg)") +
  xlab("Weight (kg)")



#### Height ####

with(dat, summary(nHeight1))
with(dat, sd(nHeight1))
with(dat, summary(nHeight2))
with(dat, sd(nHeight2))
# compare the distribution of Height 1 and 2
par(mfrow = c(2,1))
with(dat, hist(nHeight1))
with(dat, hist(nHeight2))
# generate a variable looking at the difference between the two measures
diff <- dat$nHeight1 - dat$nHeight2
summary(diff)
# check those with difference greater than -2 cm
dat[which(diff < -2), c("cISubjectID", "nHeight1", "nHeight2", "nWeight1", "nWeight2")]
# nHeight2 is too high compared to the child weight. Set to missing
dat[dat$cISubjectID == "EMP232W","nHeight2"] <- NA
rm(diff)

# average the height measures and create an average height variable
dat$avgHeight <- rowMeans(dat[,c("nHeight1","nHeight2")], na.rm = T)
summary(dat$avgHeight)
sd(dat$avgHeight)
par(mfrow = c(1,1))
hist(dat$avgWeight)

# check for outliers
cutoffU <- median(dat$avgHeight) + 5*mad(dat$avgHeight)
cutoffL <- median(dat$avgHeight) - 5*mad(dat$avgHeight)
sum(dat$avgHeight > cutoffU)
sum(dat$avgHeight < cutoffL)
# no outliers based on the criteria established in the Analysis Plan
rm(cutoffU); rm(cutoffL)

ggplot(dat, aes(x = avgHeight)) + geom_histogram(colour="black", fill="white", binwidth = 1.5) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 45)) + ggtitle("PMMST: Height (cm)") +
  xlab("Height (cm)")



#### Sitting Height ####

with(dat, summary(nSittingHeight1))
with(dat, sd(nSittingHeight1))
with(dat, summary(nSittingHeight2))
with(dat, sd(nSittingHeight2))
# compare the distribution of sitting height 1 and 2
par(mfrow = c(2,1))
with(dat, hist(nSittingHeight1))
with(dat, hist(nSittingHeight2))
# generate a variable looking at the difference between the two measures
diff <- dat$nSittingHeight1 - dat$nSittingHeight2
summary(diff)
rm(diff)

# average the sitting height measures and create an average sitting height variable
dat$avgSitHeight <- (dat$nSittingHeight1 +  dat$nSittingHeight2)/2
summary(dat$avgSitHeight)
sd(dat$avgSitHeight)
par(mfrow = c(1,1))
hist(dat$avgSitHeight)

# set sitting height as missing
dat[dat$cISubjectID == "EMP091Z", "avgSitHeight"] <- NA


# check for outliers
cutoffU <- median(dat$avgSitHeight) + 5*mad(dat$avgSitHeight)
cutoffL <- median(dat$avgSitHeight) - 5*mad(dat$avgSitHeight)
sum(dat$avgSitHeight > cutoffU)
sum(dat$avgSitHeight < cutoffL)
# no outliers based on the criteria established in the Analysis Plan
rm(cutoffU); rm(cutoffL)

ggplot(dat, aes(x = avgSitHeight)) + geom_histogram(colour="black", fill="white", binwidth = 1.5) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 60)) + ggtitle("PMMST: Sitting Height (cm)") +
  xlab("Sitting Height (cm)")



#### Waist Circumference ####

with(dat, summary(nWaist1))
with(dat, sd(nWaist1))
with(dat, summary(nWaist2))
with(dat, sd(nWaist2))
# compare the distribution of waist 1 and 2
par(mfrow = c(2,1))
with(dat, hist(nWaist1))
with(dat, hist(nWaist2))
# generate a variable looking at the difference between the two measures
diff <- dat$nWaist1 - dat$nWaist2
summary(diff)
# check those with difference greater than -2 cm
dat[which(diff < -2), c("cISubjectID", "nWaist1", "nWaist2")]
# measures are correct. Checked with Ayden on 30/06/2017
rm(diff)

# average the waist measures and create an average waist variable
dat$avgWaist <- (dat$nWaist1 +  dat$nWaist2)/2
summary(dat$avgWaist)
sd(dat$avgWaist)
par(mfrow = c(1,1))
hist(dat$avgWaist)

# check for outliers
cutoffU <- median(dat$avgWaist) + 5*mad(dat$avgWaist)
cutoffL <- median(dat$avgWaist) - 5*mad(dat$avgWaist)
sum(dat$avgWaist > cutoffU)
sum(dat$avgWaist < cutoffL)
# no outliers based on the criteria established in the Analysis Plan
rm(cutoffU); rm(cutoffL)

ggplot(dat, aes(x = avgWaist)) + geom_histogram(colour ="black", fill="white", binwidth = 1.5) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 65)) + ggtitle("PMMST: Waist circumference (cm)") +
  xlab("Waist circumference (cm)")



#### Head Circumferece ####

with(dat, summary(nHC1))
with(dat, sd(nHC1))
with(dat, summary(nHC2))
with(dat, sd(nHC2))
# compare the distribution of HC 1 and 2
par(mfrow = c(2,1))
with(dat, hist(nHC1))
with(dat, hist(nHC2))
# generate a variable looking at the difference between the two measures
diff <- dat$nHC1 - dat$nHC2
summary(diff)
rm(diff)

# average the HC measures and create an average HC variable
dat$avgHC <- (dat$nHC1 +  dat$nHC2)/2
summary(dat$avgHC)
sd(dat$avgHC)
par(mfrow = c(1,1))
hist(dat$avgHC)
attributes(dat$avgHC)$label <- "pmmst - average HC (cm)"

# check for outliers
cutoffU <- median(dat$avgHC) + 5*mad(dat$avgHC)
cutoffL <- median(dat$avgHC) - 5*mad(dat$avgHC)
sum(dat$avgHC > cutoffU)
sum(dat$avgHC < cutoffL)
# 1 outlier based on the criteria established in the Analysis Plan

# create a flag variable for HC
dat$flagHC <- ifelse(dat$avgHC < cutoffL, 1, 0)
dat$flagHC <- factor(dat$flagHC, labels = c("No Outlier", "Outlier"))
table(dat$flagHC)
rm(cutoffU); rm(cutoffL)
attributes(dat$flagHC)$label <- "pmmst - HC outlier?"

# summary and plot (no outlier removed)
summary(dat$avgHC)
ggplot(dat, aes(x = avgHC)) + geom_histogram(colour="black", fill="white", binwidth = 1) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 85)) + ggtitle("PMMST: Head circumference (cm)") +
  xlab("Head circumference (cm)")
# excluding the outliers
dat %>% filter(flagHC == "No Outlier") %>% 
  summarise(n = n(), Min = min(avgHC), Mean = mean(avgHC), SD = 
              sd(avgHC), Median = median(avgHC), Q1 = 
              quantile(avgHC, probs = 0.25), Q3 = 
              quantile(avgHC, probs = 0.75), Max = max(avgHC))
ggplot(subset(dat, flagHC == "No Outlier"), aes(x = avgHC)) + 
  geom_histogram(colour="black", fill="white", binwidth = 1) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 80)) + 
  ggtitle("PMMST (No Outliers): Head Circumference (cm)") +
  xlab("HC (cm)")



#### MUAC ####

with(dat, summary(nMUAC1))
with(dat, sd(nMUAC1))
with(dat, summary(nMUAC2))
with(dat, sd(nMUAC2))
# compare the distribution of MUAC 1 and 2
par(mfrow = c(2,1))
with(dat, hist(nMUAC1))
with(dat, hist(nMUAC2))
# generate a variable looking at the difference between the two measures
diff <- dat$nMUAC1 - dat$nMUAC2
summary(diff)
rm(diff)

# average the MUAC measures and create an average MUAC variable
dat$avgMUAC <- (dat$nMUAC1 +  dat$nMUAC2)/2
summary(dat$avgMUAC)
sd(dat$avgMUAC)
par(mfrow = c(1,1))
hist(dat$avgMUAC)

# check for outliers
cutoffU <- median(dat$avgMUAC) + 5*mad(dat$avgMUAC)
cutoffL <- median(dat$avgMUAC) - 5*mad(dat$avgMUAC)
sum(dat$avgMUAC > cutoffU)
sum(dat$avgMUAC < cutoffL)
# no outliers based on the criteria established in the Analysis Plan
rm(cutoffU); rm(cutoffL)

ggplot(dat, aes(x = avgMUAC)) + geom_histogram(colour="black", fill="white", binwidth = 0.8) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 95)) + ggtitle("PMMST: MUAC (cm)") +
  xlab("MUAC (cm)")



#### Chest Circumference ####

with(dat, summary(nChest1))
with(dat, sd(nChest1))
with(dat, summary(nChest2))
with(dat, sd(nChest2))
# compare the distribution of chest circumference
par(mfrow = c(2,1))
with(dat, hist(nChest1))
with(dat, hist(nChest2))
# generate a variable looking at the difference between the two measures
diff <- dat$nChest1 - dat$nChest2
summary(diff)
rm(diff)

# average the chest measures and create an average chest circumference variable
dat$avgChest <- (dat$nChest1 +  dat$nChest2)/2
summary(dat$avgChest)
sd(dat$avgChest)
par(mfrow = c(1,1))
hist(dat$avgChest)

# check for outliers
cutoffU <- median(dat$avgChest) + 5*mad(dat$avgChest)
cutoffL <- median(dat$avgChest) - 5*mad(dat$avgChest)
sum(dat$avgChest > cutoffU)
sum(dat$avgChest < cutoffL)
# no outliers based on the criteria established in the Analysis Plan
rm(cutoffU); rm(cutoffL)

ggplot(dat, aes(x = avgChest)) + geom_histogram(colour="black", fill="white", binwidth = 1.5) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 75)) + ggtitle("PMMST: Chest Circumference (cm)") +
  xlab("Chest Circumference (cm)")


#### Hip Circumference ####

with(dat, summary(nHip1))
with(dat, sd(nHip1))
with(dat, summary(nHip2))
with(dat, sd(nHip2))
# compare the distribution of hip 1 and 2
par(mfrow = c(2,1))
with(dat, hist(nHip1))
with(dat, hist(nHip2))
# generate a variable looking at the difference between the two measures
diff <- dat$nHip1 - dat$nHip2
summary(diff)
rm(diff)

# average the hip measures and create an average hip circumference variable
dat$avgHip <- (dat$nHip1 +  dat$nHip2)/2
summary(dat$avgHip)
sd(dat$avgHip)
par(mfrow = c(1,1))
hist(dat$avgHip)

# check for outliers
cutoffU <- median(dat$avgHip) + 5*mad(dat$avgHip)
cutoffL <- median(dat$avgHip) - 5*mad(dat$avgHip)
sum(dat$avgHip > cutoffU)
sum(dat$avgHip < cutoffL)
# 1 outliers based on the criteria established in the Analysis Plan

# create a flag variable for Hip
dat$flagHip <- ifelse(dat$avgHip > cutoffU, 1, 0)
dat$flagHip <- factor(dat$flagHip, labels = c("No Outlier", "Outlier"))
table(dat$flagHip)
rm(cutoffU); rm(cutoffL)
attributes(dat$flagHip)$label <- "pmmst - Hip circumference outlier?"

# summary and plot (no outlier removed)
summary(dat$avgHip)
ggplot(dat, aes(x = avgHip)) + geom_histogram(colour="black", fill="white", binwidth = 1) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 50)) + ggtitle("PMMST: Hip circumference (cm)") +
  xlab("Hip circumference (cm)")
# excluding the outliers
dat %>% filter(flagHip == "No Outlier") %>% 
  summarise(n = n(), Min = min(avgHip), Mean = mean(avgHip), SD = 
              sd(avgHip), Median = median(avgHip), Q1 = 
              quantile(avgHip, probs = 0.25), Q3 = 
              quantile(avgHip, probs = 0.75), Max = max(avgHip))
ggplot(subset(dat, flagHip == "No Outlier"), aes(x = avgHip)) + 
  geom_histogram(colour="black", fill="white", binwidth = 1) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 50)) + 
  ggtitle("PMMST (No Outliers): Hip circumference (cm)") +
  xlab("Hip Circumference (cm)")



#### Triceps ####

with(dat, summary(nTricepsSFT1))
with(dat, sd(nTricepsSFT1))
with(dat, summary(nTricepsSFT2))
with(dat, sd(nTricepsSFT2))
with(dat, summary(nTricepsSFT3))
with(dat, sd(nTricepsSFT3))

# compare the distribution of triceps 1, 2 and 3
par(mfrow = c(3,1))
with(dat, hist(nTricepsSFT1))
with(dat, hist(nTricepsSFT2))
with(dat, hist(nTricepsSFT3))
# generate a variable looking at the difference between the first and second measures
diff <- dat$nTricepsSFT1 - dat$nTricepsSFT2
summary(diff)
rm(diff)
# generate a variable looking at the difference between the first and third measures
diff <- dat$nTricepsSFT1 - dat$nTricepsSFT3
summary(diff)
rm(diff)
# generate a variable looking at the difference between the second and third measures
diff <- dat$nTricepsSFT2 - dat$nTricepsSFT3
summary(diff)
rm(diff)

# average the triceps measures and create an average triceps variable
dat$avgTriceps <- (dat$nTricepsSFT1 +  dat$nTricepsSFT2 + dat$nTricepsSFT3)/3
summary(dat$avgTriceps)
sd(dat$avgTriceps)
par(mfrow = c(1,1))
hist(dat$avgTriceps)

# check for outliers
cutoffU <- median(dat$avgTriceps) + 5*mad(dat$avgTriceps)
cutoffL <- median(dat$avgTriceps) - 5*mad(dat$avgTriceps)
sum(dat$avgTriceps > cutoffU)
sum(dat$avgTriceps < cutoffL)
# 0 outliers based on the criteria established in the Analysis Plan
rm(cutoffU); rm(cutoffL)

ggplot(dat, aes(x = avgTriceps)) + geom_histogram(colour="black", fill="white", binwidth = 1) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 80)) + ggtitle("PMMST: Triceps (mm)") +
  xlab("Triceps (mm)")



#### Biceps ####

with(dat, summary(nBicepsSFT1))
with(dat, sd(nBicepsSFT1))
with(dat, summary(nBicepsSFT2))
with(dat, sd(nBicepsSFT2))
with(dat, summary(nBicepsSFT3))
with(dat, sd(nBicepsSFT3))
# compare the distribution of biceps 1, 2 and 3
par(mfrow = c(3,1))
with(dat, hist(nBicepsSFT1))
with(dat, hist(nBicepsSFT2))
with(dat, hist(nBicepsSFT3))
# generate a variable looking at the difference between the first and second measures
diff <- dat$nBicepsSFT1 - dat$nBicepsSFT2
summary(diff)
rm(diff)
# generate a variable looking at the difference between the first and third measures
diff <- dat$nBicepsSFT1 - dat$nBicepsSFT3
summary(diff)
rm(diff)

# average the biceps measures and create an average triceps variable
dat$avgBiceps <- (dat$nBicepsSFT1 +  dat$nBicepsSFT2 + dat$nBicepsSFT3)/3
summary(dat$avgBiceps)
sd(dat$avgBiceps)
par(mfrow = c(1,1))
hist(dat$avgBiceps)

# check for outliers
cutoffU <- median(dat$avgBiceps) + 5*mad(dat$avgBiceps)
cutoffL <- median(dat$avgBiceps) - 5*mad(dat$avgBiceps)
sum(dat$avgBiceps > cutoffU)
sum(dat$avgBiceps < cutoffL)
# 0 outliers based on the criteria established in the Analysis Plan
rm(cutoffU); rm(cutoffL)

ggplot(dat, aes(x = avgBiceps)) + geom_histogram(colour="black", fill="white", binwidth = 0.5) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 80)) + ggtitle("PMMST: Biceps (mm)") +
  xlab("Biceps (mm)")



#### Subscapular ####

with(dat, summary(nSubscapSFT1))
with(dat, sd(nSubscapSFT1))
with(dat, summary(nSubscapSFT2))
with(dat, sd(nSubscapSFT2))
with(dat, summary(nSubscapSFT3))
with(dat, sd(nSubscapSFT3))
# compare the distribution of subscapular 1, 2 and 3
par(mfrow = c(3,1))
with(dat, hist(nSubscapSFT1))
with(dat, hist(nSubscapSFT2))
with(dat, hist(nSubscapSFT3))
# generate a variable looking at the difference between the first and second measures
diff <- dat$nSubscapSFT1 - dat$nSubscapSFT2
summary(diff)
rm(diff)
# generate a variable looking at the difference between the first and third measures
diff <- dat$nSubscapSFT1 - dat$nSubscapSFT3
summary(diff)
rm(diff)
# generate a variable looking at the difference between the second and third measures
diff <- dat$nSubscapSFT2 - dat$nSubscapSFT3
summary(diff)
rm(diff)

# average the subscapular measures and create an average subscapular variable
dat$avgSubscap <- (dat$nSubscapSFT1 +  dat$nSubscapSFT2 + dat$nSubscapSFT3)/3
summary(dat$avgSubscap)
sd(dat$avgSubscap)
par(mfrow = c(1,1))
hist(dat$avgSubscap)

# check for outliers
cutoffU <- median(dat$avgSubscap) + 5*mad(dat$avgSubscap)
cutoffL <- median(dat$avgSubscap) - 5*mad(dat$avgSubscap)
sum(dat$avgSubscap > cutoffU)
sum(dat$avgSubscap < cutoffL)
# 0 outliers based on the criteria established in the Analysis Plan
rm(cutoffU); rm(cutoffL)

ggplot(dat, aes(x = avgSubscap)) + geom_histogram(colour="black", fill="white", binwidth = 0.5) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 80)) + ggtitle("PMMST: Subscapular Skinfold (mm)") +
  xlab("Subscapular Skinfold (mm)")



#### Suprailiac ####

with(dat, summary(nSupraSFT1))
with(dat, sd(nSupraSFT1))
with(dat, summary(nSupraSFT2))
with(dat, sd(nSupraSFT2))
with(dat, summary(nSupraSFT3))
with(dat, sd(nSupraSFT3))
# compare the distribution of suprailiac 1, 2 and 3
par(mfrow = c(3,1))
with(dat, hist(nSupraSFT1))
with(dat, hist(nSupraSFT2))
with(dat, hist(nSupraSFT3))
# generate a variable looking at the difference between the first and second measures
diff <- dat$nSupraSFT1 - dat$nSupraSFT2
summary(diff)
rm(diff)
# generate a variable looking at the difference between the first and third measures
diff <- dat$nSupraSFT1 - dat$nSupraSFT3
summary(diff)
rm(diff)
# generate a variable looking at the difference between the second and third measures
diff <- dat$nSupraSFT2 - dat$nSupraSFT3
summary(diff)
rm(diff)

# average the suprailiac measures and create an average subscapular variable
dat$avgSupra <- (dat$nSupraSFT1 +  dat$nSupraSFT2 + dat$nSupraSFT3)/3
summary(dat$avgSupra)
sd(dat$avgSupra)
par(mfrow = c(1,1))
hist(dat$avgSupra)

# check for outliers
cutoffU <- median(dat$avgSupra) + 5*mad(dat$avgSupra)
cutoffL <- median(dat$avgSupra) - 5*mad(dat$avgSupra)
sum(dat$avgSupra > cutoffU)
sum(dat$avgSupra < cutoffL)
# 1 outliers based on the criteria established in the Analysis Plan

# create a flag variable for suprailiac
dat$flagSupra <- ifelse(dat$avgSupra > cutoffU, 1, 0)
dat$flagSupra <- factor(dat$flagSupra, labels = c("No Outlier", "Outlier"))
table(dat$flagSupra)
rm(cutoffU); rm(cutoffL)
attributes(dat$flagSupra)$label <- "pmmst - Suprailiac outlier?"

# summary and plot (no outlier removed)
summary(dat$avgSupra)
ggplot(dat, aes(x = avgSupra)) + geom_histogram(colour="black", fill="white", binwidth = 0.5) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 120)) + ggtitle("PMMST: Suprailiac (mm)") +
  xlab("Suprailiac (mm)")
# excluding the outliers
dat %>% filter(flagSupra == "No Outlier") %>% 
  summarise(n = n(), Min = min(avgSupra), Mean = mean(avgSupra), SD = 
              sd(avgSupra), Median = median(avgSupra), Q1 = 
              quantile(avgSupra, probs = 0.25), Q3 = 
              quantile(avgSupra, probs = 0.75), Max = max(avgSupra))
ggplot(subset(dat, flagSupra == "No Outlier"), aes(x = avgSupra)) + 
  geom_histogram(colour="black", fill="white", binwidth = 0.5) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 120)) + 
  ggtitle("PMMST (No Outliers): Suprailiac (mm)") +
  xlab("Suprailiac (mm)")



#########################################################################################################
#### Multivariate Analysis ##############################################################################
#########################################################################################################

# 1. Check each variable in relation with others
# 2. Check if same measures are consistent



#### Height, Weight and Sitting Height ####

# height and sitting height should have a positive relation. Check
with(dat, plot(avgHeight, avgSitHeight, xlab = "Height (cm)", ylab = "Sitting Height (cm)",
               col = ifelse(avgSitHeight == 52.4, "red", "black"), pch = 19))

# one child with questionable values
dat[(dat$avgHeight < 140 & dat$avgSitHeight < 55),c("cISubjectID", "avgHeight","nSittingHeight1", "nSittingHeight2", 
                                                    "avgSitHeight", "nHeight1", "nHeight2")]
# height and weight should have a positive relation. Check
with(dat, plot(avgWeight, avgHeight))

#### Skinfolds ####
pairs(~ avgTriceps + avgBiceps + avgSubscap + avgSupra + avgWeight,
      data = dat, main = "Skinfolds and weight measures")


#########################################################################################################
#### Derive Measures ####################################################################################
#########################################################################################################

# 1. Compute derived measures listed in the analysis plan
# 2. Check if same measures are consistent
# 3. Check for outlier using the criteria established in the analysis plan (median +/- 5MAD)
# 4. Retain the outliers but create a variable with 0/1 (named flagvariablename) where 1 indicated
#    the outlier based on the analysis plan.


#### Compute BMI ####
dat$bmi <- with(dat, avgWeight/(avgHeight*0.01)^2)
summary(dat$bmi)
hist(dat$bmi)

# check for outliers
cutoffU <- median(dat$bmi, na.rm = T) + 5*mad(dat$bmi, na.rm = T)
cutoffL <- median(dat$bmi, na.rm = T) - 5*mad(dat$bmi, na.rm = T)
sum(dat$bmi > cutoffU, na.rm = T)
sum(dat$bmi < cutoffL, na.rm = T)
# 1 outlier identified based on the criteria in the analysis plan

# create a flag variable for bmi
dat$flagbmi <- ifelse(dat$bmi > cutoffU | dat$bmi < cutoffL, 1, 0)
dat$flagbmi <- factor(dat$flagbmi, labels = c("No Outlier", "Outlier"))
table(dat$flagbmi)
rm(cutoffU); rm(cutoffL)

# summary and plot (no outlier removed)
summary(dat$bmi)
ggplot(dat, aes(x = bmi)) + geom_histogram(colour="black", fill="white", binwidth = 0.8) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 90)) + 
  ggtitle(expression(paste("SARAS Kids: BMI (", kg/m^{2},")"))) +
  xlab(expression(paste("BMI (", kg/m^{2}, ")")))
# excluding the outliers
dat %>% filter(flagbmi == "No Outlier") %>% 
  summarise(n = n(), Min = min(bmi), Mean = mean(bmi), SD = 
              sd(bmi), Median = median(bmi), Q1 = 
              quantile(bmi, probs = 0.25), Q3 = 
              quantile(bmi, probs = 0.75), Max = max(bmi))
ggplot(subset(dat, flagbmi == "No Outlier"), aes(x = bmi)) + 
  geom_histogram(colour="black", fill="white", binwidth = 0.8) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) + 
  ggtitle(expression(paste("SARAS Kids (No Outliers): BMI (", kg/m^{2},")"))) +
  xlab(expression(paste("BMI (", kg/m^{2}, ")")))



#### Compute Leg Length ####
dat$leg <- with(dat, (avgHeight - avgSitHeight))
summary(dat$leg)
hist(dat$leg)
attributes(dat$leg)$label <- "pmmst - leg length (cm)"

# check for outliers
cutoffU <- median(dat$leg, na.rm = T) + 5*mad(dat$leg, na.rm = T)
cutoffL <- median(dat$leg, na.rm = T) - 5*mad(dat$leg, na.rm = T)
sum(dat$leg > cutoffU, na.rm = T)
sum(dat$leg < cutoffL, na.rm = T)
# no outliers were identified

# summary and plot
summary(dat$leg)
ggplot(dat, aes(x = leg)) + geom_histogram(colour="black", fill="white", binwidth = 1.5) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 60)) + 
  ggtitle("PMMST: Leg Length (cm)") +
  xlab("Leg length(cm)")



#### Compute Leg Length to Sitting Height Ratio ####
with(dat, plot(leg, avgSitHeight))
dat$legshtratio <- with(dat, (leg/avgSitHeight))
summary(dat$legshtratio)
hist(dat$legshtratio)

# check for outliers
cutoffU <- median(dat$legshtratio, na.rm = T) + 5*mad(dat$legshtratio, na.rm = T)
cutoffL <- median(dat$legshtratio, na.rm = T) - 5*mad(dat$legshtratio, na.rm = T)
sum(dat$legshtratio > cutoffU, na.rm = T)
sum(dat$legshtratio < cutoffL, na.rm = T)
# no outliers were identified

# summary and plot
summary(dat$legshtratio)
ggplot(dat, aes(x = legshtratio)) + geom_histogram(colour="black", fill="white", binwidth = 0.02) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 45)) + 
  ggtitle("PMMST: Leg Length to Sitting Height") +
  xlab("Leg Length to Sitting Height")


#### Total sum of skinfold ####
dat$totsskin <- with(dat, (avgTriceps + avgBiceps + avgSubscap + avgSupra))
attributes(dat$totsskin)$label <- "pmmst - total sum of skinfolds (mm)"
summary(dat$totsskin)
hist(dat$totsskin)

# check for outliers
cutoffU <- median(dat$totsskin, na.rm = T) + 5*mad(dat$totsskin, na.rm = T)
cutoffL <- median(dat$totsskin, na.rm = T) - 5*mad(dat$totsskin, na.rm = T)
sum(dat$totsskin > cutoffU, na.rm = T)
sum(dat$totsskin < cutoffL, na.rm = T)
# no outliers have been identified

ggplot(dat, aes(x = totsskin)) + geom_histogram(colour="black", fill="white", binwidth = 1.5) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 60)) +
  ggtitle("PMMST Kids: Total sum of skinfolds (mm)") + xlab("Total sum of skinfolds (mm)")



#### Compute SD scores based on WHO criteria ####

# SD score are computed based on the WHO criteria
# The WHO provides on their website an already-made macro for computing z-scores.
# The macro was downloaded (18/05/2017) and saved as who2007_R. Instruction for using the macro can be
# found on Readme.pdf

# Restore reference data sets:
wfawho2007 <- 
  read.table("T:/Chiara/EMPHASIS/Phenotype Data Processing/who2007_R/wfawho2007.txt",
             header=T,sep="",skip=0) 
hfawho2007 <- 
  read.table("T:/Chiara/EMPHASIS/Phenotype Data Processing/who2007_R/hfawho2007.txt",
             header=T,sep="",skip=0) 
bfawho2007 <- 
  read.table("T:/Chiara/EMPHASIS/Phenotype Data Processing/who2007_R/bfawho2007.txt",
             header=T,sep="",skip=0) 

# source the who2007.R file
source("T:/Chiara/EMPHASIS/Phenotype Data Processing/who2007_R/who2007.r")

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
dat$cISubjectID <- toupper(dat$cISubjectID)
dat$cISubjectID <- as.character(dat$cISubjectID)
datsex_age$age <- as.numeric(datsex_age$age)

# merge
dat <- inner_join(dat, datsex_age)
rm(datsex_age)

# take only the variable needed
datZscore <- dat %>% select(cISubjectID, age, sex, avgWeight, avgHeight, bmi)
head(datZscore)

# If sex is character, it must be “m” or “M” for males and “f” or “F” for females.
# Transform sex in character
datZscore$sex <- as.character(datZscore$sex)
# age must be in months not rounded. Transform age in month
datZscore$age <- datZscore$age*12
# height and weight must be in centimeter and kg respectively

# transfrom Zscore in dataframe
datZscore <- data.frame(datZscore)
str(datZscore)

# strip variable of their attributes
attributes(datZscore$age) <- NULL
attributes(datZscore$avgWeight) <- NULL
attributes(datZscore$avgHeight) <- NULL
attributes(datZscore$bmi) <- NULL
str(datZscore)

who2007(FileLab = "PMMSTZscore",
        FilePath = "T:/Chiara/EMPHASIS/Phenotype Data Processing",
        mydf = datZscore, sex = sex, age = age, weight = avgWeight, height = avgHeight)


# remove unneed functions and dataset
rm(bfawho2007); rm(datZscore); rm(hfawho2007); rm(matprev); rm(matz); rm(wfawho2007) 
rm(calc.zbmi); rm(calc.zhfa); rm(calc.zwei); rm(prevnh); rm(prevnh.L); rm(prevph)
rm(prevph.L); rm(rounde); rm(who2007); rm(wmean); rm(wsd)

# Merge dataset with WHO z-score 
datZscore <- read.csv("PMMSTZscore_z.csv")
keep <- c("cISubjectID", "zhfa", "zwfa", "zbfa")
datZscore <- datZscore[, names(datZscore) %in% keep]
datZscore$cISubjectID <- as.character(datZscore$cISubjectID)
dat <- inner_join(dat, datZscore)
dat <- data.frame(dat)
# remove unused data
rm(datZscore); rm(keep)




#### Height for age (WHO z-score) ####

summary(dat$zhfa)
sd(dat$zhfa)
ggplot(dat, aes(x = zhfa)) + geom_histogram(colour="black", fill="white", binwidth = 0.5) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 80)) +
  ggtitle("PMMST: Height-for-age (WHO z-score)") + xlab("Height-for-age (WHO z-score)")

#### Weight for age (WHO z-score) ####

summary(dat$zwfa)
sd(dat$zwfa, na.rm = T)
ggplot(dat, aes(x = zwfa)) + geom_histogram(colour="black", fill="white", binwidth = 0.5) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 75)) +
  ggtitle("PMMST: Weight-for-age (WHO z-score)") + xlab("Weight-for-age (WHO z-score)")


# None of the wfa and hfa have been identified by WHO as outliers


#### Compute the proportion of stunting ####

# stunting is defined as height-for-age less than -2 z-scores.
dat <- data.frame(dat)
dat$stunting <- rep(NA, nrow(dat))
dat[!is.na(dat$zhfa) & dat$zhfa <= -2, "stunting"] <- "Yes"
dat[!is.na(dat$zhfa) & dat$zhfa > -2, "stunting"] <- "No"
addmargins(table(dat$stunting))


#### Compute the proportion of wasting ####

# wasting is defined as bmi-for-age less than -2 z-scores.
dat$wasting <- rep(NA, nrow(dat))
dat[!is.na(dat$zbfa) & dat$zbfa <= -2, "wasting"] <- "Yes"
dat[!is.na(dat$zbfa) & dat$zbfa > -2, "wasting"] <- "No"
addmargins(table(dat$wasting))


#### Compute the proportion of underweight  ####

# wasting is defined as weight-for-age less -2 z-scores.
dat$underweight <- rep(NA, nrow(dat))
dat[!is.na(dat$zwfa) & dat$zwfa <= -2, "underweight"] <- "Yes"
dat[!is.na(dat$zwfa) & dat$zwfa > -2, "underweight"] <- "No"
addmargins(table(dat$underweight))



#########################################################################################################
#### Save Dataset #######################################################################################
#########################################################################################################

# take only relevant variables
dat <- dat %>% dplyr::select(cISubjectID, sex, age, MasterGroupNo, avgWeight, avgHeight, avgSitHeight, 
                             avgWaist, avgHC, avgMUAC, avgChest, avgHip, avgTriceps, 
                             avgBiceps, avgSubscap, avgSupra, bmi, leg, legshtratio, 
                             totsskin, zhfa, zwfa, stunting, wasting, underweight, flagHC, flagHip, 
                             flagSupra, flagbmi)

# rename variables
names(dat) <- c("serno", "csex", "cage", "allocation", "cwt", "cht", "csitht", "cavgWaist", "cavgHC",
                "cavgMUAC", "cavgCC", "cavgHip", "cavgTric", "cavgBic", "cavgSubsc", "cavgSup", "cbmi", "cleg", 
                "clegsithtratio", "ctotsskin", "chfa", "cwfa", "cstunting", "cwasting", "cunderweight",
                "flagcavgHC", "flagcavgHip", "flagcavgSup", "flagcbmi") 

# create a visit variable
dat$cvisit <- 1
attributes(dat$cvisit)$label <- "pmmst - child clinic visit"

# add visit 2 data
datV2 <- read_dta("./Cleaned Dataset/PMMSTGrowthBodyCompCV2v1_20180126.dta")
# transform atomic in factors to ease the merge
datV2$flagcwt <- as_factor(datV2$flagcwt)
datV2$flagcavgCC <- as_factor(datV2$flagcavgCC)
datV2$flagcavgTric <- as_factor(datV2$flagcavgTric)
datV2$flagcavgSubsc <- as_factor(datV2$flagcavgSubsc)
datV2$flagctotsskin <- as_factor(datV2$flagctotsskin)
datV2$flagcavgHC <- as_factor(datV2$flagcavgHC)
datV2$flagcavgSup <- as_factor(datV2$flagcavgSup)
datV2$flagcbmi <- as_factor(datV2$flagcbmi)
# join data together
dat <- full_join(dat, datV2)
rm(datV2)
# sort data by ID and visit
dat <- dat[order(dat$serno, dat$cvisit),]
# for the variables starting with flag, fill NA as no outliers
flagvar <- names(dat)[grep("^flag", names(dat))]
dat[, names(dat) %in% flagvar][is.na(dat[, names(dat) %in% flagvar])] <-  "No Outlier"
# add body composition as measured by DXA (those will be available fro V2 only) and should be NA for V1

# add variable description
attributes(dat$serno)$label <- "pmmst - Subject ID"
attributes(dat$csex)$label <- "pmmst - child sex"
attributes(dat$cage)$label <- "pmmst - child age (years)"
attributes(dat$cstunting)$label <- "pmmst - stunting"
attributes(dat$cwasting)$label <- "pmmst - wasting"
attributes(dat$cunderweight)$label <- "pmmst - underweight"
attributes(dat$allocation)$label <- "pmmst - allocation group"
attributes(dat$cwt)$label <- "pmmst - child weight (kg)"
attributes(dat$cht)$label <- "pmmst - child height (cm)"
attributes(dat$csitht)$label <- "pmmst - child sitting height (cm)"
attributes(dat$cavgWaist)$label <- "pmmst - child average waist (cm)"
attributes(dat$cavgHC)$label <- "pmmst - child average HC (cm)"
attributes(dat$cavgMUAC)$label <- "pmmst - child average MUAC (cm)"
attributes(dat$cavgCC)$label <- "pmmst - child average CC (cm)"
attributes(dat$cavgHip)$label <- "pmmst - child average Hip (cm)"
attributes(dat$cavgTric)$label <- "pmmst - child average triceps (mm)"
attributes(dat$cavgBic)$label <- "pmmst - child average biceps (mm)"
attributes(dat$cavgSubsc)$label <- "pmmst - child average subscapular (mm)"
attributes(dat$cavgSup)$label <- "pmmst - child average suprailiac (mm)"


# get the code for cleaning DXA
source("PMMSTPhenotypeCleaningDXA.R")

# look if the fat dataframe has been imported correctly
str(datFat)
# change ID name to match that in the antro data
names(datFat)[names(datFat) == "ID"] <- "serno"
# take only the variable we need
datFat <- datFat %>% dplyr::select(serno, totfatmass_kg, totleanmass_kg, androidfatmass_kg, gynoidfatmass_kg,  
                                   fmi, lmi, fatper, flagtotfatmass_kg, flagandroidfatmass_kg,  
                                   flaggynoidfatmass_kg, flagfmi)
names(datFat) <- c("serno", "cfatmass", "cleanmass", "candroidfat", "cgynoidfat", "cfmi", "clmi", "cfatper",
                   "flagcfatmass", "flagcandroidfat", "flagcgynoidfat", "flagcfmi")

# create a visit variable
datFat$cvisit <- 2
# merge with the antros
dat <- full_join(dat, datFat)
rm(datFat, datBone)

attributes(dat$cbmi)$label <- "pmmst - bmi (kg/m2)"
attributes(dat$cwfa)$label <- "pmmst - weight for age (WHO z-score)"
attributes(dat$chfa)$label <- "pmmst - height for age (WHO z-score)"
attributes(dat$clegsithtratio)$label <- "pmmst - leg length to sitting height ratio"
attributes(dat$flagcbmi)$label <- "pmmst - bmi outlier?"
attributes(dat$cvisit)$label <- "pmmst - clinic visit"
attributes(dat$ctotsskin)$label <- "pmmst - total sum of skinfolds (mm)"
attributes(dat$cleg)$label <- "pmmst - leg length (cm)"
attributes(dat$flagcavgHC)$label <- "pmmst - HC outlier?"
attributes(dat$flagcavgHip)$label <- "pmmst - hip circumference outlier?"
attributes(dat$flagcavgSup)$label <- "pmmst - suprailiac skinfold outlier?"
attributes(dat$flagcwt)$label <- "pmmst - weight outlier?"
attributes(dat$flagcavgCC)$label <- "pmmst - CC skinfold outlier?"
attributes(dat$flagcavgTric)$label <- "pmmst - triceps outlier?"
attributes(dat$flagcavgSubsc)$label <- "pmmst - subscapular skinfold outlier?"
attributes(dat$flagctotsskin)$label <- "pmmst - sum of skinfolds outlier?"

# change level of csex as Boy and Girl
table(dat$csex)
dat$csex[dat$csex == "M"] <- "Boy"
dat$csex[dat$csex== "F"] <- "Girl"
table(dat$csex)

# save as stata, SPSS and csv
write.table(dat, "./Cleaned Dataset/PMMSTGrowthBodyCompv1_20180131.csv", row.names = FALSE)
write_dta(dat, "./Cleaned Dataset/PMMSTGrowthBodyCompv1_20180131.dta")
write_sav(dat, "./Cleaned Dataset/PMMSTGrowthBodyCompv1_20180131.sav")
