##############################################################################
##### PMMST: Phenotype Data Processing Clinic Visit 2 Anthropometry  #########
##############################################################################

# 12/07/2017

# load packages
library(haven)
library(readxl)
library(data.table)
library(dplyr)
library(tidyr)

# Data were sent by Ayden on 31/05/2017. 
# File wasin a zipped folder. Folder name is Emphasis_GMB_CV_data.zip. File name is
# Emphasis_CV1_AnthropsBP_170531 (excel file)
# On June 7, 2017 queries were sent to Ayden. Ayden came back with answers to most of those queries on
# 30/60/2017. A new version of the data was sent (Emphasis_CV1_AnthropsBP_170630.xlsx)

# load datasets with antrhopometry and close all connections
dat <- read_excel("./PMMST Outcome Data/Emphasis_CV2_Anthrops_170630.xlsx",
                  na = "NULL", col_names = TRUE)

# sanity check on imported dataset
dim(dat)
str(dat)

# Check all subject have a visit date and all have anthropometrics done
addmargins(table(dat$dDateMeasurement))
# visit date format: YYYY-MM-DD
table(dat$cAnthrops_done)

# Transform every ID in upper case
dat$cStudyNo <- toupper(dat$cStudyNo)

#############################
#### Univariate Analysis ####
#############################

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
with(dat, summary(nWeight3))
with(dat, sd(nWeight3))
# compare the distribution of weight
par(mfrow = c(3,1))
with(dat, hist(nWeight1))
with(dat, hist(nWeight2))
with(dat, hist(nWeight3))
# generate a variable looking at the difference between measures
diff <- dat$nWeight1 - dat$nWeight2
summary(diff)
# generate a variable looking at the difference between measures
diff <- dat$nWeight1 - dat$nWeight3
summary(diff)
# generate a variable looking at the difference between measures
diff <- dat$nWeight2 - dat$nWeight3
summary(diff)
rm(diff)

# average the weight measures and create an average weight variable
dat$avgWeight <- (dat$nWeight1 +  dat$nWeight2 + dat$nWeight3)/3
summary(dat$avgWeight)
sd(dat$avgWeight)
par(mfrow = c(1,1))
hist(dat$avgWeight)
attributes(dat$avgWeight)$label <- "pmmst - average weight (kg)"

# check for outliers
cutoffU <- median(dat$avgWeight) + 5*mad(dat$avgWeight)
cutoffL <- median(dat$avgWeight) - 5*mad(dat$avgWeight)
sum(dat$avgWeight > cutoffU)
sum(dat$avgWeight < cutoffL)
# 1 outliers based on the criteria established in the Analysis Plan
(dat[dat$avgWeight > cutoffU, c("cStudyNo", "nWeight1", "nWeight2", "nWeight3")])  

# create a flag variable for Weight
dat$flagWeight <- ifelse(dat$avgWeight < cutoffL | dat$avgWeight > cutoffU, 1, 0)
dat$flagWeight <- factor(dat$flagWeight, labels = c("No Outlier", "Outlier"))
table(dat$flagWeight)
rm(cutoffU); rm(cutoffL)
attributes(dat$flagWeight)$label <- "pmmst - weight outlier?"

# summary and plot (no outlier removed)
summary(dat$avgWeight)
ggplot(dat, aes(x = avgWeight)) + geom_histogram(colour="black", fill="white", binwidth = 1) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 40)) + ggtitle("PMMST: Weight (kg)") +
  xlab("Weight (kg)")
# excluding the outliers
dat %>% filter(flagWeight == "No Outlier") %>% 
  summarise(n = n(), Min = min(avgWeight), Mean = mean(avgWeight), SD = 
              sd(avgWeight), Median = median(avgWeight), Q1 = 
              quantile(avgWeight, probs = 0.25), Q3 = 
              quantile(avgWeight, probs = 0.75), Max = max(avgWeight))
ggplot(subset(dat, flagWeight == "No Outlier"), aes(x = avgWeight)) + 
  geom_histogram(colour="black", fill="white", binwidth = 1) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 40)) + 
  ggtitle("PMMST (No Outliers): Weight (kg)") +
  xlab("Weight (kg)")


#### Height ####

with(dat, summary(nStandheight1))
with(dat, sd(nStandheight1))
with(dat, summary(nStandheight2))
with(dat, sd(nStandheight2))
with(dat, summary(nStandheight3))
with(dat, sd(nStandheight3))
# compare the distribution of weight
par(mfrow = c(3,1))
with(dat, hist(nStandheight1))
with(dat, hist(nStandheight2))
with(dat, hist(nStandheight3))
# generate a variable looking at the difference between measures
diff <- dat$nStandheight1 - dat$nStandheight2
summary(diff)
# check those with difference greater than -2 cm
dat[which(diff < -2), c("cStudyNo", "nStandheight1", "nStandheight2", "nStandheight3")]
# differences were check by the Gambian data group. Those were the values recorded.
# generate a variable looking at the difference between measures
diff <- dat$nStandheight1 - dat$nStandheight3
summary(diff)
# check those with difference greater than -2 cm
dat[which(diff < -2), c("cStudyNo", "nStandheight1", "nStandheight2", "nStandheight3")]
dat[which(diff > 2), c("cStudyNo", "nStandheight1", "nStandheight2", "nStandheight3")]
# differences were check by the Gambian data group. Those were the values recorded.
rm(diff)
diff <- dat$nStandheight2 - dat$nStandheight3
summary(diff)
dat[which(diff > 2), c("cStudyNo", "nStandheight1", "nStandheight2", "nStandheight3")]

# remove those values that are different from the others and from height at V1
dat[dat$cStudyNo == "EMP017D", "nStandheight2"] <- NA
dat[dat$cStudyNo == "EMP099V", "nStandheight2"] <- NA
dat[dat$cStudyNo == "EMP303R", "nStandheight3"] <- NA

# average the height measures and create an average height variable
dat$avgHeight <- rowMeans(dat[,c("nStandheight1", "nStandheight2", "nStandheight3")], na.rm = T)
summary(dat$avgHeight)
sd(dat$avgHeight)
par(mfrow = c(1,1))
hist(dat$avgHeight)
attributes(dat$avgHeight)$label <- "pmmst - child height (cm)"

# check for outliers
cutoffU <- median(dat$avgHeight) + 5*mad(dat$avgHeight)
cutoffL <- median(dat$avgHeight) - 5*mad(dat$avgHeight)
sum(dat$avgHeight > cutoffU)
sum(dat$avgHeight < cutoffL)
# no outliers based on the criteria established in the Analysis Plan
rm(cutoffU); rm(cutoffL)

summary(dat$avgHeight)
ggplot(dat, aes(x = avgHeight)) + geom_histogram(colour="black", fill="white", binwidth = 1.5) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 45)) + ggtitle("PMMST: Height (cm)") +
  xlab("Height (cm)")


#### Sitting Height ####

with(dat, summary(nSittingHeight1))
with(dat, sd(nSittingHeight1))
with(dat, summary(nSittingHeight2))
with(dat, sd(nSittingHeight2))
with(dat, summary(nSittingHeight3))
with(dat, sd(nSittingHeight3))
# compare the distribution of sitting height 
par(mfrow = c(3,1))
with(dat, hist(nSittingHeight1))
with(dat, hist(nSittingHeight2))
with(dat, hist(nSittingHeight3))
dat[which(dat$nSittingHeight1 < 45), c("cStudyNo", "nSittingHeight1", "nSittingHeight2", "nSittingHeight3")]
# values were on the form.
# generate a variable looking at the difference between measures
diff <- dat$nSittingHeight1 - dat$nSittingHeight2
summary(diff)
# generate a variable looking at the difference between measures
diff <- dat$nSittingHeight1 - dat$nSittingHeight3
summary(diff)
rm(diff)
diff <- dat$nSittingHeight2 - dat$nSittingHeight3
summary(diff)


# average the sitting height measures and create an average height variable
dat$avgSitHeight <- (dat$nSittingHeight1 +  dat$nSittingHeight2 + dat$nSittingHeight3)/3
summary(dat$avgSitHeight)
sd(dat$avgSitHeight)
par(mfrow = c(1,1))
hist(dat$avgSitHeight)
attributes(dat$avgSitHeight)$label <- "pmmst - child sitting height (cm)"
dat[which.min(dat$avgSitHeight),]

# sitting height too low (was 62.4 at CV1) -> remove
dat[dat$cStudyNo == "EMP138F","avgSitHeight"] <- NA

# after looking at V1. Set sitting height as missing 
dat[dat$cStudyNo == "EMP192T", c("avgSitHeight")] <- NA

# EMP116L average sitting height is 11 cm lower than visit 1 -> remove visit 2  
dat[dat$cStudyNo == "EMP116L", "avgSitHeight"] <- NA

# check for outliers
cutoffU <- median(dat$avgSitHeight, na.rm = T) + 5*mad(dat$avgSitHeight, na.rm = T)
cutoffL <- median(dat$avgSitHeight, na.rm = T) - 5*mad(dat$avgSitHeight, na.rm = T)
sum(dat$avgSitHeight > cutoffU, na.rm = T)
sum(dat$avgSitHeight < cutoffL, na.rm = T)
# 0 outliers based on the criteria established in the Analysis Plan

# summary and plot 
summary(dat$avgSitHeight)
ggplot(dat, aes(x = avgSitHeight)) + geom_histogram(colour="black", fill="white", binwidth = 1) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 50)) + ggtitle("PMMST: Sitting Height (cm)") +
  xlab("Sitting Height (cm)")


#### Waist Circumference ####

with(dat, summary(nWaistCircumference1))
with(dat, sd(nWaistCircumference1))
with(dat, summary(nWaistCircumference2))
with(dat, sd(nWaistCircumference2))
with(dat, summary(nWaistCircumference3))
with(dat, sd(nWaistCircumference3))
# compare the distribution of waist circumference
par(mfrow = c(3,1))
with(dat, hist(nWaistCircumference1))
with(dat, hist(nWaistCircumference2))
with(dat, hist(nWaistCircumference3))

# generate a variable looking at the difference between measures
diff <- dat$nWaistCircumference1 - dat$nWaistCircumference2
summary(diff)
rm(diff)
diff <- dat$nWaistCircumference1 - dat$nWaistCircumference3
summary(diff)
rm(diff)
diff <- dat$nWaistCircumference2 - dat$nWaistCircumference3
summary(diff)
rm(diff)

# average the waist circumference measures and create an average waist variable
dat$avgWaist <- (dat$nWaistCircumference1 +  dat$nWaistCircumference2 + dat$nWaistCircumference3)/3
summary(dat$avgWaist)
sd(dat$avgWaist)
par(mfrow = c(1,1))
hist(dat$avgWaist)
attributes(dat$avgWaist)$label <- "pmmst - waist circumference (cm)"

# check for outliers
cutoffU <- median(dat$avgWaist) + 5*mad(dat$avgWaist)
cutoffL <- median(dat$avgWaist) - 5*mad(dat$avgWaist)
sum(dat$avgWaist > cutoffU)
sum(dat$avgWaist < cutoffL)
# No outliers were identified

# summary and plot 
summary(dat$avgWaist)
ggplot(dat, aes(x = avgWaist)) + geom_histogram(colour="black", fill="white", binwidth = 1) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 55)) + ggtitle("PMMST: Waist Circumference (cm)") +
  xlab("Waist Circumference (cm)")



#### Head Circumferece ####

with(dat, summary(nHeadCircumference1))
with(dat, sd(nHeadCircumference1))
with(dat, summary(nHeadCircumference2))
with(dat, sd(nHeadCircumference2))
with(dat, summary(nHeadCircumference3))
with(dat, sd(nHeadCircumference3))
# compare the distribution of HC 
par(mfrow = c(3,1))
with(dat, hist(nHeadCircumference1))
with(dat, hist(nHeadCircumference2))
with(dat, hist(nHeadCircumference3))
# generate a variable looking at the difference between measures
diff <- dat$nHeadCircumference1 - dat$nHeadCircumference2
summary(diff)
dat[which(diff < -5), c("cStudyNo", "nHeadCircumference1", 
                                       "nHeadCircumference2", "nHeadCircumference3")]

# Gambia team confirmed the values.
# second measure of HC does not agree with first and third or with HC measured at CV1 -> remove!
dat[dat$cStudyNo == "EMP012M", "nHeadCircumference2"] <- NA

# generate a variable looking at the difference between measures
diff <- dat$nHeadCircumference1 - dat$nHeadCircumference3
summary(diff)
dat[which(diff > 5), c("cStudyNo", "nHeadCircumference1", 
                        "nHeadCircumference2", "nHeadCircumference3")]
# Gambia team confirmed the values.
# third measure of HC does not agree with first and second or with HC measured at CV1 -> remove!
dat[dat$cStudyNo == "EMP255J", "nHeadCircumference3"] <- NA

# generate a variable looking at the difference between measures
diff <- dat$nHeadCircumference2 - dat$nHeadCircumference3
summary(diff)
rm(diff)

# average the head circumference measures and create an average head variable
dat$avgHC <- rowMeans(dat[,c("nHeadCircumference1", "nHeadCircumference2", "nHeadCircumference3")], na.rm = T)
summary(dat$avgHC)
sd(dat$avgHC)
par(mfrow = c(1,1))
hist(dat$avgHC)
attributes(dat$avgHC)$label <- "pmmst - head circumference (cm)"

# check for outliers
cutoffU <- median(dat$avgHC) + 5*mad(dat$avgHC)
cutoffL <- median(dat$avgHC) - 5*mad(dat$avgHC)
sum(dat$avgHC > cutoffU)
sum(dat$avgHC < cutoffL)
# 2 outliers were identified

# create a flag variable for head circumference
dat$flagHC <- ifelse(dat$avgHC < cutoffL | dat$avgHC > cutoffU, 1, 0)
dat$flagHC <- factor(dat$flagHC, labels = c("No Outlier", "Outlier"))
table(dat$flagHC)
rm(cutoffU); rm(cutoffL)
attributes(dat$flagHC)$label <- "pmmst - HC outlier?"

# summary and plot (no outlier removed)
summary(dat$avgHC)
ggplot(dat, aes(x = avgHC)) + geom_histogram(colour="black", fill="white", binwidth = 1) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 80)) + ggtitle("PMMST: Head Circumference (cm)") +
  xlab("Head Circumference (cm)")
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
  xlab("Head Circumference (cm)")


#### MUAC ####

with(dat, summary(nMUAC1))
with(dat, sd(nMUAC1))
with(dat, summary(nMUAC2))
with(dat, sd(nMUAC2))
with(dat, summary(nMUAC3))
with(dat, sd(nMUAC3))
# compare the distribution of MUAC
par(mfrow = c(3,1))
with(dat, hist(nMUAC1))
with(dat, hist(nMUAC2))
with(dat, hist(nMUAC3))
# value was confirmed by the Gambia team.
# set MUAC = 60 to NA (implausible values)
dat[which(dat$nMUAC1 > 30), c("nMUAC1", "nMUAC2", "nMUAC3")] <- NA
dat[dat$cStudyNo == "EMP192T", c("nMUAC1", "nMUAC2", "nMUAC3")]

# generate a variable looking at the difference between measures
diff <- dat$nMUAC1 - dat$nMUAC2
summary(diff)
rm(diff)
diff <- dat$nMUAC1 - dat$nMUAC3
summary(diff)
diff <- dat$nMUAC2 - dat$nMUAC3
summary(diff)
rm(diff)

# average the MUAC measures and create an average MUAC variable
dat$avgMUAC <- (dat$nMUAC1 +  dat$nMUAC2 + dat$nMUAC3)/3
summary(dat$avgMUAC)
sd(dat$avgMUAC, na.rm = T)
par(mfrow = c(1,1))
hist(dat$avgMUAC)
attributes(dat$avgMUAC)$label <- "pmmst - child MUAC (cm)"

# check for outliers
cutoffU <- median(dat$avgMUAC, na.rm = T) + 5*mad(dat$avgMUAC, na.rm = T)
cutoffL <- median(dat$avgMUAC, na.rm = T) - 5*mad(dat$avgMUAC, na.rm = T)
sum(dat$avgMUAC > cutoffU, na.rm = T)
sum(dat$avgMUAC < cutoffL, na.rm = T)
# 0 outlier identified based on the criteria in the Analysis Plan

# summary and plot 
summary(dat$avgMUAC)
ggplot(dat, aes(x = avgMUAC)) + geom_histogram(colour="black", fill="white", binwidth = 1) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 80)) + ggtitle("PMMST: MUAC (cm)") +
  xlab("MUAC (cm)")


#### Chest Circumference ####

with(dat, summary(nChestCircumference1))
with(dat, sd(nChestCircumference1))
with(dat, summary(nChestCircumference2))
with(dat, sd(nChestCircumference2))
with(dat, summary(nChestCircumference3))
with(dat, sd(nChestCircumference3))
# compare the distribution of chest circumference
par(mfrow = c(3,1))
with(dat, hist(nChestCircumference1))
with(dat, hist(nChestCircumference2))
with(dat, hist(nChestCircumference3))
dat[which(dat$nChestCircumference1 > 80), 
    c("cStudyNo", "nChestCircumference1", "nChestCircumference2", "nChestCircumference3", "avgWeight")]
# data were confirmed by the Gambia team
# CC too high also compared to V1 -> set to mising
dat[dat$cStudyNo == "EMP010A", c("nChestCircumference1", "nChestCircumference2", "nChestCircumference3")] <- NA
dat[dat$cStudyNo == "EMP200L", c("nChestCircumference1", "nChestCircumference2", "nChestCircumference3")] <- NA
dat[dat$cStudyNo == "EMP228P", c("nChestCircumference1", "nChestCircumference2", "nChestCircumference3")] <- NA

#generate a variable looking at the difference between measures
diff <- dat$nChestCircumference1 - dat$nChestCircumference2
summary(diff)
rm(diff)
diff <- dat$nChestCircumference1 - dat$nChestCircumference3
summary(diff)
dat[which(diff < -8), 
    c("cStudyNo", "nChestCircumference1", "nChestCircumference2", "nChestCircumference3")]
# third measure too high -> remove
dat[dat$cStudyNo == "EMP158N", c("nChestCircumference3")] <- NA
diff <- dat$nChestCircumference2 - dat$nChestCircumference3
summary(diff)
rm(diff)

# average the chest measures and create an average chest circumference variable
dat$avgChest <- rowMeans(dat[,c("nChestCircumference1", "nChestCircumference2", "nChestCircumference3")], na.rm = T)
summary(dat$avgChest)
sd(dat$avgChest, na.rm = T)
par(mfrow = c(1,1))
hist(dat$avgChest)

# check for outliers
cutoffU <- median(dat$avgChest, na.rm = T) + 5*mad(dat$avgChest, na.rm = T)
cutoffL <- median(dat$avgChest, na.rm = T) - 5*mad(dat$avgChest, na.rm = T)
sum(dat$avgChest > cutoffU, na.rm = T)
sum(dat$avgChest < cutoffL, na.rm = T)
# 1 outliers were identified based on the criteria in the analysis plan

# create a flag variable for chest circumference
dat$flagChest <- ifelse(dat$avgChest < cutoffL | dat$avgChest > cutoffU, 1, 0)
dat$flagChest <- factor(dat$flagChest, labels = c("No Outlier", "Outlier"))
table(dat$flagChest)
rm(cutoffU); rm(cutoffL)
attributes(dat$flagChest)$label <- "pmmst - CC outlier?"

# summary and plot (no outlier removed)
summary(dat$avgChest)
ggplot(dat, aes(x = avgChest)) + geom_histogram(colour="black", fill="white", binwidth = 1) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 55)) + ggtitle("PMMST: Chest Circumference (cm)") +
  xlab("Chest Circumference (cm)")
# excluding the outliers
dat %>% filter(flagChest == "No Outlier") %>% 
  summarise(n = n(), Min = min(avgChest), Mean = mean(avgChest), SD = 
              sd(avgChest), Median = median(avgChest), Q1 = 
              quantile(avgChest, probs = 0.25), Q3 = 
              quantile(avgChest, probs = 0.75), Max = max(avgChest))
ggplot(subset(dat, flagChest == "No Outlier"), aes(x = avgChest)) + 
  geom_histogram(colour="black", fill="white", binwidth = 1) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 50)) + 
  ggtitle("PMMST (No Outliers): Chest Circumference (cm)") +
  xlab("Chest Circumference (cm)")



#### Hip Circumference ####

with(dat, summary(nHipCircumference1))
with(dat, sd(nHipCircumference1))
with(dat, summary(nHipCircumference2))
with(dat, sd(nHipCircumference2))
with(dat, summary(nHipCircumference3))
with(dat, sd(nHipCircumference3))
# compare the distribution of hip 
par(mfrow = c(3,1))
with(dat, hist(nHipCircumference1))
with(dat, hist(nHipCircumference2))
with(dat, hist(nHipCircumference3))
dat[which(dat$nHipCircumference1 > 80), c("cStudyNo", "nHipCircumference1", 
                                          "nHipCircumference2", "nHipCircumference3")]
# data were confirmed by the Gambia team.
# values too high -> remove!
dat[which(dat$nHipCircumference1 > 80), c("nHipCircumference1", 
                                          "nHipCircumference2", "nHipCircumference3")] <- NA
#generate a variable looking at the difference between measures
diff <- dat$nHipCircumference1 - dat$nHipCircumference2
summary(diff)
rm(diff)
diff <- dat$nHipCircumference1 - dat$nHipCircumference3
summary(diff)
diff <- dat$nHipCircumference2 - dat$nHipCircumference3
summary(diff)
rm(diff)

# average the hip measures and create an average hip circumference variable
dat$avgHip <- rowMeans(dat[,c("nHipCircumference1", "nHipCircumference2", "nHipCircumference3")], na.rm = T)
summary(dat$avgHip)
sd(dat$avgHip, na.rm = T)
par(mfrow = c(1,1))
hist(dat$avgHip)

# check for outliers
cutoffU <- median(dat$avgHip, na.rm = T) + 5*mad(dat$avgHip, na.rm = T)
cutoffL <- median(dat$avgHip, na.rm = T) - 5*mad(dat$avgHip, na.rm = T)
sum(dat$avgHip > cutoffU, na.rm = T)
sum(dat$avgHip < cutoffL, na.rm = T)
# 0 outliers identified based on the criteria in the analysis plan

# summary and plot 
summary(dat$avgHip)
ggplot(dat, aes(x = avgHip)) + geom_histogram(colour="black", fill="white", binwidth = 1) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 40)) + ggtitle("PMMST: Hip Circumference (cm)") +
  xlab("Hip Circumference (cm)")




#### Triceps ####

with(dat, summary(nTricepSkinfold1))
with(dat, sd(nTricepSkinfold1))
with(dat, summary(nTricepSkinfold2))
with(dat, sd(nTricepSkinfold2))
with(dat, summary(nTricepSkinfold3))
with(dat, sd(nTricepSkinfold3))
# compare the distribution of triceps
par(mfrow = c(3,1))
with(dat, hist(nTricepSkinfold1))
with(dat, hist(nTricepSkinfold2))
with(dat, hist(nTricepSkinfold3))
# generate a variable looking at the difference between measures
diff <- dat$nTricepSkinfold1 - dat$nTricepSkinfold2
summary(diff)
# generate a variable looking at the difference between measures
diff <- dat$nTricepSkinfold1 - dat$nTricepSkinfold3
summary(diff)
# generate a variable looking at the difference between measures
diff <- dat$nTricepSkinfold2 - dat$nTricepSkinfold3
summary(diff)
rm(diff)

# average the triceps measures and create an average triceps variable
dat$avgTriceps <- (dat$nTricepSkinfold1 +  dat$nTricepSkinfold2 + dat$nTricepSkinfold3)/3
summary(dat$avgTriceps)
sd(dat$avgTriceps)
par(mfrow = c(1,1))
hist(dat$avgTriceps)

# check for outliers
cutoffU <- median(dat$avgTriceps) + 5*mad(dat$avgTriceps)
cutoffL <- median(dat$avgTriceps) - 5*mad(dat$avgTriceps)
sum(dat$avgTriceps > cutoffU)
sum(dat$avgTriceps < cutoffL)
# 2 outliers based on the criteria established in the Analysis Plan

# create a flag variable for triceps
dat$flagTriceps <- ifelse(dat$avgTriceps < cutoffL | dat$avgTriceps > cutoffU, 1, 0)
dat$flagTriceps <- factor(dat$flagTriceps, labels = c("No Outlier", "Outlier"))
table(dat$flagTriceps)
rm(cutoffU); rm(cutoffL)
attributes(dat$flagTriceps)$label <- "pmmst - triceps outlier?"

# summary and plot (no outlier removed)
summary(dat$avgTriceps)
ggplot(dat, aes(x = avgTriceps)) + geom_histogram(colour="black", fill="white", binwidth = 1) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 60)) + ggtitle("PMMST: Triceps (mm)") +
  xlab("Triceps (mm)")
# excluding the outliers
dat %>% filter(flagTriceps == "No Outlier") %>% 
  summarise(n = n(), Min = min(avgTriceps), Mean = mean(avgTriceps), SD = 
              sd(avgTriceps), Median = median(avgTriceps), Q1 = 
              quantile(avgTriceps, probs = 0.25), Q3 = 
              quantile(avgTriceps, probs = 0.75), Max = max(avgTriceps))
ggplot(subset(dat, flagTriceps == "No Outlier"), aes(x = avgTriceps)) + 
  geom_histogram(colour="black", fill="white", binwidth = 1) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 60)) + 
  ggtitle("PMMST (No Outliers): Triceps (mm)") +
  xlab("Triceps (mm)")



#### Biceps ####

with(dat, summary(nBicepSkinfold1))
with(dat, sd(nBicepSkinfold1))
with(dat, summary(nBicepSkinfold2))
with(dat, sd(nBicepSkinfold2))
with(dat, summary(nBicepSkinfold3))
with(dat, sd(nBicepSkinfold3))
# compare the distribution of biceps
par(mfrow = c(3,1))
with(dat, hist(nBicepSkinfold1))
with(dat, hist(nBicepSkinfold2))
with(dat, hist(nBicepSkinfold3))
# generate a variable looking at the difference between measures
diff <- dat$nBicepSkinfold1 - dat$nBicepSkinfold2
summary(diff)
# generate a variable looking at the difference between measures
diff <- dat$nBicepSkinfold1 - dat$nBicepSkinfold3
summary(diff)
# generate a variable looking at the difference between measures
diff <- dat$nBicepSkinfold2 - dat$nBicepSkinfold3
summary(diff)
rm(diff)

# average the biceps measures and create an average biceps variable
dat$avgBiceps <- (dat$nBicepSkinfold1 +  dat$nBicepSkinfold2 + dat$nBicepSkinfold3)/3
summary(dat$avgBiceps)
sd(dat$avgBiceps)
par(mfrow = c(1,1))
hist(dat$avgBiceps)

# check for outliers
cutoffU <- median(dat$avgBiceps) + 5*mad(dat$avgBiceps)
cutoffL <- median(dat$avgBiceps) - 5*mad(dat$avgBiceps)
sum(dat$avgBiceps > cutoffU)
sum(dat$avgBiceps < cutoffL)
# 0 outlier based on the criteria established in the Analysis Plan
rm(cutoffU); rm(cutoffL)

# summary and plot
summary(dat$avgBiceps)
ggplot(dat, aes(x = avgBiceps)) + geom_histogram(colour="black", fill="white", binwidth = 1) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 80)) + ggtitle("PMMST: Biceps (mm)") +
  xlab("Biceps (mm)")


#### Subscapular ####

with(dat, summary(nSubscapularSkinfold1))
with(dat, sd(nSubscapularSkinfold1))
with(dat, summary(nSubscapularSkinfold2))
with(dat, sd(nSubscapularSkinfold2))
with(dat, summary(nSubscapularSkinfold3))
with(dat, sd(nSubscapularSkinfold3))
# compare the distribution of subscapular 1, 2 and 3
par(mfrow = c(3,1))
with(dat, hist(nSubscapularSkinfold1))
with(dat, hist(nSubscapularSkinfold2))
with(dat, hist(nSubscapularSkinfold3))
# generate a variable looking at the difference between measures
diff <- dat$nSubscapularSkinfold1 - dat$nSubscapularSkinfold2
summary(diff)
# generate a variable looking at the difference between measures
diff <- dat$nSubscapularSkinfold1 - dat$nSubscapularSkinfold3
summary(diff)
# generate a variable looking at the difference between measures
diff <- dat$nSubscapularSkinfold2 - dat$nSubscapularSkinfold3
summary(diff)
rm(diff)

# average the subscapular measures and create an average subscapular variable
dat$avgSubscap <- (dat$nSubscapularSkinfold1 +  dat$nSubscapularSkinfold2 + dat$nSubscapularSkinfold3)/3
summary(dat$avgSubscap)
sd(dat$avgSubscap)
par(mfrow = c(1,1))
hist(dat$avgSubscap)

# check for outliers
cutoffU <- median(dat$avgSubscap) + 5*mad(dat$avgSubscap)
cutoffL <- median(dat$avgSubscap) - 5*mad(dat$avgSubscap)
sum(dat$avgSubscap > cutoffU)
sum(dat$avgSubscap < cutoffL)
# 1 outlier based on the criteria established in the Analysis Plan

# create a flag variable for subscapular
dat$flagSubscap <- ifelse(dat$avgSubscap < cutoffL | dat$avgSubscap > cutoffU, 1, 0)
dat$flagSubscap <- factor(dat$flagSubscap, labels = c("No Outlier", "Outlier"))
table(dat$flagSubscap)
rm(cutoffU); rm(cutoffL)
attributes(dat$flagSubscap)$label <- "pmmst - subscapular outlier?"

# summary and plot (no outlier removed)
summary(dat$avgSubscap)
ggplot(dat, aes(x = avgSubscap)) + geom_histogram(colour="black", fill="white", binwidth = 1) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 120)) + ggtitle("PMMST: Subscapular Skinfold (mm)") +
  xlab("Subscapular Skinfold (mm)")
# excluding the outliers
dat %>% filter(flagSubscap == "No Outlier") %>% 
  summarise(n = n(), Min = min(avgSubscap), Mean = mean(avgSubscap), SD = 
              sd(avgTriceps), Median = median(avgSubscap), Q1 = 
              quantile(avgSubscap, probs = 0.25), Q3 = 
              quantile(avgSubscap, probs = 0.75), Max = max(avgSubscap))
ggplot(subset(dat, flagSubscap == "No Outlier"), aes(x = avgSubscap)) + 
  geom_histogram(colour="black", fill="white", binwidth = 1) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 120)) + 
  ggtitle("PMMST (No Outliers):  Subscapular Skinfold (mm)") +
  xlab(" Subscapular Skinfold (mm)")


#### Suprailiac ####

with(dat, summary(nSuprailiacSkinfold1))
with(dat, sd(nSuprailiacSkinfold1))
with(dat, summary(nSuprailiacSkinfold2))
with(dat, sd(nSuprailiacSkinfold2))
with(dat, summary(nSuprailiacSkinfold3))
with(dat, sd(nSuprailiacSkinfold3))
# compare the distribution of suprailiac 1, 2 and 3
par(mfrow = c(3,1))
with(dat, hist(nSuprailiacSkinfold1))
with(dat, hist(nSuprailiacSkinfold2))
with(dat, hist(nSuprailiacSkinfold3))
# generate a variable looking at the difference between measures
diff <- dat$nSuprailiacSkinfold1 - dat$nSuprailiacSkinfold2
summary(diff)
# generate a variable looking at the difference between measures
diff <- dat$nSuprailiacSkinfold1 - dat$nSuprailiacSkinfold3
summary(diff)
# generate a variable looking at the difference between measures
diff <- dat$nSuprailiacSkinfold2 - dat$nSuprailiacSkinfold3
summary(diff)
rm(diff)

# average the suprailiac measures and create an average suprailiac variable
dat$avgSupra <- (dat$nSuprailiacSkinfold1 +  dat$nSuprailiacSkinfold2 + dat$nSuprailiacSkinfold3)/3
summary(dat$avgSupra)
sd(dat$avgSupra)
par(mfrow = c(1,1))
hist(dat$avgSupra)

# check for outliers
cutoffU <- median(dat$avgSupra) + 5*mad(dat$avgSupra)
cutoffL <- median(dat$avgSupra) - 5*mad(dat$avgSupra)
sum(dat$avgSupra > cutoffU)
sum(dat$avgSupra < cutoffL)
# 4 outliers based on the criteria established in the Analysis Plan

# create a flag variable for suprailiac
dat$flagSupra <- ifelse(dat$avgSupra > cutoffU, 1, 0)
dat$flagSupra <- factor(dat$flagSupra, labels = c("No Outlier", "Outlier"))
table(dat$flagSupra)
rm(cutoffU); rm(cutoffL)
attributes(dat$flagSupra)$label <- "pmmst - Suprailiac outlier?"

# summary and plot (no outlier removed)
summary(dat$avgSupra)
ggplot(dat, aes(x = avgSupra)) + geom_histogram(colour="black", fill="white", binwidth = 0.5) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 90)) + ggtitle("PMMST: Suprailiac (mm)") +
  xlab("Suprailiac (mm)")
# excluding the outliers
dat %>% filter(flagSupra == "No Outlier") %>% 
  summarise(n = n(), Min = min(avgSupra), Mean = mean(avgSupra), SD = 
              sd(avgSupra), Median = median(avgSupra), Q1 = 
              quantile(avgSupra, probs = 0.25), Q3 = 
              quantile(avgSupra, probs = 0.75), Max = max(avgSupra))
ggplot(subset(dat, flagSupra == "No Outlier"), aes(x = avgSupra)) + 
  geom_histogram(colour="black", fill="white", binwidth = 0.5) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 90)) + 
  ggtitle("PMMST (No Outliers): Suprailiac (mm)") +
  xlab("Suprailiac (mm)")



#########################################################################################################
#### Multivariate Analysis ##############################################################################
#########################################################################################################

# 1. Check each variable in relation with others
# 2. Check if same measures are consistent



#### Height, Weight and Sitting Height ####

# height and sitting height should have a positive relation. Check
with(dat, plot(avgHeight, avgSitHeight, xlab = "Height (cm)", ylab = "Sitting Height (cm)"))

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
attributes(dat$bmi)$label <- "pmmst - bmi (kg/m2)"

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
attributes(dat$flagbmi)$label <- "pmmst - bmi outlier?"


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
  scale_y_continuous(expand = c(0, 0), limits = c(0, 80)) + 
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
attributes(dat$legshtratio)$label <- "pmmst - leg length to sitting height ratio"

# check for outliers
cutoffU <- median(dat$legshtratio, na.rm = T) + 5*mad(dat$legshtratio, na.rm = T)
cutoffL <- median(dat$legshtratio, na.rm = T) - 5*mad(dat$legshtratio, na.rm = T)
sum(dat$legshtratio > cutoffU, na.rm = T)
sum(dat$legshtratio < cutoffL, na.rm = T)
# no outliers were identified

# summary and plot
summary(dat$legshtratio)
ggplot(dat, aes(x = legshtratio)) + geom_histogram(colour="black", fill="white", binwidth = 0.02) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 60)) + 
  ggtitle("SARAS Kids: Leg Length to Sitting Height") +
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
# 2 outliers have been identified

# create a flag variable for skinfolds
dat$flagtotsskin <- ifelse(dat$totsskin > cutoffU | dat$totsskin < cutoffL, 1, 0)
dat$flagtotsskin <- factor(dat$flagtotsskin, labels = c("No Outlier", "Outlier"))
table(dat$flagtotsskin)
rm(cutoffU); rm(cutoffL)
attributes(dat$flagtotsskin)$label <- "pmmst - total sum of skinfolds outlier?"

# summary and plot (no outlier removed)
summary(dat$totsskin)
ggplot(dat, aes(x = totsskin)) + geom_histogram(colour="black", fill="white", binwidth = 2) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 60)) + 
  ggtitle("PMMST: Total Sum of Skinfolds (mm)") +
  xlab("Total Sum of Skinfold (mm)")
# excluding the outliers
dat %>% filter(flagtotsskin == "No Outlier") %>% 
  summarise(n = n(), Min = min(totsskin), Mean = mean(totsskin), SD = 
              sd(totsskin), Median = median(totsskin), Q1 = 
              quantile(totsskin, probs = 0.25), Q3 = 
              quantile(totsskin, probs = 0.75), Max = max(totsskin))
ggplot(subset(dat, flagtotsskin == "No Outlier"), aes(x = totsskin)) + 
  geom_histogram(colour="black", fill="white", binwidth = 2) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 60)) + 
  ggtitle("PMMST (No Outliers): Total Sum of Skinfolds (mm)") +
  xlab("Total Sum of Skinfold (mm)")



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
datsex_age <- datsex_age %>% dplyr::select(cSubjectID, Sex, DateOfBirth, MasterGroupNo)
# rename cSubjectID to match the name in dat
names(datsex_age)[names(datsex_age) == "cSubjectID"] <- "cStudyNo"
names(datsex_age)[names(datsex_age) == "Sex"] <- "sex"
dat$cStudyNo <- toupper(dat$cStudyNo)
dat$cStudyNo <- as.character(dat$cStudyNo)

# merge
dat <- inner_join(dat, datsex_age)
rm(datsex_age)

# compute age 
dat$age <- with(dat, (dDateMeasurement - DateOfBirth)/365.25)
dat$age <- as.numeric(dat$age)
summary(dat$age)

# take only the variable needed
datZscore <- dat %>% select(cStudyNo, age, sex, avgWeight, avgHeight, bmi)
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

who2007(FileLab = "PMMSTZscoreV2",
        FilePath = "T:/Chiara/EMPHASIS/Phenotype Data Processing",
        mydf = datZscore, sex = sex, age = age, weight = avgWeight, height = avgHeight)


# remove unneed functions and dataset
rm(bfawho2007); rm(datZscore); rm(hfawho2007); rm(matprev); rm(matz); rm(wfawho2007) 
rm(calc.zbmi); rm(calc.zhfa); rm(calc.zwei); rm(prevnh); rm(prevnh.L); rm(prevph)
rm(prevph.L); rm(rounde); rm(who2007); rm(wmean); rm(wsd)

# Merge dataset with WHO z-score 
datZscore <- read.csv("PMMSTZscoreV2_z.csv")
keep <- c("cStudyNo", "zhfa", "zwfa", "zbfa")
datZscore <- datZscore[, names(datZscore) %in% keep]
datZscore$cISubjectID <- as.character(datZscore$cStudyNo)
dat <- inner_join(dat, datZscore)
dat <- data.frame(dat)
# remove unused data
rm(datZscore); rm(keep)

attributes(dat$zwfa)$label <- "pmmst - weight for age (WHO z-score)"
attributes(dat$zhfa)$label <- "pmmst - height for age (WHO z-score)"
attributes(dat$zbfa)$label <- "pmmst - bmi for age (WHO z-score)"


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
dat <- dat %>% dplyr::select(cStudyNo, sex, age, MasterGroupNo, avgWeight, avgHeight, avgSitHeight, 
                             avgWaist, avgHC, avgMUAC, avgChest, avgHip, avgTriceps, 
                             avgBiceps, avgSubscap, avgSupra, bmi, leg, legshtratio, 
                             totsskin, zhfa, zwfa, stunting, wasting, underweight, flagWeight, flagHC, 
                             flagChest, flagTriceps, flagSubscap, flagSupra, flagbmi, flagtotsskin)

# rename variables
names(dat) <- c("serno", "csex", "cage", "allocation", "cwt", "cht", "csitht", "cavgWaist", "cavgHC",
                "cavgMUAC", "cavgCC", "cavgHip", "cavgTric", "cavgBic", "cavgSubsc", "cavgSup", "cbmi", "cleg", 
                "clegsithtratio", "ctotsskin", "chfa", "cwfa", "cstunting", "cwasting", "cunderweight",
                "flagcwt", "flagcavgHC", "flagcavgCC", "flagcavgTric", "flagcavgSubsc", "flagcavgSup",
                "flagcbmi", "flagctotsskin") 

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


# create a visit variable
dat$cvisit <- 2
attributes(dat$cvisit)$label <- "pmmst - child clinic visit"

# save as stata
write_dta(dat, "./Cleaned Dataset/PMMSTGrowthBodyCompCV2v1_20180126.dta")