###############################################################################
##### SARAS Kids: Phenotype Data Processing Cardio-metabolic risk markers #####
###############################################################################

# 2018/04/09

###############################################################################
##### Blood Pressure data #####################################################
###############################################################################


# Dataset needed for the cleaning process:
# sarascp01.dta (contains blood pressure data)
# sarascp01extra.dta (contains blood pressure data)
# SARASKidsBloodCountPP.dta

# load packages
library(haven)
library(dplyr)
library(data.table)
library(ggplot2)
library(broom)
library(Hmisc)
library(readxl)

# Data on blood pressure data were sent by Patsy on 05/06/2017. 
# File name is sarascp01.dta (STATA v14 file). 
# Only children in the per-protocol analysis are processed.

dat <- read_dta("sarascp01.dta")
attributes(dat$serno) <- NULL
dat <- subset(dat, select = -sex)
# check for duplicates
length(unique(dat$serno))

# add information on sex, age, allocation group
datList <- read_dta("SARASKidsBloodCountPP.dta")
datList <- datList %>% dplyr::select(serno, allocation, DateTest, BloodTestDate, age, sex)
datList$sex <- factor(datList$sex, labels = c("Boy", "Girl"))
datList$allocation <- factor(datList$allocation, labels = c("B", "A"))
# check if all sernos in datList are in dat
datList$serno[!(datList$serno %in% dat$serno)] 
# 11 sernos are not in dat:
# 981, 1536, 1548, 1939, 2587, 4250, 4500, 4513, 4580, 4701, 5218

# Vanessa and Harshad looked for the 7 sernos and create an extra dataset with their 
# anthropometry and blood pressure measures

datextra <- read_dta("sarascp01extra.dta")
attributes(datextra$serno) <- NULL
datextra <- subset(datextra, select = -sex)
sum(names(dat) %in% names(datextra))
# since columns are in the same order, add the extra 7 to the main data
dat <- rbind(dat, datextra)
# remove datextra
rm(datextra)
# join information on allocation group, sex and age with anthropometry
dat <- inner_join(x = dat, y = datList)
str(dat)
datList$serno[!(datList$serno %in% dat$serno)]
# 4 serial number still missing.  
# 981 1536 1939 2587
rm(datList)

# take only blood pressure variables
datBP <- dat %>% dplyr::select(serno, mcdob, allocation, DateTest, BloodTestDate, age, sex, 
                               mstdhgtc, msbpuarm1c, msbpuarm2c, 
                               msbpuarm3c, mdbpuarm1c, mdbpuarm2c, mdbpuarm3c) 

# remove unused dataset
rm(dat)

# Harshad sent the available data by email on 09 June 2017. Extra data are in document Anthro.xlsx.
# Add data manually
more.rows <- data.frame(serno = c(981, 1536, 1939, 2587), stringsAsFactors=F)
datBP[(nrow(datBP) + 1):(nrow(datBP) + nrow(more.rows)), names(more.rows)] <- more.rows
rm(more.rows)
tail(datBP, 4)
# transform dat to data.table before adding extra information
datBP <- data.table(datBP)
# serno 981 and 2587 does not have bloop pressure measured
# serno 1536
datBP[serno == 1536, c("mstdhgtc", "msbpuarm1c", "msbpuarm2c", "msbpuarm3c", "mdbpuarm1c",
                    "mdbpuarm2c","mdbpuarm3c")  := 
      list(113.6, 85, 94, 91, 62, 53, 61), ]
# serno 1939
datBP[serno == 1939, c("mstdhgtc", "msbpuarm1c", "msbpuarm2c", "msbpuarm3c", "mdbpuarm1c",
                     "mdbpuarm2c","mdbpuarm3c")  := 
      list(111.2, 95,	93,	90,	63,	51, 55), ]

# sernos 981, 1536, 1939, 2587,have no sex recorded
# use the SARAS Masterlist v9 sent to CCMB to look for sex and add to the dataset.
datBP[serno == 981, sex := "Boy", ]
datBP[serno == 1536, sex := "Boy", ]
datBP[serno == 1939, sex := "Boy", ]
datBP[serno == 2587, sex := "Boy", ]
# sernos 981 1536 1939 2587 do not have any information od date of birth, allocation group, date of blood 
# test and date of first testing. Add the information using the SARAS Masterlist version 9
datBP[serno == 981, c("mcdob", "allocation", "DateTest", "BloodTestDate") := 
      list(as.Date("2007-11-28"), "B", as.Date("2014-04-23"), as.Date("2014-09-13")), ]
datBP[serno == 1536, c("mcdob", "allocation", "DateTest", "BloodTestDate") := 
      list(as.Date("2009-02-18"), "A", as.Date("2016-01-20"), as.Date("2016-01-21")), ]
datBP[serno == 1939, c("mcdob", "allocation", "DateTest", "BloodTestDate") := 
      list(as.Date("2009-01-07"), "A", as.Date("2016-02-17"), as.Date("2016-02-18")), ]
datBP[serno == 2587, c("mcdob", "allocation", "DateTest", "BloodTestDate") := 
      list(as.Date("2011-06-08"), "A", as.Date("2016-12-28"), as.Date("2016-12-29")), ]

# add a column for comments
datBP$Comments <- ""
# serno 2587 has Down syndrome. Add to the comment column
datBP[serno == 2587, "Comments" := "Child has Down Syndrome",]

# create the age variable and check children ange is between 5 and 7
datBP$age <- as.numeric(with(datBP, (BloodTestDate - mcdob)/365.25))
summary(datBP$age)
datBP[age < 5, ,]
ggplot(datBP, aes(x = age)) + geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 150)) + ggtitle("SARAS Kids: Age (years)") +
  xlab("Age (years)")

# check if sex agrees with the one sent to CCMB (v9 of the masterlist)
table(datBP$sex)

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
with(datBP, summary(msbpuarm1c))
with(datBP, sd(msbpuarm1c, na.rm = T))
with(datBP, summary(msbpuarm2c))
with(datBP, sd(msbpuarm2c, na.rm = T))
with(datBP, summary(msbpuarm3c))
with(datBP, sd(msbpuarm3c, na.rm = T))
# compare the distribution of SBP 1, 2 and 3
par(mfrow = c(3,1))
with(datBP, hist(msbpuarm1c))
with(datBP, hist(msbpuarm2c))
with(datBP, hist(msbpuarm3c))

# differences in the 3 measures of BP are possible. Only check unplausible values
datBP[datBP$msbpuarm1c > 130,]
datBP[datBP$msbpuarm1c < 60,]
datBP[datBP$msbpuarm2c > 130,]
datBP[datBP$msbpuarm2c < 60,]
datBP[datBP$msbpuarm3c > 130,]
datBP[datBP$msbpuarm3c < 60,]

# entries were checked with Harshad. On June 27, 2017 Harshad sent the following corrections
# serno 3099
datBP[serno == 3099, c("msbpuarm1c", "msbpuarm2c", "msbpuarm3c"):= list(94, 95, 99), ]
# serno 1095
datBP[serno == 1095, c("msbpuarm1c", "msbpuarm2c", "msbpuarm3c"):= list(101, 94, 93), ]
# serno 4270
datBP[serno == 4270, c("msbpuarm1c", "msbpuarm2c", "msbpuarm3c"):= list(104, 109, 106), ]

# compute average SBP
datBP$cavgsbp <- with(datBP, (msbpuarm1c + msbpuarm2c + msbpuarm3c)/3)
summary(datBP$cavgsbp)
sd(datBP$cavgsbp, na.rm = T)
par(mfrow = c(1,1))
hist(datBP$cavgsbp)
attributes(datBP$cavgsbp)$label <- "saras - average SBP (mmHg)"
# 4 children do not ha SBP measured

# check for outliers
cutoffU <- median(datBP$cavgsbp, na.rm = T) + 5*mad(datBP$cavgsbp, na.rm = T)
cutoffL <- median(datBP$cavgsbp, na.rm = T) - 5*mad(datBP$cavgsbp, na.rm = T)
sum(datBP$cavgsbp > cutoffU, na.rm = T)
sum(datBP$cavgsbp < cutoffL, na.rm = T)
# no outliers

# plot
ggplot(datBP, aes(x = cavgsbp)) + geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 75)) + 
  ggtitle("SARAS Kids: Systolic Blood Pressure (mmHg)") +
  xlab("SBP (mmHg)")


#### Diastolic Blood Pressure ####
with(datBP, summary(mdbpuarm1c))
with(datBP, sd(mdbpuarm1c, na.rm = T))
with(datBP, summary(mdbpuarm2c))
with(datBP, sd(mdbpuarm2c, na.rm = T))
with(datBP, summary(mdbpuarm3c))
with(datBP, sd(mdbpuarm3c, na.rm = T))
# compare the distribution of SBP 1, 2 and 3
par(mfrow = c(3,1))
with(datBP, hist(mdbpuarm1c))
with(datBP, hist(mdbpuarm2c))
with(datBP, hist(mdbpuarm3c))

# differences in the 3 measures of BP are possible. Only check unplausible values
datBP[datBP$mdbpuarm1c > 100,]
datBP[datBP$mdbpuarm1c < 30,]
datBP[datBP$mdbpuarm2c > 100,]
datBP[datBP$mdbpuarm2c < 30,]
datBP[datBP$mdbpuarm3c > 100,]
datBP[datBP$mdbpuarm3c < 30,]

# entries were checked with Harshad. On June 27, 2017 Harshad sent the following corrections
# serno 4270
datBP[serno == 4270, c("mdbpuarm1c", "mdbpuarm2c", "mdbpuarm3c")  := list(66, 63, 66), ]
# serno 5019
datBP[serno == 5019, c("mdbpuarm1c", "mdbpuarm2c", "mdbpuarm3c")  := list(76, 70, 68), ]


# compute average DBP
datBP$cavgdbp <- with(datBP, (mdbpuarm1c + mdbpuarm2c + mdbpuarm3c)/3)
summary(datBP$cavgdbp)
sd(datBP$cavgdbp, na.rm = T)
par(mfrow = c(1,1))
hist(datBP$cavgdbp)
datBP[datBP$cavgdbp < 40,]
attributes(datBP$cavgdbp)$label <- "saras - average DBP (mmHg)"
# 4 children do not ha DBP measured

# check for outliers
cutoffU <- median(datBP$cavgdbp, na.rm = T) + 5*mad(datBP$cavgdbp, na.rm = T)
cutoffL <- median(datBP$cavgdbp, na.rm = T) - 5*mad(datBP$cavgdbp, na.rm = T)
sum(datBP$cavgdbp > cutoffU, na.rm = T)
sum(datBP$cavgdbp < cutoffL, na.rm = T)
# 1 outlier based on the criteria outlined in the analysis plan
datBP[datBP$cavgdbp > cutoffU,,]

# create a flag variable for DBP
datBP$flagcavgdbp <- ifelse(datBP$cavgdbp > cutoffU, 1, 0)
datBP$flagcavgdbp <- factor(datBP$flagcavgdbp, labels = c("No Outlier", "Outlier"))
table(datBP$flagcavgdbp)
#rm(cutoffU); rm(cutoffL)
attributes(datBP$flagcavgdbp)$label <- "saras - DBP outlier?"

# summary and plot (no outlier removed)
summary(datBP$cavgdbp)
ggplot(datBP, aes(x = cavgdbp)) + geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 120)) + 
  ggtitle("SARAS Kids: Diastolic Blood Pressure (mmHg)") +
  xlab("DBP (mmHg)") + geom_vline(xintercept = cutoffU, col = "red")
# excluding the outliers
datBP %>% filter(flagcavgdbp == "No Outlier") %>% 
  summarise(n = n(), Min = min(cavgdbp), Mean = mean(cavgdbp), SD = 
              sd(cavgdbp), Median = median(cavgdbp), Q1 = 
              quantile(cavgdbp, probs = 0.25), Q3 = 
              quantile(cavgdbp, probs = 0.75), Max = max(cavgdbp))
ggplot(subset(datBP, flagcavgdbp == "No Outlier"), aes(x = cavgdbp)) + 
  geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) + 
  ggtitle("SARAS Kids (No Outliers): Diastolic Blood Pressure (mmHg)") +
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

# zhfa was already computed (stored in SARASKidsGrowthBodyCompv2_20170713.csv) import the data and take the variable we need
datHFA <- read.csv("./Cleaned Dataset/SARASKidsGrowthBodyCompv2_20170713.csv", sep = "")
datHFA <- datHFA %>% dplyr::select(serno, chfa)
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
  if(datBP$sex[i] == "Boy" & !is.na(datBP$chfa[i])){
    meanSBP[i] <- meanBP(age = datBP$age[i], height = datBP$chfa[i], 
                         alpha = 102.19768, b1 = 1.82416, b2 = 0.12776, b3 = 0.12776, b4 = -0.00135, 
                         g1 = 2.73157, g2 = -0.19618, g3 = -0.04659, g4 = 0.00947)}
  if(datBP$sex[i] == "Girl"  & !is.na(datBP$chfa[i])){
    meanSBP[i] <- meanBP(age = datBP$age[i], height = datBP$chfa[i], 
                         alpha = 102.01027, b1 = 1.94397, b2 = 0.00598, b3 = -0.000789, b4 = -0.00059, 
                         g1 = 2.03526, g2 = 0.02534, g3 = -0.01884, g4 = 0.00121)}
}
# add to the dataset
datBP$meanSBP <- meanSBP

# Diastolic Blood Pressure 
meanDBP <- rep(NA, nrow(datBP))
for(i in 1:nrow(datBP)){
  if(datBP$sex[i] == "Boy" & !is.na(datBP$chfa[i])){
    meanDBP[i] <- meanBP(age = datBP$age[i], height = datBP$chfa[i], 
                         alpha = 61.01217, b1 = 0.68314, b2 = -0.09835, b3 = 0.01711, b4 = 0.00045, 
                         g1 = 1.46993, g2 = -0.07849, g3 = -0.03144, g4 = 0.00967)}
  if(datBP$sex[i] == "Girl"  & !is.na(datBP$chfa[i])){
    meanDBP[i] <- meanBP(age = datBP$age[i], height = datBP$chfa[i], 
                         alpha = 60.50510, b1 = 1.01301, b2 = 0.01157, b3 = 0.00424, b4 = -0.00137, 
                         g1 = 1.16641, g2 = 0.12795, g3 = -0.03869, g4 = -0.00079)}
}
# add to the dataset
datBP$meaDBP <- meanDBP


# Convert the observed BP to a Z-score (Zbp)

# add a sigma in the dataset based on the child sex
for(i in 1:nrow(datBP)){
  if(datBP$sex[i] == "Boy"){
    datBP$sigmaSBP[i] <- 10.7128
    datBP$sigmaDBP[i] <- 11.6032}
  if(datBP$sex[i] == "Girl"){
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


###############################################################################
##### Cardio Metabolic Data ###################################################
###############################################################################

# all the cardio metabolic data (with the exception of 120min glucose) were processed in the Pune lab.
# data were sent in 2 batches and put together by Patsy in saraskem01.dta (23/01/2018). Load the data
# total cholesterol is not available for this cohort
datCardio <- read_dta("saraskem01.dta")
# merge with child blood data to check if all are in the PP
datList <- read_dta("SARASKidsBloodCountPP.dta")
datList <- datList %>% dplyr::select(serno, allocation, DateTest, BloodTestDate, age, sex)
datList$sex <- factor(datList$sex, labels = c("Boy", "Girl"))
datList$allocation <- factor(datList$allocation, labels = c("B", "A"))
# check if all sernos in datList are in dat
datCardio <- inner_join(datCardio, datList) 
# out of 713 done, 700 are in the PP analysis
datList$serno[!(datList$serno %in% datCardio$serno)] # missing sernos
rm(datList)

# The glucose values were measured in Mumbai and Pune. The Mumbai values are more credible than 
# the Pune ones (some of the samples analysed in Pune were likely to have degraded)
# Use the Mumbai glucose values throughout (0, 30 and 120 minutes).

datglu <- read_dta("sarasparag01.dta")
# keep only serno and pglu120
datglu <- datglu %>% dplyr::select(serno, pglu0, pglu30, pglu120)
# merge with the data from the Pune lab
datCardio <- left_join(datCardio, datglu)
# remove datglu
rm(datglu)
# drop the Pune glucose value from the dataset
datCardio <- subset(datCardio, select = -c(glu0, glu30))
# rename pglu0, pglu30 and pglu120 
names(datCardio)[names(datCardio) == "pglu0"] <- "glu0"
names(datCardio)[names(datCardio) == "pglu30"] <- "glu30"
names(datCardio)[names(datCardio) == "pglu120"] <- "glu120"

# convert in right units: glucose, LDL-cholesterol, HDL-cholesterol and triglycerides are mg/dl, 
# and would need to be converted into mmol/l. 
datCardio$glu0 <- 0.0555*datCardio$glu0 
datCardio$glu30 <- 0.0555*datCardio$glu30
datCardio$glu120 <- 0.0555*datCardio$glu120
# other units of conversion can be found on the following website:
# https://www.ncbi.nlm.nih.gov/books/NBK33478/
# cholesterol: To get from mg/dL mmol/L multiply by 0.02586
datCardio$hdl <- 0.02586*datCardio$hdl
datCardio$ldl <- 0.02586*datCardio$ldl
# tryglicerides: to get from mg/dL to mmol/l units, multiply by 0.01129
datCardio$tg <- 0.01129*datCardio$tg
# Insulin: as per Caroline email (18/12/2017), the values need to be multiplied by 6
datCardio$ins0 <- 6*datCardio$ins0
datCardio$ins30 <- 6*datCardio$ins30


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

#### Fasting Glucose ####
summary(datCardio$glu0)
sd(datCardio$glu0, na.rm = T)
hist(datCardio$glu0, main = "Fasting Glucose", xlab = "Fasting Glucose")
# 11 missing
# check for outliers
cutoffU <- median(datCardio$glu0, na.rm = T) + 5*mad(datCardio$glu0, na.rm = T)
cutoffL <- median(datCardio$glu0, na.rm = T) - 5*mad(datCardio$glu0, na.rm = T)
sum(datCardio$glu0 > cutoffU, na.rm = T)
sum(datCardio$glu0 < cutoffL, na.rm = T)
# 3 outliers based on criteria established in the analysis plan

# create a flag variable for fasting glucose
datCardio$flagglu0 <- ifelse(datCardio$glu0 > cutoffU, 1, 0)
datCardio$flagglu0 <- factor(datCardio$flagglu0, labels = c("No Outlier", "Outlier"))
table(datCardio$flagglu0)
#rm(cutoffU); rm(cutoffL)
attributes(datCardio$flagglu0)$label <- "saras - fasting glucose outleir?"

# summary and plot (no outlier removed)
summary(datCardio$glu0)
ggplot(datCardio, aes(x = glu0)) + geom_histogram(colour="black", fill="white", binwidth = 0.26) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 200)) + ggtitle("SARAS Kids: Fasting Glucose (mmol/l)") +
  xlab("Fasting Glucose (mmol/l)") + geom_vline(xintercept = cutoffU, col = "red")
# excluding the outliers
datCardio %>% filter(flagglu0 == "No Outlier") %>% 
  summarise(n = n(), Min = min(glu0), Mean = mean(glu0), SD = 
              sd(glu0), Median = median(glu0), Q1 = 
              quantile(glu0, probs = 0.25), Q3 = 
              quantile(glu0, probs = 0.75), Max = max(glu0))
ggplot(subset(datCardio, flagglu0 == "No Outlier"), aes(x = glu0)) + 
  geom_histogram(colour="black", fill="white", binwidth = 0.26) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 200)) + 
  ggtitle("SARAS Kids (No Outliers): Fasting Glucose (mmol/l)") +
  xlab("Fasting Glucose (mmol/l)")


#### 30 min Glucose ####
summary(datCardio$glu30)
sd(datCardio$glu30, na.rm = T)
hist(datCardio$glu30, main = "30min Glucose", xlab = "30min Glucose")
# check for outliers
cutoffU <- median(datCardio$glu30, na.rm = T) + 5*mad(datCardio$glu30, na.rm = T)
cutoffL <- median(datCardio$glu30, na.rm = T) - 5*mad(datCardio$glu30, na.rm = T)
sum(datCardio$glu30 > cutoffU, na.rm = T)
sum(datCardio$glu30 < cutoffL, na.rm = T)
# 0 outliers based on criteria established in the analysis plan

# summary and plot (no outlier removed)
summary(datCardio$glu30)
ggplot(datCardio, aes(x = glu30)) + geom_histogram(colour="black", fill="white", binwidth = 0.8) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 150)) + ggtitle("SARAS Kids: 30min Glucose (mmol/l)") +
  xlab("30min Glucose (mmol/l)")


#### 120 min Glucose ####
summary(datCardio$glu120)
sd(datCardio$glu120, na.rm = T)
hist(datCardio$glu120, main = "120min Glucose", xlab = "120min Glucose")
# check for outliers
cutoffU <- median(datCardio$glu120, na.rm = T) + 5*mad(datCardio$glu120, na.rm = T)
cutoffL <- median(datCardio$glu120, na.rm = T) - 5*mad(datCardio$glu120, na.rm = T)
sum(datCardio$glu120 > cutoffU, na.rm = T)
sum(datCardio$glu120 < cutoffL, na.rm = T)
# 0 outliers based on criteria established in the analysis plan

# summary and plot (no outlier removed)
summary(datCardio$glu120)
ggplot(datCardio, aes(x = glu120)) + geom_histogram(colour="black", fill="white", binwidth = 0.45) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 150)) + ggtitle("SARAS Kids: 120min Glucose (mmol/l)") +
  xlab("120min Glucose (mmol/l)")


#### HDL Cholesterol ####
summary(datCardio$hdl)
sd(datCardio$hdl, na.rm = T)
hist(datCardio$hdl, main = "HDL cholesterol", xlab = "HDL cholesterol")
# check for outliers
cutoffU <- median(datCardio$hdl, na.rm = T) + 5*mad(datCardio$hdl, na.rm = T)
cutoffL <- median(datCardio$hdl, na.rm = T) - 5*mad(datCardio$hdl, na.rm = T)
sum(datCardio$hdl > cutoffU, na.rm = T)
sum(datCardio$hdl < cutoffL, na.rm = T)
# 0 outliers based on criteria established in the analysis plan

# summary and plot (no outlier removed)
summary(datCardio$hdl)
ggplot(datCardio, aes(x = hdl)) + geom_histogram(colour="black", fill="white", binwidth = 0.10) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 150)) + ggtitle("SARAS Kids: HDL Cholesterol (mmol/l)") +
  xlab("HDL Cholesterol (mmol/l)")



#### LDL Cholesterol ####
summary(datCardio$ldl)
sd(datCardio$ldl, na.rm = T)
hist(datCardio$ldl, main = "LDL cholesterol", xlab = "LDL cholesterol")
# check for outliers
cutoffU <- median(datCardio$ldl, na.rm = T) + 5*mad(datCardio$ldl, na.rm = T)
cutoffL <- median(datCardio$ldl, na.rm = T) - 5*mad(datCardio$ldl, na.rm = T)
sum(datCardio$ldl > cutoffU, na.rm = T)
sum(datCardio$ldl < cutoffL, na.rm = T)
# 1 outliers based on criteria established in the analysis plan
# after talking with Caroline, delete the LDL cholesterol, as it might be related to some genetic
# problem
datCardio[datCardio$ldl > cutoffU & !is.na(datCardio$ldl), "ldl"] <- NA
# serno 1545

# summary and plot
summary(datCardio$ldl)
ggplot(datCardio, aes(x = ldl)) + geom_histogram(colour="black", fill="white", binwidth = 0.25) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 150)) + ggtitle("SARAS Kids: LDL Cholesterol (mmol/l)") +
  xlab("LDL Cholesterol (mmol/l)") + geom_vline(xintercept = cutoffU, col = "red")


#### Fasting Insulin ####
summary(datCardio$ins0)
sd(datCardio$ins0, na.rm = T)
hist(datCardio$ins0, main = "Fasting Insulin", xlab = "Fasting Insulin")
# check for outliers
cutoffU <- median(datCardio$ins0, na.rm = T) + 5*mad(datCardio$ins0, na.rm = T)
cutoffL <- median(datCardio$ins0, na.rm = T) - 5*mad(datCardio$ins0, na.rm = T)
sum(datCardio$ins0 > cutoffU, na.rm = T)
sum(datCardio$ins0 < cutoffL, na.rm = T)
# 31 outliers based on criteria established in the analysis plan
# line with outliers
temp <- datCardio[!is.na(datCardio$ins0) & datCardio$ins0 > cutoffU, c("glu0", "flagglu0")]
rm(temp)
# 2 of those 31 were identified to have outlier fasting glucose

# create a flag variable for fastin insulin
datCardio$flagins0 <- ifelse(datCardio$ins0 > cutoffU, 1, 0)
datCardio$flagins0 <- factor(datCardio$flagins0, labels = c("No Outlier", "Outlier"))
table(datCardio$flagins0)
#rm(cutoffU); rm(cutoffL)
attributes(datCardio$flagins0)$label <- "saras - fasting insulin outlier?"

# summary and plot (no outlier removed)
summary(datCardio$ins0)
ggplot(datCardio, aes(x = ins0)) + geom_histogram(colour="black", fill="white", binwidth = 25) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 350)) + ggtitle("SARAS Kids: Fasting Insulin (pmol/l)") +
  xlab("Fasting Insulin (pmol/l)") + geom_vline(xintercept = cutoffU, col = "red")
# excluding the outliers
datCardio %>% filter(flagins0 == "No Outlier") %>% 
  summarise(n = n(), Min = min(ins0), Mean = mean(ins0), SD = 
              sd(ins0), Median = median(ins0), Q1 = 
              quantile(ins0, probs = 0.25), Q3 = 
              quantile(ins0, probs = 0.75), Max = max(ins0))
ggplot(subset(datCardio, flagins0 == "No Outlier"), aes(x = ins0)) + 
  geom_histogram(colour="black", fill="white", binwidth = 5) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 175)) + 
  ggtitle("SARAS Kids (No Outliers): Fasting Insulin (pmol/l)") +
  xlab("Fasting Insulin (pmol/l)")


#### 30 min Insulin ####
summary(datCardio$ins30)
sd(datCardio$ins30, na.rm = T)
hist(datCardio$ins30, main = "30min Insulin", xlab = "30min Insulin")
# check for outliers
cutoffU <- median(datCardio$ins30, na.rm = T) + 5*mad(datCardio$ins30, na.rm = T)
cutoffL <- median(datCardio$ins30, na.rm = T) - 5*mad(datCardio$ins30, na.rm = T)
sum(datCardio$ins30 > cutoffU, na.rm = T)
sum(datCardio$ins30 < cutoffL, na.rm = T)
# 7 outliers based on criteria established in the analysis plan
datCardio[!is.na(datCardio$ins30) & datCardio$ins30 > cutoffU, c("ins0", "flagins0")]
# only one had outlier fasting insulin, none had outlier fasting glucose

# create a flag variable for 30min insulin
datCardio$flagins30 <- ifelse(datCardio$ins30 > cutoffU, 1, 0)
datCardio$flagins30 <- factor(datCardio$flagins30, labels = c("No Outlier", "Outlier"))
table(datCardio$flagins30)
#rm(cutoffU); rm(cutoffL)
attributes(datCardio$flagins30)$label <- "saras - 30min insulin outlier?"

# summary and plot (no outlier removed)
summary(datCardio$ins30)
ggplot(datCardio, aes(x = ins30)) + geom_histogram(colour="black", fill="white", binwidth = 50) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 150)) + ggtitle("SARAS Kids: 30min Insulin (pmol/l)") +
  xlab("30min Insulin (pmol/l)") + geom_vline(xintercept = cutoffU, col = "red")
# excluding the outliers
datCardio %>% filter(flagins30 == "No Outlier") %>% 
  summarise(n = n(), Min = min(ins30), Mean = mean(ins30), SD = 
              sd(ins30), Median = median(ins30), Q1 = 
              quantile(ins30, probs = 0.25), Q3 = 
              quantile(ins30, probs = 0.75), Max = max(ins30))
ggplot(subset(datCardio, flagins30 == "No Outlier"), aes(x = ins30)) + 
  geom_histogram(colour="black", fill="white", binwidth = 50) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 175)) + 
  ggtitle("SARAS Kids (No Outliers): 30min Insulin (pmol/l)") +
  xlab("30min Insulin (pmol/l)")

# caroline asked for the oulier high fasting insulin
#datOutlier <- datCardio %>% filter(flagins0 == "Outlier") %>%
#  dplyr::select(serno, glu0, glu30, glu120, ins0, ins30)
# save as csv
#write.table(datOutlier, file = "IsulinOutliers.csv", row.names = FALSE)

#### Triglycerides ####
summary(datCardio$tg)
sd(datCardio$tg, na.rm = T)
hist(datCardio$tg, main = "Triglycerides", xlab = "Triglycerides")
# check for outliers
cutoffU <- median(datCardio$tg, na.rm = T) + 5*mad(datCardio$tg, na.rm = T)
cutoffL <- median(datCardio$tg, na.rm = T) - 5*mad(datCardio$tg, na.rm = T)
sum(datCardio$tg > cutoffU, na.rm = T)
sum(datCardio$tg < cutoffL, na.rm = T)
# 7 outleirs identified based on the analysis plan
datCardio[!is.na(datCardio$tg) & datCardio$tg > cutoffU, c("ins0", "flagins0")]
datCardio[!is.na(datCardio$tg) & datCardio$tg > cutoffU, c("glu0", "flagglu0")]
# none had outleirs fasting insulin and glucose

# create a flag variable for triglycerides
datCardio$flagtg <- ifelse(datCardio$tg > cutoffU, 1, 0)
datCardio$flagtg <- factor(datCardio$flagtg, labels = c("No Outlier", "Outlier"))
table(datCardio$flagtg)
#rm(cutoffU); rm(cutoffL)
attributes(datCardio$flagtg)$label <- "saras - triglycerides outlier?"

# summary and plot (no outlier removed)
summary(datCardio$tg)
ggplot(datCardio, aes(x = tg)) + geom_histogram(colour="black", fill="white", binwidth = 0.13) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 150)) + ggtitle("SARAS Kids: Triglycerides (mmol/l)") +
  xlab("Triglycerides (mmol/l)") + geom_vline(xintercept = cutoffU, col = "red")
# excluding the outliers
datCardio %>% filter(flagtg == "No Outlier") %>% 
  summarise(n = n(), Min = min(tg), Mean = mean(tg), SD = 
              sd(tg), Median = median(tg), Q1 = 
              quantile(tg, probs = 0.25), Q3 = 
              quantile(tg, probs = 0.75), Max = max(tg))
ggplot(subset(datCardio, flagtg == "No Outlier"), aes(x = tg)) + 
  geom_histogram(colour="black", fill="white", binwidth = 0.2) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 250)) + 
  ggtitle("SARAS Kids (No Outliers): Triglycerides (mmol/l)") +
  xlab("Triglycerides (mmol/l)")


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
# fileForHOMA <- datCardio %>% dplyr::select(serno, glu0, ins0)
# fasting insulin needs to be at least 20 for the calculator to work. Transform all the values <20 to 20
# fileForHOMA$ins0[fileForHOMA$ins0 < 20] <- 20
# write.table(fileForHOMA, "SARASfileForHOMA.csv", sep = ",", row.names = F)
# rm(fileForHOMA)
# an excel spreadsheet named childHOMA.xlxs was created with the HOMA values
# load the HOMA values and add them to the dataset

# start again from here
datHOMA <-  read_excel("ChildHOMA.xlsx", sheet = "India", na = "NA")
datHOMA <- datHOMA %>% dplyr::select(serno, chomas, chomair)
# merge datHOMA with datCardio
datCardio <- inner_join(datCardio, datHOMA)
rm(datHOMA)
datCardio$chomair <- as.numeric(datCardio$chomair)
datCardio$chomas <- as.numeric(datCardio$chomas)

# summary of HOMA-IR
summary(datCardio$chomair)
sd(datCardio$chomair, na.rm = T)
hist(datCardio$chomair)
datCardio %>% filter(chomair > 3) %>% select(glu0, ins0, flagins0, flagglu0, chomair)
# those with HOMA-IR greater than 3 have outlier values for fasting insulin. 

# check for outliers
cutoffU <- median(datCardio$chomair, na.rm = T) + 5*mad(datCardio$chomair, na.rm = T)
cutoffL <- median(datCardio$chomair, na.rm = T) - 5*mad(datCardio$chomair, na.rm = T)
sum(datCardio$chomair > cutoffU, na.rm = T)
sum(datCardio$chomair < cutoffL, na.rm = T)
# based on the analysis plan there are 188 outliers. All with highe HOMA-IR values

# create a flag variable for triglycerides
datCardio$flagchomair <- ifelse(datCardio$chomair > cutoffU, 1, 0)
datCardio$flagchomair <- factor(datCardio$flagchomair, labels = c("No Outlier", "Outlier"))
table(datCardio$flagchomair)
#rm(cutoffU); rm(cutoffL)
attributes(datCardio$flagchomair)$label <- "saras - HOMA-IR outlier?"

# summary and plot (no outlier removed)
summary(datCardio$chomair)
ggplot(datCardio, aes(x = chomair)) + geom_histogram(colour="black", fill="white", binwidth = 0.35) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 500)) + ggtitle("SARAS Kids: HOMA-IR") +
  xlab("HOMA-IR") + geom_vline(xintercept = cutoffU, col = "red")
# excluding the outliers
datCardio %>% filter(flagchomair == "No Outlier") %>% 
  summarise(n = n(), Min = min(chomair), Mean = mean(chomair), SD = 
              sd(chomair), Median = median(chomair), Q1 = 
              quantile(chomair, probs = 0.25), Q3 = 
              quantile(chomair, probs = 0.75), Max = max(chomair))
ggplot(subset(datCardio, flagchomair == "No Outlier"), aes(x = chomair)) + 
  geom_histogram(colour="black", fill="white", binwidth = 0.02) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 250)) + 
  ggtitle("SARAS Kids (No Outliers): HOMA-IR") +
  xlab("HOMA-IR")

# HOMA-IR is the reciprocal of HOMA-S (100/%S). Check:
plot(datCardio$chomair, 1/datCardio$chomas)
cor(datCardio$chomair, 1/datCardio$chomas, use = "complete.obs")

#### Compute Insulinogenic Index ####
# computed as (Insulin30 - Insulin0)/(Glucose30-Glucose0)
# There are still quite a lot of children whose 30 minute glucose is lower than their fasting value. 
# Clive suggests using the following to calculate the insulogenic index, to deal with this problem:
# (log [insulin30] – log [insulin0]) / (log [glucose30] – log [glucose0])

datCardio$cinsindex <- with(datCardio, (log(ins30) - log(ins0))/(log(glu30) - log(glu0)))
summary(datCardio$cinsindex)
sum(!is.na(datCardio$cinsindex)) 
# 662 have insulinogenic index
sd(datCardio$cinsindex, na.rm = T)
# how many have insulinogenic index negative?
datCardio[datCardio$cinsindex <0 & !is.na(datCardio$cinsindex), c("serno", "cinsindex")]
# 70 out of 662
(70/662)*100
hist(datCardio$cinsindex)

# check for outliers
cutoffU <- median(datCardio$cinsindex, na.rm = T) + 5*mad(datCardio$cinsindex, na.rm = T)
cutoffL <- median(datCardio$cinsindex, na.rm = T) - 5*mad(datCardio$cinsindex, na.rm = T)
sum(datCardio$cinsindex > cutoffU, na.rm = T)
sum(datCardio$cinsindex < cutoffL, na.rm = T)
# based on the analysis plan there are 56 outliers

# create a flag variable for insulinogenic index
datCardio$flagcinsindex <- ifelse(datCardio$cinsindex > cutoffU | datCardio$cinsindex < cutoffL, 1, 0)
datCardio$flagcinsindex <- factor(datCardio$flagcinsindex, labels = c("No Outlier", "Outlier"))
table(datCardio$flagcinsindex)
#rm(cutoffU); rm(cutoffL)
attributes(datCardio$flagcinsindex)$label <- "saras - Insulinogenic index outlier?"

# summary and plot (no outlier removed)
summary(datCardio$cinsindex)
ggplot(datCardio, aes(x = cinsindex)) + 
  geom_histogram(colour="black", fill="white", binwidth = 30) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 700)) + 
  ggtitle("SARAS Kids: Insulinogenic Index") +
  xlab("Insulinogenic Index") + geom_vline(xintercept = cutoffU, col = "red") + 
  geom_vline(xintercept = cutoffL, col = "red")
# excluding the outliers
datCardio %>% filter(flagcinsindex == "No Outlier") %>% 
  summarise(n = n(), Min = min(cinsindex), Mean = mean(cinsindex), SD = 
              sd(cinsindex), Median = median(cinsindex), Q1 = 
              quantile(cinsindex, probs = 0.25), Q3 = 
              quantile(cinsindex, probs = 0.75), Max = max(cinsindex))
ggplot(subset(datCardio, flagcinsindex == "No Outlier"), aes(x = cinsindex)) + 
  geom_histogram(colour="black", fill="white", binwidth = 2) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 200)) + 
  ggtitle("SARAS Kids (No Outliers): Insulinogenic Index") +
  xlab("Insulinogenic Index")


#### Compute Disposition Index ####

# disposition index is calculated as Insulogenic index/HOMA-IR
datCardio$cdispindex <- with(datCardio, cinsindex/chomair)
summary(datCardio$cdispindex)
sd(datCardio$cdispindex, na.rm = T)
sum(!is.na(datCardio$cdispindex)) # 662

# check for outliers
cutoffU <- median(datCardio$cdispindex, na.rm = T) + 5*mad(datCardio$cdispindex, na.rm = T)
cutoffL <- median(datCardio$cdispindex, na.rm = T) - 5*mad(datCardio$cdispindex, na.rm = T)
sum(datCardio$cdispindex > cutoffU, na.rm = T)
sum(datCardio$cdispindex < cutoffL, na.rm = T)
# based on the analysis plan there are 46 outliers

# create a flag variable for triglycerides
datCardio$flagcdispindex <- ifelse(datCardio$cdispindex > cutoffU | datCardio$cdispindex < cutoffL, 1, 0)
datCardio$flagcdispindex <- factor(datCardio$flagcdispindex, labels = c("No Outlier", "Outlier"))
table(datCardio$flagcdispindex)
rm(cutoffU); rm(cutoffL)
attributes(datCardio$flagcdispindex)$label <- "saras - Disposition index outlier?"

# summary and plot (no outlier removed)
summary(datCardio$cdispindex)
hist(datCardio$cdispindex, main = "SARAD Kinds: Disposition Index", xlab = "Disposition Index")
# excluding the outliers
datCardio %>% filter(flagcdispindex == "No Outlier") %>% 
  summarise(n = n(), Min = min(cdispindex), Mean = mean(cdispindex), SD = 
              sd(cdispindex), Median = median(cdispindex), Q1 = 
              quantile(cdispindex, probs = 0.25), Q3 = 
              quantile(cdispindex, probs = 0.75), Max = max(cdispindex))
ggplot(subset(datCardio, flagcdispindex == "No Outlier"), aes(x = cdispindex)) + 
  geom_histogram(colour="black", fill="white", binwidth = 5) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 150)) + 
  ggtitle("SARAS Kids (No Outliers): Disposition Index") +
  xlab("Disposition Index")

# check insulinogenic index X HOMA-%S
plot(datCardio$cdispindex, datCardio$cinsindex*datCardio$chomas)
# same!

#### Compute Metabolic Syndrome ####

# Create a 0/1 variable where 1 is for children who are above the sex-specific within-population upper quartile 
# for android fat on DXA and SBP and triglycerides and HOMA-IR, and below the lower quartile for HDL-cholesterol.

# to do this we need to merge the blood pressure and the cardio metabolic data, and add the DXA android fat
datBP <- datBP %>% dplyr::select(serno, sex, allocation, age, cavgsbp, cavgdbp, flagcavgdbp, cSBPper99, 
                                 cDBPper99, Comments)
# delete sex, allocation and age from datCardio
datCardio <- datCardio[, !(names(datCardio) %in% c("sex", "age", "allocation"))]
dat <- left_join(datBP, datCardio)
rm(datBP, datCardio)
rm(meanDBP, meanSBP)
# add android fat from DXA
datDXA <- read_dta("./Cleaned Dataset/SARASKidsGrowthBodyCompv2_20170713.dta")
datDXA <- datDXA %>% dplyr::select(serno, candroidfat)
# add android fat to dataset
dat <- left_join(dat, datDXA)
rm(datDXA)

# sex-specific within-population higher than the upper quartile for android fat on DXA
DXAquart <- dat %>% group_by(sex) %>% summarise(Q3 = quantile(candroidfat, prob = 0.75, na.rm = T))
qDXA <- c()
for(i in 1:nrow(dat)){
  if(dat$sex[i] == "Boy" & !is.na(dat$candroidfat[i]) & dat$candroidfat[i] > DXAquart[1,"Q3"]){qDXA[i] <- 1
  } else if(dat$sex[i] == "Boy" & !is.na(dat$candroidfat[i]) & dat$candroidfat[i] <= DXAquart[1,"Q3"]){
       qDXA[i] <- 0
  } else if(dat$sex[i] == "Girl" & !is.na(dat$candroidfat[i]) & dat$candroidfat[i] > DXAquart[2,"Q3"]){
         qDXA[i] <- 1
  } else if(dat$sex[i] == "Girl" & !is.na(dat$candroidfat[i]) & dat$candroidfat[i] <= DXAquart[2,"Q3"]){
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
  if(dat$sex[i] == "Boy" & !is.na(dat$cavgsbp[i]) & dat$cavgsbp[i] > SBPquart[1,"Q3"]){qSBP[i] <- 1
  } else if(dat$sex[i] == "Boy" & !is.na(dat$cavgsbp[i]) & dat$cavgsbp[i] <= SBPquart[1,"Q3"]){
    qSBP[i] <- 0
  } else if(dat$sex[i] == "Girl" & !is.na(dat$cavgsbp[i]) & dat$cavgsbp[i] > SBPquart[2,"Q3"]){
    qSBP[i] <- 1
  } else if(dat$sex[i] == "Girl" & !is.na(dat$cavgsbp[i]) & dat$cavgsbp[i] <= SBPquart[2,"Q3"]){
    qSBP[i] <- 0
  } else{ qSBP[i] <- NA}
}
dat$qSBP <- qSBP
rm(qSBP)
rm(SBPquart)

# sex-specific within-population upper quartile for triglycerides
tgquart <- dat %>% group_by(sex) %>% summarise(Q3 = quantile(tg, prob = 0.75, na.rm = T))
qtg <- c()
for(i in 1:nrow(dat)){
  if(dat$sex[i] == "Boy" & !is.na(dat$tg[i]) & dat$tg[i] > tgquart[1,"Q3"]){qtg[i] <- 1
  } else if(dat$sex[i] == "Boy" & !is.na(dat$tg[i]) & dat$tg[i] <= tgquart[1,"Q3"]){
    qtg[i] <- 0
  } else if(dat$sex[i] == "Girl" & !is.na(dat$tg[i]) & dat$tg[i] > tgquart[2,"Q3"]){
    qtg[i] <- 1
  } else if(dat$sex[i] == "Girl" & !is.na(dat$tg[i]) & dat$tg[i] <= tgquart[2,"Q3"]){
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
  if(dat$sex[i] == "Boy" & !is.na(dat$chomair[i]) & dat$chomair[i] > HOMAquart[1,"Q3"]){
    qHOMA[i] <- 1
  } else if(dat$sex[i] == "Boy" & !is.na(dat$chomair[i]) & dat$chomair[i] <= HOMAquart[1,"Q3"]){
    qHOMA[i] <- 0
  } else if(dat$sex[i] == "Girl" & !is.na(dat$chomair[i]) & dat$chomair[i] > HOMAquart[2,"Q3"]){
    qHOMA[i] <- 1
  } else if(dat$sex[i] == "Girl" & !is.na(dat$chomair[i]) & dat$chomair[i] <= HOMAquart[2,"Q3"]){
    qHOMA[i] <- 0
  } else{qHOMA[i] <- NA}
}
dat$qHOMA <- qHOMA
rm(qHOMA)
rm(HOMAquart)

# sex-specific within-population below the lower quartile for HDL
HDLquart <- dat %>% group_by(sex) %>% summarise(Q1 = quantile(hdl, prob = 0.25, na.rm = T))
qHDL <- c()
for(i in 1:nrow(dat)){
  if(dat$sex[i] == "Boy" & !is.na(dat$hdl[i]) & dat$hdl[i] < HDLquart[1,"Q1"]){qHDL[i] <- 1
  } else if(dat$sex[i] == "Boy" & !is.na(dat$hdl[i]) & dat$hdl[i] >= HDLquart[1,"Q1"]){
    qHDL[i] <- 0
  } else if(dat$sex[i] == "Girl" & !is.na(dat$hdl[i]) & dat$hdl[i] < HDLquart[2,"Q1"]){
    qHDL[i] <- 1
  } else if(dat$sex[i] == "Girl" & !is.na(dat$hdl[i]) & dat$hdl[i] >= HDLquart[2,"Q1"]){
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


#######################################################################################################################
#### Save Dataset #####################################################################################################
#######################################################################################################################

# select only relevant variable (based on the analysis plan)
dat <- dat %>% dplyr::select(serno, sex, age, allocation, cavgsbp, cavgdbp, cSBPper99, cDBPper99, glu0, glu30, glu120,
                             tg, hdl, ldl, ins0, ins30, chomair, cinsindex, cdispindex, cmetsynd,
                             flagcavgdbp, flagglu0, flagins0, flagins30, flagtg, flagchomair,
                             flagcinsindex, flagcdispindex, Comments)  
names(dat) <- c("serno", "csex", "cage", "allocation", "cavgSBP", "cavgDBP", "cSBPper99", "cDBPper99", "cglu0",
                "cglu30", "cglu120", "ctg", "chdl", "cldl", "cins0", "cins30", "choma", "cinsindex", "cdispindex",
                "cmetsynd", "flagcavgDBP", "flagcglu0", "flagcins0", "flagcins30", "flagctg", "flagchoma",
                "flagcinsindex", "flagcdispindex", "comments")

# add attributes
attributes(dat$serno)$label <- "saras - subject ID"
attributes(dat$csex)$label <- "saras - sex"
attributes(dat$cage)$label <- "saras - age (years)"
attributes(dat$allocation)$label <- "saras - allocation group"
attributes(dat$cSBPper99)$label <- "saras - SBP greater than 99th percentile"
attributes(dat$cDBPper99)$label <- "saras - DBP greater than 99th percentile"
attributes(dat$choma)$label <- "saras - HOMA-IR"
attributes(dat$cmetsynd)$label <- "saras - metabolic syndrome"
attributes(dat$cglu0)$label <- "saras - fasting glucose (mmol/l)"
attributes(dat$cglu30)$label <- "saras - 30min glucose (mmol/l)"
attributes(dat$cglu120)$label <- "saras - 120 glucose (mmol/l)"
attributes(dat$ctg)$label <- "saras - triglycerides (mmol/l)"
attributes(dat$chdl)$label <- "saras - HDL cholesterol (mmol/l)"
attributes(dat$cldl)$label <- "saras - LDL cholesterol (mmol/l)"
attributes(dat$cins0)$label <- "saras - fasting insulin (pmol/l)"
attributes(dat$cins30)$label <- "saras - 30min insulin (pmol/l)"
attributes(dat$cinsindex)$label <- "saras - insulinogenic index"
attributes(dat$cdispindex)$label <- "saras - disposition index"
attributes(dat$comments)$label <- "saras - comments"


# save as csv as stata ans as spss file
write.table(dat, "./Cleaned Dataset/SARASKidsCardioMetv4_20180409.csv", row.names = FALSE)
write_dta(dat, "./Cleaned Dataset/SARASKidsCardioMetv4_20180409.dta")
write_sav(dat, "./Cleaned Dataset/SARASKidsCardioMetv4_20180409.sav")
