###############################################################
##### SARAS Kids: Phenotype Data Processing Birth Outcome #####
###############################################################

# 20/06/2017

# Dataset needed for the cleaning process:
# MMNPanalysis_birthoutcomes.dta
# SARASKidsBloodCountPP.dta
# MMNPanalysis_birthoutcomes.dta was created by Ella and used for the analysis reported in th
# the Potdar et al paper (2014)
# Only children in the per-protocol analysis are processed.

# load packages
library(haven)
library(dplyr)
library(ggplot2)
library(readxl)

# import data and check import is correct 
dat <- read_dta("MMNPanalysis_birthoutcomes.dta")
dim(dat)
# delect sex, as the categories are not updated
dat <- subset(dat, select = -sex)

# add information on blood counts, allocation group, sex and age. PP analysis only
datList <- read_dta("SARASKidsBloodCountPP.dta")
datList$sex <- as_factor(datList$sex)
datList <- datList %>% select(serno, AnalysisType, allocation, BloodTestDate, age, sex)
# check if all sernos in datList are in dat
datList$serno[!(datList$serno %in% dat$serno)] 
dat <- inner_join(x = dat, y = datList)
str(dat)
rm(datList)


###########################################################################################################
#### Univariate Analysis ##################################################################################
###########################################################################################################

# 1. Check each variable singularly
# 2. Check for outliers using the criteria established in the analysis plan (median +/- 5MAD)
# 5. Retain the outliers but create a variable with 0/1 (named flagvariablename) where 1 indicated
#    the outlier based on the analysis plan.


##########################################################################################################
#### Measured outcome variables ##########################################################################
##########################################################################################################

#### Birth weight ####
summary(dat$bwt1)
# 156 (22%) do not have any birth weight recorded
sd(dat$bwt1)

# check for outliers
cutoffU <- median(dat$bwt1, na.rm = T) + 5*mad(dat$bwt1, na.rm = T)
cutoffL <- median(dat$bwt1, na.rm = T) - 5*mad(dat$bwt1, na.rm = T)
sum(dat$bwt1 > cutoffU, na.rm = T)
sum(dat$bwt1 < cutoffL, na.rm = T)
# no outliers based on the criteria established in the Analysis Plan
rm(cutoffU); rm(cutoffL)

ggplot(dat, aes(x = bwt1)) + geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 70)) + 
  ggtitle("SARAS Kids: Birthweight (g)") +
  xlab("Birthweight (g)")



#### Birth length ####
summary(dat$blen1)
# 164 (23%) do not have any birth length recorded
sd(dat$blen1, na.rm = T)

# check for outliers
cutoffU <- median(dat$blen1, na.rm = T) + 5*mad(dat$blen1, na.rm = T)
cutoffL <- median(dat$blen1, na.rm = T) - 5*mad(dat$blen1, na.rm = T)
sum(dat$blen1 > cutoffU, na.rm = T)
sum(dat$blen1 < cutoffL, na.rm = T)
# no outliers based on the criteria established in the Analysis Plan
rm(cutoffU); rm(cutoffL)

ggplot(dat, aes(x = blen1)) + geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 70)) + 
  ggtitle("SARAS Kids: Birth length (cm)") +
  xlab("Length (cm)")

  
  
#### Head Circumference ####

summary(dat$bhead1)
# 167 (28%) do not have any birth head circumference recorded
sd(dat$bhead1, na.rm = T)
  
# check for outliers
cutoffU <- median(dat$bhead1, na.rm = T) + 5*mad(dat$bhead1, na.rm = T)
cutoffL <- median(dat$bhead1, na.rm = T) - 5*mad(dat$bhead1, na.rm = T)
sum(dat$bhead1 > cutoffU, na.rm = T)
sum(dat$bhead1 < cutoffL, na.rm = T)
# no outliers based on the criteria established in the Analysis Plan
rm(cutoffU); rm(cutoffL)

ggplot(dat, aes(x = bhead1)) + geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 60)) + 
  ggtitle("SARAS Kids: Head circumference (cm)") +
  xlab("Head circumference (cm)")



#### Chest circumference ####
summary(dat$bchst1)
# 164 (23%) do not have any birth chest circumference recorded
sd(dat$bchst1, na.rm = T)

# check for outliers
cutoffU <- median(dat$bchst1, na.rm = T) + 5*mad(dat$bchst1, na.rm = T)
cutoffL <- median(dat$bchst1, na.rm = T) - 5*mad(dat$bchst1, na.rm = T)
sum(dat$bchst1 > cutoffU, na.rm = T)
sum(dat$bchst1 < cutoffL, na.rm = T)
# 5 outliers based on the criteria established in the Analysis Plan

# create a variable to flag outliers for chest circumference
dat$flagbirthCC <- ifelse(dat$bchst1 > cutoffU | dat$bchst1 < cutoffL, 1, 0)
dat$flagbirthCC <- factor(dat$flagbirthCC, labels = c("No Outlier", "Outlier"))
attributes(dat$flagbirthCC)$label <- "saras - birth CC outlier?"
rm(cutoffU); rm(cutoffL)

# summary and plot (no outlier removal)
summary(dat$bchst1)
ggplot(dat, aes(x = bchst1)) + geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) + 
  ggtitle("SARAS Kids: Chest circumference (cm)") +
  xlab("Chest circumference (cm)")
# excluding the outliers
dat %>% filter(flagbirthCC == "No Outlier") %>% 
  summarise(n = n(), Min = min(bchst1), Mean = mean(bchst1), SD = 
              sd(bchst1), Median = median(bchst1), Q1 = 
              quantile(bchst1, probs = 0.25), Q3 = 
              quantile(bchst1, probs = 0.75))
ggplot(subset(dat, flagbirthCC == "No Outlier"), aes(x = bchst1)) + 
  geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 65)) + 
  ggtitle("SARAS Kids (No Outliers): Chest circumference (cm)") +
  xlab("Chest circumference (cm)")



#### Abdominal circumference ####
summary(dat$bab1)
# 164 (23%) do not have any birth abdominal circumference recorded
sd(dat$bab1, na.rm = T)

# check for outliers
cutoffU <- median(dat$bab1, na.rm = T) + 5*mad(dat$bab1, na.rm = T)
cutoffL <- median(dat$bab1, na.rm = T) - 5*mad(dat$bab1, na.rm = T)
sum(dat$bab1 > cutoffU, na.rm = T)
sum(dat$bab1 < cutoffL, na.rm = T)
# 2 outliers based on the criteria established in the Analysis Plan

# create a variable to flag outliers for abdominal circumference
dat$flagbirthAC <- ifelse(dat$bab1 > cutoffU | dat$bab1 < cutoffL, 1, 0)
dat$flagbirthAC <- factor(dat$flagbirthAC, labels = c("No Outlier", "Outlier"))
attributes(dat$flagbirthAC)$label <- "saras - birth AC outlier?"
rm(cutoffU); rm(cutoffL)

# summary and plot (no outlier removal)
summary(dat$bab1)
ggplot(dat, aes(x = bab1)) + geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) + 
  ggtitle("SARAS Kids: Abdominal circumference (cm)") +
  xlab("Abdominal circumference (cm)")
# excluding the outliers
dat %>% filter(flagbirthAC == "No Outlier") %>% 
  summarise(n = n(), Min = min(bab1), Mean = mean(bab1), SD = 
              sd(bab1), Median = median(bab1), Q1 = 
              quantile(bab1, probs = 0.25), Q3 = 
              quantile(bab1, probs = 0.75))
ggplot(subset(dat, flagbirthAC == "No Outlier"), aes(x = bab1)) + 
  geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 70)) + 
  ggtitle("SARAS Kids (No Outliers): Abdominal circumference (cm)") +
  xlab("Abdominal circumference (cm)")



#### Triceps skinfolds ####
summary(dat$btri1)
# 163 (23%) do not have any birth triceps recorded
sd(dat$btri1, na.rm = T)

# check for outliers
cutoffU <- median(dat$btri1, na.rm = T) + 5*mad(dat$btri1, na.rm = T)
cutoffL <- median(dat$btri1, na.rm = T) - 5*mad(dat$btri1, na.rm = T)
sum(dat$btri1 > cutoffU, na.rm = T)
sum(dat$btri1 < cutoffL, na.rm = T)
# 2 outliers based on the criteria established in the Analysis Plan

# create a variable to flag outliers for triceps
dat$flagbirthTri <- ifelse(dat$btri1 > cutoffU | dat$btri1 < cutoffL, 1, 0)
dat$flagbirthTri <- factor(dat$flagbirthTri, labels = c("No Outlier", "Outlier"))
attributes(dat$flagbirthTri)$label <- "saras - birth triceps outlier?"
rm(cutoffU); rm(cutoffL)

# summary and plot (no outlier removal)
summary(dat$btri1)
ggplot(dat, aes(x = btri1)) + geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) + 
  ggtitle("SARAS Kids: Triceps (mm)") +
  xlab("Triceps (mm)")
# excluding the outliers
dat %>% filter(flagbirthTri == "No Outlier") %>% 
  summarise(n = n(), Min = min(btri1), Mean = mean(btri1), SD = 
              sd(btri1), Median = median(btri1), Q1 = 
              quantile(btri1, probs = 0.25), Q3 = 
              quantile(btri1, probs = 0.75))
ggplot(subset(dat, flagbirthTri == "No Outlier"), aes(x = btri1)) + 
  geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 70)) + 
  ggtitle("SARAS Kids (No Outliers): Triceps (mm)") +
  xlab("Triceps (mm)")



#### Subscapular skinfold ####
summary(dat$bsub1)
# 164 (23%) do not have any birth subscapular skinfold recorded
sd(dat$bsub1, na.rm = T)

# check for outliers
cutoffU <- median(dat$bsub1, na.rm = T) + 5*mad(dat$bsub1, na.rm = T)
cutoffL <- median(dat$bsub1, na.rm = T) - 5*mad(dat$bsub1, na.rm = T)
sum(dat$bsub1 > cutoffU, na.rm = T)
sum(dat$bsub1 < cutoffL, na.rm = T)
# 2 outliers based on the criteria established in the Analysis Plan

# create a variable to flag outliers for subscapular
dat$flagbirthSub <- ifelse(dat$bsub1 > cutoffU | dat$bsub1 < cutoffL, 1, 0)
dat$flagbirthSub <- factor(dat$flagbirthSub, labels = c("No Outlier", "Outlier"))
attributes(dat$flagbirthSub)$label <- "saras - birth subscapular outlier?"
rm(cutoffU); rm(cutoffL)

# summary and plot (no outlier removal)
summary(dat$bsub1)
ggplot(dat, aes(x = bsub1)) + geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 80)) + 
  ggtitle("SARAS Kids: Subscapular skinfold (mm)") +
  xlab("Subscapular skinfold (mm)")
# excluding the outliers
dat %>% filter(flagbirthSub == "No Outlier") %>% 
  summarise(n = n(), Min = min(bsub1), Mean = mean(bsub1), SD = 
              sd(bsub1), Median = median(bsub1), Q1 = 
              quantile(bsub1, probs = 0.25), Q3 = 
              quantile(bsub1, probs = 0.75))
ggplot(subset(dat, flagbirthSub == "No Outlier"), aes(x = bsub1)) + 
  geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 70)) + 
  ggtitle("SARAS Kids (No Outliers): Subscapular skinfold (mm)") +
  xlab("Subscapular skinfold (mm)")



#### MUAC ####
summary(dat$bmuac1)
# 164 (23%) do not have any birth MUAC recorded
sd(dat$bmuac1, na.rm = T)

# check for outliers
cutoffU <- median(dat$bmuac1, na.rm = T) + 5*mad(dat$bmuac1, na.rm = T)
cutoffL <- median(dat$bmuac1, na.rm = T) - 5*mad(dat$bmuac1, na.rm = T)
sum(dat$bmuac1 > cutoffU, na.rm = T)
sum(dat$bmuac1 < cutoffL, na.rm = T)
# no outliers based on the criteria established in the Analysis Plan
rm(cutoffU); rm(cutoffL)

ggplot(dat, aes(x = bmuac1)) + geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 70)) + 
  ggtitle("SARAS Kids: MUAC (cm)") +
  xlab("MUAC (cm)")



#### Gestational age ####

summary(dat$bestgest)
# 22 (3%) do not have gestational age recorded
sd(dat$bestgest, na.rm = T)

# 48.29 is too high. Set as missing
dat[which(dat$bestgest > 45), "bestgest"] <- NA

ggplot(dat, aes(x = bestgest)) + geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 170)) + 
  ggtitle("SARAS Kids: Gestational age (weeks)") +
  xlab("Gestational age (weeks)")



#########################################################################################################
#### Multivariate Analysis ##############################################################################
#########################################################################################################

# 1. Check each variable in relation with others
# 2. Check if same measures are consistent


#### Birth weight and length ####
with(dat, plot(bwt1, blen1))

#### Birth abdominal circumference and skinfolds ####
with(dat, plot(btri1, bab1))
with(dat, plot(btri1, bsub1))
with(dat, plot(bab1, bsub1))

#### Birth abdominal circumference and chest circumference ####
with(dat, plot(bchst1, bab1))



#########################################################################################################
#### Derive Measures ####################################################################################
#########################################################################################################

# 1. Compute derived measures listed in the analysis plan
# 2. Check if same measures are consistent
# If applicable:
# 3. Check for outlier using the criteria established in the analysis plan (median +/- 5MAD)
# 4. Retain the outliers but create a variable with 0/1 (named flagvariablename) where 1 indicated
#    the outlier based on the analysis plan.



#### Preterm ####
# preterm is defined as gestational age < 37 weeks
dat$preterm <- ifelse(dat$bestgest < 37, 1, 0)
table(dat$preterm)
dat$preterm <- factor(dat$preterm, labels = c("No", "Yes"))
addmargins(table(dat$preterm)) # available only for 686 people

attributes(dat$preterm)$label <- "saras - Preterm"


#### Low birth weight ####
# low birth weight is defined as weight < 2500 g
dat$lbw <- ifelse(dat$bwt1 < 2500, 1, 0)
table(dat$lbw)
dat$lbw <- factor(dat$lbw, labels = c("No", "Yes"))
addmargins(table(dat$lbw)) # available only for 553 people

attributes(dat$lbw)$label <- "saras - Low birth weight"


#### SGA ####
# Small for gestational age defined using Intergrowth data: José Villar et al. for the 
# International Fetal and Newborn Growth Consortium for the 21st Century (INTERGROWTH-21st)

# SGA was defined as birth weight <10th percentile for gestation age and sex using INTERGROWTH standards

# downlaod the INTERGROWTH Tables (centiles of weight by sex)
# and use them to find the 10th centiles. Table downloaded on June 20, 2017.

# Import INTERGROWTH Standard
datIntergrowth <- read_excel("Intergrowth10CentileBW_20160620.xlsx")
# change name of gestational age to bestgest
names(datIntergrowth)[names(datIntergrowth) == "gestation(week)"] <- "bestgest"
names(datIntergrowth)[names(datIntergrowth) == "10thCentileBoys"] <- "CentileBoys"
names(datIntergrowth)[names(datIntergrowth) == "10thCentileGirls"] <- "CentileGirls"
datIntergrowth <- subset(datIntergrowth, 
                         select = c("bestgest", "CentileBoys", "CentileGirls"))
dat$bestgest <- round(dat$bestgest, 3)
datIntergrowth$bestgest <- round(datIntergrowth$bestgest, 3)
# join the two dataset
dat <- left_join(dat, datIntergrowth, by = 'bestgest', type = "all")

# create birth weight in kg
dat$bwt_kg <- (dat$bwt1)/1000
summary(dat$bwt_kg)
attributes(dat$bwt_kg)$label <- "saras - Birt hweight (kg)"

# create sga variable
for(i in 1:nrow(dat)){
  # set to missing all the boys centiles if subject is a girl
  dat$CentileBoys[i] <- ifelse(dat$sex[i] == "Girl", NA, dat$CentileBoys[i])
  # set to missing all the girl centiles if subject is a boy
  dat$CentileGirls[i] <- ifelse(dat$sex[i] == "Boy", NA, dat$CentileGirls[i])
  
}
# merge the centiles boys and girls in a unique variable
dat <- data.frame(dat)
dat$Centile <- dat$CentileBoys
dat$Centile[is.na(dat$CentileBoys)] <- dat$CentileGirls[is.na(dat$CentileBoys)]
# create sga
dat$sga <- ifelse(dat$bwt_kg < dat$Centile, 1, 0)
table(dat$sga)
dat$sga <- factor(dat$sga, labels = c("No", "Yes"))
addmargins(table(dat$sga)) # available only for 535 people
sum(!is.na(dat$bestgest) & !is.na(dat$bwt1))
# 537 people have gestational age and borth weight recorded.
# 2 where assigned as missing as their gestation is greater 
# than 42.714 days (maximum allowed in Intergrowth)
# Look at those 7
dat[!is.na(dat$bestgest) & !is.na(dat$bwt1) & is.na(dat$sga) & is.na(dat$Centile), 
    c("serno","sex", "bestgest", "CentileGirls", "CentileBoys","bwt1")]
# two have gestation greater than 42.714 weeks (serno 409 boy and serno 1885 girl). 
dat[dat$serno == 49, "Centile"] <- 3.21
dat[dat$serno == 1885, "Centile"] <- 3.03
# recompute small for gestational age
dat$sga <- ifelse(dat$bwt_kg < dat$Centile, 1, 0)
table(dat$sga)
dat$sga <- factor(dat$sga, labels = c("No SGA", "SGA"))
addmargins(table(dat$sga)) # available only for 535 people

attributes(dat$sga)$label <- "saras - SGA (Intergrowth standard)"

#########################################################################################################
#### Save Dataset #######################################################################################
#########################################################################################################

# select only relevant variable (based on the analysis plan)
dat <- dat %>% select(serno, bwt1, blen1, bmuac1, bhead1, bab1, bchst1, bsub1, btri1, 
                      allocation, sex, flagbirthCC, flagbirthAC, flagbirthTri, 
                      flagbirthSub, lbw, bestgest, preterm, sga)  
names(dat) <- c("serno", "birthWeight", "birthLength", "birthMUAC", "birthHC", "birthAC",
                "birthCC", "birthSub", "birthTri", "allocation", "csex", "flagbirthCC", "flagbirthAC", 
                "flagbirthTri", "flagbirthSub", "lbw", "ga", "preterm", "sga")
# add attributes
attributes(dat$serno)$label <- "saras - Subject ID"
attributes(dat$birthWeight)$label <- "saras - Birth weight (g)"
attributes(dat$birthLength)$label <- "saras - Birth length (cm)"
attributes(dat$birthMUAC)$label <- "saras - Birth MUAC (cm)"
attributes(dat$birthHC)$label <- "saras - Birth HC (cm)"
attributes(dat$birthAC)$label <- "saras - Birth AC (cm)"
attributes(dat$birthCC)$label <- "saras - Birth CC (cm)"
attributes(dat$birthSub)$label <- "saras - Birth Subscapular Skinfold (mm)"
attributes(dat$birthTri)$label <- "saras - Birth Triceps (mm)"
attributes(dat$ga)$label <- "saras - Gestational Age (weeks)"
attributes(dat$allocation)$label <- "saras - allocation group"
attributes(dat$csex)$label <- "saras - newborn sex"

# save as csv and as stata file
write.table(dat, "./CleanedDataset/SARASKidsBirthv1_20170622.csv", row.names = FALSE)
write_dta(dat, "./CleanedDataset/SARASKidsBirthv1_20170622.dta")
