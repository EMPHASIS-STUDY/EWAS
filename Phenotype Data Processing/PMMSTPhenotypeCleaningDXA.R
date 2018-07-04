###################################################################################################################
##### PMMST: Phenotype Data Processing Clinic Visit 2 DXA  ########################################################
###################################################################################################################

# 28/01/2018

# load packages
library(haven)
library(readxl)
library(ggplot2)
library(dplyr)

# Data were sent by Gail on 23/01/2018. 

########################################################################################################
#### DXA BODY COMPOSITION ############################################################################## 
########################################################################################################

# load dataset
# load sheet with total fat mass and lean mass first
datFat <- read_excel("./PMMST Outcome Data/CDBH Emphasis data January 24 2018.xlsx", 
                  sheet = " All WB data 19 Jan",
                  # add column names and specify column types
                  col_names = TRUE, col_types = c("text", "date", "numeric", 
                                                  "text", "text", rep("numeric", 6)))
str(datFat)
# 3 children did not have DXA scan as they were afraid of the machine (EMP101D, EMP159G, EMP290I). 
# Set their day of the analysis to missing
datFat[datFat$ID %in% c("EMP101D", "EMP159G", "EMP290I"), "Measure Date"] <- NA
# remove white space from columns name.
colnames(datFat) <- gsub(" ", "", colnames(datFat))
# remove (g) from fat and lean measure to helps with the analsyis
colnames(datFat) <- gsub("(g)", "", colnames(datFat), fixed=TRUE)
# remove cm and kg from height and weight
datFat[,c("HeightatExam", "WeightatExam")] <- lapply(datFat[,c("HeightatExam", "WeightatExam")], 
                                                  function(x) as.numeric(sub("\\s+\\D+$", "", x)))

# add information on android and gynoid fat
datGA <- read_excel("./PMMST Outcome Data/CDBH Emphasis data January 24 2018.xlsx", 
                  sheet = "Other DXA data 23 Jan",
                  # add column names and specify column types
                  col_names = TRUE, col_types = c("text", "date", "numeric", 
                                                  "text", "text", rep("numeric", 2))) 
# remove white space from columns name.
colnames(datGA) <- gsub(" ", "", colnames(datGA))
# remove (g) from fat and lean measure to helps with the analsyis
colnames(datGA) <- gsub("(g)", "", colnames(datGA), fixed=TRUE)
# take only variables that are needed for the analysis
datGA <- datGA %>% dplyr::select(ID, AndroidFatMass, GynoidFatMass)
# merge with information on fat and lean mass
datFat <- inner_join(datFat, datGA)
# remove datGA
rm(datGA)

##################################################################################################################

# quickly check height and weight
summary(datFat$HeightatExam)
hist(datFat$HeightatExam)
summary(datFat$WeightatExam)
hist(datFat$WeightatExam)
with(datFat, plot(WeightatExam, HeightatExam))

# Check if the weight is similar to the total total mass
temp <- data.frame(ID = datFat$ID, HDXA = datFat$HeightatExam, 
                   WDXA = datFat$WeightatExam, tocheck = datFat$TotalTotalMass/1000)
temp$diff <- temp$WDXA - temp$tocheck
summary(temp$diff)
hist(temp$diff, main = "Difference between Weight and Total Total Mass")
# take those with a difference greater than 1 kg
temp[temp$diff > 1 & !(is.na(temp$diff)),]
# EMP010A and EMP050Q. Color code them in the height vs weight plot
# Create column filled with default colour
temp$col <- "black"
# Set new column values to appropriate colours
temp$col[temp$ID %in% c("EMP010A", "EMP050Q")]="red"
with(temp, plot(WDXA, HDXA, pch = 19, col = col)) 
# one is one the the heaviest, the other one is the cloud of points.
# Queries those children back to Gail.
rm(temp)

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

#########################################################################################################
#### Measured outcome variables #########################################################################
#########################################################################################################

#### Total fat mass ####

# transform fat mass in kg
datFat$totfatmass_kg <- datFat$TotalFatMass/1000
attributes(datFat$totfatmass_kg)$label <- "pmmst - fat mass (kg)"
# summary statistics and histogram
summary(datFat$totfatmass_kg)
hist(datFat$totfatmass_kg)

# graph fat and weight
with(datFat, plot(totfatmass_kg, WeightatExam))

# check for outliers using median +/- 5 mad
cutoffU <- median(datFat$totfatmass_kg, na.rm = T) + 5*mad(datFat$totfatmass_kg, na.rm = T)
cutoffL <- median(datFat$totfatmass_kg, na.rm = T) - 5*mad(datFat$totfatmass_kg, na.rm = T)
sum(datFat$totfatmass_kg > cutoffU, na.rm = T)
sum(datFat$totfatmass_kg < cutoffL, na.rm = T)
# 2 outliers based on the analysis plan

# create a flag variable for fat mass and check the outliers
datFat$flagtotfatmass_kg <- ifelse(datFat$totfatmass_kg > cutoffU, 1, 0)
datFat$flagtotfatmass_kg <- factor(datFat$flagtotfatmass_kg, labels = c("No Outlier", "Outlier"))
attributes(datFat$flagtotfatmass_kg)$label <- "pmmst - fat mass outlier?"
datFat[datFat$flagtotfatmass_kg == "Outlier", ] # 2 of the heaviest children in the study

# summary and plot
# including the outliers
summary(datFat$totfatmass_kg)
ggplot(datFat, aes(x = totfatmass_kg)) + geom_histogram(colour="black", fill="white", binwidth = 0.5) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 60)) + ggtitle("PMMST: Total Fat Mass (kg)") +
  xlab("Total Fat Mass (kg)")
# excluding the outliers
datFat %>% filter(flagtotfatmass_kg == "No Outlier") %>% 
  summarise(n = n(), Min = min(totfatmass_kg), Mean = mean(totfatmass_kg), SD = 
              sd(totfatmass_kg), Median = median(totfatmass_kg), Q1 = 
              quantile(totfatmass_kg, probs = 0.25), Q3 = 
              quantile(totfatmass_kg, probs = 0.75))
ggplot(subset(datFat, flagtotfatmass_kg == "No Outlier"), aes(x = totfatmass_kg)) + 
  geom_histogram(colour="black", fill="white", binwidth = 0.5) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 60)) + 
  ggtitle("PMMST (No Outliers): Total Fat Mass (kg)") +
  xlab("Total Fat Mass (kg)")


#### Total lean mass ####

# transform lean mass in kg
datFat$totleanmass_kg <- datFat$TotalLeanMass/1000
attributes(datFat$totleanmass_kg)$label <- "pmmst - lean mass (kg)"
# summary statistics and histogram
summary(datFat$totleanmass_kg)
hist(datFat$totleanmass_kg)

# check for outliers using median +/- 5 mad
cutoffU <- median(datFat$totleanmass_kg, na.rm = T) + 5*mad(datFat$totleanmass_kg, na.rm = T)
cutoffL <- median(datFat$totleanmass_kg, na.rm = T) - 5*mad(datFat$totleanmass_kg, na.rm = T)
sum(datFat$totleanmass_kg > cutoffU, na.rm = T)
sum(datFat$totleanmass_kg < cutoffL, na.rm = T)
# 0 outliers are identified

# graph lean and weight
with(datFat, plot(totleanmass_kg, WeightatExam))

# plot
ggplot(datFat, aes(x = totleanmass_kg)) + 
  geom_histogram(colour="black", fill="white", binwidth = 1) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 60)) + ggtitle("PMMST: Total Lean Mass (kg)") +
  xlab("Total Lean Mass (kg)")


#### Android Fat ####

# transform android fat in kg
datFat$androidfatmass_kg <- datFat$AndroidFatMass/1000
attributes(datFat$androidfatmass_kg)$label <- "pmmst - android fat mass (kg)"
# summary statistics and histogram
summary(datFat$androidfatmass_kg)
hist(datFat$androidfatmass_kg)

# check for outliers using median +/- 5 mad
cutoffU <- median(datFat$androidfatmass_kg, na.rm = T) + 5*mad(datFat$androidfatmass_kg, na.rm = T)
cutoffL <- median(datFat$androidfatmass_kg, na.rm = T) - 5*mad(datFat$androidfatmass_kg, na.rm = T)
sum(datFat$androidfatmass_kg > cutoffU, na.rm = T)
sum(datFat$androidfatmass_kg < cutoffL, na.rm = T)
# 5 were identified as outliers based on the analysis plan.

# create a flag variable for android mass and check the outliers
datFat$flagandroidfatmass_kg <- ifelse(datFat$androidfatmass_kg > cutoffU, 1, 0)
datFat$flagandroidfatmass_kg <- factor(datFat$flagandroidfatmass_kg, labels = c("No Outlier", "Outlier"))
attributes(datFat$flagandroidfatmass_kg)$label <- "pmmst - android fat outlier?"
datFat[datFat$flagandroidfatmass_kg == "Outlier" & !is.na(datFat$androidfatmass_kg), ] 
# 2 of the 5 had also outlier fat mass.

# summary and plot
# including the outliers
summary(datFat$androidfatmass_kg)
ggplot(datFat, aes(x = androidfatmass_kg)) + geom_histogram(colour="black", fill="white", binwidth = 0.02) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 60)) + ggtitle("PMMST: Android Fat (kg)") +
  xlab("Android Fat (kg)")
# excluding the outliers
datFat %>% filter(flagandroidfatmass_kg == "No Outlier") %>% 
  summarise(n = n(), Min = min(androidfatmass_kg), Mean = mean(androidfatmass_kg), SD = 
              sd(androidfatmass_kg), Median = median(androidfatmass_kg), Q1 = 
              quantile(androidfatmass_kg, probs = 0.25), Q3 = 
              quantile(androidfatmass_kg, probs = 0.75))
ggplot(subset(datFat, flagandroidfatmass_kg == "No Outlier"), aes(x = androidfatmass_kg)) + 
  geom_histogram(colour="black", fill="white", binwidth = 0.018) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 60)) + 
  ggtitle("PMMST (No Outliers): Android Fat (kg)") +
  xlab("Android Fat (kg)")


#### Gynoid Fat ####

# transform gynoid fat in kg
datFat$gynoidfatmass_kg <- datFat$GynoidFatMass/1000
attributes(datFat$gynoidfatmass_kg)$label <- "pmmst - gynoid fat mass (kg)"
# summary statistics and histogram
summary(datFat$gynoidfatmass_kg)
hist(datFat$gynoidfatmass_kg)

# check for outliers using median +/- 5 mad
cutoffU <- median(datFat$gynoidfatmass_kg, na.rm = T) + 5*mad(datFat$gynoidfatmass_kg, na.rm = T)
cutoffL <- median(datFat$gynoidfatmass_kg, na.rm = T) - 5*mad(datFat$gynoidfatmass_kg, na.rm = T)
sum(datFat$gynoidfatmass_kg > cutoffU, na.rm = T)
sum(datFat$gynoidfatmass_kg < cutoffL, na.rm = T)
# 2 were identified as outliers based on the analysis plan.

# create a flag variable for gynoid fat and check the outliers
datFat$flaggynoidfatmass_kg <- ifelse(datFat$gynoidfatmass_kg > cutoffU, 1, 0)
datFat$flaggynoidfatmass_kg <- factor(datFat$flaggynoidfatmass_kg, labels = c("No Outlier", "Outlier"))
attributes(datFat$flaggynoidfatmass_kg)$label <- "pmmst - gynoid fat outlier?"
datFat[datFat$flaggynoidfatmass_kg == "Outlier" & !is.na(datFat$gynoidfatmass_kg), ] 
# 2 that also have high fat mass.

# summary and plot
# including the outliers
summary(datFat$gynoidfatmass_kg)
ggplot(datFat, aes(x = gynoidfatmass_kg)) + geom_histogram(colour="black", fill="white", binwidth = 0.1) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 60)) + ggtitle("PMMST: Gynoid Fat (kg)") +
  xlab("Gynoid Fat (kg)")
# excluding the outliers
datFat %>% filter(flaggynoidfatmass_kg == "No Outlier") %>% 
  summarise(n = n(), Min = min(gynoidfatmass_kg), Mean = mean(gynoidfatmass_kg), SD = 
              sd(gynoidfatmass_kg), Median = median(gynoidfatmass_kg), Q1 = 
              quantile(gynoidfatmass_kg, probs = 0.25), Q3 = 
              quantile(gynoidfatmass_kg, probs = 0.75))
ggplot(subset(datFat, flaggynoidfatmass_kg == "No Outlier"), aes(x = gynoidfatmass_kg)) + 
  geom_histogram(colour="black", fill="white", binwidth = 0.1) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 60)) + 
  ggtitle("PMMST (No Outliers): Gynoid Fat (kg)") +
  xlab("Gynoid Fat (kg)")


#########################################################################################################
#### Derived outcome variables ##########################################################################
#########################################################################################################

#### fat mass index (FMI): fat mass (kg)/height (m)^2 ####
# height is in cm, need to transform in m. Can do directly in the FMI formula
datFat$fmi <- datFat$totfatmass_kg/((0.01*datFat$HeightatExam)^2)
attributes(datFat$fmi)$label <- "pmmst - fat mass index (kg/m2)" 

# summary and plot 
summary(datFat$fmi)
hist(datFat$fmi)

# check for outliers
cutoffU <- median(datFat$fmi, na.rm = T) + 5*mad(datFat$fmi, na.rm = T)
cutoffL <- median(datFat$fmi, na.rm = T) - 5*mad(datFat$fmi, na.rm = T)
sum(datFat$fmi > cutoffU, na.rm = T)
sum(datFat$fmi < cutoffL, na.rm = T)
# two outliers have been identified

# create a flag variable for fat mass and check the outliers
datFat$flagfmi <- ifelse(datFat$fmi > cutoffU, 1, 0)
datFat$flagfmi <- factor(datFat$flagfmi, labels = c("No Outlier", "Outlier"))
attributes(datFat$flagfmi)$label <- "pmmst - fmi outlier?"
datFat[datFat$flagfmi == "Outlier", ] 
# EMP278W had also outlier fat mass. EMP168E had fat mass (altough almost 2 times
# higher than the median) within range but high FMI

# plot: including the outliers
ggplot(datFat, aes(x = fmi)) + geom_histogram(colour="black", fill="white", binwidth = 0.4) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 90)) + 
  ggtitle(expression(paste("PMMST: Fat Mass Index (", kg/m^{2},")"))) +
  xlab(expression(paste("Fat Mass Index (", kg/m^{2}, ")")))

# excluding the outliers
datFat %>% filter(flagfmi == "No Outlier") %>% 
  summarise(n = n(), Min = min(fmi), Mean = mean(fmi), SD = 
              sd(fmi), Median = median(fmi), Q1 = 
              quantile(fmi, probs = 0.25), Q3 = 
              quantile(fmi, probs = 0.75))
ggplot(subset(datFat, flagfmi == "No Outlier"), aes(x = fmi)) + 
  geom_histogram(colour="black", fill="white", binwidth = 0.4) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 90)) + 
  ggtitle(expression(paste("PMMST (No Outliers): Fat Mass Index (", kg/m^{2},")"))) +
  xlab(expression(paste("Fat Mass Index (", kg/m^{2}, ")")))


#### lean mass index (FMI): fat mass (kg)/height (m)^2 ####

datFat$lmi <- datFat$totleanmass_kg/((0.01*datFat$HeightatExam)^2)
attributes(datFat$lmi)$label <- "pmmst - lean mass index (kg/m2)"
summary(datFat$lmi)
hist(datFat$lmi)

# check for outliers
cutoffU <- median(datFat$lmi, na.rm = T) + 5*mad(datFat$lmi, na.rm = T)
cutoffL <- median(datFat$lmi, na.rm = T) - 5*mad(datFat$lmi, na.rm = T)
sum(datFat$lmi > cutoffU, na.rm = T)
sum(datFat$lmi < cutoffL, na.rm = T)
# no outliers have been identified

# plot
ggplot(datFat, aes(x = lmi)) + 
  geom_histogram(colour="black", fill="white", binwidth = 0.4) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 50)) + 
  ggtitle(expression(paste("PMMST: Lean Mass Index (", kg/m^{2},")"))) +
  xlab(expression(paste("Lean Mass Index (", kg/m^{2}, ")")))


#### Body Fat % ((Fat Mass x 100) / weight) ####

datFat$fatper <- with(datFat, 100*(totfatmass_kg/WeightatExam))
attributes(datFat$fatper)$label <- "pmmst - body fat %"

# summary and plot (including outliers)
summary(datFat$fatper)
hist(datFat$fatper)

# check for outliers
cutoffU <- median(datFat$fatper, na.rm = T) + 5*mad(datFat$fatper, na.rm = T)
cutoffL <- median(datFat$fatper, na.rm = T) - 5*mad(datFat$fatper, na.rm = T)
sum(datFat$fatper > cutoffU, na.rm = T)
sum(datFat$fatper < cutoffL, na.rm = T)
# 0 outliers identified
# take the one with the highest fat %
datFat[which.max(datFat$fatper),] 
# was identifed as outlier for fat mass.

# plot
ggplot(datFat, aes(x = fatper)) + geom_histogram(colour="black", fill="white", binwidth = 1.5) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 60)) + 
  ggtitle("PMMST: Body fat %") +
  xlab("Body fat %")


########################################################################################################
#### DXA BONE DATA ##################################################################################### 
########################################################################################################

# import data on BMD, BMC and BA
datBone <- read_excel("./PMMST Outcome Data/CDBH Emphasis data January 24 2018.xlsx", 
                  sheet = "WBLH data 23 Jan",
                  # add column names and specify column types
                  col_names = TRUE, col_types = c("text", "date", "numeric", 
                                                  "text", "text", rep("numeric", 6)))
str(datBone)
# 3 children did not have DXA scan as they were afraid of the machine (EMP101D, EMP159G, EMP290I). 
# Set their day of the analysis to missing
# set to missing also EMP182C. Gail is investigating why the bone measure is messing
datBone[datBone$ID %in% c("EMP101D", "EMP159G", "EMP290I", "EMP182C"), "Measure Date"] <- NA
# remove white space from columns name.
colnames(datBone) <- gsub(" ", "", colnames(datBone))
# remove cm and kg from height and weight
datBone[,c("HeightatExam", "WeightatExam")] <- lapply(datBone[,c("HeightatExam", "WeightatExam")], 
                                                  function(x) as.numeric(sub("\\s+\\D+$", "", x)))

# Add data on spine BMC, BMD and BA
datSpine <- read_excel("./PMMST Outcome Data/CDBH Emphasis data January 24 2018.xlsx", 
                      sheet = "Spine data 23 Jan",
                      # add column names and specify column types
                      col_names = TRUE, col_types = c("text", "date", "numeric", 
                                                      "text", "text", rep("numeric", 15)))
str(datSpine)
# 3 children did not have DXA scan as they were afraid of the machine (EMP101D, EMP159G, EMP290I). 
# Set their day of the analysis to missing
# set to missing also EMP182C. Gail is investigating why the bone measure is messing
datSpine[datSpine$ID %in% c("EMP101D", "EMP159G", "EMP290I", "EMP182C"), "Measurement Date"] <- NA
# remove white space from columns name.
colnames(datSpine) <- gsub(" ", "", colnames(datSpine))
# remove -  from columns name.
colnames(datSpine) <- gsub("-", "", colnames(datSpine))
# delete variable we do not need
datSpine <- datSpine %>% select(-one_of(c("MeasurementDate", "Age", "HeightatExam", "WeightatExam")))
datBone <- inner_join(datBone, datSpine)
# remove datSpine
rm(datSpine)

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

#### Total body bone mineral content (BMC, whole body minus the head) ####

summary(datBone$TBLHBMC)
hist(datBone$TBLHBMC) # EMP009Y has the highest BMC
attributes(datBone$TBLHBMC)$label <- "pmmst - total BMC excluding head (g)"

# check for outliers
cutoffU <- median(datBone$TBLHBMC, na.rm = T) + 5*mad(datBone$TBLHBMC, na.rm = T)
cutoffL <- median(datBone$TBLHBMC, na.rm = T) - 5*mad(datBone$TBLHBMC, na.rm = T)
sum(datBone$TBLHBMC > cutoffU, na.rm = T)
sum(datBone$TBLHBMC < cutoffL, na.rm = T)
# 0 outliers identified

# plot
ggplot(datBone, aes(x = TBLHBMC)) + geom_histogram(colour="black", fill="white", binwidth = 1.5) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 60)) + 
  ggtitle("PMMST: Total body BMC excluding head (g)") +
  xlab("Total body BMC excluding head (g)")


#### Total body bone mineral area (BA, whole body minus the head) ####

summary(datBone$TBLHArea)
hist(datBone$TBLHArea) # EMP009Y has the highest BMC
attributes(datBone$TBLHArea)$label <- "pmmst - total BA excluding head (cm^2)"

# check for outliers
cutoffU <- median(datBone$TBLHArea, na.rm = T) + 5*mad(datBone$TBLHArea, na.rm = T)
cutoffL <- median(datBone$TBLHArea, na.rm = T) - 5*mad(datBone$TBLHArea, na.rm = T)
sum(datBone$TBLHArea > cutoffU, na.rm = T)
sum(datBone$TBLHArea < cutoffL, na.rm = T)
# 0 outliers identified

# plot
ggplot(datBone, aes(x = TBLHArea)) + geom_histogram(colour="black", fill="white", binwidth = 45) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 60))  + 
  ggtitle(expression(paste("PMMST: Total body BA (", cm^{2},")"))) +
  xlab(expression(paste("Total body BA (", cm^{2},")")))


#### Total body bone mineral density (BMD, whole body minus the head) ####

summary(datBone$TBLHBMD)
hist(datBone$TBLHBMD) 
attributes(datBone$TBLHBMD)$label <- "pmmst - total BMD excluding head (g/cm^2)"

# check for outliers
cutoffU <- median(datBone$TBLHBMD, na.rm = T) + 5*mad(datBone$TBLHBMD, na.rm = T)
cutoffL <- median(datBone$TBLHBMD, na.rm = T) - 5*mad(datBone$TBLHBMD, na.rm = T)
sum(datBone$TBLHBMD > cutoffU, na.rm = T)
sum(datBone$TBLHBMD < cutoffL, na.rm = T)
# 0 outliers identified

# plot
ggplot(datBone, aes(x = TBLHBMD)) + geom_histogram(colour="black", fill="white", binwidth = 0.025) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 60))  + 
  ggtitle(expression(paste("PMMST: Total body BMD (g/", cm^{2},")"))) +
  xlab(expression(paste("Total body BMD (g/", cm^{2},")")))


#### Spine BMC ####

summary(datBone$L1L4BMC)
hist(datBone$L1L4BMC) 
names(datBone)[names(datBone) == "L1L4BMC"] <- "cspineBMC"
attributes(datBone$cspineBMC)$label <- "pmmst - spine BMC"

# check for outliers
cutoffU <- median(datBone$cspineBMC, na.rm = T) + 5*mad(datBone$cspineBMC, na.rm = T)
cutoffL <- median(datBone$cspineBMC, na.rm = T) - 5*mad(datBone$cspineBMC, na.rm = T)
sum(datBone$cspineBMC > cutoffU, na.rm = T)
sum(datBone$cspineBMC < cutoffL, na.rm = T)
# 1 outlier identified

# create a flag variable for spine BMC and check the outliers
datBone$flagcspineBMC <- ifelse(datBone$cspineBMC > cutoffU, 1, 0)
datBone$flagcspineBMC <- factor(datBone$flagcspineBMC, labels = c("No Outlier", "Outlier"))
attributes(datBone$flagcspineBMC)$label <- "pmmst - spine BMC outlier?"
datFat[datBone$flagcspineBMC == "Outlier", ] # EMP009Y who also has higher BMD and BMC

# plot (with outlier)
ggplot(datBone, aes(x = cspineBMC)) + geom_histogram(colour="black", fill="white", binwidth = 1.5) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 60)) + 
  ggtitle("PMMST: Spine BMC (g)") +
  xlab("Spine BMC (g)")

# summary and plot (without outlier)
datBone %>% filter(flagcspineBMC == "No Outlier") %>% 
  summarise(n = n(), Min = min(cspineBMC), Mean = mean(cspineBMC), SD = 
              sd(cspineBMC), Median = median(cspineBMC), Q1 = 
              quantile(cspineBMC, probs = 0.25), Q3 = 
              quantile(cspineBMC, probs = 0.75))
ggplot(subset(datBone, flagcspineBMC == "No Outlier"), aes(x = cspineBMC)) + 
  geom_histogram(colour="black", fill="white", binwidth = 1.5) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 60))+ 
  ggtitle("PMMST (No Outliers): Spine BMC (g)") +
  xlab("Spine BMC (g)")


#### Spine BA ####

summary(datBone$L1L4Area)
hist(datBone$L1L4Area) 
names(datBone)[names(datBone) == "L1L4Area"] <- "cspineBA"
attributes(datBone$cspineBA)$label <- "pmmst - spine BA"

# check for outliers
cutoffU <- median(datBone$cspineBA, na.rm = T) + 5*mad(datBone$cspineBA, na.rm = T)
cutoffL <- median(datBone$cspineBA, na.rm = T) - 5*mad(datBone$cspineBA, na.rm = T)
sum(datBone$cspineBA > cutoffU, na.rm = T)
sum(datBone$cspineBA < cutoffL, na.rm = T)
# 0 outliers identified

ggplot(datBone, aes(x = cspineBA)) + geom_histogram(colour="black", fill="white", binwidth = 1.35) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 60)) + 
  ggtitle(expression(paste("PMMST: Spine BA (", cm^{2},")"))) +
  xlab(expression(paste("Spine BA (", cm^{2},")")))


#### Spine BMD ####

summary(datBone$L1L4BMD)
hist(datBone$L1L4BMD) 
names(datBone)[names(datBone) == "L1L4BMD"] <- "cspineBMD"
attributes(datBone$cspineBMD)$label <- "pmmst - spine BMD"

# check for outliers
cutoffU <- median(datBone$cspineBMD, na.rm = T) + 5*mad(datBone$cspineBMD, na.rm = T)
cutoffL <- median(datBone$cspineBMD, na.rm = T) - 5*mad(datBone$cspineBMD, na.rm = T)
sum(datBone$cspineBMD > cutoffU, na.rm = T)
sum(datBone$cspineBMD < cutoffL, na.rm = T)
# 0 outliers identified

ggplot(datBone, aes(x = cspineBMD)) + geom_histogram(colour="black", fill="white", binwidth = 0.03) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 60)) + 
  ggtitle(expression(paste("PMMST: Spine BMD (", g/cm^{2},")"))) +
  xlab(expression(paste("Spine BMD (", g/cm^{2},")")))


#########################################################################################################
#### Derived outcome variables ##########################################################################
#########################################################################################################


#### Bone mineral apparent density ####

# Bone mineral apparent density was calculated based on the formula in Ward et al.
# UK reference data for the Hologic QDR Discovery dual-energy
# xray absorptiometry scanner in healthy children and young adults aged 6–17 years (2007)
datBone$cBMAD <- with(datBone,
                      (L1BMC + L2BMC + L3BMC + L4BMC)/(L1Area^(3/2) + L2Area^(3/2) + L3Area^(3/2) + L4Area^(3/2)))
summary(datBone$cBMAD)
attributes(datBone$cBMAD)$label <- "pmmst - spine BMAD"
hist(datBone$cBMAD)

# check for outliers
cutoffU <- median(datBone$cBMAD, na.rm = T) + 5*mad(datBone$cBMAD, na.rm = T)
cutoffL <- median(datBone$cBMAD, na.rm = T) - 5*mad(datBone$cBMAD, na.rm = T)
sum(datBone$cBMAD > cutoffU, na.rm = T)
sum(datBone$cBMAD < cutoffL, na.rm = T)
# 0 outliers identified

ggplot(datBone, aes(x = cBMAD)) + geom_histogram(colour="black", fill="white", binwidth = 0.012) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 60)) + 
  ggtitle(expression(paste("PMMST: Bone Mineral Apparent Density (", g/cm^{3},")"))) +
  xlab(expression(paste("Bone Mineral Apparent Density (", g/cm^{3},")")))


#########################################################################################################
#### Get DXA variables dataset ##########################################################################
#########################################################################################################

datBone <- datBone %>% dplyr::select(ID, TBLHBMD, TBLHBMC, TBLHArea, cspineBMD, cspineBMC, cspineBA, cBMAD, 
                                     flagcspineBMC)
names(datBone) <- c("serno", "cBMD", "cBMC", "cBA", "cspineBMD", "cspineBMC", "cspineBA", "cBMAD", 
                    "flagcspineBMC")
head(datBone, 4)

dat <- read_dta("./Cleaned Dataset/PMMSTGrowthBodyCompv1_20180131.dta")
dat <- dat %>% dplyr::filter(cvisit == 2) %>% dplyr::select(serno, csex, cage, allocation)

datBone <- inner_join(datBone, dat)
attributes(datBone$serno)$label <- "pmmst - Subject ID"
names(datBone)[names(datBone) == "cage"] <- "cageV2"
# add visit data
datBone$cvisit <- 2
attributes(datBone$cvisit)$label <- "pmmst - clinic visit"


#########################################################################################################
#### Add pQCT ###########################################################################################
#########################################################################################################

datPQCT <- read.csv("./PMMST Outcome Data/final_master_emp_11_06.csv")

#### Tibia 8%: total bone content ####
 
# summary statistics and histogram
summary(datPQCT$tot_cnt)
hist(datPQCT$tot_cnt)

# check for outliers using median +/- 5 mad
cutoffU <- median(datPQCT$tot_cnt, na.rm = T) + 5*mad(datPQCT$tot_cnt, na.rm = T)
cutoffL <- median(datPQCT$tot_cnt, na.rm = T) - 5*mad(datPQCT$tot_cnt, na.rm = T)
sum(datPQCT$tot_cnt > cutoffU, na.rm = T)
sum(datPQCT$tot_cnt < cutoffL, na.rm = T)

# plot 
ggplot(datPQCT, aes(x = tot_cnt)) + geom_histogram(colour="black", fill="white", binwidth = 5) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 35)) + 
  ggtitle("PMMST: Tibia 8% Total Bone Content (mg/mm)") +
  xlab("Total Bone Content (mg/mm)")

#### Tibia 8%: total bone density ####

# summary statistics and histogram
summary(datPQCT$tot_den)
hist(datPQCT$tot_den)

# check for outliers using median +/- 5 mad
cutoffU <- median(datPQCT$tot_den, na.rm = T) + 5*mad(datPQCT$tot_den, na.rm = T)
cutoffL <- median(datPQCT$tot_den, na.rm = T) - 5*mad(datPQCT$tot_den, na.rm = T)
sum(datPQCT$tot_den > cutoffU, na.rm = T)
sum(datPQCT$tot_den < cutoffL, na.rm = T)

# plot 
ggplot(datPQCT, aes(x = tot_den)) + geom_histogram(colour="black", fill="white", binwidth = 7.5) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 35)) + 
  ggtitle(expression(paste("PMMST: Tibia 8% Total Bone Density (", mg/cm^{3},")"))) +
  xlab(expression(paste("Tibia 8% Total Bone Density (", mg/cm^{3},")")))

#### Tibia 8%: tibia area ####

# summary statistics and histogram
summary(datPQCT$tot_a_8)
hist(datPQCT$tot_a_8)

# check for outliers using median +/- 5 mad
cutoffU <- median(datPQCT$tot_a_8, na.rm = T) + 5*mad(datPQCT$tot_a_8, na.rm = T)
cutoffL <- median(datPQCT$tot_a_8, na.rm = T) - 5*mad(datPQCT$tot_a_8, na.rm = T)
sum(datPQCT$tot_a_8 > cutoffU, na.rm = T)
sum(datPQCT$tot_a_8 < cutoffL, na.rm = T)

# plot 
ggplot(datPQCT, aes(x = tot_a_8)) + geom_histogram(colour="black", fill="white", binwidth = 25) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 60)) + 
  ggtitle(expression(paste("PMMST: Tibia 8% Area (", mm^{2},")"))) +
  xlab(expression(paste("Tibia 8% Area (", mm^{2},")")))

# check the one wit high tibia area
temp <- datPQCT %>% filter(tot_a_8 == max(tot_a_8, na.rm = T))
rm(temp)

#### Tibia 8%: cortical and subcortical content ####

# summary statistics and histogram
summary(datPQCT$crtsub_cnt)
hist(datPQCT$crtsub_cnt)

# check for outliers using median +/- 5 mad
cutoffU <- median(datPQCT$crtsub_cnt, na.rm = T) + 5*mad(datPQCT$crtsub_cnt, na.rm = T)
cutoffL <- median(datPQCT$crtsub_cnt, na.rm = T) - 5*mad(datPQCT$crtsub_cnt, na.rm = T)
sum(datPQCT$crtsub_cnt > cutoffU, na.rm = T)
sum(datPQCT$crtsub_cnt < cutoffL, na.rm = T)

# plot 
ggplot(datPQCT, aes(x = crtsub_cnt)) + 
  geom_histogram(colour="black", fill="white", binwidth = 5) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 60)) + 
  ggtitle("PMMST: Tibia 8% Cortical and Subcortical Content (mg/mm)") +
  xlab("Cortical and Subcortical Content (mg/mm)")


#### Tibia 8%: cortical and subcortical density ####

# summary statistics and histogram
summary(datPQCT$crtsub_den)
hist(datPQCT$crtsub_den)

# check for outliers using median +/- 5 mad
cutoffU <- median(datPQCT$crtsub_den, na.rm = T) + 5*mad(datPQCT$crtsub_den, na.rm = T)
cutoffL <- median(datPQCT$crtsub_den, na.rm = T) - 5*mad(datPQCT$crtsub_den, na.rm = T)
sum(datPQCT$crtsub_den > cutoffU, na.rm = T)
sum(datPQCT$crtsub_den < cutoffL, na.rm = T)

# plot 
ggplot(datPQCT, aes(x = crtsub_den)) + 
  geom_histogram(colour="black", fill="white", binwidth = 20) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 60)) + 
  ggtitle(expression(paste("PMMST: Tibia 8% Cortical and Subcortical Density (", mg/cm^{3},")"))) +
  xlab(expression(paste("Tibia 8% Cortical and Subcortical Density (", mg/cm^{3},")")))


#### Tibia 8%: cortical and subcortical area ####

# summary statistics and histogram
summary(datPQCT$crtsub_a)
hist(datPQCT$crtsub_a)

# check for outliers using median +/- 5 mad
cutoffU <- median(datPQCT$crtsub_a, na.rm = T) + 5*mad(datPQCT$crtsub_a, na.rm = T)
cutoffL <- median(datPQCT$crtsub_a, na.rm = T) - 5*mad(datPQCT$crtsub_a, na.rm = T)
sum(datPQCT$crtsub_a > cutoffU, na.rm = T)
sum(datPQCT$crtsub_a < cutoffL, na.rm = T)

# plot 
ggplot(datPQCT, aes(x = crtsub_a)) + 
  geom_histogram(colour="black", fill="white", binwidth = 15) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 60)) + 
  ggtitle(expression(paste("PMMST: Tibia 8% Cortical and Subcortical Area (", mm^{2},")"))) +
  xlab(expression(paste("Tibia 8% Cortical and Subcortical Area (", mm^{2},")")))


#### Tibia 8%: trabecular content ####

# summary statistics and histogram
summary(datPQCT$trab_cnt)
hist(datPQCT$trab_cnt)

# check for outliers using median +/- 5 mad
cutoffU <- median(datPQCT$trab_cnt, na.rm = T) + 5*mad(datPQCT$trab_cnt, na.rm = T)
cutoffL <- median(datPQCT$trab_cnt, na.rm = T) - 5*mad(datPQCT$trab_cnt, na.rm = T)
sum(datPQCT$trab_cnt > cutoffU, na.rm = T)
sum(datPQCT$trab_cnt < cutoffL, na.rm = T)

# plot 
ggplot(datPQCT, aes(x = trab_cnt)) + geom_histogram(colour="black", fill="white", binwidth = 1.8) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 35)) +
  ggtitle("PMMST: Tibia 8% Trabecular Content (mg/mm)") +
  xlab("Trabecular Content (mg/mm)")


#### Tibia 8%: trabecular density ####

# summary statistics and histogram
summary(datPQCT$trab_den)
hist(datPQCT$trab_den)

# check for outliers using median +/- 5 mad
cutoffU <- median(datPQCT$trab_den, na.rm = T) + 5*mad(datPQCT$trab_den, na.rm = T)
cutoffL <- median(datPQCT$trab_den, na.rm = T) - 5*mad(datPQCT$trab_den, na.rm = T)
sum(datPQCT$trab_den > cutoffU, na.rm = T)
sum(datPQCT$trab_den < cutoffL, na.rm = T)

# plot 
ggplot(datPQCT, aes(x = trab_den)) + geom_histogram(colour="black", fill="white", binwidth = 6.5) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 35)) + 
  ggtitle(expression(paste("PMMST: Tibia 8% Trabecular Density (", mg/cm^{3},")"))) +
  xlab(expression(paste("Tibia 8% Trabecular Density (", mg/cm^{3},")")))


#### Tibia 8%: trabecular area ####

# summary statistics and histogram
summary(datPQCT$trab_a)
hist(datPQCT$trab_a)

# check for outliers using median +/- 5 mad
cutoffU <- median(datPQCT$trab_a, na.rm = T) + 5*mad(datPQCT$trab_a, na.rm = T)
cutoffL <- median(datPQCT$trab_a, na.rm = T) - 5*mad(datPQCT$trab_a, na.rm = T)
sum(datPQCT$trab_a > cutoffU, na.rm = T)
sum(datPQCT$trab_a < cutoffL, na.rm = T)

# plot 
ggplot(datPQCT, aes(x = trab_a)) + geom_histogram(colour="black", fill="white", binwidth = 10) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 60)) + 
  ggtitle(expression(paste("PMMST: Tibia 8% Trabecular Area (", mm^{2},")"))) +
  xlab(expression(paste("Tibia Trabecular 8% Area (", mm^{2},")")))


#### Tibia 38%: Cortical content ####

# summary statistics and histogram
summary(datPQCT$tib_38_crt_cnt)
hist(datPQCT$tib_38_crt_cnt)

# check for outliers using median +/- 5 mad
cutoffU <- median(datPQCT$tib_38_crt_cnt, na.rm = T) + 5*mad(datPQCT$tib_38_crt_cnt, na.rm = T)
cutoffL <- median(datPQCT$tib_38_crt_cnt, na.rm = T) - 5*mad(datPQCT$tib_38_crt_cnt, na.rm = T)
sum(datPQCT$tib_38_crt_cnt > cutoffU, na.rm = T)
sum(datPQCT$tib_38_crt_cnt < cutoffL, na.rm = T)

# plot 
ggplot(datPQCT, aes(x = tib_38_crt_cnt)) + 
  geom_histogram(colour="black", fill="white", binwidth = 9) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 60)) + 
  ggtitle("PMMST: Tibia 38% Cortical Content (mg/mm)") +
  xlab("Tibia 38% Cortical Content (mg/mm)")


#### Tibia 38%: Cortical density ####

# summary statistics and histogram
summary(datPQCT$tib_38_crt_den)
hist(datPQCT$tib_38_crt_den)

# check for outliers using median +/- 5 mad
cutoffU <- median(datPQCT$tib_38_crt_den, na.rm = T) + 5*mad(datPQCT$tib_38_crt_den, na.rm = T)
cutoffL <- median(datPQCT$tib_38_crt_den, na.rm = T) - 5*mad(datPQCT$tib_38_crt_den, na.rm = T)
sum(datPQCT$tib_38_crt_den > cutoffU, na.rm = T)
sum(datPQCT$tib_38_crt_den < cutoffL, na.rm = T)

# plot 
ggplot(datPQCT, aes(x = tib_38_crt_den)) + 
  geom_histogram(colour="black", fill="white", binwidth = 15) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 60)) + 
  ggtitle(expression(paste("PMMST: Tibia 38% Cortical Density (", mg/cm^{3},")"))) +
  xlab(expression(paste("Tibia 38% Cortical Density (", mg/cm^{3},")")))

#### Tibia 38%: Cortical area ####

# summary statistics and histogram
summary(datPQCT$tib_38_crt_a)
hist(datPQCT$tib_38_crt_a)

# check for outliers using median +/- 5 mad
cutoffU <- median(datPQCT$tib_38_crt_a, na.rm = T) + 5*mad(datPQCT$tib_38_crt_a, na.rm = T)
cutoffL <- median(datPQCT$tib_38_crt_a, na.rm = T) - 5*mad(datPQCT$tib_38_crt_a, na.rm = T)
sum(datPQCT$tib_38_crt_a > cutoffU, na.rm = T)
sum(datPQCT$tib_38_crt_a < cutoffL, na.rm = T)

# plot 
ggplot(datPQCT, aes(x = tib_38_crt_a)) + 
  geom_histogram(colour="black", fill="white", binwidth = 9) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 60)) + 
  ggtitle(expression(paste("PMMST: Tibia 38% Cortical Area (", mm^{2},")"))) +
  xlab(expression(paste("Tibia 38% Cortical Area (", mm^{2},")")))

#### Tibia 38%: Cortical thickness ####

# summary statistics and histogram
summary(datPQCT$tib_38_crt_thk_c)
hist(datPQCT$tib_38_crt_thk_c)

# check for outliers using median +/- 5 mad
cutoffU <- median(datPQCT$tib_38_crt_thk_c, na.rm = T) + 5*mad(datPQCT$tib_38_crt_thk_c, na.rm = T)
cutoffL <- median(datPQCT$tib_38_crt_thk_c, na.rm = T) - 5*mad(datPQCT$tib_38_crt_thk_c, na.rm = T)
sum(datPQCT$tib_38_crt_thk_c > cutoffU, na.rm = T)
sum(datPQCT$tib_38_crt_thk_c < cutoffL, na.rm = T)

# plot 
ggplot(datPQCT, aes(x = tib_38_crt_thk_c)) + 
  geom_histogram(colour="black", fill="white", binwidth = 0.15) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 60)) + 
  ggtitle("PMMST: Tibia 38% Cortical Thickness (mm)") +
  xlab("Tibia 38% Cortical Thickness (mm)")


#### Tibia 38%: Axial area moment of inertia ####

# summary statistics and histogram
summary(datPQCT$tib_38_i_circ)
hist(datPQCT$tib_38_i_circ)

# check for outliers using median +/- 5 mad
cutoffU <- median(datPQCT$tib_38_i_circ, na.rm = T) + 5*mad(datPQCT$tib_38_i_circ, na.rm = T)
cutoffL <- median(datPQCT$tib_38_i_circ, na.rm = T) - 5*mad(datPQCT$tib_38_i_circ, na.rm = T)
sum(datPQCT$tib_38_i_circ > cutoffU, na.rm = T)
sum(datPQCT$tib_38_i_circ < cutoffL, na.rm = T)

# plot 
ggplot(datPQCT, aes(x = tib_38_i_circ)) + 
  geom_histogram(colour="black", fill="white", binwidth = 550) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 60)) + 
  ggtitle(expression(paste("PMMST: Tibia 38% Axial area moment of inertia (", mm^{4},")"))) +
  xlab(expression(paste("Tibia 38% Axial area moment of inertia (", mm^{4},")")))


#### Tibia 38%: Tibia area ####

# summary statistics and histogram
summary(datPQCT$tib_38_tot_a)
hist(datPQCT$tib_38_tot_a)

# check for outliers using median +/- 5 mad
cutoffU <- median(datPQCT$tib_38_tot_a, na.rm = T) + 5*mad(datPQCT$tib_38_tot_a, na.rm = T)
cutoffL <- median(datPQCT$tib_38_tot_a, na.rm = T) - 5*mad(datPQCT$tib_38_tot_a, na.rm = T)
sum(datPQCT$tib_38_tot_a > cutoffU, na.rm = T)
sum(datPQCT$tib_38_tot_a < cutoffL, na.rm = T)

# plot 
ggplot(datPQCT, aes(x = tib_38_tot_a)) + 
  geom_histogram(colour="black", fill="white", binwidth = 15) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 60)) + 
  ggtitle(expression(paste("PMMST: Tibia 38% Area (", mm^{2},")"))) +
  xlab(expression(paste("Tibia 38% Area (", mm^{2},")")))


#### Tibia 50%: Cortical content ####

# summary statistics and histogram
summary(datPQCT$tib_50_crt_cnt)
hist(datPQCT$tib_50_crt_cnt)

# check for outliers using median +/- 5 mad
cutoffU <- median(datPQCT$tib_50_crt_cnt, na.rm = T) + 5*mad(datPQCT$tib_50_crt_cnt, na.rm = T)
cutoffL <- median(datPQCT$tib_50_crt_cnt, na.rm = T) - 5*mad(datPQCT$tib_50_crt_cnt, na.rm = T)
sum(datPQCT$tib_50_crt_cnt > cutoffU, na.rm = T)
sum(datPQCT$tib_50_crt_cnt < cutoffL, na.rm = T)

# plot 
ggplot(datPQCT, aes(x = tib_50_crt_cnt)) + 
  geom_histogram(colour="black", fill="white", binwidth = 9) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 60)) + 
  ggtitle("PMMST: Tibia 50% Cortical Content (mg/mm)") +
  xlab("Tibia 50% Cortical Content (mg/mm)")


#### Tibia 50%: Cortical density ####

# summary statistics and histogram
summary(datPQCT$tib_50_crt_den)
hist(datPQCT$tib_50_crt_den)

# check for outliers using median +/- 5 mad
cutoffU <- median(datPQCT$tib_50_crt_den, na.rm = T) + 5*mad(datPQCT$tib_50_crt_den, na.rm = T)
cutoffL <- median(datPQCT$tib_50_crt_den, na.rm = T) - 5*mad(datPQCT$tib_50_crt_den, na.rm = T)
sum(datPQCT$tib_50_crt_den > cutoffU, na.rm = T)
sum(datPQCT$tib_50_crt_den < cutoffL, na.rm = T)

# plot 
ggplot(datPQCT, aes(x = tib_50_crt_den)) + 
  geom_histogram(colour="black", fill="white", binwidth = 15) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 60)) + 
  ggtitle(expression(paste("PMMST: Tibia 50% Cortical Density (", mg/cm^{3},")"))) +
  xlab(expression(paste("Tibia 50% Cortical Density (", mg/cm^{3},")")))

#### Tibia 50%: Cortical area ####

# summary statistics and histogram
summary(datPQCT$tib_50_crt_a)
hist(datPQCT$tib_50_crt_a)

# check for outliers using median +/- 5 mad
cutoffU <- median(datPQCT$tib_50_crt_a, na.rm = T) + 5*mad(datPQCT$tib_50_crt_a, na.rm = T)
cutoffL <- median(datPQCT$tib_50_crt_a, na.rm = T) - 5*mad(datPQCT$tib_50_crt_a, na.rm = T)
sum(datPQCT$tib_50_crt_a > cutoffU, na.rm = T)
sum(datPQCT$tib_50_crt_a < cutoffL, na.rm = T)

# plot 
ggplot(datPQCT, aes(x = tib_50_crt_a)) + 
  geom_histogram(colour="black", fill="white", binwidth = 9) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 60)) + 
  ggtitle(expression(paste("PMMST: Tibia 50% Cortical Area (", mm^{2},")"))) +
  xlab(expression(paste("Tibia 50% Cortical Area (", mm^{2},")")))

#### Tibia 50%: Cortical thickness ####

# summary statistics and histogram
summary(datPQCT$tib_50_crt_thk_c)
hist(datPQCT$tib_50_crt_thk_c)

# check for outliers using median +/- 5 mad
cutoffU <- median(datPQCT$tib_50_crt_thk_c, na.rm = T) + 5*mad(datPQCT$tib_50_crt_thk_c, na.rm = T)
cutoffL <- median(datPQCT$tib_50_crt_thk_c, na.rm = T) - 5*mad(datPQCT$tib_50_crt_thk_c, na.rm = T)
sum(datPQCT$tib_50_crt_thk_c > cutoffU, na.rm = T)
sum(datPQCT$tib_50_crt_thk_c < cutoffL, na.rm = T)

# plot 
ggplot(datPQCT, aes(x = tib_50_crt_thk_c)) + 
  geom_histogram(colour="black", fill="white", binwidth = 0.15) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 60)) + 
  ggtitle("PMMST: Tibia 50% Cortical Thickness (mm)") +
  xlab("Tibia 50% Cortical Thickness (mm)")


#### Tibia 50%: Axial area moment of inertia ####

# summary statistics and histogram
summary(datPQCT$tib_50_i_circ)
hist(datPQCT$tib_50_i_circ)

# check for outliers using median +/- 5 mad
cutoffU <- median(datPQCT$tib_50_i_circ, na.rm = T) + 5*mad(datPQCT$tib_50_i_circ, na.rm = T)
cutoffL <- median(datPQCT$tib_50_i_circ, na.rm = T) - 5*mad(datPQCT$tib_50_i_circ, na.rm = T)
sum(datPQCT$tib_50_i_circ > cutoffU, na.rm = T)
sum(datPQCT$tib_50_i_circ < cutoffL, na.rm = T)

# plot 
ggplot(datPQCT, aes(x = tib_50_i_circ)) + 
  geom_histogram(colour="black", fill="white", binwidth = 550) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 60)) + 
  ggtitle(expression(paste("PMMST: Tibia 50% Axial area moment of inertia (", mm^{4},")"))) +
  xlab(expression(paste("Tibia 50% Axial area moment of inertia (", mm^{4},")")))


#### Tibia 50%: Tibia area ####

# summary statistics and histogram
summary(datPQCT$tib_50_tot_a)
hist(datPQCT$tib_50_tot_a)

# check for outliers using median +/- 5 mad
cutoffU <- median(datPQCT$tib_50_tot_a, na.rm = T) + 5*mad(datPQCT$tib_50_tot_a, na.rm = T)
cutoffL <- median(datPQCT$tib_50_tot_a, na.rm = T) - 5*mad(datPQCT$tib_50_tot_a, na.rm = T)
sum(datPQCT$tib_50_tot_a > cutoffU, na.rm = T)
sum(datPQCT$tib_50_tot_a < cutoffL, na.rm = T)

# plot 
ggplot(datPQCT, aes(x = tib_50_tot_a)) + 
  geom_histogram(colour="black", fill="white", binwidth = 15) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 60)) + 
  ggtitle(expression(paste("PMMST: Tibia 50% Area (", mm^{2},")"))) +
  xlab(expression(paste("Tibia 50% Area (", mm^{2},")")))


#### Tibia 66%: Muscle + bone area ####

summary(datPQCT$mu.ba)
hist(datPQCT$mu.ba)

# check for outliers using median +/- 5 mad
cutoffU <- median(datPQCT$mu.ba, na.rm = T) + 5*mad(datPQCT$mu.ba, na.rm = T)
cutoffL <- median(datPQCT$mu.ba, na.rm = T) - 5*mad(datPQCT$mu.ba, na.rm = T)
sum(datPQCT$mu.ba > cutoffU, na.rm = T)
sum(datPQCT$mu.ba < cutoffL, na.rm = T)

# plot 
ggplot(datPQCT, aes(x = mu.ba)) + 
  geom_histogram(colour="black", fill="white", binwidth = 250) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 60)) + 
  ggtitle(expression(paste("PMMST: Tibia 66% Muscle + Bone Area (", mm^{2},")"))) +
  xlab(expression(paste("Tibia 66% Muscle + Bone Area (", mm^{2},")")))

#### Tibia 66%: Bone area ####

summary(datPQCT$bone_csa)
hist(datPQCT$bone_csa)

# check for outliers using median +/- 5 mad
cutoffU <- median(datPQCT$bone_csa, na.rm = T) + 5*mad(datPQCT$bone_csa, na.rm = T)
cutoffL <- median(datPQCT$bone_csa, na.rm = T) - 5*mad(datPQCT$bone_csa, na.rm = T)
sum(datPQCT$bone_csa > cutoffU, na.rm = T)
sum(datPQCT$bone_csa < cutoffL, na.rm = T)

# plot 
ggplot(datPQCT, aes(x = bone_csa)) + 
  geom_histogram(colour="black", fill="white", binwidth = 20) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 60)) + 
  ggtitle(expression(paste("PMMST: Tibia 66% Bone Area (", mm^{2},")"))) +
  xlab(expression(paste("Tibia 66% Bone Area (", mm^{2},")")))

#### Tibia 66%: Total area ####

summary(datPQCT$total_area)
hist(datPQCT$total_area)

# check for outliers using median +/- 5 mad
cutoffU <- median(datPQCT$total_area, na.rm = T) + 5*mad(datPQCT$total_area, na.rm = T)
cutoffL <- median(datPQCT$total_area, na.rm = T) - 5*mad(datPQCT$total_area, na.rm = T)
sum(datPQCT$total_area > cutoffU, na.rm = T)
sum(datPQCT$total_area < cutoffL, na.rm = T)
# one outlier. Check the outlier
temp <- subset(datPQCT, datPQCT$total_area > cutoffU)

# flag as outliers
datPQCT$flagtotal_area <- ifelse(datPQCT$total_area > cutoffU, 1, 0)
datPQCT$flagtotal_area <- factor(datPQCT$flagtotal_area, labels = c("No Outlier", "Outlier"))
attributes(datPQCT$flagtotal_area)$label <- "pmmst - tibia 66% total area outlier?"
table(datPQCT$flagtotal_area)

# summary and plot
# including the outliers
summary(datPQCT$total_area)
ggplot(datPQCT, aes(x = total_area)) + 
  geom_histogram(colour="black", fill="white", binwidth = 250) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 60)) + 
  ggtitle(expression(paste("PMMST: Tibia 66% Total Area (", mm^{2},")"))) +
  xlab(expression(paste("Tibia 66% Total Area (", mm^{2},")")))

# excluding the outliers
datPQCT %>% filter(flagtotal_area == "No Outlier") %>% 
  summarise(n = n(), Min = min(total_area), Mean = mean(total_area), SD = 
              sd(total_area), Median = median(total_area), Q1 = 
              quantile(total_area, probs = 0.25), Q3 = 
              quantile(total_area, probs = 0.75))
ggplot(subset(datPQCT, flagtotal_area == "No Outlier"), aes(x = total_area)) + 
  geom_histogram(colour="black", fill="white", binwidth = 250) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 60))  + 
  ggtitle(expression(paste("PMMST (No Outliers): Tibia 66% Total Area (", mm^{2},")"))) +
  xlab(expression(paste("Tibia 66% Total Area (", mm^{2},")")))



#### Check trends/missing/scanner ####

# where were the scanners made?
addmargins(table(datPQCT$scanner))
# 2 did not have a scanner? need to check!

# check the two children with no scan whereabouts
temp <- datPQCT[is.na(datPQCT$scanner), ]
temp$patid
rm(temp)
# EMP101D and EMP276K not have any PQCT data

# check those who had missing tibia 8%
datPQCT[is.na(datPQCT$tot_a_8), ]
# do have the data at other sides of the tibia
datPQCT[is.na(datPQCT$tib_38_crt_a), ] # only the two missing children
datPQCT[is.na(datPQCT$tib_50_crt_a), ] # the two missing children plus 2 more.


#########################################################################################################
#### Save dataset #######################################################################################
#########################################################################################################


# save in csv, stata and SPSS (version 1 contained DXA data only)
# write.table(datBone, "./Cleaned Dataset/PMMSTBoneDXAv1_20180202.csv", row.names = FALSE)
# write_dta(datBone, "./Cleaned Dataset/PMMSTBoneDXAv1_20180202.dta")
# write_sav(datBone, "./Cleaned Dataset/PMMSTBoneDXAv1_20180202.sav")

# version 2 has both DXA and pQCT

datPQCT <- datPQCT %>% 
  dplyr::select(patid, scanner, tot_cnt, tot_den, crtsub_cnt, crtsub_den, trab_cnt, 
                trab_den, tot_a_8, trab_a, crtsub_a, tib_38_crt_cnt, tib_38_crt_den, 
                tib_38_crt_a, tib_38_crt_thk_c, tib_38_i_circ, tib_38_tot_a, tib_50_crt_cnt,
                tib_50_crt_den, tib_50_crt_a, tib_50_crt_thk_c, tib_50_i_circ, tib_50_tot_a, mu.ba, bone_csa,
                total_area, flagtotal_area)

names(datPQCT) <- c("serno", "cscannerPQCT", "ctotbonecon_8", "ctotboneden_8", "ccorticalcon_8", 
                    "ccorticalden_8", "ctrabecularcon_8", "ctrabecularden_8", "ctibiaarea_8",
                    "ctrabeculararea_8", "ccorticalarea_8", "ccorticalcon_38", "ccorticalden_38", 
                    "ccorticalarea_38", "ccorticalthick_38", "cinertia_38", "ctibiaarea_38",
                    "ccorticalcon_50", "ccorticalden_50", 
                    "ccorticalarea_50", "ccorticalthick_50", "cinertia_50", "ctibiaarea_50", 
                    "cmusclebone_66",
                    "cbone_66", "ctotalarea_66", "flagctotalarea_66")
head(datPQCT, 4)

dat <- left_join(datPQCT, datBone)
attributes(dat$serno)$label <- "pmmst - Subject ID"
attributes(dat$cvisit)$label <- "pmmst - clinic visit"
# add other attributes
attributes(dat$cscannerPQCT)$label <- "pmmst - location of pQCT scan"
attributes(dat$ctotbonecon_8)$label <- 
  "pmmst - tibia 8%: total bone content (mg/mm)"
attributes(dat$ctotboneden_8)$label <- 
  "pmmst - tibia 8%: total bone density (mg/cm^3)"
attributes(dat$ccorticalcon_8)$label <- 
  "pmmst - tibia 8%: cortical and subcortical content (mg/mm)"
attributes(dat$ccorticalden_8)$label <- 
  "pmmst - tibia 8%: cortical and subcortical density (mg/cm^3)"
attributes(dat$ctrabecularcon_8)$label <- "pmmst - tibia 8%: trabecular content (mg/mm)"
attributes(dat$ctrabecularden_8)$label <- "pmmst - tibia 8%: trabecular density (mg/cm^3)"
attributes(dat$ctibiaarea_8)$label <- "pmmst - tibia 8%: tibia area (mm^2)"
attributes(dat$ctrabeculararea_8)$label <- "pmmst - tibia 8%: trabecular area (mm^2)"
attributes(dat$ccorticalarea_8)$label <- 
  "pmmst - tibia 8%: cortical and subcortical area (mm^2)"
attributes(dat$ccorticalcon_38)$label <- 
  "pmmst - tibia 38%: cortical and subcortical content (mm^2)"
attributes(dat$ccorticalden_38)$label <- 
  "pmmst - tibia 38%: cortical and subcortical denisty (mg/mm^3)"
attributes(dat$ccorticalarea_38)$label <- 
  "pmmst - tibia 38%: cortical and subcortical area (mm^2)"
attributes(dat$ccorticalthick_38)$label <- "pmmst - tibia 38%: cortical thickness (mm)"
attributes(dat$cinertia_38)$label <- "pmmst - tibia 38%: axial moment of inertia (mm^4)"
attributes(dat$ctibiaarea_38)$label <- "pmmst - tibia 38%: tibia area (mm^2)"
attributes(dat$ccorticalcon_50)$label <- 
  "pmmst - tibia 50%: cortical and subcortical content (mm^2)"
attributes(dat$ccorticalden_50)$label <- 
  "pmmst - tibia 50%: cortical and subcortical denisty (mg/mm^3)"
attributes(dat$ccorticalarea_50)$label <- 
  "pmmst - tibia 50%: cortical and subcortical area (mm^2)"
attributes(dat$ccorticalthick_50)$label <- 
  "pmmst - tibia 50%: cortical thickness (mm)"
attributes(dat$cinertia_50)$label <- 
  "pmmst - tibia 50%: axial moment of inertia (mm^4)"
attributes(dat$ctibiaarea_50)$label <- "pmmst - tibia 50%: tibia area (mm^2)"
attributes(dat$cmusclebone_66)$label <- "pmmst - tibia 66%: muscle + bone area (mm^2)"
attributes(dat$cbone_66)$label <- "pmmst - tibia 66%: bone area (mm^2)"
attributes(dat$ctotalarea_66)$label <- "pmmst - tibia 66%: total area (mm^2)"
attributes(dat$flagctotalarea_66)$label <- "pmmst - tibia 66%: is total area outlier?"

# save in csv, stata and SPSS (version 1 contained DXA data only)
write.table(dat, "./Cleaned Dataset/PMMSTBoneDXApQCTv2_20180628.csv", row.names = FALSE)
write_dta(dat, "./Cleaned Dataset/PMMSTBoneDXApQCTv2_20180628.dta")
write_sav(dat, "./Cleaned Dataset/PMMSTBoneDXApQCTv2_20180628.sav")
