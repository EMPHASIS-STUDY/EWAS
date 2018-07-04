##############################################
##### SARAS Kids: Phenotype Data DXA #########
##############################################

# 16/06/2017

# load packages
library(haven)
library(dplyr)
library(ggplot2)

# Data were sent by Vanessa on 08/06/2017. 
# File name is sarasdxa04.dta (STATA v14 file). 
# Only children in the per-protocol analysis are processed.

datDXA <- read_dta("sarasdxa04.dta")
# check for duplicates
length(unique(datDXA$serno)) # there no duplicates
# remove sex from DXA data as sex was already checked in SARASKidsBloodCountPP.dta
datDXA <- subset(datDXA, select = -sex)

# add information on allocation group, sex and age. PP analysis only
datList <- read_dta("SARASKidsBloodCountPP.dta")
datList <- datList %>% dplyr::select(serno, AnalysisType, allocation, age, sex)
# check if all sernos in datList are in dat
length(datList$serno[!(datList$serno %in% datDXA$serno)]) 
# 3 sernos do not have DXA data
datList$serno[!(datList$serno %in% datDXA$serno)]
# 1762 (Down syndrome), 2587, 2785 did not have DXA data 
datDXA <- inner_join(datDXA, datList)
# remove dataset not longer useful
rm(datList)


# sanity check
str(datDXA)
dim(datDXA)

#############################
#### Analysis Steps #########
#############################

# 1. Check each variable singularly
# 2. Check each variable in relation with other relevant variables collected.
# 3. Check for outliers using the average computed in 3. and the criteria established in the analysis plan (median +/- 5MAD)
# 4. Retain the outliers but create a variable with 0/1 (named flag_variable name) where 1 indicated
#    the outlier based on the analysis plan.


########################################################################################################
#### DXA BODY COMPOSITION ############################################################################## 
########################################################################################################

# total fat mass 
summary(datDXA$totfatmass)
hist(datDXA$totfatmass, main = "Total fat mass (g)")
# no one is missing total fat mass values

# check for outliers using median +/- 5 mad
cutoffU <- median(datDXA$totfatmass) + 5*mad(datDXA$totfatmass)
cutoffL <- median(datDXA$totfatmass) - 5*mad(datDXA$totfatmass)
sum(datDXA$totfatmass > cutoffU)
sum(datDXA$totfatmass < cutoffL)
# 9 outliers as also identified by histogram
# check the outliers before removal
datDXA[datDXA$totfatmass > cutoffU , c("serno", "ht", "wt", "totfatmass", "totleanmass", "totbonemass")]
# graph fat mass vs weight and highlight the outliers
ggplot(datDXA, aes(x = wt, y = totfatmass)) + geom_point() +
  scale_x_continuous(breaks = round(seq(min(datDXA$wt), max(datDXA$wt), by = 5), 0)) +
  scale_y_continuous(breaks = round(seq(700, 15000, by = 2000), 0)) +
  geom_point(data = datDXA[datDXA$totfatmass > cutoffU,], aes(x = wt, y = totfatmass), 
             colour = "red") + xlab("Weight (kg)") + ylab("Total Fat Mass (g)") 

# create a flag variable for fat mass
datDXA$flagfatmass <- ifelse(datDXA$totfatmass > cutoffU, 1, 0)
datDXA$flagfatmass <- factor(datDXA$flagfatmass, labels = c("No Outlier", "Outlier"))
attributes(datDXA$flagfatmass)$label <- "saras - fat mass outlier?"

# transform fat mass in kg
datDXA$totfatmass_kg <- datDXA$totfatmass/1000
attributes(datDXA$totfatmass_kg)$label <- "saras - fat mass (kg)"

# summary and plot
# including the outliers
summary(datDXA$totfatmass_kg)
ggplot(datDXA, aes(x = totfatmass_kg)) + geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 200)) + ggtitle("SARAS Kids: Total Fat Mass (kg)") +
  xlab("Total Fat Mass (kg)")
# excluding the outliers
datDXA %>% filter(flagfatmass == "No Outlier") %>% 
  summarise(n = n(), Min = min(totfatmass_kg), Mean = mean(totfatmass_kg), SD = 
              sd(totfatmass_kg), Median = median(totfatmass_kg), Q1 = 
              quantile(totfatmass_kg, probs = 0.25), Q3 = 
              quantile(totfatmass_kg, probs = 0.75))
ggplot(subset(datDXA, flagfatmass == "No Outlier"), aes(x = totfatmass_kg)) + 
  geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 80)) + 
  ggtitle("SARAS Kids (No Outliers): Total Fat Mass (kg)") +
  xlab("Total Fat Mass (kg)")


######################################################################################################

# total lean mass 
summary(datDXA$totleanmass)
hist(datDXA$totleanmass, main = "Total lean mass (g)")
# no one is missing total lean mass values

# check for outliers using median +/- 5 mad
cutoffU <- median(datDXA$totleanmass) + 5*mad(datDXA$totleanmass)
cutoffL <- median(datDXA$totleanmass) - 5*mad(datDXA$totleanmass)
sum(datDXA$totleanmass > cutoffU)
sum(datDXA$totleanmass < cutoffL)
# one outlier as also identified by histogram
# check the outlier before removal
datDXA[datDXA$totleanmass > cutoffU , c("serno", "ht", "wt", "totfatmass", "totleanmass", "totbonemass")]
# graph lean mass vs weight and highlight the outliers
ggplot(datDXA, aes(x = wt, y = totleanmass)) + geom_point() +
  scale_x_continuous(breaks = round(seq(min(datDXA$wt), max(datDXA$wt), by = 5), 0)) +
  scale_y_continuous(breaks = round(seq(8000, 22000, by = 2000), 0)) +
  geom_point(data = datDXA[datDXA$totleanmass > cutoffU,], aes(x = wt, y = totleanmass), 
             colour = "red") + xlab("Weight (kg)") + ylab("Total Lean Mass (g)") 

# create a flag variable for lean mass
datDXA$flagleanmass <- ifelse(datDXA$totleanmass > cutoffU, 1, 0)
datDXA$flagleanmass <- factor(datDXA$flagleanmass, labels = c("No Outlier", "Outlier"))
attributes(datDXA$flagleanmass)$label <- "saras - lean mass outlier?"

# transform lean mass in kg
datDXA$totleanmass_kg <- datDXA$totleanmass/1000
attributes(datDXA$totleanmass_kg)$label <- "saras - lean mass (kg)"

# summary and plot
# including the outliers
summary(datDXA$totleanmass_kg)
ggplot(datDXA, aes(x = totleanmass_kg)) + geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 200)) + ggtitle("SARAS Kids: Total Lean Mass (kg)") +
  xlab("Total Lean Mass (kg)")
# excluding the outliers
datDXA %>% filter(flagleanmass == "No Outlier") %>% 
  summarise(n = n(), Min = min(totleanmass_kg), Mean = mean(totleanmass_kg), SD = 
              sd(totleanmass_kg), Median = median(totleanmass_kg), Q1 = 
              quantile(totleanmass_kg, probs = 0.25), Q3 = 
              quantile(totleanmass_kg, probs = 0.75))
ggplot(subset(datDXA, flagleanmass == "No Outlier"), aes(x = totleanmass_kg)) + 
  geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 80)) + 
  ggtitle("SARAS Kids (No Outliers): Total Lean Mass (kg)") +
  xlab("Total Lean Mass (kg)")

########################################################################################################

# android fat 
summary(datDXA$androidfatmass)
hist(datDXA$androidfatmass, main = "Total android fat (g)")
# no one is missing android fat

# check for outliers using median +/- 5 mad
cutoffU <- median(datDXA$androidfatmass) + 5*mad(datDXA$androidfatmass)
cutoffL <- median(datDXA$androidfatmass) - 5*mad(datDXA$androidfatmass)
sum(datDXA$androidfatmass > cutoffU)
sum(datDXA$androidfatmass < cutoffL)
# 12 as also identified by histogram
# check the outlier before removal
datDXA[datDXA$androidfatmass > cutoffU , c("serno", "ht", "wt", "totfatmass", "androidfatmass")]
ggplot(datDXA, aes(x = wt, y = androidfatmass)) + geom_point() +
  geom_point(data = datDXA[datDXA$androidfatmass > cutoffU,], aes(x = wt, y = androidfatmass), 
             colour = "red") + xlab("Weight (kg)") + ylab("Total Android fat (g)") 

# create a flag variable for android fat
datDXA$flagandroidfatmass <- ifelse(datDXA$androidfatmass > cutoffU, 1, 0)
datDXA$flagandroidfatmass <- factor(datDXA$flagandroidfatmass, labels = c("No Outlier", "Outlier"))
attributes(datDXA$flagandroidfatmass)$label <- "saras - android fat outlier?"

# transform android fat in kg
datDXA$androidfatmass_kg <- datDXA$androidfatmass/1000
attributes(datDXA$androidfatmass_kg)$label <- "android fat mass (kg)"


# summary and plot
# including the outliers
summary(datDXA$androidfatmass_kg)
ggplot(datDXA, aes(x = androidfatmass_kg)) + geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 200)) + ggtitle("SARAS Kids: Total Android Fat (kg)") +
  xlab("Android Fat (kg)")
# excluding the outliers
datDXA %>% filter(flagandroidfatmass == "No Outlier") %>% 
  summarise(n = n(), Min = min(androidfatmass_kg), Mean = mean(androidfatmass_kg), SD = 
              sd(androidfatmass_kg), Median = median(androidfatmass_kg), Q1 = 
              quantile(androidfatmass_kg, probs = 0.25), Q3 = 
              quantile(androidfatmass_kg, probs = 0.75))
ggplot(subset(datDXA, flagandroidfatmass == "No Outlier"), aes(x = androidfatmass_kg)) + 
  geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 85)) + 
  ggtitle("SARAS Kids (No Outliers): Android Fat (kg)") +
  xlab("Android Fat (kg)")


########################################################################################################

# gynoid fat 
summary(datDXA$gynoidfatmass)
hist(datDXA$gynoidfatmass, main = "Total gynoid fat (g)")
# no one is missing android fat

# check for outliers using median +/- 5 mad
cutoffU <- median(datDXA$gynoidfatmass) + 5*mad(datDXA$gynoidfatmass)
cutoffL <- median(datDXA$gynoidfatmass) - 5*mad(datDXA$gynoidfatmass)
sum(datDXA$gynoidfatmass > cutoffU)
sum(datDXA$gynoidfatmass < cutoffL)
# 5 as also identified by histogram
# check the outlier before removal
datDXA[datDXA$gynoidfatmass > cutoffU , c("serno", "ht", "wt", "totfatmass", "gynoidfatmass")]
ggplot(datDXA, aes(x = wt, y = gynoidfatmass)) + geom_point() +
  geom_point(data = datDXA[datDXA$gynoidfatmass > cutoffU,], aes(x = wt, y = gynoidfatmass), 
             colour = "red") + xlab("Weight (kg)") + ylab("Total Gynoid fat (g)") 

# create a flag variable for gynoid fat
datDXA$flaggynoidfatmass <- ifelse(datDXA$gynoidfatmass > cutoffU, 1, 0)
datDXA$flaggynoidfatmass <- factor(datDXA$flaggynoidfatmass, labels = c("No Outlier", "Outlier"))
attributes(datDXA$flaggynoidfatmass)$label <- "saras - gynoid fat outlier?"

# transform gynoid fat in kg
datDXA$gynoidfatmass_kg <- datDXA$gynoidfatmass/1000
attributes(datDXA$gynoidfatmass_kg)$label <- "saras - gynoid fat mass (kg)"

# summary and plot
# including the outliers
summary(datDXA$gynoidfatmass_kg)
ggplot(datDXA, aes(x = gynoidfatmass_kg)) + geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 200)) + ggtitle("SARAS Kids: Total Gynoid Fat (kg)") +
  xlab("Gynoid Fat (kg)")
# excluding the outliers
datDXA %>% filter(flaggynoidfatmass == "No Outlier") %>% 
  summarise(n = n(), Min = min(gynoidfatmass_kg), Mean = mean(gynoidfatmass_kg), SD = 
              sd(gynoidfatmass_kg), Median = median(gynoidfatmass_kg), Q1 = 
              quantile(gynoidfatmass_kg, probs = 0.25), Q3 = 
              quantile(gynoidfatmass_kg, probs = 0.75))
ggplot(subset(datDXA, flaggynoidfatmass == "No Outlier"), aes(x = gynoidfatmass_kg)) + 
  geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 85)) + 
  ggtitle("SARAS Kids (No Outliers): Gynoid Fat (kg)") +
  xlab("Gynoid Fat (kg)")



########################################################################################################

# fat mass index (FMI): fat mass (kg)/height (m)^2
summary(datDXA$ht)
hist(datDXA$ht)
# check height values that are too high
datDXA[datDXA$ht > 150 , c("serno", "ht", "wt")]
# serno 780 height was recorded wrong. Real height is 109.8 
datDXA[datDXA$serno == 780, "ht"] <- 109.8
# height is in cm, need to transform in m. Can do directly in the FMI formula
datDXA$fmi <- datDXA$totfatmass_kg/((0.01*datDXA$ht)^2)
attributes(datDXA$fmi)$label <- "saras - fat mass index (kg/m2)" 

# summary and plot (no outliers removed)
summary(datDXA$fmi)
ggplot(datDXA, aes(x = fmi)) + geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 180)) + 
  ggtitle(expression(paste("SARAS Kids: Fat Mass Index (", kg/m^{2},")"))) +
  xlab(expression(paste("Fat Mass Index (", kg/m^{2}, ")")))
datDXA[datDXA$fmi > cutoffU,]

# check for outliers
cutoffU <- median(datDXA$fmi) + 5*mad(datDXA$fmi)
cutoffL <- median(datDXA$fmi) - 5*mad(datDXA$fmi)
sum(datDXA$fmi > cutoffU)
sum(datDXA$fmi < cutoffL)
# 6 outliers identified. Check those six
ggplot(datDXA, aes(x = ht, y = totfatmass_kg)) + geom_point() +
  scale_x_continuous(breaks = round(seq(min(datDXA$ht), max(datDXA$ht), by = 5), 0)) +
  geom_point(data = datDXA[datDXA$fmi > cutoffU,], aes(x = ht, y = totfatmass_kg), 
             colour = "red") + xlab("Height (cm)") + ylab("Fat Mass (kg)") 

# create a flag variable for fmi
datDXA$flagfmi <- ifelse(datDXA$fmi > cutoffU, 1, 0)
datDXA$flagfmi <- factor(datDXA$flagfmi, labels = c("No Outlier", "Outlier"))
attributes(datDXA$flagfmi)$label <- "saras - fat mass index outlier?"

# excluding the outliers
datDXA %>% filter(flagfmi == "No Outlier") %>% 
  summarise(n = n(), Min = min(fmi), Mean = mean(fmi), SD = 
              sd(fmi), Median = median(fmi), Q1 = 
              quantile(fmi, probs = 0.25), Q3 = 
              quantile(fmi, probs = 0.75))
ggplot(subset(datDXA, flagfmi == "No Outlier"), aes(x = fmi)) + 
  geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 85)) + 
  ggtitle(expression(paste("SARAS Kids (No Outliers): Fat Mass Index (", kg/m^{2},")"))) +
  xlab(expression(paste("Fat Mass Index (", kg/m^{2}, ")")))


# lean mass index (LMI): lean mass (kg)/height (m)^2
# height is in cm, need to transform in m. Can do directly in the FMI formula
datDXA$lmi <- datDXA$totleanmass_kg/((0.01*datDXA$ht)^2)
attributes(datDXA$lmi)$label <- "saras - lean mass index (kg/m2)"

# summary and plot (including outliers)
summary(datDXA$lmi)
ggplot(datDXA, aes(x = lmi)) + geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 90)) + 
  ggtitle(expression(paste("SARAS Kids: Lean Mass Index (", kg/m^{2},")"))) +
  xlab(expression(paste("Lean Mass Index (", kg/m^{2}, ")")))

# check for outliers
cutoffU <- median(datDXA$lmi) + 5*mad(datDXA$lmi)
cutoffL <- median(datDXA$lmi) - 5*mad(datDXA$lmi)
sum(datDXA$lmi > cutoffU)
sum(datDXA$lmi < cutoffL)
# 1 outliers identified. 
datDXA[datDXA$lmi > cutoffU,]

# create a flag variable for fmi
datDXA$flaglmi <- ifelse(datDXA$lmi > cutoffU, 1, 0)
datDXA$flaglmi <- factor(datDXA$flaglmi, labels = c("No Outlier", "Outlier"))
attributes(datDXA$flaglmi)$label <- "saras - lean mass index outlier?"


# excluding the outliers
datDXA %>% filter(flaglmi == "No Outlier") %>% 
  summarise(n = n(), Min = min(lmi), Mean = mean(lmi), SD = 
              sd(lmi), Median = median(lmi), Q1 = 
              quantile(lmi, probs = 0.25), Q3 = 
              quantile(lmi, probs = 0.75))
ggplot(subset(datDXA, flaglmi == "No Outlier"), aes(x = lmi)) + 
  geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 70)) + 
  ggtitle(expression(paste("SARAS Kids (No Outliers): Lean Mass Index (", kg/m^{2},")"))) +
  xlab(expression(paste("Lean Mass Index (", kg/m^{2}, ")")))


########################################################################################################


# Body Fat % ((Fat Mass x 100) / weight)
# check the variable weight
summary(datDXA$wt)
hist(datDXA$wt)
# children with high observations: weight is real. Those were heavier children.
# create a variable for body fat % (fatpr)
datDXA$fatper <- with(datDXA, 100*(totfatmass_kg/wt))
attributes(datDXA$fatper)$label <- "saras - body fat %"

# summary and plot (including outliers)
summary(datDXA$fatper)
ggplot(datDXA, aes(x = fatper)) + geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 95)) + 
  ggtitle("SARAS Kids: Body fat %") +
  xlab("Body fat %")

# check for outliers
cutoffU <- median(datDXA$fatper) + 5*mad(datDXA$fatper)
cutoffL <- median(datDXA$fatper) - 5*mad(datDXA$fatper)
sum(datDXA$fatper > cutoffU)
sum(datDXA$fatper < cutoffL)
# 2 outliers identified. Check those 
datDXA[datDXA$fatper > cutoffU,]

# create a flag variable for fmi
datDXA$flagfatper <- ifelse(datDXA$fatper > cutoffU, 1, 0)
datDXA$flagfatper <- factor(datDXA$flagfatper, labels = c("No Outlier", "Outlier"))
attributes(datDXA$flagfatper)$label <- "saras - body fat % outlier?"


# excluding the outliers
datDXA %>% filter(flagfatper == "No Outlier") %>% 
  summarise(n = n(), Min = min(fatper), Mean = mean(fatper), SD = 
              sd(fatper), Median = median(fatper), Q1 = 
              quantile(fatper, probs = 0.25), Q3 = 
              quantile(fatper, probs = 0.75), Max = max(fatper))
ggplot(subset(datDXA, flagfatper == "No Outlier"), aes(x = fatper)) + 
  geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 80)) + 
  ggtitle("SARAS Kids: Body fat %") +
  xlab("Body fat %")




########################################################################################################
#### DXA BONE MEASURES ################################################################################# 
########################################################################################################

# Total body bone mineral content (BMC, whole body minus the head) (g)
summary(datDXA$totbmc)
# Total BMC values also have head incorporated. Need to subtract head
datDXA$totbmcNohead <- with(datDXA, totbmc - headbmc)
summary(datDXA$totbmcNohead)
hist(datDXA$totbmcNohead)
attributes(datDXA$totbmcNohead)$label <- "saras - total BMC excluding head"

# summary and plot (no outlier removed)
summary(datDXA$totbmcNohead)
ggplot(datDXA, aes(x = totbmcNohead)) + geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 90)) + 
  ggtitle("SARAS Kids: Total body BMC excluding head (g)") +
  xlab("Total body BMC excluding head (g)")


# check for outliers in BMC using median +/- 5 mad
cutoffU <- median(datDXA$totbmcNohead) + 5*mad(datDXA$totbmcNohead)
cutoffL <- median(datDXA$totbmcNohead) - 5*mad(datDXA$totbmcNohead)
sum(datDXA$totbmcNohead > cutoffU)
sum(datDXA$totbmcNohead < cutoffL)
# 2 as also identified by histogram
# check the outlier
datDXA[datDXA$totbmcNohead > cutoffU , c("serno", "ht", "wt", "headbmc", "headbmd", "totbmc")]

# create a flag variable for BMC
datDXA$flagtotbmcNohead <- ifelse(datDXA$totbmcNohead > cutoffU, 1, 0)
datDXA$flagtotbmcNohead <- factor(datDXA$flagtotbmcNohead, labels = c("No Outlier", "Outlier"))
rm(cutoffU); rm(cutoffL)
attributes(datDXA$flagtotbmcNohead)$label <- "saras - total BMC outlier?"

# excluding the outliers
datDXA %>% filter(flagtotbmcNohead == "No Outlier") %>% 
  summarise(n = n(), Min = min(totbmcNohead), Mean = mean(totbmcNohead), SD = 
              sd(totbmcNohead), Median = median(totbmcNohead), Q1 = 
              quantile(totbmcNohead, probs = 0.25), Q3 = 
              quantile(totbmcNohead, probs = 0.75))
ggplot(subset(datDXA, flagtotbmcNohead == "No Outlier"), aes(x = totbmcNohead)) + 
  geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 90)) + 
  ggtitle("SARAS Kids (No Outliers): Total body BMC excluding head (g)") +
  xlab("Total body BMC excluding head (g)")




########################################################################################################

# Total body bone area (BA, whole body minus the head) (cm2)
summary(datDXA$totarea)
# Total BA values also have head incorporated. Need to subtract head
datDXA$totareaNohead <- with(datDXA, totarea - headarea)
attributes(datDXA$totareaNohead)$label <- "saras - total BA excluding head"
summary(datDXA$totareaNohead)
hist(datDXA$totareaNohead)

# summary and plot (no outliers removed)
summary(datDXA$totareaNohead)
ggplot(datDXA, aes(x = totareaNohead)) + geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 90)) + 
  ggtitle(expression(paste("SARAS Kids: Total body bone area (", cm^{2},")"))) +
  xlab(expression(paste("Total body bone area (", cm^{2},")")))

# check for outliers in BA using median +/- 5 mad
cutoffU <- median(datDXA$totareaNohead) + 5*mad(datDXA$totareaNohead)
cutoffL <- median(datDXA$totareaNohead) - 5*mad(datDXA$totareaNohead)
sum(datDXA$totareaNohead > cutoffU)
sum(datDXA$totareaNohead < cutoffL)
# 1 as also identified by histogram
# check the outlier before removal
datDXA[datDXA$totareaNohead > cutoffU , c("serno", "ht", "wt", "totarea", "headarea")]

# create a flag variable for BA
datDXA$flagtotareaNohead <- ifelse(datDXA$totareaNohead > cutoffU, 1, 0)
datDXA$flagtotareaNohead <- factor(datDXA$flagtotareaNohead, labels = c("No Outlier", "Outlier"))
rm(cutoffU); rm(cutoffL)
attributes(datDXA$flagtotareaNohead)$label <- "saras - total BA outlier?"

# excluding the outliers
datDXA %>% filter(flagtotareaNohead == "No Outlier") %>% 
  summarise(n = n(), Min = min(totareaNohead), Mean = mean(totareaNohead), SD = 
              sd(totareaNohead), Median = median(totareaNohead), Q1 = 
              quantile(totareaNohead, probs = 0.25), Q3 = 
              quantile(totareaNohead, probs = 0.75))
ggplot(subset(datDXA, flagtotareaNohead == "No Outlier"), aes(x = totareaNohead)) + 
  geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 90)) + 
  ggtitle(expression(paste("SARAS Kids (No Outliers): Total body bone area (", cm^{2},")"))) +
  xlab(expression(paste("Total body bone area (", cm^{2},")")))



########################################################################################################

# Total body bone mineral density (BMD, whole body minus the head) (g/cm^2)
summary(datDXA$totbmd)
# Total BMD values also have head incorporated. BMD is computed as BMC/BA. Recompute the 
# ratio considering BMC and BA without the head
datDXA$totbmdNohead <- with(datDXA, totbmcNohead/totareaNohead)

# summary and plots (no outlier removed)
summary(datDXA$totbmdNohead)
ggplot(datDXA, aes(x = totbmdNohead)) + geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 90)) + 
  ggtitle(expression(paste("SARAS Kids: Total body bone density (", g/cm^{2},")"))) +
  xlab(expression(paste("Total body bone density (", g/cm^{2},")")))


# check for outliers in BMD using median +/- 5 mad
cutoffU <- median(datDXA$totbmdNohead) + 5*mad(datDXA$totbmdNohead)
cutoffL <- median(datDXA$totbmdNohead) - 5*mad(datDXA$totbmdNohead)
sum(datDXA$totbmdNohead > cutoffU)
sum(datDXA$totbmdNohead < cutoffL)
# no additional outliers



########################################################################################################

# As they are spine data are unsuable in the Indian cohort

# Spine BMC, BMD and BA

# summary and plots (no outliers removal)
summary(datDXA$spbmc)
summary(datDXA$spbmd)
summary(datDXA$sparea)

ggplot(datDXA, aes(x = spbmc)) + geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 95)) + 
  ggtitle("SARAS Kids: Spine BMC (g)") +
  xlab("Spine BMC (g)")
ggplot(datDXA, aes(x = spbmd)) + geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 75)) + 
  ggtitle(expression(paste("SARAS Kids: Spine BMD (", g/cm^{2},")"))) +
  xlab(expression(paste("Spine BMD (", g/cm^{2},")")))
ggplot(datDXA, aes(x = sparea)) + geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 90)) + 
  ggtitle(expression(paste("SARAS Kids: Spine BA (", cm^{2},")"))) +
  xlab(expression(paste("Spine BA (", cm^{2},")")))


# check for outliers 
cutoffU <- median(datDXA$spbmc) + 5*mad(datDXA$spbmc)
cutoffL <- median(datDXA$spbmc) - 5*mad(datDXA$spbmc)
sum(datDXA$spbmc > cutoffU)
sum(datDXA$spbmc < cutoffL)
# 1 spine BMC outlier
# check the outlier before removal
datDXA[datDXA$spbmc > cutoffU , c("serno", "ht", "wt", "spbmc", "spbmd", "sparea")]

# create a flag variable for spine BMC
datDXA$flagspbmc <- ifelse(datDXA$spbmc > cutoffU, 1, 0)
datDXA$flagspbmc <- factor(datDXA$flagspbmc, labels = c("No Outlier", "Outlier"))
rm(cutoffU); rm(cutoffL)
attributes(datDXA$flagspbmc)$label <- "saras - spine BMC outlier?"


# excluding the outliers
# spine BMC
datDXA %>% filter(flagspbmc == "No Outlier") %>% 
  summarise(n = n(), Min = min(spbmc), Mean = mean(spbmc), SD = 
              sd(spbmc), Median = median(spbmc), Q1 = 
              quantile(spbmc, probs = 0.25), Q3 = 
              quantile(spbmc, probs = 0.75))
ggplot(subset(datDXA, flagspbmc == "No Outlier"), aes(x = spbmc)) + 
  geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 75))+ 
  ggtitle("SARAS Kids (No Outliers): Spine BMC (g)") +
  xlab("Spine BMC (g)")


cutoffU <- median(datDXA$spbmd) + 5*mad(datDXA$spbmd)
cutoffL <- median(datDXA$spbmd) - 5*mad(datDXA$spbmd)
sum(datDXA$spbmd > cutoffU)
sum(datDXA$spbmd < cutoffL)
# 0 spine BMD outliers
rm(cutoffU); rm(cutoffL)

cutoffU <- median(datDXA$sparea) + 5*mad(datDXA$sparea)
cutoffL <- median(datDXA$sparea) - 5*mad(datDXA$sparea)
sum(datDXA$sparea > cutoffU)
sum(datDXA$sparea < cutoffL)
# 1 spine BA outlier
# check the outlier before removal
datDXA[datDXA$sparea > cutoffU , c("serno", "ht", "wt", "spbmc", "spbmd", "sparea")]
# create a flag variable for spine BA
datDXA$flagsparea <- ifelse(datDXA$sparea > cutoffU, 1, 0)
datDXA$flagsparea <- factor(datDXA$flagsparea, labels = c("No Outlier", "Outlier"))
rm(cutoffU); rm(cutoffL)
attributes(datDXA$flagsparea)$label <- "saras - spine BA outlier?"

# excluding the outliers
# spine BA
datDXA %>% filter(flagsparea == "No Outlier") %>% 
  summarise(n = n(), Min = min(sparea), Mean = mean(sparea), SD = 
              sd(sparea), Median = median(sparea), Q1 = 
              quantile(sparea, probs = 0.25), Q3 = 
              quantile(sparea, probs = 0.75))
ggplot(subset(datDXA, flagsparea == "No Outlier"), aes(x = sparea)) + 
  geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 75)) + 
  ggtitle(expression(paste("SARAS Kids (No Outliers): Spine BA (", g/cm^{2},")"))) +
  xlab(expression(paste("Spine BA (", cm^{2},")")))

# we do not have individual vertebrae, hence bone mineral apparent density cannot be calculated as of
# now


########################################################################################################
#### DXA SAVE BONE DATA ################################################################################ 
########################################################################################################

names(datDXA)
datDXA <- datDXA %>% dplyr::select(serno, age, sex, allocation, totbmdNohead,
                                   totbmcNohead, totareaNohead,
                                   flagtotbmcNohead, flagtotareaNohead) %>% arrange(-desc(serno))
datDXA$sex <- as_factor(datDXA$sex)
datDXA$allocation <- as_factor(datDXA$allocation)

names(datDXA) <- c("serno", "cage", "csex", "allocation", "cBMD", "cBMC", "cBA", "flagcBMC", 
                   "flagcBA")

attributes(datDXA$serno)$label <- "saras - subject ID"

# save as csv as stata ans as spss file
write.table(datDXA, "./Cleaned Dataset/SARASKidsBoneDXAv1_20180302.csv", row.names = FALSE)
write_dta(datDXA, "./Cleaned Dataset/SARASKidsBoneDXAv1_20180302.dta")
write_sav(datDXA, "./Cleaned Dataset/SARASKidsBoneDXAv1_20180302.sav")
