###############################################################
##### SARAS Kids: Phenotype Data Processing Anthropometry #####
###############################################################

# 19/06/2017

# Dataset needed for the cleaning process:
# sarascp01.dta
# sarascp01extra.dta
# SARASKidsBloodCountPP.dta

# load packages
library(haven)
library(dplyr)
library(data.table)
library(ggplot2)

# Data were sent by Patsy on 05/06/2017. 
# File name is sarascp01.dta (STATA v14 file). 
# Only children in the per-protocol analysis are processed.
 
dat <- read_dta("sarascp01.dta")
attributes(dat$serno) <- NULL
# check for duplicates
length(unique(dat$serno))
dat <- subset(dat, select = -sex)

# add information on blood counts, allocation group, sex and age. PP analysis only
datList <- read_dta("SARASKidsBloodCountPP.dta")
datList <- datList %>% select(serno, allocation, DateTest, BloodTestDate, age, sex)
datList$sex <- factor(datList$sex, labels = c("Boy", "Girl"))
datList$allocation <- factor(datList$allocation, labels = c("B", "A"))
# check if all sernos in datList are in dat
datList$serno[!(datList$serno %in% dat$serno)] 
# 11 sernos are not in dat:
# 981, 1536, 1548, 1939, 2587, 4250, 4500, 4513, 4580, 4701, 5218

# Vanessa and Harshad looked for the 7 sernos and create an extra dataset with their 
# anthropometry measures

datextra <- read_dta("sarascp01extra.dta")
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

# Harshad sent the available data by email on 09 June 2017. Extra data are in document Anthro.xlsx.
# Add data manually
more.rows <- data.frame(serno = c(981, 1536, 1939, 2587), stringsAsFactors=F)
dat[(nrow(dat) + 1):(nrow(dat) + nrow(more.rows)), names(more.rows)] <- more.rows
rm(more.rows)
tail(dat)
# transform dat to data.table before adding extra information
dat <- data.table(dat)
# serno 981
dat[serno == 981, c("mweightc", "mstdhgtc", "msithgtc", "mmuac1c", "mmuac2c", "mmuac3c", "mhdcir1c", 
                   "mhdcir2c", "mhdcir3c", "mchest1c", "mchest2c", "mchest3c", "mabdo1c", "mabdo2c", 
                   "mabdo3c", "mtricp1c", "mtricp2c", "mtricp3c", "mbicep1c", "mbicep2c", "mbicep3c",     
                   "msubscp1c", "msubscp2c", "msubscp3c", "msupri1c", "msupri2c", "msupri3c")  := 
     list(16.7, 110.9, 83.6, 15.5, 16, 15.8, 49.4, 49, 50.9, 51, 51.2, 51.6, 52.2, 52.6, 52.6, 8.2, 8.4,
          8.2, 5.8, 5.6, 5.6, 6.2, 6.4, 6.8, 5.4, 5.4, 5.6), ]
# serno 1536
dat[serno == 1536, c("mweightc", "mstdhgtc", "msithgtc", "mmuac1c", "mmuac2c", "mmuac3c", "mhdcir1c", 
                    "mhdcir2c", "mhdcir3c", "mchest1c", "mchest2c", "mchest3c", "mabdo1c", "mabdo2c", 
                    "mabdo3c", "mtricp1c", "mtricp2c", "mtricp3c", "mbicep1c", "mbicep2c", "mbicep3c",     
                    "msubscp1c", "msubscp2c", "msubscp3c", "msupri1c", "msupri2c", "msupri3c")  := 
      list(16.9,	113.6, 88.1, 16.2, 16.4, 16.2, 49.8, 49.3, 49.7, 52.2, 52.4, 52.5, 50.9, 50.3, 49.8, 
           7.2, 6.8, 7, 4.8, 5, 4.8, 7.8, 8, 8, 5.2, 5.4, 5.4), ]
# serno 1939
dat[serno == 1939, c("mweightc", "mstdhgtc", "msithgtc", "mmuac1c", "mmuac2c", "mmuac3c", "mhdcir1c", 
                     "mhdcir2c", "mhdcir3c", "mchest1c", "mchest2c", "mchest3c", "mabdo1c", "mabdo2c", 
                     "mabdo3c", "mtricp1c", "mtricp2c", "mtricp3c", "mbicep1c", "mbicep2c", "mbicep3c",     
                     "msubscp1c", "msubscp2c", "msubscp3c", "msupri1c", "msupri2c", "msupri3c")  := 
      list(16, 111.2, 86.6, 14.8, 14.6, 14.5, 47.9, 47.8, 47.6, 53.2, 53.2, 53.1, 47.2, 47.4, 47.6, 6.2, 
           6, 6, 4, 4.2, 4.2, 5.2, 5.4, 5.8, 3.2, 3.6, 3.4), ]
# serno 2587
dat[serno == 2587, c("mweightc", "mstdhgtc", "msithgtc", "mmuac1c", "mmuac2c", "mmuac3c", "mhdcir1c", 
                     "mhdcir2c", "mhdcir3c", "mchest1c", "mchest2c", "mchest3c", "mabdo1c", "mabdo2c", 
                     "mabdo3c", "mtricp1c", "mtricp2c", "mtricp3c", "mbicep1c", "mbicep2c", "mbicep3c",     
                     "msubscp1c", "msubscp2c", "msubscp3c", "msupri1c", "msupri2c", "msupri3c")  := 
      list(26.25, NA, NA, 20.6, 20.4, 20.7, 46.7, 46.6, 46.7, 66.6, 66.2, 66.2, 63.4, 63.2, 63.4, NA, NA, 
           NA, NA, NA, NA, NA, NA, NA, NA, NA, NA), ]

# take only the variable needed to create and clean anthropometry outcome
keep <- c("serno", "sex", "mcdob", "mweightc", "mstdhgtc", "msithgtc", "mmuac1c", "mmuac2c", "mmuac3c", 
          "mhdcir1c", "allocation", "BloodTestDate", "DateTest", "mhdcir2c", "mhdcir3c", "mchest1c", 
          "mchest2c", "mchest3c", "mabdo1c", "mabdo2c", "mabdo3c", "mtricp1c", "mtricp2c", "mtricp3c", 
          "mbicep1c", "mbicep2c", "mbicep3c", "msubscp1c", "msubscp2c", "msubscp3c", "msupri1c", 
          "msupri2c", "msupri3c") 
dat <- dat[, keep, with = FALSE]
rm(keep)

sum(is.na(dat$sex))
dat$serno[is.na(dat$sex)]
# sernos 981, 1536, 1939, 2587,have no sex recorded
# use the SARAS Masterlist v9 sent to CCMB to look for sex and add to the dataset.
dat[serno == 981, sex := "Boy", ]
dat[serno == 1536, sex := "Boy", ]
dat[serno == 1939, sex := "Boy", ]
dat[serno == 1548, sex := "Girl", ]
dat[serno == 2587, sex := "Boy", ]
# sernos 981 1536 1939 2587 do not have any information od date of birth, allocation group, date of blood 
# test and date of first testing. Add the information using the SARAS Masterlist version 9
dat[serno == 981, c("mcdob", "allocation", "DateTest", "BloodTestDate") := 
      list(as.Date("2007-11-28"), "B", as.Date("2014-04-23"), as.Date("2014-09-13")), ]
dat[serno == 1536, c("mcdob", "allocation", "DateTest", "BloodTestDate") := 
      list(as.Date("2009-02-18"), "A", as.Date("2016-01-20"), as.Date("2016-01-21")), ]
dat[serno == 1939, c("mcdob", "allocation", "DateTest", "BloodTestDate") := 
      list(as.Date("2009-01-07"), "A", as.Date("2016-02-17"), as.Date("2016-02-18")), ]
dat[serno == 2587, c("mcdob", "allocation", "DateTest", "BloodTestDate") := 
      list(as.Date("2011-06-08"), "A", as.Date("2016-12-28"), as.Date("2016-12-29")), ]

# add a column for comments
dat$Comments <- ""
# serno 2587 has Down syndrome. Add to the comment column
dat[serno == 2587, "Comments" := "Child has Down Syndrome",]

# create the age variable and check children ange is between 5 and 7
dat$age <- as.numeric(with(dat, (BloodTestDate - mcdob)/365.25))
summary(dat$age)
dat[age < 5, ,]
ggplot(dat, aes(x = age)) + geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 150)) + ggtitle("SARAS Kids: Age (years)") +
  xlab("Age (years)")

# checck if sex agrees with the one sent to CCMB (v9 of the masterlist)
table(dat$sex)

# Data is ready for cleaning and processing



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


#### Weight ####
with(dat, summary(mweightc))
with(dat, sd(mweightc, na.rm = T))
with(dat, hist(mweightc))

# check for outliers
cutoffU <- median(dat$mweightc) + 5*mad(dat$mweightc)
cutoffL <- median(dat$mweightc) - 5*mad(dat$mweightc)
sum(dat$mweightc > cutoffU)
sum(dat$mweightc < cutoffL)
# 6 outliers based on the criteria established in the Analysis Plan
# plot height vs weight to check whether those are plausible
ggplot(dat, aes(x = mweightc, y = mstdhgtc)) + geom_point() +
  geom_point(data = dat[dat$mweightc > cutoffU,], aes(x = mweightc, y = mstdhgtc), 
             colour = "red") + xlab("Weight (kg)") + ylab("Height (cm)")

# create a flag variable for weight
dat$flagweight <- ifelse(dat$mweightc > cutoffU, 1, 0)
dat$flagweight <- factor(dat$flagweight, labels = c("No Outlier", "Outlier"))
table(dat$flagweight)
rm(cutoffU); rm(cutoffL)
attributes(dat$flagweight)$label <- "saras - weight outlier?"

# summary and plot (no outlier removed)
summary(dat$mweightc)
ggplot(dat, aes(x = mweightc)) + geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 150)) + ggtitle("SARAS Kids: Weight (kg)") +
  xlab("Weight (kg)")
# excluding the outliers
dat %>% filter(flagweight == "No Outlier") %>% 
  summarise(n = n(), Min = min(mweightc), Mean = mean(mweightc), SD = 
              sd(mweightc), Median = median(mweightc), Q1 = 
              quantile(mweightc, probs = 0.25), Q3 = 
              quantile(mweightc, probs = 0.75), Max = max(mweightc))
ggplot(subset(dat, flagweight == "No Outlier"), aes(x = mweightc)) + 
  geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 90)) + 
  ggtitle("SARAS Kids (No Outliers): Weight (kg)") +
  xlab("Weight (kg)")



#### Height ####
with(dat, summary(mstdhgtc))
with(dat, sd(mstdhgtc, na.rm = T))
with(dat, hist(mstdhgtc))

# check for outliers
cutoffU <- median(dat$mstdhgtc, na.rm = T) + 5*mad(dat$mstdhgtc, na.rm = T)
cutoffL <- median(dat$mstdhgtc, na.rm = T) - 5*mad(dat$mstdhgtc, na.rm = T)
# no outliers based on criteria set in the analysis plan
rm(cutoffU); rm(cutoffL)

ggplot(dat, aes(x = mstdhgtc)) + geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 75)) + ggtitle("SARAS Kids: Height (cm)") +
  xlab("Height (cm)")


#### Sitting Height ####

# sitting height is recorded including the height of the stool.
# Harshad email about the height of the stool on 22/06/2017:
# Height of Stool (Mother and Father) =31.54
# Height of Stool (Child) = 25.50

# replace the standing height variable with one with stool height subtracted

dat$msithgtc <- dat$msithgtc - 25.5

with(dat, summary(msithgtc))
with(dat, sd(msithgtc, na.rm = T))
with(dat, hist(msithgtc))

# check for outliers
cutoffU <- median(dat$msithgtc, na.rm = T) + 5*mad(dat$msithgtc, na.rm = T)
cutoffL <- median(dat$msithgtc, na.rm = T) - 5*mad(dat$msithgtc, na.rm = T)
# no outliers based on criteria set in the analysis plan
rm(cutoffU); rm(cutoffL)

ggplot(dat, aes(x = msithgtc)) + geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 85)) + ggtitle("SARAS Kids: Sitting height (cm)") +
  xlab("Sitting height (cm)")


#### MUAC ####
with(dat, summary(mmuac1c))
with(dat, sd(mmuac1c, na.rm = T))
with(dat, summary(mmuac2c))
with(dat, sd(mmuac2c, na.rm = T))
with(dat, summary(mmuac3c))
with(dat, sd(mmuac3c, na.rm = T))
# compare the distribution of MUAC 1, 2 and 3
par(mfrow = c(3,1))
with(dat, hist(mmuac1c))
with(dat, hist(mmuac2c))
with(dat, hist(mmuac3c))
# generate a variable looking at the difference between the measures
diff <- dat$mmuac2c - dat$mmuac1c
summary(diff)
rm(diff)
diff <- dat$mmuac3c - dat$mmuac1c
summary(diff)
# check those that are more than 10 cm difference
dat[diff > 10, c("serno", "mmuac1c", "mmuac2c", "mmuac3c"), with = FALSE]
rm(diff)
diff <- dat$mmuac3c - dat$mmuac2c
summary(diff)
# check those that are more than 10 cm difference
dat[diff > 10, c("serno", "mmuac1c", "mmuac2c", "mmuac3c"), with = FALSE]
rm(diff)
# serno 18 values of mmuac3c should be 16.9 (checked with Harshad on 17/06/2017)
dat[serno == 18, "mmuac3c" := 16.9,]

# compute average MUAC
dat$avgMUAC <- with(dat, (mmuac1c + mmuac2c + mmuac3c)/3)
summary(dat$avgMUAC)
attributes(dat$avgMUAC)$label <- "saras - average muac (cm)"

# check for outliers
cutoffU <- median(dat$avgMUAC, na.rm = T) + 5*mad(dat$avgMUAC, na.rm = T)
cutoffL <- median(dat$avgMUAC, na.rm = T) - 5*mad(dat$avgMUAC, na.rm = T)
sum(dat$avgMUAC > cutoffU, na.rm = T)
sum(dat$avgMUAC < cutoffL, na.rm = T)
# 2 outliers based on the criteria established in the Analysis Plan (heaviest children)
dat[dat$avgMUAC > cutoffU,,]
attributes(dat$avgMUAC)$label <- "saras - average MUAC (cm)"

# create a flag variable for MUAC
dat$flagMUAC <- ifelse(dat$avgMUAC > cutoffU, 1, 0)
dat$flagMUAC <- factor(dat$flagMUAC, labels = c("No Outlier", "Outlier"))
table(dat$flagMUAC)
rm(cutoffU); rm(cutoffL)
attributes(dat$flagMUAC)$label <- "saras - MUAC outlier?"


# summary and plot (no outlier removed)
summary(dat$avgMUAC)
ggplot(dat, aes(x = avgMUAC)) + geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 120)) + ggtitle("SARAS Kids: MUAC (cm)") +
  xlab("MUAC (cm)")
# excluding the outliers
dat %>% filter(flagMUAC == "No Outlier") %>% 
  summarise(n = n(), Min = min(avgMUAC), Mean = mean(avgMUAC), SD = 
              sd(avgMUAC), Median = median(avgMUAC), Q1 = 
              quantile(avgMUAC, probs = 0.25), Q3 = 
              quantile(avgMUAC, probs = 0.75), Max = max(avgMUAC))
ggplot(subset(dat, flagMUAC == "No Outlier"), aes(x = avgMUAC)) + 
  geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 90)) + 
  ggtitle("SARAS Kids (No Outliers): MUAC (cm)") +
  xlab("MUAC (cm)")



#### Head Circumferece ####
with(dat, summary(mhdcir1c))
with(dat, sd(mhdcir1c))
with(dat, summary(mhdcir2c))
with(dat, sd(mhdcir2c))
with(dat, summary(mhdcir3c))
with(dat, sd(mhdcir3c))
# compare the distribution of HC 1, 2 and 3
par(mfrow = c(3,1))
with(dat, hist(mhdcir1c))
with(dat, hist(mhdcir2c))
with(dat, hist(mhdcir3c))
# generate a variable looking at the difference between the measures
diff <- dat$mhdcir2c - dat$mhdcir1c
summary(diff)
rm(diff)
diff <- dat$mhdcir3c - dat$mhdcir1c
summary(diff)
rm(diff)
diff <- dat$mhdcir3c - dat$mhdcir2c
summary(diff)
rm(diff)

# average the HC measures and create an average HC variable
dat$avgHC <- (dat$mhdcir1c +  dat$mhdcir2c + dat$mhdcir3c)/3
summary(dat$avgHC)
sd(dat$avgHC)
par(mfrow = c(1,1))
hist(dat$avgHC)
attributes(dat$avgHC)$label <- "saras - average HC (cm)"

# check for outliers
cutoffU <- median(dat$avgHC) + 5*mad(dat$avgHC)
cutoffL <- median(dat$avgHC) - 5*mad(dat$avgHC)
sum(dat$avgHC > cutoffU)
sum(dat$avgHC < cutoffL)
# 1 outlier based on criteria set in the analysis plan
dat[dat$avgHC < cutoffL, , ]

# create a flag variable for HC
dat$flagHC <- ifelse(dat$avgHC < cutoffL, 1, 0)
dat$flagHC <- factor(dat$flagHC, labels = c("No Outlier", "Outlier"))
table(dat$flagHC)
rm(cutoffU); rm(cutoffL)
attributes(dat$flagHC)$label <- "saras - HC outlier?"

# summary and plot (no outlier removed)
summary(dat$avgHC)
ggplot(dat, aes(x = avgHC)) + geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) + ggtitle("SARAS Kids: Head circumference (cm)") +
  xlab("Head circumference (cm)")
# excluding the outliers
dat %>% filter(flagHC == "No Outlier") %>% 
  summarise(n = n(), Min = min(avgHC), Mean = mean(avgHC), SD = 
              sd(avgHC), Median = median(avgHC), Q1 = 
              quantile(avgHC, probs = 0.25), Q3 = 
              quantile(avgHC, probs = 0.75), Max = max(avgHC))
ggplot(subset(dat, flagHC == "No Outlier"), aes(x = avgHC)) + 
  geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 90)) + 
  ggtitle("SARAS Kids (No Outliers): HC (cm)") +
  xlab("HC (cm)")



#### Waist Circumference ####
with(dat, summary(mabdo1c))
with(dat, sd(mabdo1c, na.rm = T))
with(dat, summary(mabdo2c))
with(dat, sd(mabdo2c, na.rm = T))
with(dat, summary(mabdo3c))
with(dat, sd(mabdo3c, na.rm = T))
# compare the distribution of waist 1, 2 and 3
par(mfrow = c(3,1))
with(dat, hist(mabdo1c))
with(dat, hist(mabdo2c))
with(dat, hist(mabdo3c))
# generate a variable looking at the difference between the measures
diff <- dat$mabdo2c - dat$mabdo1c
summary(diff)
rm(diff)
diff <- dat$mabdo3c - dat$mabdo1c
summary(diff)
# check those that are more than 4 cm difference
dat[diff > 4, c("serno", "mabdo1c", "mabdo2c", "mabdo3c"), with = FALSE]
rm(diff)
diff <- dat$mabdo3c - dat$mabdo2c
summary(diff)
rm(diff)
# measures checked with Harshad (17/06/2017). All measures are correct.

# average the AC measures and create an average AC variable
dat$avgAC <- (dat$mabdo1c +  dat$mabdo2c + dat$mabdo3c)/3
summary(dat$avgAC)
sd(dat$avgAC, na.rm = T)
par(mfrow = c(1,1))
hist(dat$avgAC)
attributes(dat$avgAC)$label <- "saras - average AC (cm)"

# check for outliers
cutoffU <- median(dat$avgAC, na.rm = T) + 5*mad(dat$avgAC, na.rm = T)
cutoffL <- median(dat$avgAC, na.rm = T) - 5*mad(dat$avgAC, na.rm = T)
sum(dat$avgAC > cutoffU, na.rm = T)
sum(dat$avgAC < cutoffL, na.rm = T)
# 4 outlier based on criteria set in the analysis plan
dat[dat$avgAC > cutoffU, , ]
# 4 heaviest kids

# create a flag variable for AC
dat$flagAC <- ifelse(dat$avgAC > cutoffU, 1, 0)
dat$flagAC <- factor(dat$flagAC, labels = c("No Outlier", "Outlier"))
table(dat$flagAC)
rm(cutoffU); rm(cutoffL)
attributes(dat$flagAC)$label <- "saras - AC outlier?"

# summary and plot (no outlier removed)
summary(dat$avgAC)
ggplot(dat, aes(x = avgAC)) + geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 160)) + 
  ggtitle("SARAS Kids: Abdominal circumference (cm)") +
  xlab("Abdominal circumference (cm)")
# excluding the outliers
dat %>% filter(flagAC == "No Outlier") %>% 
  summarise(n = n(), Min = min(avgAC), Mean = mean(avgAC), SD = 
              sd(avgAC), Median = median(avgAC), Q1 = 
              quantile(avgAC, probs = 0.25), Q3 = 
              quantile(avgAC, probs = 0.75), Max = max(avgAC))
ggplot(subset(dat, flagAC == "No Outlier"), aes(x = avgAC)) + 
  geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) + 
  ggtitle("SARAS Kids (No Outliers): AC (cm)") +
  xlab("AC (cm)")



#### Chest circumference ####
with(dat, summary(mchest1c))
with(dat, sd(mchest1c, na.rm = T))
with(dat, summary(mchest2c))
with(dat, sd(mchest2c, na.rm = T))
with(dat, summary(mchest3c))
with(dat, sd(mchest3c, na.rm = T))
# compare the distribution of waist 1, 2 and 3
par(mfrow = c(3,1))
with(dat, hist(mchest1c))
with(dat, hist(mchest2c))
with(dat, hist(mchest3c))
# generate a variable looking at the difference between the measures
diff <- dat$mchest2c - dat$mchest1c
summary(diff)
rm(diff)
diff <- dat$mchest3c - dat$mchest1c
summary(diff)
# check those that are more than 4 cm difference
dat[diff < -4, c("serno", "mchest1c", "mchest2c", "mchest3c"), with = FALSE]
rm(diff)
diff <- dat$mabdo3c - dat$mabdo2c
summary(diff)
rm(diff)
# measures checked with Harshad (17/06/2017). All measures are correct.

# average the CC measures and create an average CC variable
dat$avgCC <- (dat$mchest1c +  dat$mchest2c + dat$mchest3c)/3
summary(dat$avgCC)
sd(dat$avgCC)
par(mfrow = c(1,1))
hist(dat$avgCC)
attributes(dat$avgCC)$label <- "saras - average CC (cm)"

# check for outliers
cutoffU <- median(dat$avgCC) + 5*mad(dat$avgCC)
cutoffL <- median(dat$avgCC) - 5*mad(dat$avgCC)
sum(dat$avgCC > cutoffU)
sum(dat$avgCC < cutoffL)
# 4 outlier based on criteria set in the analysis plan
dat[dat$avgCC > cutoffU, , ]
# 4 heavies kids

# create a flag variable for CC
dat$flagCC <- ifelse(dat$avgCC > cutoffU, 1, 0)
dat$flagCC <- factor(dat$flagCC, labels = c("No Outlier", "Outlier"))
table(dat$flagCC)
rm(cutoffU); rm(cutoffL)
attributes(dat$flagCC)$label <- "saras - CC outlier?"

# summary and plot (no outlier removed)
summary(dat$avgCC)
ggplot(dat, aes(x = avgCC)) + geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 150)) + ggtitle("SARAS Kids: Chest circumference (cm)") +
  xlab("Chest circumference (cm)")
# excluding the outliers
dat %>% filter(flagCC == "No Outlier") %>% 
  summarise(n = n(), Min = min(avgCC), Mean = mean(avgCC), SD = 
              sd(avgCC), Median = median(avgCC), Q1 = 
              quantile(avgCC, probs = 0.25), Q3 = 
              quantile(avgCC, probs = 0.75), Max = max(avgCC))
ggplot(subset(dat, flagCC == "No Outlier"), aes(x = avgCC)) + 
  geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) + 
  ggtitle("SARAS Kids (No Outliers): CC (cm)") +
  xlab("CC (cm)")

#### Hip circumference ####
# NOT MEASURED IN MUMBAI


#### Triceps ####
with(dat, summary(mtricp1c))
with(dat, sd(mtricp1c, na.rm = T))
with(dat, summary(mtricp2c))
with(dat, sd(mtricp2c, na.rm = T))
with(dat, summary(mtricp3c))
with(dat, sd(mtricp3c, na.rm = T))
# compare the distribution of triceps 1, 2 and 3
par(mfrow = c(3,1))
with(dat, hist(mtricp1c))
with(dat, hist(mtricp2c))
with(dat, hist(mtricp3c))
# generate a variable looking at the difference between the measures
diff <- dat$mtricp2c - dat$mtricp1c
summary(diff)
rm(diff)
diff <- dat$mtricp3c - dat$mtricp1c
summary(diff)
rm(diff)
diff <- dat$mtricp3c - dat$mtricp2c
summary(diff)
rm(diff)

# average the triceps measures and create an average triceps variable
dat$avgTriceps <- (dat$mtricp1c +  dat$mtricp3c + dat$mtricp3c)/3
summary(dat$avgTriceps)
sd(dat$avgTriceps, na.rm = T)
par(mfrow = c(1,1))
hist(dat$avgTriceps)
attributes(dat$avgTriceps)$label <- "saras - average triceps (mm)"

# check for outliers
cutoffU <- median(dat$avgTriceps, na.rm = T) + 5*mad(dat$avgTriceps, na.rm = T)
cutoffL <- median(dat$avgTriceps, na.rm = T) - 5*mad(dat$avgTriceps, na.rm = T)
sum(dat$avgTriceps > cutoffU, na.rm = T)
sum(dat$avgTriceps < cutoffL, na.rm = T)
# 8 outliers have been identified using the criteria on the analysis plan
# check the outliers before removing
dat[dat$avgTriceps > cutoffU, , ]

# create a flag variable for triceps
dat$flagTriceps <- ifelse(dat$avgTriceps > cutoffU, 1, 0)
dat$flagTriceps <- factor(dat$flagTriceps, labels = c("No Outlier", "Outlier"))
table(dat$flagTriceps)
rm(cutoffU); rm(cutoffL)
attributes(dat$flagTriceps)$label <- "saras - Triceps outlier?"

# summary and plot (no outlier removed)
summary(dat$avgTriceps)
ggplot(dat, aes(x = avgTriceps)) + geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 120)) + ggtitle("SARAS Kids: Triceps (mm)") +
  xlab("Triceps (mm)")
# excluding the outliers
dat %>% filter(flagTriceps == "No Outlier") %>% 
  summarise(n = n(), Min = min(avgTriceps), Mean = mean(avgTriceps), SD = 
              sd(avgTriceps), Median = median(avgTriceps), Q1 = 
              quantile(avgTriceps, probs = 0.25), Q3 = 
              quantile(avgTriceps, probs = 0.75), Max = max(avgTriceps))
ggplot(subset(dat, flagTriceps == "No Outlier"), aes(x = avgTriceps)) + 
  geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) + 
  ggtitle("SARAS Kids (No Outliers): Triceps (mm)") +
  xlab("Triceps (mm)")




#### Biceps ####
with(dat, summary(mbicep1c))
with(dat, sd(mbicep1c, na.rm = T))
with(dat, summary(mbicep2c))
with(dat, sd(mbicep2c, na.rm = T))
with(dat, summary(mbicep3c))
with(dat, sd(mbicep3c, na.rm = T))
# compare the distribution of biceps 1, 2 and 3
par(mfrow = c(3,1))
with(dat, hist(mbicep1c))
with(dat, hist(mbicep2c))
with(dat, hist(mbicep3c))
# generate a variable looking at the difference between the measures
diff <- dat$mbicep2c - dat$mbicep1c
summary(diff)
rm(diff)
diff <- dat$mbicep3c - dat$mbicep1c
summary(diff)
rm(diff)
diff <- dat$mbicep3c - dat$mbicep2c
summary(diff)
rm(diff)

# average the biceps measures and create an average biceps variable
dat$avgBiceps <- (dat$mbicep1c +  dat$mbicep2c + dat$mbicep3c)/3
summary(dat$avgBiceps)
sd(dat$avgBiceps, na.rm = T)
par(mfrow = c(1,1))
hist(dat$avgBiceps)
attributes(dat$avgBiceps)$label <- "saras - average biceps (mm)"

# check for outliers
cutoffU <- median(dat$avgBiceps, na.rm = T) + 5*mad(dat$avgBiceps, na.rm = T)
cutoffL <- median(dat$avgBiceps, na.rm = T) - 5*mad(dat$avgBiceps, na.rm = T)
sum(dat$avgBiceps > cutoffU, na.rm = T)
sum(dat$avgBiceps < cutoffL, na.rm = T)
# 5 outliers have been identified using the criteria on the analysis plan
# check the outliers before removing
dat[dat$avgBiceps > cutoffU, , ]

# create a flag variable for biceps
dat$flagBiceps <- ifelse(dat$avgBiceps > cutoffU, 1, 0)
dat$flagBiceps <- factor(dat$flagBiceps, labels = c("No Outlier", "Outlier"))
table(dat$flagBiceps)
rm(cutoffU); rm(cutoffL)
attributes(dat$flagBiceps)$label <- "saras - Biceps outlier?"

# summary and plot (no outlier removed)
summary(dat$avgBiceps)
ggplot(dat, aes(x = avgBiceps)) + geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) + ggtitle("SARAS Kids: Biceps (mm)") +
  xlab("Biceps (mm)")
# excluding the outliers
dat %>% filter(flagBiceps == "No Outlier") %>% 
  summarise(n = n(), Min = min(avgBiceps), Mean = mean(avgBiceps), SD = 
              sd(avgBiceps), Median = median(avgBiceps), Q1 = 
              quantile(avgBiceps, probs = 0.25), Q3 = 
              quantile(avgBiceps, probs = 0.75), Max = max(avgBiceps))
ggplot(subset(dat, flagBiceps == "No Outlier"), aes(x = avgTriceps)) + 
  geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 120)) + 
  ggtitle("SARAS Kids (No Outliers): Biceps (mm)") +
  xlab("Biceps (mm)")



#### Subscapular ####
with(dat, summary(msubscp1c))
with(dat, sd(msubscp1c, na.rm = T))
with(dat, summary(msubscp2c))
with(dat, sd(msubscp2c, na.rm = T))
with(dat, summary(msubscp3c))
with(dat, sd(msubscp3c, na.rm = T))
# compare the distribution of subscapular 1, 2 and 3
par(mfrow = c(3,1))
with(dat, hist(msubscp1c))
with(dat, hist(msubscp2c))
with(dat, hist(msubscp3c))
# generate a variable looking at the difference between the measures
diff <- dat$msubscp2c - dat$msubscp1c
summary(diff)
rm(diff)
diff <- dat$msubscp3c - dat$msubscp2c
summary(diff)
rm(diff)
diff <- dat$msubscp3c - dat$msubscp1c
summary(diff)
rm(diff)

# average the subscapular measures and create an average subscapular variable
dat$avgSubsc <- (dat$msubscp1c +  dat$msubscp2c + dat$msubscp3c)/3
summary(dat$avgSubsc)
sd(dat$avgSubsc, na.rm = T)
par(mfrow = c(1,1))
hist(dat$avgSubsc)
attributes(dat$avgSubsc)$label <- "saras - average subscapular (mm)"

# check for outliers
cutoffU <- median(dat$avgSubsc, na.rm = T) + 5*mad(dat$avgSubsc, na.rm = T)
cutoffL <- median(dat$avgSubsc, na.rm = T) - 5*mad(dat$avgSubsc, na.rm = T)
sum(dat$avgSubsc > cutoffU, na.rm = T)
sum(dat$avgSubsc < cutoffL, na.rm = T)
# 10 outliers have been identified using the criteria on the analysis plan
# check the outliers before removing
dat[dat$avgSubsc > cutoffU, , ]

# create a flag variable for subscapular
dat$flagSubsc <- ifelse(dat$avgSubsc > cutoffU, 1, 0)
dat$flagSubsc <- factor(dat$flagSubsc, labels = c("No Outlier", "Outlier"))
table(dat$flagSubsc)
rm(cutoffU); rm(cutoffL)
attributes(dat$flagSubsc)$label <- "saras - Subscapular outlier?"

# summary and plot (no outlier removed)
summary(dat$avgSubsc)
ggplot(dat, aes(x = avgSubsc)) + geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 200)) + ggtitle("SARAS Kids: Subscapular (mm)") +
  xlab("Subscapular (mm)")
# excluding the outliers
dat %>% filter(flagSubsc == "No Outlier") %>% 
  summarise(n = n(), Min = min(avgSubsc), Mean = mean(avgSubsc), SD = 
              sd(avgSubsc), Median = median(avgSubsc), Q1 = 
              quantile(avgSubsc, probs = 0.25), Q3 = 
              quantile(avgSubsc, probs = 0.75), Max = max(avgSubsc))
ggplot(subset(dat, flagSubsc == "No Outlier"), aes(x = avgSubsc)) + 
  geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 90)) + 
  ggtitle("SARAS Kids (No Outliers): Subscapular (mm)") +
  xlab("Subscapular (mm)")



#### Suprailiac ####
with(dat, summary(msupri1c))
with(dat, sd(msupri1c, na.rm = T))
with(dat, summary(msupri2c))
with(dat, sd(msupri2c, na.rm = T))
with(dat, summary(msupri3c))
with(dat, sd(msupri3c, na.rm = T))
# compare the distribution of suprailiac 1, 2 and 3
par(mfrow = c(3,1))
with(dat, hist(msupri1c))
with(dat, hist(msupri2c))
with(dat, hist(msupri3c))
# generate a variable looking at the difference between the measures
diff <- dat$msupri2c - dat$msupri1c
summary(diff)
rm(diff)
diff <- dat$msupri3c - dat$msupri2c
summary(diff)
rm(diff)
diff <- dat$msupri3c - dat$msupri1c
summary(diff)
rm(diff)

# average the suprailiac measures and create an average suprailiac variable
dat$avgSup <- (dat$msupri1c +  dat$msupri2c + dat$msupri3c)/3
summary(dat$avgSup)
sd(dat$avgSup, na.rm = T)
par(mfrow = c(1,1))
hist(dat$avgSup)
attributes(dat$avgSup)$label <- "saras - average suprailiac (mm)"

# check for outliers
cutoffU <- median(dat$avgSup, na.rm = T) + 5*mad(dat$avgSup, na.rm = T)
cutoffL <- median(dat$avgSup, na.rm = T) - 5*mad(dat$avgSup, na.rm = T)
sum(dat$avgSup > cutoffU, na.rm = T)
sum(dat$avgSup < cutoffL, na.rm = T)
# 6 outliers have been identified using the criteria on the analysis plan
# check the outliers before removing
dat[dat$avgSup > cutoffU, , ]


# create a flag variable for Suprailiac
dat$flagSup <- ifelse(dat$avgSup > cutoffU, 1, 0)
dat$flagSup <- factor(dat$flagSup, labels = c("No Outlier", "Outlier"))
table(dat$flagSup)
rm(cutoffU); rm(cutoffL)
attributes(dat$flagSup)$label <- "saras - Suprailiac outlier?"

# summary and plot (no outlier removed)
summary(dat$avgSup)
ggplot(dat, aes(x = avgSup)) + geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 200)) + ggtitle("SARAS Kids: Suprailiac (mm)") +
  xlab("Suprailiac (mm)")
# excluding the outliers
dat %>% filter(flagSup == "No Outlier") %>% 
  summarise(n = n(), Min = min(avgSup), Mean = mean(avgSup), SD = 
              sd(avgSup), Median = median(avgSubsc), Q1 = 
              quantile(avgSup, probs = 0.25), Q3 = 
              quantile(avgSup, probs = 0.75), Max = max(avgSup))
ggplot(subset(dat, flagSup == "No Outlier"), aes(x = avgSup)) + 
  geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 90)) + 
  ggtitle("SARAS Kids (No Outliers): Suprailiac (mm)") +
  xlab("Suprailiac (mm)")





#########################################################################################################
#### Multivariate Analysis ##############################################################################
#########################################################################################################

# 1. Check each variable in relation with others
# 2. Check if same measures are consistent



#### Height, Weight and Sitting Height ####

# height and sitting height should have a positive relation. Check
with(dat, plot(mstdhgtc, msithgtc))
# one child with questionable values
dat[(mstdhgtc < 110 & msithgtc < 50),,]
# set sitting height to missing
dat[serno == 1762, "msithgtc" := NA, with = FALSE]

# height and weight should have a positive relation. Check
with(dat, plot(mstdhgtc, mweightc))



#### Skinfolds ####
pairs(~ avgTriceps + avgBiceps + avgSubsc + avgSup + mweightc,
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
dat$bmi <- with(dat, mweightc/(mstdhgtc*0.01)^2)
summary(dat$bmi)
hist(dat$bmi)
attributes(dat$bmi)$label <- "saras - bmi (kg/m2)"

# check for outliers
cutoffU <- median(dat$bmi, na.rm = T) + 5*mad(dat$bmi, na.rm = T)
cutoffL <- median(dat$bmi, na.rm = T) - 5*mad(dat$bmi, na.rm = T)
sum(dat$bmi > cutoffU, na.rm = T)
sum(dat$bmi < cutoffL, na.rm = T)
# 9 outliers have been identified using the criteria on the analysis plan
# check the outliers
dat[dat$bmi > cutoffU | dat$bmi < cutoffL, , ]

# create a flag variable for bmi
dat$flagbmi <- ifelse(dat$bmi > cutoffU | dat$bmi < cutoffL, 1, 0)
dat$flagbmi <- factor(dat$flagbmi, labels = c("No Outlier", "Outlier"))
table(dat$flagbmi)
rm(cutoffU); rm(cutoffL)
attributes(dat$flagbmi)$label <- "saras - bmi outlier?"

# summary and plot (no outlier removed)
summary(dat$bmi)
ggplot(dat, aes(x = bmi)) + geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 200)) + 
  ggtitle(expression(paste("SARAS Kids: BMI (", kg/m^{2},")"))) +
  xlab(expression(paste("BMI (", kg/m^{2}, ")")))
# excluding the outliers
dat %>% filter(flagbmi == "No Outlier") %>% 
  summarise(n = n(), Min = min(bmi), Mean = mean(bmi), SD = 
              sd(bmi), Median = median(bmi), Q1 = 
              quantile(bmi, probs = 0.25), Q3 = 
              quantile(bmi, probs = 0.75), Max = max(bmi))
ggplot(subset(dat, flagbmi == "No Outlier"), aes(x = bmi)) + 
  geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) + 
  ggtitle(expression(paste("SARAS Kids (No Outliers): BMI (", kg/m^{2},")"))) +
  xlab(expression(paste("BMI (", kg/m^{2}, ")")))



#### Compute Leg Length ####
with(dat, plot(mstdhgtc, msithgtc))
dat$leg <- with(dat, (mstdhgtc - msithgtc))
summary(dat$leg)
hist(dat$leg)
attributes(dat$leg)$label <- "saras - leg length (cm)"

# check for outliers
cutoffU <- median(dat$leg, na.rm = T) + 5*mad(dat$leg, na.rm = T)
cutoffL <- median(dat$leg, na.rm = T) - 5*mad(dat$leg, na.rm = T)
sum(dat$leg > cutoffU, na.rm = T)
sum(dat$leg < cutoffL, na.rm = T)

# no outlier were identified
ggplot(dat, aes(x = leg)) + geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 200)) + 
  ggtitle("SARAS Kids: Leg length (cm)") +
  xlab("Leg length(cm)")


#### Compute Leg Length to Sitting Height Ratio ####
with(dat, plot(leg, msithgtc))
dat$legshtratio <- with(dat, (leg/msithgtc))
summary(dat$legshtratio)
hist(dat$legshtratio)
attributes(dat$legshtratio)$label <- "saras - leg length to sitting height ratio"

# check for outliers
cutoffU <- median(dat$legshtratio, na.rm = T) + 5*mad(dat$legshtratio, na.rm = T)
cutoffL <- median(dat$legshtratio, na.rm = T) - 5*mad(dat$legshtratio, na.rm = T)
sum(dat$legshtratio > cutoffU, na.rm = T)
sum(dat$legshtratio < cutoffL, na.rm = T)


# no outlier were identified
ggplot(dat, aes(x = legshtratio)) + geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 200)) + 
  ggtitle("SARAS Kids: Leg length to Sitting Height") +
  xlab("Leg Length to Sitting Height")



#### Total sum of skinfold ####
dat$totsskin <- with(dat, (avgTriceps + avgBiceps + avgSubsc + avgSup))
attributes(dat$totsskin)$label <- "saras - total sum of skinfolds"
summary(dat$totsskin)
hist(dat$totsskin)

# check for outliers
cutoffU <- median(dat$totsskin, na.rm = T) + 5*mad(dat$totsskin, na.rm = T)
cutoffL <- median(dat$totsskin, na.rm = T) - 5*mad(dat$totsskin, na.rm = T)
sum(dat$totsskin > cutoffU, na.rm = T)
sum(dat$totsskin < cutoffL, na.rm = T)
# 7 outliers have been identified using the criteria on the analysis plan
# check the outliers
dat[dat$totsskin > cutoffU | dat$bmi < cutoffL, , ]
# most of them have at least one or two skinfold measures that were deemed as outliers in
# the univariate analysis

# create a flag variable for total sum of skinfold
dat$flagtotsskin <- ifelse(dat$totsskin > cutoffU | dat$totsskin < cutoffL, 1, 0)
dat$flagtotsskin <- factor(dat$flagtotsskin, labels = c("No Outlier", "Outlier"))
table(dat$flagtotsskin)
rm(cutoffU); rm(cutoffL)
attributes(dat$flagtotsskin)$label <- "saras - total skinfolds outlier?"

# summary and plot (no outlier removed)
summary(dat$totsskin)
ggplot(dat, aes(x = totsskin)) + geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 150)) +
  ggtitle("SARAS Kids: Total sum of skinfolds (mm)") + xlab("Total sum of skinfolds (mm)")
# excluding the outliers
dat %>% filter(flagbmi == "No Outlier") %>% 
  summarise(n = n(), Min = min(totsskin), Mean = mean(totsskin), SD = 
              sd(totsskin), Median = median(totsskin), Q1 = 
              quantile(totsskin, probs = 0.25), Q3 = 
              quantile(totsskin, probs = 0.75), Max = max(totsskin))
ggplot(subset(dat, flagtotsskin == "No Outlier"), aes(x = totsskin)) + 
  geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) +
  ggtitle("SARAS Kids (No Outliers): Total sum of skinfolds (mm)") + 
  xlab("Total sum of skinfolds (mm)")




#### Compute SD scores based on WHO criteria ####

# SD score are computed based on the WHO criteria
# The WHO provides on their website an already-made macro for computing z-scores.
# The macro was downloaded (18/05/2017) and saved as who2007_R. Instruction for using the macro can be
# found on Readme.pdf

# Restore reference data sets:
wfawho2007 <- 
  read.table("C:/Users/chiara.digravio/Desktop/EMPHASIS/PhenotypeDataProcessing/who2007_R/wfawho2007.txt",
             header=T,sep="",skip=0) 
hfawho2007 <- 
  read.table("C:/Users/chiara.digravio/Desktop/EMPHASIS/PhenotypeDataProcessing/who2007_R/hfawho2007.txt",
             header=T,sep="",skip=0) 
bfawho2007 <- 
  read.table("C:/Users/chiara.digravio/Desktop/EMPHASIS/PhenotypeDataProcessing/who2007_R/bfawho2007.txt",
             header=T,sep="",skip=0) 

# source the who2007.R file
source("C:/Users/chiara.digravio/Desktop/EMPHASIS/PhenotypeDataProcessing/who2007_R/who2007.r")

# take only the variable needed
datZscore <- dat %>% select(serno, age, sex, mweightc, mstdhgtc, bmi)

# If sex is character, it must be “m” or “M” for males and “f” or “F” for females.
# Transform sex in character
table(datZscore$sex)
levels(datZscore$sex) <- c("M", "F")
datZscore$sex <- as.character(datZscore$sex)
# age must be in months not rounded. Transform age in month
datZscore$age <- datZscore$age*12
# height and weight must be in centimeter and kg respectively

# transfrom Zscore in dataframe
datZscore <- data.frame(datZscore)
str(datZscore)

# strip variable of their attributes
attributes(datZscore$age) <- NULL
attributes(datZscore$mweightc) <- NULL
attributes(datZscore$mstdhgtc) <- NULL
attributes(datZscore$bmi) <- NULL
str(datZscore)


who2007(FileLab = "SARASKidsZscore",
        FilePath = "C:/Users/chiara.digravio/Desktop/EMPHASIS/PhenotypeDataProcessing",
        mydf = datZscore, sex = sex, age = age, weight = mweightc, height = mstdhgtc)

# remove unneed functions and dataset
rm(bfawho2007); rm(datZscore); rm(hfawho2007); rm(matprev); rm(matz); rm(wfawho2007) 
rm(calc.zbmi); rm(calc.zhfa); rm(calc.zwei); rm(prevnh); rm(prevnh.L); rm(prevph)
rm(prevph.L); rm(rounde); rm(who2007); rm(wmean); rm(wsd)

# Merge dataset with WHO z-score 
datZscore <- read.csv("SARASKidsZscore_z.csv")
keep <- c("serno", "zhfa", "zwfa", "zbfa")
datZscore <- datZscore[, names(datZscore) %in% keep]
dat <- inner_join(dat, datZscore)
dat <- data.frame(dat)

# remove unused data
rm(datZscore); rm(keep)

attributes(dat$zwfa)$label <- "saras - weight for age (WHO z-score)"
attributes(dat$zhfa)$label <- "saras - height for age (WHO z-score)"


#### Height for age (WHO z-score) ####

summary(dat$zhfa)
sd(dat$zhfa, na.rm = T)
# there are 2 NAs
# serno 2587 (height not measured) and 1459 (age < 5)
# we can add the z-score for 1459 by looking at the WHO table for the corresponding age
dat <- data.table(dat)
dat[serno == 1459, c("zhfa", "zwfa", "zbfa") := list(-0.93, -2.16, -2.13), ]
ggplot(dat, aes(x = zhfa)) + geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 80)) +
  ggtitle("SARAS Kids: Height-for-age (WHO z-score)") + xlab("Height-for-age (WHO z-score)")


#### Weight for age (WHO z-score) ####

summary(dat$zwfa)
sd(dat$zwfa, na.rm = T)
ggplot(dat, aes(x = zwfa)) + geom_histogram(colour="black", fill="white") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 110)) +
  ggtitle("SARAS Kids: Weight-for-age (WHO z-score)") + xlab("Weight-for-age (WHO z-score)")

# check possible outliers 
cutoffU <- median(dat$zwfa, na.rm = T) + 5*mad(dat$zwfa, na.rm = T)
cutoffL <- median(dat$zwfa, na.rm = T) - 5*mad(dat$zwfa, na.rm = T)
sum(dat$zwfa > cutoffU, na.rm = T)
sum(dat$zwfa < cutoffL, na.rm = T)
# 1 outliers based on the criteria established in the analysis plan
dat[dat$zwfa > cutoffU, ,]
# heaviest child in the cohort

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

# select only relevant variable (based on the analysis plan)
dat <- dat %>% select(serno, sex, age, allocation, mweightc, mstdhgtc, msithgtc, avgMUAC, 
                      avgHC, avgAC,  avgCC,  avgTriceps, avgBiceps,  avgSubsc,  avgSup, bmi,
                      leg, legshtratio, totsskin, zhfa, zwfa, stunting, wasting, underweight,
                      flagweight, flagMUAC, flagHC, flagAC, flagCC, flagTriceps, flagBiceps, flagSubsc,
                      flagSup,  flagbmi, flagtotsskin,  Comments)
                      
                      
  
names(dat) <- c("serno", "csex", "cage", "allocation", "cwt", "cht", "csitht", "cavgMUAC", "cavgHC", 
                "cavgWaist", "cavgCC", "cavgTric", "cavgBic", "cavgSubsc", "cavgSup", "cbmi",
                "cleg", "clegsithtratio", "ctotsskin", "chfa", "cwfa", "cstunting", "cwasting",
                "cunderweight", "flagcwt", "flagcavgMUAC", "flagcavgHC", "flagcavgWaist", 
                "flagcavgCC", "flagcavgTric", "flagcavgBic", "flagcavgSubsc",
                "flagcavgSup",  "flagcbmi", "flagctotsskin", "Comments")
# add attributes
attributes(dat$serno)$label <- "saras - Subject ID"
attributes(dat$csex)$label <- "saras - child sex"
attributes(dat$cage)$label <- "saras - child age (years)"
attributes(dat$cstunting)$label <- "saras - stunting"
attributes(dat$cwasting)$label <- "saras - wasting"
attributes(dat$cunderweight)$label <- "saras - underweight"
attributes(dat$allocation)$label <- "saras - allocation group"
attributes(dat$cwt)$label <- "saras - child weight (kg)"
attributes(dat$cht)$label <- "saras - child height (cm)"
attributes(dat$csitht)$label <- "saras - child sitting height (cm)"


# add body composition information from DXA data
source("SARASKidsPhenotypeCleaningDXA.R")
datDXA <- datDXA %>% select(serno, totfatmass_kg, totleanmass_kg, androidfatmass_kg, gynoidfatmass_kg,
                            fmi, lmi, fatper, flagfatmass, flagleanmass, flagandroidfatmass, 
                            flaggynoidfatmass, flagfmi, flaglmi, flagfatper)
names(datDXA)[names(datDXA) == "totfatmass_kg"] <- "cfatmass"
names(datDXA)[names(datDXA) == "totleanmass_kg"] <- "cleanmass"
names(datDXA)[names(datDXA) == "androidfatmass_kg"] <- "candroidfat"
names(datDXA)[names(datDXA) == "gynoidfatmass_kg"] <- "cgynoiddfat"
names(datDXA)[names(datDXA) == "fmi"] <- "cfmi"
names(datDXA)[names(datDXA) == "lmi"] <- "clmi"
names(datDXA)[names(datDXA) == "fatper"] <- "cfatper"
names(datDXA)[names(datDXA) == "flagfatmass"] <- "flagcfatmass"
names(datDXA)[names(datDXA) == "flagleanmass"] <- "flagcleanmass"
names(datDXA)[names(datDXA) == "flagandroidfatmass"] <- "flagcandroidfat"
names(datDXA)[names(datDXA) == "flaggynoidfatmass"] <- "flagcgynoidfat"
names(datDXA)[names(datDXA) == "flagfmi"] <- "flagcfmi"
names(datDXA)[names(datDXA) == "flaglmi"] <- "flagclmi"
names(datDXA)[names(datDXA) == "flagfatper"] <- "flagcfatper"

dat <- left_join(dat, datDXA, all = T)
  
  
# save as csv, stata and spss file
write.table(dat, "./CleanedDataset/SARASKidsGrowthBodyCompv2_20170713.csv", row.names = FALSE)
write_dta(dat, "./CleanedDataset/SARASKidsGrowthBodyCompv2_20170713.dta")
write_sav(dat, "./CleanedDataset/SARASKidsGrowthBodyCompv2_20170713.sav")
