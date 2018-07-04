############################################################################################################################
##### PMMST Phenotype Data Processing Birth Outcome ########################################################################
############################################################################################################################

# 03/02/2018

# load packages
library(readxl)
library(dplyr)
library(ggplot2)
library(haven)


# load datasets with birth outcomes
dat <- read_excel("./PMMST Outcome Data/PMMST_pdata_redux_170718.xlsx", na = "NA", col_names = TRUE)
# sanity check on imported dataset
dim(dat)
str(dat)
# Transform every ID in upper case
dat$cSubjectID <- toupper(dat$cSubjectID)
# live births?
table(dat$Live_bth) 
# some have comments. Check those with comments
dat %>% filter(!is.na(comments)) %>% dplyr::select(cSubjectID, Bth_wt, Bth_len, hc, LMP, comments)
# some children where deliverd at the coast; hence, they do not have any information
# those have values for birth measure that is negative. Set to missing
dat$Bth_wt[dat$Bth_wt < 0] <- NA
dat$Bth_len[dat$Bth_len < 0] <- NA
dat$hc[dat$hc < 0] <- NA
# in the comments we could retreive weight from EMP099V:
dat$Bth_wt[dat$cSubjectID == "EMP099V"] <- 3.46
# check these coastal delivery with information on children computed based on confounders
datcoast <- read_dta("./Cleaned Dataset/PMMSTConfoundersv1_20180212.dta") %>%
  dplyr::select(serno, ccoast) %>% filter(ccoast == 2)
# only one remained at the coast at visit 2
rm(datcoast)

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
summary(dat$Bth_wt)
# 23 (7%) do not have any birth weight recorded
sd(dat$Bth_wt, na.rm = T)
hist(dat$Bth_wt)

# check for outliers
cutoffU <- median(dat$Bth_wt, na.rm = T) + 5*mad(dat$Bth_wt, na.rm = T)
cutoffL <- median(dat$Bth_wt, na.rm = T) - 5*mad(dat$Bth_wt, na.rm = T)
sum(dat$Bth_wt > cutoffU, na.rm = T)
sum(dat$Bth_wt < cutoffL, na.rm = T)
# no outliers based on the criteria established in the Analysis Plan
rm(cutoffU); rm(cutoffL)

# transform birth weight in g
dat$Bth_wt <- dat$Bth_wt*1000

ggplot(dat, aes(x = Bth_wt)) + geom_histogram(colour="black", fill="white", binwidth = 150) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 55)) + 
  ggtitle("PMMST: Birthweight (g)") +
  xlab("Birthweight (g)")


#### Birth length ####
summary(dat$Bth_len)
# 38 (12%) do not have any birth length recorded
sd(dat$Bth_len, na.rm = T)
hist(dat$Bth_len)

# check for outliers
cutoffU <- median(dat$Bth_len, na.rm = T) + 5*mad(dat$Bth_len, na.rm = T)
cutoffL <- median(dat$Bth_len, na.rm = T) - 5*mad(dat$Bth_len, na.rm = T)
sum(dat$Bth_len > cutoffU, na.rm = T)
sum(dat$Bth_len < cutoffL, na.rm = T)
# 1 outlier based on the criteria established in the Analysis Plan
# check outlier
dat[dat$Bth_len < cutoffL & !is.na(dat$Bth_len), c("Bth_len", "Bth_wt")]
with(dat, plot(Bth_wt, Bth_len))
# birth length too low. Set to missing
dat[dat$Bth_len < cutoffL & !is.na(dat$Bth_len), c("Bth_len")] <- NA
summary(dat$Bth_len)

# check for outliers
cutoffU <- median(dat$Bth_len, na.rm = T) + 5*mad(dat$Bth_len, na.rm = T)
cutoffL <- median(dat$Bth_len, na.rm = T) - 5*mad(dat$Bth_len, na.rm = T)
sum(dat$Bth_len > cutoffU, na.rm = T)
sum(dat$Bth_len < cutoffL, na.rm = T)
# 0 outliers based on the criteria established in the Analysis Plan

# plot 
ggplot(dat, aes(x = Bth_len)) + geom_histogram(colour="black", fill="white", binwidth = 0.9) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 55)) + 
  ggtitle("PMMST: Birth length (cm)") +
  xlab("PMMST (cm)")


#### Head Circumference ####

summary(dat$hc)
# 37 (12%) do not have any birth head circumference recorded
sd(dat$hc, na.rm = T)
hist(dat$hc)

# check for outliers
cutoffU <- median(dat$hc, na.rm = T) + 5*mad(dat$hc, na.rm = T)
cutoffL <- median(dat$hc, na.rm = T) - 5*mad(dat$hc, na.rm = T)
sum(dat$hc > cutoffU, na.rm = T)
sum(dat$hc < cutoffL, na.rm = T)
# 1 outlier based on the criteria established in the Analysis Plan

# create a variable to flag outliers for chest circumference
dat$flagbirthHC <- ifelse(dat$hc > cutoffU | dat$hc < cutoffL, 1, 0)
dat$flagbirthHC <- factor(dat$flagbirthHC, labels = c("No Outlier", "Outlier"))
attributes(dat$flagbirthHC)$label <- "pmmst - birth HC outlier?"
rm(cutoffU); rm(cutoffL)

# plot (with outlier)
ggplot(dat, aes(x = hc)) + geom_histogram(colour="black", fill="white", binwidth = 1) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 90)) + 
  ggtitle("PMMST: Head circumference (cm)") +
  xlab("Head circumference (cm)")

# remove the outlier
dat %>% filter(flagbirthHC == "No Outlier") %>% 
  summarise(n = n(), Min = min(hc), Mean = mean(hc), SD = 
              sd(hc), Median = median(hc), Q1 = 
              quantile(hc, probs = 0.25), Q3 = 
              quantile(hc, probs = 0.75))
ggplot(subset(dat, flagbirthHC == "No Outlier"), aes(x = hc)) + 
  geom_histogram(colour="black", fill="white", binwidth = 0.4) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 60)) + 
  ggtitle("PMMST (No Outliers): Head circumference (cm)") +
  xlab("Head circumference (cm)")

# chest circumference, abdominal circumference, skinfolds and MUAC not available.


#### Gestational age ####

summary(dat$BthGAwkscombi)
# 6 have missing gestational age. 
dat[dat$BthGAwkscombi < 0 & !is.na(dat$BthGAwkscombi), ]
# One has negative value. Need to check. LMP and conception date are incorrect. Set to NA
# together with gestational age
dat[dat$cSubjectID == "EMP202X", c("LMP", "ConceptionDate", "BthGAwkscombi")] <- NA
summary(dat$BthGAwkscombi)
hist(dat$BthGAwkscombi)
# 18.29 seems to little. Set to missing
dat[which.min(dat$BthGAwkscombi), c("LMP", "BthGAwkscombi")] <- NA
# 8 missing gestational age


#######################################################################################################################
#### Multivariate Analysis ############################################################################################
#######################################################################################################################

# 1. Check each variable in relation with others
# 2. Check if same measures are consistent


#### Birth weight and length ####
with(dat, plot(Bth_wt, Bth_len))


#######################################################################################################################
#### Derived Measures #################################################################################################
#######################################################################################################################

# 1. Compute derived measures listed in the analysis plan
# 2. Check if same measures are consistent
# If applicable:
# 3. Check for outlier using the criteria established in the analysis plan (median +/- 5MAD)
# 4. Retain the outliers but create a variable with 0/1 (named flagvariablename) where 1 indicated
#    the outlier based on the analysis plan.


#### Preterm ####
# preterm is defined as gestational age < 37 weeks
dat$preterm <- ifelse(dat$BthGAwkscombi < 37, 1, 0)
table(dat$preterm)
dat$preterm <- factor(dat$preterm, labels = c("No", "Yes"))
addmargins(table(dat$preterm)) # available only for 686 people
attributes(dat$preterm)$label <- "pmmst - Preterm"


#### Low birth weight ####
# low birth weight is defined as weight < 2500 g
dat$lbw <- ifelse(dat$Bth_wt < 2500, 1, 0)
table(dat$lbw)
dat$lbw <- factor(dat$lbw, labels = c("No", "Yes"))
addmargins(table(dat$lbw))
attributes(dat$lbw)$label <- "pmmst - Low birth weight"


#### SGA ####
# Small for gestational age defined using Intergrowth data: José Villar et al. for the 
# International Fetal and Newborn Growth Consortium for the 21st Century (INTERGROWTH-21st)

# SGA was defined as birth weight <10th percentile for gestation age and sex using INTERGROWTH standards

# downlaod the INTERGROWTH Tables (centiles of weight by sex)
# and use them to find the 10th centiles. Table downloaded on June 20, 2017.

# Import INTERGROWTH Standard
datIntergrowth <- read_excel("Intergrowth10CentileBW_20160620.xlsx")
# change name of gestational age to BthGAwkscombi
names(datIntergrowth)[names(datIntergrowth) == "gestation(week)"] <- "BthGAwkscombi"
names(datIntergrowth)[names(datIntergrowth) == "10thCentileBoys"] <- "CentileBoys"
names(datIntergrowth)[names(datIntergrowth) == "10thCentileGirls"] <- "CentileGirls"
datIntergrowth <- subset(datIntergrowth, 
                         select = c("BthGAwkscombi", "CentileBoys", "CentileGirls"))
dat$BthGAwkscombi <- round(dat$BthGAwkscombi, 3)
datIntergrowth$BthGAwkscombi <- round(datIntergrowth$BthGAwkscombi, 3)
# join the two dataset
dat <- left_join(dat, datIntergrowth, by = 'BthGAwkscombi', type = "all")

# create birth weight in kg
dat$bwt_kg <- (dat$Bth_wt)/1000
summary(dat$bwt_kg)
attributes(dat$bwt_kg)$label <- "pmmst - Birth weight (kg)"

# create sga variable
for(i in 1:nrow(dat)){
  # set to missing all the boys centiles if subject is a girl
  dat$CentileBoys[i] <- ifelse(dat$Sex[i] == "F", NA, dat$CentileBoys[i])
  # set to missing all the girl centiles if subject is a boy
  dat$CentileGirls[i] <- ifelse(dat$Sex[i] == "M", NA, dat$CentileGirls[i])
  
}
# merge the centiles boys and girls in a unique variable
dat <- data.frame(dat)
dat$Centile <- dat$CentileBoys
dat$Centile[is.na(dat$CentileBoys)] <- dat$CentileGirls[is.na(dat$CentileBoys)]
# create sga
dat$sga <- ifelse(dat$bwt_kg < dat$Centile, 1, 0)
table(dat$sga)
dat$sga <- factor(dat$sga, labels = c("No", "Yes"))
addmargins(table(dat$sga))
sum(!is.na(dat$BthGAwkscombi) & !is.na(dat$Bth_wt)) # 254
# 273 people have gestational age and birth weight recorded.
dat[!is.na(dat$BthGAwkscombi) & !is.na(dat$Bth_wt) & is.na(dat$sga) & is.na(dat$Centile), 
    c("cSubjectID","Sex", "BthGAwkscombi", "CentileGirls", "CentileBoys","Bth_wt")]
# 18 were assigned as missing as their gestation is greater 
# than 42.714 days (maximum allowed in Intergrowth). One has GA less than 27. 
# most have gestation greater than 42.714 weeks  
dat[!is.na(dat$BthGAwkscombi) & dat$BthGAwkscombi > 42.714 & 
      is.na(dat$sga) & dat$Sex == "M", "Centile"] <- 3.21
dat[!is.na(dat$BthGAwkscombi) & dat$BthGAwkscombi > 42.714 & 
      is.na(dat$sga) & dat$Sex == "F", "Centile"] <- 3.03
dat[dat$cSubjectID == "EMP175Q", "Centile"] <- 0.74
# recompute small for gestational age
dat$sga <- ifelse(dat$bwt_kg < dat$Centile, 1, 0)
table(dat$sga)
dat$sga <- factor(dat$sga, labels = c("No SGA", "SGA"))
addmargins(table(dat$sga))
attributes(dat$sga)$label <- "pmmst - SGA (Intergrowth standard)"

#######################################################################################################################
#### Save Dataset #####################################################################################################
#######################################################################################################################

dat <- dat %>% dplyr::select(cSubjectID, Sex, MasterGroupNo, Bth_wt, Bth_len, hc, BthGAwkscombi, 
                            preterm, lbw, sga, flagbirthHC)
names(dat) <- c("serno", "sex", "allocation", "birthWeight", "birthLength", "birthHC", "ga", 
                "preterm", "lbw", "sga", "flagbirthHC")

# change newborn sex in Boy and Girl
table(dat$sex)
dat$sex[dat$sex=="M"] <- "Boy"
dat$sex[dat$sex=="F"] <- "Girl"

# add attributes
attributes(dat$serno)$label <- "pmmst - Subject ID"
attributes(dat$birthWeight)$label <- "pmmst - Birth weight (g)"
attributes(dat$birthLength)$label <- "pmmst - Birth length (cm)"
attributes(dat$birthHC)$label <- "pmmst - Birth HC (cm)"
attributes(dat$ga)$label <- "pmmst - Gestational Age (weeks)"
attributes(dat$allocation)$label <- "pmmst - allocation group"
attributes(dat$sex)$label <- "pmmst - newborn sex"

# save as csv, stata and SPSS file
write.table(dat, "./Cleaned Dataset/PMMSTBirthv1_20180305.csv", row.names = FALSE)
write_dta(dat, "./Cleaned Dataset/PMMSTBirthv1_20180305.dta")
write_sav(dat, "./Cleaned Dataset/PMMSTBirthv1_20180305.sav")
