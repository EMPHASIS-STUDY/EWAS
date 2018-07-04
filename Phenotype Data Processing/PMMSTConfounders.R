#####################################################################################################################
##### PMMST: Phenotype Data Processing Possible Confounders #########################################################
#####################################################################################################################

# 2018/02/06
# Chiara Di Gravio

library(readxl)
library(dplyr)
library(ggplot2)
library(haven)

dat <- read_excel("./Confounders/EMPHASIS_DSS_Biobank_Info_06022018_BS.xlsx", sheet = "Sheet1")

# take only the variable needed in the analysis
dat <- dat %>% dplyr::select(cSubjectID, DateOfBirth, mDoB, mAge, mEventType, mEduLevel, mEduYrs, 
                             mEndDate, fDoB, fAge, fEventType,
                             fEndDate, fEduLevel, fEduYrs, dDateMeasurement, EndDate, TerminatingEventType)

# take only those children that are included in the EMHASIS study and have data for clinic visit 1
datCV1 <- read_excel("./PMMST Outcome Data/Emphasis_CV1_AnthropsBP_170630.xlsx", 
                  sheet = "AnthropBP", na = "NULL", col_names = TRUE)
datCV1 <- datCV1 %>% dplyr::select(cISubjectID, dVisitDate)

# change the subject ID variable name to match the one in the confounders dataset info variables.
names(datCV1)[names(datCV1) == "cISubjectID"] <- "cSubjectID"
# Transform every ID in upper case
datCV1$cSubjectID <- toupper(datCV1$cSubjectID)

# EMP180Q has only visit 2 data so we need to add information on this child separately
datCV2 <- anti_join(dat, datCV1) %>% dplyr::filter(cSubjectID == "EMP180Q")
# merge the two dataset together
dat <- inner_join(datCV1, dat)
dat <- dat[!duplicated(dat), ]
# 298 as expected
dat <- rbind(dat, datCV2)
datCV2
# add an extra row to the dataset
dat[nrow(dat) + 1,] <- c(datCV2$cSubjectID, NA, "2007-12-23", "1965-07-02", datCV2$mAge, datCV2$mEventType, 
                         datCV2$mEduLevel, datCV2$mEduYrs, NA, "1967-07-02", datCV2$fAge, datCV2$fEventType, 
                         NA, datCV2$fEduLevel, datCV2$fEduYrs,  "2016-10-27", NA, datCV2$TerminatingEventType)
rm(datCV1, datCV2)

# measurement date should be clinic visit 2 check
temp <- data.frame(s = dat$cSubjectID, v1 = dat$dVisitDate, v2 = dat$dDateMeasurement)
temp$diff <- as.numeric(temp$v2 - temp$v1)
summary(temp$diff)
# check those who have NA
temp[is.na(temp$diff),]
# those are the three who only have v1 and the one who only have visit 2
rm(temp)

# add birth measurements
datBirth <- read_excel("./PMMST Outcome Data/PMMST_pdata_redux_170718.xlsx", na = "NA", col_names = TRUE)
# take maternal BMI 
datBirth <- datBirth %>% dplyr::select(cSubjectID, MotherBMI)
dat <- left_join(dat, datBirth)
rm(datBirth)

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

#### Maternal Age ####

# maternal age at clinic visit 2
summary(dat$mAge)
sd(dat$mAge, na.rm = T)
names(dat)[names(dat) == "mAge"] <- "mageCV2"
dat$mageCV2 <- as.numeric(dat$mageCV2)
attributes(dat$mageCV2)$label <- "pmmst - maternal age at CV2 (years)" 
hist(dat$mageCV2, main = "Maternal age at CV2", xlab = "Age (years)")
     
# compute maternal age (in years) at clinic visit 1.
dat$mageCV1 <- as.numeric((dat$dVisitDate - dat$mDoB)/365.25)
summary(dat$mageCV1)
sd(dat$mageCV1, na.rm = T)
attributes(dat$mageCV1)$label <- "pmmst - maternal age at CV1 (years)" 
hist(dat$mageCV1, main = "Maternal age at CV1", xlab = "Age (years)")
# compare this with the one in the confouders dataset
summary(dat$mageCV2)
par(mfrow = c(2,1))
hist(dat$mageCV1)
hist(dat$mageCV2)
par(mfrow = c(1,1))
# check difference between age at visit 1 and 2 is positive
diff <- dat$mageCV2 - dat$mageCV1 
sum(diff < 0, na.rm = T)
rm(diff)

# compute maternal age at delivery
dat$mageDel <- as.numeric((dat$DateOfBirth - dat$mDoB)/365.25)
summary(dat$mageDel) # no NA
sd(dat$mageDel)
attributes(dat$mageDel)$label <- "pmmst - maternal age at delivery (years)" 
hist(dat$mageDel, main = "Maternal age at delivery", xlab = "Age (years)")

# the difference between maternal age at delivery and that at visit 1 should give information on the child age
dat$cageCV1 <- as.numeric((dat$mageCV1 - dat$mageDel))
summary(dat$cageCV1)
sd(dat$cageCV1)
attributes(dat$cageCV1)$label <- "pmmst - child age at CV1 (years)" 
hist(dat$cageCV1, main = "Children age at CV1", xlab = "Age (years)")


#### Maternal education level ####
dat$mEduLevel <- factor(dat$mEduLevel)
addmargins(table(dat$mEduLevel))
# CU (COLLEGE/UNIVERSITY), JS (JUNIOR SECONDARY), LB (LOWER BASIC), MA (MADRASSA), NO (NONE), NS (NURSERY SCHOOL),
# OT (OTHER), QE (QURANIC EDUCATION), SS (SENIOR SECONDARY), UB (UPPER BASIC)
sum(is.na(dat$mEduLevel))
# 44 women do not have any information on education level
temp <- dat[is.na(dat$mEduLevel), ]
table(temp$mEventType)
rm(temp)
# hows the parents residency status; EXT for OutMigration, DTH for death and Blank for resident. 
# The enddates indicate dates for these events.
# There are two maternal deaths


#### Create a dummy variable to indicate maternal death ####
dat$mdeath <- ifelse(dat$mEventType == "DTH", 1, 0)
dat$mdeath[is.na(dat$mdeath)] <- 0
table(dat$mdeath)
dat$mdeath <- factor(dat$mdeath, labels = c("No", "Yes"))
attributes(dat$mdeath)$label <- "pmmst - maternal death?" 
# 4 mothers are death in the study check which one and when
temp <- dat[dat$mdeath == "Yes" & !is.na(dat$mdeath),]
temp$diff <- temp$dVisitDate - temp$mEndDate
# for all of those mothers that died before CV1, we need to set age at visit 1 and 2 as missing
dat$mageCV1[(dat$cSubjectID %in% subset(temp, diff > 0)$cSubjectID)] <- NA
dat$mageCV2[(dat$cSubjectID %in% subset(temp, diff > 0)$cSubjectID)] <- NA
rm(temp)
# of the 4 who died only 1 died after the visits 1 and 2 EMP046J
summary(dat$mageCV1)
summary(dat$mageCV2)


#### Maternal BMI ####
summary(dat$MotherBMI)
sd(dat$MotherBMI)
attributes(dat$MotherBMI)$label <- "pmmst - maternal BMI" 
hist(dat$MotherBMI)


#### Paternal Age ####
summary(dat$fAge)
sd(dat$fAge, na.rm = T)
# 32 NA's for age  
dat[is.na(dat$fAge),]
names(dat)[names(dat) == "fAge"] <- "fageCV2"
dat$fageCV2 <- as.numeric(dat$fageCV2)
attributes(dat$fageCV2)$label <- "pmmst - paternal age at CV2 (years)" 
hist(dat$fageCV2, main = "Paternal age at CV2", xlab = "Age (years)")


# compute paternal age (in years) at clinic visit 1.
dat$fageCV1 <- as.numeric((dat$dVisitDate - dat$fDoB)/365.25)
summary(dat$fageCV1) # 29 NA
sd(dat$fageCV1, na.rm = T)
attributes(dat$fageCV1)$label <- "pmmst - paternal age at CV1 (years)" 
hist(dat$fageCV1, main = "Paternal age at CV1", xlab = "Age (years)")
# check those with high paternal age
dat[dat$fageCV1 > 70 & !(is.na(dat$fageCV1)), ]

# check visit 1 and 2 ages
sum(dat$fageCV2 - dat$fageCV1 < 0, na.rm = T)

# compute paternal age at delivery
dat$fageDel <- as.numeric((dat$DateOfBirth - dat$fDoB)/365.25)
summary(dat$fageDel) # 29 NA
sd(dat$fageDel, na.rm = T)
attributes(dat$fageDel)$label <- "pmmst - paternal age at delivery (years)" 
hist(dat$fageDel, main = "Paternal age at delivery", xlab = "Age (years)")

# check visit 1 and delivery ages
sum(dat$fageCV1 - dat$fageDel < 0, na.rm = T)

#### Paternal education level ####
dat$fEduLevel <- factor(dat$fEduLevel)
addmargins(table(dat$fEduLevel))
# CU (COLLEGE/UNIVERSITY), JS (JUNIOR SECONDARY), LB (LOWER BASIC), MA (MADRASSA), NO (NONE), NS (NURSERY SCHOOL),
# OT (OTHER), QE (QURANIC EDUCATION), SS (SENIOR SECONDARY), UB (UPPER BASIC)
sum(is.na(dat$fEduLevel))
# 81 men do not have any information on education level
temp <- dat[is.na(dat$fEduLevel), ]
table(temp$fEventType)
rm(temp)
# There are 23 paternal deaths among those with no education information


#### Create a dummy variable to indicate paternal death ####
dat$fdeath <- ifelse(dat$fEventType == "DTH", 1, 0)
dat$fdeath[is.na(dat$fdeath)] <- 0
table(dat$fdeath)
dat$fdeath <- factor(dat$fdeath, labels = c("No", "Yes"))
attributes(dat$fdeath)$label <- "pmmst - paternal death?" 
# 36 fathers are death in the study check which one and when
temp <- dat[dat$fdeath == "Yes" & !is.na(dat$fdeath),]
temp$diff <- temp$dVisitDate - temp$fEndDate
# for all of those fathers that died before CV1 (32), we need to set age at visit 1 and 2 as missing
dat$fageCV1[(dat$cSubjectID %in% subset(temp, diff > 0)$cSubjectID)] <- NA
dat$fageCV2[(dat$cSubjectID %in% subset(temp, diff > 0)$cSubjectID)] <- NA
# check those who died after clinic visit 1 to see whether they where alive at CV2
temp2 <- temp %>% dplyr::filter(diff < 0)
rm(temp2) # they all died after
# check if they have died before child's birth
temp$diff <- temp$fEndDate - temp$DateOfBirth
temp <- temp %>% dplyr::filter(diff < 0)
# these should have age at delivery as missing
dat$fageDel[(dat$cSubjectID %in% temp$cSubjectID)] <- NA
rm(temp)


#### Create a variable looking at whether the child lives on the coast or internal area? ####

table(dat$TerminatingEventType) # 90 children should live on the coast.
dat$ccoast <- ifelse(dat$TerminatingEventType == "EXT", 1, 0)
dat$ccoast[is.na(dat$ccoast)] <- 0
dat$ccoast <- factor(dat$ccoast, labels = c("No", "Yes"))
table(dat$ccoast)

# check whether they moved to the coast before birth
temp <- dat[dat$ccoast == "Yes",]
temp$diff <- temp$EndDate - temp$DateOfBirth
# no one moved before birth
# check if they moved before cv2
temp <- dat[dat$ccoast == "Yes",]
temp$diff <- temp$EndDate - temp$dDateMeasurement
temp <- temp %>% dplyr::filter(diff > 0)

# 7 move to the coast after visit 2. These should have a no in the ccosat variable
dat$ccoast[(dat$cSubjectID %in% temp$cSubjectID)] <- "No"
table(dat$ccoast)
rm(temp)

#### Add other information on child medium of study ####

datAtlantis <- read_excel("./PMMST Outcome Data/Emphasis_CV2_Cognitive_170630.xlsx", 
                          sheet = "Atlantis", na = "null", skip = 1, col_names = TRUE)

# change every ID number in upper case
datAtlantis$cISubjectID <- toupper(datAtlantis$cISubjectID)
str(datAtlantis)
# transform school attended in a factor. Levels are as following:
# 1. English
# 2. Arabic
# 3. Both
# 4. None
addmargins(table(datAtlantis$cSchoolAttended))
datAtlantis$cSchoolAttended <- factor(datAtlantis$cSchoolAttended,
                                      labels=c("English","Arabic","Both", "None"))
addmargins(table(datAtlantis$cSchoolAttended))
# take only the variable needed for the analysis
datAtlantis <- datAtlantis %>% dplyr::select(cISubjectID, cSchoolAttended)
names(datAtlantis)[names(datAtlantis) == "cISubjectID"] <- "cSubjectID"
# merge with the other confounders data
dat <- left_join(dat, datAtlantis)
rm(datAtlantis)

########################################################################################################################
#### Save Dataset ######################################################################################################
########################################################################################################################

dat <- dat %>% arrange(-desc(cSubjectID)) %>% 
  dplyr::select(cSubjectID, dVisitDate, dDateMeasurement, mDoB, mageDel, mageCV1, mageCV2, mEduLevel, MotherBMI, mdeath,
                fDoB, fageDel, fageCV1, fageCV2, fEduLevel, fdeath, ccoast, cSchoolAttended) 
names(dat) <- c("serno", "cdateCV1", "cdateCV2", "mDoB", "mageDel", "mageCV1", "mageCV2", "mEduLevel", "mbmi", "mdeath",
                "fDoB", "fageDel", "fageCV1", "fageCV2", "fEduLevel", "fdeath", "ccoast", "cschool")
str(dat)
# delete times from date
dat$cdateCV2 <- as.POSIXct(format(dat$cdateCV2, format='%Y-%m-%d'))
dat$mDoB <- as.POSIXct(format(dat$mDoB, format='%Y-%m-%d'))
dat$fDoB <- as.POSIXct(format(dat$fDoB, format='%Y-%m-%d'))

# label variables
attributes(dat$serno)$label <- "pmmst - Subject ID"
attributes(dat$cdateCV1)$label <- "pmmst - clinic visit 1 (date)"
attributes(dat$cdateCV2)$label <- "pmmst - clinic visit 2 (date)"
attributes(dat$mDoB)$label <- "pmmst - maternal date of birth"
attributes(dat$mEduLevel)$label <- "pmmst - maternal education level"
attributes(dat$mbmi)$label <- "pmmst - maternal BMI"
attributes(dat$fDoB)$label <- "pmmst - paternal date of birth"
attributes(dat$fEduLevel)$label <- "pmmst - paternal education level"
attributes(dat$ccoast)$label <- "pmmst - child on the coast?"
attributes(dat$cschool)$label <- "pmmst - language of school attended by the child"

# save as STATA, SPSS and csv
write.table(dat, "./Cleaned Dataset/PMMSTConfoundersv1_20180215.csv", row.names = FALSE)
write_dta(dat, "./Cleaned Dataset/PMMSTConfoundersv1_20180215.dta")
write_sav(dat, "./Cleaned Dataset/PMMSTConfoundersv1_20180215.sav")
