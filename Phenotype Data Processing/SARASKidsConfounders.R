#######################################################################################################################
##### SARAS Kids: Phenotype Data Processing Possible Confouners #######################################################
#######################################################################################################################

# March 14, 2018

# load packages
library(haven)
library(dplyr)
library(ggplot2)

# the first file emphasis01.dta was sent by Patsy. After queries to both Harshad and Patsy. Patsy recreated the 
# paerental data file. One file for the mothers (sarasmother01.dta), one file for fathers (sarasfather01.dta).
# The data was created on March 8, 2018.

# import maternal data
datM <- read_dta("./Confounders/sarasmother01.dta")
# check for duplicates
length(unique(datM$serno))
# import paternal data
datF <- read_dta("./Confounders/sarasfather01.dta")
# check for duplicates
length(unique(datF$serno))
# merge
dat <- left_join(datM, datF)
rm(datM, datF)

# add information on allocation group and child age. PP analysis only
datList <- read_dta("SARASKidsBloodCountPP.dta")
datList <- datList %>% dplyr::select(serno, allocation, age)
# check if all sernos in datList are in dat
datList$allocation <- as_factor(datList$allocation)
# merge the two dataset
dat <- left_join(datList, dat)
rm(datList)
dim(dat) # 709 

# delete some of not useful variables from the dataset
drop <- c("mmdatereg", "mmcno", "mmmdmsdy", "mmmdmspc", "mmoccspc", "mmtbcodly", "mmtbco", "mmmfccigrt", "mmbedi",
          "mmganja", "mmmashri", "mmsnuff", "mmcwgtbco", "mmpanzrda", "mmgutka", "mmtbcoth", "mmtbcothsp", "mmaudit1", 
          "mmaudit2", "mmaudit3", "mmaudit4", "mmaudit5", "mmaudit6", "mmaudit7", "mmaudit8", "mmaudit9", "mmaudit10",
          "mmaudits", "mmdateant", "mmpaname", "mmhbp", "mmhbpy", "mmhbpm", "mmhbpme", "mmhbpmed", "mmtp2db", "mmtp2dby", 
          "mmtp2dbm", "mmtp2dbme", "mmtp2dbmed", "mmhrtds", "mmhrtdsy", "mmhrtdsm", "mmhrtdsme", "mmhrtdsmed", 
          "mmstrok", "mmstroky", "mmstrokm", "mmstrokme", "mmstrokmed", "mmhicho", "mmhichoy", "mmhichom", 
          "mmhichome", "mmhichomed", "mmcancr", "mmcancry", "mmcancrm", "mmcancrme", "mmcancrmed", "mmchlgdis", 
          "mmchlgdisy", "mmchlgdism", "mmchlgdisme", "mmchlgdismed", "mmchkddis", "mmchkddisy", "mmchkddism", 
          "mmchkddisme", "mmchkddismed", "mmothspc", "mmothspcc", "mmothspcy", "mmothspcm", "mmothspcme", "mmothspcmed",
          "mmcursup", "mmpnam1", "mmpnam1d", "mmpnam1o", "mmpnam1l", "mmpnam2", "mmpnam2d", "mmpnam2o", "mmpnam2l",
          "mmpnam3", "mmpnam3d", "mmpnam3o", "mmpnam3l", "mmpnam4", "mmpnam4d", "mmpnam4o", "mmpnam4l", "mmdiasaras",
          "mmcurpreg", "mmnurse", "mmpulse1", "mmpulse2", "mmpulse3", "mmdategrp", 
          "mmrtrdg1", "mmrtrdg2", "mmrtrdg3", "mmltrdg1", "mmltrdg2", "mmltrdg3", "mmhand", "datesli", 
          "mfcno", "mfcname", "mfmdmsdy", "mfmdmspc", "mfoccspc", "mftbco", "mftbcodly", "mfmfccigrt", "mfbedi",
          "mfganja", "mfmashri", "mfsnuff", "mfcwgtbco", "mfpanzrda", "mfgutka", "mftbcoth", "mftbcothsp", "mfaudit1",
          "mfaudit2", "mfaudit3", "mfaudit4", "mfaudit5", "mfaudit6", "mfaudit7", "mfaudit9", "mfaudit10", "mfaudits",
          "mfmdate", "mfpaname", "mfhbp", "mfhbpy", "mfhbpm", "mfhbpme", "mfhbpmed", "mftp2db", "mftp2dby", 
          "mftp2dbm", "mftp2dbme", "mftp2dbmed", "mfhrtds", "mfhrtdsy", "mfhrtdsm", "mfhrtdsme", "mfhrtdsmed", 
          "mfstrok", "mfstroky", "mfstrokm", "mfstrokme", "mfstrokmed", "mfhicho", "mfhichoy", "mfhichom", 
          "mfhichome", "mfhichomed", "mfcancr", "mfcancry", "mfcancrm", "mfcancrme", "mfcancrmed", "mfchlgdis",
          "mfchlgdisy", "mfchlgdism", "mfchlgdisme", "mfchlgdismed", "mfchkddis", "mfchkddisy", "mfchkddism", 
          "mfchkddisme", "mfchkddismed", "mfothsp", "mfothspc", "mfothspcy", "mfothspcm", "mfothspcme", "mfothspcmed",
          "mfcursup", "mfpnam1", "mfpnam1d", "mfpnam1o", "mfpnam1l", "mfpnam2", "mfpnam2d", "mfpnam2o", "mfpnam2l",
          "mfpnam3", "mfpnam3d", "mfpnam3o", "mfpnam3l", "mfpnam4", "mfpnam4d", "mfpnam4o", "mfpnam4l", "mfdateant",
          "mfnurse", "mfpulse1", "mfpulse2", "mfpulse3", "mfdategrp", "mfrtrdg1", "mfrtrdg2", "mfrtrdg3", "mfltrdg1", 
          "mfltrdg2", "mfltrdg3", "mfhand")
dat <- dat[,!(names(dat) %in% drop)]
str(dat)
rm(drop)

#########################################################################################################################
#### Data Cleaning ######################################################################################################
#########################################################################################################################

# transform relevant variable in factors
dat$allocation <- as_factor(dat$allocation)
# mother alive or death?
# add information based on what collected in emphasis01.dta
# sernos 1407, 3846, 4075 and 4743 died. Add a variable with the information
dat$mdeath <- 0
dat[dat$serno %in% c(1407, 3846, 4075, 4743), "mdeath"] <- 1
table(dat$mdeath)
dat$mdeath <- factor(dat$mdeath, labels = c("No", "Yes"))
# father alive or death
# add information based on what collected in emphasis01.dta
dat$fdeath <- 0
dat[dat$serno %in% c(589, 1762, 2947, 3303, 3618, 4371, 5307), "fdeath"] <- 1
table(dat$fdeath)
dat$fdeath <- factor(dat$fdeath, labels = c("No", "Yes"))

#######################################################################################################################
#### Maternal Data ####################################################################################################
#######################################################################################################################

# maternal education level 
table(dat$mmedulev)
# create 3 category: primary or less, secondary, graduate and above
dat$mmedulev[dat$mmedulev == 5] <- 0
dat$mmedulev[dat$mmedulev == 6] <- 0
dat$mmedulev[dat$mmedulev == 7] <- 0
dat$mmedulev[dat$mmedulev == 1] <- 2
dat$mmedulev[dat$mmedulev == 3] <- 1
dat$mmedulev[dat$mmedulev == 4] <- 1
addmargins(table(dat$mmedulev))
# 6 have maternal education level missing. Check with registration data if there are available information
temp <- dat %>% dplyr::filter(is.na(mmedulev))
# add only information for the 2 that are still alive at time of the study (4 are death).
dat[dat$serno %in% c(499, 1443),"mmedulev"] <- 1
rm(temp)
# set to missing the 4 women who died
dat$mmedulev <- ifelse(dat$mdeath == "Yes", NA, dat$mmedulev)
dat$mmedulev <- factor(dat$mmedulev, labels = c("Primary or less", "Secondary", "Graduate and above"))
addmargins(table(dat$mmedulev))


# maternal age
summary(dat$mmagey)
# there are no all the values of maternal age. 
# create a variable by looking at difference between date of visit and dob
dat$mage <-  with(dat, (mmdate - mmdob)/365.25) 
# for those who do not have age, we take the registration age, add the time from registration to delivery and add 
# the child age.
# first add age at registration and delivery date
agereg <- read_dta("T:/Chiara/EMPHASIS/MMNPanalysis.dta")
# keep only variables that are needed
agereg <- agereg %>% dplyr::select(serno, deldate, visitdate, calcage)
# merge with cuttern dataset
dat <- inner_join(dat, agereg)
rm(agereg)
# difference in year between visit date and delivery date:
dat$diff <- with(dat, deldate - visitdate)/365.25
# add this difference as well as child age to the registration calcage
dat$calculatedage <- as.numeric(dat$calcage + dat$diff + dat$age)
dat$mage[is.na(dat$mage)] <- dat$calculatedage[is.na(dat$mage)]
dat$mage <- as.numeric(dat$mage)
hist(dat$mage)
dat[dat$mage < 20 & !is.na(dat$mage), c("serno", "mmdob","mmantdate", "mmagey", "mage")]
# 1 mothers have age < 20. This mother has a wrong date of birth recorded. 
# based on Harshad's email (March 6, 2018), serno 3898 DOB is 02-03-1972. 
dat$mmdob[dat$serno == 3898] <- "1977-03-02"
dat$mage[dat$serno == 3898] <- (dat$mmdate[dat$serno == 3898] - dat$mmdob[dat$serno == 3898])/365.25
# remove mmagey and calculated age
dat <- dat[, !(names(dat) %in% c("mmagey", "calculatedage"))]
# sumamry statistics of age
summary(dat$mage)
hist(dat$mage)
# set to missing the 4 womenthat are deatf
dat$mage <- ifelse(dat$mdeath == "Yes", NA, dat$mage)
summary(dat$mage)



# maternal occupation
addmargins(table(dat$mmocc))
# those with 0 or with missing occupation are currently unemployed. Transform NA to 0
dat$mmocc[is.na(dat$mmocc)] <- 0
# recategorise women occupation based on the Potdar 2014 paper
dat$mmocc[dat$mmocc == 1 | dat$mmocc == 2] <- 8 
dat$mmocc[dat$mmocc == 5 | dat$mmocc == 6] <- 1
dat$mmocc[dat$mmocc == 3 | dat$mmocc == 4] <- 2
dat$mmocc[dat$mmocc == 8] <- 3
addmargins(table(dat$mmocc))
# set to missing the 4 women who died
dat$mmocc <- ifelse(dat$mdeath == "Yes", NA, dat$mmocc)
dat$mmocc <- factor(dat$mmocc, labels = c("Not working", "unskilled/semi-skilled", 
                                          "skilled/self-employed", "semi-professional/professional"))
addmargins(table(dat$mmocc))



# maternal anthropometry

# weight and height
summary(dat$mmweight)
hist(dat$mmweight)
# check weight greater than 100
dat[dat$mmweight > 100 & !is.na(dat$mmweight), c("serno", "mage", "mmweight")]
# check the minimum weight
dat[dat$mmweight == min(dat$mmweight, na.rm = T) & !is.na(dat$mmweight), 
    c("serno", "mmweight", "mmstdhgt")]
# check weight with height
with(dat, plot(mmstdhgt, mmweight))

# check those with outlier weight/height combination
dat[dat$mmweight > 60 & dat$mmstdhgt < 120 & !is.na(dat$mmweight), 
    c("serno", "mmweight", "mmstdhgt", "mmsithgt", "mmmuac", "mmwaist", "mmtricp1")] 
# standing height for 4513 is too low: remove the point after checking with Harshad:
dat[dat$serno == 4513, "mmstdhgt"] <- NA
dat[dat$mmweight < 40 & dat$mmstdhgt < 140 & !is.na(dat$mmweight), 
    c("serno", "mmweight", "mmstdhgt", "mmsithgt", "mmmuac", "mmwaist", "mmtricp1")]
# those were checked with Harshad and are OK

# check that sitting height is less than standing height
diff <- dat$mmstdhgt - dat$mmsithgt
summary(diff)
sum(diff < 0, na.rm = T)
rm(diff)
# check sitting height
hist(dat$mmsithgt)
dat[dat$mmsithgt > 130 & !is.na(dat$mmsithgt), 
    c("serno", "mmweight", "mmstdhgt", "mmsithgt", "mmmuac", "mmwaist", "mmtricp1")]
# serno 4371 has a wrong value of sitting height (checked with Harshad). Delete
dat[dat$serno == 4371, "mmsithgt"] <- NA
# take off the height of the stal for sitting height (31.54 cm)
dat$mmsithgt <- dat$mmsithgt - 31.54
hist(dat$mmsithgt)

# compute leg length
dat$mleg <- with(dat, mmstdhgt - mmsithgt)
summary(dat$mleg)
hist(dat$mleg)
# find the minimum of leg length
dat[which.min(dat$mleg), c("serno", "mmweight", "mmstdhgt", "mmsithgt", "mmmuac", "mmwaist", "mleg")]

# compute maternal BMI
dat$mbmi <- with(dat, mmweight/(0.01*mmstdhgt)^2)
summary(dat$mbmi)
hist(dat$mbmi)

# head circumference
summary(dat$mmhdcir)
hist(dat$mmhdcir)
dat %>% filter(mdeath == "Yes") %>% 
  select(serno, mmweight, mmstdhgt, mmsithgt, mmmuac, mmwaist, mmhdcir)

# MUAC
summary(dat$mmmuac)
hist(dat$mmmuac)
# check those MUAC greater than 40 cm
dat[dat$mmmuac > 40 & !is.na(dat$mmmuac), 
    c("serno", "mmweight", "mmstdhgt", "mmsithgt", "mmmuac", "mmwaist", "mmtricp1")]
dat %>% filter(serno == 2829) %>% 
  select(serno, mmweight, mmstdhgt, mmsithgt, mmmuac, mmwaist, mmhdcir)
# 2829 MUAC was too high (set to missing)

# waist
summary(dat$mmwaist)
hist(dat$mmwaist)

# hip
summary(dat$mmhip)
hist(dat$mmhip)

# triceps, biceps, subscapular and suprialiac skinfolds
summary(dat$mmtricp1)
summary(dat$mmtricp2)
summary(dat$mmtricp3)
# check women has similar measure
diff12 <- with(dat, mmtricp1 - mmtricp2)
summary(diff12)
rm(diff12)
diff13 <- with(dat, mmtricp1 - mmtricp3)
summary(diff13)
rm(diff13)
diff23 <- with(dat, mmtricp2 - mmtricp3)
summary(diff23)
rm(diff23)
# average triceps measures
dat$avgmTric <- rowMeans(dat[, c("mmtricp1", "mmtricp2", "mmtricp3")])
summary(dat$avgmTric)
hist(dat$avgmTric)

summary(dat$mmbicep1)
summary(dat$mmbicep2)
summary(dat$mmbicep3)
# check women has similar measure
diff12 <- with(dat, mmbicep1 - mmbicep2)
summary(diff12)
rm(diff12)
diff13 <- with(dat, mmbicep1 - mmbicep3)
summary(diff13)
rm(diff13)
diff23 <- with(dat, mmbicep2 - mmbicep3)
summary(diff23)
rm(diff23)
# average biceps measures
dat$avgmBic <- rowMeans(dat[, c("mmbicep1", "mmbicep2", "mmbicep3")])
summary(dat$avgmBic)
hist(dat$avgmBic)
dat[dat$avgmBic > 30 & !is.na(dat$avgmBic), 
    c("serno", "mmweight", "mmstdhgt", "mmsithgt", "mmmuac", "mmwaist", "mmhip", "mmtricp1",
      "avgmBic")]

summary(dat$mmsubscp1)
summary(dat$mmsubscp2)
summary(dat$mmsubscp3)
# check women has similar measure
diff12 <- with(dat, mmsubscp1 - mmsubscp2)
summary(diff12)
dat[diff12 < -2, 
    c("serno", "mmweight",
      "mmsubscp1", "mmsubscp2", "mmsubscp3")]
# after talking with Patsy (08/02/2017) set to missing those values that are different
dat[dat$serno == 1103, "mmsubscp1"] <- NA
dat[dat$serno == 1542, "mmsubscp2"] <- NA
rm(diff12)
diff13 <- with(dat, mmsubscp1 - mmsubscp3)
summary(diff13)
rm(diff13)
diff23 <- with(dat, mmsubscp2 - mmsubscp3)
summary(diff23)
rm(diff23)
# average subscapular measures
dat$avgmSubsc <- rowMeans(dat[, c("mmsubscp1", "mmsubscp2", "mmsubscp3")], na.rm = T)
summary(dat$avgmSubsc)
hist(dat$avgmSubsc)

summary(dat$mmsupri1)
summary(dat$mmsupri2)
summary(dat$mmsupri3)
# check women has similar measure
diff12 <- with(dat, mmsupri1 - mmsupri2)
summary(diff12)
rm(diff12)
diff13 <- with(dat, mmsupri1 - mmsupri3)
summary(diff13)
dat[diff13 < -2, 
    c("serno","mmweight", "mmstdhgt", "mmsithgt", "mmmuac",
      "mmsupri1", "mmsupri2", "mmsupri3")]
# suprialiac measure for 150 was 9.6
dat[dat$serno == 150, "mmsupri3"] <- 9.6
# after talking with Patsy (08/02/2017) set to missing those values that are different
dat[diff13 > 2, 
    c("serno", "mmweight", "mmstdhgt", "mmsithgt", "mmmuac",
      "mmsupri1", "mmsupri2", "mmsupri3")]
dat[dat$serno == 17, "mmsupri3"] <- NA
rm(diff13)
diff23 <- with(dat, mmsupri2 - mmsupri3)
summary(diff23)
rm(diff23)
# average suprailiac measures
dat$avgmSup <- rowMeans(dat[, c("mmsupri1", "mmsupri2", "mmsupri3")], na.rm = T)
summary(dat$avgmSup)
hist(dat$avgmSup)


# blood pressure

# systolic
summary(dat$mmsbp1)
summary(dat$mmsbp2)
summary(dat$mmsbp3)
# average blood pressure values
dat$avgmsys <- rowMeans(dat[, c("mmsbp1", "mmsbp2", "mmsbp3")], na.rm = T)
summary(dat$avgmsys)
hist(dat$avgmsys)

# diastolic
summary(dat$mmdbp1)
summary(dat$mmdbp2)
summary(dat$mmdbp3)
# average blood pressure values
dat$avgmdia <- rowMeans(dat[, c("mmdbp1", "mmdbp2", "mmdbp3")], na.rm = T)
summary(dat$avgmdia)
hist(dat$avgmdia)


#######################################################################################################################
#### Paternal Data ####################################################################################################
#######################################################################################################################

# father age
summary(dat$mfagey)
# there are no all the values of maternal age. 
# create a variable by looking at difference between date of visit and dob
dat$fage <-  with(dat, (mfdate - mfdob)/365.25) 
temp <- data.frame(calcage = dat$fage, age = dat$mfagey)
rm(temp)
dat$fage[is.na(dat$fage)] <- dat$mfagey[is.na(dat$fage)]
dat$fage <- as.numeric(dat$fage)
summary(dat$fage)
hist(dat$fage)
# 83 are missing. No available information
# look at those with father death and check whether they are missing.
dat %>% filter(fdeath == "Yes") %>% select(serno, fage)
# all 7 with father death have missing age
# remove mfagey
dat <- dat[, !(names(dat) %in% c("mfagey"))]

# education level
# create 3 category: primary or less, secondary, graduate and above
dat$mfedulev[dat$mfedulev == 5] <- 0
dat$mfedulev[dat$mfedulev == 6] <- 0
dat$mfedulev[dat$mfedulev == 7] <- 0
dat$mfedulev[dat$mfedulev == 1] <- 2
dat$mfedulev[dat$mfedulev == 3] <- 1
dat$mfedulev[dat$mfedulev == 4] <- 1
addmargins(table(dat$mfedulev))
# add the missing available from woman registration data
temp <- dat[is.na(dat$mfedulev), c("serno", "mfedulev")]
# change missing based on those observed at registration in the MMNP
dat$mfedulev[dat$serno %in% c(121, 1370, 1548, 3245, 3465, 4060)] <- 0
dat$mfedulev[dat$serno %in% c(777, 2970, 3199, 4371)] <- 2
dat$mfedulev[is.na(dat$mfedulev)] <- 1
# set to missing the 7 men who died
dat$mfedulev <- ifelse(dat$fdeath == "Yes", NA, dat$mfedulev)
dat$mfedulev <- factor(dat$mfedulev, labels = c("Primary or less", "Secondary", 
                                                "Graduate and above"))
addmargins(table(dat$mfedulev))

# paternal occupation
addmargins(table(dat$mfocc))
# those with 0 or with missing occupation are currently unemployed. Transform NA to 0
dat$mfocc[is.na(dat$mfocc)] <- 0
# recategorise men occupation based on the Potdar 2014 paper
dat$mfocc[dat$mfocc == 1 | dat$mfocc == 2] <- 8 
dat$mfocc[dat$mfocc == 5 | dat$mfocc == 6] <- 1
dat$mfocc[dat$mfocc == 3 | dat$mfocc == 4] <- 2
dat$mfocc[dat$mfocc == 8] <- 3
addmargins(table(dat$mfocc))
# set to missing the 7 men who died
dat$mfocc <- ifelse(dat$fdeath == "Yes", NA, dat$mfocc)
dat$mfocc <- factor(dat$mfocc, labels = c("Not working", "unskilled/semi-skilled", 
                                          "skilled/self-employed", "semi-professional/professional"))
addmargins(table(dat$mfocc))


# paternal anthropometry

# weight and height
summary(dat$mfwt)
hist(dat$mfwt)
dat[dat$serno == 4007, "mfwt"] <- NA # weight too low compared to other measures
# check weight with height
with(dat, plot(mfstdht, mfwt))
# check those with outlier weight/height combination
dat[dat$mfwt < 60 & dat$mfstdht < 120 & !is.na(dat$mfwt), 
    c("serno", "mfwt", "mfstdht", "mfsitht", "mfmuac", "mfwaist", "mftricp1")]
# serno 1707 height is 161.5 (based on Harshad's email). Need to correct
dat[dat$serno == 1707, "mfstdht"] <- 161.5

# check that sitting height is less than standing height (leg length)
diff <- dat$mfstdht - dat$mfsitht
summary(diff)
sum(diff < 0, na.rm = T)
rm(diff)

# sitting height
hist(dat$mfsitht)
dat[dat$mfsitht > 140 & !is.na(dat$mfsitht), 
    c("serno", "mfwt", "mfstdht", "mfsitht", "mfmuac", "mfwaist", "mftricp1")]
# correct sitting height for serno 1109 is 115 as confirmed by Harshad's email
dat[dat$serno == 1109, "mfsitht"] <- 115
dat[dat$mfsitht < 110 & !is.na(dat$mfsitht), 
    c("serno", "mfwt", "mfstdht", "mfsitht", "mfmuac", "mfwaist", "mftricp1")]
# take off the height of the stal for sitting height (31.54 cm)
dat$mfsitht <- dat$mfsitht - 31.54
hist(dat$mfsitht)
summary(dat$mfsitht)


# compute leg length
dat$fleg <- with(dat, mfstdht - mfsitht)
summary(dat$fleg)
hist(dat$fleg)

# compute paternal BMI
dat$fbmi <- with(dat, mfwt/(0.01*mfstdht)^2)
summary(dat$fbmi)
hist(dat$fbmi)

# head circumference
summary(dat$mfhdcir)
hist(dat$mfhdcir)

# MUAC
summary(dat$mfmuac)
hist(dat$mfmuac)
# check those MUAC greater than 50 cm
dat[dat$mfmuac > 40 & !is.na(dat$mfmuac), 
    c("serno", "mfwt", "mfstdht", "mfsitht", "mfmuac", "mfwaist", "mftricp1")]
# MUAC for 1161 is too high compared to other measures. This was double checked with Harshad. Delete
dat[dat$serno == 1161, "mfmuac"] <- NA

# waist
summary(dat$mfwaist)
hist(dat$mfwaist)
# check waist that are less than 40
dat[dat$mfwaist < 40 & !is.na(dat$mfwaist), 
    c("serno", "mfwt", "mfstdht", "mfsitht", "mfmuac", "mfwaist", "mftricp1")]
# the waist of those 2 are wrong. Set them to missing
dat[dat$mfwaist < 40 & !is.na(dat$mfwaist), c("mfwaist")] <- NA

# hip
summary(dat$mfhip)
hist(dat$mfhip)

# triceps, biceps, subscapular and suprialiac skinfolds
summary(dat$mftricp1)
summary(dat$mftricp2)
summary(dat$mftricp3)
# check similar measure
diff12 <- with(dat, mftricp1 - mftricp2)
summary(diff12)
rm(diff12)
diff13 <- with(dat, mftricp1 - mftricp3)
summary(diff13)
rm(diff13)
diff23 <- with(dat, mftricp2 - mftricp3)
summary(diff23)
rm(diff23)
# average triceps measures
dat$avgfTric <- rowMeans(dat[, c("mftricp1", "mftricp2", "mftricp3")])
summary(dat$avgfTric)
hist(dat$avgfTric)

summary(dat$mfbicep1)
summary(dat$mfbicep2)
summary(dat$mfbicep3)
# check women has similar measure
diff12 <- with(dat, mfbicep1 - mfbicep2)
summary(diff12)
rm(diff12)
diff13 <- with(dat, mfbicep1 - mfbicep3)
summary(diff13)
rm(diff13)
diff23 <- with(dat, mfbicep2 - mfbicep3)
summary(diff23)
rm(diff23)
# average biceps measures
dat$avgfBic <- rowMeans(dat[, c("mfbicep1", "mfbicep2", "mfbicep3")])
summary(dat$avgfBic)
hist(dat$avgfBic)
dat[dat$avgfBic > 25 & !is.na(dat$avgfBic), 
    c("serno", "mfwt", "mfstdht", "mfsitht", "mfmuac", "mfwaist", "mfhip", "mftricp1",
      "avgfBic")]

summary(dat$mfsubscp1)
summary(dat$mfsubscp2)
summary(dat$mfsubscp3)
# check women has similar measure
diff12 <- with(dat, mfsubscp1 - mfsubscp2)
summary(diff12)
rm(diff12)
diff13 <- with(dat, mfsubscp1 - mfsubscp3)
summary(diff13)
rm(diff13)
diff23 <- with(dat, mfsubscp2 - mfsubscp3)
summary(diff23)
rm(diff23)
# average subscapular measures
dat$avgfSubsc <- rowMeans(dat[, c("mfsubscp1", "mfsubscp2", "mfsubscp3")], na.rm = T)
summary(dat$avgfSubsc)
hist(dat$avgfSubsc)

summary(dat$mfsupri1)
summary(dat$mfsupri2)
summary(dat$mfsupri3)
# check women has similar measure
diff12 <- with(dat, mfsupri1 - mfsupri2)
summary(diff12)
rm(diff12)
diff13 <- with(dat, mfsupri1 - mfsupri3)
summary(diff13)
dat[diff13 > 2 & !is.na(dat$mfsupri1), 
    c("serno", "mfwt", "mfstdht", "mfsitht",
      "mfsupri1", "mfsupri2", "mfsupri3")]
dat[dat$serno == 18, "mfsupri3"] <- NA
rm(diff13)
diff23 <- with(dat, mfsupri2 - mfsupri3)
summary(diff23)
rm(diff23)
# average suprailiac measures
dat$avgfSup <- rowMeans(dat[, c("mfsupri1", "mfsupri2", "mfsupri3")], na.rm = T)
summary(dat$avgfSup)
hist(dat$avgfSup)


# blood pressure

# systolic
summary(dat$mfsys1)
hist(dat$mfsys1)
dat[dat$mfsys1 < 30 & !is.na(dat$mfsys1), 
    c("serno", "mcalcagef", "mfwt", "mfstdht", "mfsitht", "mfmuac",
      "mfsys1", "mfsys2", "mfsys3")]
dat[dat$serno == 508, "mfsys1"] <- NA
summary(dat$mfsys2)
summary(dat$mfsys3)
# average blood pressure values
dat$avgfsys <- rowMeans(dat[, c("mfsys1", "mfsys2", "mfsys3")], na.rm = T)
summary(dat$avgfsys)
hist(dat$avgfsys)

# diastolic
summary(dat$mfdia1)
summary(dat$mfdia2)
summary(dat$mfdia3)
# average blood pressure values
dat$avgfdia <- rowMeans(dat[, c("mfdia1", "mfdia2", "mfdia3")])
summary(dat$avgfdia)
hist(dat$avgfdia)
rm(temp)

#######################################################################################################################
##### Other variables #################################################################################################
#######################################################################################################################

# sli
hist(dat$sliscore)
summary(dat$sliscore)

# add other variables using emphasis01.dta
emph01 <- read_dta("./Confounders/emphasis01.dta")
# take only the variable we need
emph01 <- emph01 %>% select(serno, mrelig, mcschool, mcmedstdy, hiscore)
# merge the data
dat <- inner_join(dat, emph01)
rm(emph01)

# religion
# create 3 categoris: Hindu, Muslim, Other
dat$mrelig[dat$mrelig == 5] <- 3
dat$mrelig[dat$mrelig == 7] <- 3
addmargins(table(dat$mrelig))
# check those with missing religion
dat[is.na(dat$mrelig),c("serno", "mdeath")]
# take them from the registration data
dat[dat$serno %in% c(1164, 1744, 3813, 4752, 4885, 5104), "mrelig"] <- 1
dat[dat$serno == 2014, "mrelig"] <- 3
# those 4 who are death set religion to missing
dat$mrelig <- ifelse(dat$mdeath == "Yes", NA, dat$mrelig)
dat$mrelig <- factor(dat$mrelig, labels = c("Hindu", "Muslim", "Other"))
addmargins(table(dat$mrelig))

# child goes to school?
table(dat$mcschool) 
addmargins(table(dat$mcschool))
# check missing values with Harshad. Only 9 children did not go to school (6 of them are in the EMPHASIS study)
dat$mcschool <- ifelse(dat$serno %in% c(1185, 1356, 1547, 350, 1761, 1762, 5034, 5119, 4885), 0, 1)
dat$mcschool <- factor(dat$mcschool, labels = c("No", "Yes"))
addmargins(table(dat$mcschool))

# child medium of study
dat$mcmedstdy <- factor(dat$mcmedstdy, labels = c("English", "Marathi", "Hindi", "Urdu", "Other"))
addmargins(table(dat$mcmedstdy))

# home inventory score
summary(dat$hiscore)

# Add GDM using the 1999 cutoff (use the GDM computed by Ella in the paper)

gdmdat <- read_dta("T:/Chiara/MMNP_SARASKids/MMNP Glucose and Diabetes (Ella)/MMNP_glucoseanalysis.dta")
# take only the variables of interest
gdmdat <- gdmdat %>% dplyr::select(serno, gdmwho_2006)
names(gdmdat) <- c("serno", "mGDM")
dat <- inner_join(dat, gdmdat)
addmargins(table(dat$mGDM))
rm(gdmdat)
dat$mGDM <- factor(dat$mGDM, labels = c("No GDM", "GDM"))

#######################################################################################################################
##### Save Dataset ####################################################################################################
#######################################################################################################################

# take only the variables needed
dat <- dat %>% dplyr::select(serno, allocation, mage, mrelig, sliscore, mmedulev, mmocc, mmweight, mmstdhgt, 
                             mmsithgt, mmmuac, mmhdcir, mmwaist, mleg, avgmTric, avgmBic, avgmSubsc, avgmSup, 
                             avgmsys, avgmdia, mbmi, mGDM, mdeath,
                             fage, mfedulev, mfocc, mfwt, mfstdht, mfsitht, mfmuac, mfhdcir, mfwaist, fleg, avgfTric,
                             avgfBic, avgfSubsc, avgfSup, avgfsys, avgfdia, fbmi, fdeath, 
                             mcschool, mcmedstdy, hiscore)

names(dat) <- c("serno", "allocation", "mage","mrel", "msliscore", "mEduLevel", "mocc", "mwt", "mht", "msitht", "mMUAC", 
                "mHC","mWaist", "mleg", "mavgTric", "mavgBic", "mavgSubsc", "mavgSup", "mavgSBP", "mavgDBP", "mbmi", 
                "mGDM", "mdeath",
                "fage", "fEduLevel", "focc", "fwt", "fht", "fsitht", "fMUAC", "fHC",
                "fWaist", "fleg", "favgTric", "favgBic", "favgSubsc", "favgSup", "favgSBP", "favgDBP", "fbmi", "fdeath",
                "cschool", "cschoolLang", "chiscore")

# transform NaN in NA
dat$mavgSubsc[is.nan(dat$mavgSubsc)] <- NA
dat$mavgSup[is.nan(dat$mavgSup)] <- NA
dat$mavgSBP[is.nan(dat$mavgSBP)] <- NA
dat$mavgDBP[is.nan(dat$mavgDBP)] <- NA

dat$favgSubsc[is.nan(dat$favgSubsc)] <- NA
dat$favgSup[is.nan(dat$favgSup)] <- NA
dat$favgSBP[is.nan(dat$favgSBP)] <- NA
dat$favgDBP[is.nan(dat$favgDBP)] <- NA

# add attributes label to each variable
attributes(dat$serno)$label <- "saras - Subject ID"
attributes(dat$mage)$label <- "saras - maternal age (years)"
attributes(dat$fage)$label <- "saras - paternal age (years)"
attributes(dat$allocation)$label <- "saras - allocation group"
attributes(dat$mrel)$label <- "saras - religion"
attributes(dat$msliscore)$label <- "saras - SLI score"
attributes(dat$mEduLevel)$label <- "saras - maternal education"
attributes(dat$mocc)$label <- "saras - maternal occupation"
attributes(dat$fEduLevel)$label <- "saras - paternal education"
attributes(dat$focc)$label <- "saras - paternal occupation"
attributes(dat$fdeath)$label <- "saras - paternal death?"
attributes(dat$mdeath)$label <- "saras - maternal death?"

attributes(dat$mwt)$label <- "saras - maternal weight (kg)"
attributes(dat$mht)$label <- "saras - maternal height (cm)"
attributes(dat$msitht)$label <- "saras - maternal sitting height (cm)"
attributes(dat$fwt)$label <- "saras - paternal weight (kg)"
attributes(dat$fht)$label <- "saras - paternal height (cm)"
attributes(dat$fsitht)$label <- "saras - paternal sitting height (cm)"
attributes(dat$mMUAC)$label <- "saras - maternal MUAC (cm)"
attributes(dat$mHC)$label <- "saras - maternal HC (cm)"
attributes(dat$mWaist)$label <- "saras - maternal Waist (cm)"
attributes(dat$fMUAC)$label <- "saras - paternal MUAC (cm)"
attributes(dat$fHC)$label <- "saras - paternal HC (cm)"
attributes(dat$fWaist)$label <- "saras - paternal Waist (cm)"

attributes(dat$mleg)$label <- "saras - maternal leg length (cm)"
attributes(dat$mbmi)$label <- "saras - maternal BMI (kg/m2)"
attributes(dat$fleg)$label <- "saras - paternal leg length (cm)"
attributes(dat$fbmi)$label <- "saras - paternal BMI (kg/m2)"
attributes(dat$mGDM)$label <- "saras - maternal GDM"

attributes(dat$mavgTric)$label <- "saras - maternal triceps (mm)"
attributes(dat$mavgBic)$label <- "saras - maternal biceps (mm)"
attributes(dat$mavgSubsc)$label <- "saras - maternal subscapular skinfold (mm)"
attributes(dat$mavgSup)$label <- "saras - maternal suprailiac (mm)"
attributes(dat$favgTric)$label <- "saras - paternal triceps (mm)"
attributes(dat$favgBic)$label <- "saras - paternal biceps (mm)"
attributes(dat$favgSubsc)$label <- "saras - paternal subscapular skinfold (mm)"
attributes(dat$favgSup)$label <- "saras - paternal suprailiac (mm)"

attributes(dat$mavgSBP)$label <- "saras - maternal SBP"
attributes(dat$mavgDBP)$label <- "saras - maternal DBP"
attributes(dat$favgSBP)$label <- "saras - paternal SBP"
attributes(dat$favgDBP)$label <- "saras - paternal DBP"

attributes(dat$cschool)$label <- "saras - child attend school?"
attributes(dat$cschoolLang)$label <- "saras - child language of education"
attributes(dat$chiscore)$label <- "saras - child home inventory score"

# save as csv, stata and spss
write.table(dat, "./Cleaned Dataset/SARASKidsConfoundersv1_20180314.csv", row.names = FALSE)
write_dta(dat, "./Cleaned Dataset/SARASKidsConfoundersv1_20180314.dta")
write_sav(dat, "./Cleaned Dataset/SARASKidsConfoundersv1_20180314.sav")
