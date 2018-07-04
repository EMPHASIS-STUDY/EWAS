###################################################################################################
#### OUTCOME EWAS ANALYSIS: Cardio-metabolic Data #################################################
###################################################################################################

## 2018/04/11
# load packages
library(haven)
library(dplyr)
library(ggplot2)
library(limma)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(skimr)
library(data.table)

###################################################################################################
#### QQ-PLOT: FUNCTION #####################à######################################################
###################################################################################################

gg_qqplot <- function(p) {
  n  <- length(p)
  df <- data.frame(observed = -log10(sort(p)), expected = -log10(ppoints(n)))
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  ggplot(df) + theme_bw() +
    geom_point(aes(expected, observed), size = 2) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    xlab(log10Pe) + ylab(log10Po)
}

###################################################################################################
#### THE GAMBIA ###################################################################################
###################################################################################################

# load methylation data
load("./Methylation Data/EMPH_EWAS_redux.RData")

names(DMPs_pcs) # differentially methylated position found using PC
DMPs_pcs$coefficients # estimated coefficient for each PC

# remove unused data
rm(DMPs_pcs)
rm(DMRs_pcs)

# add relevant outcomes
datSBP <- read_dta("PMMST_BloodPressureTemp_20180522.dta")
datSBP <- datSBP %>% select(cISubjectID, cavgsbp, cavgdbp) 
# change variable name in datSBP to match that in pdata
names(datSBP)[names(datSBP) == "cISubjectID"] <- "Subject_ID"
pdata <- inner_join(pdata, datSBP)
rm(datSBP)

# get annotation 
anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)


##################### UNADJUSTED ANALYSIS: SBP ########################################

# remove the 2 subject_ID that do not have SBP
# from the design matrix
todropID <- pdata[is.na(pdata$cavgsbp), c("Subject_ID")]
pdata_redux <- pdata[!(pdata$Subject_ID %in% todropID),]
rm(todropID)
# from the M-value matrix
# use the Slide and the Array
todropSlide <- pdata[is.na(pdata$cavgsbp), c("Slide", "Array")]
# merge slide and array together
todropSlide <- paste(todropSlide$Slide, todropSlide$Array, sep = "_")
norm_mval_fil_redux <- norm_mval_fil[, !(colnames(norm_mval_fil) %in% todropSlide)]

# get annotation 
anno_sub <- anno[match(rownames(norm_mval_fil_redux),anno$Name),
                 c(1,2,4,12:17,22,23,26,35)]

# look at SBP
with(pdata, hist(cavgsbp))

#### All ####
design_all <- model.matrix(~ Slide + Array + Bcell + CD4T + CD8T + Eos + Mono + NK + 
                             Neu + MooreSoC +
                             Sex + Age + cavgsbp, pdata_redux)
DMPs_all <- eBayes(lmFit(norm_mval_fil_redux, design_all))
res_DMPs_all <- topTable(DMPs_all, coef = "cavgsbp", number = Inf, 
                         genelist = anno_sub, sort.by="B")
head(res_DMPs_all, 5)
# save all
write.csv(res_DMPs_all[1:100,],file="SBP_EWAS_DMPs_all_Top100.csv")

with(res_DMPs_all, hist(P.Value))

# QQplot
lambda <- signif(median(qchisq(1-res_DMPs_all$P.Value, 1))/qchisq(0.5, 1), 5)
lambda <- paste0(paste0(expression(lambda)," = "),lambda)
gg_qqplot(res_DMPs_all$P.Value) + 
  ggtitle("M ~ Slide + Array + BCCs + Sex + MooreSoC + Age + SBP") + 
  annotate("text", x = 1, y = 6, label = lambda) 
ggsave("SBP_EWAS_DMPs_all_QQ.png")


# remove dataset not used anymore
rm(design_all, DMPs_all, res_DMPs_all)


#### All but no BCC ####
design_all_no_BCC <-  model.matrix(~ Slide + Array + Sex + 
                                     MooreSoC +
                                     Age + cavgsbp, pdata_redux)
DMPs_all_no_BCC <- eBayes(lmFit(norm_mval_fil_redux,design_all_no_BCC))
res_DMPs_all_no_BCC <- topTable(DMPs_all_no_BCC, coef = "cavgsbp", number = Inf, 
                                genelist = anno_sub, sort.by="B")
head(res_DMPs_all_no_BCC, 5)

# save all no bcc
write.csv(res_DMPs_all_no_BCC[1:100,],file="SBP_EWAS_DMPs_all_no_BCC_Top100.csv")

with(res_DMPs_all_no_BCC, hist(P.Value))

# QQplot
lambda <- signif(median(qchisq(1-res_DMPs_all_no_BCC$P.Value, 1))/qchisq(0.5, 1), 5)
lambda <- paste0(paste0(expression(lambda)," = "),lambda)
gg_qqplot(res_DMPs_all_no_BCC$P.Value) + 
  ggtitle("M ~ Slide + Array + Sex + MooreSoC + Age + SBP") + 
  annotate("text", x = 1, y = 6, label = lambda) 
ggsave("SBP_EWAS_DMPs_all_no_BCC_QQ.png")


# remove dataset not used anymore
rm(design_all_no_BCC, DMPs_all_no_BCC, res_DMPs_all_no_BCC)


#### PCs ####

# need to first drop the missing from design_pcs
todropID <- pdata[is.na(pdata$cavgsbp), c("Subject_ID")]
pdata_redux <- pdata[!(pdata$Subject_ID %in% todropID),]
rm(todropID)
# and from the pcs dataset
todropSlide <- pdata[is.na(pdata$cavgsbp), c("Slide", "Array")]
# merge slide and array together
todropSlide <- paste(todropSlide$Slide, todropSlide$Array, sep = "_")
pcs_redux <- pcs[!(rownames(pcs) %in% todropSlide), ]
dim(pcs_redux)

# design matrix and analysis
design_pcs <- model.matrix(~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + 
                             PC10 +
                             PC11 + PC12 + PC13 + PC15 + Age + MooreSoC + cavgsbp , 
                           cbind(pdata_redux,pcs_redux))
DMPs_pcs <- eBayes(lmFit(norm_mval_fil_redux,design_pcs))
res_DMPs_pcs <- topTable(DMPs_pcs, coef = "cavgsbp", number = Inf, 
                         genelist = anno_sub, sort.by="B")
head(res_DMPs_pcs)

# save all PCs
write.csv(res_DMPs_pcs[1:100,],file="SBP_EWAS_DMPs_PCs_Top100.csv")

with(res_DMPs_pcs, hist(P.Value))

# QQplot
lambda <- signif(median(qchisq(1-res_DMPs_pcs$P.Value, 1))/qchisq(0.5, 1), 5)
lambda <- paste0(paste0(expression(lambda)," = "),lambda)
gg_qqplot(res_DMPs_pcs$P.Value) + 
  ggtitle("M ~ PCs + Age + MooreSoC + SBP") + 
  annotate("text", x = 1, y = 6, label = lambda) 
ggsave("SBP_EWAS_DMPs_pcs_QQ.png")


# volcano plot
# Volcano plots of the effect sizes plotted against ???log10(p value) for epigenome-wide 
# discovery meta-analysis associations with systolic and diastolic BP. 
# Effect size units are percent change in DNA methylation per 1-unit change in blood pressure.
ggplot(res_DMPs_pcs, aes(x = logFC, y = -log10(res_DMPs_pcs$P.Value))) +
  geom_point(size=2.5) + theme_bw()
ggsave("PCs_SBP_EWAS_Volcano.png")


##################### ADJUSTED ANALYSIS: SBP ###########################################

# first look at all possible combination of missing data 
head(pdata_redux)
# find missing data
apply(pdata_redux, 2, function(x){sum(is.na(x))})
# no missing heoght, weight, and BMI

# add data on fat mass
datFatMass <- read_dta("PMMSTGrowthBodyCompv1_20180131.dta") %>% dplyr::filter(cvisit == 2) %>%
  dplyr::select(serno, cfatmass, flagcfatmass) 
names(datFatMass)[names(datFatMass) == "serno"] <-  "Subject_ID"  
pdata <- left_join(pdata, datFatMass)
rm(datFatMass)
# 5 have missing Fat Mass. 2 have outlier fat mass
with(pdata, hist(cfatmass))
# check the ones with missing fat mass or outlying fat mass
pdata %>% dplyr::filter(is.na(cfatmass) | flagcfatmass == 2)
# the missing one have blood pressure measured (and so are the outliers)

skim(pdata)
# 2 missing SBP, 5 missing variables fat mass

# SBP by Sex
pdata %>% group_by(Sex) %>% skim(cavgsbp)
# both missing are in female
t.test(pdata$cavgsbp ~ pdata$Sex, var.equal = F)
# SBP by Age
ggplot(pdata, aes(x = Age, y = cavgsbp)) + geom_point()
cor.test(pdata$cavgsbp, pdata$Age, method = "spearman", exact = F)
# SBP by heigth
ggplot(pdata, aes(x = height, y = cavgsbp)) + geom_point()
cor.test(pdata$cavgsbp, pdata$height, method = "pearson", exact = F)
# SBP by MasterGroup
pdata %>% group_by(MasterGroupNo) %>% skim(cavgsbp)
t.test(pdata$cavgsbp ~ pdata$MasterGroupNo, var.equal = F)
# SBP by Season of Conception
pdata %>% group_by(MooreSoC) %>% skim(cavgsbp)
t.test(pdata$cavgsbp ~ pdata$MooreSoC, var.equal = F)
# SBP by fat mass
ggplot(pdata, aes(x = cfatmass, y = cavgsbp)) + geom_point()
cor.test(pdata$cfatmass, pdata$cavgsbp, method = "pearson", exact = F)

# covariates: use the PCs and look and adjustments
# design matrix and analysis

# adjustment 1: Interaction between Intervention and SoC
design_pcs <- model.matrix(~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + 
                             PC9 + PC10 +
                             PC11 + PC12 + PC13 + PC15 + Age + cavgsbp + 
                             MooreSoC*MasterGroupNo, 
                           cbind(pdata_redux,pcs_redux))
DMPs_pcs <- eBayes(lmFit(norm_mval_fil_redux,design_pcs))
res_DMPs_pcs <- topTable(DMPs_pcs, coef = "cavgsbp", number = Inf, 
                         genelist = anno_sub, sort.by="B")
head(res_DMPs_pcs)
with(res_DMPs_pcs, hist(P.Value))
# save
write.csv(res_DMPs_pcs[1:100,],file="SBP_EWAS_DMPs_PCs_Top100_Adj1_Interaction.csv")

# QQplot
lambda <- signif(median(qchisq(1-res_DMPs_pcs$P.Value, 1))/qchisq(0.5, 1), 5)
lambda <- paste0(paste0(expression(lambda)," = "),lambda)
gg_qqplot(res_DMPs_pcs$P.Value) + 
  ggtitle("M ~ PCs + Age + MooreSoC*Group + SBP") + 
  annotate("text", x = 1, y = 6, label = lambda) 
ggsave("SBP_EWAS_DMPs_pcs_QQ_Adj1_Interaction.png")


# adjustment 2: add height and fat mass

# first we need to exclude the 5 children with no fat mass recorded.
# from the design matrix
todropID <- pdata_redux[is.na(pdata_redux$cfatmass), c("Subject_ID")]
pdata_redux <- pdata_redux[!(pdata_redux$Subject_ID %in% todropID),]
rm(todropID)
# from the M-value matrix
# use the Slide and the Array
todropSlide <- pdata[is.na(pdata$cfatmass), c("Slide", "Array")]
# merge slide and array together
todropSlide <- paste(todropSlide$Slide, todropSlide$Array, sep = "_")
norm_mval_fil_redux <- 
  norm_mval_fil_redux[, !(colnames(norm_mval_fil_redux) %in% todropSlide)]
# reduce pcs
pcs_redux <- pcs_redux[!(rownames(pcs_redux) %in% todropSlide), ]
dim(pcs_redux)

# get annotation 
anno_sub <- anno[match(rownames(norm_mval_fil_redux),anno$Name),
                 c(1,2,4,12:17,22,23,26,35)]

design_pcs <- model.matrix(~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + 
                             PC9 + PC10 +
                             PC11 + PC12 + PC13 + PC15 + Age + cavgsbp +
                             MooreSoC*MasterGroupNo + cfatmass + height, 
                           cbind(pdata_redux,pcs_redux))

DMPs_pcs <- eBayes(lmFit(norm_mval_fil_redux,design_pcs))
res_DMPs_pcs <- topTable(DMPs_pcs, coef = "cavgsbp", number = Inf, 
                         genelist = anno_sub, sort.by="B")
head(res_DMPs_pcs)
with(res_DMPs_pcs, hist(P.Value))
# save
write.csv(res_DMPs_pcs[1:100,],file="SBP_EWAS_DMPs_PCs_Top100_Adj2_FatMass.csv")

# QQplot
lambda <- signif(median(qchisq(1-res_DMPs_pcs$P.Value, 1))/qchisq(0.5, 1), 5)
lambda <- paste0(paste0(expression(lambda)," = "),lambda)
gg_qqplot(res_DMPs_pcs$P.Value) + 
  ggtitle("M ~ PCs + Age + MooreSoC*Group + Fat Mass + Height + SBP") + 
  annotate("text", x = 1, y = 6, label = lambda) 
ggsave("SBP_EWAS_DMPs_PCs_Top100_Adj2_FatMass.png")


#######################################################################################
##################### UNADJUSTED ANALYSIS: DBP ########################################
#######################################################################################

# No child has missing DBP

# get annotation 
anno_sub <- anno[match(rownames(norm_mval_fil),anno$Name),
                 c(1,2,4,12:17,22,23,26,35)]

# look at DBP
with(pdata, hist(cavgdbp))

#### All ####
design_all <- model.matrix(~ Slide + Array + Bcell + CD4T + CD8T + Eos + Mono + NK + 
                             Neu + MooreSoC +
                             Sex + Age + cavgdbp, pdata)
DMPs_all <- eBayes(lmFit(norm_mval_fil, design_all))
res_DMPs_all <- topTable(DMPs_all, coef = "cavgdbp", genelist = anno_sub,
                         number = Inf, sort.by="B")
head(res_DMPs_all, 5)
# save all
write.csv(res_DMPs_all[1:100,],file="DBP_EWAS_DMPs_all_Top100.csv")

with(res_DMPs_all, hist(P.Value))

# QQplot
lambda <- signif(median(qchisq(1-res_DMPs_all$P.Value, 1))/qchisq(0.5, 1), 5)
lambda <- paste0(paste0(expression(lambda)," = "),lambda)
gg_qqplot(res_DMPs_all$P.Value) + 
  ggtitle("M ~ Slide + Array + BCCs + MooreSoc + Sex + Age + DBP") + 
  annotate("text", x = 1, y = 6, label = lambda) 
ggsave("DBP_EWAS_DMPs_all_QQ.png")

# remove dataset not used anymore
rm(design_all, DMPs_all, res_DMPs_all)


#### All but no BCC ####
design_all_no_BCC <-  model.matrix(~ Slide + Array + Sex + Age +
                                     MooreSoC + cavgdbp, pdata)

DMPs_all_no_BCC <- eBayes(lmFit(norm_mval_fil,design_all_no_BCC))
res_DMPs_all_no_BCC <- topTable(DMPs_all_no_BCC, coef = "cavgdbp", genelist = anno_sub,
                                number = Inf, sort.by="B")
head(res_DMPs_all_no_BCC, 5)

# save all no bcc
write.csv(res_DMPs_all_no_BCC[1:100,],file="DBP_EWAS_DMPs_all_no_BCC_Top100.csv")

with(res_DMPs_all_no_BCC, hist(P.Value))

# QQplot
lambda <- signif(median(qchisq(1-res_DMPs_all_no_BCC$P.Value, 1))/qchisq(0.5, 1), 5)
lambda <- paste0(paste0(expression(lambda)," = "),lambda)
gg_qqplot(res_DMPs_all_no_BCC$P.Value) + 
  ggtitle("M ~ Slide + Array + Sex + Age + DBP") + 
  annotate("text", x = 1, y = 6, label = lambda) 
ggsave("DBP_EWAS_DMPs_all_no_BCC_QQ.png")

# remove dataset not used anymore
rm(design_all_no_BCC, DMPs_all_no_BCC, res_DMPs_all_no_BCC)


#### PCs ####

# design matrix and analysis
design_pcs <- model.matrix(~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + 
                             PC10 +
                             PC11 + PC12 + PC13 + PC15 + Age + cavgdbp + MooreSoC, 
                           cbind(pdata,pcs))
DMPs_pcs <- eBayes(lmFit(norm_mval_fil,design_pcs))
res_DMPs_pcs <- topTable(DMPs_pcs, coef = "cavgdbp", number = Inf, sort.by="B",
                         genelist = anno_sub)
head(res_DMPs_pcs)

# save all PCs
write.csv(res_DMPs_pcs[1:100,],file="DBP_EWAS_DMPs_PCs_Top100.csv")

with(res_DMPs_pcs, hist(P.Value))

# QQplot
lambda <- signif(median(qchisq(1-res_DMPs_pcs$P.Value, 1))/qchisq(0.5, 1), 5)
lambda <- paste0(paste0(expression(lambda)," = "),lambda)
gg_qqplot(res_DMPs_pcs$P.Value) + 
  ggtitle("M ~ PCs + Age + MooreSoC + DBP") + 
  annotate("text", x = 1, y = 6, label = lambda) 
ggsave("DBP_EWAS_DMPs_PCs_QQ.png")


# volcano plot
# Volcano plots of the effect sizes plotted against ???log10(p value) for epigenome-wide 
# discovery meta-analysis associations with systolic and diastolic BP. 
# Effect size units are percent change in DNA methylation per 1-unit change in blood pressure.
ggplot(res_DMPs_pcs, aes(x = logFC, y = -log10(res_DMPs_pcs$P.Value))) +
  geom_point(size=2.5) + theme_bw()
ggsave("PCs_DBP_EWAS_Volcano.png")



##################### ADJUSTED ANALYSIS: DBP ###########################################

# DBP by Sex
pdata %>% group_by(Sex) %>% skim(cavgdbp)
t.test(pdata$cavgdbp ~ pdata$Sex, var.equal = F)
# DBP by Age
ggplot(pdata, aes(x = Age, y = cavgdbp)) + geom_point()
cor.test(pdata$cavgdbp, pdata$Age, method = "spearman", exact = F)
# DBP by heigth
ggplot(pdata, aes(x = height, y = cavgdbp)) + geom_point()
cor.test(pdata$cavgdbp, pdata$height, method = "pearson", exact = F)
# DBP by MasterGroup
pdata %>% group_by(MasterGroupNo) %>% skim(cavgdbp)
t.test(pdata$cavgdbp ~ pdata$MasterGroupNo, var.equal = F)
# DBP is lower in those in group 2.
summary(lm(cavgdbp ~ MasterGroupNo + Sex + height, data = pdata))
summary(lm(cavgdbp ~ MasterGroupNo*MooreSoC + Sex + height, data = pdata))
# DBP by Season of Conception
pdata %>% group_by(MooreSoC) %>% skim(cavgdbp)
t.test(pdata$cavgdbp ~ pdata$MooreSoC, var.equal = F)
# DBP by fat mass
ggplot(pdata, aes(x = cfatmass, y = cavgdbp)) + geom_point()
cor.test(pdata$cfatmass, pdata$cavgdbp, method = "pearson", exact = F)

# adjustment 1: Interaction between Intervention and SoC
design_pcs <- model.matrix(~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + 
                             PC9 + PC10 +
                             PC11 + PC12 + PC13 + PC15 + Age + cavgdbp + 
                             MooreSoC*MasterGroupNo, 
                           cbind(pdata,pcs))
DMPs_pcs <- eBayes(lmFit(norm_mval_fil,design_pcs))
res_DMPs_pcs <- topTable(DMPs_pcs, coef = "cavgdbp", number = Inf, 
                         genelist = anno_sub, sort.by="B")
head(res_DMPs_pcs)
with(res_DMPs_pcs, hist(P.Value))
# save
write.csv(res_DMPs_pcs[1:100,],file="DBP_EWAS_DMPs_PCs_Top100_Adj1_Interaction.csv")

# QQplot
lambda <- signif(median(qchisq(1-res_DMPs_pcs$P.Value, 1))/qchisq(0.5, 1), 5)
lambda <- paste0(paste0(expression(lambda)," = "),lambda)
gg_qqplot(res_DMPs_pcs$P.Value) + 
  ggtitle("M ~ PCs + Age + MooreSoC*Group + DBP") + 
  annotate("text", x = 1, y = 6, label = lambda) 
ggsave("DBP_EWAS_DMPs_pcs_QQ_Adj1_Interaction.png")

# adjustment 2: add fat mass and height
# first we'll need to remove the 5 with missing fat mass

todropID <- pdata[is.na(pdata$cfatmass), c("Subject_ID")]
pdata_redux <- pdata[!(pdata$Subject_ID %in% todropID),]
rm(todropID)
# from the M-value matrix
# use the Slide and the Array
todropSlide <- pdata[is.na(pdata$cfatmass), c("Slide", "Array")]
# merge slide and array together
todropSlide <- paste(todropSlide$Slide, todropSlide$Array, sep = "_")
norm_mval_fil_redux <- norm_mval_fil[, !(colnames(norm_mval_fil) %in% todropSlide)]
# reduce pcs
pcs_redux <- pcs[!(rownames(pcs) %in% todropSlide), ]
dim(pcs_redux)

# get annotation 
anno_sub <- anno[match(rownames(norm_mval_fil_redux),anno$Name),
                 c(1,2,4,12:17,22,23,26,35)]

design_pcs <- model.matrix(~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + 
                             PC9 + PC10 +
                             PC11 + PC12 + PC13 + PC15 + Age + cavgdbp +
                             MooreSoC*MasterGroupNo + cfatmass + height, 
                           cbind(pdata_redux,pcs_redux))

DMPs_pcs <- eBayes(lmFit(norm_mval_fil_redux,design_pcs))
res_DMPs_pcs <- topTable(DMPs_pcs, coef = "cavgdbp", number = Inf, 
                         genelist = anno_sub, sort.by="B")
head(res_DMPs_pcs)
with(res_DMPs_pcs, hist(P.Value))
# save
write.csv(res_DMPs_pcs[1:100,],file="DBP_EWAS_DMPs_PCs_Top100_Adj2_FatMass.csv")

# QQplot
lambda <- signif(median(qchisq(1-res_DMPs_pcs$P.Value, 1))/qchisq(0.5, 1), 5)
lambda <- paste0(paste0(expression(lambda)," = "),lambda)
gg_qqplot(res_DMPs_pcs$P.Value) + 
  ggtitle("M ~ PCs + Age + MooreSoC*Group + Fat Mass + Height + DBP") + 
  annotate("text", x = 1, y = 6, label = lambda) 
ggsave("DBP_EWAS_DMPs_PCs_Top100_Adj2_FatMass.png")


###################################################################################################
#### INDIA #######################################################################################
###################################################################################################

# To run just once (no need after M-values are created)
# load methylation data
# load("./Methylation Data/norm.beta.pc20.Robj")

# the methylation data are b-value. Transform in M values
# norm.beta[norm.beta == 0] <- 0.00001
# for(i in 1:ncol(norm.beta)){
# norm.beta[, i] <- logit2(norm.beta[, i])
# }
# save m-values
# write.csv(norm.beta, "./Methylation Data/norm_mval", row.names = F)
# rm(norm.beta, i)

# load the M-values
norm_mval <- fread("./Methylation Data/norm_mval")

# load blood pressure data
datBP <- read_dta("SARASKidsCardioMetv4_20180409.dta")
datBP <- datBP %>% select(serno, cage, cavgSBP, cavgDBP, flagcDBP)
head(datBP)
datBP$serno <- paste("SARAS", datBP$serno, sep = "_")

# load IVS and techical variables
isv <- read.csv("SARAS_Int_EWAS_ISVs.csv")
head(isv)
# rename X as serno so it matches datBP
names(isv)[names(isv) == "X"] <- "serno"
# merge isv with datBP
isvBP <- inner_join(isv, datBP)

# remove dataset not needed
rm(datBP, isv)

# 4 misssing SBP and DBP
sum(is.na(isvBP$cavgSBP))
hist(isvBP$cavgSBP)
sum(is.na(isvBP$cavgSBP))
hist(isvBP$cavgDBP)


### ISVs ####

# need to first drop the missing from design_pcs
todropID <- isvBP[is.na(isvBP$cavgSBP), c("serno")]
isvBP_redux <- isvBP[!(isvBP$serno %in% todropID),]
# and from the mvalues dataset
norm_mval <- data.frame(norm_mval)
norm_mval_redux <- norm_mval[, !(colnames(norm_mval) %in% todropID), ]
rm(todropID)

# design matrix and analysis

design_isv <- model.matrix(~ ISV1 + ISV2 + ISV4 + ISV5 + ISV6 + ISV7 + ISV8 + ISV9 + 
                             ISV10 + ISV11 + ISV12 + ISV13 + ISV14  + cage + cavgSBP, 
                           isvBP_redux)
DMPs_pcs <- eBayes(lmFit(norm_mval_redux,design_isv))
res_DMPs_pcs <- topTable(DMPs_pcs, coef = "cavgsbp", number = Inf, 
                         genelist = anno_sub, sort.by="B")
head(res_DMPs_pcs)
