#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
#R script for running mQTL analysis for EWAS using GEM/matrixEQTL 
#
#inputs: matrix of methylation beta values (EPIC), matrix of SNP genotypes (GSA),
#        phenotype data
#
# authors. Ayden Saffari <ayden.saffari@lshtm.ac.uk> (MRC ING, LSHTM)
#          Ashutosh Singh Tomar (CSIR, CCMB)
#          Prachand Issarapu (CSIR, CCMB)
#     
# NOT FOR DISTRIBUTION/ PUBLIC CONSUMPTION
#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/

library("GEM")
library("plyr")
library("dplyr")
library("reshape2")
library("ggplot2")
library("gamplotlib")

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#initialization
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
##########
#load data
##########
pdata <- readRDS("../R_objects/pdata.rds")
res_DMPs_pcs <- readRDS("../R_objects/res_DMPs_pcs.rds")
DMRs_CpGs <- readRDS("../R_objects/EMPH_GMB_DMRs_CpGs.rds")
norm_beta_fil <- readRDS("../R_objects/norm_beta_fil.rds")
GMB_CpGs <- norm_beta_fil[which(rownames(norm_beta_fil) %in% 
                          unique(c(res_DMPs_pcs$Name[res_DMPs_pcs$adj.P.Val < 0.1],DMRs_CpGs))),]
pcs <- readRDS("../R_objects/pcs.rds")
pdata <- cbind(pdata,pcs[,1:15])
GSA_sample_sheet <- read.csv("/data/GSA/emphasis/EMPHASIS_GMB_GSA_Samplesheet.csv")

GMB_SNPs <- read.table("../data/GSA_GMB_PLINKfiltered_a_hwe_geno_maf_recodeA_t.traw", sep="\t", head=T)
GMB_SNPs <- GMB_SNPs[,-c(1,3,4,5,6)]

#############################
#produce summary stats for genotypes
#############################
summary_table <- apply(GMB_SNPs[,-1],1,function(x){summary(as.factor(x))})
summary_table <- ldply(summary_table,function(s){t(data.frame(unlist(s)))})
summary_table_fil <- summary_table[!(summary_table$`NA's` > 30),]
#NAs
dim(summary_table) - dim(summary_table_fil)

#monomorphic AA/BB
#plink A coding chooses minor/major allele based on freq, so monos always coded as 0
length(which(is.na(summary_table_fil$`0`) & is.na(summary_table_fil$`1`)))

#polymorphic - any minor alleles (1 or 2)
length(which(apply(summary_table_fil[,-1],1,function(x){length(which(is.na(x)))}) <= 1))

#heterzygote only
length(which((is.na(summary_table_fil$`0`)) & (is.na(summary_table_fil$`2`)) & 
               !(is.na(summary_table_fil$`1`))))

##########################
#reshape GSA data for GEM
###########################
SNPs <- GMB_SNPs[,1,drop=F] 
GMB_SNPs <- GMB_SNPs[,-1]

#edit sample names to match those in sample sheet/beta matrix
colnames(GMB_SNPs) <- paste0("2",sapply(colnames(GMB_SNPs),
                                            function(x){strsplit(x,"X?[0:9]*_2")[[1]][2]}))

#match sample order in sample sheet to GSA data
GSA_sample_sheet <- GSA_sample_sheet[match(colnames(GMB_SNPs),
                                           GSA_sample_sheet$Array.info),]
all(GSA_sample_sheet$Array.info == colnames(GMB_SNPs))

#replace arrays with sample IDs for GSA
colnames(GMB_SNPs) <- GSA_sample_sheet$Sample.ID[
                          GSA_sample_sheet$Array.info == colnames(GMB_SNPs)]

#change genotype coding to 1,2,3
#REMOVE -shouldn't need to do this, 0,1,2 coding *should* work
#GMB_SNPs <- GMB_SNPs + 1

#add column for probe names back in
GMB_SNPs <- cbind(SNPs,GMB_SNPs)
rownames(GMB_SNPs) <- GMB_SNPs$SNP
GMB_SNPs <- GMB_SNPs[,-1]

###########################
#reshape EPIC data for GEM
############################
#replace arrays with sample IDs for EPIC array
GMB_CpGs <- GMB_CpGs[,match(rownames(pdata),colnames(GMB_CpGs))]
colnames(GMB_CpGs) <- pdata$Subject_ID[colnames(GMB_CpGs) == rownames(pdata)]
GMB_CpGs <- as.data.frame(GMB_CpGs)

#dont need to do - add ID column and move to 1
#GMB_CpGs$ID <- rownames(GMB_CpGs)
#GMB_CpGs <- GMB_CpGs[,c(ncol(GMB_CpGs),1:(ncol(GMB_CpGs) - 1))]

GMB_SNPs <- GMB_SNPs[,colnames(GMB_SNPs) %in% colnames(GMB_CpGs)]
GMB_CpGs <- GMB_CpGs[,match(colnames(GMB_SNPs),colnames(GMB_CpGs))]

#rownames(GMB_CpGs) <- GMB_CpGs$ID
#GMB_CpGs <- GMB_CpGs[,-1]

all(colnames(GMB_SNPs) == colnames(GMB_CpGs))
dim(GMB_CpGs)
dim(GMB_SNPs)

#match pdata sample order to GSA and EPIC
pdata <- pdata[match(as.factor(colnames(GMB_SNPs)),pdata$Subject_ID),]
all(colnames(GMB_SNPs) == pdata$Subject_ID)
all(colnames(GMB_CpGs) == pdata$Subject_ID)
dim(pdata)

#try with just ESM1
ESM1_DMR <- res_DMPs_pcs[res_DMPs_pcs$Name %in% DMRs_CpGs,13,drop=F]
ESM1_DMR <- ESM1_DMR[grep("ESM1",ESM1_DMR$GencodeCompV12_NAME),,drop=F]
GMB_CpGs <- GMB_CpGs[rownames(GMB_CpGs) %in% rownames(ESM1_DMR),]
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# run GEM mQTL analysis
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#create env and cov objects
env <- t(pdata[,colnames(pdata) == "MasterGroupNo",drop=F])

cov <- pdata[,colnames(pdata) %in% 
               c("Subject_ID","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10",
                 "PC11","PC12","PC13","PC15","Age","MooreSoC","MasterGroupNo"),drop=F]
cov <- dcast(melt(cov, id.var = "Subject_ID"), ... ~ Subject_ID )
cov[cov=="dry"] <- 1
cov[cov=="rainy"] <- 2
cov <- cov[-2,]

rownames(cov) <- cov$variable
cov <- cov[,-1]
dim(cov)

cov_num <- sapply(cov[,], as.numeric)
rownames(cov_num) <- rownames(cov)

env_num <- t(as.data.frame(sapply(env[,,drop=F], as.numeric)))
rownames(env_num) <- rownames(env)
colnames(env_num) <- colnames(env)
#save as text files
write.table(GMB_SNPs,"../data/GMB_SNPs.txt",sep="\t")
write.table(GMB_CpGs,"../data/GMB_CpGs.txt",sep="\t")
write.table(env_num,"../data/GMB_env.txt",sep="\t")
write.table(cov_num,"../data/GMB_cov.txt",sep="\t")
write.table(rbind(cov_num,env_num),"../data/GMB_gxe.txt",sep="\t")

#run GEM models
GEM_Emodel("../data/GMB_env.txt", "../data/GMB_cov.txt", "../data/GMB_CpGs.txt",
           1,"../results/GEM/Result_Emodel.txt", "../results/GEM/Emodel_QQ.jpg", 
           savePlot=T)
GEM_Gmodel("../data/GMB_SNPs.txt","../data/GMB_cov.txt","../data/GMB_CpGs.txt",
          1e-04, "../results/GEM/Result_Gmodel.txt")
GEM_GxEmodel("../data/GMB_SNPs.txt", "../data/GMB_gxe.txt", "../data/GMB_CpGs.txt", 
             1, "../results/GEM/Result_GEmodel.txt", topKplot = 1, savePlot=T)

#Run regression with genotype and interaction
GxE_reg_top <- cbind(t(GMB_SNPs[rownames(GMB_SNPs) %in% c("rs1423249","rs10239100","rs278368"),]),
                     t(GMB_CpGs[rownames(GMB_CpGs) %in% c("cg06837426","cg20673840","cg20451680",
                                                        "cg14972155","cg21180956","cg20059697","cg13106512"),]))

GxE_reg_top <- merge(GxE_reg_top,pdata,by.x='row.names',by.y='Subject_ID')
#GxE_reg_top$rs10239100 <- as.factor(GxE_reg_top$rs10239100)
#GxE_reg_top$rs1423249 <- as.factor(GxE_reg_top$rs1423249)

summary(lm(cg20673840 ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + 
           PC11 + PC12 + PC13 + PC15 + Age + MooreSoC + 
           rs10239100 * MasterGroupNo + 
           rs1423249 * MasterGroupNo,GxE_reg_top))

#frequency plots
histo_geno_inter <- as.data.frame(table(GxE_reg_top$rs10239100:GxE_reg_top$MasterGroupNo))
colnames(histo_geno_inter)[1] <- "Genotype:MasterGroupNo" 
ggplot(histo_geno_inter, aes(x=`Genotype:MasterGroupNo`,y=Freq)) + geom_histogram(stat="identity")
ggsave("GMB_mQTL_rs10239100_geno_inter_histogram.pdf",height=7,width=7)

histo_geno_inter <- as.data.frame(table(GxE_reg_top$rs1423249:GxE_reg_top$MasterGroupNo))
colnames(histo_geno_inter)[1] <- "Genotype:MasterGroupNo" 
ggplot(histo_geno_inter, aes(x=`Genotype:MasterGroupNo`,y=Freq)) + geom_histogram(stat="identity")
ggsave("GMB_mQTL_rs1423249_geno_inter_histogram.pdf",height=7,width=7)

#meth x inter x geno plot
GxE_reg_top_fil <- GxE_reg_top[,colnames(GxE_reg_top) %in% 
       c("Row.names","rs10239100","rs1423249","rs278368","cg20673840","cg14972155","MasterGroupNo")]
GxE_reg_top_fil <- na.omit(GxE_reg_top_fil) 
colnames(GxE_reg_top_fil)[7] <- "intervention"
GxE_reg_top_fil$intervention <- revalue(GxE_reg_top_fil$intervention, c("1"="intervention","2"="control"))
GxE_reg_top_fil$intervention <- relevel(GxE_reg_top_fil$intervention,"control")
GxE_reg_top_fil$rs10239100 <- as.factor(GxE_reg_top_fil$rs10239100)
GxE_reg_top_fil$rs1423249 <- as.factor(GxE_reg_top_fil$rs1423249)
GxE_reg_top_fil$rs278368 <- as.factor(GxE_reg_top_fil$rs278368)
GxE_reg_top_fil$rs10239100 <- revalue(GxE_reg_top_fil$rs10239100, c("0"="CC","1"="CA","2"="AA"))
GxE_reg_top_fil$rs1423249 <- revalue(GxE_reg_top_fil$rs1423249, c("0"="GG","1"="GA","2"="AA"))
GxE_reg_top_fil$rs278368 <- revalue(GxE_reg_top_fil$rs278368, c("0"="GG","1"="GA","2"="AA"))

ggplot(GxE_reg_top_fil, aes(intervention,cg20673840)) + 
       geom_point(color="#264049") + stat_summary(aes(y = cg20673840,group=intervention), fun.y=mean,
       colour="#30F3B0", geom="line",group=1) + facet_wrap( ~ rs10239100) + 
       theme_gamplotlib() + theme(strip.background = element_blank()) + 
       ggtitle("cg20673840 ~ rs10239100:intervention")
ggsave("../results/GMB_mQTL_cg20673840_rs10239100_GxE_scatter.pdf",height=7,width=7)

ggplot(GxE_reg_top_fil, aes(intervention,cg14972155)) + 
       geom_point(color="#264049") + stat_summary(aes(y = cg14972155,group=intervention), fun.y=mean,
       colour="#30F3B0", geom="line",group=1) + facet_wrap( ~ rs278368) + 
       theme_gamplotlib() + theme(strip.background = element_blank()) + 
       ggtitle("cg14972155 ~ rs278368:intervention")
ggsave("../results/GMB_mQTL_cg14972155_rs278368_GxE_scatter.pdf",height=7,width=7)



#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# run additional analyses
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#G only mQTLs
######
#rs1423249

mQTL_cpgs <- c("cg06837426","cg20673840","cg20451680","cg14972155","cg21180956","cg20059697","cg13106512")
res_mQTL_cpgs  <- lapply(mQTL_cpgs, function(x) {
                         lm(substitute(cpg ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + 
                                       PC11 + PC12 + PC13 + PC15 + Age + MooreSoC + rs1423249 + MasterGroupNo,
                                       list(cpg = as.name(x))), data = GxE_reg_top)})
names(res_mQTL_cpgs) <- mQTL_cpgs

#AIC
lapply(res_mQTL_cpgs,AIC)

#coeffs and adj R sqrd
res_mQTL_cpgs <- lapply(res_mQTL_cpgs,summary)
lapply(res_mQTL_cpgs,function(x){x$coefficients[c(18,19),]})
lapply(res_mQTL_cpgs,function(x){x$adj.r.squared})

#GxEs
#####
GxE_1 <- lm(cg20673840 ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + 
           PC11 + PC12 + PC13 + PC15 + Age + MooreSoC + 
           rs10239100 * MasterGroupNo + 
           rs1423249,GxE_reg_top)

GxE_2 <- lm(cg14972155 ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + 
           PC11 + PC12 + PC13 + PC15 + Age + MooreSoC + 
           rs278368 * MasterGroupNo + 
           rs1423249 ,GxE_reg_top)

summary(GxE_1)$adj.r.squared
AIC(GxE_1)

summary(GxE_2)$adj.r.squared
AIC(GxE_2)

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# phenoscanner GWAS lookup
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#rs1423249, rs10239100, rs278368

#summarise by trait, give min p
pheno_GWAS <- read.csv("../results/EWAS/mQTL/PhenoScanner/GMB_EPIC_mQTLs_11855_PhenoScanner_GWAS.csv")
pheno_GWAS <- pheno_GWAS[pheno_GWAS$SNP %in% c("rs1423249", "rs10239100", "rs278368"),]

pheno_GWAS_summ <- inner_join(pheno_GWAS %>% group_by(Trait) %>% tally, 
                         pheno_GWAS %>% group_by(Trait) %>% summarise(minP=min(P)), by="Trait")
pheno_GWAS_summ <- arrange(pheno_GWAS_summ, -n)

write.csv(pheno_GWAS_summ,"../results/EWAS/mQTL/pheno_GWAS_summ.csv")

