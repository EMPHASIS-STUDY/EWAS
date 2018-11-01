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
res_DMPs_pcs <- readRDS("../R_objects/res_DMPs_pcs.rds")
DMRs_CpGs <- readRDS("../R_objects/EMPH_GMB_DMRs_CpGs.rds")
norm_beta_fil <- readRDS("../R_objects/norm_beta_fil.rds")
#GMB_CpGs <- norm_beta_fil[which(rownames(norm_beta_fil) %in%
#                          unique(c(res_DMPs_pcs$Name[
#                          res_DMPs_pcs$adj.P.Val < 0.1],DMRs_CpGs))),]

GMB_CpGs <- norm_beta_fil[which(rownames(norm_beta_fil) %in%
                          DMRs_CpGs),]

pcs <- readRDS("../R_objects/pcs.rds")
pdata <- cbind(pdata,pcs[,1:15])
GSA_sample_sheet <- read.csv("/data/GSA/emphasis/EMPHASIS_GMB_GSA_Samplesheet.csv")

GMB_SNPs <- read.table("../data/GSA_GMB_PLINKfiltered_a_hwe_geno_maf_recodeA_t.traw",
                       sep="\t", head=T)
GMB_SNPs <- GMB_SNPs[,-c(1,3,4,5,6)]

#############################
#produce genotype summary stat table
#############################
summary_table <- apply(GMB_SNPs[,-1],1,function(x){summary(as.factor(x))})
summary_table <- ldply(summary_table,function(s){t(data.frame(unlist(s)))})
summary_table_fil <- summary_table[!(summary_table$`NA's` > 30),]
summary_table_fil$SNP <- rownames(GMB_SNPs)
summary_table_fil <- summary_table_fil[,match(c("SNP","0", "1","2","NA's"),
                                        colnames(summary_table_fil))]
##########################
#reshape GSA data for GEM
###########################
rownames(GMB_SNPs) <- GMB_SNPs[,1]
GMB_SNPs <- GMB_SNPs[,-1]

pdata <- readRDS("../R_objects/pdata.rds")
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
#REMOVE -don't need to do this, 0,1,2 coding is equivalent 
#GMB_SNPs_f <- t(as.data.frame(apply(GMB_SNPs,1,factor)))
#GMB_SNPs_f <- revalue(GMB_SNPs_f, c("0"="1", "1"="2", "2"="3"))

###########################
#reshape EPIC data for GEM
############################
#replace arrays with sample IDs for EPIC array
GMB_CpGs <- GMB_CpGs[,match(rownames(pdata),colnames(GMB_CpGs))]
all(colnames(GMB_CpGs) == rownames(pdata))
colnames(GMB_CpGs) <- pdata$Subject_ID
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

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# run GEM mQTL analysis
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#create env and cov objects
env <- pdata[,colnames(pdata) %in% c("Subject_ID","MasterGroupNo")]
rownames(env) <- env$Subject_ID
env <- t(env[,colnames(env) == "MasterGroupNo",drop=F])
rownames(env)

cov <- pdata[,colnames(pdata) %in% c("Subject_ID","PC1","PC2","PC3","PC4",
"PC5","PC6","PC7","PC8","PC9","PC10", "PC11","PC12","PC13","PC15","Age",
"MooreSoC","MasterGroupNo"),drop=F]
cov$MooreSoC <- recode(cov$MooreSoC,dry="1",rainy="2")
cov$MooreSoC <- relevel(cov$MooreSoC,"1")

cov <- dcast(melt(cov, id.var = "Subject_ID"), ... ~ Subject_ID )
rownames(cov) <- cov$variable
cov <- cov[,-1]
cov <- cov[,match(pdata$Subject_ID,colnames(cov))]
dim(cov)

#create combined cov file
cov_env <- rbind(cov[rownames(cov) != "MasterGroupNo",],
                 cov[rownames(cov) == "MasterGroupNo",])

#remove mastergroup from cov
cov <- cov[rownames(cov) != "MasterGroupNo",]

#convert to numeric
cov_num <- sapply(cov[,], as.numeric)
rownames(cov_num) <- rownames(cov)

cov_env_num <- sapply(cov_env[,], as.numeric)
rownames(cov_env_num) <- rownames(cov_env)

env_num <- t(as.data.frame(sapply(env[,,drop=F], as.numeric)))
rownames(env_num) <- rownames(env)
colnames(env_num) <- colnames(env)

dim(cov_num)
rownames(cov_num)
dim(env_num)
rownames(env_num)
dim(cov_env_num)
rownames(cov_env_num)
all.equal(colnames(cov_num),as.character(pdata$Subject_ID))
all.equal(colnames(cov_num),colnames(GMB_SNPs))
all.equal(colnames(cov_num),colnames(GMB_CpGs))
all.equal(colnames(cov_num),colnames(env_num))
all.equal(colnames(env_num),as.character(pdata$Subject_ID))
all.equal(colnames(env_num),colnames(GMB_SNPs))
all.equal(colnames(env_num),colnames(GMB_CpGs))
all.equal(colnames(cov_env_num),colnames(cov_num))

#save as text files
write.table(GMB_SNPs,"../data/GMB_SNPs.txt",sep="\t")
write.table(GMB_CpGs,"../data/GMB_CpGs.txt",sep="\t")
write.table(env_num,"../data/GMB_env.txt",sep="\t")
write.table(cov_num,"../data/GMB_cov.txt",sep="\t")
write.table(cov_env_num,"../data/GMB_gxe.txt",sep="\t")

#run GEM models
GEM_Emodel("../data/GMB_env.txt", "../data/GMB_cov.txt", "../data/GMB_CpGs.txt",
           1,"../results/GEM/Result_Emodel.txt", "../results/GEM/Emodel_QQ.jpg",
           savePlot=T)
GEM_Gmodel("../data/GMB_SNPs.txt","../data/GMB_cov.txt","../data/GMB_CpGs.txt",
          1e-04, "../results/GEM/Result_Gmodel.txt")
GEM_GxEmodel("../data/GMB_SNPs.txt", "../data/GMB_gxe.txt", "../data/GMB_CpGs.txt",
             1, "../results/GEM/Result_GEmodel.txt", topKplot = 1, savePlot=T)

#Run regression with additive genotype and interaction
GxE_reg_top <- merge(t(GMB_SNPs[rownames(GMB_SNPs) %in%
                     c("rs1423249"),]),
                     t(GMB_CpGs[rownames(GMB_CpGs) %in% c("cg06837426","cg20673840","cg20451680",
                        "cg14972155","cg20059697",
                        "cg13106512","cg21180956"),]),by="row.names")
rownames(GxE_reg_top) <- GxE_reg_top$Row.names
GxE_reg_top <- GxE_reg_top[,-1]
GxE_reg_top <- GxE_reg_top[match(pdata$Subject_ID,
                           rownames(GxE_reg_top)),]
GxE_reg_top <- merge(GxE_reg_top,pdata,by.x='row.names',by.y='Subject_ID')

#alternative with genotype as factor
#not run

summary(lm(cg14972155 ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 +
           PC7 + PC8 + PC9 + PC10 +
           PC11 + PC12 + PC13 + PC15 + Age + MooreSoC +
           MasterGroupNo +
           rs1423249,GxE_reg_top))

#meth x inter x geno plot
GxE_reg_top_fil <- GxE_reg_top[,colnames(GxE_reg_top) %in%
       c("Row.names","rs10239100","rs1423249","rs278368","cg20673840",
         "cg06837426","cg14972155","MasterGroupNo")]
GxE_reg_top_fil <- na.omit(GxE_reg_top_fil)
colnames(GxE_reg_top_fil)[which(colnames(GxE_reg_top_fil) == "MasterGroupNo")] <- "intervention"
GxE_reg_top_fil$intervention <- revalue(GxE_reg_top_fil$intervention, c("1"="intervention","2"="control"))
GxE_reg_top_fil$intervention <- relevel(GxE_reg_top_fil$intervention,"control")
GxE_reg_top_fil$rs1423249 <- as.factor(GxE_reg_top_fil$rs1423249)
GxE_reg_top_fil$rs1423249 <- revalue(GxE_reg_top_fil$rs1423249, c("0"="GG","1"="GA","2"="AA"))

#cg06837426 ~ rs1423249
ggplot(GxE_reg_top_fil, aes(rs1423249,cg06837426),color=rs1423249) + 
       geom_point(aes(color = rs1423249)) + 
       scale_color_manual(values=c("#C04B8E","#C04B8E","#C04B8E")) + 
       stat_summary(aes(y = cg06837426,group=rs1423249),fun.y=mean,colour="#252997",
       geom="line",group=1) + 
       theme_gamplotlib() + theme(strip.background = element_blank(),
       panel.grid.major.x = element_blank(),legend.position="none",
       aspect.ratio=1) +
       scale_x_discrete() +
       ylab("methylation Beta value") + xlab("genotype") +
       ggtitle("cg06837426 ~ rs1423249")
ggsave("../results/GMB_mQTL_cg06837426_rs1423249_G_scatter.pdf",width=(4),
       height=3.5, units="in", dpi=300)

#cg06837426 ~ rs1423249:intervention
ggplot(GxE_reg_top_fil, aes(intervention,cg06837426),color=intervention) +
       geom_point(aes(color = intervention)) + 
       scale_color_manual(values=c("#46617A","#00B8A2")) + 
       stat_summary(aes(y = cg06837426,group=intervention), fun.y=mean,
       colour="#252997", geom="line",group=1) + facet_wrap( ~ rs1423249) +
       theme_gamplotlib() + theme(strip.background = element_blank(),
       panel.grid.major.x = element_blank(),legend.position="none",
       aspect.ratio=1) +
       scale_x_discrete(labels=c("control","inter.")) +
       ylab("methylation Beta value") + xlab("group") +
       ggtitle("cg06837426 ~ rs1423249:intervention")
ggsave("../results/GMB_mQTL_cg06837426_rs1423249_GxE_scatter.pdf",width=(4),
       height=3.5, units="in", dpi=300)

#cg20673840 ~ rs1423249
ggplot(GxE_reg_top_fil, aes(rs1423249,cg20673840),color=rs1423249) +
       geom_point(aes(color = rs1423249)) + scale_color_manual(values=c("#C04B8E",
       "#C04B8E","#C04B8E")) + stat_summary(aes(y = cg20673840,group=rs1423249), 
       fun.y=mean,colour="#252997", geom="line",group=1) +
       theme_gamplotlib() + theme(strip.background = element_blank(),
       panel.grid.major.x = element_blank(),legend.position="none",
       aspect.ratio=1) +
       scale_x_discrete() +
       ylab("methylation Beta value") + xlab("genotype") +
       ggtitle("cg20673840 ~ rs1423249")
ggsave("../results/GMB_mQTL_cg20673840_rs1423249_G_scatter.pdf",width=(4),
height=3.5, units="in", dpi=300)

#cg20673840 ~ rs1423249:intervention
ggplot(GxE_reg_top_fil, aes(intervention,cg20673840),color=intervention) +
       geom_point(aes(color = intervention)) + scale_color_manual(values=c("#46617A",
       "#00B8A2")) + stat_summary(aes(y = cg20673840,group=intervention), fun.y=mean,
       colour="#252997", geom="line",group=1) + facet_wrap( ~ rs1423249) + 
       theme_gamplotlib() + theme(strip.background = element_blank(),
       panel.grid.major.x = element_blank(),legend.position="none",
       aspect.ratio=1) +
       scale_x_discrete(labels=c("control","inter.")) +
       ylab("methylation Beta value") + xlab("group") + 
       ggtitle("cg20673840 ~ rs1423249:intervention")
ggsave("GMB_mQTL_cg20673840_rs1423249_GxE_scatter.pdf",width=(4),
       height=3.5, units="in", dpi=300)

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# run additional analyses
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#G only mQTLs
######

#factor for intervention vs control
GxE_reg_top$intervention <- relevel(GxE_reg_top$MasterGroupNo,"2")

#rs1423249
mQTL_cpgs <- c("cg06837426","cg20673840","cg20451680",
                "cg14972155","cg20059697",
                "cg13106512","cg21180956")

#G
###
print("G")
res_mQTL_cpgs  <- lapply(mQTL_cpgs, function(x) {
                         lm(substitute(cpg ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
                                       PC11 + PC12 + PC13 + PC15 + Age + MooreSoC + rs1423249,
                                       list(cpg = as.name(x))), data = GxE_reg_top)})
names(res_mQTL_cpgs) <- mQTL_cpgs

#coeffs and adj R sqrd
res_mQTL_cpgs_summ <- lapply(res_mQTL_cpgs,summary)
lapply(res_mQTL_cpgs_summ,function(x){x$coefficients[c(18),]})
lapply(res_mQTL_cpgs_summ,function(x){x$adj.r.squared})
#AIC
lapply(res_mQTL_cpgs,AIC)

#E
###
print("E")
res_mQTL_cpgs  <- lapply(mQTL_cpgs, function(x) {
                         lm(substitute(cpg ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
                                       PC11 + PC12 + PC13 + PC15 + Age + MooreSoC + intervention,
                                       list(cpg = as.name(x))), data = GxE_reg_top)})
names(res_mQTL_cpgs) <- mQTL_cpgs

#coeffs and adj R sqrd
res_mQTL_cpgs_summ <- lapply(res_mQTL_cpgs,summary)
lapply(res_mQTL_cpgs_summ,function(x){x$coefficients[c(18),]})
lapply(res_mQTL_cpgs_summ,function(x){x$adj.r.squared})
#AIC
lapply(res_mQTL_cpgs,AIC)

#G+E
####
print("G + E")
res_mQTL_cpgs  <- lapply(mQTL_cpgs, function(x) {
                         lm(substitute(cpg ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
                                       PC11 + PC12 + PC13 + PC15 + Age + MooreSoC + rs1423249 + intervention,
                                       list(cpg = as.name(x))), data = GxE_reg_top)})
names(res_mQTL_cpgs) <- mQTL_cpgs

#coeffs and adj R sqrd
res_mQTL_cpgs_summ <- lapply(res_mQTL_cpgs,summary)
lapply(res_mQTL_cpgs_summ,function(x){x$coefficients[c(18,19),]})
lapply(res_mQTL_cpgs_summ,function(x){x$adj.r.squared})
#AIC
lapply(res_mQTL_cpgs,AIC)



#GxE
####
print("G x E")
res_mQTL_cpgs  <- lapply(mQTL_cpgs, function(x) {
  lm(substitute(cpg ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
                  PC11 + PC12 + PC13 + PC15 + Age + MooreSoC + rs1423249*MasterGroupNo,
                list(cpg = as.name(x))), data = GxE_reg_top)})
names(res_mQTL_cpgs) <- mQTL_cpgs

#coeffs and adj R sqrd
res_mQTL_cpgs_summ <- lapply(res_mQTL_cpgs,summary)
lapply(res_mQTL_cpgs_summ,function(x){x$coefficients[c(18,19,20),]})
lapply(res_mQTL_cpgs_summ,function(x){x$adj.r.squared})
#AIC
lapply(res_mQTL_cpgs,AIC)

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# redo mQTL with imputed data for chr 5 and 8
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# chr 5
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

GMB_SNPs_chr5_3 <- read.table("/data/GSA/emphasis/imputed_geno_chr5_chr8/CHROMOSOME_5/PED_FORMAT/plink/fcgene_plink_chr5_info03_a_hwe_geno_maf_recodeA_t.traw", sep="\t", head=T)

GMB_SNPs_chr5_9 <- read.table("/data/GSA/emphasis/imputed_geno_chr5_chr8/CHROMOSOME_5/PED_FORMAT/plink/fcgene_plink_chr5_info09_a_hwe_geno_maf_recodeA_t.traw", sep="\t", head=T)

GMB_SNPs_chr5 <- GMB_SNPs_chr5_9
GMB_SNPs_chr5 <- GMB_SNPs_chr5[,-c(1,3,4,5,6)]

##########################
#reshape GSA data for GEM
###########################
SNPs_chr5 <- GMB_SNPs_chr5[,1,drop=F]
GMB_SNPs_chr5 <- GMB_SNPs_chr5[,-1]

#edit sample names to match those in sample sheet/beta matrix
colnames(GMB_SNPs_chr5) <- paste0("",sapply(colnames(GMB_SNPs_chr5),
                                           function(x){strsplit(x,"X?[0:9]*_")[[1]][2]}))

#match sample order in sample sheet to GSA data
GSA_sample_sheet_chr5 <- GSA_sample_sheet[match(colnames(GMB_SNPs_chr5),
                                           GSA_sample_sheet$Sample.ID),]
all(GSA_sample_sheet_chr5$Sample.ID == colnames(GMB_SNPs_chr5))

#add probe names back in
#rownames(GMB_SNPs_chr5) <- make.names(t(SNPs_chr5), unique=TRUE)
rownames(GMB_SNPs_chr5) <- t(SNPs_chr5)
###########################
#reshape EPIC data for GEM
############################
GMB_SNPs_chr5 <- GMB_SNPs_chr5[,colnames(GMB_SNPs_chr5) %in% colnames(GMB_CpGs)]
GMB_CpGs_chr5 <- GMB_CpGs[,match(colnames(GMB_SNPs_chr5),colnames(GMB_CpGs))]

all(colnames(GMB_SNPs_chr5) == colnames(GMB_CpGs_chr5))
dim(GMB_CpGs_chr5)
dim(GMB_SNPs_chr5)

#match pdata sample order to GSA and EPIC
pdata_chr5 <- pdata[match(as.factor(colnames(GMB_SNPs_chr5)),pdata$Subject_ID),]
all(colnames(GMB_SNPs_chr5) == pdata_chr5$Subject_ID)
all(colnames(GMB_CpGs_chr5) == pdata_chr5$Subject_ID)
dim(pdata_chr5)

################
#run GEM models
################
#create env and cov objects
env <- pdata_chr5[,colnames(pdata_chr5) %in% c("Subject_ID","MasterGroupNo")]
rownames(env) <- env$Subject_ID
env <- t(env[,colnames(env) == "MasterGroupNo",drop=F])
rownames(env)

cov <- pdata_chr5[,colnames(pdata_chr5) %in% c("Subject_ID","PC1","PC2","PC3","PC4",
                 "PC5","PC6","PC7","PC8","PC9","PC10", "PC11","PC12","PC13","PC15","Age",
                 "MooreSoC","MasterGroupNo"),drop=F]
cov$MooreSoC <- recode(cov$MooreSoC,dry="1",rainy="2")
cov$MooreSoC <- relevel(cov$MooreSoC,"1")

cov <- dcast(melt(cov, id.var = "Subject_ID"), ... ~ Subject_ID )
rownames(cov) <- cov$variable
cov <- cov[,-1]
cov <- cov[,match(pdata_chr5$Subject_ID,colnames(cov))]
dim(cov)

#create combined cov file
cov_env <- rbind(cov[rownames(cov) != "MasterGroupNo",],
                 cov[rownames(cov) == "MasterGroupNo",])

#remove mastergroup from cov
cov <- cov[rownames(cov) != "MasterGroupNo",]

cov_num <- sapply(cov[,], as.numeric)
rownames(cov_num) <- rownames(cov)

cov_env_num <- sapply(cov_env[,], as.numeric)
rownames(cov_env_num) <- rownames(cov_env)

env_num <- t(as.data.frame(sapply(env[,,drop=F], as.numeric)))
rownames(env_num) <- rownames(env)
colnames(env_num) <- colnames(env)

dim(cov_num)
rownames(cov_num)
dim(env_num)
rownames(env_num)
dim(cov_env_num)
rownames(cov_env_num)
all.equal(colnames(cov_num),as.character(pdata_chr5$Subject_ID))
all.equal(colnames(cov_num),colnames(GMB_SNPs_chr5))
all.equal(colnames(cov_num),colnames(GMB_CpGs_chr5))
all.equal(colnames(cov_num),colnames(env_num))
all.equal(colnames(env_num),as.character(pdata_chr5$Subject_ID))
all.equal(colnames(env_num),colnames(GMB_SNPs_chr5))
all.equal(colnames(env_num),colnames(GMB_CpGs_chr5))
all.equal(colnames(cov_env_num),colnames(cov_num))

#save as text files
write.table(GMB_SNPs_chr5,"../data/GMB_SNPs_chr5.txt",sep="\t")
write.table(GMB_CpGs_chr5,"../data/GMB_CpGs_chr5.txt",sep="\t")
write.table(env_num,"../data/GMB_env_chr5.txt",sep="\t")
write.table(cov_num,"../data/GMB_cov_chr5.txt",sep="\t")
write.table(cov_env_num,"../data/GMB_gxe_chr5.txt",sep="\t")

#run GEM models
GEM_Gmodel("../data/GMB_SNPs_chr5.txt",
			"../data/GMB_cov_chr5.txt",
			"../data/GMB_CpGs_chr5.txt",1e-04,
			"../results/GEM/Result_Gmodel_chr5_imputed_9.txt")
GEM_GxEmodel("../data/GMB_SNPs_chr5.txt",
			"../data/GMB_gxe_chr5.txt",
			"../data/GMB_CpGs_chr5.txt", 1,
			 "../results/GEM/Result_GEmodel_chr5_imputed_9.txt",
			  topKplot = 1, savePlot=T)


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# chr 8
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

GMB_SNPs_chr8_3 <- read.table("/data/GSA/emphasis/imputed_geno_chr5_chr8/CHROMOSOME_8/PED_FORMAT/plink/fcgene_plink_chr8_info03_a_hwe_geno_maf_recodeA_t.traw", sep="\t", head=T)

GMB_SNPs_chr8_9 <- read.table("/data/GSA/emphasis/imputed_geno_chr5_chr8/CHROMOSOME_8/PED_FORMAT/plink/fcgene_plink_chr8_info09_a_hwe_geno_maf_recodeA_t.traw", sep="\t", head=T)

GMB_SNPs_chr8 <- GMB_SNPs_chr8_9
GMB_SNPs_chr8 <- GMB_SNPs_chr8[,-c(1,3,4,5,6)]

##########################
#reshape GSA data for GEM
###########################
SNPs_chr8 <- GMB_SNPs_chr8[,1,drop=F]
GMB_SNPs_chr8 <- GMB_SNPs_chr8[,-1]

#edit sample names to match those in sample sheet/beta matrix
colnames(GMB_SNPs_chr8) <- paste0("",sapply(colnames(GMB_SNPs_chr8),
                                           function(x){strsplit(x,"X?[0:9]*_")[[1]][2]}))

#match sample order in sample sheet to GSA data
GSA_sample_sheet_chr8 <- GSA_sample_sheet[match(colnames(GMB_SNPs_chr8),
                                           GSA_sample_sheet$Sample.ID),]
all(GSA_sample_sheet_chr8$Sample.ID == colnames(GMB_SNPs_chr8))

#add probe names back in
rownames(GMB_SNPs_chr8) <- t(SNPs_chr8)

###########################
#reshape EPIC data for GEM
############################
GMB_SNPs_chr8 <- GMB_SNPs_chr8[,colnames(GMB_SNPs_chr8) %in% colnames(GMB_CpGs)]
GMB_CpGs_chr8 <- GMB_CpGs[,match(colnames(GMB_SNPs_chr8),colnames(GMB_CpGs))]

all(colnames(GMB_SNPs_chr8) == colnames(GMB_CpGs_chr8))
dim(GMB_CpGs_chr8)
dim(GMB_SNPs_chr8)

#match pdata sample order to GSA and EPIC
pdata_chr8 <- pdata[match(as.factor(colnames(GMB_SNPs_chr8)),pdata$Subject_ID),]
all(colnames(GMB_SNPs_chr8) == pdata_chr8$Subject_ID)
all(colnames(GMB_CpGs_chr8) == pdata_chr8$Subject_ID)
dim(pdata_chr8)

#create env and cov objects
env <- pdata_chr8[,colnames(pdata_chr8) %in% c("Subject_ID","MasterGroupNo")]
rownames(env) <- env$Subject_ID
env <- t(env[,colnames(env) == "MasterGroupNo",drop=F])
rownames(env)

cov <- pdata_chr8[,colnames(pdata_chr8) %in% c("Subject_ID","PC1","PC2","PC3","PC4",
"PC5","PC6","PC7","PC8","PC9","PC10", "PC11","PC12","PC13","PC15","Age","MooreSoC","MasterGroupNo"),drop=F]
cov$MooreSoC <- recode(cov$MooreSoC,dry="1",rainy="2")
cov$MooreSoC <- relevel(cov$MooreSoC,"1")

cov <- dcast(melt(cov, id.var = "Subject_ID"), ... ~ Subject_ID )
rownames(cov) <- cov$variable
cov <- cov[,-1]
cov <- cov[,match(pdata_chr8$Subject_ID,colnames(cov))]
dim(cov)

#create combined cov file
cov_env <- rbind(cov[rownames(cov) != "MasterGroupNo",],
                 cov[rownames(cov) == "MasterGroupNo",])

#remove mastergroup from cov
cov <- cov[rownames(cov) != "MasterGroupNo",]

cov_num <- sapply(cov[,], as.numeric)
rownames(cov_num) <- rownames(cov)

cov_env_num <- sapply(cov_env[,], as.numeric)
rownames(cov_env_num) <- rownames(cov_env)

env_num <- t(as.data.frame(sapply(env[,,drop=F], as.numeric)))
rownames(env_num) <- rownames(env)
colnames(env_num) <- colnames(env)

dim(cov_num)
rownames(cov_num)
dim(env_num)
rownames(env_num)
dim(cov_env_num)
rownames(cov_env_num)
all.equal(colnames(cov_num),as.character(pdata_chr8$Subject_ID))
all.equal(colnames(cov_num),colnames(GMB_SNPs_chr8))
all.equal(colnames(cov_num),colnames(GMB_CpGs_chr8))
all.equal(colnames(cov_num),colnames(env_num))
all.equal(colnames(env_num),as.character(pdata_chr8$Subject_ID))
all.equal(colnames(env_num),colnames(GMB_SNPs_chr8))
all.equal(colnames(env_num),colnames(GMB_CpGs_chr8))
all.equal(colnames(cov_env_num),colnames(cov_num))

#save as text files
write.table(GMB_SNPs_chr8,"../data/GMB_SNPs_chr8.txt",sep="\t")
write.table(GMB_CpGs_chr8,"../data/GMB_CpGs_chr8.txt",sep="\t")
write.table(env_num,"../data/GMB_env_chr8.txt",sep="\t")
write.table(cov_num,"../data/GMB_cov_chr8.txt",sep="\t")
write.table(cov_env_num,"../data/GMB_gxe_chr8.txt",sep="\t")

#run GEM models
GEM_Gmodel("../data/GMB_SNPs_chr8.txt",
			"../data/GMB_cov_chr8.txt",
			"../data/GMB_CpGs_chr8.txt",1e-04,
			"../results/GEM/Result_Gmodel_chr8_imputed_9.txt")
GEM_GxEmodel("../data/GMB_SNPs_chr8.txt",
			"../data/GMB_gxe_chr8.txt",
			"../data/GMB_CpGs_chr8.txt", 1,
			 "../results/GEM/Result_GEmodel_chr8_imputed_9.txt",
			  topKplot = 1, savePlot=F)
