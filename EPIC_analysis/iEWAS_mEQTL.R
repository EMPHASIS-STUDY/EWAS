#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
#R script for running mQTL analysis for EWAS using GEM/matrixEQTL 
#
#inputs: matrix of methylation beta values (EPIC), matrix of SNP genotypes (GSA),
#        phenotype data
#
# author. Ayden Saffari <ayden.saffari@lshtm.ac.uk>
# affiliations. MRC ING, LSHTM
#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
library("plyr")
library("GEM")
library("reshape2")

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#initialization
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

##########
#load data
##########
pdata <- readRDS("../R_objects/pdata.rds")
res_DMPs_pcs <- readRDS("../R_objects/res_DMPs_pcs.rds")
norm_beta_fil <- readRDS("../R_objects/norm_beta_fil.rds")
sig_norm_beta_fil <- norm_beta_fil[which(rownames(norm_beta_fil) %in% 
                                  res_DMPs_pcs$Name[res_DMPs_pcs$adj.P.Val < 0.1]),]

pcs <- readRDS("../R_objects/pcs.rds")
pdata <- cbind(pdata,pcs[,1:15])
GSA_sample_sheet <- read.csv("/data/GSA/emphasis/EMPHASIS_GMB_GSA_Samplesheet.csv")

GMB_SNPs <- read.table("../data/GSA_GMB_PLINKfiltered_geno_t.traw", sep="\t", head=T)
GMB_SNPs <- GMB_SNPs[,-c(1,3,4,5,6)]

#############################
#summary stats for genotypes
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
############################
SNPs <- GMB_SNPs[,1,drop=F] 
GMB_SNPs <- GMB_SNPs[,-1]

#edit sample names to match those in sample sheet/beta matrix
colnames(GMB_SNPs)[-1] <- paste0("2",sapply(colnames(GMB_SNPs)[-1],
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
rownames(GMB_SNPs) <- GMB_SNPs$ID
GMB_SNPs <- GMB_SNPs[,-1]

###########################
#reshape EPIC data for GEM
############################
#replace arrays with sample IDs for EPIC array
sig_norm_beta_fil <- sig_norm_beta_fil[,match(rownames(pdata),
                                           colnames(sig_norm_beta_fil))]
colnames(sig_norm_beta_fil) <- pdata$Subject_ID[colnames(
                                              sig_norm_beta_fil) == rownames(pdata)]
sig_norm_beta_fil <- as.data.frame(sig_norm_beta_fil)

#add ID column and move to 1
sig_norm_beta_fil$ID <- rownames(sig_norm_beta_fil)
sig_norm_beta_fil <- sig_norm_beta_fil[,c(ncol(sig_norm_beta_fil),
                                          1:(ncol(sig_norm_beta_fil) - 1))]
GMB_CpGs <- sig_norm_beta_fil

GMB_SNPs <- GMB_SNPs[,colnames(GMB_SNPs) %in% colnames(GMB_CpGs)]


GMB_CpGs <- GMB_CpGs[,match(colnames(GMB_SNPs),colnames(GMB_CpGs))]

rownames(GMB_CpGs) <- GMB_CpGs$ID
GMB_CpGs <- GMB_CpGs[,-1]

all(colnames(GMB_SNPs) == colnames(GMB_CpGs))
dim(GMB_CpGs)
dim(GMB_SNPs)

#match pdata sample order to GSA and EPIC
pdata <- pdata[match(as.factor(colnames(GMB_SNPs)[-1]),pdata$Subject_ID),]
dim(pdata)

#create env and cov object
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

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# run GEM
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
GEM_Emodel("../data/GMB_env.txt", "../data/GMB_cov.txt", "../data/GMB_CpGs.txt",
           1,"../results/GEM/Result_Emodel.txt", "../results/GEM/Emodel_QQ.txt", 
           savePlot=T)
GEM_Gmodel("../data/GMB_SNPs.txt","../data/GMB_cov.txt","../data/GMB_CpGs.txt",
          1e-04, "../results/GEM/Result_Gmodel.txt")
GEM_GxEmodel("../data/GMB_SNPs.txt", "../data/GMB_gxe.txt", "../data/GMB_CpGs.txt", 
             1, "../results/GEM/Result_GEmodel.txt", topKplot = 1, savePlot=T)


#TODOs


#Run through analysis with example dataset
#Limit to cis within 5 Kb
#Rerun with all probes in DMRs
#Check GE number of tests
#Run regression with genotype and interaction
#Try with season of conception as exposure
