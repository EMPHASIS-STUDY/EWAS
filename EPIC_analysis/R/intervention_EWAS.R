#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
#R script for running EWAS using limma with batch adjustment using ISVA/SV/PCs 
#
#inputs: normalized beta values, phenotypic data
#
# author. Ayden Saffari <ayden.saffari@lshtm.ac.uk>
# affiliations. MRC ING, LSHTM
#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/

library("meffil")
library("minfi")
library("missMethyl")
library("sva")
library("limma")
library("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
library("gridExtra")
library("ggplot2")

source("gamplotlib/ggQQplot.R")
source("gamplotlib/ggVplot.R")
source("gamplotlib/ggPCAplot.R")
source("ilEPICfilter.R")
source("DMR.plot.R")
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#initialization
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
load("../R_objects/norm_beta.RData")
load("../R_objects/pdata.RData")

#ensure samples in beta matrix match those in pdata
norm_beta <- norm_beta[,colnames(norm_beta) %in% rownames(pdata)]
norm_beta <- norm_beta[,match(colnames(norm_beta),rownames(pdata))]
all(colnames(norm_beta) == rownames(pdata))

#convert to m vals
#deal with zero values by adding offset (if needed)
norm_beta_noinf <- norm_beta
norm_beta_noinf[norm_beta_noinf == 0] <- 0.00001
norm_mval <- logit2(norm_beta_noinf)
rm(norm_beta_noinf)

#filter out cross-hybridising probes and chrX and chrY probes
#NB. deal with SNP probes at later stage
norm_beta_fil <- ilEPICfilter(norm_beta, cross=T, sex=T, snpfil=F)
#norm_beta_fil[norm_beta_fil == 0] <- 0.00001
norm_mval_fil <- logit2(norm_beta_fil)

#get annotation
anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
anno_sub <- anno[match(rownames(norm_mval_fil),anno$Name),
c(1,2,4,12:17,22,23,26,35)]

#subset of 200k most variable probes to use for PCA/SVA/ISVA 
#(idea from meffil)
variable <- pdata[,colnames(pdata) %in% c("MasterGroupNo")]
covariates <- pdata[,colnames(pdata) %in% c("Age","MooreSoC")]
var_idx <- order(rowVars(norm_mval_fil, na.rm=T), decreasing=T)[1:200000]
norm_mval_200kvar <- meffil:::impute.matrix(norm_mval_fil[var_idx,,drop=F])

random_seed <- 20170817
set.seed(random_seed)

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#Batch adjustment methods
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

#####################
#Principal components
#####################
pca <- prcomp(t(norm_mval_200kvar),center=T, scale.=F)
pcs <- pca$x
pcpvar <- data.frame(PC = as.factor(seq(1:length(pca$sdev))),
                          per_var = ((pca$sdev^2)/sum(pca$sdev^2)) * 100)
#get variance explained
summary(pca)

#plot scree
ggplot(pcpvar[1:30,], aes(x=PC,y=per_var)) +  
  geom_bar(stat = "identity") + 
  ylab("percentage variance explained") +
  theme_bw() + 
  theme(panel.grid.major = element_line(colour = "#F6F7F7"),
        panel.grid.minor = element_line(colour = "#F6F7F7"))

#correlations/associations with batch variables
pcs_pdata_assoc <- batchPCAcorr(pcs[,1:15],pdata[,-c(8)],15)
write.csv(pcs_pdata_assoc,"../results/EMPH_EPIC_PCA_assoc_200k_fil.csv")

#####################
#SVA
#####################
#cant have sex + mastergroup in same model!
mod <- model.matrix(~.,cbind(covariates,variable))
mod0 <- model.matrix(~.,covariates)
sva <- sva(norm_mval_200kvar, mod=mod, mod0=mod0)
svs <- as.data.frame(sva$sv)
colnames(svs) <- sapply(colnames(svs),function(x){paste0("S",x)})
rownames(svs) <- colnames(norm_mval_200kvar)

svs_pdata_assoc <- batchPCAcorr(svs,pdata[,-c(8)],ncol(svs))
write.csv(svs_pdata_assoc,"../results/EMPH_EPIC_SVA_assoc_200k_fil.csv")


#####################
# ISVA
#####################
#ISVA unsupervised
isva_un <- meffil:::isva(norm_mval_200kvar, as.numeric(variable), verbose=T)
isva_un <- as.data.frame(isva_un$isv)
colnames(isva_un) <- sapply(colnames(isva_un),function(x){paste0("IS",x)})
rownames(isva_un) <- colnames(norm_mval_200kvar)

isvs_un_pdata_assoc <- batchPCAcorr(isva_un,pdata[,-c(8)],ncol(isva_un))
write.csv(isvs_un_pdata_assoc,"../results/EMPH_EPIC_ISV_un_assoc_200k_fil.csv")


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#EWAS 
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

####
#All
####
design_all <- model.matrix(~ Slide + Array + Bcell + CD4T + CD8T + Eos + Mono + NK + Neu + 
                           Sex + Age + MooreSoC + MasterGroupNo, pdata)

DMPs_all <- lmFit(norm_mval_fil,design_all)
DMPs_all <- eBayes(DMPs_all)
res_DMPs_all <- topTable(DMPs_all, coef = "MasterGroupNo2", number = Inf,
                         genelist=anno_sub, sort.by="B")

#write results
write.csv(as.data.frame(res_DMPs_all)[1:100,],file="../results/EPIC_EWAS_top_100_DMPs_all.csv")

#QQplot
lambda <- signif(median(qchisq(1-res_DMPs_all$P.Value,1))/qchisq(0.5,1),5)
lambda <- paste0(paste0(expression(lambda)," = "),lambda)
ggQQplot(res_DMPs_all$P.Value,ylim = c(0,8)) + theme_bw() + 
ggtitle("M ~ Slide + Array + BCCs + Sex + Age + MooreSoC + MasterGroupNo") + 
theme(plot.title = element_text(colour="black", size=10)) +
annotate("text", x = 4, y = 8, label = lambda)
ggsave("../results/EPIC_EWAS_DMPs_all_QQ.png")

############
#All no BCC
############
design_all_no_BCC <-  model.matrix(~ Slide + Array + Sex + Age + MooreSoC + MasterGroupNo, pdata)

DMPs_all_no_BCC <- lmFit(norm_mval_fil,design_all_no_BCC)
DMPs_all_no_BCC  <- eBayes(DMPs_all_no_BCC )
res_DMPs_all_no_BCC <- topTable(DMPs_all_no_BCC, coef = "MasterGroupNo2", number = Inf,
                         genelist=anno_sub, sort.by="B")

write.csv(as.data.frame(res_DMPs_all_no_BCC)[1:100,],file="../results/EPIC_EWAS_top_100_DMPs_all_no_BCC.csv")

#QQplot
lambda <- signif(median(qchisq(1-res_DMPs_all_no_BCC$P.Value,1))/qchisq(0.5,1),5)
lambda <- paste0(paste0(expression(lambda)," = "),lambda)
ggQQplot(res_DMPs_all_no_BCC$P.Value,ylim = c(0,8)) + theme_bw() + 
ggtitle("M ~ Slide + Array + Sex + Age + MooreSoC + MasterGroupNo") + 
theme(plot.title = element_text(colour="black", size=10)) +
annotate("text", x = 4, y = 8, label = lambda)
ggsave("../results/EPIC_EWAS_DMPs_all_NO_BCC_QQ.png")

#####
#SVs
#####
design_svs <- model.matrix(~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + SV11 + 
                           SV12 + SV13 + SV15 + SV17 + SV19 + SV20 + Age + MooreSoC + MasterGroupNo, cbind(pdata,svs)

DMPs_svs <- lmFit(norm_mval_fil,design_svs)
DMPs_svs <- eBayes(DMPs_svs)
res_DMPs_svs <- topTable(DMPs_svs, coef = "MasterGroupNo2", number = Inf,
                         genelist=anno_sub, sort.by="B")

#write results
write.csv(as.data.frame(res_DMPs_svs)[1:100,],file="../results/EPIC_EWAS_top_100_DMPs_svs.csv")

#QQplot
lambda <- signif(median(qchisq(1-res_DMPs_svs$P.Value,1))/qchisq(0.5,1),5)
lambda <- paste0(paste0(expression(lambda)," = "),lambda)
ggQQplot(res_DMPs_svs$P.Value,ylim = c(0,8)) + theme_bw() + 
ggtitle("M ~ SVs + Age + MooreSoC + MasterGroupNo") + 
theme(plot.title = element_text(colour="black", size=10)) +
annotate("text", x = 4, y = 8, label = lambda)
ggsave("../results/EPIC_EWAS_DMPs_svs_QQ.png")

###############
#ISVs (unsup) 
##############
#+ bio vars (excluding sex) 
design_isvs_un <- model.matrix(~ ISV1 + ISV3 + ISV4 + ISV5 + ISV6 + Sex + 
                               Age + MooreSoC + MasterGroupNo,cbind(pdata,isva_un))

DMPs_isvs <- lmFit(norm_mval_fil,design_isvs_un)
DMPs_isvs <- eBayes(DMPs_isvs)
res_DMPs_isvs <- topTable(DMPs_isvs, coef = "MasterGroupNo2", number = Inf,
                         genelist=anno_sub, sort.by="B")

#write results
write.csv(as.data.frame(res_DMPs_isvs)[1:100,],file="../results/EPIC_EWAS_top_100_DMPs_isvs.csv")

#QQplot
lambda <- signif(median(qchisq(1-res_DMPs_isvs$P.Value,1))/qchisq(0.5,1),5)
lambda <- paste0(paste0(expression(lambda)," = "),lambda)
ggQQplot(res_DMPs_isvs$P.Value,ylim = c(0,8)) + theme_bw() + 
ggtitle("M ~ ISVs + Sex + Age + MooreSoC + MasterGroupNo") + 
theme(plot.title = element_text(colour="black", size=10)) +
annotate("text", x = 4, y = 8, label = lambda)
ggsave("../results/EPIC_EWAS_DMPs_isvs_QQ.png")

#########
#PCs
#########
design_pcs <- model.matrix(~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
                          PC11 + PC12 + PC13 + PC15 + Age + MooreSoC + MasterGroupNo, cbind(pdata,pcs))

DMPs_pcs <- lmFit(norm_mval_fil,design_pcs)
DMPs_pcs <- eBayes(DMPs_pcs)
res_DMPs_pcs <- topTable(DMPs_pcs, coef = "MasterGroupNo2", number = Inf,
                         genelist=anno_sub, sort.by="B")

#write results
write.csv(as.data.frame(res_DMPs_pcs)[1:100,],file="../results/EPIC_EWAS_top_100_DMPs_pcs.csv")

#QQplot
lambda <- signif(median(qchisq(1-res_DMPs_pcs$P.Value,1))/qchisq(0.5,1),5)
lambda <- paste0(paste0(expression(lambda)," = "),lambda)
ggQQplot(res_DMPs_pcs$P.Value,ylim = c(0,8)) + theme_bw() + 
ggtitle("M ~ PCs + Age + MooreSoC + MasterGroupNo") + 
theme(plot.title = element_text(colour="black", size=10)) +
annotate("text", x = 4, y = 8, label = lambda)
ggsave("../results/EPIC_EWAS_DMPs_pcs_QQ.png")
 
                           
#add column for delta Beta
add_delta_B <- function(res,beta,inter){
  beta_res <- beta[match(res$Name,rownames(beta)),]
  delta_B <- apply(beta_res,1, function(x){
                      mean(x[which(inter == 2)]) -  mean(x[which(inter == 1)])})
  delta_B <- cbind(delta_B)
  res <- merge(res,delta_B,by.x='row.names',by.y='row.names') 
  res <- res[with(res,order(-B)),]
  return(res) 
}

res_DMPs_pcs_deltab <- add_delta_B(as.data.frame(res_DMPs_pcs)[1:100,],norm_beta_fil,pdata$MasterGroupNo)  
write.csv(res_DMPs_pcs_deltab,file="../results/EPIC_EWAS_top_100_DMPs_pcs_deltab.csv") 
                           
#######################################
#PCs with MasterGroup x SoC interaction 
#######################################
design_pcs_inter <- model.matrix(~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
                          PC11 + PC12 + PC13 + PC15 + Age + MooreSoC * MasterGroupNo, cbind(pdata,pcs))

DMPs_pcs_inter <- lmFit(norm_mval_fil,design_pcs_inter)
DMPs_pcs_inter <- eBayes(DMPs_pcs_inter)
res_DMPs_pcs_inter <- topTable(DMPs_pcs_inter, coef = "MooreSoCrainy:MasterGroupNo2", number = Inf,
                         genelist=anno_sub, sort.by="B")

inter_main_group <- topTable(DMPs_pcs_inter, coef = "MasterGroupNo2", number = Inf, genelist=anno_sub, sort.by="B")
inter_main_season <- topTable(DMPs_pcs_inter, coef = "MooreSoCrainy", number = Inf, genelist=anno_sub, sort.by="B")
inter_group_season <- topTable(DMPs_pcs_inter, 
                               coef = "MooreSoCrainy:MasterGroupNo2", number = Inf, genelist=anno_sub, sort.by="B")
                           
#QQplot
lambda <- signif(median(qchisq(1-inter_main_group$P.Value,1))/qchisq(0.5,1),5)
lambda <- paste0(paste0(expression(lambda)," = "),lambda)
ggQQplot(inter_main_group$P.Value,ylim = c(0,8)) + theme_bw() + 
ggtitle("M ~ PCs + Age + MooreSoC * MasterGroupNo") + 
theme(plot.title = element_text(colour="black", size=10)) +
annotate("text", x = 4, y = 8, label = lambda)
ggsave("../results/EPIC_EWAS_DMPs_pcs_inter_QQ.png")

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#Pathway analysis
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
###
#GO
###
gst_all <- gometh(sig.cpg=res_DMPs_all$Name[1:1000],collection="GO",array.type="EPIC")
gst_all_no_bcc <- gometh(sig.cpg=res_DMPs_all_no_BCC$Name[1:1000],collection="GO",array.type="EPIC")
gst_svs <- gometh(sig.cpg=res_DMPs_svs$Name[1:1000],collection="GO",array.type="EPIC")
gst_isvs <- gometh(sig.cpg=res_DMPs_isvs$Name[1:1000],collection="GO",array.type="EPIC")
gst_pcs <- gometh(sig.cpg=res_DMPs_pcs$Name[1:1000],collection="GO",array.type="EPIC")

#write results                           
write.csv(topGO(gst_all),file="../results/EPIC_EWAS_GO_pathways_all.csv")
write.csv(topGO(gst_all_no_bcc),file="../results/EPIC_EWAS_GO_pathways_all_no_BCC.csv")
write.csv(topGO(gst_svs),file="../results/EPIC_EWAS_GO_pathways_SVs.csv")
write.csv(topGO(gst_isvs),file="../results/EPIC_EWAS_GO_pathways_ISVs.csv")
write.csv(topGO(gst_pcs),file="../results/EPIC_EWAS_GO_pathways_PCs.csv")

#####
#KEGG
#####
gst_all_KEGG <- gometh(sig.cpg=res_DMPs_all$Name[1:1000],collection="KEGG",array.type="EPIC")
gst_all_no_bcc_KEGG <- gometh(sig.cpg=res_DMPs_all_no_BCC$Name[1:1000],collection="KEGG",array.type="EPIC")
gst_svs_KEGG <- gometh(sig.cpg=res_DMPs_svs$Name[1:1000],collection="KEGG",array.type="EPIC")
gst_isvs_KEGG <- gometh(sig.cpg=res_DMPs_isvs$Name[1:1000],collection="KEGG",array.type="EPIC")
gst_pcs_KEGG <- gometh(sig.cpg=res_DMPs_pcs$Name[1:1000],collection="KEGG",array.type="EPIC")

#write results
write.csv(topKEGG(gst_all_KEGG),file="../results/EPIC_EWAS_KEGG_pathways_all.csv")
write.csv(topKEGG(gst_all_no_bcc_KEGG),file="../results/EPIC_EWAS_KEGG_pathways_all_no_BCC.csv")
write.csv(topKEGG(gst_svs_KEGG),file="../results/EPIC_EWAS_KEGG_pathways_SVs.csv")
write.csv(topKEGG(gst_isvs_KEGG),file="../results/EPIC_EWAS_KEGG_pathways_ISVs.csv")
write.csv(topKEGG(gst_pcs_KEGG),file="../results/EPIC_EWAS_KEGG_pathways_PCs.csv")

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#Regional analysis
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

DMRs_pcs <- fastEPICdmrcate(norm_mval_fil,design_pcs,coef="MasterGroupNo2",pcutoff=0.1, mc.cores=4)
DMRs_svs <- fastEPICdmrcate(norm_mval_fil,design_svs,coef="MasterGroupNo2",pcutoff=0.1, mc.cores=4)
DMRs_isvs <- fastEPICdmrcate(norm_mval_fil,design_isvs,coef="MasterGroupNo2",pcutoff=0.1, mc.cores=4)

#write results
write.csv(DMRs_pcs$results,file="../results/EPIC_EWAS_DMRs_PCs.csv")
write.csv(DMRs_svs$results,file="../results/EPIC_EWAS_DMRs_SVs.csv")
write.csv(DMRs_isvs$results,file="../results/EPIC_EWAS_DMRs_ISVs.csv")

#annotated regions
DMRs_pcs_annot <- extractRanges(DMRs_pcs, genome = "hg19")
write.csv(DMRs_pcs_annot,file="../results/EPIC_EWAS_DMRs_annot_PCs.csv")
        
#visualizations
#DMR plots                           
type <- pdata$MasterGroupNo
groups <- c("1"="#F97E86", "2"="#0F7597")
cols <- groups[as.character(type)]
DMR.plot(ranges=DMRs_pcs_annot, dmr=1, CpGs=norm_beta_fil, what="Beta", arraytype = "EPIC",
         phen.col=cols, genome="hg19")
DMR.plot(ranges=DMRs_pcs_annot, dmr=277, CpGs=norm_beta_fil, what="Beta", arraytype = "EPIC",
         phen.col=cols, genome="hg19")
DMR.plot(ranges=DMRs_pcs_annot, dmr=2, CpGs=norm_beta_fil, what="Beta", arraytype = "EPIC",
         phen.col=cols, genome="hg19")
                           
#vplot
ggVplot(DMRs_pcs$results[,c(3,4,5)],fdrcut=0.1, lfccut=0.05, xlim=c(-0.1,0.1),n=DMRs_pcs$results[,2]) + theme_gamplotlib()
 + scale_fill_gamplotlib() + scale_color_gamplotlib() + xlab("fold change")    
ggsave("../results/DMR_pcs_fc_vplot.pdf")          
                           
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#Additional Plots
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

#plot PCs  
pcs_df <- cbind(sample_sheet, pcs[,1:4])

pal <- c("#F87E86","#0E7597","#2E7B39","#FDA80A",
"#B31918","#9CD893","#E6BD9B","#B7BAAA","#683A2E","#211318")

#Sex
pdf(file="../results/ENID_PCS_Sex.pdf", onefile=FALSE)
p1 <- ggPCAplot(pcs_df, pcs=c(1,2), batch="Sex",pal=pal) + theme_bw() + coord_fixed(ratio=1) + theme(legend.position="none")
p2 <- ggPCAplot(pcs_df, pcs=c(1,3), batch="Sex",pal=pal) + theme_bw() + coord_fixed(ratio=1) + theme(legend.position="none")
p3 <- ggPCAplot(pcs_df, pcs=c(2,3), batch="Sex",pal=pal) + theme_bw() + coord_fixed(ratio=1) + theme(legend.position="none")
p3l <- ggPCAplot(pcs_df, pcs=c(2,3), batch="Sex",pal=pal) + theme_bw() + coord_fixed(ratio=1) 

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

legr1 <- g_legend(p3l)
grid.arrange(p1,p2,p3,legr1,nrow=1,heights=c(2),widths=c(2,2,2,1.1))
dev.off()

#sentrix_row
pdf(file="../results/ENID_PCS_sentrix_row.pdf", onefile=FALSE)
p1 <- ggPCAplot(pcs_df, pcs=c(1,2), batch="sentrix_row",pal=pal) + theme_bw() + coord_fixed(ratio=1) + theme(legend.position="none")
p2 <- ggPCAplot(pcs_df, pcs=c(1,3), batch="sentrix_row",pal=pal) + theme_bw() + coord_fixed(ratio=1) + theme(legend.position="none")
p3 <- ggPCAplot(pcs_df, pcs=c(2,3), batch="sentrix_row",pal=pal) + theme_bw() + coord_fixed(ratio=1) + theme(legend.position="none")
p3l <- ggPCAplot(pcs_df, pcs=c(2,3), batch="sentrix_row",pal=pal) + theme_bw() + coord_fixed(ratio=1) 

legr1 <-g_legend(p3l)
grid.arrange(p1,p2,p3,legr1,nrow=1,heights=c(2),widths=c(2,2,2,1.1))
dev.off()

#Sample plate
pdf(file="../results/ENID_PCS_Sample_Plate.pdf", onefile=FALSE)
p1 <- ggPCAplot(pcs_df, pcs=c(1,2), batch="Sample_Plate",pal=pal) + theme_bw() + coord_fixed(ratio=1) + theme(legend.position="none")
p2 <- ggPCAplot(pcs_df, pcs=c(1,3), batch="Sample_Plate",pal=pal) + theme_bw() + coord_fixed(ratio=1) + theme(legend.position="none")
p3 <- ggPCAplot(pcs_df, pcs=c(2,3), batch="Sample_Plate",pal=pal) + theme_bw() + coord_fixed(ratio=1) + theme(legend.position="none")
p3l <- ggPCAplot(pcs_df, pcs=c(2,3), batch="Sample_Plate",pal=pal) + theme_bw() + coord_fixed(ratio=1) 
legr1 <-g_legend(p3l)
grid.arrange(p1,p2,p3,legr1,nrow=1,heights=c(2),widths=c(2,2,2,1.1))
dev.off()

#TODO - remaining vars
