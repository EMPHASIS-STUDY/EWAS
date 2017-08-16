#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
#R script for running EWAS using limma with bacth adjustment using ISVA/SV/PCs 
#
#inputs: normalized beta values, phenotypic data
#
# author. Ayden Saffari <ayden.saffari@lshtm.ac.uk>
# affiliations. MRC ING, LSHTM
#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/

library("meffil")
library("minfi")
library("sva")
library("limma")
library("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
library("gridExtra")
library("ggplot2")

source("ggPCAplot.R")
source("ilEPICfilter.R")

#^^^^^^^^^^^^^^^^^^^^^^^^
#initialization
#^^^^^^^^^^^^^^^^^^^^^^^^
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

#^^^^^^^^^^^^^^^^^^^^^^^^
#Principal components
#^^^^^^^^^^^^^^^^^^^^^^^^
#NB - perfom using unfiltered dataset
pca <- prcomp(t(norm_mval),center=T, scale.=F)
pcs <- pca$x
pcpvar <- data.frame(PC = as.factor(seq(1:length(pca$sdev))),
                          per_var = ((pca$sdev^2)/sum(pca$sdev^2)) * 100)
#get variance explained
summary(pca)

#correlations/associations with batch variables
pcs_pdata_assoc <- batchPCAcorr(pcs[,1:6],pdata[,-c(8)],6)
write.csv(pcs_pdata_assoc,"../results/EMPH_EPIC_PCA_assoc.csv")

#^^^^^^^^^^^^^^^^^^^^^^^^
#SVA
#^^^^^^^^^^^^^^^^^^^^^^^^
#NB - perfom using unfiltered dataset


#^^^^^^^^^^^^^^^^^^^^^^^^
# ISVA
#^^^^^^^^^^^^^^^^^^^^^^^^
#NB - perfom using unfiltered dataset

#TODO: replace this with custom script
featureset <- meffil:::guess.featureset(rownames(norm_mval))
features <- meffil.get.features(featureset)
autosomal.sites <- meffil.get.autosomal.sites(featureset)
autosomal.sites <- intersect(autosomal.sites, rownames(norm_mval))

#ISVA0
variable <- pheno_m_PC_wbc[,colnames(pheno_m_PC_wbc) %in% c("Case")]
var_idx <- order(rowVars(norm_beta_m_fil, na.rm=T), decreasing=T)[1:50000]
beta_isva <- meffil:::impute.matrix(norm_beta_m_fil[var_idx,,drop=F])
isva0 <- meffil:::DoISVA(beta_isva, variable, verbose=T)
isva0 <- as.data.frame(isva0$isv)

        
#ISVA1
covs <- pheno_m_PC_wbc[,colnames(pheno_m_PC_wbc) %in% c("Slide","Sample_Plate","sentrix_row","sentrix_col","Sex","Bcell","CD4T","CD8T","Eos","Mono","NK")]
factor_log <- sapply(cbind(isva0, covs), is.factor)
isva1 <- meffil:::DoISVA(beta_isva, variable, 
                        cf.m=cbind(isva0, covs), 
                        factor.log=factor_log,
                        verbose=T)
isva1 <- as.data.frame(isva1$isv)

design_isva0 <- model.matrix(~ .,cbind(isva0,variable))
design_isva1 <- model.matrix(~ .,cbind(isva1,variable))



#^^^^^^^^^^^^^^^^^^^^^^^^
#create design matrices
#^^^^^^^^^^^^^^^^^^^^^^^^

#PCs 1,2,3 + Sex + WBCs
design_PC <- model.matrix(~ PC1 + PC2 + PC3 + Sex + Bcell + CD4T + CD8T + Eos + Mono + Neu + NK + Case, sample_sheet)
all(colnames(norm_mval_fil) == rownames(design_PC))

#PCs (as above) + leave out a blood cell type
design_PC_trim_neu <- model.matrix(~ PC1 + PC2 + PC3 + Sex + Bcell + CD4T + CD8T + Eos + Mono + NK + Case, sample_sheet)

#all and no PCs
design_all <- model.matrix(~ Sample_Plate + Slide + sentrix_row + sentrix_col + Sex + Bcell + CD4T + CD8T + Eos + Mono + Neu + NK + Case, sample_sheet)

#^^^^^^^^^^^^^^^^^^^^^^^^
#run limma using PCs
#^^^^^^^^^^^^^^^^^^^^^^^^
#PCs 1,2,3 + Sex + WBCs
DMPs_PC <- lmFit(norm_mval_fil,design_PC)
DMPs_PC <- eBayes(DMPs_PC)
res_DMPs_PC <- topTable(DMPs_PC, coef = "Case", number = Inf, genelist=anno_sub, sort.by="B")

#PCs (as above) + leave out a blood cell type
#TODO

#all and no PCs
#TODO




#^^^^^^^^^^^^^^^^^^^^^^^^
#Plots
#^^^^^^^^^^^^^^^^^^^^^^^^

#scree
ggplot(pcpvar[1:10,], aes(x=PC,y=per_var)) +  
  geom_bar(stat = "identity") + 
  ylab("percentage variance explained") +
  theme_bw() + 
  theme(panel.grid.major = element_line(colour = "#F6F7F7"),
        panel.grid.minor = element_line(colour = "#F6F7F7"))

#plot PCs  
pcs_df <- cbind(sample_sheet, pcs$x[,1:4])

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
