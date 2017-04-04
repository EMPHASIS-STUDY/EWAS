#UNFINISHED
#TEST SCRIPT

#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
# TEST R script to perform EWAS using limma 
#
# #norm_betas = normalized beta values
# #sample_sheet = list of covariates
#
#
# author. Ayden Saffari <ayden.saffari@lshtm.ac.uk>
# affiliations. MRC ING, LSHTM
#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/

library("limma")
library("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
library("gridExtra")
library("ggplot2")
source("ggPCAplot.R")
source("ilEPICfilter.R")
#^^^^^^^^^^^^^^^^^^^^^^^^
#initialization
#^^^^^^^^^^^^^^^^^^^^^^^^

all(colnames(norm_betas) == rownames(sample_sheet))

#filter probes
#TODO
norm_betas_fil <- ilEPICfilter(norm_betas, cross=T, multi=T, sex=T) 

#convert to m vals
norm_mval_fil <- logit2(norm_betas_fil)

#get annotation
anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
anno_sub <- anno[match(rownames(norm_mval_m_fil),anno$Name),
c(1:4,12:19,24:ncol(anno))]

#^^^^^^^^^^^^^^^^^^^^^^^^
#Principal components
#^^^^^^^^^^^^^^^^^^^^^^^^

pcs <- prcomp(t(norm_mval_m_fil),retx=T, center=T, scale=T )
#get variance explained
summary(pcs)

pal <- c("#F87E86","#0E7597","#2E7B39","#FDA80A",
"#B31918","#9CD893","#E6BD9B","#B7BAAA","#683A2E","#211318")

#Sex
pdf(file="../results/ENID_PCS_Sex.pdf", onefile=FALSE)
p1 <- ggpcaplot(pcs_df, pcs=c(1,2), batch="Sex",pal=pal) + theme_bw() + coord_fixed(ratio=1) + theme(legend.position="none")
p2 <- ggpcaplot(pcs_df, pcs=c(1,3), batch="Sex",pal=pal) + theme_bw() + coord_fixed(ratio=1) + theme(legend.position="none")
p3 <- ggpcaplot(pcs_df, pcs=c(2,3), batch="Sex",pal=pal) + theme_bw() + coord_fixed(ratio=1) + theme(legend.position="none")
p3l <- ggpcaplot(pcs_df, pcs=c(2,3), batch="Sex",pal=pal) + theme_bw() + coord_fixed(ratio=1) 

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
p1 <- ggpcaplot(pcs_df, pcs=c(1,2), batch="sentrix_row",pal=pal) + theme_bw() + coord_fixed(ratio=1) + theme(legend.position="none")
p2 <- ggpcaplot(pcs_df, pcs=c(1,3), batch="sentrix_row",pal=pal) + theme_bw() + coord_fixed(ratio=1) + theme(legend.position="none")
p3 <- ggpcaplot(pcs_df, pcs=c(2,3), batch="sentrix_row",pal=pal) + theme_bw() + coord_fixed(ratio=1) + theme(legend.position="none")
p3l <- ggpcaplot(pcs_df, pcs=c(2,3), batch="sentrix_row",pal=pal) + theme_bw() + coord_fixed(ratio=1) 

legr1 <-g_legend(p3l)
grid.arrange(p1,p2,p3,legr1,nrow=1,heights=c(2),widths=c(2,2,2,1.1))
dev.off()

#Sample plate
pdf(file="../results/ENID_PCS_Sample_Plate.pdf", onefile=FALSE)
p1 <- ggpcaplot(pcs_df, pcs=c(1,2), batch="Sample_Plate",pal=pal) + theme_bw() + coord_fixed(ratio=1) + theme(legend.position="none")
p2 <- ggpcaplot(pcs_df, pcs=c(1,3), batch="Sample_Plate",pal=pal) + theme_bw() + coord_fixed(ratio=1) + theme(legend.position="none")
p3 <- ggpcaplot(pcs_df, pcs=c(2,3), batch="Sample_Plate",pal=pal) + theme_bw() + coord_fixed(ratio=1) + theme(legend.position="none")
p3l <- ggpcaplot(pcs_df, pcs=c(2,3), batch="Sample_Plate",pal=pal) + theme_bw() + coord_fixed(ratio=1) 
legr1 <-g_legend(p3l)
grid.arrange(p1,p2,p3,legr1,nrow=1,heights=c(2),widths=c(2,2,2,1.1))
dev.off()

#TODO - remaining vars

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
# TODO - try with ISVA
#^^^^^^^^^^^^^^^^^^^^^^^^

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
# TODO - try with SVs
#^^^^^^^^^^^^^^^^^^^^^^^^
#TODO
