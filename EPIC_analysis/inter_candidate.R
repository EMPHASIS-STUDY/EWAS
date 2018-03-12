#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
# R script for performing pyroseq candidate genes analysis
#
# author. Ayden Saffari <ayden.saffari@lshtm.ac.uk> (MRC ING, LSHTM)
#               
#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/

library("ggplot2")
library("reshape")
library("plyr")
library("gtools")
library("ggplot2")
library("ggdendro")
library("gplots")
library("ggthemes")
library("RColorBrewer")
library("limma")

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# RXRA
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

####
#init
####

#read in data
RXRA_SARAS <- read.csv("RXRA_SARAS_cg1-14_clean07.03.18.csv")
RXRA_SARAS$plate <- as.factor(RXRA_SARAS$plate)

#subset : processed cpgs only
RXRA_SARAS_cg_c <- RXRA_SARAS[,grep("cg[0-9]{1,}_c",colnames(RXRA_SARAS))]
rownames(RXRA_SARAS_cg_c) <- RXRA_SARAS$sample

#subset : processed cpgs with outliers removed
RXRA_SARAS_cg_all_o <- RXRA_SARAS[,grep("cg[0-9]{1,}_all_o",colnames(RXRA_SARAS))]
RXRA_SARAS_cg_all_o <- cbind(RXRA_SARAS_cg_c[,!(colnames(RXRA_SARAS_cg_c) %in% c("cg2_c","cg11_c","cg12_c"))],RXRA_SARAS_cg_all_o)
RXRA_SARAS_cg_all_o <- RXRA_SARAS_cg_all_o[,
						match(mixedsort(levels(melt_RXRA_SARAS_cg_all_o$CpG)),
						colnames(RXRA_SARAS_cg_all_o))]
												
#batch adjust data
design <-  model.matrix(~ sex + age + group + plate, RXRA_SARAS,
                          contrasts.arg = list(sex="contr.sum",
                          group="contr.treatment",plate="contr.sum"))

RXRA_SARAS_cg_all_o_adj <- lmFit(t(RXRA_SARAS_cg_all_o),design)
betas <- RXRA_SARAS_cg_all_o_adj$coefficients[,-c(1,2,3,4)]
RXRA_SARAS_cg_all_o_adj <- t(RXRA_SARAS_cg_all_o) - (betas %*% t(design[,-c(1,2,3,4)]))
RXRA_SARAS_cg_all_o_adj <- t(RXRA_SARAS_cg_all_o_adj)


#check batch adjustment has worked
summary(lm(cg9_c ~ sex + age + group + plate, data=cbind(RXRA_SARAS_cg_all_o,RXRA_SARAS[,colnames(RXRA_SARAS) %in% c("sex","age","group","plate")])))

summary(lm(cg9_c ~ sex + age + group + plate, data=cbind(RXRA_SARAS_cg_all_o_adj,RXRA_SARAS[,colnames(RXRA_SARAS) %in% c("sex","age","group","plate")])))

##############
#density plots
##############

#processed
melt_RXRA_SARAS_cg_c <- melt(RXRA_SARAS_cg_c)
colnames(melt_RXRA_SARAS_cg_c)[1] <- "CpG"

ggplot(melt_RXRA_SARAS_cg_c,aes(value,color=CpG)) + geom_density() + scale_color_manual(values=c(brewer.pal(8, "Dark2"),"#1f78b4","#33a02c",
"#e31a1c","#ff7f00")) + theme_gamplotlib()
ggsave("RXRA_SARAS_cg_c_density.pdf",width=7,height=7)

#processed with additional outlier removal
melt_RXRA_SARAS_cg_all_o <- melt(RXRA_SARAS_cg_all_o)
colnames(melt_RXRA_SARAS_cg_all_o)[1] <- "CpG"

ggplot(melt_RXRA_SARAS_cg_all_o,aes(value,color=CpG)) + geom_density() + scale_color_manual(values=c(brewer.pal(8, "Dark2"),"#1f78b4","#33a02c",
"#e31a1c","#ff7f00")) + theme_gamplotlib()
ggsave("RXRA_SARAS_cg_all_o_density.pdf",width=7,height=7)


#####################
#correlation heatmaps
#####################

col_heat <- c("#88419d","#8c96c6","#b3cde3","#edf8fb","#edf8fb","#b2e2e2","#66c2a4","#238b45")

#processed with additional outlier removal
pdf(file="RXRA_SARAS_corr_heat_cg_all_o.pdf")
heatmap.2(cor(RXRA_SARAS_cg_all_o,use="complete.obs"), trace="none", col = col_heat, margin=c(8, 8),keysize=1)
dev.off()

#processed with additional outlier removal and batch adjustment
pdf(file="RXRA_SARAS_corr_heat_cg_all_o_batch_adjusted.pdf")
heatmap.2(cor(RXRA_SARAS_cg_all_o_adj_2,use="complete.obs"), trace="none", col = col_heat, margin=c(8, 8),keysize=1)
dev.off()


#####################################
#run regressions (exclude cgs 13 and 14)
#####################################

#######################
#cluster-based analysis
#######################
