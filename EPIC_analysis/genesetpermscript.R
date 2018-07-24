#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
# Emphasis DMP comparison:
# Goeman Buhlmann self contained geneset permutation test
#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
require("doParallel")
registerDoParallel(cores=4)


##########
#load data
##########
pdata <- readRDS("../R_objects/pdata.rds")
res_DMPs_pcs <- readRDS("../R_objects/res_DMPs_pcs.rds")
DMRs_CpGs <- readRDS("../R_objects/EMPH_GMB_DMRs_CpGs.rds")
norm_beta_fil <- readRDS("../R_objects/norm_beta_fil.rds")
GMB_CpGs <- norm_beta_fil[which(rownames(norm_beta_fil) %in%
                          unique(c(res_DMPs_pcs$Name[
                          res_DMPs_pcs$adj.P.Val < 0.1],DMRs_CpGs))),]
pcs <- readRDS("../R_objects/pcs.rds")

###################
# helper functions
##################
#function to add column for delta Beta
add_delta_B <- function(res,beta,inter){
  beta_res <- beta[match(res$Name,rownames(beta)),]
  delta_B <- apply(beta_res,1, function(x){
                      mean(x[which(inter == 2)]) -  mean(x[which(inter == 1)])})
  delta_B <- cbind(delta_B)
  res <- merge(res,delta_B,by.x='row.names',by.y='row.names')
  res <- res[with(res,order(-B)),]
  return(res)
}

genesetpermtest <- function(perm){
	val <- tryCatch(
	{
	design_pcs <- model.matrix(formula(paste0(
	"~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 +
	PC15 + Age + MooreSoC + ",names(perm))),cbind(cbind(pdata,pcs),perm))
	DMPs_pcs <- lmFit(mvals,design_pcs)
	DMPs_pcs <- eBayes(DMPs_pcs)
	res_DMPs_pcs <- topTable(DMPs_pcs, coef = paste0(names(perm),"2"), number = Inf,
                         genelist=anno, sort.by="B")
	res_DMPs_pcs <- res_DMPs_pcs[res_DMPs_pcs$P.Value < 0.05,]
	res_DMPs_pcs <- add_delta_B(as.data.frame(res_DMPs_pcs),
								betas,perm)
	res_DMPs_pcs <- res_DMPs_pcs[abs(res_DMPs_pcs$delta_B) >= 0.02,]
	#return size of overlap
	val <- nrow(res_DMPs_pcs)
	},
	error=function(x){return(0)},
	warning=function(x){return(0)}
	)
	return(val)
}


############################
# geneset permutation test
#######################

#generate permutations
nperms <- 9999
permMasterGroupNo <- as.data.frame(matrix(nrow=length(
                                   pdata$MasterGroupNo),ncol=nperms))
for(i in 1:nperms){
	permMasterGroupNo[,i] <- sample(pdata$MasterGroupNo,nrow(permMasterGroupNo))
}
permMasterGroupNo$obs <- pdata[,14]

#####
#MEs
#####
betas <- norm_beta_fil[rownames(norm_beta_fil) %in% MEs,]
mvals <- norm_mval_fil[rownames(norm_mval_fil) %in% MEs,]
anno <- anno_sub[anno_sub$Name %in% rownames(mvals),]

perm_res <- foreach(i=1:ncol(permMasterGroupNo),.combine=c, .verbose=T) %dopar% {
        			genesetpermtest(permMasterGroupNo[,i,drop=F])
        			}
(length(which(perm_res[1:9999] >= perm_res[10000]))) / 10000

#####
#ICRs
#####

betas <- norm_beta_fil[rownames(norm_beta_fil) %in% ICR_cpgs,]
mvals <- norm_mval_fil[rownames(norm_mval_fil) %in% ICR_cpgs,]
anno <- anno_sub[anno_sub$Name %in% rownames(mvals),]

perm_res <- foreach(i=1:ncol(permMasterGroupNo),.combine=c, .verbose=T) %dopar% {
        			genesetpermtest(permMasterGroupNo[,i,drop=F])
        			}
perm_res[10000]
(length(which(perm_res[1:9999] >= perm_res[10000]))) / 10000

#########
#nutripgs
#########

betas <- norm_beta_fil[rownames(norm_beta_fil) %in% nutricpgs,]
mvals <- norm_mval_fil[rownames(norm_mval_fil) %in% nutricpgs,]
anno <- anno_sub[anno_sub$Name %in% rownames(mvals),]

perm_res <- foreach(i=1:ncol(permMasterGroupNo),.combine=c, .verbose=T) %dopar% {
        			genesetpermtest(permMasterGroupNo[,i,drop=F])
        			}
perm_res[10000]
(length(which(perm_res[1:9999] >= perm_res[10000]))) / 10000


#var controls
betas <- norm_beta_fil[rownames(norm_beta_fil) %in% var_controls,]
mvals <- norm_mval_fil[rownames(norm_mval_fil) %in% var_controls,]
anno <- anno_sub[anno_sub$Name %in% rownames(mvals),]

perm_res <- foreach(i=1:ncol(permMasterGroupNo),.combine=c, .verbose=T) %dopar% {
        			genesetpermtest(permMasterGroupNo[,i,drop=F])
        			}
(length(which(perm_res[1:9999] >= perm_res[10000]))) / 10000
