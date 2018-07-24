#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
#
# R script to perform self contained geneset permutation test
# author. Ayden Saffari <ayden.saffari@lshtm.ac.uk> (MRC ING, LSHTM)
#
# NOT FOR DISTRIBUTION/ PUBLIC CONSUMPTION
#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/

library("limma")
library("doParallel")
registerDoParallel(cores=4)

##########
#load data
##########
pdata <- readRDS("../R_objects/pdata.rds")
pcs <- readRDS("../R_objects/pcs.rds")
norm_beta_fil <- readRDS("../R_objects/norm_beta_fil.rds")
norm_mval_fil <- readRDS("../R_objects/norm_mval_fil.rds")

MEs <- readRDS("../R_objects/MEs.rds")
ICR_cpgs <- readRDS("../R_objects/ICR_cpgs.rds")
nutricpgs <- readRDS("../R_objects/nutricpgs.rds")
var_controls <- readRDS("../R_objects/var_controls.rds")

#merge pcs into pdata
pdata <- cbind(pdata,pcs)

###################
# helper functions
##################
#function to add column for delta Beta
add_delta_B <- function(res,beta,inter){
  beta_res <- beta[match(rownames(res),rownames(beta)),]
  delta_B <- apply(beta_res,1, function(x){
                      mean(x[which(inter == 2)]) -  mean(x[which(inter == 1)])})
  delta_B <- cbind(delta_B)
  res <- merge(res,delta_B,by.x='row.names',by.y='row.names')
  res <- res[with(res,order(-B)),]
  return(res)
}

#geneset permutation testing function
genesetpermtest <- function(perm,form){
  val <- tryCatch(
	{
	design <- model.matrix(formula(paste0(form,names(perm))),
                             cbind(pdata,perm))
	DMPs <- lmFit(mvals,design)
	DMPs <- eBayes(DMPs)
	res_DMPs <- topTable(DMPs, coef = paste0(names(perm),"2"),
                           number = Inf, sort.by="B")
	res_DMPs <- res_DMPs[res_DMPs$P.Value < 0.05,]
	res_DMPs <- add_delta_B(as.data.frame(res_DMPs),
								betas,perm)
	res_DMPs <- res_DMPs[abs(res_DMPs$delta_B) >= 0.02,]

	#return size of overlap
	val <- nrow(res_DMPs)
	},
	error=function(x){return(0)},
	warning=function(x){return(0)}
	)
	return(val)
}

###################################
# perform geneset permutation tests
###################################

#model formula for limma, should be same as used in the EWAS except with
#mastergroupnumber/samplegroup/intervention removed (leave the trailing
# + and space)
modform <- "~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC15 + Age + MooreSoC + "

#generate permutations
nperms <- 9999
permMasterGroupNo <- as.data.frame(matrix(nrow=length(
                                   pdata$MasterGroupNo),ncol=nperms))
for(i in 1:nperms){
	permMasterGroupNo[,i] <- sample(pdata$MasterGroupNo,nrow(permMasterGroupNo))
}
permMasterGroupNo$obs <- pdata[,"MasterGroupNo"]

#####
#MEs
#####
betas <- norm_beta_fil[rownames(norm_beta_fil) %in% MEs,]
mvals <- norm_mval_fil[rownames(norm_mval_fil) %in% MEs,]

perm_res <- foreach(i=1:ncol(permMasterGroupNo),.combine=c, .verbose=T) %dopar%
                   {
        			          genesetpermtest(permMasterGroupNo[,i,drop=F],modform)
        			     }
message("##### MEs #####")
message(paste0("overlap: ",perm_res[10000]))
message(paste0("p-value overlap: ",
        (length(which(perm_res[1:9999] >= perm_res[10000]))) / 10000))

#####
#ICRs
#####
betas <- norm_beta_fil[rownames(norm_beta_fil) %in% ICR_cpgs,]
mvals <- norm_mval_fil[rownames(norm_mval_fil) %in% ICR_cpgs,]

perm_res <- foreach(i=1:ncol(permMasterGroupNo),.combine=c, .verbose=T) %dopar%
                   {
        			          genesetpermtest(permMasterGroupNo[,i,drop=F],modform)
        			     }
message("##### ICRs #####")
message(paste0("overlap: ",perm_res[10000]))
message(paste0("p-value overlap: ",
        (length(which(perm_res[1:9999] >= perm_res[10000]))) / 10000))

#########
#nutripgs
#########
betas <- norm_beta_fil[rownames(norm_beta_fil) %in% nutricpgs,]
mvals <- norm_mval_fil[rownames(norm_mval_fil) %in% nutricpgs,]

perm_res <- foreach(i=1:ncol(permMasterGroupNo),.combine=c, .verbose=T) %dopar%
                    {
        			           genesetpermtest(permMasterGroupNo[,i,drop=F],modform)
        			      }

message("##### nutricpgs #####")
message(paste0("overlap: ",perm_res[10000]))
message(paste0("p-value overlap: ",(length(
        which(perm_res[1:9999] >= perm_res[10000]))) / 10000))

##############
#var controls
##############
betas <- norm_beta_fil[rownames(norm_beta_fil) %in% var_controls,]
mvals <- norm_mval_fil[rownames(norm_mval_fil) %in% var_controls,]

perm_res <- foreach(i=1:ncol(permMasterGroupNo),.combine=c, .verbose=T) %dopar%
                    {
        			         genesetpermtest(permMasterGroupNo[,i,drop=F],modform)
        			      }

message("##### variable controls #####")
message(paste0("overlap: ",perm_res[10000]))
message(paste0("p-value overlap: ",(length(
        which(perm_res[1:9999] >= perm_res[10000]))) / 10000))
