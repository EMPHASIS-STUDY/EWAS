#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
# R script for emphasis candidate genes analysis - robust regressions and regional method
#
# author. Ayden Saffari <ayden.saffari@lshtm.ac.uk> (MRC ING, LSHTM)
#               
#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/

library("ggplot2")
library("reshape")
library("dplyr")
library("plyr")
library("broom")
library("gtools")
library("ggplot2")
library("ggfortify")
library("ggdendro")
library("gplots")
library("ggthemes")
library("RColorBrewer")
library("limma")
library("devtools")
install_github("asaffa/gamplotlib")
library("gamplotlib")
library("MASS")
library("sandwich")
library("lmtest")
library("car")
library("EmpiricalBrownsMethod")

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#init
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#read in pdata
pdata <- readRDS("../../EPIC/R_objects/EMPH_EPIC_sample_sheet.rds")
load("../../EPIC/R_objects/cell_counts.RData")
pdata <- merge(pdata,cell_counts,by.x="Sample_Name",by.y="IID",all.x=T)
pdata <- pdata[-grep("ILC",pdata$Subject_ID),]
pdata_redux <- pdata[,colnames(pdata) %in% c("Subject_ID","sex","BMI","age",
					 "MasterGroupNo","MooreSoC","MaternalBMI","Bcell","CD4T",
					 "CD8T","Eos","Mono","Neu","NK")]

#read in pyro data
csv_files <- list.files("../data/GMB/",pattern="*.csv",full.names=T)
all_candi <- lapply(csv_files,read.csv)
names(all_candi)  <- gsub(".csv","",list.files("../data/GMB",pattern="*.csv"),fixed=T)
all_candi <- lapply(all_candi,function(x){x$batch <- as.factor(x$batch)
                    x})

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#run robust regressions 
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

all_res_rlm <- list()

for(i in 1:length(all_candi)){
    message("#################")
	message(names(all_candi)[i])
	message("#################")
    
	combi_pyro_pdata <- merge(pdata_redux, all_candi[[i]], 
	                        by.x="Subject_ID",by.y="sample_id")
	combi_pyro_pdata <- combi_pyro_pdata %>% na.omit	  
    combi_pyro_pdata <- droplevels(combi_pyro_pdata)                        
    print(dim(combi_pyro_pdata))
	res <- list()

	for (probe in                          
		colnames(combi_pyro_pdata)[grep("cg",colnames(combi_pyro_pdata))]){
		design <- as.formula(paste0(probe,
		                     " ~ batch + age + sex + 
		                     Bcell + CD4T + CD8T + Eos + Mono + NK +
		                     MooreSoC + MasterGroupNo"))
		rlmmod <-  rlm(design,data=combi_pyro_pdata)
		print(summary(rlmmod))
		res[[probe]] <- coeftest(rlmmod, vcov=vcovHC(rlmmod, type="HC0"))
		 
	}

	all_res_rlm[names(all_candi)[i]] <- list(res)    
	
	         
}

lapply(all_res_rlm,function(x){lapply(x,function(y){y["MasterGroupNo2",]})})

all_res_rlm_tidy <- lapply(all_res_rlm,function(x){lapply(x,function(y){tidy(y)})})

#write results
#############
rlm_res_file <- "../results/all_candis_rlm_res.txt"
for (i in 1:length(all_res_rlm_tidy)){
	cat(names(all_res_rlm_tidy)[i],"\n",
		file=rlm_res_file, append=TRUE)
	for (j in 1:length(all_res_rlm_tidy[[i]])){
		cat(names(all_res_rlm_tidy[[i]])[j],"\n",
			file=rlm_res_file, append=TRUE)
		write.table(all_res_rlm_tidy[[i]][[j]],
					file=rlm_res_file, append=TRUE,
           			quote=F, row.names=F)
		cat("\n",file=rlm_res_file, append=TRUE)
	}
}

#QQ plot for all candis
rlm_pvals <- unlist(lapply(all_res_rlm,function(x){lapply(x,function(y){
                    y["MasterGroupNo2","Pr(>|z|)"]})}))
ggQQplot(rlm_pvals,ylim=c(0,2)) + theme_gamplotlib()
ggsave(filename="../results/rlm_all_candis_qqplot.png")



#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#regional pvalue combination
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

all_pvals_rlm <- lapply(all_res_rlm,function(x){unlist(lapply(x,function(y){
                   y["MasterGroupNo2",4]}))})
eb_combp_res_rlm <- vector(mode="numeric")
for(i in 1:length(all_res_rlm)){
	eb_combp_res_rlm <- c(eb_combp_res_rlm,
	                  empiricalBrownsMethod(t(all_candi[[i]][,
	                  colnames(all_candi[[i]]) %in% 
	                  names(all_pvals_rlm[[i]])] %>% na.omit), all_pvals_rlm[[i]],
	                  extra_info=T)$P_test)
}

names(eb_combp_res_rlm) <- names(all_candi)
p.adjust(eb_combp_res_rlm,method="fdr")

combi_reg_p_rlm <- cbind(as.data.frame(eb_combp_res_rlm),
			as.data.frame(p.adjust(eb_combp_res_rlm,method="fdr")))
colnames(combi_reg_p_rlm)[1] <- "p.value"
colnames(combi_reg_p_rlm)[2] <- "fdr"		
write.table(combi_reg_p_rlm,file="../results/regional_p_fdr.txt",quote=F)
