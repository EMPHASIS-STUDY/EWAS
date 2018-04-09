
#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
# R script to select subsample for technical validation by pyroseq
# so as to cover entire range of methylation values for selected hits
#
# author. Ayden Saffari <ayden.saffari@lshtm.ac.uk> (MRC ING, LSHTM)
#           
# NOT FOR DISTRIBUTION/ PUBLIC CONSUMPTION
#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/

library("reshape")
library("ggplot2")
library("gamplotlib")
library("dplyr")
library("lumi")

colpalricht1 <- c("#FA980D","#0E7624","#F97E86","#0F7597","#B477AB",
                  "#834041","#6D6F72","#2C191F")
      
colpalmin <- c("#303030","#F0516F")
                  
pdata <- readRDS("../../EPIC/R_objects/pdata.rds")
norm_mval_fil <- readRDS("../../EPIC/R_objects/norm_mval_fil.rds")
ESM1_loc <- c("cg14972155","cg06837426","cg20451680","cg20673840")
LZTS1_loc <- c("cg27655507")

norm_mval_fil_pyro <- norm_mval_fil[rownames(norm_mval_fil) %in% c(ESM1_loc,LZTS1_loc),]
norm_beta_fil_pyro <- m2beta(norm_mval_fil_pyro)

########################
#plot starting densities
########################
pyro_melt <- melt(norm_mval_fil_pyro)
colnames(pyro_melt) <- c("CpG","sample","M")
ggplot(pyro_melt, aes(M, colour = CpG)) + ylim(0,2) + geom_density() + 
	  theme_gamplotlib() + scale_color_gamplotlib()
ggsave("../results/EMPH_pyro_valid_start_density.pdf",height=7,width=7)

###################################################
#select subsample of 96 samples (48 from each arm)
# covering range of methylation values
###################################################
#select samples with in and max values across all cpgs
minmaxsamps <- unique(c(apply(norm_mval_fil_pyro,1,which.min),
                      apply(norm_mval_fil_pyro,1,which.max)))
                      
remMaster1 <- which(pdata$MasterGroupNo[minmaxsamps] == 1)
remMaster2 <- which(pdata$MasterGroupNo[minmaxsamps] == 2)

#split remaining sample into cases and controls
cases <- which(pdata$MasterGroupNo == 1)
cases <- cases[!(cases %in% minmaxsamps[remMaster1])]
ctrls <- which(pdata$MasterGroupNo == 2)
ctrls <- ctrls[!(ctrls %in% minmaxsamps[remMaster2])]
                    
#generating function
swapsampleoptim <- function(subs){
    cases <- cases
    controls <- ctrls 
    n1 <- sample(seq(0:2),1)
    n2 <- sample(seq(0:2),1)
    add_rand_case <- sample(cases[!(cases %in% subs)],n1)
    add_rand_cont <- sample(controls[!(controls %in% subs)],n2)
    rem_rand_case <- sample(cases[(cases %in% subs)],n1)
    rem_rand_cont <- sample(controls[(controls %in% subs)],n2)
    subs <- c(subs[!(subs %in% c(rem_rand_case,rem_rand_cont))],
                  add_rand_case,add_rand_cont)
    return(subs)
}

#target uniform distribution
#target normal
targetuni <- apply(norm_mval_fil_pyro,1,function(x){runif(289,min(x),max(x))})
targetnorm <- apply(norm_mval_fil_pyro,1,function(x){rnorm(289, mean = mean(x), sd = sd(x))})
plot(density(targetnorm),xlim=c(-5,5))
lines(density(targetuni))

#scoring function (sum squared error from deciles of uniform)
subsgetdistoptimuni <- function(subs){
    x <- norm_mval_fil_pyro
    y <- targetuni
    dist <- vector()
    qsx <- list()
    qsy <- list()
    for(i in rownames(x)){
    	xsel <- x[rownames(x) == i,subs]
    	ysel <- y[,colnames(y) == i]
		qsx[[i]] <- quantile(xsel,c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9))
		qsy[[i]] <- quantile(ysel,c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9))
		dist <- c(dist,(qsx[[i]] - qsy[[i]])^2)
		
	}
	return(sum(dist))
}
                   
#simulated annealing to optimize subsample selection
#based on matching deciles of uniform distribution
temp <- optim(subsgetdistoptimuni, par=c(sample(cases,48 - length(remMaster1)),
              sample(ctrls,48 - length(remMaster2))), gr=swapsampleoptim, 
              method="SANN", control=list(maxit=50000, temp=10000,tmax=1000,
              trace=TRUE,REPORT=1))
                   
##############
#plot results
###############

#quick QC check
head(norm_mval_fil[rownames(norm_mval_fil) %in% c(ESM1_loc,LZTS1_loc),
     c(temp$par,minmaxsamps)])
plot(density(norm_mval_fil[rownames(norm_mval_fil) %in% c(ESM1_loc,LZTS1_loc),
     c(temp$par,minmaxsamps)]),xlim=c(-1,3),ylim=c(0,2))
apply(norm_mval_fil[rownames(norm_mval_fil) %in% c(ESM1_loc,LZTS1_loc),
      c(temp$par,minmaxsamps)],1,function(x){lines(density(x))})

#plot densities
pyro_subs_melt <- melt(norm_mval_fil_pyro[,c(temp$par,minmaxsamps)])
colnames(pyro_subs_melt) <- c("CpG","sample","M")
ggplot(pyro_subs_melt, aes(M, colour = CpG)) + ylim(0,2) + geom_density() + 
	  theme_gamplotlib() + scale_color_gamplotlib()
ggsave("../results/EMPH_pyro_valid_subs_density.pdf",height=7,width=7)

#plot histograms
pyro_melt_beta <- melt(norm_beta_fil_pyro)
colnames(pyro_melt_beta) <- c("CpG","sample","Beta")

pyro_subs_melt_beta <- melt(norm_beta_fil_pyro[,c(temp$par,minmaxsamps)])
colnames(pyro_subs_melt_beta) <- c("CpG","sample","Beta")
cpgdistohisto <- list()
for (i in 1:length(levels(pyro_melt_beta$CpG))){
	p1 <- eval(substitute(ggplot() + geom_histogram(data=pyro_melt_beta[
  						  pyro_melt_beta$CpG == levels(pyro_melt_beta$CpG)[i],],
  						  aes(x=Beta,fill="all"),alpha=0.25, binwidth=0.025,
  						  color=alpha(colpalmin[1],0.5)) +
 						  geom_histogram(data=pyro_subs_melt_beta[
 						  pyro_subs_melt_beta$CpG ==levels(pyro_melt_beta$CpG)[i],],
 						  aes(x=Beta,fill="subsample"),
 						  alpha=0.5, binwidth=0.025, color=alpha(colpalmin[2],0.75)) +
 						  theme_gamplotlib() + 
 						  scale_color_manual(values=c(all=colpalmin[1],
 						  subsample=colpalmin[2]),name="Dataset") +
 						  scale_fill_manual(values=c(all=colpalmin[1],
 						  subsample=colpalmin[2]),name="Dataset") +
 						  xlim(0,1) + ylim(0,80) +
 						  ggtitle(levels(pyro_melt_beta$CpG)[i])
 						  ,list(i = i)))
    cpgdistohisto[[i]] <- p1
}  

pdf(file="../results/EMPH_pyro_valid_subs_histo.pdf",height=10,width=10)
multiplot(plotlist = cpgdistohisto, layout = matrix(seq(1:6),ncol=2,nrow=3,byrow=T))
dev.off()

####################
# write summary stats
####################
write.csv(apply(norm_mval_fil_pyro,1,summary),
		   "../results/EMPH_pyro_valid_summstat_all.csv")
write.csv(apply(norm_mval_fil_pyro[,c(temp$par,minmaxsamps)],1,summary),
		  "../results/EMPH_pyro_valid_summstat_subs.csv")
##############
#write results
##############
write.csv(pdata[c(temp$par,minmaxsamps),],quote=T,
          file="../results/EMPH_pyro_valid_subsselect_pdata.csv")
