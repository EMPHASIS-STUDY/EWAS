#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
# batchPCAcorr <- R function to check for batch variable correlation with principal  
# components constructed from 450k methylation data
#
# input : pcs = matrix of derived principal components, batch = matrix of batch variables to be tested against using:
#
# 1. Wilcoxon test for categorical variables with 2 levels
# 2. Kruskal Wallace test for categorical variables with > 2 levels
# 3. Spearman's correlation for continuous variables
#
# output : matrix of association p-values
#
# author. Ayden Saffari <ayden.saffari@lshtm.ac.uk>
# affiliations. MRC ING, LSHTM
#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
  
batchPCAcorr <- function(pcs=NULL,batch=NULL,npcs=10){  
	res <- matrix(nrow=ncol(batch),ncol=npcs,dimnames=list(colnames(batch),colnames(pcs)))
	
	#ensure samples are in the same order
	if(!(all.equal(rownames(pcs),rownames(batch)))){
	stop("samples must the same in both PC and batch variable data (i.e. rows should be equivalent)")
	}
	
	#iterate through batch vars and PCs running the relevant test for association
	for(m in rownames(res)){
		cat <- ifelse(length(levels(batch[,m]))==2,TRUE,FALSE)
		for(n in colnames(res)){
			if(cat){
				res[m,n] <- wilcox.test(pcs[,n] ~ batch[,m])$p.value
			}else if(!(cat)){
				if(typeof(levels(batch[,m]))!= "NULL"){
					res[m,n] <- kruskal.test(pcs[,n] ~ 	batch[,m])$p.value
				}else{
				    #fit <- lm(pcs[,n] ~ batch[,m])
					#res[m,n] <- summary.lm(fit)$coefficients[2,4]
					res[m,n] <- cor.test(pcs[,n],batch[,m],method='spearman',
				                exact=F,continuity=T)$p.value	
				}
			}
		}
	}
	return(res)
}  


