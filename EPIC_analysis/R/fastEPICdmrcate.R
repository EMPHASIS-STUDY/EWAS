#TODO - make compatible with meffil.ewas object
#TODO - test
#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
# mefdmrcate <- R function to perform regional analysis on methylation data processed 
# by meffil using the dmrcate method
#
# input : matrix of beta values, meffil.ewas object
# output : list of differentially methylated regions from DMRcate
#
#          Annotation argument for EPIC arrays (default):
#          c(array = "IlluminaHumanMethylationEPIC", annotation ="ilm10b2.hg19")
#
#		   Annotation argument for 450K arrays:
#		   c(array = "IlluminaHumanMethylation450k", annotation ="ilmn12.hg19")
#
# author. Ayden Saffari <ayden.saffari@lshtm.ac.uk>
# affiliations. MRC ING, LSHTM
#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/


fastEPICdmrcate <- function(mvals,design, 
 					   annotation = c(array = "IlluminaHumanMethylationEPIC",
 					   annotation ="ilm10b2.hg19"),analysis.type="differential",
 					   coef="variable",fdr=0.05, lambda=500, C=2, pcutoff=0.05,
 					   min.cpgs=2, mc.cores=1,...) {
	
	require("DMRcate")	
	dm <- cpg.annotate(datatype="array", mvals, annotation=annotation,design=design,
	 				  analysis.type=analysis.type,fdr=fdr,coef=coef,...)
	DMRs <- dmrcate(dm,lambda=lambda, C=C, pcutoff=pcutoff, min.cpgs=min.cpgs,
	mc.cores=mc.cores)
	return(DMRs)
}
