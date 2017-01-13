#TODO - adapt for EPIC
#TODO - test
#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
# il450Kfilter <- R function to filter 450K methylation data in order to remove SNP,
# cross-reactive/non-specific, multi-mapping, and sex specific probes
#
# adapted from the SNPfilter function in minfi
#
# input : matrix of beta or m values, with rownames = probe identfiers (cgXXXXXXXX)
# output : filtered matrix
#
# author. Ayden Saffari <ayden.saffari@lshtm.ac.uk>
# affiliations. MRC ING, LSHTM
#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/


il450Kfilter <- function (object, cross = F, multi = F, sex = F, snps = c("CpG", "SBE"), maf = 0) 
{
  
  require("IlluminaHumanMethylation450kanno.ilmn12.hg19")
  require("RCurl")
  
  #remove rs probes
  object <- object[!(grepl('rs',rownames(object))),]

  #load 450 annotation 
  snpAnno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  
  #get list of cross reactive probes 
  #ref - Chen et al., Epigenetics 2013
  cross_react <- getURL( 		"https://raw.githubusercontent.com/sirselim/illumina450k_filtering/master/48639-non-specific-probes-Illumina450k.csv")
  cross_react <- read.csv(text = cross_react, head = T, as.is = T)
  cross_react_probes <- as.character(cross_react$TargetID)

  #get list of multimapping probes 
  #ref - Benton et al., Genome Bio 2015
  multi_map <- getURL(   
"https://raw.githubusercontent.com/sirselim/illumina450k_filtering/master/HumanMethylation450_15017482_v.1.1_hg19_bowtie_multimap.txt")
  multi_map <- read.csv(text = multi_map, head = F, as.is = T)
  multi_map_probes <- as.character(multi_map$V1)

  sex_probes <- snpAnno$Name[snpAnno$chr %in% c("chrX","chrY")]
  ######################################################
  #1. remove SNP probes (as per original minfi function)
  ######################################################
  maf_cols <- paste0(snps, "_maf")
  snpDF <- snpAnno[rownames(object),]
  choices <- c("Probe_maf", "CpG_maf", "SBE_maf")
  if (!all(choices %in% colnames(snpDF))) 
    stop("The specificed 'snpAnno' is not supported by this function")
  if (sum(!(maf_cols %in% choices)) > 0) {
    stop(
      "snps vector argument must be a combination of  \"Probe\", \"CpG\" and \"SBE\"")
  }
  if (!is.numeric(maf) || maf < 0 || maf > 1) {
    stop("maf argument must be a numeric value between 0 and 1")
  }
  wh <- Reduce(union, lapply(maf_cols, function(xx) {
    which(snpDF[, xx] >= maf)
  }))
  
  #############################################
  #2. remove cross-reactive/non-specific probes
  #############################################
  if(cross){
    cr <- which(snpDF$Name %in% cross_react_probes)
  }else{
    cr <- NA
  }
  
  ##############################
  #3. remove multi mapped probes 
  ##############################
  if(multi){
    mm <- which(snpDF$Name %in% multi_map_probes)
  }else{
    mm <- NA
  }
  
  #####################################
  #4. remove probes on sex chromosomes
  #####################################
  if(sex){
    sx <- which(snpDF$Name %in% sex_probes)
  }else{
    sx <- NA
  }
  
  fil <- c(wh,cr,mm,sx)
  fil <- fil[!(is.na(fil))]
  fil <- sort(fil)
  fil <- unique(fil)
  
  snpDF <- snpDF[-fil, ]
  object[which(rownames(object) %in% rownames(snpDF)),]
}
