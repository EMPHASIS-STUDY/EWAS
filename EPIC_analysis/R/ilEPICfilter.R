#TODO - filter SNP in probe (based on Pidsley paper)
#TODO - filter SNP at single base extension (based on Pidsley paper)
#TODO - test

#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
# ilEPICfilter <- R function to filter EPIC methylation data in order to remove SNP,
# cross-reactive/non-specific, and sex chromosome probes
#
# based on the SNPfilter function from minfi
#
# input : matrix of beta or m values, with rownames = probe identfiers (cgXXXXXXXX)
# output : filtered matrix
#
# author. Ayden Saffari <ayden.saffari@lshtm.ac.uk>
# affiliations. MRC ING, LSHTM
#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/

ilEPICfilter <- function (object, cross = F, sex = F, snpfil = F, snps = c("CpG", "SBE"), maf = 0) 
{
  
  require("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
  require("RCurl")
  
  #remove rs genotyping probes
  object <- object[!(grepl('rs',rownames(object))),]

  #load EPIC annotation 
  snpAnno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
  
  #get list of cross reactive probes 
  #ref - Pidsley et al., Genome Biology 2016
  cross_react <- getURL("https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5055731/bin/13059_2016_1066_MOESM1_ESM.csv")
  cross_react <- read.csv(text = cross_react, head = T, as.is = T)
  cross_react_probes <- as.character(cross_react[,1])

  sex_probes <- snpAnno$Name[snpAnno$chr %in% c("chrX","chrY")]
  
  ######################################################
  #1. remove SNP probes (as per original minfi SNPfilter function)
  ######################################################
  maf_cols <- paste0(snps, "_maf")
  snpDF <- snpAnno[rownames(object),]
  choices <- c("Probe_maf", "CpG_maf", "SBE_maf")
  if (!all(choices %in% colnames(snpDF))) 
    stop("The specificed 'snpAnno' is not supported by this function")
  if (sum(!(maf_cols %in% choices)) > 0) {
    stop("snps vector argument must be a combination of  \"Probe\", \"CpG\" and \"SBE\"")
  }
  if (!is.numeric(maf) || maf < 0 || maf > 1) {
    stop("maf argument must be a numeric value between 0 and 1")
  }
  
  if(snpfil){
  wh <- Reduce(union, lapply(maf_cols, function(xx) {
    which(snpDF[, xx] >= maf)
  }))
  }else{
    wh <- NA
  }
  
  #############################################
  #2. remove cross-reactive/non-specific probes
  #############################################
  if(cross){
    cr <- which(snpDF$Name %in% cross_react_probes)
  }else{
    cr <- NA
  }
    
  #####################################
  #3. remove probes on sex chromosomes
  #####################################
  if(sex){
    sx <- which(snpDF$Name %in% sex_probes)
  }else{
    sx <- NA
  }
  
  fil <- c(wh,cr,sx)
  fil <- fil[!(is.na(fil))]
  fil <- sort(fil)
  fil <- unique(fil)
  
  if(length(fil) > 0){
    snpDF <- snpDF[-fil, ]
  }
  object[which(rownames(object) %in% rownames(snpDF)),]
}
