#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
# R batch script to perform QC on EPIC methylation data
# using the meffil package and following the goDMC protocol
#
# arg 1: path to data directory 
# arg 2: number of cores to use
#
# e.g. Rscript 01_meffil_QC.R /home/meffil_batch_test 8
#
# output : meffil QC object, QC summary, QC report
#
# author. Ayden Saffari <ayden.saffari@lshtm.ac.uk>
# affiliations. MRC ING, LSHTM
#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# initialization
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#install and/or load required packages 

cran_p <- c("devtools","markdown","knitr","ggplot2")
bioc_p <- c("illuminaio","limma","IlluminaHumanMethylation450kmanifest",
            "IlluminaHumanMethylation450kanno.ilmn12.hg19",
            "IlluminaHumanMethylationEPICanno.ilm10b2.hg19","DNAcopy","meffil")

test_cran <- lapply(cran_p,require,character.only = TRUE)
lapply(cran_p[which(test_cran == FALSE)],function(x){
       install.packages(x,character.only = TRUE,
       repos="http://cran.ma.imperial.ac.uk/")})
source('http://www.bioconductor.org/biocLite.R')
test_bioc <- lapply(bioc_p,require,character.only = TRUE)
lapply(bioc_p[which(test_bioc == FALSE)],
	   function(x){biocLite(x,character.only = TRUE)})

lapply(cran_p,library,character.only = TRUE)
lapply(bioc_p,library,character.only = TRUE)

#parse arguments
args <- commandArgs(trailingOnly=TRUE)

#set up directory structure
res_path <- paste0(args[1],"/results/")
robj_path <- paste0(args[1],"/R_objects/")
rwork_path <- paste0(args[1],"/R_workspaces/")
dir.create(res_path)
dir.create(robj_path)
dir.create(rwork_path)

#reserve multiprocessor cores
ncores <- ifelse(!(is.na(args[2])),args[2],8)
options(mc.cores=ncores)

#read in sample sheet
sample_sheet <- meffil.read.samplesheet(base=args[1],recursive=T,
					ignore.case=T,verbose=T)

#sometimes samplesheet function does not guess filenames correctly, 
#if so, use list of idats from basenames function (seems to be more reliable)
if(any(sample_sheet$Basename == "character(0)")){
    sample_sheet$Basename <- as.character(sample_sheet$Basename)
    idats <- meffil.basenames(args[1], recursive = T)
    idats <- data.frame(barcode=basename(idats),Basename=idats)
    idats$Basename <- as.character(idats$Basename)
    sample_sheet$Basename[sample_sheet$barcode %in% idats$barcode] <- 
    idats[idats$barcode %in% sample_sheet$barcode,2]
}

sample_sheet <- sample_sheet[-which(sample_sheet$Basename == "character(0)"),]	

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# create QC objects, perform background correction, dye bias correction, sex prediction,
# cell count estimates
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

qc_obj <- meffil.qc(sample_sheet, cell.type.reference="blood gse35069 complete",
				    verbose=T)
out_file <- paste0(robj_path,"qc_obj.RData")
message(paste0("saving to: ",out_file))
save(qc_obj,file=out_file)

message(paste0("samples processed: ",length(qc_obj)))
message(c("sample ids:\n",paste(names(qc_obj),"\n")))


qc_param <- meffil.qc.parameters(beadnum.samples.threshold = 0.1,
				 detectionp.samples.threshold = 0.1,
				 detectionp.cpgs.threshold = 0.1,
				 beadnum.cpgs.threshold = 0.1,
				 sex.outlier.sd = 5)
								 
qc_summ <- meffil.qc.summary(qc_obj,parameters = qc_param)
out_file <- paste0(robj_path,"qc_summ.RData")
message(paste0("saving to: ",out_file))
save(qc_summ, file=out_file)

out_file <- paste0(res_path,"qc_report.html")
message(paste0("saving to: ",out_file))
meffil.qc.report(qc_summ, output.file=out_file)

outliers <- qc_summ$bad.samples
message(paste0("total number of sample outliers: ", length(unique(outliers$sample.name))))

idx <- outliers$issue %in% c("Control probe (dye.bias)", "Methylated vs Unmethylated",
                              "X-Y ratio outlier","Low bead numbers",
                              "Detection p-value","Genotype mismatch",
                              "Control probe (bisulfite1)","Control probe (bisulfite2)")
message("outliers to remove: \n")
message(c("Control probe (dye.bias)\n","Methylated vs Unmethylated\n","X-Y ratio outlier\n",
	"Low bead numbers\n","Detection p-value\n","Genotype mismatch\n",
	"Control probe (bisulfite1)\n","Control probe(bisulfite2)\n"))                              
outliers <- outliers[idx,]
message(paste0("removing this many samples: ", length(unique(outliers$sample.name))))
qc_obj <- meffil.remove.samples(qc_obj, outliers$sample.name)
message(paste0("samples remaining: ",length(qc_obj)))

out_file <- paste0(robj_path,"qc_obj_clean.RData")
message(paste0("saving to: ",out_file))
save(qc_obj,file=out_file)

qc_summ <- meffil.qc.summary(qc_obj,parameters = qc_param)
out_file <- paste0(robj_path,"qc_summ_clean.RData")
message(paste0("saving to: ",out_file))
save(qc_summ, file=out_file)
out_file <- paste0(res_path,"qc_report_clean.html")
meffil.qc.report(qc_summ, output.file=out_file)

pc_plot <- meffil.plot.pc.fit(qc_obj)
out_file <- paste0(res_path,"pc_fit.pdf")
message(paste0("saving to: ",out_file))
ggsave(pc_plot$plot,filename=out_file,height=6,width=6)

out_file <- paste0("meffil","_QC.RData")
out_file <- paste0(rwork_path,out_file)
message(paste0("saving to: ",out_file))
save.image(out_file)


