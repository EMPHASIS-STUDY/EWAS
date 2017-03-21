#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
# R batch script to perform normalisation on EPIC methylation data
# using the meffil package and following the goDMC protocol
#
# arg 1: path to data directory 
# arg 2: number of PCs to adjust for
# arg 3: list of batch variables (as single row .tsv file) (assumed to be in same dir)
# arg 4: number of cores to use
#
# e.g. Rscript 02_meffil_norm.R /home/meffil_batch_test 14 batch.tsv 8
#
# output : meffil normilized object, normilzation report, matrix of normalized beta 
# values, list of bad CPGs, matrix of estimated cell counts
#
# author. Ayden Saffari <ayden.saffari@lshtm.ac.uk>
# affiliations. MRC ING, LSHTM
#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/

library("meffil")
library("ggplot2")

#parse arguments
args <- commandArgs(trailingOnly=TRUE)

#load data
setwd(args[1])
in_file <- paste0(args[1],"/R_workspaces/meffil_QC.RData")
load(in_file)

#number of pcs
pcs <- as.integer(args[2])

#batch variables
batch_var <- if(!(is.na(args[3]))) scan(file=args[3],what="character") else 
					c("Sex","Slide","Array")

#reserve multiprocessor cores
ncores <- ifelse(!(is.na(args[4])),args[4],8)
options(mc.cores=ncores)

norm_obj <- meffil.normalize.quantiles(qc_obj, number.pcs=pcs)
out_file <- paste0(robj_path,"norm_obj.RData")
message(paste0("saving to: ",out_file))
save(norm_obj,file=out_file)
norm_beta <- meffil.normalize.samples(norm_obj, cpglist.remove=qc_summ$bad.cpgs$name)
message(paste(paste0("beta matrix dimensions: ",dim(norm_beta)[1]),dim(norm_beta)[2]))
out_file <- paste0(robj_path,"norm_beta.RData")
message(paste0("saving to: ",out_file))
save(norm_beta,file=out_file)


out_file <- paste0(robj_path,"bad_cpgs.RData")
badcpgs <- qc_summ$bad.cpgs$name
message(paste0("number of cpgs failing QC: ",length(badcpgs)))
message(paste0("saving to: ",out_file))
save(badcpgs,file=out_file)


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# generate normalization report
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
norm_param <- meffil.normalization.parameters(
    norm_obj,
    variables=batch_var,
    control.pcs=1:pcs,
    batch.pcs=1:pcs,
    batch.threshold=0.01
)

meff_pcs <- meffil.methylation.pcs(norm_beta,probe.range=20000)
out_file <- paste0(robj_path,"meff_pcs.RData")
message(paste0("saving to: ",out_file))
save(meff_pcs,file=out_file)
norm_summ <- meffil.normalization.summary(norm_obj, pcs=meff_pcs,parameters=norm_param)
out_file <- paste0(robj_path,"norm_summ.RData")
message(paste0("saving to: ",out_file))
save(norm_summ,file=out_file)
out_file <- paste0(res_path,"normalization-report.html")
message(paste0("saving to: ",out_file))
meffil.normalization.report(norm_summ, output.file=out_file)

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# export cell counts
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
cell_counts <- t(sapply(qc_obj, function(x) x$cell.counts$counts))
cell_counts <- data.frame(IID=row.names(cell_counts),cell_counts)
out_file <- paste0(res_path,"cell_counts.txt")
message(paste0("saving to: ",out_file))
write.table(cell_counts,out_file,sep="\t",row.names=F,col.names=T,quote=F)

out_file <- paste0(robj_path,"cell_counts.RData")
message(paste0("saving to: ",out_file))
save(cell_counts,file=out_file)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

out_file <- paste0("meffil","_QC.RData")
out_file <- paste0(rwork_path,out_file)
message(paste0("saving to: ",out_file))
save.image(out_file)