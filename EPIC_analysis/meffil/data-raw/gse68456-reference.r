


#' Defines cell type reference "cord blood gse68456" 
#' for estimating cord blood cell counts.
#' The first is based on
#' six cell types: CD4T, CD8T, Mono, Bcell, NK, Gran.
#' The second is based on 
#' the same cell types but with Gran replaced by Neu and Eos.

retrieve.gse68456 <- function(dir) {
    wd <- getwd()
    on.exit(setwd(wd))
    setwd(dir)
    
    cat("Downloading data ...\n")
    filename <- "gse68456.tar"
    download.file("http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE68456&format=file", filename)
    
    cat("Unzipping data ...\n")
    system(paste("tar xvf", filename))

    filenames <- list.files(path=".", pattern="Red.idat.gz$")
    basenames <- sub("_Red.idat.gz$", "", filenames)
    samples <- data.frame(Basename=basenames,
                          gsm=sub("([^_]+)_.*", "\\1", basenames),
                          participant=sub(".*_([^_]_)_.*", "\\1", basenames),
                          cell.type=sub(".*_.*_(.+)$", "\\1", basenames),
                          stringsAsFactors=F)

    samples$cell.type[which(samples$cell.type == "B")] <- "Bcell"
    samples$cell.type[which(samples$cell.type == "G")] <- "Gran"
    samples$cell.type[which(samples$cell.type == "Mo")] <- "Mono"
    samples$Sex <- NA
    samples$Sample_Name <- samples$gsm

    samples
}

create.gse68456.reference <- function() {
    number.pcs <- 5
    verbose <- T
    
    dir.create(temp.dir <- tempfile(tmpdir="."))
    on.exit(unlink(temp.dir, recursive=TRUE))
    
    samplesheet <- retrieve.gse68456(temp.dir)
    samplesheet$Basename <- file.path(temp.dir, samplesheet$Basename)

    ds <- meffil.normalize.dataset(samplesheet,
                                   just.beta=F,
                                   qc.file="gse68456-qc-report.html",
                                   author="Olivia de Goede, et al.",
                                   study="Purified cord blood cell type methylation (GEO:GSE68456)",
                                   number.pcs=number.pcs,
                                   norm.file="gse68456-normalization-report.html",
                                   chip="450k",
                                   featureset="common",
                                   verbose=verbose)

    samplesheet <- samplesheet[match(colnames(ds$M), samplesheet$Sample_Name),]
    
    cell.types <- c("CD3T", "CD4T","CD8T","Mono","Bcell","NK","Gran")
    selected <- samplesheet$cell.type %in% cell.types
    meffil.add.cell.type.reference("cord blood gse68456",
                                   ds$M[,selected], ds$U[,selected],
                                   cell.types=samplesheet$cell.type[selected],
                                   chip="450k",
                                   featureset="common",
                                   verbose=verbose)
}
