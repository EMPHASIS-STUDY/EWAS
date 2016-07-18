#' Dye bias correction
#'
#' Adjusts dye bias by scaling each color channel to have the same mean intensity.
#'
#' Infinium HumanMethylation450 BeadChip
#' @param rg Cy5/Cy3 signal generated by \code{\link{read.rg}()}.
#' @param probes Output from \code{\link{meffil.probe.info}()} compatible with \code{rg}.
#' @param intensity Intensity of both color channels after correct (Default: 5000).
#' @return Cy5/Cy3 signals with average intensity equal to \code{intensity}.
dye.bias.correct <- function(rg, probes, intensity=5000, verbose=F) {
    msg(verbose=verbose)
    stopifnot(is.rg(rg))
    stopifnot(intensity > 100)

    rg$R[,"Mean"] <- rg$R[,"Mean"] * intensity/calculate.intensity.R(rg, probes)
    rg$G[,"Mean"] <- rg$G[,"Mean"] * intensity/calculate.intensity.G(rg, probes)
    rg
}

calculate.intensity.R <- function(rg, probes) {
    addresses <- probes$address[which(probes$target %in% c("NORM_A", "NORM_T")
                                      & probes$dye == "R")]
    idx <- match(addresses, rownames(rg$R))
    idx <- na.omit(idx)
    if (length(idx) < 10) ## there are 93 in total
        stop("Seems like this IDAT file does not match the supplied chip annotation")

    mean(rg$R[idx,"Mean"], na.rm=T)
}

calculate.intensity.G <- function(rg, probes) {
    addresses <- probes$address[which(probes$target %in% c("NORM_G", "NORM_C")
                                      & probes$dye == "G")]
    idx <- match(addresses, rownames(rg$G))
    idx <- na.omit(idx)
    if (length(idx) < 10) ## there are 93 in total
        stop("Seems like this IDAT file does not match the supplied chip annotation")
        
    mean(rg$G[idx,"Mean"], na.rm=T)
}
