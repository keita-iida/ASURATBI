#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
#' Add metadata of variables and samples.
#'
#' This function adds metadata of variables and samples.
#'
#' @param sce A SingleCellExperiment object.
#' @param mitochondria_symbol A string representing for mitochondrial genes.
#'   This function computes percents of reads that map to the mitochondrial
#'   genes. Examples are `^MT-`, `^mt-`, etc.
#'
#' @return A SingleCellExperiment object.
#' @import SingleCellExperiment
#' @import SummarizedExperiment
#' @export
#'
#' @examples
#' data(pbmc_eg)
#' pbmc <- add_metadata(sce = pbmc_eg, mitochondria_symbol = "^MT-")
#'
add_metadata <- function(sce = NULL, mitochondria_symbol = NULL){
  mat <- assay(sce, "counts")
  #--------------------------------------------------
  # Variable metadata
  #--------------------------------------------------
  rowData(sce)$nSamples <- apply(mat, 1, function(x) sum(x > 0))
  #--------------------------------------------------
  # Sample metadata
  #--------------------------------------------------
  sce$nReads <- as.integer(apply(mat, 2, function(x) sum(x)))
  sce$nGenes <- as.integer(apply(mat, 2, function(x) sum(x > 0)))
  sce$percMT <- apply(mat[grepl(mitochondria_symbol, rownames(mat)), ],
                      2, function(x) 100 * sum(x)) /
                  apply(mat, 2, function(x) sum(x))

  return(sce)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
#' Remove variables based on expression profiles across samples.
#'
#' This function removes low expressed variable data.
#'
#' @param sce A SingleCellExperiment object.
#' @param min_nsamples An integer. This function removes variables for which
#'   the numbers of non-zero expressing samples are less than this value.
#'
#' @return A SingleCellExperiment object.
#' @import SingleCellExperiment
#' @import SummarizedExperiment
#' @export
#'
#' @examples
#' data(pbmc_eg)
#' pbmc <- add_metadata(sce = pbmc_eg, mitochondria_symbol = "^MT-")
#' pbmc <- remove_variables(sce = pbmc, min_nsamples = 10)
#'
remove_variables <- function(sce = NULL, min_nsamples = 0){
  mat <- assay(sce, "counts")
  inds <- which(apply(mat, 1, function(x) sum(x > 0)) >= min_nsamples)

  return(sce[inds, ])
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
#' Remove samples based on expression profiles across variables.
#'
#' This function removes sample data by setting minimum and maximum threshold
#'   values for the metadata.
#'
#' @param sce A SingleCellExperiment object.
#' @param min_nReads A minimum threshold value of the number of reads.
#' @param max_nReads A maximum threshold value of the number of reads.
#' @param min_nGenes A minimum threshold value of the number of non-zero
#'   expressed genes.
#' @param max_nGenes A maximum threshold value of the number of non-zero
#'   expressed genes.
#' @param min_percMT A minimum threshold value of the percent of reads that map
#'   to mitochondrial genes.
#' @param max_percMT A maximum threshold value of the percent of reads that map
#'   to mitochondrial genes.
#'
#' @return A SingleCellExperiment object.
#' @import SingleCellExperiment
#' @import SummarizedExperiment
#' @export
#'
#' @examples
#' data(pbmc_eg)
#' pbmc <- add_metadata(sce = pbmc_eg, mitochondria_symbol = "^MT-")
#' pbmc <- remove_samples(sce = pbmc, min_nReads = 0, max_nReads = 1e+10,
#'                        min_nGenes = 0, max_nGenes = 1e+10,
#'                        min_percMT = NULL, max_percMT = NULL)
#'
remove_samples <- function(
  sce = NULL, min_nReads = NULL, max_nReads = NULL,
  min_nGenes = NULL, max_nGenes = NULL, min_percMT = NULL, max_percMT = NULL
){
  mat <- assay(sce, "counts")
  if(!((is.null(min_nReads)) || (is.null(max_nReads)))){
    inds_1 <- which((sce$nReads >= min_nReads) & (sce$nReads <= max_nReads))
  }else{
    inds_1 <- seq_len(dim(sce)[2])
  }
  if(!((is.null(min_nGenes)) || (is.null(max_nGenes)))){
    inds_2 <- which((sce$nGenes >= min_nGenes) & (sce$nGenes <= max_nGenes))
  }else{
    inds_2 <- seq_len(dim(sce)[2])
  }
  if(!((is.null(min_percMT)) || (is.null(max_percMT)))){
    inds_3 <- which((sce$percMT >= min_percMT) & (sce$percMT <= max_percMT))
  }else{
    inds_3 <- seq_len(dim(sce)[2])
  }
  inds <- intersect(intersect(inds_1, inds_2), inds_3)

  return(sce[, inds])
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
#' Remove variables based on the mean expression levels across samples.
#'
#' This function removes variable data such that the mean expression levels
#'   across samples are less than `min_meannReads`.
#'
#' @param sce A SingleCellExperiment object.
#' @param min_meannReads An integer. This function removes variables for which
#'   the mean read counts are less than this value.
#'
#' @return A SingleCellExperiment object.
#' @import SingleCellExperiment
#' @import SummarizedExperiment
#' @export
#'
#' @examples
#' data(pbmc_eg)
#' pbmc <- remove_variables_second(sce = pbmc_eg, min_meannReads = 0.01)
#'
remove_variables_second <- function(sce = NULL, min_meannReads = 0){
  mat <- assay(sce, "counts")
  inds <- which(apply(mat, 1, mean) >= min_meannReads)

  return(sce[inds, ])
}
