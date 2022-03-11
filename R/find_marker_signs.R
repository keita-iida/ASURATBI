#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
#' Compute separation indices of sign scores for given two clusters.
#'
#' This function computes separation indices of sign scores for given two
#'   clusters.
#'
#' @param sce A SingleCellExperiment object.
#' @param labels A vector of labels of all the samples.
#' @param nrand_samples An integer for the number of samples used for
#'   random sampling, which samples at least one sample per cluster.
#' @param ident_1 Label names identifying cluster numbers,
#'   e.g., ident_1 = 1, ident_1 = c(1, 3).
#' @param ident_2 Label names identifying cluster numbers,
#'   e.g., ident_2 = 2, ident_2 = c(2, 4).
#'
#' @return A SingleCellExperiment object.
#' @import SingleCellExperiment
#' @import SummarizedExperiment
#' @import S4Vectors
#' @export
#'
#' @examples
#' data(pbmcs_eg)
#' labels <- SummarizedExperiment::colData(pbmcs_eg$GO)$seurat_clusters
#' pbmcs_eg$GO <- compute_sepI_clusters(sce = pbmcs_eg$GO, labels = labels,
#'                                      nrand_samples = 10, ident_1 = 1,
#'                                      ident_2 = c(0, 2))
#' # The results are stored in `metadata(pbmcs_eg$GO)$marker_signs`.
#'
compute_sepI_clusters <- function(
  sce = NULL, labels = NULL, nrand_samples = NULL,
  ident_1 = NULL, ident_2 = NULL
){
  #--------------------------------------------------
  # Error handling
  #--------------------------------------------------
  if((length(unique(sort(labels))) <= 1) || (is.null(labels))){
    stop("Insufficient labels for calculating separation index.")
  }
  if(length(intersect(ident_1, ident_2)) != 0){
    stop("ident_1 and ident_2 must be disjoint.")
  }
  if((!is.element(ident_1, labels)) || (!is.element(ident_2, labels))){
    stop("labels must include both ident_1 and ident_2.")
  }
  if(dim(sce)[2] != length(labels)){
    stop("dim(sce)[2] must be equal to length(labels).")
  }
  if(!is.null(nrand_samples)){
    tmp <- sce[, colnames(sce)[which(labels %in% union(ident_1, ident_2))]]
    if(dim(tmp)[2] < nrand_samples){
      stop("nrand_samples must be<= the number of identified samples.")
    }
  }
  #--------------------------------------------------
  # Preparation
  #--------------------------------------------------
  idents <- sort(union(ident_1, ident_2))
  labels <- as.character(labels)
  #--------------------------------------------------
  # Random sampling
  #--------------------------------------------------
  if(!(is.null(nrand_samples))){
    index_list <- list()
    baseindex_list <- list()
    nonbaseindex_list <- list()
    for(i in seq_len(length(idents))){
      index_list[[i]] <- which(labels == idents[i])
      if(length(index_list[[i]]) == 1){
        baseindex_list[[i]] <- index_list[[i]]
      }else{
        baseindex_list[[i]] <- sample(index_list[[i]], 1, prob = NULL)
      }
      nonbaseindex_list[[i]] <- setdiff(index_list[[i]], baseindex_list[[i]])
    }
    nonbaseindices <- unlist(nonbaseindex_list)
    n <- nrand_samples - length(idents)
    if(length(nonbaseindices) >= 2){
      nonbaseindices <- sample(nonbaseindices, n, prob = NULL)
    }
    final_idents <- union(unlist(baseindex_list), nonbaseindices)
    subsce <- sce[, final_idents]
    popu_1 <- colnames(sce)[final_idents[labels[final_idents] %in% ident_1]]
  }else{
    subsce <- sce[, colnames(sce)[which(labels %in% union(ident_1, ident_2))]]
    popu_1 <- colnames(sce)[which(labels %in% ident_1)]
  }
  submat <- as.matrix(assay(subsce, "counts"))
  #--------------------------------------------------
  # Compute separation indices.
  #--------------------------------------------------
  res <- data.frame(
    Ident_1 = paste(ident_1, collapse = "/"),
    Ident_2 = paste(ident_2, collapse = "/"),
    SignID = NA,
    Description = NA,
    CorrGene = NA,
    WeakCorrGene = NA,
    sepI = NA,
    Rank = rep(NA, dim(subsce)[1])
  )
  for(i in seq_len(dim(subsce)[1])){
    res$SignID[i] <- rownames(subsce)[i]
    res$Description[i] <- rowData(subsce)$Description[i]
    res$CorrGene[i] <- rowData(subsce)$CorrGene[i]
    res$WeakCorrGene[i] <- rowData(subsce)$WeakCorrGene[i]
    #--------------------------------------------------
    # Preparation
    #--------------------------------------------------
    data <- submat[which(rownames(submat) == res$SignID[i]), ]
    vec_1 <- sort(data, decreasing = FALSE)
    vec_1 <- ifelse(names(vec_1) %in% popu_1, 1, 0)
    #--------------------------------------------------
    # Count the number of steps in bubble sort.
    #--------------------------------------------------
    dist1 <- bubble_sort(list(vec_1, 0))[[2]]     # dist(vec_1, (0,...,1)).
    dist2 <- bubble_sort(list(1 - vec_1, 0))[[2]] # dist(vec_1, (1,...,0)).
    res$sepI[i] <- round((dist2 - dist1) / (dist2 + dist1), digits = 6)
  }
  #--------------------------------------------------
  # Arrange the data frame in order of res$sepI
  #--------------------------------------------------
  inds <- order(res$sepI, decreasing = TRUE)
  res <- res[inds, ]
  res$Rank <- seq_len(length(inds))
  rownames(res) <- seq_len(nrow(res))
  #--------------------------------------------------
  # Output
  #--------------------------------------------------
  slot_name <- paste("Label_", paste(ident_1, collapse = "/"), "_vs_",
                     paste(ident_2, collapse = "/"), sep = "")
  metadata(sce)$marker_signs[[slot_name]] <- res

  return(sce)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
#' Compute separation indices for each cluster against the others.
#'
#' This function computes separation indices for each cluster versus the others.
#'
#' @param sce A SingleCellExperiment object.
#' @param labels A vector of labels of all the samples (cells).
#' @param nrand_samples An integer for the number of samples used for
#'   random sampling, which samples at least one sample per cluster.
#'
#' @return A SingleCellExperiment object.
#' @import SingleCellExperiment
#' @import SummarizedExperiment
#' @import S4Vectors
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
#'
#' @examples
#' data(pbmcs_eg)
#' labels <- SummarizedExperiment::colData(pbmcs_eg$GO)$seurat_clusters
#' pbmcs_eg$GO <- compute_sepI_all(sce = pbmcs_eg$GO, labels = labels,
#'                                 nrand_samples = 10)
#' # The results are stored in `metadata(pbmcs_eg$GO)$marker_signs`.
#'
compute_sepI_all <- function(sce = NULL, labels = NULL, nrand_samples = NULL){
  #--------------------------------------------------
  # Error handling
  #--------------------------------------------------
  if((length(unique(sort(labels))) <= 1) || (is.null(labels))){
    stop("Insufficient labels for calculating separation index.")
  }
  #--------------------------------------------------
  # Preparation
  #--------------------------------------------------
  res <- list()
  tmp <- sce
  metadata(tmp)$marker_signs <- NULL
  #--------------------------------------------------
  # Loop
  #--------------------------------------------------
  idents <- unique(sort(labels))

  # Initializes the progress bar
  pb <- txtProgressBar(min = 0, max = length(idents), style = 3,
                       width = 50, char = "=")
  for(i in seq_len(length(idents))){
    setTxtProgressBar(pb, i)
    ident_1 <- idents[i]
    ident_2 <- setdiff(idents, ident_1)
    tmp <- compute_sepI_clusters(sce = tmp, labels = labels,
                                 nrand_samples = nrand_samples,
                                 ident_1 = ident_1, ident_2 = ident_2)
    res[[i]] <- metadata(tmp)$marker_signs[[i]]
    slot_name <- paste("Label_", paste(ident_1, collapse = "/"), "_vs_",
                       paste(ident_2, collapse = "/"), sep = "")
    metadata(sce)$marker_signs[[slot_name]] <- res[[i]]
  }
  close(pb)
  res_all <- c()
  for(i in seq_len(length(res))){
    res_all <- rbind(res_all, res[[i]])
  }
  rownames(res_all) <- seq_len(nrow(res_all))
  metadata(sce)$marker_signs$all <- res_all

  return(sce)
}

