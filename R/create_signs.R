#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
#' Remove signs including too few or too many genes.
#'
#' This function removes signs including too few or too many genes.
#'
#' @param sce A SingleCellExperiment object.
#' @param min_ngenes Minimum number of genes, which must be greater than one.
#' @param max_ngenes Maximum number of genes, which must be greater than one.
#'
#' @return A SingleCellExperiment object.
#' @import SingleCellExperiment
#' @import SummarizedExperiment
#' @import S4Vectors
#' @export
#'
#' @examples
#' data(pbmc_eg)
#' data(human_GO_eg)
#' pbmcs <- list(GO = pbmc_eg)
#' S4Vectors::metadata(pbmcs$GO) <- list(sign = human_GO_eg[["BP"]])
#' pbmcs$GO <- remove_signs(sce = pbmcs$GO, min_ngenes = 2, max_ngenes = 1000)
#' # The results are stored in `metadata(pbmcs$GO)$sign`.
#'
remove_signs <- function(sce = NULL, min_ngenes = 2, max_ngenes = 1000){
  #--------------------------------------------------
  # Error handling
  #--------------------------------------------------
  if(is.null(metadata(sce)$sign)){
    stop("metadata(sce)$sign must have a database.")
  }
  if(min_ngenes <= 1){ stop("min_ngenes should be >= 2.") }
  if(max_ngenes <= 1){ stop("max_ngenes should be >= 2.") }
  #--------------------------------------------------
  # Select innate genes.
  #--------------------------------------------------
  innate_g <- data.frame(gene = rownames(sce), geneID = rowData(sce)$geneID)
  tmp <- metadata(sce)$sign
  for(i in seq_len(nrow(tmp))){
    #------------------------------
    # GeneID: ENTREZ Gene ID
    #------------------------------
    geneIDs <- unlist(strsplit(tmp$GeneID[i], "/"))
    geneIDs <- geneIDs[which(geneIDs %in% innate_g$geneID)]
    tmp$GeneID[i] <- paste(geneIDs, collapse = "/")
    #------------------------------
    # Gene: gene symbol
    #------------------------------
    genes <- c()
    for(j in seq_len(length(geneIDs))){
      genes <- c(genes, innate_g[which(innate_g$geneID == geneIDs[j]), ]$gene)
    }
    tmp$Gene[i] <- paste(genes, collapse = "/")
    #------------------------------
    # Count
    #------------------------------
    tmp$Count[i] <- as.integer(length(geneIDs))
  }
  metadata(sce)$sign <- tmp
  #--------------------------------------------------
  # Remove IDs including too few or too many genes.
  #--------------------------------------------------
  tmp <- metadata(sce)$sign
  res <- tmp[which((tmp$Count >= min_ngenes) & (tmp$Count <= max_ngenes)), ]
  rownames(res) <- seq_len(nrow(res))
  metadata(sce)$sign <- as.data.frame(res)

  return(sce)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
#' Cluster each functional gene set into three groups.
#'
#' This function clusters each functional gene set into strongly, variably, and
#'   weakly correlated gene sets.
#'
#' @param sce A SingleCellExperiment object.
#' @param cormat A correlation matrix of gene expressions.
#' @param th_posi A threshold of positive correlation coefficient.
#' @param th_nega A threshold of negative correlation coefficient.
#'
#' @return A SingleCellExperiment object.
#' @import SingleCellExperiment
#' @import SummarizedExperiment
#' @import S4Vectors
#' @importFrom cluster pam
#' @export
#'
#' @examples
#' data(pbmc_eg)
#' data(human_GO_eg)
#' mat <- t(as.matrix(SummarizedExperiment::assay(pbmc_eg, "centered")))
#' pbmc_cormat <- cor(mat, method = "spearman")
#' pbmcs <- list(GO = pbmc_eg)
#' S4Vectors::metadata(pbmcs$GO) <- list(sign = human_GO_eg[["BP"]])
#' pbmcs$GO <- remove_signs(sce = pbmcs$GO, min_ngenes = 2, max_ngenes = 1000)
#' pbmcs$GO <- cluster_genesets(sce = pbmcs$GO, cormat = pbmc_cormat,
#'                              th_posi = 0.24, th_nega = -0.20)
#' # The results are stored in `metadata(pbmcs$GO)$sign`.
#'
cluster_genesets <- function(
  sce = NULL, cormat = NULL, th_posi = NULL, th_nega = NULL
){
  #--------------------------------------------------
  # Error handling
  #--------------------------------------------------
  if(th_posi < 0){ stop("th_posi should be >= 0.") }
  if(th_nega > 0){ stop("th_nega should be <= 0.") }
  #--------------------------------------------------
  # Correlation graph-based decomposition of FGSs
  #--------------------------------------------------
  df <- metadata(sce)$sign
  res <- c("SignID", "Description", "IC",
           "CountStrgCorrGene", "CountVariCorrGene", "CountWeakCorrGene",
           "CorrStrg", "CorrVari", "CorrWeak",
           "StrgCorrGene", "VariCorrGene", "WeakCorrGene",
           "StrgCorrGeneID", "VariCorrGeneID", "WeakCorrGeneID")
  res <- data.frame(matrix(ncol = 15, nrow = 0, dimnames = list(NULL, res)))
  for(i in seq_len(nrow(df))){
    genes <- unlist(strsplit(df$Gene[i], "/"))
    if(length(genes) <= 1){
      next
    }
    inds <- which(rownames(cormat) %in% genes)
    if(length(inds) <= 2){
      next
    }
    mat <- cormat[inds, inds]
    #--------------------------------------------------
    # Select genes having strong correlations with others.
    #--------------------------------------------------
    diag(mat) <- -99
    inds_posi <- which(apply(mat, 2, function(x) sum(x >= th_posi)) > 0)
    diag(mat) <- 99
    inds_nega <- which(apply(mat, 2, function(x) sum(x <= th_nega)) > 0)
      
    if(length(inds_posi) == 0){
      #--------------------------------------------------
      # Define strongly, variably, and weakly correlated gene sets.
      #--------------------------------------------------
      genes_strg <- NA ; genes_vari <- NA ; genes_weak <- genes
      mean_strg <- NA ; mean_vari <- NA
      mean_weak <- ifelse(length(genes_weak) <= 1, NA, mean(mat[mat <= 1]))
    }else{
      if(length(inds_nega) == 0){
        #--------------------------------------------------
        # Define strongly, variably, and weakly correlated gene sets.
        #--------------------------------------------------
        genes_vari <- NA ; mean_vari <- NA
          
        submat <- mat[inds_posi, inds_posi]
        genes_strg <- rownames(submat)
        mean_strg <- ifelse(length(genes_strg) <= 1, NA,
                            mean(submat[submat <= 1]))
        if(mean_strg < th_posi){
          genes_strg <- NA ; mean_strg <- NA
          genes_weak <- genes
          mean_weak <- ifelse(length(genes_weak) <= 1, NA, mean(mat[mat <= 1]))
        }else{
          genes_weak <- setdiff(genes, genes_strg)  
          inds <- which(rownames(mat) %in% genes_weak)
          tmp <- mat[inds, inds]
          mean_weak <- ifelse(length(inds) <= 1, NA, mean(tmp[tmp <= 1]))
        }
      }else{
        inds <- union(inds_posi, inds_nega)
        submat <- mat[inds, inds] ; diag(submat) <- 1
        #--------------------------------------------------
        # PAM clustering using cluster package
        #--------------------------------------------------
        c <- pam(submat, k = 2)
        diag(submat) <- 99
        gset <- list() ; mu <- c()
        for(j in seq_len(2)){
          gset[[j]] <- names(c[["clustering"]])[which(c[["clustering"]] == j)]
          inds <- which(rownames(submat) %in% gset[[j]])
          tmp <- submat[inds, inds]
          mu[j] <- ifelse(length(inds) <= 1, NA, mean(tmp[tmp <= 1]))
        }
        #--------------------------------------------------
        # Define strongly, variably, and weakly correlated gene sets.
        #--------------------------------------------------
        inds <- order(mu, decreasing = TRUE)
        genes_strg <- gset[[inds[1]]]
        mean_strg <- mu[inds[1]]
        if(mean_strg < th_posi){
          genes_strg <- NA ; mean_strg <- NA
          genes_vari <- NA ; mean_vari <- NA
          genes_weak <- genes
          mean_weak <- ifelse(length(genes_weak) <= 1, NA, mean(mat[mat <= 1]))
        }else{
          genes_vari <- gset[[inds[2]]]
          mean_vari <- mu[inds[2]]
          genes_weak <- setdiff(genes, union(gset[[1]], gset[[2]]))  
          inds <- which(rownames(mat) %in% genes_weak)
          tmp <- mat[inds, inds]
          mean_weak <- ifelse(length(inds) <= 1, NA, mean(tmp[tmp <= 1]))
        }
      }
    }
    #--------------------------------------------------
    # Result
    #--------------------------------------------------
    tmp <- as.data.frame(rowData(sce))
    if(is.element(NA, genes_strg)){
      genes_strg <- NA ; geneIDs_strg <- NA ; count_strg <- 0
    }else{
      count_strg <- length(genes_strg)
      gs <- c()
      for(j in seq_len(length(genes_strg))){
        inds <- which(rownames(tmp) == genes_strg[j])
        if(ncol(tmp) == 1){
          if(colnames(tmp) == "geneID"){
            gs <- c(gs, tmp[inds, ])
          }else{
            stop("rowData() has only one column; the name must be \"geneID\".")
          }
        }else{
          gs <- c(gs, tmp[inds, ]$geneID)
        }
      }
      geneIDs_strg <- paste(gs, collapse = "/")
      genes_strg <- paste(genes_strg, collapse = "/")
    }
    if(is.element(NA, genes_vari)){
      genes_vari <- NA ; geneIDs_vari <- NA ; count_vari <- 0
    }else{
      count_vari <- length(genes_vari)
      gs <- c()
      for(j in seq_len(length(genes_vari))){
        inds <- which(rownames(tmp) == genes_vari[j])
        if(ncol(tmp) == 1){
          if(colnames(tmp) == "geneID"){
            gs <- c(gs, tmp[inds, ])
          }else{
            stop("rowData() has only one column; the name must be \"geneID\".")
          }
        }else{
          gs <- c(gs, tmp[inds, ]$geneID)
        }
      }
      geneIDs_vari <- paste(gs, collapse = "/")
      genes_vari <- paste(genes_vari, collapse = "/")
    }
    if(is.element(NA, genes_weak)){
      genes_weak <- NA ; geneIDs_weak <- NA ; count_weak <- 0
    }else{
      count_weak <- length(genes_weak)
      gs <- c()
      for(j in seq_len(length(genes_weak))){
        inds <- which(rownames(tmp) == genes_weak[j])
        if(ncol(tmp) == 1){
          if(colnames(tmp) == "geneID"){
            gs <- c(gs, tmp[inds, ])
          }else{
            stop("rowData() has only one column; the name must be \"geneID\".")
          }
        }else{
          gs <- c(gs, tmp[inds, ]$geneID)
        }
      }
      geneIDs_weak <- paste(gs, collapse = "/")
      genes_weak <- paste(genes_weak, collapse = "/")
    }
    dg <- data.frame(SignID = df$ID[i], Description = df$Description[i],
                     IC = df$IC[i],
                     CountStrgCorrGene = count_strg,
                     CountVariCorrGene = count_vari,
                     CountWeakCorrGene = count_weak,
                     CorrStrg = mean_strg,
                     CorrVari = mean_vari,
                     CorrWeak = mean_weak,
                     StrgCorrGene = genes_strg,
                     VariCorrGene = genes_vari,
                     WeakCorrGene = genes_weak,
                     StrgCorrGeneID = geneIDs_strg,
                     VariCorrGeneID = geneIDs_vari,
                     WeakCorrGeneID = geneIDs_weak)
    res <- rbind(res, dg)
  }
  metadata(sce)$sign <- res

  return(sce)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
#' Define signs for strongly and variably correlated gene sets.
#'
#' This function define signs for strongly and variably correlated gene sets.
#'
#' @param sce A SingleCellExperiment object.
#' @param min_cnt_strg An integer for the cutoff value for strongly correlated
#'   gene sets.
#' @param min_cnt_vari An integer for the cutoff value for variably correlated
#'   gene sets.
#'
#' @return A SingleCellExperiment object.
#' @import SingleCellExperiment
#' @import SummarizedExperiment
#' @import S4Vectors
#' @export
#'
#' @examples
#' data(pbmc_eg)
#' data(human_GO_eg)
#' mat <- t(as.matrix(SummarizedExperiment::assay(pbmc_eg, "centered")))
#' pbmc_cormat <- cor(mat, method = "spearman")
#' pbmcs <- list(GO = pbmc_eg)
#' S4Vectors::metadata(pbmcs$GO) <- list(sign = human_GO_eg[["BP"]])
#' pbmcs$GO <- remove_signs(sce = pbmcs$GO, min_ngenes = 2, max_ngenes = 1000)
#' pbmcs$GO <- cluster_genesets(sce = pbmcs$GO, cormat = pbmc_cormat,
#'                              th_posi = 0.24, th_nega = -0.20)
#' pbmcs$GO <- create_signs(sce = pbmcs$GO, min_cnt_strg = 2, min_cnt_vari = 2)
#' # The results are stored in `metadata(pbmcs$GO)$sign_all`.
#'
create_signs <- function(sce = NULL, min_cnt_strg = 2, min_cnt_vari = 2){
  #--------------------------------------------------
  # Error handling
  #--------------------------------------------------
  if(min_cnt_strg <= 0){ stop("min_cnt_strg should be >= 1.") }
  if(min_cnt_vari <= 0){ stop("min_cnt_vari should be >= 1.") }

  df <- metadata(sce)$sign
  #--------------------------------------------------
  # Strongly correlated gene sets
  #--------------------------------------------------
  df_strg <- df[which(df$CountStrgCorrGene >= min_cnt_strg), ]
  if(nrow(df_strg) >= 1){
    df_strg <- data.frame(
      SignID = df_strg$SignID, Description = df_strg$Description,
      IC = df_strg$IC,
      CountStrgCorrGene = df_strg$CountStrgCorrGene,
      CountWeakCorrGene = df_strg$CountWeakCorrGene,
      CorrStrg = df_strg$CorrStrg, CorrWeak = df_strg$CorrWeak,
      StrgCorrGene = df_strg$StrgCorrGene, WeakCorrGene = df_strg$WeakCorrGene,
      StrgCorrGeneID = df_strg$StrgCorrGeneID,
      WeakCorrGeneID = df_strg$WeakCorrGeneID)
  }
  #--------------------------------------------------
  # Variably correlated gene sets
  #--------------------------------------------------
  df_vari <- df[which(df$CountVariCorrGene >= min_cnt_vari), ]
  if(nrow(df_vari) >= 1){
    df_vari <- data.frame(
      SignID = df_vari$SignID, Description = df_vari$Description,
      IC = df_vari$IC,
      CountVariCorrGene = df_vari$CountVariCorrGene,
      CountWeakCorrGene = df_vari$CountWeakCorrGene,
      CorrVari = df_vari$CorrVari, CorrWeak = df_vari$CorrWeak,
      VariCorrGene = df_vari$VariCorrGene, WeakCorrGene = df_vari$WeakCorrGene,
      VariCorrGeneID = df_vari$VariCorrGeneID,
      WeakCorrGeneID = df_vari$WeakCorrGeneID)
  }
  #--------------------------------------------------
  # Store the results.
  #--------------------------------------------------
  metadata(sce)$sign_SCG <- df_strg
  metadata(sce)$sign_VCG <- df_vari
  #--------------------------------------------------
  # Add metadata(sce)$sign_all.
  #--------------------------------------------------
  metadata(sce)$sign <- NULL
  tmp <- c("SignID", "Description", "IC",
           "CountCorrGene", "CountWeakCorrGene", "Corr", "CorrWeak",
           "CorrGene", "WeakCorrGene", "CorrGeneID", "WeakCorrGeneID")
  colnames(df_strg) <- tmp
  colnames(df_vari) <- tmp
  if(nrow(df_strg) == 0){
    stop("SCG is not defined.")
  }else{
    if(nrow(df_vari) == 0){
      metadata(sce)$sign_all <- cbind(CorrType = "SCG", df_strg)
    }else{
      metadata(sce)$sign_all <- rbind(cbind(CorrType = "SCG", df_strg),
                                      cbind(CorrType = "VCG", df_vari))
    }
  }

  return(sce)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
#' Remove redundant signs using semantic similarity matrices.
#'
#' This function removes redundant signs using semantic similarity matrices.
#'
#' @param sce A SingleCellExperiment object.
#' @param similarity_matrix A semantic similarity matrix.
#' @param threshold A threshold value of semantic similarity, used for
#'   regarding biological terms as similar ones 
#' @param keep_rareID If TRUE, biological terms with the larger ICs are kept.
#'
#' @return A SingleCellExperiment object.
#' @import SingleCellExperiment
#' @import SummarizedExperiment
#' @import S4Vectors
#' @export
#'
#' @examples
#' data(pbmc_eg)
#' data(human_GO_eg)
#' mat <- t(as.matrix(SummarizedExperiment::assay(pbmc_eg, "centered")))
#' pbmc_cormat <- cor(mat, method = "spearman")
#' pbmcs <- list(GO = pbmc_eg)
#' S4Vectors::metadata(pbmcs$GO) <- list(sign = human_GO_eg[["BP"]])
#' pbmcs$GO <- remove_signs(sce = pbmcs$GO, min_ngenes = 2, max_ngenes = 1000)
#' pbmcs$GO <- cluster_genesets(sce = pbmcs$GO, cormat = pbmc_cormat,
#'                              th_posi = 0.24, th_nega = -0.20)
#' pbmcs$GO <- create_signs(sce = pbmcs$GO, min_cnt_strg = 2, min_cnt_vari = 2)
#' pbmcs$GO <- remove_signs_redundant(
#'   sce = pbmcs$GO, similarity_matrix = human_GO_eg$similarity_matrix$BP,
#'   threshold = 0.80, keep_rareID = TRUE)
#' # The results are stored in `metadata(pbmcs$GO)$sign_SCG`,
#' # `metadata(pbmcs$GO)$sign_VCG`, `metadata(pbmcs$GO)$sign_all`,
#' # and if there exist, `metadata(pbmcs$GO)$sign_SCG_redundant` and
#' # `metadata(pbmcs$GO)$sign_VCG_redundant`.
#'
remove_signs_redundant <- function(
  sce = NULL, similarity_matrix = NULL, threshold = NULL, keep_rareID = NULL
){
  #--------------------------------------------------
  # Error handling
  #--------------------------------------------------
  if(is.null(similarity_matrix)){
    stop("Semantic similarity is not defined. Skip this process.")
  }

  report <- list()
  s_or_v <- c("sign_SCG", "sign_VCG")
  for(k in seq_len(length(s_or_v))){
    #--------------------------------------------------
    # Report format
    #--------------------------------------------------
    df <- metadata(sce)[[s_or_v[k]]]
    if(nrow(df) == 0){
      next
    }
    inds <- which(rownames(similarity_matrix) %in% df$SignID)
    if(length(inds) < 2){
      next
    }
    submat <- similarity_matrix[inds, inds]
    inds <- order(df$IC, decreasing = FALSE)
    submat <- submat[inds, inds]
    report[[k]] <- c("Similarity", "Kept_SignID", "Removed_SignID",
                     "Kept_Description", "Removed_Description",
                     "Kept_GeneCount", "Removed_GeneCount",
                     "Kept_Gene", "Removed_Gene", "Kept_IC", "Removed_IC")
    report[[k]] <- data.frame(matrix(ncol = 11, nrow = 0,
                                     dimnames = list(NULL, report[[k]])))
    flag <- rep(0, dim(submat)[1])
    if(keep_rareID == TRUE){
      #--------------------------------------------------
      # case I
      #--------------------------------------------------
      I <- dim(submat)[1]
      for(i in seq(I, 2)){
        J <- dim(submat)[2] - (I - i + 1) # Be careful
        if(flag[i] == 1){                 # Be careful
          next
        }
        for(j in seq_len(J)){
          if(is.na(submat[i, j])){
            next
          }
          if(abs(submat[i, j]) >= threshold){
            ind_kept <- which(df$SignID == rownames(submat)[i])
            ind_remo <- which(df$SignID == rownames(submat)[j])
            if(k == 1){
              kept_count <- df[ind_kept, ]$CountStrgCorrGene +
                            df[ind_kept, ]$CountWeakCorrGene
              remo_count <- df[ind_remo, ]$CountStrgCorrGene +
                            df[ind_remo, ]$CountWeakCorrGene
              kept_genes <- paste(df[ind_kept, ]$StrgCorrGene,
                                  df[ind_kept, ]$WeakCorrGene, sep = "/")
              remo_genes <- paste(df[ind_remo, ]$StrgCorrGene,
                                  df[ind_remo, ]$WeakCorrGene, sep = "/")
            }else if(k == 2){
              kept_count <- df[ind_kept, ]$CountVariCorrGene +
                            df[ind_kept, ]$CountWeakCorrGene
              remo_count <- df[ind_remo, ]$CountVariCorrGene +
                            df[ind_remo, ]$CountWeakCorrGene
              kept_genes <- paste(df[ind_kept, ]$VariCorrGene,
                                  df[ind_kept, ]$WeakCorrGene, sep = "/")
              remo_genes <- paste(df[ind_remo, ]$VariCorrGene,
                                  df[ind_remo, ]$WeakCorrGene, sep = "/")
            }
            report[[k]] <- rbind(report[[k]], data.frame(
              Similarity = submat[i, j],
              Kept_SignID = df[ind_kept, ]$SignID,
              Removed_SignID = df[ind_remo, ]$SignID,
              Kept_Description = df[ind_kept, ]$Description,
              Removed_Description = df[ind_remo, ]$Description,
              Kept_GeneCount = kept_count,
              Removed_GeneCount = remo_count,
              Kept_Gene = kept_genes,
              Removed_Gene = remo_genes,
              Kept_IC = df[ind_kept, ]$IC,
              Removed_IC = df[ind_remo, ]$IC
            ))
            flag[j] <- 1
          }
        }
      }
    }else{
      #--------------------------------------------------
      # case II
      #--------------------------------------------------
      I <- dim(submat)[1] - 1
      J <- dim(submat)[2]
      for(i in seq_len(I)){
        j0 <- i + 1       # Be careful
        if(flag[i] == 1){ # Be careful
          next
        }
        for(j in seq(j0, J)){
          if(is.na(submat[i, j])){
            next
          }
          if(abs(submat[i, j]) >= threshold){
            ind_kept <- which(df$SignID == rownames(submat)[i])
            ind_remo <- which(df$SignID == rownames(submat)[j])
            if(k == 1){
              kept_count <- df[ind_kept, ]$CountStrgCorrGene +
                            df[ind_kept, ]$CountWeakCorrGene
              remo_count <- df[ind_remo, ]$CountStrgCorrGene +
                            df[ind_remo, ]$CountWeakCorrGene
              kept_genes <- paste(df[ind_kept, ]$StrgCorrGene,
                                  df[ind_kept, ]$WeakCorrGene, sep = "/")
              remo_genes <- paste(df[ind_remo, ]$StrgCorrGene,
                                  df[ind_remo, ]$WeakCorrGene, sep = "/")
            }else if(k == 2){
              kept_count <- df[ind_kept, ]$CountVariCorrGene +
                            df[ind_kept, ]$CountWeakCorrGene
              remo_count <- df[ind_remo, ]$CountVariCorrGene +
                            df[ind_remo, ]$CountWeakCorrGene
              kept_genes <- paste(df[ind_kept, ]$VariCorrGene,
                                  df[ind_kept, ]$WeakCorrGene, sep = "/")
              remo_genes <- paste(df[ind_remo, ]$VariCorrGene,
                                  df[ind_remo, ]$WeakCorrGene, sep = "/")
            }
            report[[k]] <- rbind(report[[k]], data.frame(
              Similarity = submat[i, j],
              Kept_SignID = df[ind_kept, ]$SignID,
              Removed_SignID = df[ind_remo, ]$SignID,
              Kept_Description = df[ind_kept, ]$Description,
              Removed_Description = df[ind_remo, ]$Description,
              Kept_GeneCount = kept_count,
              Removed_GeneCount = remo_count,
              Kept_Gene = kept_genes,
              Removed_Gene = remo_genes,
              Kept_IC = df[ind_kept, ]$IC,
              Removed_IC = df[ind_remo, ]$IC
            ))
            flag[j] <- 1
          }
        }
      }
    }
    #--------------------------------------------------
    # Remove redundant signs.
    #--------------------------------------------------
    if(length(report[[k]]$Removed_SignID) >= 1){
      kept_IDs <- setdiff(df$SignID, report[[k]]$Removed_SignID)
      res <- df[which(df$SignID %in% kept_IDs), ]
      rownames(res) <- seq_len(nrow(res))
    }else{
      res <- df
    }
    metadata(sce)[[s_or_v[k]]] <- res
  }
  #--------------------------------------------------
  # Modify metadata(sce)$sign_all.
  #--------------------------------------------------
  df_strg <- as.data.frame(metadata(sce)$sign_SCG)
  df_vari <- as.data.frame(metadata(sce)$sign_VCG)
  tmp <- c("SignID", "Description", "IC",
           "CountCorrGene", "CountWeakCorrGene", "Corr", "CorrWeak",
           "CorrGene", "WeakCorrGene", "CorrGeneID", "WeakCorrGeneID")
  colnames(df_strg) <- tmp
  colnames(df_vari) <- tmp
  if(nrow(df_strg) == 0){
    stop("SCG is not defined.")
  }else{
    if(nrow(df_vari) == 0){
      metadata(sce)$sign_all <- cbind(CorrType = "SCG", df_strg)
    }else{
      metadata(sce)$sign_all <- rbind(cbind(CorrType = "SCG", df_strg),
                                      cbind(CorrType = "VCG", df_vari))
    }
  }
  #--------------------------------------------------
  # Store the report.
  #--------------------------------------------------
  for(k in seq_len(length(s_or_v))){
    if(nrow(metadata(sce)[[s_or_v[k]]]) < 2){
      next
    }
    text <- paste(s_or_v[k], "redundant", sep = "_")
    metadata(sce)[[text]] <- report[[k]]
  }

  return(sce)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
#' Remove signs by specifying keywords.
#'
#' This function removes signs by specifying keywords.
#'
#' @param sce A SingleCellExperiment object.
#' @param keywords keywords separated by pipes `|`.
#'
#' @return A SingleCellExperiment object.
#' @import SingleCellExperiment
#' @import SummarizedExperiment
#' @import S4Vectors
#' @export
#'
#' @examples
#' data(pbmc_eg)
#' data(human_GO_eg)
#' mat <- t(as.matrix(SummarizedExperiment::assay(pbmc_eg, "centered")))
#' pbmc_cormat <- cor(mat, method = "spearman")
#' pbmcs <- list(GO = pbmc_eg)
#' S4Vectors::metadata(pbmcs$GO) <- list(sign = human_GO_eg[["BP"]])
#' pbmcs$GO <- remove_signs(sce = pbmcs$GO, min_ngenes = 2, max_ngenes = 1000)
#' pbmcs$GO <- cluster_genesets(sce = pbmcs$GO, cormat = pbmc_cormat,
#'                              th_posi = 0.24, th_nega = -0.20)
#' pbmcs$GO <- create_signs(sce = pbmcs$GO, min_cnt_strg = 2, min_cnt_vari = 2)
#' keywords <- "Covid19|foofoo|hogehoge"
#' pbmcs$GO <- remove_signs_manually(sce = pbmcs$GO, keywords = keywords)
#' # The results are stored in `metadata(pbmcs$GO)$sign_SCG`,
#' # `metadata(pbmcs$GO)$sign_VCG`, and `metadata(pbmcs$GO)$sign_all`.
#'
remove_signs_manually <- function(sce = NULL, keywords = NULL){
  if(is.null(keywords)){
    return(sce)
  }
  s_or_v <- c("sign_SCG", "sign_VCG")
  for(k in seq_len(length(s_or_v))){
    df <- metadata(sce)[[s_or_v[k]]]
    if(nrow(df) != 0){
      tmp <- ((grepl(keywords, df$SignID)) | (grepl(keywords, df$Description)))
      df <- df[!tmp, ]
      rownames(df) <- seq_len(nrow(df))
    }
    metadata(sce)[[s_or_v[k]]] <- df
  }
  #--------------------------------------------------
  # Modify metadata(sce)$sign_all.
  #--------------------------------------------------
  df_strg <- as.data.frame(metadata(sce)$sign_SCG)
  df_vari <- as.data.frame(metadata(sce)$sign_VCG)
  tmp <- c("SignID", "Description", "IC",
           "CountCorrGene", "CountWeakCorrGene", "Corr", "CorrWeak",
           "CorrGene", "WeakCorrGene", "CorrGeneID", "WeakCorrGeneID")
  colnames(df_strg) <- tmp
  colnames(df_vari) <- tmp
  if(nrow(df_strg) == 0){
    stop("SCG is not defined.")
  }else{
    if(nrow(df_vari) == 0){
      metadata(sce)$sign_all <- cbind(CorrType = "SCG", df_strg)
    }else{
      metadata(sce)$sign_all <- rbind(cbind(CorrType = "SCG", df_strg),
                                      cbind(CorrType = "VCG", df_vari))
    }
  }

  return(sce)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
#' Select signs by specifying keywords.
#'
#' This function selects signs by specifying keywords.
#'
#' @param sce An ASURAT object.
#' @param keywords Keywords separated by a pipe.
#'
#' @return An ASURAT object.
#' @import SingleCellExperiment
#' @import SummarizedExperiment
#' @import S4Vectors
#' @export
#'
#' @examples
#' data(pbmc_eg)
#' data(human_GO_eg)
#' mat <- t(as.matrix(SummarizedExperiment::assay(pbmc_eg, "centered")))
#' pbmc_cormat <- cor(mat, method = "spearman")
#' pbmcs <- list(GO = pbmc_eg)
#' S4Vectors::metadata(pbmcs$GO) <- list(sign = human_GO_eg[["BP"]])
#' pbmcs$GO <- remove_signs(sce = pbmcs$GO, min_ngenes = 2, max_ngenes = 1000)
#' pbmcs$GO <- cluster_genesets(sce = pbmcs$GO, cormat = pbmc_cormat,
#'                              th_posi = 0.24, th_nega = -0.20)
#' pbmcs$GO <- create_signs(sce = pbmcs$GO, min_cnt_strg = 2, min_cnt_vari = 2)
#' keywords <- "cell|process"
#' pbmcs$GO <- select_signs_manually(sce = pbmcs$GO, keywords = keywords)
#' # The results are stored in `metadata(pbmcs$GO)$sign_SCG`,
#' # `metadata(pbmcs$GO)$sign_VCG`, and `metadata(pbmcs$GO)$sign_all`.
#'
select_signs_manually <- function(sce = NULL, keywords = NULL){
  if(is.null(keywords)){
    return(sce)
  }
  s_or_v <- c("sign_SCG", "sign_VCG")
  for(k in seq_len(length(s_or_v))){
    df <- metadata(sce)[[s_or_v[k]]]
    if(nrow(df) != 0){
      tmp <- ((grepl(keywords, df$SignID)) | (grepl(keywords, df$Description)))
      df <- df[tmp, ]
      rownames(df) <- seq_len(nrow(df))
    }
    metadata(sce)[[s_or_v[k]]] <- df
  }
  #--------------------------------------------------
  # Modify metadata(sce)$sign_all.
  #--------------------------------------------------
  df_strg <- as.data.frame(metadata(sce)$sign_SCG)
  df_vari <- as.data.frame(metadata(sce)$sign_VCG)
  tmp <- c("SignID", "Description", "IC",
           "CountCorrGene", "CountWeakCorrGene", "Corr", "CorrWeak",
           "CorrGene", "WeakCorrGene", "CorrGeneID", "WeakCorrGeneID")
  colnames(df_strg) <- tmp
  colnames(df_vari) <- tmp
  if(nrow(df_strg) == 0){
    stop("SCG is not defined. The keywords may not be included in the data.")
  }else{
    if(nrow(df_vari) == 0){
      metadata(sce)$sign_all <- cbind(CorrType = "SCG", df_strg)
    }else{
      metadata(sce)$sign_all <- rbind(cbind(CorrType = "SCG", df_strg),
                                      cbind(CorrType = "VCG", df_vari))
    }
  }

  return(sce)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
#' Create a new SingleCellExperiment object for sign-by-sample matrices.
#'
#' This function creates a new SingleCellExperiment object for sign-by-sample
#'   matrices (SSM) by concatenating SSMs for strongly and variably correlated
#'   gene sets.
#'
#' @param sce A SingleCellExperiment object.
#' @param weight_strg A weight parameter for strongly correlated gene sets.
#' @param weight_vari A weight parameter for variably correlated gene sets.
#'
#' @return A SingleCellExperiment object.
#' @import SingleCellExperiment
#' @import SummarizedExperiment
#' @import S4Vectors
#' @export
#'
#' @examples
#' data(pbmc_eg)
#' data(human_GO_eg)
#' mat <- t(as.matrix(SummarizedExperiment::assay(pbmc_eg, "centered")))
#' pbmc_cormat <- cor(mat, method = "spearman")
#' pbmcs <- list(GO = pbmc_eg)
#' S4Vectors::metadata(pbmcs$GO) <- list(sign = human_GO_eg[["BP"]])
#' pbmcs$GO <- remove_signs(sce = pbmcs$GO, min_ngenes = 2, max_ngenes = 1000)
#' pbmcs$GO <- cluster_genesets(sce = pbmcs$GO, cormat = pbmc_cormat,
#'                              th_posi = 0.24, th_nega = -0.20)
#' pbmcs$GO <- create_signs(sce = pbmcs$GO, min_cnt_strg = 2, min_cnt_vari = 2)
#' pbmcs$GO <- makeSignMatrix(sce = pbmcs$GO, weight_strg = 0.5,
#'                            weight_vari = 0.5)
#' # The resutls can be check by, e.g., assay(pbmcs$GO, "counts").
#'
makeSignMatrix <- function(sce = NULL, weight_strg = 0.5, weight_vari = 0.5){
  #--------------------------------------------------
  # Error handling
  #--------------------------------------------------
  if((weight_strg == 99) & (weight_strg == 99)){
    message("Search a way to bring harmony.")
  }
  if((weight_strg < 0) | (weight_strg > 1) |
     (weight_vari < 0) | (weight_vari > 1)){
    stop("weight_* must be values in [0, 1].")
  }

  mat <- assay(sce, "centered")
  res <- matrix(ncol = dim(mat)[2], nrow = 0,
                dimnames = list(NULL, colnames(mat)))
  rowdata <- c("ParentSignID", "Description", "CorrGene", "WeakCorrGene", "IC")
  rowdata <- matrix(ncol = 5, nrow = 0, dimnames = list(NULL, rowdata))

  s_or_v <- c("sign_SCG", "sign_VCG")
  for(k in seq_len(length(s_or_v))){
    df <- metadata(sce)[[s_or_v[k]]]
    if(nrow(df) != 0){
      for(i in seq_len(nrow(df))){
        #--------------------------------------------------
        # Create a sign-by-sample matrix for weakly correlated gene sets
        #--------------------------------------------------
        if(df$CountWeakCorrGene[i] == 0){
          vector_weak <- matrix(0, ncol = dim(mat)[2], nrow = 1,
                                dimnames = list(NULL, colnames(mat)))
        }else{
          genes_weak <- unlist(strsplit(df$WeakCorrGene[i], "/"))
          submat_weak <- matrix(ncol = dim(mat)[2], nrow = 0,
                                dimnames = list(NULL, colnames(mat)))
          inds <- which(rownames(mat) %in% genes_weak)
          submat_weak <- rbind(submat_weak, mat[inds, ])
          vector_weak <- apply(submat_weak, 2, mean)
        }
        #--------------------------------------------------
        # Create a sign-by-sample matrix for strongly correlated gene sets
        #--------------------------------------------------
        if(is.null(df$CountStrgCorrGene[i])){
          ssm_strg <- NA
        }else{
          if(df$CountStrgCorrGene[i] == 0){
            vector_strg <- matrix(0, ncol = dim(mat)[2], nrow = 1,
                                  dimnames = list(NULL, colnames(mat)))
          }else{
            genes_strg <- unlist(strsplit(df$StrgCorrGene[i], "/"))
            submat_strg <- matrix(ncol = dim(mat)[2], nrow = 0,
                                  dimnames = list(NULL, colnames(mat)))
            inds <- which(rownames(mat) %in% genes_strg)
            submat_strg <- rbind(submat_strg, mat[inds, ])
            vector_strg <- apply(submat_strg, 2, mean)
          }
          ssm_strg <- matrix(ncol = dim(mat)[2], nrow = 0,
                             dimnames = list(NULL, colnames(mat)))
          tmp <- weight_strg * vector_strg + (1.0 - weight_strg) * vector_weak
          ssm_strg <- rbind(ssm_strg, tmp)
          rownames(ssm_strg) <- paste(df$SignID[i], "-", "S", sep = "")
          res <- rbind(res, ssm_strg)
          tmp <- data.frame(ParentSignID = df$SignID[i],
                            Description = df$Description[i],
                            CorrGene = df$StrgCorrGene[i],
                            WeakCorrGene = df$WeakCorrGene[i],
                            IC = df$IC[i])
          rowdata <- rbind(rowdata, tmp)
        }
        #--------------------------------------------------
        # Create a sign-by-sample matrix for variably correlated gene sets
        #--------------------------------------------------
        if(is.null(df$CountVariCorrGene[i])){
          ssm_vari <- NA
        }else{
          if(df$CountVariCorrGene[i] == 0){
            vector_vari <- matrix(0, ncol = dim(mat)[2], nrow = 1,
                                  dimnames = list(NULL, colnames(mat)))
          }else{
            genes_vari <- unlist(strsplit(df$VariCorrGene[i], "/"))
            submat_vari <- matrix(ncol = dim(mat)[2], nrow = 0,
                                  dimnames = list(NULL, colnames(mat)))
            inds <- which(rownames(mat) %in% genes_vari)
            submat_vari <- rbind(submat_vari, mat[inds, ])
            vector_vari <- apply(submat_vari, 2, mean)
          }
          ssm_vari <- matrix(ncol = dim(mat)[2], nrow = 0,
                             dimnames = list(NULL, colnames(mat)))
          tmp <- weight_vari * vector_vari + (1.0 - weight_vari) * vector_weak
          ssm_vari <- rbind(ssm_vari, tmp)
          rownames(ssm_vari) <- paste(df$SignID[i], "-", "V", sep = "")
          res <- rbind(res, ssm_vari)
          tmp <- data.frame(ParentSignID = df$SignID[i],
                            Description = df$Description[i],
                            CorrGene = df$VariCorrGene[i],
                            WeakCorrGene = df$WeakCorrGene[i],
                            IC = df$IC[i])
          rowdata <- rbind(rowdata, tmp)
        }
      }
    }
  }
  new_sce <- SingleCellExperiment(assays = list(counts = res),
                                  rowData = rowdata, colData = colData(sce))
  for(k in seq_len(length(s_or_v))){
    metadata(new_sce)[[s_or_v[k]]] <- metadata(sce)[[s_or_v[k]]]
  }
  if(!is.null(metadata(sce)[["sign_all"]])){
    metadata(new_sce)[["sign_all"]] <- metadata(sce)[["sign_all"]]
  }
  if(length(altExpNames(sce)) != 0){
    for(i in seq_len(length(altExpNames(sce)))){
      altExp(new_sce, altExpNames(sce)[i]) <- altExp(sce, altExpNames(sce)[i])
    }
  }

  return(new_sce)
}

