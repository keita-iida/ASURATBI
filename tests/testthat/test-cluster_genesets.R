#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
data(pbmc_eg)
data(human_GO_eg)
mat <- t(as.matrix(SummarizedExperiment::assay(pbmc_eg, "centered")))
pbmc_cormat <- cor(mat, method = "spearman")
pbmcs <- list(GO = pbmc_eg)
S4Vectors::metadata(pbmcs$GO) <- list(sign = human_GO_eg[["BP"]])
pbmcs$GO <- remove_signs(sce = pbmcs$GO, min_ngenes = 2, max_ngenes = 1000)

test_that("cluster_genesets() runs properly.", {
  set.seed(1)
  pbmcs$GO <- cluster_genesets(sce = pbmcs$GO, cormat = pbmc_cormat,
                               th_posi = 0.24, th_nega = -0.20)
  expect_equal(nrow(S4Vectors::metadata(pbmcs$GO)$sign), 10)
})
