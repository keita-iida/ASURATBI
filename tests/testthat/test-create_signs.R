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
pbmcs$GO <- cluster_genesets(sce = pbmcs$GO, cormat = pbmc_cormat,
                             th_posi = 0.24, th_nega = -0.20)

test_that("create_signs() runs properly.", {
  pbmcs$GO <- create_signs(sce = pbmcs$GO, min_cnt_strg = 2, min_cnt_vari = 2)
  expect_equal(nrow(S4Vectors::metadata(pbmcs$GO)$sign_SCG), 8)
  expect_equal(nrow(S4Vectors::metadata(pbmcs$GO)$sign_VCG), 2)
  expect_equal(nrow(S4Vectors::metadata(pbmcs$GO)$sign_all), 10)
})
