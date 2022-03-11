#-----------------------------------------------------------------------------80
# 
#-----------------------------------------------------------------------------80
#' A SingleCellExperiment object made from a gene expression table.
#'
#' A SingleCellExperiment object, including 50 genes and 50 cells.
#'   The original data "4k PBMCs from a Healthy Donor" was downloaded from 10x
#'   Genomics database.
#'
#' @format SingleCellExperiment object.
#' @source \url{https://support.10xgenomics.com/single-cell-gene-expression}
#'
"pbmc_eg"
#-----------------------------------------------------------------------------80
# 
#-----------------------------------------------------------------------------80
#' A list of SingleCellExperiment objects made from sign-sample matrices.
#'
#' A list of SingleCellExperiment objects, consisting of small sign-by-sample
#'   matrices, pbmcs_eg$CM (using Cell Ontology and MSigDB databases),
#'   pbmcs_eg$GO (using Gene Ontology database), and pbmcs_eg$KG (KEGG).
#'   Here, pbmcs_eg$CM, pbmcs_eg$GO, and pbmcs_eg$KG include 87, 72, and 64
#'   signs, respectively, and 50 cells.
#'
#' @format A list of SingleCellExperiment objects.
#'
"pbmcs_eg"
#-----------------------------------------------------------------------------80
# 
#-----------------------------------------------------------------------------80
#' A list of small Cell Ontology and MSigDB databases for human.
#'
#' A list of small Cell Ontology and MSigDB databases for human.
#'
#' @format A list of dataframe.
#'
"human_COMSig_eg"
#-----------------------------------------------------------------------------80
# 
#-----------------------------------------------------------------------------80
#' A list of small Gene Ontology database for human.
#'
#' A list of small Gene Ontology database for human.
#'
#' @format A list of dataframe.
#'
"human_GO_eg"
#-----------------------------------------------------------------------------80
# 
#-----------------------------------------------------------------------------80
#' A list of small KEGG database for human.
#'
#' A list of small KEGG database for human.
#'
#' @format A list of dataframe.
#'
"human_KEGG_eg"
