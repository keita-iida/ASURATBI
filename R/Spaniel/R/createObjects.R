#'@rawNamespace import(Seurat, except = "Assays")
#'@import SingleCellExperiment
#'@import SummarizedExperiment

NULL

#' Create a Seurat Object From Spatial Transcriptomics Data
#'
#' This function converts a count matrix into a Seurat object. The barcodes for
#' each spot are added to the metadata of the Seurat object and are used in
#' plotting the data.
#'
#' @param counts Raw count matrix or data frame where each row represents a
#'               gene and each column represents barcoded location on a spatial
#'               transcriptomics slide. The columns should be named using the
#'               spot barcode (eg "GTCCGATATGATTGCCGC")
#' @param barcodeFile a tab seperated barcode file supplied by Spatial
#'               Trancscriptomics. The file should contains three column:
#'               The first column contains the Spatial Transcriptomics barcode,
#'               the second and third column equate to the x and y location
#' @param projectName The name of the project which is stored in the Seurat
#'                    Object.
#' @param sectionNumber The location of the sample on the slide
#' @return A Seurat Object
#' 
#' @usage  createSeurat(counts, barcodeFile, projectName = projectName,
#'                         sectionNumber = sectionNo)
#'  
#' 
#' @export
#' @examples
#' ## Data is taken from DOI: 10.1126/science.aaf2403
#' examplecounts <- readRDS(file.path(system.file(package = "Spaniel"),
#'                             "extdata/counts.rds"))
#' exampleBarcodes <- file.path(system.file(package = "Spaniel"),
#'                             "1000L2_barcodes.txt")
#' SeuratObj <- createSeurat(examplecounts,
#'                         exampleBarcodes,
#'                         projectName = "TestProj",
#'                         sectionNumber = 1
#'                         )


createSeurat <- function(counts,
                        barcodeFile,
                        projectName = "projectName",
                        sectionNumber = "sectionNo") {
    barcodes <- utils::read.csv(barcodeFile, sep = "\t", header = FALSE)
    rownames(barcodes) <- barcodes$V1
    barcodes <- barcodes[colnames(counts),]
    colnames(barcodes) <- c("spot", "x", "y")
    seuratObj <- Seurat::CreateSeuratObject(counts = counts,
                                            project = paste0(projectName,
                                                                sectionNumber,
                                                                sep = "_"))
    seuratObj@meta.data[, c("spot",
                            "x", "y")] <- barcodes[, c("spot",
                                                        "x", "y")]
    return(seuratObj)
    
}


#' Create a SingleCellExperiment Object From Spatial Transcriptomics Data
#'
#' This function converts a count matrix into a SingleCellExperiment object. 
#' The barcodes for each spot are added to the coldata of the 
#' SingleCellExperiment object and are used in  plotting the data.
#'
#' @param counts Raw count matrix or data frame where each row represents a
#'               gene and each column represents barcoded location on a spatial
#'               transcriptomics slide. The columns should be named using the
#'               spot barcode (eg "GTCCGATATGATTGCCGC")
#' @param barcodeFile a tab seperated barcode file supplied by Spatial
#'               Trancscriptomics. The file should contains three column:
#'               The first column contains the Spatial Transcriptomics barcode,
#'               the second and third column equate to the x and y location
#' @param projectName The name of the project which is stored in the Seurat
#'                    Object.
#' @param sectionNumber The location of the sample on the slide
#' @usage  createSCE(counts, barcodeFile, projectName=projectName, 
#'                 sectionNumber=sectionNo)


#' @return A SingleCellExeriment Object
#' @export
#' @examples
#' ## Data is taken from DOI: 10.1126/science.aaf2403
#' examplecounts <- readRDS(file.path(system.file(package = "Spaniel"),
#'                             "extdata/counts.rds"))
#' exampleBarcodes <- file.path(system.file(package = "Spaniel"),
#'                             "1000L2_barcodes.txt")
#' seuratOb <- createSCE(examplecounts,
#'                         exampleBarcodes,
#'                         projectName = "TestProj",
#'                         sectionNumber = 1)

createSCE <- function(counts,
                    barcodeFile,
                    projectName="projectName",
                    sectionNumber="sectionNo"){
    
    barcodes <- utils::read.csv(barcodeFile, sep="\t", header=FALSE)
    rownames(barcodes) <- barcodes$V1
    colnames(barcodes) <- c("spot", "x", "y")
    
    counts <- counts[, intersect(colnames(counts), barcodes$spot)]
    
    barcodes <- barcodes[colnames(counts), ]
    
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = as.matrix(counts)), 
        colData = barcodes)
    SummarizedExperiment::colData(sce)$project <- paste0(projectName,
                                                        sectionNumber,
                                                        sep = "_")
    
    ### calculate QC metrics
    # sce <- scater::calculateQCMetrics(
    #     sce
    # )
    
    sce <- scater::addPerCellQC(
        sce
    )
    
    return(sce)
    
}
