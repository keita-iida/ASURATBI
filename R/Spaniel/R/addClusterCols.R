#' @include utilities.R
#'
NULL



#' markClusterCol
#' 
#' A function to mark the columns containing cluster information in the metadata
#' or colData of a Seurat or SCE object. Columns are marked with "cluster_" 
#' prefix.
#' @param object Either a Seurat or SCE object containing clustering information
#' @param pattern pattern indicating which columns contain cluster information
#' @return A Seurat or SCE object
#' @examples 
#' sceObj <- readRDS(file.path(system.file(package = "Spaniel"),
#'                         "extdata/sceData.rds"))
#' sceObj <- markClusterCol(sceObj, "res")
#' @export
markClusterCol <- function(object, pattern) {
    metaData <- getMetadata(object)
    whichCols <- grep(pattern, colnames(metaData))
    colnames(metaData)[whichCols] <- paste0("cluster_",
                                            colnames(metaData)[whichCols])
    object <- updateMetadata(metaData, object)
    return(object)
}



