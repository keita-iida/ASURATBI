#' Add coordinates to Object
#' Adds output of Spot Detector coordinates to 
#' Seurat object or SCE object created by createSeurat/createSCE. Details about 
#' how to use Spot Detector can be found:
#' https://github.com/SpatialTranscriptomicsResearch/st_spot_detector
#' 
#' @param object either a Seurat object or SCE 
#' @param coordinatesFile path to coordinates file exported from Spot Detector
#' @param scaleFactor a scaling factor which can be used if the image file has 
#' been reduced in size after the coordinates were generated. For example if 
#' the image to be used is 10 percent the size of the original factor
#' scaleFactor = 10
#'
#' @return object
#' @export
#'
#' @examples
#' 
#' ### load a SingleCellExperiment Object 
#' sceObj <- readRDS(file.path(system.file(package = "Spaniel"),
#'                         "extdata/sceData.rds"))
#' ### path to coordinates file exported from spot detector                        
#' coordinatesFile <-  file.path(system.file(package = "Spaniel"),
#'                         "spot_positions.tsv")                        
#' sceObj <- addCoordinates(sceObj, coordinatesFile)  
#' 
addCoordinates <- function(object, coordinatesFile, scaleFactor = NULL){
    testObject(object)
    
    ## read in coordinates
    coordinates <- read.csv(coordinatesFile,
                            sep = "\t")
    
    coordinates$x_y <- paste0(coordinates$x, "_", coordinates$y)
    coordinates <- coordinates[!duplicated(coordinates$x_y), ]

    if (is(object, "Seurat")){
            ### add to metadata
            object$x_y <- paste0(object$x, "_", object$y)
            md <- object[[]] %>% left_join(coordinates)
            rownames(md) <- colnames(object)
            object <- AddMetaData(object, md)
            object$selected[is.na(object$selected)] <- 0
    
    ### scale coordinates
        if (!is.null(scaleFactor)){
            object$pixel_y <- object$pixel_y / scaleFactor
            object$pixel_x <- object$pixel_x / scaleFactor
        }
    
    }
    
    
    if (is(object, "SingleCellExperiment")){
        ### add to metadata
        object$x_y <- paste0(colData(object)$x, "_", colData(object)$y)
        md <- colData(object) %>% data.frame %>% left_join(coordinates)
        rownames(md) <- colnames(object)
        md <-S4Vectors::DataFrame(md)
        colData(object) <- md
        colData(object)$selected[is.na(colData(object)$selected)] <- 0
        
        ### scale coordinates
        
        if (!is.null(scaleFactor)){
            colData(object)$pixel_y <- colData(object)$pixel_y / scaleFactor
            colData(object)$pixel_x <- colData(object)$pixel_x / scaleFactor
        }
        
    }
    
    return(object)  
}
