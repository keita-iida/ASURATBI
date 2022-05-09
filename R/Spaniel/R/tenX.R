## get image dimensions
getImage <- function(spatialDir, resolution){
    if(resolution == "High"){
        img <- png::readPNG(file.path(spatialDir, "tissue_hires_image.png"))
    }else if(resolution == "Low"){
        img <- png::readPNG(file.path(spatialDir, "tissue_lowres_image.png"))
    }
    return(img)
}

## read 10X barcodes
read10XBarcodes <- function(spatialDir, imageDims, resolution){
    barcodes <- read.csv(file.path(spatialDir, "tissue_positions_list.csv"),
                         stringsAsFactors = FALSE, header = FALSE)
    colnames(barcodes) <- c("Barcode", "Section", "Spot_Y", 
                            "Spot_X", "Image_Y", "Image_X")
    #barcodes$barcode <- gsub("-.*", "", barcodes$barcode)
    scaleFactors <- jsonlite::fromJSON(txt = file.path(spatialDir, 
                                             "scalefactors_json.json"))
    
    if(resolution == "High"){
        barcodes$pixel_x <- barcodes$Image_X * scaleFactors$tissue_hires_scalef
        barcodes$pixel_y <- barcodes$Image_Y * scaleFactors$tissue_hires_scalef
    }else if(resolution == "Low"){
        barcodes$pixel_x <- barcodes$Image_X * scaleFactors$tissue_lowres_scalef
        barcodes$pixel_y <- barcodes$Image_Y * scaleFactors$tissue_lowres_scalef
    }
    ## invert Y coordinates 
    barcodes$pixel_y <- imageDims[2] - barcodes$pixel_y
    return(barcodes)
}



parseVisumImage <-  function(imgFile) {
    img <- png::readPNG(imgFile)
    g <-grid::rasterGrob(
        img,
        interpolate = FALSE,
        width = grid::unit(1, "npc"),
        height = grid::unit(1, "npc")
    )
    return(g)
}


#' createVisiumSCE
#' 
#' A function to select to import 10X data into an SCE object
#' @param tenXDir The path to Space Ranger outs directory containing spatial 
#' directory and filtered_feature_bc_matric
#' @param resolution Resolution of the tissue image to be used for plotting. 
#' Can be either "High", or "Low". Default is "Low".
#' @return SingleCellExperimentObject
#' @examples
#' tenXDir <- file.path(system.file(package = "Spaniel"), "extdata/outs")
#' sce <- createVisiumSCE(tenXDir, resolution = "Low")
#' 
#' 
#' @export
createVisiumSCE <- function(tenXDir="../outs", 
                            resolution="Low"){

    filteredDir <- file.path(tenXDir, "filtered_feature_bc_matrix")
    
    sce <- DropletUtils::read10xCounts(filteredDir, 
                                       version = "3")
    
    
    ### parse image, create grob and get image dimensions
    spatialDir =  file.path(tenXDir, "spatial")
    img <- getImage(spatialDir, resolution)
    imgageDims <- dim(img)[c(1,2)]
    grob <-grid::rasterGrob(
        img,
        interpolate = FALSE,
        width = grid::unit(1, "npc"),
        height = grid::unit(1, "npc")
    )
    
    ### add image dimensions to sce
    sce@metadata$ImgDims <- c(imgageDims, resolution)
    
    ### add grob to sce
    sce@metadata$Grob <- grob
    
    ### read barcodes and add to sce
    barcodes <- read10XBarcodes(spatialDir, imgageDims, resolution)
    
    sce@colData <- sce@colData %>% 
        data.frame() %>% 
        left_join(barcodes) %>% 
      S4Vectors::DataFrame()
    
    ## add QC metrics to sce
    per.cell <- scater::perCellQCMetrics(sce)
    colData(sce) <- cbind(colData(sce), per.cell)
    
    
    
    return(sce)
    }




