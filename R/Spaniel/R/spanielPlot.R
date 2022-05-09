#' @include spanielPlotInternals.R
#' @include utilities.R
#' 
NULL

#' Spatial Transcriptomics Plot
#'
#' This function overlays information from a Seurat object or 
#' SingleCellExperiment object containing barcodes  onto a H & E image. 
#' There are 4 plots available showing a) the number of
#' genes detected per spot, b) the number of reads detected per spot,
#' c) clustering results, d) the gene expression of a selected gene.
#'
#' @param object Either a Seurat object (version 3) or a SingleCellExperiment 
#' object containing barcode coordinates in the metadata (Seurat) or
#' colData (SingleCellExperiment). 
#' @param techType Either 1) "Original" (default) for the original Spatial 
#' Transcriptomics slides where the image has been cropped to the edge of the 
#' spots 2) "Visium" for 10X slides.
#' @param byCoord TRUE/FALSE option to plot original Spatial Transcriptomics data using pixel 
#' coordinates instead of by spot coordinates. Not required if techType = "Visium".
#' Default is FALSE.
#' @param imgDims pixel dimensions of histological image. Required when 
#' byCoord parameter is set to TRUE, Not required if techType = "Visium".
#' @param plotType There are 5 types of plots avaiable:
#'                       1) NoGenes - This shows the number of genes per spot 
#'                       and uses information from "nFeature_RNA" column of 
#'                       Seurat object or "detected" from a 
#'                       SingleCellExperiment object.
#'                       2) CountsPerSpot - This shows the number of counts per
#'                       spot. It uses information from "nCount_RNA" column of 
#'                       Seurat object or "sum" from a 
#'                       singleCellExperiment object.
#'                       3) Cluster - This plot is designed to show clustering
#'                       results stored in the meta.data or colData of an object
#'                       4) Gene- This plot shows the expression of a single 
#'                       gene. This plot uses scaled/normalised expressin data 
#'                       from the scale.data slot of Seurat object or logcounts 
#'                       of a SingleCellExperiment object.
#'                       5) Other - A generic plot to plot any column from the
#'                       meta.data or colData of an object.
#' @param gene Gene to plot
#' @param clusterRes which cluster resolution to plot 
#' @param ptSize Point size used for cluster plot default is 2
#' @param ptSizeMin Minimum point size used for QC and Gene Expression plots 
#'        default is 0
#' @param ptSizeMax Maximum point size used for QC and Gene Expression plots
#'         default is 5
#' @param customTitle Specify plot title (optional)
#' @param scaleData Show scaled data on plot (default is TRUE)
#' @param showFilter Logical filter showing pass/fail for spots
#' @param grob an grob to be used as the backgound image see(parseImage). This 
#' is used for original Spatial Transcriptomics objects but not Visium                      
#' @return A ggplot spatial transcriptomics plot
#' @export
#' @examples
#' 
#' pathToTenXOuts <- file.path(system.file(package = "Spaniel"), "extdata/outs")
#' sceObj <- createVisiumSCE(tenXDir=pathToTenXOuts, 
#'                        resolution="Low") 
#' filter <- sceObj$detected > 0             
#' spanielPlot(object = sceObj, 
#' plotType = "NoGenes", 
#' showFilter = filter, 
#' techType = "Visium", 
#' ptSizeMax = 3)
#' @usage spanielPlot(object, 
#' grob = NULL, 
#' techType = 
#' "Original", 
#' byCoord = FALSE, 
#' imgDims = NULL, 
#' plotType = c("NoGenes", "CountsPerSpot", "Cluster", "Gene"), 
#' gene= NULL, 
#' clusterRes = NULL, 
#' customTitle = NULL, 
#' scaleData = TRUE,  
#' showFilter = NULL, 
#' ptSize = 2, 
#' ptSizeMin = 0, 
#' ptSizeMax = 5)
#' 
# Main Spaniel Plot Function
# ------------------------------------------------------------------------------
spanielPlot <- function (object, 
                         grob = NULL,
                         techType = "Original",
                         byCoord = FALSE,
                         imgDims = NULL,
                         plotType = c("NoGenes", 
                                      "CountsPerSpot", 
                                      "Cluster", 
                                      "Gene"),
                         gene= NULL, 
                         clusterRes = NULL, 
                         customTitle = NULL,
                         scaleData = TRUE, 
                         showFilter = NULL, 
                         ptSize = 2,
                         ptSizeMin = 0, 
                         ptSizeMax = 5)

{
    ### Validate object is either a Seurat or SCE object
    testObject(object)
    
    if (techType == "Original" & byCoord == TRUE){
        try(if(is.null(imgDims)) 
            stop("image dimensions must be specified to plot by coordinate")
        )
    }
    
    ### set title, colour column, and size column according to plotType
    colPlot = NULL
    ungroupVars(plotTitle,cl,sz,shp, showSizeLegend, colPlot) %=% 
        setVars(object, plotType, ptSize, gene, clusterRes)
    
    shp <- convertIfNULL(shp)
    showSizeLegend <- convertIfTRUE(showSizeLegend)
    sz <- convertSize(sz)
   
    ### create data.frame for ggplot
    tmp <- makeGGDF(object, 
                    plotType, 
                    colPlot, 
                    cl, 
                    techType, 
                    byCoord, 
                    imgDims)
    
    
    ### Update tmp optional arguments if supplied
    if (!is.null(showFilter)){
        tmp$Filter = showFilter
        cl = "Filter"
    }
    
    if (!is.null(customTitle)) {
        plotTitle = customTitle
    }
    
    ### TO DO! add in option for seurat object
    if (techType == "Visium"){
        grob <- S4Vectors::metadata(object)$Grob
        imgDims <- S4Vectors::metadata(object)$ImgDims
    }
    
    ### Create plot
    p <- plotImage(grob, 
                   tmp,
                   pointColour=cl, 
                   pointSize=sz, 
                   plotTitle, 
                   sizeLegend=showSizeLegend,
                   techType = techType,
                   byCoord = byCoord, 
                   imgDims = imgDims)
    
    
    ### FOR QC PLOTS
    if (!is.null(showFilter)){
        p <- p +
            ggplot2::scale_size_continuous(range=c(ptSizeMin,
                                                   ptSizeMax)) +
            ggplot2::guides(color= ggplot2::guide_legend(), 
                            size = ggplot2::guide_legend())
    }
    
    ### FOR GENE PLOTS AND NUMBER READS/COUNTS
    if (plotType != "Cluster" & is.null(showFilter)){
        p <- p +
            ggplot2::scale_size_continuous(range=c(ptSizeMin,
                                                   ptSizeMax))
    }

    ### FOR CLUSTER PLOTS
    if (plotType == "Cluster"){
        p <- p + ggplot2::guides(size=FALSE) +
            ggplot2::scale_size_continuous(range=c(ptSize,
                                                   ptSize))
    }

    return(p)
}

