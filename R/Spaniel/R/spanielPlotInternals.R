#'@importFrom magrittr %>%



# general Functions
# ------------------------------------------------------------------------------

# Grouping the left hand side
ungroupVars <- function(...) {
    List <- as.list(substitute(list(...)))[-1L]
    class(List) <- 'varGroup'
    return(List)
}

# Binary Operator
'%=%' <- function(l, r) {
    envir <- as.environment(-1)
    for (i in seq_len(length(l))) {
        do.call('<-', list(l[[i]], r[[i]]), envir=envir)
    }
}


# Set the variables for plot
# ------------------------------------------------------------------------------
setVars <- function(object, 
                    plotType, 
                    pointpointSize = NULL, 
                    gene = NULL, 
                    clusterRes = NULL){
    if (plotType == "NoGenes") {
        plotTitle <- "Number of Genes Per Spot"
        cl <- "No_Of_Genes"
        sz <- "No_Of_Genes"
        shp <- "NULL"
        colPlot <- ifelse(is(object, "Seurat"),
                            "nFeature_RNA",
                            "detected")
    }
    
    if (plotType == "CountsPerSpot") {
        plotTitle <- "Total Counts Per Spot"
        cl <- "Exprs"
        sz <- "Exprs"
        shp <- "NULL"
        colPlot <- ifelse(is(object, "Seurat"),
                            "nCount_RNA",
                            "sum")
    }
    
    if (plotType == "Cluster") {
        plotTitle <- "Spot Clusters"
        cl <- "Cluster"
        sz <- pointpointSize
        shp <- "Cluster"
        colPlot <- clusterRes
    }
    
    if (plotType == "Gene") {
        plotTitle <- paste("Expression of", gene)
        cl <- "Sign_score"
        sz <- "Sign_score"
        shp <- "NULL"
        colPlot <- gene
    }
    
    showpointSizeLegend <- ifelse(plotType == "Cluster", FALSE, TRUE)
    
    return(c(plotTitle,
                cl,
                sz,
                shp,
                showpointSizeLegend,
                colPlot))
}


# Convert Shape and Size if Null

convertIfNULL <- function(x){
  # convert shp NULL
  if(x == "NULL"){x = NULL}
  return(x)}

# Convert ShowSizeFactor if "TRUE"
convertIfTRUE <- function(y){
  y = ifelse(y == "TRUE", TRUE, FALSE)
  return(y)
}


# Convert size if not equal to character
convertSize <- function(sz){
  if (suppressWarnings(!is.na(as.numeric(sz)))){
  sz = as.numeric(sz)
  }
  return(sz)
}



# Make Coordinates
# ------------------------------------------------------------------------------
###
makeCoordinates <- function(metaData, techType, byCoord, imgDims){
    
    if (techType == "Original" & byCoord == FALSE){
        
        coordinates <- metaData[, c("x", "y")]
        ## reverse the order of the y coordinates
        coordinates$y <- 36 - coordinates$y
    }
    if (techType == "Original" & byCoord == TRUE){
        coordinates <- metaData[, c("pixel_x", "pixel_y")]
        ## reverse the order of the pixel coordinates
        coordinates$pixel_y <- imgDims[2] - coordinates$pixel_y
        colnames(coordinates) <- c("x", "y")
    }
    if (techType == "Visium"){
        coordinates <- metaData[, c("pixel_x", "pixel_y")]
        colnames(coordinates) <- c("x", "y")
    }
    return(coordinates)
}


# Point Range
# ------------------------------------------------------------------------------
###
pointRange <- function(techType, byCoord, imgDims){
    
    if (techType == "Original" & byCoord == FALSE){
            x_min <- y_min <- 1
            x_max <- 33
            y_max <- 35
            
    }

    if (techType == "Visium" | byCoord == TRUE){
        x_min <- y_min <- 0
        x_max <- as.numeric(imgDims[2])
        y_max <- as.numeric(imgDims[1])
    
    }
    return(c(x_min, x_max, y_min, y_max))
}

utils::globalVariables(c("x_min", "x_max", "y_min", "y_max"))

# Create ggplot df for each plot type
# ------------------------------------------------------------------------------
### Make a generic function for all 4 plot types
makeGGDF <- function(object, plotType, colPlot, cl, techType, byCoord, imgDims){
    
    ### Get Metadata and Coodinates
    metaData <- getMetadata(object)
    coordinates <- makeCoordinates(metaData, techType, byCoord, imgDims)

    coordinates$spot <- rownames(metaData)
    
    if (plotType == "Gene"){
        try(if(!colPlot %in% rownames(object))  
            stop(paste0(colPlot,  "not found in rownames(object"))
        )
        # get expression data
        tmp <- getExprs(object)
        tmp <- data.frame(toPlot = tmp[colPlot, ], 
                            spot = rownames(metaData))
    } else {
        try(if(!colPlot %in% colnames(metaData))  
            stop(paste0(colPlot, "not found in colames(object"))
        )
        tmp <- data.frame(toPlot = metaData[,colPlot],
                            spot = rownames(metaData))
    }
    
    # join by spot
    tmp$spot <- as.character(tmp$spot)
    tmp <- dplyr::inner_join(tmp, coordinates, by = "spot") 
    # set colname for the to plot column
    colnames(tmp)[1] <- cl
    
   
    
    return(tmp)
}


# Plot image
# ------------------------------------------------------------------------------
plotImage <- function(grob, tmp, pointColour, pointSize, plotTitle = NULL, 
                      sizeLegend = TRUE, plotType, techType, byCoord, imgDims 
                        ){
    
    p <- ggplot2::ggplot(tmp ,ggplot2::aes_string("x", "y", 
                                                  color = pointColour, 
                                                  size = pointSize)) 
    ## get min and max coordinates for plotting
    plotDims <- pointRange(techType, 
                           byCoord, 
                           imgDims)
    
    ungroupVars(x_min,x_max,y_min,y_max) %=% 
        plotDims
    
    p <- p + ggplot2::xlim(x_min, x_max) +
        ggplot2::ylim(y_min, y_max) +
        ggplot2::annotation_custom(grob, xmin = x_min, xmax = x_max,
                                   ymin = y_min, ymax = y_max)
      ## add theme to plot
    p <- p + ggplot2::geom_point(alpha = 0.6)  +
                ggplot2::labs(title = plotTitle) +
                ggplot2::theme(axis.title.x=ggplot2::element_blank(),
                        axis.text.x=ggplot2::element_blank(),
                        axis.ticks.x=ggplot2::element_blank(),
                        axis.title.y=ggplot2::element_blank(),
                        axis.text.y=ggplot2::element_blank(),
                        axis.ticks.y=ggplot2::element_blank())
    NULL

    ### if show size false
    if (sizeLegend == FALSE){
        p <- p + ggplot2::guides(size=FALSE)}
    ### show plot
     p <- p + ggplot2::guides(color = ggplot2::guide_legend(),
                        size = ggplot2::guide_legend())
    
   
    return(p)
    
}




