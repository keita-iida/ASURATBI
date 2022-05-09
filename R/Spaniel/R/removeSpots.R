#' selectSpots
#' 
#' A function to select spots to remove from analysis
#' @param sObj Either a Seurat object (version 3) or a SingleCellExperiment 
#' object containing barcode coordinates in the metadata (Seurat) or
#' colData (SingleCellExperiment). 
#' @param imgObj a ggplot grob (see parseImage function)
#' @return Runs a shiny application
#' @examples
#' ## Run the shiny app (Not run):
#' # selectSpots(sObj, imgObj)
#' 
#' # Click on the spots to remove from downstream analysis. Once all the spots 
#' # have been selected close the shiny app window. A list of spots is 
#' # stored in a text file called points_to_remove.txt in the working directory.
#' 
#' # Once this step has been run a filtered Seurat or SCE object can be 
#' # created using removeSpots (see removeSpots for more details)
#' 
#' 
#' @export
#' 

selectSpots <- function(sObj, imgObj){
    options(shiny.maxRequestSize=100*1024^2)
    ############## UI #################################    
    ui <- fluidPage(
        # App title ----
        headerPanel("Exclude Spots"),
        fluidRow(plotOutput("plotImage", 
                            click = "plotImage_click")), 
        
        fluidRow(column(width = 6,
                        h4("Points near click"),
                        verbatimTextOutput("click_info")
        )
        )
    )
    
    ############ Server #################################
    server <- function(input, output, session) {
        output$plotImage <- renderPlot({
            ### create plot
            spanielPlot(object = sObj, 
                    grob = imgObj, 
                    plotType =  "NoGenes",
                    customTitle = NULL, 
                    scaleData = TRUE)
        })
        
        output$click_info <- renderPrint({
            metadata <- getMetadata(sObj)
            coords <-getCoordinates(metadata)
            coords$y <- 36 - coords$y
            pts <- nearPoints(coords, input$plotImage_click)
            write(as.character(pts$spot), 
                    "points_to_remove.txt", append = TRUE)
            
        })
        
        
        
        # stop App on closing the browser
        session$onSessionEnded(function() {
            stopApp()
        })
    }
    shinyApp(ui, server)
}


#' removeSpots
#' 
#' A function to filter spots from analysis. It requires selectSpots to be 
#' run first.
#' @param sObj Either a Seurat object (version 3) or a SingleCellExperiment 
#' object containing barcode coordinates in the metadata (Seurat) or
#' colData (SingleCellExperiment).  
#' @param pointsToRemove path to points to remove file. Default is 
#' "points_to_remove.txt" 
#' @return A filtered Seurat or SingleCellExperiment Object
#' @examples
#' sceObj <- readRDS(file.path(system.file(package = "Spaniel"),
#'                         "extdata/sceData.rds"))
#' toRemove <- file.path(system.file(package = "Spaniel"),
#'                         "points_to_remove.txt")
#' sceObj_filtered <- removeSpots(sObj = sceObj, pointsToRemove = toRemove)
#' @export
#' @usage  removeSpots(sObj, pointsToRemove = "points_to_remove.txt")

removeSpots <-function(sObj, 
                        pointsToRemove = "points_to_remove.txt"){
    toRemove = read.csv(pointsToRemove, 
                        header = FALSE, 
                        stringsAsFactors = FALSE)$V1
    toKeep = setdiff(colnames(sObj), toRemove)
    filter = colnames(sObj) %in% toRemove
    
    objFiltered <- sObj[,toKeep]
    
    return(objFiltered)
    
}












