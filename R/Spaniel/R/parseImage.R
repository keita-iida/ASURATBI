#### PARSE IMAGE ####
#' This function parses a HE image to use as the background for plots
#' @param imgFile Path to the image file
#' @param imgType Type of image options jpg (default), png
#' @return A rasterized grob
#' @export
#' @examples
#' imgFile <- file.path(system.file(package = "Spaniel"),
#'                     "extdata/outs/spatial/tissue_lowres_image.png")
#' img <- parseImage(imgFile, imgType = "png")
parseImage <-  function(imgFile, imgType ="jpg") {
    if (imgType == "jpg"){
        img <- jpeg::readJPEG(imgFile)}
    if (imgType == "png"){
        img <- png::readPNG(imgFile)
    }
    g <-grid::rasterGrob(
        img,
        interpolate = FALSE,
        width = grid::unit(1, "npc"),
        height = grid::unit(1, "npc")
    )
    return(g)
}



