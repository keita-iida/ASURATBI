## ----global_options, include=FALSE---------------------------------------
knitr::opts_chunk$set(warning=FALSE, message=FALSE, include = TRUE,
                        fig.height = 8, fig.width = 8, fig.align = "center",
                        echo=TRUE
                        )



## ----install, eval = FALSE-----------------------------------------------
## 
## BiocManager::install('Spaniel')
## 
## 


## ----install_dev, eval = FALSE-------------------------------------------
## 
## devtools::install_github("RachelQueen1/Spaniel", ref = "Dev" )
## 
## 


## ---- load libraries-----------------------------------------------------
library(Spaniel)
library(DropletUtils)
library(scater)



## ----counts--------------------------------------------------------------
pathToTenXOuts <- file.path(system.file(package = "Spaniel"), "extdata/outs")

sce <- createVisiumSCE(tenXDir=pathToTenXOuts, 
                            resolution="Low")



## ----barcodes------------------------------------------------------------
colData(sce)[, c("Barcode", "pixel_x", "pixel_y")]


## ----image_dimensions----------------------------------------------------
metadata(sce)$ImgDims


## ----grob----------------------------------------------------------------
metadata(sce)$Grob


## ---- qcplotting,  results = "hide"--------------------------------------

filter <- sce$detected > 0
spanielPlot(object = sce,
        plotType = "NoGenes", 
        showFilter = filter, 
        techType = "Visium", 
        ptSizeMax = 3)





## ------------------------------------------------------------------------
sce <- sce[, filter]


## ---- gene plot,  results = "hide"---------------------------------------


sce <- logNormCounts(sce)

gene <- "ENSMUSG00000024843"
p2 <- spanielPlot(object = sce,
        plotType = "Gene", 
        gene = "ENSMUSG00000024843",
        techType = "Visium", 
        ptSizeMax = 3)

p2


