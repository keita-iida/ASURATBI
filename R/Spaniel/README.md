## About Spaniel

Spaniel is an R package designed to visualise results of Spatial Transcriptomics
experiments. The current stable version of Spaniel (version 1.1.0) is available from Bioconductor:

```{r install, eval = FALSE}

BiocManager::install('Spaniel')


```

### Spaniel - with 10X import option

This vignette refers to a development version of Spaniel (version 1.2) designed to import data from a 10X Genomics Visium experiment.  
This version will be tested and pushed to Bioconductor. In the meantime, if you would like to test the features described in this vignette you can install a development Spaniel (version 1.2) using the following command:



```{r install_dev, eval = FALSE}

devtools::install_github("RachelQueen1/Spaniel", ref = "Development" )


```



```{r, load libraries}
library(Spaniel)
library(DropletUtils)
library(scater)
```

## Data

This vignette will show how to load the results of 10X Visium spatial 
transcriptomics experiment which has been run through the Space Ranger pipeline.
The data is distributed as part of the Space Ranger software package which can 
be downloaded here:

https://support.10xgenomics.com/spatial-gene-expression/software/overview/welcome

The output from from the "spaceranger testrun" is used as an example here. 


### Import the expression data

Spaniel can be load the output directly from SpaceRanger output using the createVisiumSCE function. This which imports the gene expression data, spatial barcodes, and image dimensions into the SingleCellExperiment object. 

```{r counts}
pathToTenXOuts <- file.path(system.file(package = "Spaniel"), "extData/outs")
sce <- createVisiumSCE(tenXDir=pathToTenXOuts, 
                            resolution="Low")

```

### SCE Object

The pixel coordinates are added to the colData of the SCE object shown below:


```{r barcodes}
colData(sce)[, c("Barcode", "pixel_x", "pixel_y")]
```

The image dimensions are added to the metadata of the SCE object:

```{r image_dimensions}
metadata(sce)$ImgDims
```

The image is stored as a rasterised grob.

```{r grob}
metadata(sce)$Grob
```

## Quality Control

Assessing the number of genes and number of counts per spot is a useful quality 
control step. Spaniel allows QC metrics to be viewed on top of the 
histological image so that any quality issues can be pinpointed. 
Spots within the tissue region which have a low number of genes or counts may 
be due to experimental problems which should be addressed. Conversely spots 
which lie  outside of the tissue and have a high number of counts or large 
number of genes may indicate that there is background contamination. 


### Visualisation

The plotting function allows the use of a binary filter to visualise which 
spots pass filtering thresholds. We create a filter to show spots at 1 gene is detected. Spots where no genes are detected will be removed from the remainder of the analysis.

__NOTE:__ The parameters are set for subset of counts used in this dataset.  
The filter thresholds will be experiment specific and should be adjusted as 
necessary.

```{r, qcplotting,  results = "hide" }

filter <- sce$detected > 0
spanielPlot(object = sce,
        plotType = "NoGenes", 
        showFilter = filter, 
        techType = "Visium", 
        ptSizeMax = 3)



```


Spots where no genes are detected can be removed from the remainder of the analysis.

```{r}
sce <- sce[, filter]
```

The filtered data can then be normalised using the the "normalize" function from scater and the expression of selected genes can be viewed on the histological image.

```{r, gene plot,  results = "hide" }


sce <- logNormCounts(sce)

gene <- "ENSMUSG00000024843"
p2 <- spanielPlot(object = sce,
        plotType = "Gene", gene = gene,
        showFilter = NULL, 
        techType = "Visium", 
        ptSizeMax = 3)

p2

```

