Analysis of pancreatic ductal adenocarcinoma datasets
================
Keita Iida
2022-02-20

-   [1 Install libraries](#install-libraries)
-   [2 Introduction](#introduction)
-   [3 Prepare scRNA-seq and ST data (Moncada et
    al., 2020)](#prepare-scrna-seq-and-st-data-moncada-et-al-2020)
    -   [3.1 scRNA-seq data](#scrna-seq-data)
    -   [3.2 ST data](#st-data)
-   [4 Preprocessing](#preprocessing)
    -   [4.1 Control data quality](#control-data-quality)
    -   [4.2 Normalize data](#normalize-data)

# 1 Install libraries

Attach necessary libraries:

``` r
library(ASURAT)
library(SingleCellExperiment)
library(SummarizedExperiment)
```

<br>

# 2 Introduction

In this vignette, we analyze single-cell RNA sequencing (scRNA-seq) and
spatial transcriptome (ST) data obtained from primary tumors of
pancreatic ductal adenocarcinoma (PDAC) patients (Moncada et al., Nat.
Biotechnol. 38, 2020).

<br>

# 3 Prepare scRNA-seq and ST data (Moncada et al., 2020)

## 3.1 scRNA-seq data

The data can be loaded by the following code:

``` r
pdacrna <- readRDS(url("https://figshare.com/ndownloader/files/34112468"))
```

The data are stored in
[DOI:10.6084/m9.figshare.19200254](https://figshare.com/account/projects/132986/articles/19200254)
and the generating process is described below.

<br>

From
[GSE111672](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111672),
we downloaded inDrop data with sample accession numbers GSM3036909,
GSM3036910, GSM3405527, GSM3405528, GSM3405529, and GSM3405530 (PDAC-A
inDrop1-inDrop6).

``` r
fn <- c("rawdata/2020_001_Moncada/pdac_indrop/PDACA_1/results/gene_expression.tsv",
        "rawdata/2020_001_Moncada/pdac_indrop/PDACA_2/results/gene_expression.tsv",
        "rawdata/2020_001_Moncada/pdac_indrop/PDACA_3/results/gene_expression.tsv",
        "rawdata/2020_001_Moncada/pdac_indrop/PDACA_4/results/gene_expression.tsv",
        "rawdata/2020_001_Moncada/pdac_indrop/PDACA_5/results/gene_expression.tsv",
        "rawdata/2020_001_Moncada/pdac_indrop/PDACA_6/results/gene_expression.tsv")
pdacrna <- list()
for(i in seq_len(length(fn))){
  d <- read.table(fn[i], header = TRUE, stringsAsFactors = FALSE, row.names = 1)
  colnames(d) <- paste0("PDAC-A-", i, "-", colnames(d))
  pdacrna[[i]] <- SingleCellExperiment(assays = list(counts = as.matrix(d)),
                                       rowData = data.frame(gene = rownames(d)),
                                       colData = data.frame(cell = colnames(d)))
}
```

          [,1]  [,2]
    [1,] 19811 10000
    [2,] 19811 10000
    [3,] 19811 10000
    [4,] 19811 10000
    [5,] 19811 10000
    [6,] 19811 10000

<br>

## 3.2 ST data

The data can be loaded by the following code:

``` r
pdacst <- readRDS(url("https://figshare.com/ndownloader/files/34112471"))
```

The data are stored in
[DOI:10.6084/m9.figshare.19200254](https://figshare.com/account/projects/132986/articles/19200254)
and the generating process is described below.

<br>

Load a raw read count table, convert Ensembl IDs into gene symbols, and
change the column names.

``` r
fn <- "rawdata/2020_001_Moncada/pdac_st/SRR6825057_stdata.tsv"
pdacst <- read.table(fn, header = TRUE, stringsAsFactors = FALSE, row.names = 1)
pdacst <- t(pdacst)
ensembl <- rownames(pdacst)
dictionary <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, key = ensembl,
                                    columns = c("SYMBOL", "ENTREZID"),
                                    keytype = "ENSEMBL")
dictionary <- dictionary[!duplicated(dictionary$ENSEMBL), ]
dictionary[which(is.na(dictionary$SYMBOL)),]$SYMBOL <- as.character("NA")
rownames(pdacst) <- make.unique(as.character(dictionary$SYMBOL))
colnames(pdacst) <- paste0("PDAC-A-ST1_", colnames(pdacst))
```

Create a SingleCellExperiment object by inputting the read count table.

``` r
pdacst <- SingleCellExperiment(assays = list(counts = as.matrix(pdacst)),
                               rowData = data.frame(gene = rownames(pdacst)),
                               colData = data.frame(cell = colnames(pdacst)))
```

A Seurat object, including PDAC tissue images, was obtained from the
authors of
[DOI:10.1038/s41587-019-0392-8](https://doi.org/10.1038/s41587-019-0392-8)
and
[DOI:10.1093/nar/gkab043](https://academic.oup.com/nar/article/49/9/e50/6129341),
and set the tissue image data into the metadata slot.

``` r
fn <- "rawdata/2020_001_Moncada/pdac_st/PDAC-A_ST_list.RDS"
pdacst_surt <- readRDS(file = fn)
pdacst_surt <- pdacst_surt$GSM3036911
metadata(pdacst)$images <- pdacst_surt@images
```

Since the above SingleCellExperiment object includes spatial coordinates
outside of tissues, remove such spots.

``` r
pdacst <- pdacst[, colnames(pdacst_surt)]
identical(colnames(pdacst), colnames(pdacst_surt))
```

    [1] TRUE

``` r
dim(pdacst)
```

    [1] 25807   428

<br>

# 4 Preprocessing

## 4.1 Control data quality

Remove variables (genes) and samples (cells) with low quality, by
processing the following three steps:

1.  remove variables based on expression profiles across samples,
2.  remove samples based on the numbers of reads and nonzero expressed
    variables,
3.  remove variables based on the mean read counts across samples.

First of all, add metadata for both variables and samples using ASURAT
function `add_metadata()`.

The arguments are

-   `sce`: SingleCellExperiment object, and
-   `mitochondria_symbol`: a string representing for mitochondrial
    genes.

``` r
for(i in seq_len(length(pdacrna))){
  pdacrna[[i]] <- add_metadata(sce = pdacrna[[i]], mitochondria_symbol = "^MT-")
}
pdacst <- add_metadata(sce = pdacst, mitochondria_symbol = "^MT-")
```

One can check the results in `rowData(sce)` and `colData(sce)` slots.

<br>

### 4.1.1 Remove variables based on expression profiles

ASURAT function `remove_variables()` removes variable (gene) data such
that the numbers of non-zero expressing samples (cells) are less than
`min_nsamples`.

``` r
for(i in seq_len(length(pdacrna))){
  pdacrna[[i]] <- remove_variables(sce = pdacrna[[i]], min_nsamples = 10)
}
pdacst <- remove_variables(sce = pdacst, min_nsamples = 10)
```

<br>

### 4.1.2 Remove samples based on expression profiles

Qualities of sample (cell) data are confirmed based on proper
visualization of `colData(sce)`. ASURAT function `plot_dataframe2D()`
shows scatter plots of two-dimensional data (see
[here](#visualization_lowdim) for details).

``` r
for(i in seq_len(length(pdacrna))){
  dataframe2D <- data.frame(x = colData(pdacrna[[i]])$nReads,
                            y = colData(pdacrna[[i]])$nGenes)
  p <- plot_dataframe2D(dataframe2D = dataframe2D) +
    ggplot2::labs(title = paste0("PDAC-A inDrop ", i),
                  x = "Number of reads", y = "Number of genes") +
    ggplot2::theme_classic(base_size = 20) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  filename <- paste0("figures/figure_08_0010_", i, ".png")
  ggplot2::ggsave(file = filename, plot = p, dpi = 100, width = 5, height = 5)
}
```

<img src="figures/figure_08_0010_1.png" width="200px">
<img src="figures/figure_08_0010_2.png" width="200px">
<img src="figures/figure_08_0010_3.png" width="200px">
<img src="figures/figure_08_0010_4.png" width="200px">
<img src="figures/figure_08_0010_5.png" width="200px">
<img src="figures/figure_08_0010_6.png" width="200px">

<img src="figures/figure_09_0010.png" width="200px">

``` r
for(i in seq_len(length(pdacrna))){
  dataframe2D <- data.frame(x = colData(pdacrna[[i]])$nReads,
                            y = colData(pdacrna[[i]])$percMT)
  p <- plot_dataframe2D(dataframe2D = dataframe2D) +
    ggplot2::labs(title = paste0("PDAC-A inDrop ", i),
                  x = "Number of reads", y = "Perc of MT reads") +
    ggplot2::theme_classic(base_size = 20) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  filename <- paste0("figures/figure_08_0011_", i, ".png")
  ggplot2::ggsave(file = filename, plot = p, dpi = 100, width = 5, height = 5)
}
```

<img src="figures/figure_08_0011_1.png" width="200px">
<img src="figures/figure_08_0011_2.png" width="200px">
<img src="figures/figure_08_0011_3.png" width="200px">
<img src="figures/figure_08_0011_4.png" width="200px">
<img src="figures/figure_08_0011_5.png" width="200px">
<img src="figures/figure_08_0011_6.png" width="200px">

<img src="figures/figure_09_0011.png" width="200px">

ASURAT function `remove_samples()` removes sample (cell) data by setting
cutoff values for the metadata.

The arguments are

-   `sce`: SingleCellExperiment object,
-   `min_nReads` and `max_nReads`: minimum and maximum number of reads,
-   `min_nGenes` and `max_nGenes`: minimum and maximum number of
    non-zero expressed genes, and
-   `min_percMT` and `max_percMT`: minimum and maximum percent of reads
    that map to mitochondrial genes, respectively.

``` r
for(i in seq_len(length(pdacrna))){
  pdacrna[[i]] <- remove_samples(sce = pdacrna[[i]],
                                 min_nReads = 1000, max_nReads = 10000,
                                 min_nGenes = 500, max_nGenes = 1e+10,
                                 min_percMT = 0, max_percMT = 20)
}

pdacst <- remove_samples(sce = pdacst, min_nReads = 0, max_nReads = 1e+10,
                         min_nGenes = 0, max_nGenes = 1e+10,
                         min_percMT = NULL, max_percMT = NULL)
```

<br>

### 4.1.3 Remove variables based on the mean read counts

Qualities of variable (gene) data are confirmed based on proper
visualization of `rowData(sce)`. ASURAT function `plot_dataframe2D()`
shows scatter plots of two-dimensional data (see
[here](#visualization_lowdim) for details).

``` r
for(i in seq_len(length(pdacrna))){
  dataframe2D <- data.frame(x = seq_len(nrow(rowData(pdacrna[[i]]))),
                            y = sort(rowData(pdacrna[[i]])$nSamples, decreasing = TRUE))
  p <- plot_dataframe2D(dataframe2D = dataframe2D) +
    ggplot2::labs(title = paste0("PDAC-A inDrop ", i),
                  x = "Rank of genes", y = "Mean read counts") +
    ggplot2::theme_classic(base_size = 20) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  filename <- paste0("figures/figure_08_0015_", i, ".png")
  ggplot2::ggsave(file = filename, plot = p, dpi = 100, width = 5, height = 5)
}
```

<img src="figures/figure_08_0015_1.png" width="200px">
<img src="figures/figure_08_0015_2.png" width="200px">
<img src="figures/figure_08_0015_3.png" width="200px">
<img src="figures/figure_08_0015_4.png" width="200px">
<img src="figures/figure_08_0015_5.png" width="200px">
<img src="figures/figure_08_0015_6.png" width="200px">

<img src="figures/figure_09_0015.png" width="200px">

ASURAT function `remove_variables_second()` removes variable (gene) data
such that the mean read counts across samples are less than
`min_meannReads`.

``` r
for(i in seq_len(length(pdacrna))){
  pdacrna[[i]] <- remove_variables_second(sce = pdacrna[[i]], min_meannReads = 0.10)
}
pdacst <- remove_variables_second(sce = pdacst, min_meannReads = 0.10)
```

``` r
rbind(dim(pdacrna[[1]]), dim(pdacrna[[2]]), dim(pdacrna[[3]]),
      dim(pdacrna[[4]]), dim(pdacrna[[5]]), dim(pdacrna[[6]]), dim(pdacst))
```

         [,1] [,2]
    [1,] 6083  393
    [2,] 6151  395
    [3,] 4814  342
    [4,] 5413  311
    [5,] 5737  294
    [6,] 5669  281
    [7,] 4497  428

<br>

## 4.2 Normalize data

Integrate data using bayNorm functions (Tang et al., Bioinformatics,
2020).

``` r
# Take intersection of genes.
genes <- Reduce(intersect, list(rownames(pdacrna[[1]]), rownames(pdacrna[[2]]),
                                rownames(pdacrna[[3]]), rownames(pdacrna[[4]]),
                                rownames(pdacrna[[5]]), rownames(pdacrna[[6]]),
                                rownames(pdacst)))
for(i in seq_len(length(pdacrna))){
  pdacrna[[i]] <- pdacrna[[i]][genes, ]
}
pdacst <- pdacst[genes, ]
# Horizontally concatenate SingleCellExperiment objects.
pdac <- cbind(pdacrna[[1]], pdacrna[[2]], pdacrna[[3]], pdacrna[[4]],
              pdacrna[[5]], pdacrna[[6]], pdacst)
```

and log transformation with a pseudo count. Then, we center the
log-normalized data by subtracting with the mean expression levels
across samples (cells). The resulting normalized-and-centered data are
used for downstream analyses.

Perform `bayNorm()` for attenuating technical biases with respect to
zero inflation and variation of capture efficiencies between samples
(cells).

``` r
bayout <- bayNorm::bayNorm(Data = counts(sclc), mode_version = TRUE)
assay(sclc, "normalized") <- bayout$Bay_out
```

Perform log-normalization with a pseudo count.

``` r
assay(sclc, "logcounts") <- log(assay(sclc, "normalized") + 1)
```

Center row data.

``` r
mat <- assay(sclc, "logcounts")
assay(sclc, "centered") <- sweep(mat, 1, apply(mat, 1, mean), FUN = "-")
```
