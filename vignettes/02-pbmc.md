Analysis of peripheral blood mononuclear cell datasets
================
Keita Iida
2022-02-20

-   [1 Install libraries](#install-libraries)
-   [2 Introduction](#introduction)
-   [3 Prepare scRNA-seq data](#prepare-scrna-seq-data)
    -   [3.1 PBMC 4k and 6k](#pbmc-4k-and-6k)
    -   [3.2 PBMCs with and without sepsis (Reyes et
        al., 2020)](#pbmcs-with-and-without-sepsis-reyes-et-al-2020)
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

In this vignette, we analyze single-cell RNA sequencing (scRNA-seq) data
obtained from peripheral blood mononuclear cells (PBMCs) of healthy
donors, and PBMCs of donors with and without bacterial sepsis (Reyes et
al., Nat. Med. 26, 2020).

<br>

# 3 Prepare scRNA-seq data

## 3.1 PBMC 4k and 6k

The data can be loaded by the following code:

``` r
pbmc4k <- readRDS(url("https://figshare.com/ndownloader/files/34112459"))
pbmc6k <- readRDS(url("https://figshare.com/ndownloader/files/34112462"))
```

The data are stored in
[DOI:10.6084/m9.figshare.19200254](https://figshare.com/account/projects/132986/articles/19200254)
and the generating process is described below.

<br>

The data were obtained from 10x Genomics repository (PBMC 4k and 6k).
Create SingleCellExperiment objects by inputting raw read count tables.

``` r
#--------------------------------------------------
# PBMC 4k
#--------------------------------------------------
path_dir <- "rawdata/2020_001_10xgenomics/pbmc_4k/"
path_dir <- paste0(path_dir, "filtered_gene_bc_matrices/GRCh38/")
pbmc4k <- Seurat::Read10X(data.dir = path_dir, gene.column = 2,
                          unique.features = TRUE, strip.suffix = FALSE)
pbmc4k <- SingleCellExperiment(assays = list(counts = as.matrix(pbmc4k)),
                               rowData = data.frame(gene = rownames(pbmc4k)),
                               colData = data.frame(cell = colnames(pbmc4k)))
#--------------------------------------------------
# PBMC 6k
#--------------------------------------------------
path_dir <- "rawdata/2020_001_10xgenomics/pbmc_6k/"
path_dir <- paste(path_dir, "filtered_matrices_mex/hg19/", sep = "")
pbmc6k <- Seurat::Read10X(data.dir = path_dir, gene.column = 2,
                          unique.features = TRUE, strip.suffix = FALSE)
pbmc6k <- SingleCellExperiment(assays = list(counts = as.matrix(pbmc6k)),
                               rowData = data.frame(gene = rownames(pbmc6k)),
                               colData = data.frame(cell = colnames(pbmc6k)))
```

``` r
rbind(dim(pbmc4k), dim(pbmc6k))
```

          [,1] [,2]
    [1,] 33694 4340
    [2,] 32738 5419

<br>

## 3.2 PBMCs with and without sepsis (Reyes et al., 2020)

The data can be loaded by the following code:

``` r
pbmc130k <- readRDS(url("https://figshare.com/ndownloader/files/34112465"))
```

The data are stored in
[DOI:10.6084/m9.figshare.19200254](https://figshare.com/account/projects/132986/articles/19200254)
and the generating process is described below.

<br>

The data were obtained from Broad Institute Single Cell Portal:
[SCP548](https://singlecell.broadinstitute.org/single_cell?type=study&page=1&terms=SCP548).
Since the file size of the raw read count table is huge (\~5.61 GB), we
briefly removed gene and cell data such that the numbers of non-zero
expressing cells are less than 100 and the numbers of sequencing reads
are less than 2000, respectively, by using perl scripts as follows:

``` perl
perl ../perl/pg_01_add_nsamples.pl
perl ../perl/pg_02_add_nreads.pl
perl ../perl/pg_03_remove_variables.pl
perl ../perl/pg_05_remove_samples_2nd.pl
```

Create a SingleCellExperiment object by inputting a raw read count
table.

``` r
fn <- "rawdata/2020_001_Reyes/SCP548/expression/"
fn <- paste0(fn, "scp_gex_matrix_red2.csv")
pbmc130k <- read.csv(fn)
genes <- pbmc130k[-1, 2]
cells <- colnames(pbmc130k)[-seq_len(2)]
pbmc130k <- pbmc130k[-1, -seq_len(2)]
rownames(pbmc130k) <- genes
colnames(pbmc130k) <- cells
pbmc130k <- SingleCellExperiment(assays = list(counts = as.matrix(pbmc130k)),
                                 rowData = data.frame(gene = rownames(pbmc130k)),
                                 colData = data.frame(cell = colnames(pbmc130k)))
```

``` r
dim(pbmc130k)
```

    [1] 14973 60022

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
function `add_metadata()`. The arguments are

-   `sce`: SingleCellExperiment object, and
-   `mitochondria_symbol`: a string representing for mitochondrial
    genes.

``` r
pbmc4k <- add_metadata(sce = pbmc4k, mitochondria_symbol = "^MT-")
pbmc6k <- add_metadata(sce = pbmc6k, mitochondria_symbol = "^MT-")
pbmc130k <- add_metadata(sce = pbmc130k, mitochondria_symbol = "^MT-")
```

One can check the results in `rowData(sce)` and `colData(sce)` slots.

<br>

### 4.1.1 Remove variables based on expression profiles

ASURAT function `remove_variables()` removes variable (gene) data such
that the numbers of non-zero expressing samples (cells) are less than
`min_nsamples`.

``` r
pbmc4k <- remove_variables(sce = pbmc4k, min_nsamples = 10)
pbmc6k <- remove_variables(sce = pbmc6k, min_nsamples = 10)
pbmc130k <- remove_variables(sce = pbmc130k, min_nsamples = 100)
```

<br>

### 4.1.2 Remove samples based on expression profiles

Qualities of sample (cell) data are confirmed based on proper
visualization of `colData(sce)`. ASURAT function `plot_dataframe2D()`
shows scatter plots of two-dimensional data.

``` r
dataframe2D <- data.frame(x = colData(pbmc4k)$nReads, y = colData(pbmc4k)$nGenes)
p <- plot_dataframe2D(dataframe2D = dataframe2D) +
  ggplot2::labs(title = "PBMC 4k", x = "Number of reads", y = "Number of genes") +
  ggplot2::theme_classic(base_size = 20) +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
filename <- "figures/figure_04_0010.png"
ggplot2::ggsave(file = filename, plot = p, dpi = 100, width = 5, height = 5)
```

<img src="figures/figure_04_0010.png" width="200px">
<img src="figures/figure_05_0010.png" width="200px">
<img src="figures/figure_06_0010.png" width="200px">

``` r
dataframe2D <- data.frame(x = colData(pbmc4k)$nReads, y = colData(pbmc4k)$percMT)
p <- plot_dataframe2D(dataframe2D = dataframe2D) +
  ggplot2::labs(title = "PBMC 4k", x = "Number of reads", y = "Perc of MT reads") +
  ggplot2::theme_classic(base_size = 20) +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
filename <- "figures/figure_04_0011.png"
ggplot2::ggsave(file = filename, plot = p, dpi = 100, width = 5, height = 5)
```

<img src="figures/figure_04_0011.png" width="200px">
<img src="figures/figure_05_0011.png" width="200px">
<img src="figures/figure_06_0011.png" width="200px">

ASURAT function `remove_samples()` removes sample (cell) data by setting
cutoff values for the metadata. The arguments are

-   `sce`: SingleCellExperiment object,
-   `min_nReads` and `max_nReads`: minimum and maximum number of reads,
-   `min_nGenes` and `max_nGenes`: minimum and maximum number of
    non-zero expressed genes, and
-   `min_percMT` and `max_percMT`: minimum and maximum percent of reads
    that map to mitochondrial genes, respectively.

``` r
pbmc4k <- remove_samples(sce = pbmc4k, min_nReads = 2200, max_nReads = 20000,
                         min_nGenes = 900, max_nGenes = 1e+10,
                         min_percMT = 0, max_percMT = 10)

pbmc6k <- remove_samples(sce = pbmc6k, min_nReads = 2000, max_nReads = 10000,
                         min_nGenes = 900, max_nGenes = 1e+10,
                         min_percMT = 0, max_percMT = 10)

pbmc130k <- remove_samples(sce = pbmc130k, min_nReads = 2000, max_nReads = 30000,
                           min_nGenes = 1000, max_nGenes = 1e+10,
                           min_percMT = 0, max_percMT = 10)
```

<br>

### 4.1.3 Remove variables based on the mean read counts

Qualities of variable (gene) data are confirmed based on proper
visualization of `rowData(sce)`. ASURAT function `plot_dataframe2D()`
shows scatter plots of two-dimensional data.

``` r
dataframe2D <- data.frame(x = seq_len(nrow(rowData(pbmc4k))),
                          y = sort(rowData(pbmc4k)$nSamples, decreasing = TRUE))
p <- plot_dataframe2D(dataframe2D = dataframe2D) +
  ggplot2::labs(title = "PBMC 4k", x = "Rank of genes", y = "Mean read counts") +
  ggplot2::theme_classic(base_size = 20) +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
filename <- "figures/figure_04_0015.png"
ggplot2::ggsave(file = filename, plot = p, dpi = 100, width = 5, height = 5)
```

<img src="figures/figure_04_0015.png" width="200px">
<img src="figures/figure_05_0015.png" width="200px">
<img src="figures/figure_06_0015.png" width="200px">

ASURAT function `remove_variables_second()` removes variable (gene) data
such that the mean read counts across samples are less than
`min_meannReads`.

``` r
pbmc4k <- remove_variables_second(sce = pbmc4k, min_meannReads = 0.05)
pbmc6k <- remove_variables_second(sce = pbmc6k, min_meannReads = 0.05)
pbmc130k <- remove_variables_second(sce = pbmc130k, min_meannReads = 0.05)
```

``` r
rbind(dim(pbmc4k), dim(pbmc6k), dim(pbmc130k))
```

         [,1]  [,2]
    [1,] 6671  3797
    [2,] 4984  1056
    [3,] 6270 39914

<br>

## 4.2 Normalize data
