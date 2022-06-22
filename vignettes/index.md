Applying ASURAT to single-cell and spatial transcriptome datasets
================
Keita Iida
2022-06-22

# Table of contents

-   [Analyses of single-cell transcriptome of PBMC 4k from healthy
    donors (10x
    Genomics)](https://keita-iida.github.io/ASURATBI/02-pbmc4k.html)
-   [Analyses of single-cell transcriptome of PBMC 6k from healthy
    donors (10x
    Genomics)](https://keita-iida.github.io/ASURATBI/03-pbmc6k.html)
-   [Analyses of single-cell transcriptomes of PBMCs from control and
    sepsis donors (Reyes et al.,
    2020)](https://keita-iida.github.io/ASURATBI/04-pbmc130k.html)
-   [Additional computations for the PBMC
    datasets](https://keita-iida.github.io/ASURATBI/06-pbmcs.html)
-   [Analyses of single-cell transcriptomes of small cell lung cancer
    (Stewart et al.,
    2020)](https://keita-iida.github.io/ASURATBI/01-sclc.html)
-   [Analyses of single-cell and spatial transcriptomes of pancreatid
    ductal adenocarcinoma (Moncada et al.,
    2020)](https://keita-iida.github.io/ASURATBI/05-pdac.html)

<br>

# Session information

``` r
sessionInfo()
#> R version 4.0.4 (2021-02-15)
#> Platform: x86_64-apple-darwin17.0 (64-bit)
#> Running under: macOS Big Sur 10.16
#> 
#> Matrix products: default
#> BLAS:   /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRblas.dylib
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib
#> 
#> locale:
#> [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> loaded via a namespace (and not attached):
#>  [1] compiler_4.0.4  magrittr_2.0.3  fastmap_1.1.0   cli_3.3.0      
#>  [5] tools_4.0.4     htmltools_0.5.2 rstudioapi_0.13 yaml_2.3.5     
#>  [9] stringi_1.7.6   rmarkdown_2.14  knitr_1.39      stringr_1.4.0  
#> [13] xfun_0.31       digest_0.6.29   rlang_1.0.2     evaluate_0.15
```

\`\`\` R version 4.0.4 (2021-02-15) Platform: x86_64-apple-darwin17.0
(64-bit) Running under: macOS Big Sur 11.1

Matrix products: default LAPACK:
/Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib

locale: \[1\]
en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages: \[1\] parallel stats4 stats graphics grDevices
utils datasets methods base

other attached packages: \[1\] org.Hs.eg.db_3.12.0 AnnotationDbi_1.52.0
SingleCellExperiment_1.12.0 \[4\] SummarizedExperiment_1.20.0
Biobase_2.50.0 GenomicRanges_1.42.0  
\[7\] GenomeInfoDb_1.26.7 IRanges_2.24.1 S4Vectors_0.28.1  
\[10\] BiocGenerics_0.36.1 MatrixGenerics_1.2.1 matrixStats_0.62.0  
\[13\] ASURAT_0.99.15

loaded via a namespace (and not attached): \[1\] circlize_0.4.15
shape_1.4.6 GetoptLong_1.0.5  
\[4\] xfun_0.31 lattice_0.20-45 tcltk_4.0.4  
\[7\] colorspace_2.0-3 vctrs_0.4.1 htmltools_0.5.2  
\[10\] yaml_2.3.5 blob_1.2.3 rlang_1.0.2  
\[13\] DBI_1.1.3 plot3D_1.4 bit64_4.0.5  
\[16\] RColorBrewer_1.1-3 GenomeInfoDbData_1.2.4 zlibbioc_1.36.0  
\[19\] GlobalOptions_0.1.2 memoise_2.0.1 evaluate_0.15  
\[22\] misc3d_0.9-1 knitr_1.39 ComplexHeatmap_2.6.2  
\[25\] fastmap_1.1.0 Cairo_1.5-15 Rcpp_1.0.8.3  
\[28\] cachem_1.0.6 DelayedArray_0.16.3 XVector_0.30.0  
\[31\] bit_4.0.4 rjson_0.2.21 png_0.1-7  
\[34\] digest_0.6.29 grid_4.0.4 clue_0.3-61  
\[37\] cli_3.3.0 tools_4.0.4 bitops_1.0-7  
\[40\] RCurl_1.98-1.7 RSQLite_2.2.14 cluster_2.1.3  
\[43\] pkgconfig_2.0.3 crayon_1.5.1 Matrix_1.4-1  
\[46\] rmarkdown_2.14 compiler_4.0.4\`\`\`
