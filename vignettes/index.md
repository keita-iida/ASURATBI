Applying ASURAT to single-cell and spatial transcriptome datasets
================
Keita Iida
2022-06-30

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
-   [Miscellaneous](https://keita-iida.github.io/ASURATBI/07-misc.html)

<br>

# Session information

``` r
sessionInfo()
#> R version 4.2.1 (2022-06-23)
#> Platform: x86_64-apple-darwin17.0 (64-bit)
#> Running under: macOS Big Sur ... 10.16
#> 
#> Matrix products: default
#> BLAS:   /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRblas.0.dylib
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
#> 
#> locale:
#> [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> loaded via a namespace (and not attached):
#>  [1] compiler_4.2.1  magrittr_2.0.3  fastmap_1.1.0   cli_3.3.0      
#>  [5] tools_4.2.1     htmltools_0.5.2 rstudioapi_0.13 yaml_2.3.5     
#>  [9] stringi_1.7.6   rmarkdown_2.14  knitr_1.39      stringr_1.4.0  
#> [13] xfun_0.31       digest_0.6.29   rlang_1.0.3     evaluate_0.15
```
