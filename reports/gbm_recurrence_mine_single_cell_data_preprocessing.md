GBM Recurrence - Mine Single Cell data
================
Miguel Cosenza-Contreras

- <a href="#initial-data-loading-and-wrangling"
  id="toc-initial-data-loading-and-wrangling"><span
  class="toc-section-number">1</span> Initial data loading and
  wrangling</a>
- <a href="#normalize-cpm" id="toc-normalize-cpm"><span
  class="toc-section-number">2</span> Normalize (CPM)</a>
- <a href="#log2-transform" id="toc-log2-transform"><span
  class="toc-section-number">3</span> log2 transform</a>
- <a href="#wrangle" id="toc-wrangle"><span
  class="toc-section-number">4</span> Wrangle</a>

``` r
knitr::opts_chunk$set(echo = TRUE, 
                      message = FALSE, 
                      warning = FALSE)

source(here::here("scr/helper_functions.R"))

## Required packages ----
library(tidyverse)
```

    Warning: package 'ggplot2' was built under R version 4.3.1

    Warning: package 'purrr' was built under R version 4.3.1

    ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ✔ dplyr     1.1.2     ✔ readr     2.1.4
    ✔ forcats   1.0.0     ✔ stringr   1.5.0
    ✔ ggplot2   3.4.3     ✔ tibble    3.2.1
    ✔ lubridate 1.9.2     ✔ tidyr     1.3.0
    ✔ purrr     1.0.2     
    ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ✖ dplyr::filter() masks stats::filter()
    ✖ dplyr::lag()    masks stats::lag()
    ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
library(mixOmics)
```

    Loading required package: MASS

    Attaching package: 'MASS'

    The following object is masked from 'package:dplyr':

        select

    Loading required package: lattice

    Loaded mixOmics 6.24.0
    Thank you for using mixOmics!
    Tutorials: http://mixomics.org
    Bookdown vignette: https://mixomicsteam.github.io/Bookdown
    Questions, issues: Follow the prompts at http://mixomics.org/contact-us
    Cite us:  citation('mixOmics')


    Attaching package: 'mixOmics'

    The following object is masked from 'package:purrr':

        map

``` r
library(fs)
```

    Warning: package 'fs' was built under R version 4.3.1

``` r
library(kableExtra)
```


    Attaching package: 'kableExtra'

    The following object is masked from 'package:dplyr':

        group_rows

``` r
library(sva)
```

    Loading required package: mgcv
    Loading required package: nlme

    Warning: package 'nlme' was built under R version 4.3.1


    Attaching package: 'nlme'

    The following object is masked from 'package:dplyr':

        collapse

    This is mgcv 1.8-42. For overview type 'help("mgcv-package")'.
    Loading required package: genefilter

    Attaching package: 'genefilter'

    The following object is masked from 'package:MASS':

        area

    The following object is masked from 'package:readr':

        spec

    Loading required package: BiocParallel

``` r
library(limma)
library(naniar)
library(missForest)
library(clusterProfiler)
```

    Warning: package 'clusterProfiler' was built under R version 4.3.1


    Registered S3 methods overwritten by 'treeio':
      method              from    
      MRCA.phylo          tidytree
      MRCA.treedata       tidytree
      Nnode.treedata      tidytree
      Ntip.treedata       tidytree
      ancestor.phylo      tidytree
      ancestor.treedata   tidytree
      child.phylo         tidytree
      child.treedata      tidytree
      full_join.phylo     tidytree
      full_join.treedata  tidytree
      groupClade.phylo    tidytree
      groupClade.treedata tidytree
      groupOTU.phylo      tidytree
      groupOTU.treedata   tidytree
      inner_join.phylo    tidytree
      inner_join.treedata tidytree
      is.rooted.treedata  tidytree
      nodeid.phylo        tidytree
      nodeid.treedata     tidytree
      nodelab.phylo       tidytree
      nodelab.treedata    tidytree
      offspring.phylo     tidytree
      offspring.treedata  tidytree
      parent.phylo        tidytree
      parent.treedata     tidytree
      root.treedata       tidytree
      rootnode.phylo      tidytree
      sibling.phylo       tidytree
    clusterProfiler v4.8.2  For help: https://yulab-smu.top/biomedical-knowledge-mining-book/

    If you use clusterProfiler in published research, please cite:
    T Wu, E Hu, S Xu, M Chen, P Guo, Z Dai, T Feng, L Zhou, W Tang, L Zhan, X Fu, S Liu, X Bo, and G Yu. clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. The Innovation. 2021, 2(3):100141

    Attaching package: 'clusterProfiler'

    The following object is masked from 'package:lattice':

        dotplot

    The following object is masked from 'package:MASS':

        select

    The following object is masked from 'package:purrr':

        simplify

    The following object is masked from 'package:stats':

        filter

``` r
library(ReactomePA)
```

    ReactomePA v1.44.0  For help: https://yulab-smu.top/biomedical-knowledge-mining-book/

    If you use ReactomePA in published research, please cite:
    Guangchuang Yu, Qing-Yu He. ReactomePA: an R/Bioconductor package for reactome pathway analysis and visualization. Molecular BioSystems 2016, 12(2):477-479

    Attaching package: 'ReactomePA'

    The following object is masked from 'package:lattice':

        dotplot

``` r
library(org.Hs.eg.db)
```

    Loading required package: AnnotationDbi

    Warning: package 'AnnotationDbi' was built under R version 4.3.1

    Loading required package: stats4
    Loading required package: BiocGenerics

    Attaching package: 'BiocGenerics'

    The following object is masked from 'package:limma':

        plotMA

    The following object is masked from 'package:fs':

        path

    The following objects are masked from 'package:lubridate':

        intersect, setdiff, union

    The following objects are masked from 'package:dplyr':

        combine, intersect, setdiff, union

    The following objects are masked from 'package:stats':

        IQR, mad, sd, var, xtabs

    The following objects are masked from 'package:base':

        anyDuplicated, aperm, append, as.data.frame, basename, cbind,
        colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
        get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
        match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
        Position, rank, rbind, Reduce, rownames, sapply, setdiff, sort,
        table, tapply, union, unique, unsplit, which.max, which.min

    Loading required package: Biobase
    Welcome to Bioconductor

        Vignettes contain introductory material; view with
        'browseVignettes()'. To cite Bioconductor, see
        'citation("Biobase")', and for packages 'citation("pkgname")'.

    Loading required package: IRanges

    Warning: package 'IRanges' was built under R version 4.3.1

    Loading required package: S4Vectors

    Attaching package: 'S4Vectors'

    The following object is masked from 'package:clusterProfiler':

        rename

    The following objects are masked from 'package:lubridate':

        second, second<-

    The following objects are masked from 'package:dplyr':

        first, rename

    The following object is masked from 'package:tidyr':

        expand

    The following object is masked from 'package:utils':

        findMatches

    The following objects are masked from 'package:base':

        expand.grid, I, unname


    Attaching package: 'IRanges'

    The following object is masked from 'package:clusterProfiler':

        slice

    The following object is masked from 'package:nlme':

        collapse

    The following object is masked from 'package:lubridate':

        %within%

    The following objects are masked from 'package:dplyr':

        collapse, desc, slice

    The following object is masked from 'package:purrr':

        reduce

    The following object is masked from 'package:grDevices':

        windows


    Attaching package: 'AnnotationDbi'

    The following object is masked from 'package:clusterProfiler':

        select

    The following object is masked from 'package:MASS':

        select

    The following object is masked from 'package:dplyr':

        select

``` r
library(DT)
library(here)
```

    here() starts at C:/Users/migue/OneDrive/Documentos/R_Projects/5_projects/gbm/gbm_manuscript_data_analysis

``` r
library(janitor)
```


    Attaching package: 'janitor'

    The following objects are masked from 'package:stats':

        chisq.test, fisher.test

``` r
library(drawProteins)
```


    Attaching package: 'drawProteins'

    The following object is masked from 'package:clusterProfiler':

        parse_gff

``` r
library(seqinr)
```


    Attaching package: 'seqinr'

    The following object is masked from 'package:limma':

        zscore

    The following object is masked from 'package:nlme':

        gls

    The following object is masked from 'package:dplyr':

        count

``` r
library(ggpubr)
library(ggrepel)
library(extrafont)
```

    Registering fonts with R

``` r
extrafont::loadfonts(device = "win")
```

    Agency FB already registered with windowsFonts().
    Algerian already registered with windowsFonts().
    Arial Black already registered with windowsFonts().
    Arial already registered with windowsFonts().
    Arial Narrow already registered with windowsFonts().
    Arial Rounded MT Bold already registered with windowsFonts().
    Bahnschrift already registered with windowsFonts().
    Baskerville Old Face already registered with windowsFonts().
    Bauhaus 93 already registered with windowsFonts().
    Bell MT already registered with windowsFonts().
    Berlin Sans FB already registered with windowsFonts().
    Berlin Sans FB Demi already registered with windowsFonts().
    Bernard MT Condensed already registered with windowsFonts().
    Blackadder ITC already registered with windowsFonts().
    Bodoni MT already registered with windowsFonts().
    Bodoni MT Black already registered with windowsFonts().
    Bodoni MT Condensed already registered with windowsFonts().
    Bodoni MT Poster Compressed already registered with windowsFonts().
    Book Antiqua already registered with windowsFonts().
    Bookman Old Style already registered with windowsFonts().
    Bookshelf Symbol 7 already registered with windowsFonts().
    Bradley Hand ITC already registered with windowsFonts().
    Britannic Bold already registered with windowsFonts().
    Broadway already registered with windowsFonts().
    Brush Script MT already registered with windowsFonts().
    Calibri already registered with windowsFonts().
    Calibri Light already registered with windowsFonts().
    Californian FB already registered with windowsFonts().
    Calisto MT already registered with windowsFonts().
    Cambria already registered with windowsFonts().
    Candara already registered with windowsFonts().
    Candara Light already registered with windowsFonts().
    Castellar already registered with windowsFonts().
    Centaur already registered with windowsFonts().
    Century already registered with windowsFonts().
    Century Gothic already registered with windowsFonts().
    Century Schoolbook already registered with windowsFonts().
    Chiller already registered with windowsFonts().
    Colonna MT already registered with windowsFonts().
    Comic Sans MS already registered with windowsFonts().
    Consolas already registered with windowsFonts().
    Constantia already registered with windowsFonts().
    Cooper Black already registered with windowsFonts().
    Copperplate Gothic Bold already registered with windowsFonts().
    Copperplate Gothic Light already registered with windowsFonts().
    Corbel already registered with windowsFonts().
    Corbel Light already registered with windowsFonts().
    Courier New already registered with windowsFonts().
    Curlz MT already registered with windowsFonts().
    Dubai already registered with windowsFonts().
    Dubai Light already registered with windowsFonts().
    Dubai Medium already registered with windowsFonts().
    Ebrima already registered with windowsFonts().
    Edwardian Script ITC already registered with windowsFonts().
    Elephant already registered with windowsFonts().
    Engravers MT already registered with windowsFonts().
    Eras Bold ITC already registered with windowsFonts().
    Eras Demi ITC already registered with windowsFonts().
    Eras Light ITC already registered with windowsFonts().
    Eras Medium ITC already registered with windowsFonts().
    Felix Titling already registered with windowsFonts().
    Footlight MT Light already registered with windowsFonts().
    Forte already registered with windowsFonts().
    Franklin Gothic Book already registered with windowsFonts().
    Franklin Gothic Demi already registered with windowsFonts().
    Franklin Gothic Demi Cond already registered with windowsFonts().
    Franklin Gothic Heavy already registered with windowsFonts().
    Franklin Gothic Medium already registered with windowsFonts().
    Franklin Gothic Medium Cond already registered with windowsFonts().
    Freestyle Script already registered with windowsFonts().
    French Script MT already registered with windowsFonts().
    Gabriola already registered with windowsFonts().
    Gadugi already registered with windowsFonts().
    Garamond already registered with windowsFonts().
    Georgia already registered with windowsFonts().
    Gigi already registered with windowsFonts().
    Gill Sans Ultra Bold already registered with windowsFonts().
    Gill Sans Ultra Bold Condensed already registered with windowsFonts().
    Gill Sans MT already registered with windowsFonts().
    Gill Sans MT Condensed already registered with windowsFonts().
    Gill Sans MT Ext Condensed Bold already registered with windowsFonts().
    Gloucester MT Extra Condensed already registered with windowsFonts().
    Goudy Old Style already registered with windowsFonts().
    Goudy Stout already registered with windowsFonts().
    Haettenschweiler already registered with windowsFonts().
    Harlow Solid Italic already registered with windowsFonts().
    Harrington already registered with windowsFonts().
    High Tower Text already registered with windowsFonts().
    HoloLens MDL2 Assets already registered with windowsFonts().
    Impact already registered with windowsFonts().
    Imprint MT Shadow already registered with windowsFonts().
    Informal Roman already registered with windowsFonts().
    Ink Free already registered with windowsFonts().
    Javanese Text already registered with windowsFonts().
    Jokerman already registered with windowsFonts().
    Juice ITC already registered with windowsFonts().
    Kristen ITC already registered with windowsFonts().
    Kunstler Script already registered with windowsFonts().
    Wide Latin already registered with windowsFonts().
    Lato already registered with windowsFonts().
    Lato Light already registered with windowsFonts().
    Lato Semibold already registered with windowsFonts().
    Leelawadee UI already registered with windowsFonts().
    Leelawadee UI Semilight already registered with windowsFonts().
    Lucida Bright already registered with windowsFonts().
    Lucida Calligraphy already registered with windowsFonts().
    Lucida Console already registered with windowsFonts().
    Lucida Fax already registered with windowsFonts().
    Lucida Handwriting already registered with windowsFonts().
    Lucida Sans already registered with windowsFonts().
    Lucida Sans Typewriter already registered with windowsFonts().
    Lucida Sans Unicode already registered with windowsFonts().
    Magneto already registered with windowsFonts().
    Maiandra GD already registered with windowsFonts().
    Malgun Gothic already registered with windowsFonts().
    Malgun Gothic Semilight already registered with windowsFonts().
    Marlett already registered with windowsFonts().
    Matura MT Script Capitals already registered with windowsFonts().
    Microsoft Himalaya already registered with windowsFonts().
    Microsoft Yi Baiti already registered with windowsFonts().
    Microsoft New Tai Lue already registered with windowsFonts().
    Microsoft PhagsPa already registered with windowsFonts().
    Microsoft Sans Serif already registered with windowsFonts().
    Microsoft Tai Le already registered with windowsFonts().
    Mistral already registered with windowsFonts().
    Modern No. 20 already registered with windowsFonts().
    Mongolian Baiti already registered with windowsFonts().
    Monotype Corsiva already registered with windowsFonts().
    MS Outlook already registered with windowsFonts().
    MS Reference Sans Serif already registered with windowsFonts().
    MS Reference Specialty already registered with windowsFonts().
    MT Extra already registered with windowsFonts().
    MV Boli already registered with windowsFonts().
    Myanmar Text already registered with windowsFonts().
    Niagara Engraved already registered with windowsFonts().
    Niagara Solid already registered with windowsFonts().
    Nirmala UI already registered with windowsFonts().
    Nirmala UI Semilight already registered with windowsFonts().
    OCR A Extended already registered with windowsFonts().
    Old English Text MT already registered with windowsFonts().
    Onyx already registered with windowsFonts().
    Palace Script MT already registered with windowsFonts().
    Palatino Linotype already registered with windowsFonts().
    Papyrus already registered with windowsFonts().
    Parchment already registered with windowsFonts().
    Perpetua already registered with windowsFonts().
    Perpetua Titling MT already registered with windowsFonts().
    Playbill already registered with windowsFonts().
    Poor Richard already registered with windowsFonts().
    Pristina already registered with windowsFonts().
    Rage Italic already registered with windowsFonts().
    Ravie already registered with windowsFonts().
    Roboto Black already registered with windowsFonts().
    Roboto already registered with windowsFonts().
    Roboto Light already registered with windowsFonts().
    Roboto Medium already registered with windowsFonts().
    Roboto Thin already registered with windowsFonts().
    Rockwell already registered with windowsFonts().
    Rockwell Condensed already registered with windowsFonts().
    Rockwell Extra Bold already registered with windowsFonts().
    Script MT Bold already registered with windowsFonts().
    Segoe Fluent Icons already registered with windowsFonts().
    Segoe MDL2 Assets already registered with windowsFonts().
    Segoe Print already registered with windowsFonts().
    Segoe Script already registered with windowsFonts().
    Segoe UI already registered with windowsFonts().
    Segoe UI Light already registered with windowsFonts().
    Segoe UI Semibold already registered with windowsFonts().
    Segoe UI Semilight already registered with windowsFonts().
    Segoe UI Black already registered with windowsFonts().
    Segoe UI Historic already registered with windowsFonts().
    Segoe UI Symbol already registered with windowsFonts().
    Segoe UI Variable already registered with windowsFonts().
    Showcard Gothic already registered with windowsFonts().
    SimSun-ExtB already registered with windowsFonts().
    Sitka Text already registered with windowsFonts().
    Snap ITC already registered with windowsFonts().
    Stencil already registered with windowsFonts().
    Sylfaen already registered with windowsFonts().
    Symbol already registered with windowsFonts().
    Tahoma already registered with windowsFonts().
    Tempus Sans ITC already registered with windowsFonts().
    Times New Roman already registered with windowsFonts().
    Trebuchet MS already registered with windowsFonts().
    Tw Cen MT already registered with windowsFonts().
    Tw Cen MT Condensed already registered with windowsFonts().
    Tw Cen MT Condensed Extra Bold already registered with windowsFonts().
    Verdana already registered with windowsFonts().
    Viner Hand ITC already registered with windowsFonts().
    Vivaldi already registered with windowsFonts().
    Vladimir Script already registered with windowsFonts().
    Webdings already registered with windowsFonts().
    Wingdings already registered with windowsFonts().
    Wingdings 2 already registered with windowsFonts().
    Wingdings 3 already registered with windowsFonts().

``` r
theme_set(theme(axis.text.x = element_text(hjust = 0.5, 
                                           vjust = 0, 
                                           size = 6, 
                                           angle = 360),
                axis.text.y = element_text(hjust = 0.5, 
                                           vjust = 0, 
                                           size = 6),
                panel.background = element_blank(),
                panel.grid.major = element_line(color = "grey"),
                panel.border = element_rect(colour = "black", 
                                            fill = NA, 
                                            size = 1.5),
                axis.title=element_text(size = 8),
                legend.title = element_text(size = 8),
                legend.key.height = unit(3, 
                                         'mm'),
                legend.key.width = unit(3, 
                                        'mm'),
                legend.position = "bottom",
                text = element_text(family = "Helvetica")))
```

    Warning: The `size` argument of `element_rect()` is deprecated as of ggplot2 3.4.0.
    ℹ Please use the `linewidth` argument instead.

# Initial data loading and wrangling

Load the matrix and annotation

``` r
annot_1 <- read_csv(here::here("data/revision_single_cell_data/annot_Human_R_GBM_Full.csv"))

annot_2_ncbi <- read_csv(here::here("data/revision_single_cell_data/from_ncbi/GSM4972210_annot.Human.GBM.R1_2_3_4_4nc.csv.gz"))

# matrix data

matrix_ncbi <- read_csv(here::here("data/revision_single_cell_data/from_ncbi/GSM4972210_Human.GBM.R1_2_3_4_4nc.filtered.gene.bc.matrix.csv"))
```

``` r
View(annot_1)
#View(annot_2_ncbi)
#View(matrix_ncbi)
```

``` r
#head(str(matrix_ncbi))
#names(matrix_ncbi)[1]
```

# Normalize (CPM)

``` r
library(scuttle)
```

# log2 transform

``` r
matrix_ncbi <- matrix_ncbi %>%
  dplyr::rename('gene' = '...1') %>%
  column_to_rownames("gene") %>% 
  as.matrix()
```

``` r
cpm_mat <- calculateCPM(matrix_ncbi)
mat_cpm_log2 <- log2(cpm_mat + 1)
rm(cpm_mat)
```

``` r
mat_cpm_log2[1:15, 1:15]
```

                  AAACCTGAGAAGGTTT-1 AAACCTGAGCCAGGAT-1 AAACCTGCAGCTTAAC-1
    RP11-34P13.7            0.000000           0.000000           0.000000
    FO538757.2              0.000000           0.000000           0.000000
    AP006222.2              0.000000           0.000000           0.000000
    RP4-669L17.10           0.000000           0.000000           0.000000
    RP11-206L10.9           0.000000           0.000000           0.000000
    FAM87B                  0.000000           0.000000           0.000000
    LINC00115               0.000000           0.000000           0.000000
    FAM41C                  0.000000           0.000000           0.000000
    NOC2L                   0.000000           0.000000           0.000000
    KLHL17                  0.000000           0.000000           0.000000
    PLEKHN1                 0.000000           0.000000           0.000000
    HES4                    0.000000           0.000000           0.000000
    ISG15                   8.663376           9.573713           8.704111
    AGRN                    0.000000           0.000000           0.000000
    C1orf159                0.000000           8.575605           0.000000
                  AAACCTGTCATTCACT-1 AAACCTGTCGCCATAA-1 AAACGGGCAGGAATGC-1
    RP11-34P13.7            0.000000            0.00000           0.000000
    FO538757.2              0.000000            0.00000           0.000000
    AP006222.2              0.000000            0.00000           0.000000
    RP4-669L17.10           0.000000            0.00000           0.000000
    RP11-206L10.9           0.000000            0.00000           0.000000
    FAM87B                  0.000000            0.00000           0.000000
    LINC00115               0.000000            0.00000           0.000000
    FAM41C                  0.000000            0.00000           0.000000
    NOC2L                   8.315426            7.65998           0.000000
    KLHL17                  0.000000            0.00000           0.000000
    PLEKHN1                 0.000000            0.00000           0.000000
    HES4                    0.000000            0.00000           0.000000
    ISG15                   0.000000            0.00000           9.431822
    AGRN                    0.000000            0.00000           0.000000
    C1orf159                0.000000            0.00000           0.000000
                  AAACGGGCAGGTCCAC-1 AAACGGGCATATACCG-1 AAACGGGGTCCGTGAC-1
    RP11-34P13.7                   0           0.000000                  0
    FO538757.2                     0           0.000000                  0
    AP006222.2                     0           0.000000                  0
    RP4-669L17.10                  0           0.000000                  0
    RP11-206L10.9                  0           0.000000                  0
    FAM87B                         0           0.000000                  0
    LINC00115                      0           0.000000                  0
    FAM41C                         0           0.000000                  0
    NOC2L                          0           7.975861                  0
    KLHL17                         0           0.000000                  0
    PLEKHN1                        0           0.000000                  0
    HES4                           0           0.000000                  0
    ISG15                          0           7.975861                  0
    AGRN                           0           0.000000                  0
    C1orf159                       0           0.000000                  0
                  AAACGGGTCCGCATCT-1 AAACGGGTCGCGCCAA-1 AAAGATGCATATGAGA-1
    RP11-34P13.7            0.000000           0.000000                  0
    FO538757.2              0.000000           0.000000                  0
    AP006222.2              0.000000           0.000000                  0
    RP4-669L17.10           0.000000           0.000000                  0
    RP11-206L10.9           0.000000           0.000000                  0
    FAM87B                  0.000000           0.000000                  0
    LINC00115               0.000000           0.000000                  0
    FAM41C                  0.000000           0.000000                  0
    NOC2L                   8.203034           0.000000                  0
    KLHL17                  0.000000           0.000000                  0
    PLEKHN1                 0.000000           0.000000                  0
    HES4                    0.000000           0.000000                  0
    ISG15                   9.200584           8.569595                  0
    AGRN                    0.000000           0.000000                  0
    C1orf159                0.000000           0.000000                  0
                  AAAGATGGTGAGCGAT-1 AAAGCAAGTAGCCTCG-1 AAAGTAGAGACACGAC-1
    RP11-34P13.7             0.00000           0.000000            0.00000
    FO538757.2               0.00000           0.000000           10.01485
    AP006222.2               0.00000           0.000000            0.00000
    RP4-669L17.10            0.00000           0.000000            0.00000
    RP11-206L10.9            0.00000           0.000000            0.00000
    FAM87B                   0.00000           0.000000            0.00000
    LINC00115                0.00000           0.000000            0.00000
    FAM41C                   0.00000           0.000000            0.00000
    NOC2L                    0.00000           0.000000            0.00000
    KLHL17                   0.00000           0.000000            0.00000
    PLEKHN1                  0.00000           0.000000            0.00000
    HES4                     0.00000           0.000000            0.00000
    ISG15                    8.17626           8.297274           10.01485
    AGRN                     0.00000           0.000000            0.00000
    C1orf159                 0.00000           0.000000            0.00000

# Wrangle

``` r
mat_cpm_df <- mat_cpm_log2 %>% 
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  relocate(gene, everything()) 
```

``` r
# ASAH1 expression
mat_log2_df_asah1 <- mat_cpm_df %>%
  filter(gene == "ASAH1")
```

``` r
long_log2_df_asah1 <- mat_log2_df_asah1 %>%
  pivot_longer(cols = -gene, 
               names_to = "cell", 
               values_to = "asah1_expression") 
```

Merge with annotation

``` r
long_log2_df_asah1_annot <- long_log2_df_asah1 %>%
  left_join(annot_1, by = c("cell" = "cell"))
```

``` r
#View(long_log2_df_asah1_annot)
```

``` r
saveRDS(long_log2_df_asah1_annot,
        here::here("rds/long_cpm_scrnaseq_asah1_annot.rds"))
```
