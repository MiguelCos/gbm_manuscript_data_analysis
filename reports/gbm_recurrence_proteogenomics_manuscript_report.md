GBM Recurrence Proteogenomics
================
Miguel Cosenza-Contreras

- <a href="#initial-data-loading-and-wrangling"
  id="toc-initial-data-loading-and-wrangling"><span
  class="toc-section-number">1</span> Initial data loading and
  wrangling</a>
- <a
  href="#localization-of-identified-peptides-and-mapping-to-variant-regions"
  id="toc-localization-of-identified-peptides-and-mapping-to-variant-regions"><span
  class="toc-section-number">2</span> Localization of identified peptides
  and mapping to variant regions</a>
  - <a
    href="#localize-peptides-within-protein-sequence-and-confirm-identification-of-saavs"
    id="toc-localize-peptides-within-protein-sequence-and-confirm-identification-of-saavs"><span
    class="toc-section-number">2.1</span> Localize peptides within protein
    sequence and confirm identification of SAAVs</a>
  - <a href="#generate-table-of-identified-peptides-mapping-to-saavs"
    id="toc-generate-table-of-identified-peptides-mapping-to-saavs"><span
    class="toc-section-number">2.2</span> Generate table of identified
    peptides mapping to SAAVs</a>
- <a href="#visualize-quantitative-features-of-identified-saavs"
  id="toc-visualize-quantitative-features-of-identified-saavs"><span
  class="toc-section-number">3</span> Visualize quantitative features of
  identified SAAVs</a>
  - <a href="#prep-abundance-matrices"
    id="toc-prep-abundance-matrices"><span
    class="toc-section-number">3.1</span> Prep abundance matrices</a>
  - <a href="#filter-merge-and-summarize-saav-abundances"
    id="toc-filter-merge-and-summarize-saav-abundances"><span
    class="toc-section-number">3.2</span> Filter, merge and summarize SAAV
    abundances</a>
    - <a href="#generate-summarized-results-tables"
      id="toc-generate-summarized-results-tables"><span
      class="toc-section-number">3.2.1</span> Generate summarized results
      tables</a>
  - <a href="#heatmap-of-all-saavs-observed"
    id="toc-heatmap-of-all-saavs-observed"><span
    class="toc-section-number">3.3</span> Heatmap of all SAAVs observed</a>
  - <a href="#heatmap-of-saavs-observed-on-every-sample"
    id="toc-heatmap-of-saavs-observed-on-every-sample"><span
    class="toc-section-number">3.4</span> Heatmap of SAAVs observed on every
    sample</a>
- <a href="#limma-on-saavs" id="toc-limma-on-saavs"><span
  class="toc-section-number">4</span> Limma on SAAVs</a>
  - <a href="#prep-design-matrix" id="toc-prep-design-matrix"><span
    class="toc-section-number">4.1</span> Prep design matrix</a>
  - <a href="#volcano-saavs" id="toc-volcano-saavs"><span
    class="toc-section-number">4.2</span> Volcano SAAVs</a>
  - <a href="#boxplot-of-differentially-abundant-saav"
    id="toc-boxplot-of-differentially-abundant-saav"><span
    class="toc-section-number">4.3</span> Boxplot of differentially abundant
    SAAV</a>
  - <a href="#boxplot-da-saav" id="toc-boxplot-da-saav"><span
    class="toc-section-number">4.4</span> Boxplot DA SAAV</a>

``` r
knitr::opts_chunk$set(echo = TRUE, 
                      message = FALSE, 
                      warning = FALSE)

source(here::here("scr/helper_functions.R"))

## Required packages ----
library(tidyverse)
```

    ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.2 ──
    ✔ ggplot2 3.4.0      ✔ purrr   0.3.4 
    ✔ tibble  3.1.8      ✔ dplyr   1.0.10
    ✔ tidyr   1.2.1      ✔ stringr 1.4.1 
    ✔ readr   2.1.2      ✔ forcats 0.5.2 

    Warning: package 'ggplot2' was built under R version 4.2.2

    ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ✖ dplyr::filter() masks stats::filter()
    ✖ dplyr::lag()    masks stats::lag()

``` r
library(mixOmics)
```

    Loading required package: MASS

    Attaching package: 'MASS'

    The following object is masked from 'package:dplyr':

        select

    Loading required package: lattice

    Loaded mixOmics 6.20.0
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

    Attaching package: 'nlme'

    The following object is masked from 'package:dplyr':

        collapse

    This is mgcv 1.8-40. For overview type 'help("mgcv-package")'.
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


    clusterProfiler v4.4.4  For help: https://yulab-smu.top/biomedical-knowledge-mining-book/

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

    ReactomePA v1.40.0  For help: https://yulab-smu.top/biomedical-knowledge-mining-book/

    If you use ReactomePA in published research, please cite:
    Guangchuang Yu, Qing-Yu He. ReactomePA: an R/Bioconductor package for reactome pathway analysis and visualization. Molecular BioSystems 2016, 12(2):477-479

    Attaching package: 'ReactomePA'

    The following object is masked from 'package:lattice':

        dotplot

``` r
library(org.Hs.eg.db)
```

    Loading required package: AnnotationDbi
    Loading required package: stats4
    Loading required package: BiocGenerics

    Attaching package: 'BiocGenerics'

    The following object is masked from 'package:limma':

        plotMA

    The following object is masked from 'package:fs':

        path

    The following objects are masked from 'package:dplyr':

        combine, intersect, setdiff, union

    The following objects are masked from 'package:stats':

        IQR, mad, sd, var, xtabs

    The following objects are masked from 'package:base':

        anyDuplicated, append, as.data.frame, basename, cbind, colnames,
        dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
        grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
        order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
        rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
        union, unique, unsplit, which.max, which.min

    Loading required package: Biobase
    Welcome to Bioconductor

        Vignettes contain introductory material; view with
        'browseVignettes()'. To cite Bioconductor, see
        'citation("Biobase")', and for packages 'citation("pkgname")'.

    Loading required package: IRanges
    Loading required package: S4Vectors

    Attaching package: 'S4Vectors'

    The following object is masked from 'package:clusterProfiler':

        rename

    The following objects are masked from 'package:dplyr':

        first, rename

    The following object is masked from 'package:tidyr':

        expand

    The following objects are masked from 'package:base':

        expand.grid, I, unname


    Attaching package: 'IRanges'

    The following object is masked from 'package:clusterProfiler':

        slice

    The following object is masked from 'package:nlme':

        collapse

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

Load id/quant data and prep coverage information

``` r
## tmt-integrator output loading 
tmt_protdata_frag_progen <- read_tsv(here("data/proteogenomics/sras_cellreports_2021_paper_razor/tmt-report/abundance_protein_MD.tsv")) %>% 
  janitor::clean_names()

tmt_peptdata_frag_progen <- read_tsv(here("data/proteogenomics/sras_cellreports_2021_paper_razor/tmt-report/abundance_peptide_MD.tsv")) %>% 
  janitor::clean_names()

identified_proteins_progen <- read_tsv(here("data/proteogenomics/sras_cellreports_2021_paper_razor/combined_protein.tsv")) %>% 
  janitor::clean_names()

identified_peptides_progen <- read_tsv(here("data/proteogenomics/sras_cellreports_2021_paper_razor/combined_peptide.tsv")) %>% 
  janitor::clean_names()

index2genepeptide <- tmt_peptdata_frag_progen %>%
  dplyr::select(index, gene, peptide, protein_id)
```

``` r
sample_annotation <- read_csv(here("data/sample_annotation.csv"))

# correct annotation

sample_annotation2 <- sample_annotation %>%
  mutate(patient = paste("x",
                         patient, 
                         sep = ""),
         recurrence = case_when(recurrence == "initial" ~ "prim",
                                recurrence == "recurrent" ~ "rec",
                                TRUE ~ recurrence)) %>%
  mutate(paired_id = paste(patient, 
                           recurrence, 
                           sep = "_")) %>%
  filter(recurrence %in% c("prim", 
                           "rec"))

sample_annotationwo6 <- sample_annotation2 %>%
  filter(patient != "x6")
```

``` r
# wide expression matrix with NAs
expr_matrix_progen <- dplyr::select(tmt_protdata_frag_progen,
                             Protein = index, 
                             starts_with("x"))

long_mat_progen <- pivot_longer(expr_matrix_progen,
                         cols = starts_with("x"),
                         names_to = c("paired_id_semi"),
                         values_to = "Abundance") %>%
  mutate(recurrence = str_extract(paired_id_semi, 
                                  "[^0-9]+$")) %>%
  mutate(patient = str_remove_all(paired_id_semi, 
                                  "[^0-9]+$"))

long_mat_progen_2 <- long_mat_progen %>%
  mutate(paired_id = paste(patient, 
                           recurrence, 
                           sep = "_")) 

paired_annot_uniq_progen <- long_mat_progen_2 %>%
  dplyr::select(patient, 
                recurrence, 
                paired_id) %>%
  distinct()

# wide peptide abundace matrix with NAs
expr_matrix_pept_progen <- dplyr::select(tmt_peptdata_frag_progen,
                             index,
                             peptide,
                             starts_with("x")) %>%
  dplyr::select(!starts_with("x6"))

expr_matrixnona_pept_progen <- expr_matrix_pept_progen %>%
  na.omit()

long_mat_pept_progen <- pivot_longer(expr_matrix_pept_progen,
                         cols = starts_with("x"),
                         names_to = c("paired_id_semi"),
                         values_to = "Abundance") %>%
  mutate(recurrence = str_extract(paired_id_semi, 
                                  "[^0-9]+$")) %>%
  mutate(patient = str_remove_all(paired_id_semi, 
                                  "[^0-9]+$"))

long_mat_pept_progen_2 <- long_mat_pept_progen %>%
  mutate(paired_id = paste(patient, 
                           recurrence, 
                           sep = "_")) 

paired_annot_uniq_pept_progen <- long_mat_pept_progen_2 %>%
  dplyr::select(patient, 
                recurrence, 
                paired_id) %>%
  distinct()

## Annotated abundance data in long format ----

quant_annotated_progen <- left_join(long_mat_progen_2,
                             sample_annotation2) %>%
  mutate(Channel_mix = paste(mixture,
                             channel, 
                             sep = "_"))

quant_annotated_pept_progen <- left_join(long_mat_pept_progen_2,
                             sample_annotation2) %>%
  mutate(Channel_mix = paste(mixture,
                             channel, 
                             sep = "_"))


## Wide xpression matrices without NAs
expr_matnona_progen <- column_to_rownames(expr_matrix_progen,
                            "Protein") %>%
          as.matrix() %>%
          na.omit()

# peptide to protein mapping from combined_peptides.tsv 
protein2peptideall <- dplyr::select(identified_peptides_progen,
                                    protein, 
                                    Peptide = sequence) %>% 
          separate(col = protein,
                   into = c("sp", 
                            "Genes", 
                            "ensembl"),
                   sep = "\\|") %>%
          dplyr::select(Genes, 
                        Peptide) %>%
          distinct() %>%
          mutate(Peptide = str_remove_all(Peptide, 
                                          "\\[.+?\\]")) %>%
          mutate(Peptide = str_remove_all(Peptide, 
                                          "^n"))

# peptide to protein mappong from tmt abundance report (quantified and normalized)
protein2peptide_tmtrept <- tmt_peptdata_frag_progen %>%
  dplyr::select(index, 
                protein_id, 
                peptide)
```

# Localization of identified peptides and mapping to variant regions

``` r
# list of identified peptides associated with 'variant' proteins
var_prots <- protein2peptideall %>%
                    mutate(variant = if_else(str_detect(Genes, 
                                                        "\\_"),
                           true = "variant",
                           false = "canonical")) %>%
                    filter(variant == "variant") %>% 
                    pull(Genes) %>%
                    unique()
```

``` r
variant_peptides_1 <- filter(protein2peptideall,
                             Genes %in% var_prots)

# prep fasta containing only identified protein sequences
fastavar <- seqinr::read.fasta(here("data/proteogenomics/sras_cellreports_2021_paper_razor/mix_1/protein.fas"), 
                            seqtype = "AA", 
                            as.string = TRUE)

fastavar2 <- seqinr::read.fasta(here("data/proteogenomics/sras_cellreports_2021_paper_razor/mix_2/protein.fas"), 
                            seqtype = "AA", 
                            as.string = TRUE)

fastavar3 <- seqinr::read.fasta(here("data/proteogenomics/sras_cellreports_2021_paper_razor/mix_3/protein.fas"), 
                            seqtype = "AA", 
                            as.string = TRUE)

fasta1 <- append(fastavar,
                fastavar2)

fasta <- append(fasta1,
                fastavar3)

fasta <- discard(fasta,
                 duplicated(names(fasta)))
```

``` r
if(!file.exists(here("results/proteogenomics/cellreps_2021/cleavage_annoated_peptides_razor.tsv"))){
  cleavage_annoated_peptides <- annotate_peptides(expr_mat = variant_peptides_1, 
                                                  fasta = fasta,
                                                  decoy_tag = "rev_")
  
  write_tsv(cleavage_annoated_peptides, "results/proteogenomics/cellreps_2021/cleavage_annoated_peptides_razor.tsv")
} else {
  cleavage_annoated_peptides <- read_tsv(here("results/proteogenomics/cellreps_2021/cleavage_annoated_peptides_razor.tsv"))
}
```

## Localize peptides within protein sequence and confirm identification of SAAVs

``` r
# tabulate location of the variant
annotated_w_variant_location <- cleavage_annoated_peptides %>%
          as_tibble() %>%
          mutate(protein_id_var = protein_id,
                 variant_type = if_else(str_detect(protein_id, 
                                                   "\\:"), 
                                        true = "other",
                                        false = "saav")) %>%
          mutate(protein_id = str_replace_all(protein_id, 
                                              "\\.", 
                                              "\\_")) %>%
          mutate(protein_id = str_replace_all(protein_id, 
                                              "\\:", 
                                              "\\_")) %>%
          separate(protein_id, 
                   into = c("ense_id", "variant_location"), 
                   sep = "\\_",
                   extra = "merge") %>%
                    na.omit() %>%
          mutate(variant_type = if_else(condition = str_detect(variant_location, 
                                                               "[0-9]$"),
                                        true = "other",
                                        false = variant_type))

# generate a long list of variant locations 

locat2protdf <- annotated_w_variant_location %>%
  separate_rows(variant_location)

long2prot <- locat2protdf %>%
          mutate(variant_num_location = readr::parse_number(variant_location))

long2prot_worivar <- left_join(long2prot, 
                               dplyr::select(annotated_w_variant_location,
                                             protein_id_var, 
                                             ense_id)) %>%
  distinct()


annotated_wo_variant_location <- annotated_w_variant_location %>%
          dplyr::select(-variant_location)

# merge long table of variant locations with the annotation of peptide position within 
# their associated protein sequences

protein2varlocat <- left_join(annotated_wo_variant_location, 
                              long2prot_worivar) %>%
          distinct() %>%
          mutate(start_num = as.numeric(start_position),
                 end_num = as.numeric(end_position)) %>%
  # evaluate if the SAAV falls on the position of a variant sequence
          mutate(variant_peptide = variant_num_location >= as.numeric(start_position) & variant_num_location <= as.numeric(end_position)) %>% 
          dplyr::select(protein_id_var, 
                        ense_id, 
                        protein_description, 
                        Peptide, 
                        start_position, 
                        start_num,
                        end_position,
                        end_num, 
                        variant_location, 
                        variant_num_location, 
                        variant_peptide, 
                        variant_type)

pept2variant <- protein2varlocat %>%
  dplyr::select(Peptide, 
                variant_peptide, 
                variant_type,
                protein_id_var,
                ense_id,
                variant_location,
                peptide_start_position = start_position,
                peptide_end_position = end_position,
                protein_description) %>%
  distinct() %>%
  na.omit() %>%
  filter(variant_type == "saav") %>% 
  mutate(duplicated_peptide = duplicated(Peptide)) %>% 
  relocate(duplicated_peptide,
           .after = variant_type)
```

## Generate table of identified peptides mapping to SAAVs

``` r
saavs_table <- pept2variant %>%
  filter(variant_peptide == TRUE) %>% 
  mutate(duplicated_peptide = duplicated(Peptide)) %>% 
  dplyr::rename(peptide = Peptide) %>%
  left_join(., index2genepeptide) %>%
  mutate(gene_variant = paste(gene, 
                              variant_location, 
                              sep = "_")) 
```

``` r
length(unique(saavs_table$Peptide))
```

    [1] 0

``` r
length(unique(saavs_table$protein_id_var))
```

    [1] 188

``` r
saavs_suppl_table <- saavs_table %>% 
  dplyr::select(gene_variant,
                variant_type,
                variant_location,
                peptide,
                peptide_start_position,
                peptide_end_position,
                protein_id_var,
                ense_id,
                gene,
                protein_description)
```

``` r
write_tsv(x = saavs_suppl_table,
          file = here("suppl_tables/supplementary_table_5_proteogenomics_identified_saavs.tsv"))
```

# Visualize quantitative features of identified SAAVs

## Prep abundance matrices

``` r
recurrencewo6 <- sample_annotationwo6$recurrence
patientwo6 <- sample_annotationwo6$patient
```

``` r
pep_prematwo6 <- dplyr::select(expr_matrix_pept_progen,
                               index, 
                               peptide, 
                               rownames(design_limmawo6))  

pep_matwo6 <- dplyr::select(expr_matrix_pept_progen,
                            index, 
                            rownames(design_limmawo6))  %>%
  column_to_rownames("index") %>%
  as.matrix()

pep_matnona_wo6 <- pep_matwo6 %>%
  na.omit()
```

## Filter, merge and summarize SAAV abundances

``` r
# extract indexing, peptide and sequence information from normalized data
index2genepeptide <- tmt_peptdata_frag_progen %>% 
  dplyr::select(index, 
                peptide, 
                gene, 
                protein_id)

# table of interesting features (saavs)
features1 <- pept2variant %>%
  dplyr::rename(peptide = Peptide) %>%
  left_join(., index2genepeptide) %>%
  na.omit()

interesting_features1 <- features1 %>%
                    filter(variant_peptide == TRUE) 

# data frame of identified and TMT-integrator summarized SAAVs + annotation
protein2varlocat_saav <- protein2varlocat %>%
  dplyr::rename(peptide = Peptide) %>%
  filter(peptide %in% interesting_features1$peptide,
         variant_peptide == TRUE,
         variant_type == "saav") %>%
  left_join(., index2genepeptide) %>%
  mutate(gene_variant = paste(gene, 
                              variant_location, 
                              sep = "_")) 

# mapping variant to index id
genevariant2index <- protein2varlocat_saav %>%
  dplyr::select(index, 
                gene_variant)

# abundance matrix of peptides mapping to variant locations
pep_premat_variantwo6 <- filter(pep_prematwo6,
                                peptide %in% protein2varlocat_saav$peptide) %>%
  left_join(., genevariant2index) %>%
  relocate(gene_variant) %>%
  mutate(duplicated_variant = duplicated(gene_variant)) %>%
  relocate(duplicated_variant,
           .after = peptide)

# summarize SAAVs (gene+variant location)
# i.e. take the median abundance of all the peptides mapping to the same SAAV

pep_premat_sumgenarwo6 <- pep_premat_variantwo6 %>%
  group_by(gene_variant) %>%
  summarise_at(vars(starts_with("x")), 
               .funs = function(x){median(x, 
                                          na.rm = TRUE)}) %>%
  ungroup() 

full_info_sumgenvar <- pep_premat_sumgenarwo6 %>%
  left_join(., 
            protein2varlocat_saav) %>%
  filter(variant_type == "saav")

pep_matvariantswo6or <- pep_premat_sumgenarwo6 %>%
  column_to_rownames("gene_variant") %>%
  as.matrix()

pep_matvariantswo6 <- pep_premat_sumgenarwo6 %>%
  filter(str_detect(gene_variant, 
                    "CAP1_", 
                    negate = TRUE)) %>%
  column_to_rownames("gene_variant") %>%
  as.matrix()

pep_matvariantswo6nona <- pep_matvariantswo6 %>%
  na.omit()
```

### Generate summarized results tables

``` r
write_tsv(protein2varlocat_saav,
          here("results/proteogenomics/cellreps_2021/all_peptides_mapping_to_saavs_gbm.tsv"))

write_tsv(pep_premat_sumgenarwo6,
          here("results/proteogenomics/cellreps_2021/summarized_median_abundance_of_saavs.tsv"))

write_tsv(full_info_sumgenvar,
          here("results/proteogenomics/cellreps_2021/summarized_median_abundance_of_saavs_full.tsv"))
```

The variant-containing peptides identified and accurately quantified and
summarized are able to account for 214 single amino acid mutations
(SAAVs).

Among these identified SAAVs, 30 were identified in all mixtures and
samples.

## Heatmap of all SAAVs observed

``` r
library(circlize)
```

``` r
library(ComplexHeatmap)
```

``` r
split_patstage <- str_split_fixed(string = colnames(pep_matvariantswo6nona), 
                                  pattern = "[:digit:]+", 
                                  n = 2) 

patient <- str_extract(string = colnames(pep_matvariantswo6nona),
                       pattern = "[:digit:]+")

patientstage <- colnames(pep_matvariantswo6nona) %>%
  str_replace(., 
              pattern = "prim", 
              replacement = "Init") %>%
  str_replace(., 
              pattern = "rec", 
              replacement = "Rec")

stage <- split_patstage[,2] %>%
  if_else(condition = . == "prim",
          true = "Initial",
          false = "Recurrent")

pep_matvariantswo6or_hm <- pep_matvariantswo6or
pep_matvariantswo6nona_hm <- pep_matvariantswo6nona

colnames(pep_matvariantswo6or_hm) <- patientstage
colnames(pep_matvariantswo6nona_hm) <- patientstage
```

``` r
col_fun_or = colorRamp2(c(min(pep_matvariantswo6or_hm, 
                              na.rm = TRUE), 
                       median(pep_matvariantswo6or_hm, 
                              na.rm = TRUE), 
                       max(pep_matvariantswo6or_hm, 
                           na.rm = TRUE)), 
                     c("#2a9d8f", 
                       "white", 
                       "red"))

Heatmap(pep_matvariantswo6or_hm, 
        col = col_fun_or,
        column_split = stage, 
        row_names_gp = gpar(fontsize = 5),
        column_names_gp = gpar(fontsize = 5),
        cluster_rows = FALSE, cluster_columns = FALSE,
        heatmap_legend_param = list(
        title = "log2-\nNormalized\nAbundance"
    ))
```

![](gbm_recurrence_proteogenomics_manuscript_report_files/figure-gfm/unnamed-chunk-22-1.png)

All identified SAAVs

``` r
col_fun_or = colorRamp2(c(min(pep_matvariantswo6or_hm, 
                              na.rm = TRUE), 
                       median(pep_matvariantswo6or_hm, 
                              na.rm = TRUE), 
                       max(pep_matvariantswo6or_hm, 
                           na.rm = TRUE)), 
                     c("#2a9d8f", 
                       "white", 
                       "red"))

Heatmap(pep_matvariantswo6or_hm, 
        col = col_fun_or,
        column_split = patient, 
        row_names_gp = gpar(fontsize = 5.5),
        column_names_gp = gpar(fontsize = 7),
        cluster_rows = FALSE, cluster_columns = FALSE,
        heatmap_legend_param = list(
        title = "log2-\nNormalized\nAbundance"
    ))
```

![](gbm_recurrence_proteogenomics_manuscript_report_files/figure-gfm/unnamed-chunk-23-1.png)

## Heatmap of SAAVs observed on every sample

``` r
col_fun_sel = colorRamp2(c(min(pep_matvariantswo6or_hm, 
                               na.rm = TRUE), 
                       median(pep_matvariantswo6or_hm, 
                              na.rm = TRUE), 
                       max(pep_matvariantswo6or_hm, 
                           na.rm = TRUE)), 
                     c("#2a9d8f", 
                       "white", 
                       "red"))


Heatmap(pep_matvariantswo6nona_hm, 
        col = col_fun_sel,
        column_split = stage, 
        row_names_gp = gpar(fontsize = 7),
        column_names_gp = gpar(fontsize = 7),
        heatmap_legend_param = list(
        title = "log2-\nNormalized\nAbundance"
    ))
```

![](gbm_recurrence_proteogenomics_manuscript_report_files/figure-gfm/unnamed-chunk-24-1.png)

# Limma on SAAVs

## Prep design matrix

``` r
design_limmawo6_saav <- model.matrix(~patient+stage)

sample_annotationwo6 <- sample_annotationwo6 %>%
  mutate(paired_id_semi = paste(patient, 
                                recurrence, 
                                sep = ""))

rownames(design_limmawo6_saav) <- sample_annotationwo6$paired_id_semi
```

## Prep abundance matrices

``` r
colnames(pep_matvariantswo6nona)
```

     [1] "x1prim"  "x1rec"   "x2prim"  "x2rec"   "x3prim"  "x3rec"   "x4prim" 
     [8] "x4rec"   "x5prim"  "x5rec"   "x7prim"  "x7rec"   "x8prim"  "x8rec"  
    [15] "x9prim"  "x9rec"   "x10prim" "x10rec"  "x11prim" "x11rec" 

``` r
rownames(design_limmawo6_saav)
```

     [1] "x1prim"  "x1rec"   "x2prim"  "x2rec"   "x3prim"  "x3rec"   "x4prim" 
     [8] "x4rec"   "x5prim"  "x5rec"   "x7prim"  "x7rec"   "x8prim"  "x8rec"  
    [15] "x9prim"  "x9rec"   "x10prim" "x10rec"  "x11prim" "x11rec" 

## SAAV-level model fit

``` r
limmafit <- lmFit(pep_matvariantswo6nona, 
                  design = design_limmawo6_saav, 
                  method = 'robust')
  
limmafit <- eBayes(limmafit)
  
fit_mod_saavs <- topTable(limmafit, 
                        coef = "stageRecurrent", 
                        number = Inf, 
                        adjust.method = "BH") %>%
  mutate(Protein = rownames(.),
         Limma = "Robust - w Patient effect")

compar_tab_saavs <-  fit_mod_saavs %>% 
  rownames_to_column("gene_variant") 
```

``` r
increased_recwo6_saavs <- compar_tab_saavs %>%
                    filter(logFC > 0,
                           adj.P.Val < 0.05) %>% 
  pull(gene_variant)

decreased_recwo6_saavs <- compar_tab_saavs %>%
                    filter(logFC < 0,
                           adj.P.Val < 0.05) %>% 
  pull(gene_variant)
```

``` r
length(increased_recwo6_saavs)
```

    [1] 1

``` r
length(decreased_recwo6_saavs)
```

    [1] 0

Only 1 SAAV is observed as differentially abundant between Recurrent and
Initial tumors.

## Volcano SAAVs

``` r
diff_saavs <- compar_tab_saavs %>%
                    filter(adj.P.Val < 0.05)

size <- 1.5

volcano_saavs <- ggplot(data = compar_tab_saavs,
                      mapping = aes(x = logFC, 
                                    y = -log10(adj.P.Val))) +
      geom_point(data = compar_tab_saavs %>% 
                   filter(logFC > 0,
                          adj.P.Val < 0.05),
                 mapping = aes(x = logFC, 
                               y = -log10(adj.P.Val)), 
                 color = "red",
                 size = size)+
      geom_point(data = compar_tab_saavs %>% 
                   filter(logFC < -0,
                          adj.P.Val < 0.05),
                 mapping = aes(x = logFC, 
                               y = -log10(adj.P.Val)), 
                 color = "red",
                 size = size) +
      geom_point(data = compar_tab_saavs %>% 
                   filter(logFC > 0,
                          adj.P.Val > 0.05),
                 mapping = aes(x = logFC, 
                               y = -log10(adj.P.Val)), 
                 color = "#2a9d8f",
                 size = size) +
      geom_point(data = compar_tab_saavs %>% 
                   filter(logFC < -0,
                          adj.P.Val > 0.05),
                 mapping = aes(x = logFC, 
                               y = -log10(adj.P.Val)), 
                 color = "#2a9d8f",
                 size = size) +
      geom_hline(yintercept = -log10(0.05),
                 color = "red", 
                 linetype = "dashed") +
      ggrepel::geom_text_repel(data = diff_saavs,
                               aes(label = gene_variant), 
                               size = 4) +
      xlab("logFC - Recurrent / Primary") + 
      theme(axis.text.x = element_text(hjust = 0.5, 
                                       vjust = 0, 
                                       size = 6, 
                                       angle = 360),
            axis.text.y = element_text(hjust = 0.95, 
                                       vjust = 0.2, 
                                       size = 8),
            panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.border = element_rect(colour = "black", 
                                        fill = NA, 
                                        size = 0.5),
            axis.title = element_text(size = 8),
            legend.text = element_text(size = 6),
            legend.title = element_text(size = 8),
            legend.key.height = unit(3, 
                                    'mm'),
            legend.key.width = unit(3, 
                                   'mm'),
            legend.position = "bottom") 
```

``` r
volcano_saavs
```

![](gbm_recurrence_proteogenomics_manuscript_report_files/figure-gfm/unnamed-chunk-33-1.png)

## Boxplot of differentially abundant SAAV

``` r
prelongsaav <- pep_matvariantswo6nona %>%
  as.data.frame() %>%
  rownames_to_column("gene_variant") %>%
   pivot_longer(cols = starts_with("x"),
                         names_to = c("patient"),
                         values_to = "Abundance") %>%
  mutate(Stage = str_extract(patient, 
                             "[^0-9]+$")) %>%
  mutate(Patient = str_remove_all(patient, 
                                  "[^0-9]+$")) %>%
  filter(gene_variant == increased_recwo6_saavs)
```

``` r
diff_abund_saav <- ggplot(prelongsaav, 
                    aes(x = Stage, 
                        y = Abundance, 
                        fill = Stage, 
                        cex.axis = 1.5)) +
  geom_boxplot() +
  geom_dotplot(binaxis = "y", 
               stackdir = "center", 
               dotsize = 0.2) +
  # Box plot with jittered points
  # 0.2 : degree of jitter in x direction
  # geom_jitter(shape=16, position=position_jitter(0.2))
  ylab("log2-Normalized Abundance") +
  geom_signif(
    comparisons = list(c("Initial", 
                         "Recurrence")),
    map_signif_level = TRUE
  ) + 
  stat_compare_means(method="t.test") +
  ggtitle(increased_recwo6_saavs) +
  theme(axis.text.x = element_text(hjust = 0.5, 
                                   vjust = 0, 
                                   size = 6, 
                                   angle = 360),
        axis.text.y = element_text(hjust = 0.95, 
                                   vjust = 0.2, 
                                   size = 8),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "black", 
                                    fill = NA, 
                                    size = 0.5),
        axis.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 8),
        legend.key.height = unit(3, 
                                'mm'),
        legend.key.width = unit(3, 
                               'mm'),
        legend.position = "bottom") 
```

## Boxplot DA SAAV

``` r
  slim_data2saav <- pivot_wider(prelongsaav %>% dplyr::select(-patient), 
                          values_from = Abundance, 
                          names_from = Stage) 

  slim_data3saav <- slim_data2saav %>% 
    mutate(Abs_diff = rec-prim) %>%
    mutate(rank = row_number(Abs_diff)) %>% 
    dplyr::rename(Initial = prim, 
                  Recurrent = rec) %>%
    ungroup() 
```

``` r
ggpaired(slim_data3saav,
         cond1 = "Initial",
         cond2 = "Recurrent", 
         fill = "condition",
         ylab = "log2-Normalized Abundance",
         label = NULL,
         repel = TRUE, 
         point.size = 0.2, 
         line.size = 0.1) +
      theme(axis.text.x = element_text(hjust = 0.5, 
                                       vjust = 0,
                                       size = 5, 
                                       angle = 360),
            axis.text.y = element_text(hjust = 0.95, 
                                       vjust = 0.2, 
                                       size = 5),
            panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.border = element_rect(colour = "black", 
                                        fill = NA, 
                                        size = 0.5),
            axis.title = element_text(size = 8),
            legend.text = element_text(size = 6),
            legend.title = element_text(size = 8),
            legend.key.height= unit(3, 
                                    'mm'),
            legend.key.width= unit(3, 
                                   'mm'),
            legend.position = "bottom",
            strip.text = element_text(hjust = 0.5, 
                                      vjust = 0, 
                                      size = 8, 
                                      angle = 360)) + 
  ggtitle(increased_recwo6_saavs)
```

![](gbm_recurrence_proteogenomics_manuscript_report_files/figure-gfm/unnamed-chunk-38-1.png)

``` r
paired_box_saavs <- ggpaired(slim_data3saav,
         cond1 = "Initial",
         cond2 = "Recurrent", 
         fill = "condition",
         ylab = "log2-Normalized Abundance",
         label = NULL,
         repel = TRUE, 
         point.size = 0.2, 
         line.size = 0.1) +
      theme(axis.text.x = element_text(hjust = 0.5, 
                                       vjust = 0, 
                                       size = 5, 
                                       angle = 360),
            axis.text.y = element_text(hjust = 0.95, 
                                       vjust = 0.2, 
                                       size = 5),
            panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.border = element_rect(colour = "black", 
                                        fill = NA, 
                                        size = 0.5),
            axis.title = element_text(size = 8),
            legend.text = element_text(size = 6),
            legend.title = element_text(size = 8),
            legend.key.height = unit(3, 
                                    'mm'),
            legend.key.width = unit(3, 
                                   'mm'),
            legend.position = "bottom",
            strip.text = element_text(hjust = 0.5, 
                                      vjust = 0, 
                                      size = 8, 
                                      angle = 360)) 
```
