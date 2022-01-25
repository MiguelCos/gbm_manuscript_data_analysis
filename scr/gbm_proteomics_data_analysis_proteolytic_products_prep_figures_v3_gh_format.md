GBM Proteomics and Proteolytic processing data analysis - Reproducible
report
================
Miguel Cosenza
25 January, 2022

``` r
knitr::opts_chunk$set(echo = TRUE, 
                      message = FALSE, 
                      warning = FALSE)
```

``` r
## Required packages ----
library(tidyverse)
library(mixOmics)
library(fs)
library(kableExtra)
library(sva)
library(limma)
library(naniar)
library(missForest)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)
library(DT)
library(here)
library(janitor)
library(drawProteins)
library(seqinr)
library(ggpubr)

source(here("scr/helper_functions.R"))

theme_set(theme(axis.text.x = element_text(hjust = 0.5, vjust = 0, size = 6, 
                                           angle = 90),
                axis.text.y = element_text(hjust = 0.5, vjust = 0, size = 6),
                panel.background = element_blank(),
                panel.grid.major = element_line(color = "grey"),
                panel.border = element_rect(colour = "black", fill=NA, size=1.5),
                axis.title=element_text(size=8),
                legend.text = element_text(size = 7),
                legend.title = element_text(size = 8),
                legend.key.height= unit(3, 'mm'),
                legend.key.width= unit(3, 'mm'),
                legend.position="bottom"))
```

The following report holds the code, results and interpretation of the
proteomics data analysis of Initial-Recurrent glioblastoma samples of
matched patients.

# Sample annotation info

``` r
sample_annotation <- read_csv(here("data/sample_annotation.csv"))

# correct annotation

sample_annotation2 <- sample_annotation %>%
  mutate(patient = paste("x",patient, sep = ""),
         recurrence = case_when(recurrence == "initial" ~ "prim",
                                recurrence == "recurrent" ~ "rec",
                                TRUE ~ recurrence)) %>%
  mutate(paired_id = paste(patient, recurrence, sep = "_")) %>%
  filter(recurrence %in% c("prim", "rec"))

sample_annotationwo6 <- sample_annotation2 %>%
  filter(patient != "x6")
```

## Sample annotation

### Brief note on sample annotation:

Twenty-two patient matched samples (corresponding to 11 patients) were
processed for protein extraction, TMT-labelled and submitted via high-pH
fractionation and further analyzed via LC-MS/MS.

<table class="table table-striped table-hover" style="font-size: 14px; margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
sample_id
</th>
<th style="text-align:left;">
patient
</th>
<th style="text-align:left;">
recurrence
</th>
<th style="text-align:left;">
channel
</th>
<th style="text-align:left;">
mixture
</th>
<th style="text-align:right;">
irt_microl
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
empty
</td>
<td style="text-align:left;">
empty
</td>
<td style="text-align:left;">
empty
</td>
<td style="text-align:left;">
126
</td>
<td style="text-align:left;">
tmt_1
</td>
<td style="text-align:right;">
4
</td>
</tr>
<tr>
<td style="text-align:left;">
2013/210
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
initial
</td>
<td style="text-align:left;">
127N
</td>
<td style="text-align:left;">
tmt_1
</td>
<td style="text-align:right;">
4
</td>
</tr>
<tr>
<td style="text-align:left;">
2014/79
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
recurrent
</td>
<td style="text-align:left;">
127C
</td>
<td style="text-align:left;">
tmt_1
</td>
<td style="text-align:right;">
4
</td>
</tr>
<tr>
<td style="text-align:left;">
2015/18
</td>
<td style="text-align:left;">
2
</td>
<td style="text-align:left;">
initial
</td>
<td style="text-align:left;">
128N
</td>
<td style="text-align:left;">
tmt_1
</td>
<td style="text-align:right;">
6
</td>
</tr>
<tr>
<td style="text-align:left;">
2015/163
</td>
<td style="text-align:left;">
2
</td>
<td style="text-align:left;">
recurrent
</td>
<td style="text-align:left;">
128C
</td>
<td style="text-align:left;">
tmt_1
</td>
<td style="text-align:right;">
6
</td>
</tr>
<tr>
<td style="text-align:left;">
2012/49
</td>
<td style="text-align:left;">
3
</td>
<td style="text-align:left;">
initial
</td>
<td style="text-align:left;">
129N
</td>
<td style="text-align:left;">
tmt_1
</td>
<td style="text-align:right;">
6
</td>
</tr>
<tr>
<td style="text-align:left;">
2013/145
</td>
<td style="text-align:left;">
3
</td>
<td style="text-align:left;">
recurrent
</td>
<td style="text-align:left;">
129C
</td>
<td style="text-align:left;">
tmt_1
</td>
<td style="text-align:right;">
8
</td>
</tr>
<tr>
<td style="text-align:left;">
2013/84
</td>
<td style="text-align:left;">
4
</td>
<td style="text-align:left;">
initial
</td>
<td style="text-align:left;">
130N
</td>
<td style="text-align:left;">
tmt_1
</td>
<td style="text-align:right;">
8
</td>
</tr>
<tr>
<td style="text-align:left;">
2016/154
</td>
<td style="text-align:left;">
4
</td>
<td style="text-align:left;">
recurrent
</td>
<td style="text-align:left;">
130C
</td>
<td style="text-align:left;">
tmt_1
</td>
<td style="text-align:right;">
8
</td>
</tr>
<tr>
<td style="text-align:left;">
pool
</td>
<td style="text-align:left;">
pool
</td>
<td style="text-align:left;">
pool
</td>
<td style="text-align:left;">
131
</td>
<td style="text-align:left;">
tmt_1
</td>
<td style="text-align:right;">
8
</td>
</tr>
<tr>
<td style="text-align:left;">
pool
</td>
<td style="text-align:left;">
pool
</td>
<td style="text-align:left;">
pool
</td>
<td style="text-align:left;">
131C
</td>
<td style="text-align:left;">
tmt_1
</td>
<td style="text-align:right;">
8
</td>
</tr>
<tr>
<td style="text-align:left;">
pool
</td>
<td style="text-align:left;">
pool
</td>
<td style="text-align:left;">
pool
</td>
<td style="text-align:left;">
126
</td>
<td style="text-align:left;">
tmt_2
</td>
<td style="text-align:right;">
4
</td>
</tr>
<tr>
<td style="text-align:left;">
empty
</td>
<td style="text-align:left;">
empty
</td>
<td style="text-align:left;">
empty
</td>
<td style="text-align:left;">
127N
</td>
<td style="text-align:left;">
tmt_2
</td>
<td style="text-align:right;">
4
</td>
</tr>
<tr>
<td style="text-align:left;">
2015/139
</td>
<td style="text-align:left;">
5
</td>
<td style="text-align:left;">
initial
</td>
<td style="text-align:left;">
127C
</td>
<td style="text-align:left;">
tmt_2
</td>
<td style="text-align:right;">
4
</td>
</tr>
<tr>
<td style="text-align:left;">
2016/183
</td>
<td style="text-align:left;">
5
</td>
<td style="text-align:left;">
recurrent
</td>
<td style="text-align:left;">
128N
</td>
<td style="text-align:left;">
tmt_2
</td>
<td style="text-align:right;">
6
</td>
</tr>
<tr>
<td style="text-align:left;">
2015/150
</td>
<td style="text-align:left;">
6
</td>
<td style="text-align:left;">
initial
</td>
<td style="text-align:left;">
128C
</td>
<td style="text-align:left;">
tmt_2
</td>
<td style="text-align:right;">
6
</td>
</tr>
<tr>
<td style="text-align:left;">
2016/178
</td>
<td style="text-align:left;">
6
</td>
<td style="text-align:left;">
recurrent
</td>
<td style="text-align:left;">
129N
</td>
<td style="text-align:left;">
tmt_2
</td>
<td style="text-align:right;">
6
</td>
</tr>
<tr>
<td style="text-align:left;">
2015/161
</td>
<td style="text-align:left;">
7
</td>
<td style="text-align:left;">
initial
</td>
<td style="text-align:left;">
129C
</td>
<td style="text-align:left;">
tmt_2
</td>
<td style="text-align:right;">
8
</td>
</tr>
<tr>
<td style="text-align:left;">
2016/211
</td>
<td style="text-align:left;">
7
</td>
<td style="text-align:left;">
recurrent
</td>
<td style="text-align:left;">
130N
</td>
<td style="text-align:left;">
tmt_2
</td>
<td style="text-align:right;">
8
</td>
</tr>
<tr>
<td style="text-align:left;">
2016/65
</td>
<td style="text-align:left;">
8
</td>
<td style="text-align:left;">
initial
</td>
<td style="text-align:left;">
130C
</td>
<td style="text-align:left;">
tmt_2
</td>
<td style="text-align:right;">
8
</td>
</tr>
<tr>
<td style="text-align:left;">
2016/189
</td>
<td style="text-align:left;">
8
</td>
<td style="text-align:left;">
recurrent
</td>
<td style="text-align:left;">
131
</td>
<td style="text-align:left;">
tmt_2
</td>
<td style="text-align:right;">
8
</td>
</tr>
<tr>
<td style="text-align:left;">
pool
</td>
<td style="text-align:left;">
pool
</td>
<td style="text-align:left;">
pool
</td>
<td style="text-align:left;">
131C
</td>
<td style="text-align:left;">
tmt_2
</td>
<td style="text-align:right;">
8
</td>
</tr>
<tr>
<td style="text-align:left;">
pool
</td>
<td style="text-align:left;">
pool
</td>
<td style="text-align:left;">
pool
</td>
<td style="text-align:left;">
126
</td>
<td style="text-align:left;">
tmt_3
</td>
<td style="text-align:right;">
4
</td>
</tr>
<tr>
<td style="text-align:left;">
pool
</td>
<td style="text-align:left;">
pool
</td>
<td style="text-align:left;">
pool
</td>
<td style="text-align:left;">
127N
</td>
<td style="text-align:left;">
tmt_3
</td>
<td style="text-align:right;">
4
</td>
</tr>
<tr>
<td style="text-align:left;">
empty
</td>
<td style="text-align:left;">
empty
</td>
<td style="text-align:left;">
empty
</td>
<td style="text-align:left;">
127C
</td>
<td style="text-align:left;">
tmt_3
</td>
<td style="text-align:right;">
4
</td>
</tr>
<tr>
<td style="text-align:left;">
2016/70
</td>
<td style="text-align:left;">
9
</td>
<td style="text-align:left;">
initial
</td>
<td style="text-align:left;">
128N
</td>
<td style="text-align:left;">
tmt_3
</td>
<td style="text-align:right;">
6
</td>
</tr>
<tr>
<td style="text-align:left;">
2016/207
</td>
<td style="text-align:left;">
9
</td>
<td style="text-align:left;">
recurrent
</td>
<td style="text-align:left;">
128C
</td>
<td style="text-align:left;">
tmt_3
</td>
<td style="text-align:right;">
6
</td>
</tr>
<tr>
<td style="text-align:left;">
2016/145
</td>
<td style="text-align:left;">
10
</td>
<td style="text-align:left;">
initial
</td>
<td style="text-align:left;">
129N
</td>
<td style="text-align:left;">
tmt_3
</td>
<td style="text-align:right;">
6
</td>
</tr>
<tr>
<td style="text-align:left;">
2016/249
</td>
<td style="text-align:left;">
10
</td>
<td style="text-align:left;">
recurrent
</td>
<td style="text-align:left;">
129C
</td>
<td style="text-align:left;">
tmt_3
</td>
<td style="text-align:right;">
8
</td>
</tr>
<tr>
<td style="text-align:left;">
2016/231
</td>
<td style="text-align:left;">
11
</td>
<td style="text-align:left;">
initial
</td>
<td style="text-align:left;">
130N
</td>
<td style="text-align:left;">
tmt_3
</td>
<td style="text-align:right;">
8
</td>
</tr>
<tr>
<td style="text-align:left;">
2017/50
</td>
<td style="text-align:left;">
11
</td>
<td style="text-align:left;">
recurrent
</td>
<td style="text-align:left;">
130C
</td>
<td style="text-align:left;">
tmt_3
</td>
<td style="text-align:right;">
8
</td>
</tr>
<tr>
<td style="text-align:left;">
empty
</td>
<td style="text-align:left;">
empty
</td>
<td style="text-align:left;">
empty
</td>
<td style="text-align:left;">
131
</td>
<td style="text-align:left;">
tmt_3
</td>
<td style="text-align:right;">
8
</td>
</tr>
<tr>
<td style="text-align:left;">
empty
</td>
<td style="text-align:left;">
empty
</td>
<td style="text-align:left;">
empty
</td>
<td style="text-align:left;">
131C
</td>
<td style="text-align:left;">
tmt_3
</td>
<td style="text-align:right;">
8
</td>
</tr>
</tbody>
</table>

# Brief note on data analysis and processing:

Three proteomics data analysis approaches were used to explore the
spectral data and will be briefly described:

1.  General large scale proteomics.
2.  Proteolytic processing.
3.  Proteogenomics

The FragPipe bioinformatic pipeline was used for each of these
approaches with varying parameters depending on the type of explorative
approach.

**1. General large scale proteomics analyses** : In order to explore the
differential abundance of proteins between Recurrent and Initial tumor,
a fully-tryptic search was performed against the EBI Human canonical
proteome (version 2021_03), appended with common contaminants and iRT
peptides. A minimum percentage of %5 of the total summed reporter ion
intensities were required to consider a peptide and reporter ion for
quantitation.

**2. Proteolytic processing**: The large scale evaluation of
differential proteolytic activity between recurrent vs primary tumor was
executed via a semi-tryptic search against the EBI Human canonical
proteome (version 2021_03), appended with common contaminants and iRT
peptides. A minimum percentage of %5 of the total summed reporter ion
intensities were required to consider a peptide and reporter ion for
quantitation.

**3. Proteogenomics**: The explorative proteogenomics analysis started
with the generation of a GBM-specific database from publicly available
RNA-seq sequencing data from tumor specimens. Fastq files were obtain
from the SRA BioProject PRJNA627121 (Cell Rep, 2021 Mar 2;34(9):108787).
Eight samples containing paired reads were selected, based on: 1) the
higher number of reads per file; 2) having a similar proportion of males
and females and 3) not bearing isocitrate-dehyderogenase (IDH)-1
mutation (IDH WT). None of these patients undergo previous treatment,
therefore potential variants in the sequencing data are not expected to
be associated to treatment-related mutagenic processes. Spectral data
was then searched against this custom database using a fully-tryptic
approac. No minimum intensity threshold was set for reporter ion
intensity for quantitation.

Exceptionally when stated differently within each sub-approach described
above, all these searches were executed using the following parameters:
MSFragger was used as a search engine, setting a fully or semi-tryptic
specificity allowing for 1 missed cleavages. Precursor mass tolerance
was set to -20/20 and fragment mass tolerance of 20 ppm with mass
calibration and parameter optimization resulting in 5 ppm. Peptide
N-terminal acetylation and peptide N-terminal TMT labeling were set as
variable modifications. TMT labelling at K and carbamidomethylation of C
were set as fixed modifications. MSBooster was used for deep-learning
based predictions of retention time and spectra. Predicted features were
used by Percolator for post-processing scoring and false discovery rate
(FDR) control via target-decoy competition of all peptide-to-spectrum
matches (PSMs) obtained from MSFragger search. ProteinProphet was used
as a protein inference algorithm keeping all PSMs with probability score
bigger than 0.9. Final report and FDR estimation was based on the
filtered PSM and protein lists before protein scoring.

Relative quantitation of identified peptides within each sample was
performed via their reported ion intensities using TMT-integrator. Only
PSMs coming from unique peptides with a minimum probability of 0.9 and
an isotopic purity of at least 50% were considered for quantitation.
Quantitative values were normalized via median centering and summarized
by protein/gene and peptide using a virtual reference channel, to
finally generate protein and peptide normalized abundance matrices of
abundance used for further statistical processing in combination with
the sample annotation.

# General proteomics analysis

## Initial data loading and wrangling

Load id/quant data and prep coverage information

``` r
## tmt-integrator output loading 
tmt_protdata_frag <- read_tsv(here("data/specific_search_fragpipe17/specific_no_ptms_2/tmt-report/abundance_protein_MD.tsv")) %>% 
  janitor::clean_names()

identified_proteins <- read_tsv(here("data/specific_search_fragpipe17/specific_no_ptms_2/combined_protein.tsv")) %>% 
  janitor::clean_names()

identified_peptides <- read_tsv(here("data/specific_search_fragpipe17/specific_no_ptms_2/combined_peptide.tsv")) %>%
  janitor::clean_names()

fasta_ident_prots <- read.fasta(here("data/specific_search_fragpipe17/specific_no_ptms_2/mix_2/protein.fas"), seqtype = "AA", as.string = TRUE)
```

Wrangle and merge protein ID/quant with sample annotation

``` r
# wide expression matrix with NAs
expr_matrix <- dplyr::select(tmt_protdata_frag,
                             Protein = index, starts_with("x"))

long_mat <- pivot_longer(expr_matrix,
                         cols = starts_with("x"),
                         names_sep = "_",
                         names_to = c("patient", "recurrence"),
                         values_to = "Abundance")

long_mat_2 <- long_mat %>%
  mutate(paired_id = paste(patient, recurrence, sep = "_"))

paired_annot_uniq <- long_mat_2 %>%
  dplyr::select(patient, recurrence, paired_id) %>%
  distinct()

# correct annotation

sample_annotation2 <- sample_annotation %>%
  mutate(patient = paste("x",patient, sep = ""),
         recurrence = case_when(recurrence == "initial" ~ "prim",
                                recurrence == "recurrent" ~ "rec",
                                TRUE ~ recurrence)) %>%
  mutate(paired_id = paste(patient, recurrence, sep = "_")) %>%
  filter(recurrence %in% c("prim", "rec"))

## Annotated abundance data in long format ----

quant_annotated <- left_join(long_mat_2,
                             sample_annotation2, 
                             by = c("paired_id", "patient","recurrence")) %>%
  mutate(Channel_mix = paste(mixture,channel,  sep = "_"))

## Wide xpression matrices without NAs
expr_matnona <- column_to_rownames(expr_matrix,
                            "Protein") %>%
          as.matrix() %>%
          na.omit()
```

## ID coverage

### Coverage plot

``` r
ids_summary <- quant_annotated %>% 
  dplyr::select(Protein, mixture, Abundance) %>% 
  distinct(Protein, mixture, .keep_all = TRUE) %>%
  group_by(mixture) %>%
  summarise(`Nr of IDs` = sum(!is.na(Abundance))) %>% 
  bind_rows(.,
            tibble(mixture = "In all three",
                   `Nr of IDs` = nrow(expr_matnona))) %>%
  bind_rows(.,
            tibble(mixture = "Overall",
                   `Nr of IDs` = nrow(expr_matrix)))
```

-   Number of identified proteins: 5954

-   Number of identified and quantified proteins in all samples: 3228

The initial analysis with MaxQuant yielded 2090 IDs in all samples
(without missing values) and 5176 overall.

``` r
ggplot(data = ids_summary,
       aes(y = `Nr of IDs`, x = mixture)) +
  geom_col(fill = "#008683") + 
  geom_text(size = 5, position = position_stack(vjust = 0.5), aes(label = `Nr of IDs`)) +
  coord_cartesian(ylim = c(0, 5000))+
  coord_flip() +
  labs(title = "Number of Identified and quantified Proteins by Mixture") +
  ylab("Nr of IDed and Quant Proteins") + 
  xlab("TMT Mixture") + 
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.1, size = 10, angle = 0),
        axis.text.y = element_text(hjust = 0.5, vjust = 0.1, size = 10, angle = 0),
        panel.background = element_blank(),
        panel.grid.major = element_line(color = "grey"),
        panel.border = element_rect(colour = "black", fill=NA, size=1.5),
        axis.title=element_text(size=12,face="bold"))
```

![](gbm_proteomics_data_analysis_proteolytic_products_prep_figures_v3_gh_format_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->
## Exploratory analysis and visualizations

### Visualization of missing values

``` r
vis_miss(expr_matrix) + 
          ggtitle("Distribution of missing values") + 
          theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.1, size = 10, angle = 90),
                axis.text.y = element_text(hjust = 0.5, vjust = 0.1, size = 10),
                panel.background = element_blank(),
                panel.grid.major = element_line(color = "grey"),
                panel.border = element_rect(colour = "black", fill=NA, size=1.5),
                axis.title=element_text(size=12,face="bold"))
```

![](gbm_proteomics_data_analysis_proteolytic_products_prep_figures_v3_gh_format_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

### QC plots after TMT integrator pre-processing and normalization (Protein level)

#### By TMT channel/sample

##### General Abudance distribution by TMT channel/sample

``` r
qcplot <- ggplot(quant_annotated,
                 mapping = aes(x = Channel_mix, 
                               y = Abundance, 
                               fill = recurrence)) +
          geom_boxplot() + 
          #facet_grid(.~recurrence, scales = "free") +
          theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.1, size = 10, angle = 90),
                axis.text.y = element_text(hjust = 0.5, vjust = 0.1, size = 10),
                panel.background = element_blank(),
                panel.grid.major = element_line(color = "grey"),
                panel.border = element_rect(colour = "black", fill=NA, size=1.5),
                axis.title=element_text(size=12,face="bold"))

print(qcplot)
```

![](gbm_proteomics_data_analysis_proteolytic_products_prep_figures_v3_gh_format_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

The median-centering after TMT-integrator seem to show a slight
mixture-dependent centered median per sample. The paired distribution of
initial and recurrent tumor allows the biological effect to be more
important than any potential TMT-mixture-based batch effect.

The different abundance distribution can also be explained by a
different distribution of missing values.

PCA and sPLS-DA analysis allow to confirm this.

#### By Protein (General)

``` r
# select/sample N proteins to observe their abundance distribution  
which_prots <- sample_proteins(x = quant_annotated$Protein,
                               size = 50,
                               seed = 156)
```

50 proteins were randomly sampled to check for the normality of their
‘normalized’ Abundance distribution.

``` r
qcplot_prots <- ggplot(quant_annotated %>% 
                         filter(Protein %in% which_prots),
                 mapping = aes(x = reorder(Protein,Abundance,na.rm = TRUE), 
                               y = Abundance)) +
          geom_boxplot() + 
          geom_jitter(color="black", size=0.4, alpha=0.9) +
          labs(subtitle = "Protein IDs arranged by median abudance",
               x = "Protein ID ") + 
          theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.1, size = 10, angle = 90),
                axis.text.y = element_text(hjust = 0.5, vjust = 0.1, size = 10),
                panel.background = element_blank(),
                panel.grid.major = element_line(color = "grey"),
                panel.border = element_rect(colour = "black", fill=NA, size=1.5),
                axis.title=element_text(size=12,face="bold"))

print(qcplot_prots)
```

![](gbm_proteomics_data_analysis_proteolytic_products_prep_figures_v3_gh_format_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

#### By Protein and Batch

50 proteins were randomly sampled to check for the normality of their
‘normalized’ Abundance distribution.

``` r
qcplot_prots2 <- ggplot(quant_annotated %>% 
                         filter(Protein %in% which_prots),
                 mapping = aes(x = reorder(Protein,Abundance,na.rm = TRUE), 
                               y = Abundance, 
                               fill = mixture)) +
          geom_boxplot() + 
          geom_jitter(color="black", size=0.4, alpha=0.9) +
          labs(subtitle = "Protein IDs arranged by median abudance",
                       x = "Protein ID ") + 
          facet_grid(.~mixture, scales = "free") +
          theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.1, size = 10, angle = 90),
                axis.text.y = element_text(hjust = 0.5, vjust = 0.1, size = 10),
                panel.background = element_blank(),
                panel.grid.major = element_line(color = "grey"),
                panel.border = element_rect(colour = "black", fill=NA, size=1.5),
                axis.title=element_text(size=12,face="bold"))

print(qcplot_prots2)
```

![](gbm_proteomics_data_analysis_proteolytic_products_prep_figures_v3_gh_format_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

The abundance distribution of individual proteins can be considered as
normaly distributed in general, which makes then suitable for
inferential statistics via linear models.

### Exploratory PCA analisys

#### General PCA

##### Percentage of explained variance per component

The PCA was applied on the expression matrix without missing values

``` r
pca_res = pca(t_expr_matnona, ncomp = 10, center = TRUE, scale = TRUE)
```

``` r
plot(pca_res)
```

![](gbm_proteomics_data_analysis_proteolytic_products_prep_figures_v3_gh_format_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

\~45% of variance based on protein abundance can be explained by
components 1 and 2 after PCA.

##### PCA plot on samples based on protein abundance

``` r
# preprocess pca results 
pca_variates <- pca_res$variates$X %>% 
  as.data.frame() %>%
  rownames_to_column("Sample") %>%
  separate(col = "Sample",
           into = c("sample", "recurrence", "mixture"),
           sep = "\\_", extra = "merge") %>%
  mutate(paired_id = paste(sample, recurrence, sep = "_")) %>%
  dplyr::select(-c(mixture, recurrence)) %>%
  left_join(. , sample_annotation2,
            by = "paired_id") %>%
  mutate(Stage = if_else(condition = recurrence == "prim", 
                         true = "Initial", 
                         false = "Recurrent"),
         Mix = case_when(mixture == "tmt_1" ~ "Mix1",
                         mixture == "tmt_2" ~ "Mix2",
                         mixture == "tmt_3" ~ "Mix3"))

pca_variateswopat6 <- pca_variates %>% 
  filter(patient != "x6")
```

``` r
ggplot(data = pca_variates,
       aes(x = PC1, y = PC2)) +
  geom_point(size = 4, aes(shape = mixture, color = recurrence)) +  
  geom_text(aes(label = patient), position = position_nudge(x = 4.5)) +
  labs(title = "PCA plot of samples\nBased on the abundance values of proteins") +
  xlab(paste("PC1", round(pca_res$prop_expl_var$X[1]*100), "% var explained")) + 
  ylab(paste("PC2", round(pca_res$prop_expl_var$X[2]*100), "% var explained")) + 
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0, size = 10, angle = 360),
        axis.text.y = element_text(hjust = 0.95, vjust = 0.2, size = 10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8),
        legend.key.height= unit(3, 'mm'),
        legend.key.width= unit(3, 'mm'),
        legend.position="bottom")
```

![](gbm_proteomics_data_analysis_proteolytic_products_prep_figures_v3_gh_format_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

There is no evident batch effect after PCA analysis.

There is a good number of individual patients that show important
differences in their PCA variates values in recurrent vs primary.

``` r
ggplot(data = pca_variateswopat6,
       aes(x = PC1, y = PC2)) +
  geom_point(size = 4, aes(shape = mixture, color = recurrence)) +  
  geom_text(aes(label = patient), position = position_nudge(x = 4.5)) +
  labs(title = "PCA plot of samples\nBased on the abundance values of proteins (patient 6 excluded)") +
  xlab(paste("PC1", round(pca_res$prop_expl_var$X[1]*100), "% var explained")) + 
  ylab(paste("PC2", round(pca_res$prop_expl_var$X[2]*100), "% var explained")) + 
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0, size = 10, angle = 360),
        axis.text.y = element_text(hjust = 0.95, vjust = 0.2, size = 10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8),
        legend.key.height= unit(3, 'mm'),
        legend.key.width= unit(3, 'mm'),
        legend.position="bottom")
```

![](gbm_proteomics_data_analysis_proteolytic_products_prep_figures_v3_gh_format_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

``` r
pca_plot <- ggplot(data = pca_variateswopat6,
       aes(x = PC1, y = PC2)) +
  geom_point(size = 4, aes(shape = Mix, color = Stage)) +  
  geom_text(aes(label = patient), position = position_nudge(x = 4.5), size = 2) +
  xlab(paste("PC1", round(pca_res$prop_expl_var$X[1]*100), "% var explained")) + 
  ylab(paste("PC2", round(pca_res$prop_expl_var$X[2]*100), "% var explained")) +
theme(axis.text.x = element_text(hjust = 0.5, vjust = 0, size = 7, angle = 360),
      axis.text.y = element_text(hjust = 0.95, vjust = 0.2, size = 7),
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, size=0.5),
      axis.title = element_text(size = 7),
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 8),
      legend.key.height= unit(3, 'mm'),
      legend.key.width= unit(3, 'mm'),
      legend.position="bottom")
```

``` r
ggsave(plot = pca_plot, 
       filename = here::here("figures/pca_proteomics_gbm.tiff"), 
       device = "tiff",
       units = "mm",
       width = 63,
       height = 50)

ggsave(plot = pca_plot, 
       filename = here::here("figures/pca_proteomics_gbm.eps"), 
       device = "eps",
       units = "mm",
       width = 63,
       height = 50)
```

### Exploratory sPLS-DA

``` r
plsda_rest1 <- plsda(X = t_expr_matnona,
                     Y = recurrence,
                     ncomp = 3)
```

``` r
plotIndiv(plsda_rest1,
          nd.names = TRUE, ellipse = TRUE, legend = FALSE, title = "sPLS-DA_Rec-Prim", 
          size.title = rel(2.2), size.xlabel = rel(1.5), size.ylabel = rel(1.5), 
          siz.axis = rel(1.2), point.lwd = 0.8, cex = 4)  
```

![](gbm_proteomics_data_analysis_proteolytic_products_prep_figures_v3_gh_format_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

Patient 6 in recurrence behaves suspiciously similar to initial samples.
We consider this to be an outlier associated to a potential problem
during sample prep.

Therefore we decided to exclude patient 6 from further analyses.

``` r
plsda_rest2 <- plsda(X = t_expr_matnona,
                     Y = mixture,
                     ncomp = 3)
```

``` r
plotIndiv(plsda_rest2,
          nd.names = TRUE, ellipse = TRUE, legend = FALSE, title = "sPLS-DA_Rec-Prim", 
          size.title = rel(2.2), size.xlabel = rel(1.5), size.ylabel = rel(1.5), 
          siz.axis = rel(1.2), point.lwd = 0.8, cex = 4)  
```

![](gbm_proteomics_data_analysis_proteolytic_products_prep_figures_v3_gh_format_files/figure-gfm/unnamed-chunk-25-1.png)<!-- -->

We can say that the biological differences between initial and recurrent
status are stronger than the differences that could come from potential
batch effects.

### Exploratory sPLS-DA without patient 6

``` r
plsda_restwo6 <- plsda(X = t_expr_matnonawo6,
                     Y = recurrencewo6,
                     ncomp = 3)
```

``` r
plotIndiv(plsda_restwo6,
          nd.names = TRUE, ellipse = TRUE, legend = FALSE, title = "sPLS-DA_Rec-Prim wo Pat6", 
          size.title = rel(2.2), size.xlabel = rel(1.5), size.ylabel = rel(1.5), 
          siz.axis = rel(1.2), point.lwd = 0.8, cex = 4)  
```

![](gbm_proteomics_data_analysis_proteolytic_products_prep_figures_v3_gh_format_files/figure-gfm/unnamed-chunk-27-1.png)<!-- -->

#### GGplot-based visualization of Exploratory sPLS-DA

``` r
# preprocess pca results 
plsda_variates <- plsda_restwo6$variates$X %>% 
  as.data.frame() %>%
  rownames_to_column("Sample") %>%
  separate(col = "Sample",
           into = c("sample", "recurrence", "mixture"),
           sep = "\\_", extra = "merge") %>%
  mutate(paired_id = paste(sample, recurrence, sep = "_")) %>%
  dplyr::select(-c(mixture, recurrence)) %>%
  left_join(. , sample_annotation2,
            by = "paired_id") 

plsda_variateswpat6 <- plsda_rest1$variates$X %>% 
  as.data.frame() %>%
  rownames_to_column("Sample") %>%
  separate(col = "Sample",
           into = c("sample", "recurrence", "mixture"),
           sep = "\\_", extra = "merge") %>%
  mutate(paired_id = paste(sample, recurrence, sep = "_")) %>%
  dplyr::select(-c(mixture, recurrence)) %>%
  left_join(. , sample_annotation2,
            by = "paired_id") 
```

``` r
ggplot(data = plsda_variateswpat6,
       aes(x = comp1, y = comp2)) +
  geom_point(size = 4, aes(shape = mixture, color = recurrence)) +  
  geom_text(aes(label = patient), position = position_nudge(x = 2)) +
  labs(title = "sPLS-DA plot of samples\nBased on the abundance values of proteins") +
  xlab(paste("sPLS-DA Comp-1", round(plsda_restwo6$prop_expl_var$X[1]*100), "% var explained")) + 
  ylab(paste("sPLS-DA Comp-2", round(plsda_restwo6$prop_expl_var$X[2]*100), "% var explained")) + 
  ggforce::geom_mark_ellipse(mapping = aes(color = recurrence)) +
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0, size = 10, angle = 360),
        axis.text.y = element_text(hjust = 0.95, vjust = 0.2, size = 10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8),
        legend.key.height= unit(3, 'mm'),
        legend.key.width= unit(3, 'mm'),
        legend.position="bottom")
```

![](gbm_proteomics_data_analysis_proteolytic_products_prep_figures_v3_gh_format_files/figure-gfm/unnamed-chunk-29-1.png)<!-- -->

Patient 6 in recurrence behaves suspiciously similar to initial samples.
We consider this to be an outlier associated to a potential problem
during sample prep.

Therefore we decided to exclude patient 6 from further analyses.

``` r
ggplot(data = plsda_variates,
       aes(x = comp1, y = comp2)) +
  geom_point(size = 4, aes(shape = mixture, color = recurrence)) +  
  geom_text(aes(label = patient), position = position_nudge(x = 2)) +
  labs(title = "sPLS-DA plot of samples - pat6 excluded \nBased on the abundance values of proteins") +
  xlab(paste("sPLS-DA Comp-1", round(plsda_restwo6$prop_expl_var$X[1]*100), "% var explained")) + 
  ylab(paste("sPLS-DA Comp-2", round(plsda_restwo6$prop_expl_var$X[2]*100), "% var explained")) + 
  ggforce::geom_mark_ellipse(mapping = aes(color = recurrence)) +
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0, size = 10, angle = 360),
        axis.text.y = element_text(hjust = 0.95, vjust = 0.2, size = 10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8),
        legend.key.height= unit(3, 'mm'),
        legend.key.width= unit(3, 'mm'),
        legend.position="bottom")
```

![](gbm_proteomics_data_analysis_proteolytic_products_prep_figures_v3_gh_format_files/figure-gfm/unnamed-chunk-30-1.png)<!-- -->

The difference between initial and recurrent tumor is better resolved in
the sPLS-DA plot when excluding patient 6.

## Sparcity reduction

Only proteins that are identified and quantified in 2 out of 3 TMT
mixtures will be kept for the inferential analysis. The proteins missing
in 1 mixture will be imputed with `missForest`.

**Function to select proteins to exclude:**

``` r
sel_proteins_missing <- function(long_matrix,
                                 threshold) {

na_count <- group_by(long_matrix,
                     Protein, mixture) %>%
    summarise(na_count = sum(is.na(Abundance)),
              total = n()) %>% 
    ungroup() %>% 
    mutate(NA_fraction = na_count/total)

na_count_perbatch <- na_count %>%
  group_by(Protein) %>%
  summarise(na_per_batch = sum(NA_fraction)) %>%
  ungroup()

proteins2exclude <- na_count_perbatch %>%
  filter(na_per_batch > threshold) %>%
  pull(Protein)

return(proteins2exclude)
}
```

**Excluded proteins:**

``` r
proteins2exclude <- sel_proteins_missing(quant_annotated,
                                         threshold = 1)
```

Exclude proteins from initial abundance matrix.

``` r
expr_matrix_filt <- expr_matrix2 %>%
  filter(!Protein %in% proteins2exclude)
```

### Visualize missing values (before and after filtering)

``` r
nofiltnas_plot <- vis_miss(expr_matrix) + 
          ggtitle("Distribution of missing values",
                  subtitle = paste("All quant proteins. Total =", nrow(expr_matrix))) + 
          theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.1, size = 10, angle = 90),
                axis.text.y = element_text(hjust = 0.5, vjust = 0.1, size = 10),
                panel.background = element_blank(),
                panel.grid.major = element_line(color = "grey"),
                panel.border = element_rect(colour = "black", fill=NA, size=1.5),
                axis.title=element_text(size=12,face="bold"))

filtnas_plot <- vis_miss(expr_matrix_filt) + 
          ggtitle("Distribution of missing values",
                  subtitle = paste("Only present in 2/3 TMT mixtures. Total =", nrow(expr_matrix_filt))) + 
          theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.1, size = 10, angle = 90),
                axis.text.y = element_text(hjust = 0.5, vjust = 0.1, size = 10),
                panel.background = element_blank(),
                panel.grid.major = element_line(color = "grey"),
                panel.border = element_rect(colour = "black", fill=NA, size=1.5),
                axis.title=element_text(size=12,face="bold"))
```

``` r
cowplot::plot_grid(nofiltnas_plot, filtnas_plot, ncol = 2)
```

![](gbm_proteomics_data_analysis_proteolytic_products_prep_figures_v3_gh_format_files/figure-gfm/unnamed-chunk-35-1.png)<!-- -->

## `missForest` imputation:

``` r
mat_filt <- expr_matrix_filt %>%
  column_to_rownames("Protein") %>%
  as.matrix()

t_mat_filt <- t(mat_filt)

if(!file.exists(here("rds/missforest_imp_mat_filt.Rds"))){
  
 t_imp_mat_filt <- missForest::missForest(t_mat_filt)
 t_imp_filt_mat <- t_imp_mat_filt$ximp
 
 mat_filt_imp <- t(t_imp_filt_mat)
 
 saveRDS(t_imp_mat_filt,
         file = here("rds/missforest_imp_mat_filt.Rds"))
 
} else {
  
 t_imp_mat_filt <- readRDS(here("rds/missforest_imp_mat_filt.Rds"))
 t_imp_filt_mat <- t_imp_mat_filt$ximp
 
 mat_filt_imp <- t(t_imp_filt_mat)
  
}
```

## Inferential analysis

**Protein to Gene annotation**

``` r
prot2gene <- tmt_protdata_frag %>%
  dplyr::select(Protein = index, Gene = gene)
```

**Prep protein abundance matrix**

``` r
mat_filt_impwo6 <- as.data.frame(mat_filt_imp) %>%
  rownames_to_column("Protein") %>%
  dplyr::select(-starts_with("x6")) %>%
  column_to_rownames("Protein") %>%
  as.matrix()

sample_annotationwo6 <- sample_annotation2 %>%
  filter(patient != "x6")
```

### Set up design matrix

Initial vs recurrent samples are patient-matched, therefore the design
matrix is setup so the linear model can account for the patient random
effect.

``` r
design_wo6 <- model.matrix(~patientwo6+recurrencewo6)

rownames(design_wo6) <- sample_annotationwo6$paired_id
```

Check design matrix:

``` r
print(design_wo6)
```

    ##          (Intercept) patientwo6x10 patientwo6x11 patientwo6x2 patientwo6x3
    ## x1_prim            1             0             0            0            0
    ## x1_rec             1             0             0            0            0
    ## x2_prim            1             0             0            1            0
    ## x2_rec             1             0             0            1            0
    ## x3_prim            1             0             0            0            1
    ## x3_rec             1             0             0            0            1
    ## x4_prim            1             0             0            0            0
    ## x4_rec             1             0             0            0            0
    ## x5_prim            1             0             0            0            0
    ## x5_rec             1             0             0            0            0
    ## x7_prim            1             0             0            0            0
    ## x7_rec             1             0             0            0            0
    ## x8_prim            1             0             0            0            0
    ## x8_rec             1             0             0            0            0
    ## x9_prim            1             0             0            0            0
    ## x9_rec             1             0             0            0            0
    ## x10_prim           1             1             0            0            0
    ## x10_rec            1             1             0            0            0
    ## x11_prim           1             0             1            0            0
    ## x11_rec            1             0             1            0            0
    ##          patientwo6x4 patientwo6x5 patientwo6x7 patientwo6x8 patientwo6x9
    ## x1_prim             0            0            0            0            0
    ## x1_rec              0            0            0            0            0
    ## x2_prim             0            0            0            0            0
    ## x2_rec              0            0            0            0            0
    ## x3_prim             0            0            0            0            0
    ## x3_rec              0            0            0            0            0
    ## x4_prim             1            0            0            0            0
    ## x4_rec              1            0            0            0            0
    ## x5_prim             0            1            0            0            0
    ## x5_rec              0            1            0            0            0
    ## x7_prim             0            0            1            0            0
    ## x7_rec              0            0            1            0            0
    ## x8_prim             0            0            0            1            0
    ## x8_rec              0            0            0            1            0
    ## x9_prim             0            0            0            0            1
    ## x9_rec              0            0            0            0            1
    ## x10_prim            0            0            0            0            0
    ## x10_rec             0            0            0            0            0
    ## x11_prim            0            0            0            0            0
    ## x11_rec             0            0            0            0            0
    ##          recurrencewo6rec
    ## x1_prim                 0
    ## x1_rec                  1
    ## x2_prim                 0
    ## x2_rec                  1
    ## x3_prim                 0
    ## x3_rec                  1
    ## x4_prim                 0
    ## x4_rec                  1
    ## x5_prim                 0
    ## x5_rec                  1
    ## x7_prim                 0
    ## x7_rec                  1
    ## x8_prim                 0
    ## x8_rec                  1
    ## x9_prim                 0
    ## x9_rec                  1
    ## x10_prim                0
    ## x10_rec                 1
    ## x11_prim                0
    ## x11_rec                 1
    ## attr(,"assign")
    ##  [1] 0 1 1 1 1 1 1 1 1 1 2
    ## attr(,"contrasts")
    ## attr(,"contrasts")$patientwo6
    ## [1] "contr.treatment"
    ## 
    ## attr(,"contrasts")$recurrencewo6
    ## [1] "contr.treatment"

``` r
mat_filt_impwo6[1:3,1:12]
```

    ##             x1_prim   x1_rec  x2_prim   x2_rec  x3_prim   x3_rec  x4_prim
    ## A0A075B6H9 19.94963 18.06669 18.24085 18.60434 19.40524 18.81092 18.01550
    ## A0A0A0MS15 18.53463 18.59081 17.84463 18.06341 19.35137 19.16811 17.79523
    ## A0A0B4J1U7 22.16874 20.92584 20.55509 21.41375 22.41883 21.81234 21.02839
    ##              x4_rec  x5_prim   x5_rec  x7_prim   x7_rec
    ## A0A075B6H9 18.47586 19.19619 18.89962 18.39509 18.95522
    ## A0A0A0MS15 18.21430 18.58008 17.37422 18.04101 19.07056
    ## A0A0B4J1U7 21.48788 22.84955 20.74536 21.02607 21.91741

### Fit limma

``` r
source(here("scr/fit_limmawo6.R"))
```

``` r
# fit limma
limma_tab_wo6 <- fit_limmawo6(mat_filt_impwo6, 
                         design_wo6, 
                         method = 'robust', 
                         Limma = "Robust - w Patient effect wo Pat6",
                         prot2gene = prot2gene)

#get proteins increased in recurrence
increased_in_rec <- limma_tab_wo6 %>%
  filter(logFC > 0,
         adj.P.Val < 0.05) %>%
  pull(Protein)

#get proteins decreased in recurrence
decreased_in_rec <- limma_tab_wo6 %>%
  filter(logFC < 0,
         adj.P.Val < 0.05) %>%
  pull(Protein)
```

### Volcano plots

``` r
our_volcano(limma_tab_wo6, 
            FC_cutoff = 0, 
            pval_cutoff = 0.05, 
            color_diffex = "red", 
            color_nondifex = "#2a9d8f", 
            interesting_proteins = NULL, 
            vert_line_col = "red",
            hline_col = "red", 
            hline_pos = 0.05, 
            linetype = "dashed",
            increased_in = "Recurrent", 
            comparison_title = "Recurrent vs Initial") +
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0, size = 10, angle = 360),
        axis.text.y = element_text(hjust = 0.95, vjust = 0.2, size = 10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8),
        legend.key.height= unit(3, 'mm'),
        legend.key.width= unit(3, 'mm'),
        legend.position="bottom") 
```

![](gbm_proteomics_data_analysis_proteolytic_products_prep_figures_v3_gh_format_files/figure-gfm/unnamed-chunk-44-1.png)<!-- -->

### Boxplots specific proteins

``` r
quant_anwgene <- left_join(quant_annotated, prot2gene)


tab_2_paired_boxplot <-  quant_anwgene %>%
  filter(Gene %in% c("ASAH1", "SYNM", "GPNMB"))
# prep data
  slim_data <- pivot_longer(expr_matrix_filt, 
                          cols = 2:length(names(expr_matrix_filt)), 
                          names_to = c("patient","condition"),
                          values_to = "Intensity",
                          names_sep = "\\_") %>% 
  group_by(Protein)
  
  slim_data2 <- pivot_wider(slim_data, 
                          values_from = Intensity, 
                          names_from = condition) %>%
    filter(patient != "x6")

  slim_data3 <- slim_data2 %>% 
    mutate(Abs_diff = rec-prim) %>%
    mutate(rank = row_number(Abs_diff)) %>% 
    dplyr::rename(Initial = prim, Recurrent = rec) %>%
    ungroup() %>%
    filter(Protein %in% c("Q13510", "Q14956", "O15061")) %>%
    mutate(Gene = case_when(Protein == "Q13510" ~ "ASAH1",
                              Protein == "Q14956" ~ "GPNMB",
                              Protein == "O15061" ~ "SYNM"))
```

``` r
ggpaired(slim_data3,
         cond1 = "Initial",
         cond2 = "Recurrent", fill = "condition",
         ylab = "Normalized Abundance",
         label = NULL,
         repel = TRUE,
         facet.by = "Gene", point.size = 0.2, line.size = 0.1) +
      theme(axis.text.x = element_text(hjust = 0.5, vjust = 0, size = 5, angle = 360),
        axis.text.y = element_text(hjust = 0.95, vjust = 0.2, size = 5),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 8),
        legend.key.height= unit(3, 'mm'),
        legend.key.width= unit(3, 'mm'),
        legend.position="bottom",
        strip.text = element_text(hjust = 0.5, vjust = 0, size = 8, angle = 360)) 
```

![](gbm_proteomics_data_analysis_proteolytic_products_prep_figures_v3_gh_format_files/figure-gfm/unnamed-chunk-50-1.png)<!-- -->

## ID coverage (final)

### Coverage plot

``` r
ids_summary2 <- quant_annotated %>% 
  dplyr::select(Protein, mixture, Abundance) %>% 
  distinct(Protein, mixture, .keep_all = TRUE) %>%
  group_by(mixture) %>%
  summarise(`Nr of IDs` = sum(!is.na(Abundance))) %>% 
  bind_rows(.,
            tibble(mixture = "Quant_in_all",
                   `Nr of IDs` = nrow(expr_matnona))) %>%
  bind_rows(.,
            tibble(mixture = "Quant_2/3_mixtures",
                   `Nr of IDs` = nrow(mat_filt_impwo6))) %>%
  bind_rows(.,
            tibble(mixture = "Overall",
                   `Nr of IDs` = nrow(expr_matrix)))

ids_summary2 <- mutate(ids_summary2,
                       Mixture = case_when(mixture == "tmt_1" ~ "Mix_1",
                                           mixture == "tmt_2" ~ "Mix_2",
                                           mixture == "tmt_3" ~ "Mix_3",
                                           mixture == "Quant_in_all" ~ "Quant. in all",
                                           mixture == "Overall" ~ "Overall",
                                           mixture == "Quant_2/3_mixtures" ~ "Quant. in 2/3")) %>%
  mutate(Mixture = factor(Mixture, levels = c("Overall","Quant. in 2/3", "Quant. in all",
                                               "Mix_3", "Mix_2", "Mix_1")))
```

Number of identified proteins: 5954

Number of identified and quantified proteins in all samples: 3228

Number of identified and quantified proteins 2 out of 3 batches: 4464

``` r
ggplot(data = ids_summary2,
       aes(y = `Nr of IDs`, x = mixture)) +
  geom_col(fill = "#008683") + 
  
  geom_text(size = 5, position = position_stack(vjust = 0.5), aes(label = `Nr of IDs`)) +
  coord_cartesian(ylim = c(0, 5000))+
  coord_flip() +
  labs(title = "Number of Identified and quantified Proteins by Mixture") +
  ylab("Nr of IDed and Quant Proteins") + 
  xlab("TMT Mixture") + 
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.1, size = 14, angle = 0),
        axis.text.y = element_text(hjust = 0.5, vjust = 0.1, size = 12, angle = 0),
        panel.background = element_blank(),
        panel.grid.major = element_line(color = "grey"),
        panel.border = element_rect(colour = "black", fill=NA, size=1.5),
        axis.title=element_text(size=12, face="bold"))
```

![](gbm_proteomics_data_analysis_proteolytic_products_prep_figures_v3_gh_format_files/figure-gfm/unnamed-chunk-54-1.png)<!-- -->

## Median abundance distribution of filtered and inputed samples

``` r
expr_matrix_filtwo6 <- as.data.frame(mat_filt_imp) %>%
  rownames_to_column("Protein") %>%
  dplyr::select(-starts_with("x6"))

long_mat_f <- pivot_longer(expr_matrix_filtwo6,
                         cols = starts_with("x"),
                         names_sep = "_",
                         names_to = c("patient", "recurrence"),
                         values_to = "Abundance")

long_mat_3 <- long_mat_f %>%
  mutate(paired_id = paste(patient, recurrence, sep = "_"))

paired_annot_uniq2 <- long_mat_3 %>%
  dplyr::select(patient, recurrence, paired_id) %>%
  distinct()

paired_annot_uniq2$paired_id
```

    ##  [1] "x1_prim"  "x1_rec"   "x2_prim"  "x2_rec"   "x3_prim"  "x3_rec"  
    ##  [7] "x4_prim"  "x4_rec"   "x5_prim"  "x5_rec"   "x7_prim"  "x7_rec"  
    ## [13] "x8_prim"  "x8_rec"   "x9_prim"  "x9_rec"   "x10_prim" "x10_rec" 
    ## [19] "x11_prim" "x11_rec"

``` r
## Annotated abundance data in long format ----

quant_annotated_f <- left_join(long_mat_3,
                             sample_annotation2, 
                             by = c("paired_id", "patient","recurrence")) %>%
  mutate(Channel_mix = paste(mixture,channel,  sep = "_"),
         paired_id = factor(paired_id, levels = paired_annot_uniq2$paired_id)) %>%
  mutate(Sample = str_replace_all(paired_id, pattern = "prim", replacement = "init"),
         Stage = if_else(recurrence == "prim",
                         true = "Initial",
                         false = "Recurrent"))
```

### By TMT channel/sample

#### General Abudance distribution by TMT channel/sample

``` r
qcplot_f <- ggplot(quant_annotated_f,
                 mapping = aes(x = Sample, 
                               y = Abundance, 
                               fill = Stage)) +
          geom_boxplot(outlier.size=0.1, size = 0.25) + 
          theme(axis.text.x = element_text(hjust = 0.5, vjust = 0, size = 7, angle = 90),
                axis.text.y = element_text(hjust = 0.5, vjust = 0, size = 6),
                panel.background = element_blank(),
                panel.grid.major = element_blank(),
                panel.border = element_rect(colour = "black", fill=NA, size=0.5),
                axis.title=element_text(size=8),
                axis.title.x = element_blank(),
                legend.text = element_text(size = 7),
                legend.title = element_text(size = 8),
                legend.key.height= unit(3, 'mm'),
                legend.key.width= unit(3, 'mm'),
                legend.position="bottom")
```

![](gbm_proteomics_data_analysis_proteolytic_products_prep_figures_v3_gh_format_files/figure-gfm/unnamed-chunk-59-1.png)<!-- -->

## Over-representation analyses/tests (ORAs)

### Prep tabular info for ORAs

``` r
tab_sig_prots_tryptic <- bind_rows(tibble(protein = increased_in_rec,
                                  characteristic = "up-regulated"),
                           tibble(protein = decreased_in_rec,
                                  characteristic = "down-regulated"))

limma_tab_wo62_tryptic <- limma_tab_wo6


unip2symb_tryptic <- bitr(limma_tab_wo62_tryptic$Protein, 
                  fromType = "UNIPROT", 
                  toType = c("SYMBOL", "ENTREZID"),
                  OrgDb = org.Hs.eg.db) %>%
  dplyr::rename(protein = UNIPROT) %>%
  mutate(Protein = protein)

tab_sig_prots_tryptic <- left_join(tab_sig_prots_tryptic, unip2symb_tryptic)

unip2symbIDed_tryptic <- bitr(identified_proteins$protein_id, 
                  fromType = "UNIPROT", 
                  toType = c("SYMBOL", "ENTREZID"),
                  OrgDb = org.Hs.eg.db)


limma_tab_wo62_tryptic <- left_join(limma_tab_wo62_tryptic, unip2symb_tryptic)


geneList <- limma_tab_wo62_tryptic$logFC
names(geneList) <- limma_tab_wo62_tryptic$Protein
geneList <-  sort(geneList, decreasing = TRUE)

geneListr <- limma_tab_wo62_tryptic$logFC
names(geneListr) <- limma_tab_wo62_tryptic$ENTREZID
geneListr <-  sort(geneListr, decreasing = TRUE)
```

### Reactome

``` r
group_comparison_react <- compareCluster(ENTREZID~characteristic, 
                                              data=tab_sig_prots_tryptic, 
                                              fun="enrichPathway",
                                              organism = "human",
                                              pvalueCutoff = 0.05,
                                              pAdjustMethod = "BH",
                                              qvalueCutoff = 0.2,
                                              universe = unip2symbIDed_tryptic$ENTREZID,
                                              minGSSize = 10,
                                              maxGSSize = 1000,
                                              readable = TRUE)
```

#### Dotplot Reactome

``` r
enrichplot::dotplot(group_comparison_react, x = "characteristic") + 
          xlab("Quant in Recurrent status") +
  scale_color_continuous(low="red", high="blue",
            guide=guide_colorbar(reverse=TRUE),breaks = c(0.01, 0.04)) +
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0, size = 8, angle = 360),
        axis.text.y = element_text(hjust = 0.95, vjust = 0.2, size = 7),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8),
        legend.key.height= unit(3, 'mm'),
        legend.key.width= unit(3, 'mm'),
        legend.position="bottom")
```

![](gbm_proteomics_data_analysis_proteolytic_products_prep_figures_v3_gh_format_files/figure-gfm/unnamed-chunk-63-1.png)<!-- -->![](gbm_proteomics_data_analysis_proteolytic_products_prep_figures_v3_gh_format_files/figure-gfm/unnamed-chunk-63-2.png)<!-- -->

``` r
dotplot_react <- enrichplot::dotplot(group_comparison_react, x = "characteristic") + 
          xlab("Quant in Recurrent status") +
  scale_color_continuous(low="red", high="blue",
            guide=guide_colorbar(reverse=TRUE),breaks = c(0.01, 0.04)) +
  scale_size_continuous(breaks = c(0.10, 0.25)) +
  theme(axis.text.x = element_text(hjust = 0.4, vjust = 0.1, size = 6, angle = 340),
        axis.text.y = element_text(hjust = 0.95, vjust = 0.2, size = 6),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.key.height= unit(2, 'mm'),
        legend.key.width= unit(2, 'mm'),
        legend.position="bottom")
```

![](gbm_proteomics_data_analysis_proteolytic_products_prep_figures_v3_gh_format_files/figure-gfm/unnamed-chunk-64-1.png)<!-- -->

``` r
ggsave(plot = dotplot_react, 
       filename = here::here("figures/dotplot_reactome_proteomics_gbm.tiff"), 
       device = "tiff",
       units = "mm",
       width = 91,
       height = 90)

ggsave(plot = dotplot_react, 
       filename = here::here("figures/dotplot_reactome_proteomics_gbm.eps"), 
       device = "eps",
       units = "mm",
       width = 91,
       height = 90)
```

#### Cnet plot Reactome

``` r
enrichplot::cnetplot(group_comparison_react) + 
  theme(legend.text = element_text(size = 7),
                legend.title = element_text(size = 8),
                legend.key.height= unit(3, 'mm'),
                legend.key.width= unit(3, 'mm'),
                legend.position="bottom")
```

![](gbm_proteomics_data_analysis_proteolytic_products_prep_figures_v3_gh_format_files/figure-gfm/unnamed-chunk-66-1.png)<!-- -->

``` r
cnet_plot <- enrichplot::cnetplot(group_comparison_react) + 
  theme(legend.text = element_text(size = 7),
        text = element_text(size = 4),
                legend.title = element_text(size = 8),
                legend.key.height= unit(3, 'mm'),
                legend.key.width= unit(3, 'mm'),
                legend.position="bottom")
```

``` r
ggsave(plot = cnet_plot, 
       filename = here::here("figures/cnet_plot_reactome_proteomics_gbm.tiff"), 
       device = "tiff",
       units = "mm",
       width = 110,
       height = 110)

ggsave(plot = cnet_plot, 
       filename = here::here("figures/cnet_plot_reactome_proteomics_gbm.eps"), 
       device = "eps",
       units = "mm",
       width = 110,
       height = 110)

ggsave(plot = cnet_plot, 
       filename = here::here("figures/cnet_plot_reactome_proteomics_gbm.svg"), 
       device = "svg",
       units = "mm",
       width = 110,
       height = 110)
```

## Protein coverage plots

``` r
source(here("scr/protein_peptide_coverage_plot.R"))
```

### Annotate peptides from interesting proteins

``` r
source(here("scr/annotate_peptides.R"))
```

``` r
# get protein to peptide pairs
pept2prot_inter <- identified_peptides %>%
  filter(gene %in% c("ASAH1", "SYNM", "GPNMB")) %>%
  dplyr::select(Peptide = sequence, Genes = protein_id)

annotated_peptides <- annotate_peptides(expr_mat = pept2prot_inter,
                                        fasta = fasta_ident_prots) 
```

    ## [1] "1 AALEALLGR out of 75"
    ## [1] "2 AASLTMHFR out of 75"
    ## [1] "3 AELQELNAR out of 75"
    ## [1] "4 AIDCLEDEK out of 75"
    ## [1] "5 ARLEDALLR out of 75"
    ## [1] "6 ATGPAAPPPR out of 75"
    ## [1] "7 AVSSQTNVR out of 75"
    ## [1] "8 AYVPIAQVK out of 75"
    ## [1] "9 DCPDPCIGW out of 75"
    ## [1] "10 DGAVGEK out of 75"
    ## [1] "11 DKVAAGASESTR out of 75"
    ## [1] "12 DVYVVTDQIPVFVTMFQK out of 75"
    ## [1] "13 DYQDLLQVK out of 75"
    ## [1] "14 EALGLEQLR out of 75"
    ## [1] "15 EGLWAEGQAR out of 75"
    ## [1] "16 EGQGGPGSVSVDVK out of 75"
    ## [1] "17 EIIFQGPISAAGK out of 75"
    ## [1] "18 EKASEER out of 75"
    ## [1] "19 ELQEALGAR out of 75"
    ## [1] "20 ELYIPSGESEVAGGASHSSGQR out of 75"
    ## [1] "21 ENLLLEEELR out of 75"
    ## [1] "22 EQAMFDKK out of 75"
    ## [1] "23 EQGKEQAMFDKK out of 75"
    ## [1] "24 ESLDVYELDAK out of 75"
    ## [1] "25 EVPISLEVSQDRR out of 75"
    ## [1] "26 EVPVYIGEDSTIAR out of 75"
    ## [1] "27 GAVPWYTINLDLPPYK out of 75"
    ## [1] "28 GHLIHGR out of 75"
    ## [1] "29 GQFETYLR out of 75"
    ## [1] "30 GRLDAELGAQQR out of 75"
    ## [1] "31 GSQTGTSIGGDAR out of 75"
    ## [1] "32 HPFFLDDR out of 75"
    ## [1] "33 HPFFLDDRR out of 75"
    ## [1] "34 IMQVVDEK out of 75"
    ## [1] "35 KSTYPPSGPTYR out of 75"
    ## [1] "36 LDAELGAQQR out of 75"
    ## [1] "37 LEDALLR out of 75"
    ## [1] "38 LGPSEVWR out of 75"
    ## [1] "39 LGPTETETSEHIAIR out of 75"
    ## [1] "40 LPGLLGNFPGPFEEEMK out of 75"
    ## [1] "41 LQTGPEK out of 75"
    ## [1] "42 LTEDVDVSDEAGLDYLLSK out of 75"
    ## [1] "43 LTVYTTLIDVTK out of 75"
    ## [1] "44 LYDYVCR out of 75"
    ## [1] "45 MREELSALTR out of 75"
    ## [1] "46 MREEYGIQAEER out of 75"
    ## [1] "47 NDQAVGVSFK out of 75"
    ## [1] "48 NMINTFVPSGK out of 75"
    ## [1] "49 NRPETIR out of 75"
    ## [1] "50 NTEAQVK out of 75"
    ## [1] "51 QFTQSPETEASADSFPDTK out of 75"
    ## [1] "52 QQLDELSWATALAEGER out of 75"
    ## [1] "53 QQLVEVIGQLEETLPER out of 75"
    ## [1] "54 QTVMTEK out of 75"
    ## [1] "55 RGFLGSGYSSSATTQQENSYGK out of 75"
    ## [1] "56 RGLDAAHER out of 75"
    ## [1] "57 SAEQMIGDIINLGLK out of 75"
    ## [1] "58 SGEFHAEPTVIEK out of 75"
    ## [1] "59 STYPPSGPTYR out of 75"
    ## [1] "60 SVISDEK out of 75"
    ## [1] "61 SYHYTDSLLQR out of 75"
    ## [1] "62 TGLSLEVATYR out of 75"
    ## [1] "63 TKTEIVVESK out of 75"
    ## [1] "64 TPQGPVSATVEVSSPTGFAQSQVLEDVSQAAR out of 75"
    ## [1] "65 TQEAGALGVSDR out of 75"
    ## [1] "66 TQKDGAVGEK out of 75"
    ## [1] "67 TSQENISFETMYDVLSTKPVLNK out of 75"
    ## [1] "68 VAAGASESTR out of 75"
    ## [1] "69 VGDYFATEESVGTQTSVR out of 75"
    ## [1] "70 VIVNSLK out of 75"
    ## [1] "71 VTYVDRK out of 75"
    ## [1] "72 WHELMLDKAPVLK out of 75"
    ## [1] "73 WKHPFFLDDR out of 75"
    ## [1] "74 WYVVQTNYDR out of 75"
    ## [1] "75 YSWQDEIVQGTR out of 75"

``` r
annotated_peptides2 <- annotated_peptides %>%
  mutate(peptide_length = str_length(Peptide))
```

### ASAH1

``` r
feat_asah1 <- drawProteins::get_features("Q13510")
```

    ## [1] "Download has worked"

``` r
feat_asah1df <- drawProteins::feature_to_dataframe(feat_asah1)
#Get peptides as coverage features for plotting

cov_feat_asah1 <- get_coverage(annotated_peptides = annotated_peptides2, 
                               id = "Q13510")
#Merge coverage feature info with protein sequence features info:

fullfeat_asah1 <- add_peptides(feat_asah1df,
                               peptide_coverage_data = cov_feat_asah1)
# drawCanvas
canvas <- draw_canvas(fullfeat_asah1)

# draw the protein chain
wchain <- draw_chains(canvas, fullfeat_asah1)

# draw protein domains
wdomainsasah1 <- draw_domains(wchain, data = fullfeat_asah1)
#pasa <- drawProteins::draw_folding(wdomainsasah1, fullfeat_asah1)
cov_feat_asah1$order <- 1.4

pepcov_asah1 <- draw_peptides(wdomainsasah1, cov_feat_asah1)
pepcov_asah1
```

![](gbm_proteomics_data_analysis_proteolytic_products_prep_figures_v3_gh_format_files/figure-gfm/unnamed-chunk-72-1.png)<!-- -->

``` r
ggsave(plot = pepcov_asah1, 
       filename = here::here("figures/asah1_coverage_proteomics_gbm.tiff"), 
       device = "tiff",
       units = "mm",
       width = 150,
       height = 100)

ggsave(plot = pepcov_asah1, 
       filename = here::here("figures/asah1_coverage_proteomics_gbm.eps"), 
       device = "eps",
       units = "mm",
       width = 150,
       height = 100)
```

### SYNM coverage plot

Get features from Uniprot

``` r
feat_SYNM <- drawProteins::get_features("O15061")
```

    ## [1] "Download has worked"

``` r
feat_SYNMdf <- drawProteins::feature_to_dataframe(feat_SYNM)
#Get peptides as coverage features for plotting

cov_feat_SYNM <- get_coverage(annotated_peptides = annotated_peptides2, 
                               id = "O15061")
#Merge coverage feature info with protein sequence features info:

fullfeat_SYNM <- add_peptides(feat_SYNMdf,
                               peptide_coverage_data = cov_feat_SYNM)
# drawCanvas
canvasSYNM <- draw_canvas(fullfeat_SYNM)

# draw the protein chain
wchainSYNM <- draw_chains(canvasSYNM, fullfeat_SYNM)
cov_feat_SYNM$order <- 1.4

pepcov_SYNM <- draw_peptides(wchainSYNM, cov_feat_SYNM)
pepcov_SYNM
```

![](gbm_proteomics_data_analysis_proteolytic_products_prep_figures_v3_gh_format_files/figure-gfm/unnamed-chunk-74-1.png)<!-- -->

``` r
ggsave(plot = pepcov_SYNM, 
       filename = here::here("figures/SYNM_coverage_proteomics_gbm.tiff"), 
       device = "tiff",
       units = "mm",
       width = 150,
       height = 100)

ggsave(plot = pepcov_SYNM, 
       filename = here::here("figures/SYNM_coverage_proteomics_gbm.eps"), 
       device = "eps",
       units = "mm",
       width = 150,
       height = 100)
```

# Analysis of proteolytic processing

## Load semi-specific search results

``` r
tmt_peptdata_frag <- read_tsv(here("data/semi_specific_search_fragpipe17/tmt-report/abundance_peptide_MD.tsv")) %>% 
  janitor::clean_names()


pept_ids_mix1 <- read_delim(here("data/semi_specific_search_fragpipe17/mix_1/peptide.tsv")) %>%
                    clean_names() %>%
          dplyr::select(-c(starts_with("empty")))

pept_ids_mix2 <- read_delim(here("data/semi_specific_search_fragpipe17/mix_2/peptide.tsv")) %>%
                    clean_names() %>%
          dplyr::select(-c(starts_with("x"), starts_with("empty"), "master"))

pept_ids_mix3 <- read_delim(here("data/semi_specific_search_fragpipe17/mix_3/peptide.tsv")) %>%
                    clean_names() %>%
          dplyr::select(-c(starts_with("x"), starts_with("empty"), "master"))

identified_proteins_semi <- read_tsv(here("data/semi_specific_search_fragpipe17/combined_protein.tsv")) %>% 
  janitor::clean_names()
```

## Wrangle and merge peptide ID/quant with sample annotation

``` r
# wide expression matrix with NAs
pept_matrix <- dplyr::select(tmt_peptdata_frag,
                             index, starts_with("x"))

pept_long_mat <- pivot_longer(pept_matrix,
                         cols = starts_with("x"),
                         names_to = c("patient"),
                         values_to = "Abundance") %>%
  mutate(recurrence = str_extract(patient, "[^0-9]+$")) %>%
  mutate(patient = str_remove_all(patient, "[^0-9]+$"))

pept_long_mat_2 <- pept_long_mat %>%
  mutate(paired_id_semi = paste(patient, recurrence, sep = "_"))

sample_annotation2 <- sample_annotation2 %>%
  mutate(paired_id_semi = paste(patient, recurrence, sep = "_"))

## Annotated abundance data in long format ----

pept_quant_annotated <- left_join(pept_long_mat_2,
                             sample_annotation2, 
                             by = c("paired_id_semi", "patient","recurrence")) %>%
  mutate(Channel_mix = paste(mixture,channel,  sep = "_"))

## Wide xpression matrices without NAs
pept_expr_matnona <- column_to_rownames(pept_matrix,
                            "index") %>%
          as.matrix() %>%
          na.omit()
```

Wrangle and merge peptide ids for further annotation based on they
proteolytic specificity

``` r
peptide_joined1 <- full_join(pept_ids_mix1, pept_ids_mix2) %>%
                    mutate(duplicated = duplicated(peptide))

peptide_joined <- full_join(peptide_joined1, pept_ids_mix3) %>%
                    mutate(duplicated = duplicated(peptide))


peptide_ids_filtered <- filter(peptide_joined,
                               peptide %in% peptide_joined$peptide) %>%
                    distinct(peptide, assigned_modifications,observed_modifications, .keep_all = TRUE) %>%
                    dplyr::rename(Peptide = peptide) %>%
                    dplyr::mutate(Genes = protein_id)
```

### Annotate peptides with cleavage information

``` r
fasta1 <- seqinr::read.fasta(here("data/semi_specific_search_fragpipe17/mix_1/protein.fas"),
                            seqtype = "AA", as.string = TRUE)

fasta2 <- seqinr::read.fasta(here("data/semi_specific_search_fragpipe17/mix_2/protein.fas"),
                             seqtype = "AA", as.string = TRUE)

fasta3 <- seqinr::read.fasta(here("data/semi_specific_search_fragpipe17/mix_3/protein.fas"),
                             seqtype = "AA", as.string = TRUE)

fasta22 <- append(fasta1, fasta2)

fasta <- append(fasta22, fasta3)

fasta <- discard(fasta,
                 duplicated(names(fasta)))

attr(fasta[[1]], "Annot")
```

    ## [1] ">sp|A0A024RBG1|NUD4B_HUMAN Diphosphoinositol polyphosphate phosphohydrolase NUDT4B OS=Homo sapiens OX=9606 GN=NUDT4B PE=3 SV=1"

``` r
extr_annt <- function(x){
  ext <- attr(x, "Annot")
}

names <- purrr::map_chr(fasta, .f = extr_annt)

seqinr::write.fasta(fasta,
                    file.out = here("data/semi_specific_search_fragpipe17/protein_combined.fas"),
                    names = names)

if(!file.exists(here("results/semi_tryptic/cleavage_annoated_peptides_final.tsv"))){
  cleavage_annoated_peptides <- annotate_peptides(expr_mat = peptide_ids_filtered, 
                                                  fasta = fasta,
                                                  decoy_tag = "rev_")
  
  write_tsv(cleavage_annoated_peptides, here("results/semi_tryptic/cleavage_annoated_peptides_final.tsv"))
} else {
  cleavage_annoated_peptides <- read_tsv(here("results/semi_tryptic/cleavage_annoated_peptides_final.tsv"))
}

cleavage_annoated_peptides <- dplyr::mutate(cleavage_annoated_peptides,
                                            Genes = protein_id)

peptides_annotated <- left_join(peptide_ids_filtered, cleavage_annoated_peptides,
                                by = c("Peptide", "protein_id"))
```

## Sparcity reduction

Only peptides that were identified and quantified in all 3 mixtures
would be used for further inferential analyses for proteolytic
processing.

``` r
pept_matrix_nona <- pept_matrix %>% 
  na.omit()

pept_nona <- pept_matrix_nona %>%
  dplyr::select(index) %>%
  separate(col = index, 
           sep = "_", into = c("protein", "peptide"))

pept_matrix_nonawop6 <- pept_matrix %>% 
  dplyr::select(-starts_with("x6")) %>%
  na.omit()
```

## Coverage information (peptide-level)

``` r
ids_summarypep <- pept_quant_annotated %>% 
  dplyr::select(index, mixture, Abundance) %>% 
  distinct(index, mixture, .keep_all = TRUE) %>%
  group_by(mixture) %>%
  summarise(`Nr of IDs` = sum(!is.na(Abundance))) %>% 
  bind_rows(.,
            tibble(mixture = "In all three",
                   `Nr of IDs` = nrow(pept_matrix_nona)))# %>%
  #bind_rows(.,
  #          tibble(mixture = "Overall",
  #                `Nr of IDs` = nrow(pept_matrix))) 

ids_summarypep <- mutate(ids_summarypep,
                       Mixture = case_when(mixture == "tmt_1" ~ "Mix_1",
                                           mixture == "tmt_2" ~ "Mix_2",
                                           mixture == "tmt_3" ~ "Mix_3",
                                           mixture == "In all three" ~ "Quant. in all",
                                           mixture == "Overall" ~ "Overall",
                                           mixture == "Quant_2/3_mixtures" ~ "Quant. in 2/3")) %>%
  mutate(Mixture = factor(Mixture, levels = c("Overall", "Quant. in all",
                                               "Mix_3", "Mix_2", "Mix_1")))
```

``` r
ggplot(data = ids_summarypep,
       aes(y = `Nr of IDs`, x = Mixture)) +
  geom_col(fill = "#008683") + 
  geom_text(size = 5, position = position_stack(vjust = 0.5), 
            aes(label = `Nr of IDs`)) +
  coord_cartesian(ylim = c(0, 5000))+
  coord_flip() +
  labs(title = "Number of Identified and quantified Peptides by Mixture") +
  ylab("Nr of IDed and Quant Peptides") + 
  xlab("TMT Mixture") + 
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0, size = 7, angle = 360),
        axis.text.y = element_text(hjust = 0.95, vjust = 0.2, size = 7),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8),
        legend.key.height= unit(3, 'mm'),
        legend.key.width= unit(3, 'mm'),
        legend.position="bottom")
```

![](gbm_proteomics_data_analysis_proteolytic_products_prep_figures_v3_gh_format_files/figure-gfm/unnamed-chunk-82-1.png)<!-- -->

``` r
barplot_pept_cover <- ggplot(data = ids_summarypep,
       aes(y = `Nr of IDs`, x = Mixture)) +
  geom_col(fill = "#008683") + 
  geom_text(size = 2, position = position_stack(vjust = 0.5), aes(label = `Nr of IDs`)) +
  coord_flip() +
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0, size = 7, angle = 360),
        axis.text.y = element_text(hjust = 0.95, vjust = 0.2, size = 7),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8),
        legend.key.height= unit(3, 'mm'),
        legend.key.width= unit(3, 'mm'),
        legend.position="bottom")
```

``` r
ggsave(plot = barplot_pept_cover, 
       filename = here::here("figures/barplot_coverage_pept_proteomics_gbm.tiff"), 
       device = "tiff",
       units = "mm",
       width = 38,
       height = 50)

ggsave(plot = barplot_pept_cover, 
       filename = here::here("figures/barplot_coverage_pept_proteomics_gbm.eps"), 
       device = "eps",
       units = "mm",
       width = 38,
       height = 50)
```

## Label N-terminal peptides

``` r
source(here("scr/annotate_nterm.R"))
```

``` r
nterannot <- annotate_nterm(peptides_annotated,
                            tmtmass = 229.1629) %>%
  clean_names()
```

### Merge peptide annotation with TMT integrator quant report

``` r
# all quantified peptides
tmt_reprt_annot <- left_join(tmt_peptdata_frag, nterannot) 

check_na <- filter(tmt_reprt_annot, is.na(semi_type))

excluded_tmt_reprt_annot <- anti_join(nterannot, tmt_reprt_annot) %>% 
                    clean_names()

#only those present in all mixtures
tmt_peptdata_fragnona <- tmt_peptdata_frag %>%
  filter(peptide %in% pept_nona$peptide)

tmt_reprt_annotnona <- left_join(tmt_peptdata_fragnona, nterannot) 
```

## Summarize feature counts per feature

``` r
source(here("scr/summarize_peptide_counts.R"))
```

``` r
summary_count <- summarize_peptide_counts(tmt_reprt_annot) 
```

### Barplot

These are the ones that would be used for inferential statistics.

``` r
to_count_info <- tmt_reprt_annot %>% 
   filter(index %in% pept_matrix_nona$index) %>%
                      dplyr::select(index, specificity, nterm, semi_type, tmt_tag,
                                    aa_after, aa_before, following_10_resid, previous_10_resid, start_position, end_position) %>%
  mutate(specificity = if_else(str_detect(semi_type, "unspecific"), 
                               true = "semi_specific",
                               false = specificity),
         duplicated = duplicated(index)) %>%
  distinct()
  
to_count_info_semi <- to_count_info %>%
                           dplyr::select(index,specificity) %>%
                           distinct() %>%
  mutate(duplicated = duplicated(index)) 

false_duplication <- to_count_info_semi %>%
  filter(duplicated) %>%
  pull(index)

to_count_info_semi <- to_count_info_semi %>%
  mutate(specificity = if_else(index %in% false_duplication, 
                               true = "semi_specific",
                               false = specificity)) %>%
  dplyr::select(-duplicated) %>%
  distinct()

  
  n_semi <- dplyr::count(to_count_info_semi, specificity) %>% 
                      dplyr::rename(feature_type = specificity) %>%
                      dplyr::mutate(category = "specificity")
  

  n_total <- tibble(feature_type = "Total",
                    n = length(unique(to_count_info$index)),
                    category = "Total")
  
  n_term <- dplyr::count(to_count_info, nterm) %>% 
                      dplyr::rename(feature_type = nterm) %>%
                      dplyr::mutate(category = "N-term")
  
  summary_count_pept_nona <- bind_rows(n_semi,
                                       n_term,
                                       n_total)


total_sumcount <- summary_count_pept_nona %>%
  filter(feature_type == "Total") %>%
  pull(n)

summary_countnona_no_tot <- summary_count_pept_nona %>%
  dplyr::filter(category != "Total") %>%
  dplyr::rename(`Feature category` = category)
```

![](gbm_proteomics_data_analysis_proteolytic_products_prep_figures_v3_gh_format_files/figure-gfm/unnamed-chunk-91-1.png)<!-- -->

``` r
ggsave(plot = summary_pept_annot, 
       filename = here::here("figures/barplot_coverage_peptidesfeatures_gbm.tiff"), 
       device = "tiff",
       units = "mm",
       width = 76,
       height = 50)

ggsave(plot = summary_pept_annot, 
       filename = here::here("figures/barplot_coverage_peptidesfeatures_gbm.eps"), 
       device = "eps",
       units = "mm",
       width = 76,
       height = 50)
```

## Inferential analysis (peptide-level)

### Prep design matrix

``` r
design_limmawo6 <- model.matrix(~patientwo6+recurrencewo6)

sample_annotationwo6 <- sample_annotationwo6 %>%
  mutate(paired_id_semi = paste(patient, recurrence, sep = ""))

rownames(design_limmawo6) <- sample_annotationwo6$paired_id_semi
```

### Prep abundance matrices

``` r
pep_matwo6 <- dplyr::select(pept_matrix_nonawop6,
                             index, rownames(design_limmawo6))  %>%
  column_to_rownames("index") %>%
  as.matrix()
```

``` r
colnames(pep_matwo6)
```

    ##  [1] "x1prim"  "x1rec"   "x2prim"  "x2rec"   "x3prim"  "x3rec"   "x4prim" 
    ##  [8] "x4rec"   "x5prim"  "x5rec"   "x7prim"  "x7rec"   "x8prim"  "x8rec"  
    ## [15] "x9prim"  "x9rec"   "x10prim" "x10rec"  "x11prim" "x11rec"

``` r
rownames(design_limmawo6)
```

    ##  [1] "x1prim"  "x1rec"   "x2prim"  "x2rec"   "x3prim"  "x3rec"   "x4prim" 
    ##  [8] "x4rec"   "x5prim"  "x5rec"   "x7prim"  "x7rec"   "x8prim"  "x8rec"  
    ## [15] "x9prim"  "x9rec"   "x10prim" "x10rec"  "x11prim" "x11rec"

The column names in our peptide abundance matrix have the same order as
the row names in our design matrix.

### Peptide-level model fit

``` r
fit_limmapeptwo6 <- function(mat, design, method,Limma){
  limmafit <- lmFit(mat, design, method = method)
  limmafit <- eBayes(limmafit)
  limma_tab <- topTable(limmafit, coef = "recurrencewo6rec", number = Inf, adjust.method = "BH") %>%
  mutate(Protein = rownames(.),
         Limma = Limma)
  
  return(limma_tab)
}
```

``` r
fit_mod_peptwo6 <- fit_limmapeptwo6(pep_matwo6, 
                         design_limmawo6, 
                         method = 'robust', 
                         Limma = "Robust - w Patient effect") 

compar_tab_peptwo6 <-  fit_mod_peptwo6 %>%
                    separate(Protein, 
                             into = c("protein", "peptide"), 
                             remove = FALSE) %>%
  dplyr::rename(index = Protein)
```

### Extract coefficients and feature specific FDR correction

#### Interesting features

``` r
features1 <- tmt_reprt_annotnona %>%
                    dplyr::select(peptide, index, specificity, nterm, 
                                  semi_type, is_terminal)

interesting_features1 <- features1 %>%
                    filter(specificity == "semi_specific",
                           is_terminal == "not_terminal") 
```

``` r
# load function for feature specific fdr correction
source(here("scr/feature_specific_fdr_correction.R"))
```

``` r
# load function for volcano plot
source(here("scr/plot_volcano.R"))
```

``` r
# excluding patient 6
compar_tab_correctedwo6_pept <- feature_fdr_correction(compar_tab_peptwo6,
                                               features1,
                                               interesting_features1)

compar_tab_interestingwo6_pept <- compar_tab_correctedwo6_pept %>% 
  dplyr::rename(Protein = protein) %>%
  filter(fdr_correction == 'feature-specific') %>%
  left_join(.,prot2gene)

increased_recwo6_pept <- compar_tab_interestingwo6_pept %>%
                    filter(logFC>0,
                           adj.P.Val < 0.05) %>% 
  pull(peptide)

decreased_recwo6_pept <- compar_tab_interestingwo6_pept %>%
                    filter(logFC<0,
                           adj.P.Val < 0.05) %>% 
  pull(peptide)

check_increased_pept <- cleavage_annoated_peptides %>%
  filter(Peptide %in% increased_recwo6_pept)
```

``` r
length(increased_recwo6_pept)
```

    ## [1] 15

``` r
length(decreased_recwo6_pept)
```

    ## [1] 3

15 semi-specific peptides are increased in Recurrent tumor and 3 are
decreased in Recurrent tumor.

### Volcano peptides

``` r
recwo6_pept22 <- compar_tab_interestingwo6_pept %>%
  filter(adj.P.Val < 0.05)

size <- 0.4

volcano_pept <- ggplot(data = compar_tab_interestingwo6_pept,
                      mapping = aes(x = logFC, y = -log10(adj.P.Val))) +
      geom_point(data = compar_tab_interestingwo6_pept %>% filter(logFC > 0,
                                           adj.P.Val < 0.05),
                 mapping = aes(x = logFC, y = -log10(adj.P.Val)), 
                 color = "red",
                 size = size)+
      geom_point(data = compar_tab_interestingwo6_pept %>% filter(logFC < -0,
                                           adj.P.Val < 0.05),
                 mapping = aes(x = logFC, y = -log10(adj.P.Val)), color = "red",
                 size = size) +
      geom_point(data = compar_tab_interestingwo6_pept %>% filter(logFC > 0,
                                           adj.P.Val > 0.05),
                 mapping = aes(x = logFC, y = -log10(adj.P.Val)), color = "#2a9d8f",
                 size = size) +
      geom_point(data = compar_tab_interestingwo6_pept %>% filter(logFC < -0,
                                           adj.P.Val > 0.05),
                 mapping = aes(x = logFC, y = -log10(adj.P.Val)), color = "#2a9d8f",
                 size = size) +
      geom_hline(yintercept = -log10(0.05),
                 color = "red", linetype = "dashed") +
      xlab("logFC - Recurrent / Primary") + 
      theme(axis.text.x = element_text(hjust = 0.5, vjust = 0, size = 6, angle = 360),
        axis.text.y = element_text(hjust = 0.95, vjust = 0.2, size = 8),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 8),
        legend.key.height= unit(3, 'mm'),
        legend.key.width= unit(3, 'mm'),
        legend.position="bottom") 
```

``` r
ggsave(plot = volcano_pept, 
       filename = here::here("figures/volcano_proteolysis_gbm.tiff"), 
       device = "tiff",
       units = "mm",
       width = 63,
       height = 50)

ggsave(plot = volcano_pept, 
       filename = here::here("figures/volcano_proteolysis_gbm.eps"), 
       device = "eps",
       units = "mm",
       width = 63,
       height = 50)
```

## Proportional abundance of proteolytic peptides

### Prepare data

**Wrangle data**

``` r
pept_semi_annot <- peptides_annotated %>%
  dplyr::select(Peptide, specificity, is_terminal) 

pept_quant_to_summary <- pept_quant_annotated %>%
  separate(col = index, 
           into = c("Protein", "Peptide"), 
           sep = "\\_", 
           remove = FALSE) %>%
  left_join(.,pept_semi_annot)

pep_quant_presummary <- pept_quant_to_summary %>%
  dplyr::select(Protein, Peptide, recurrence, paired_id, 
                Abundance, specificity, is_terminal)

pep_quant_presummary_semi <- pep_quant_presummary %>%
  filter(specificity == "semi_specific",
         is_terminal == "not_terminal")
```

**Prepare data-frame with abudance percentage of semi-specific peptides
per sample**

``` r
pept_summary_int_all <- pep_quant_presummary %>%
  group_by(paired_id, recurrence) %>%
  summarise(Sum_All = sum(Abundance, na.rm = TRUE))

pept_summary_int_semi <- pep_quant_presummary %>%
  filter(specificity == "semi_specific",
         is_terminal == "not_terminal") %>%
  group_by(paired_id, recurrence) %>%
  summarise(Sum_Semi = sum(Abundance, na.rm = TRUE))

pept_sum_summary <- left_join(pept_summary_int_all, pept_summary_int_semi) %>%
  mutate(Percentage = Sum_Semi/Sum_All * 100,
         Stage = if_else(recurrence == "prim",
                         true = "Initial",
                         false = "Recurrence"))

# summary of sum of log2 semi-specific peptides per sample

pept_summary_semi_1 <- pep_quant_presummary_semi %>% 
  group_by(paired_id, recurrence) %>%
  summarise(Sum_All = sum(Abundance, na.rm = TRUE)) %>%
  mutate(Stage = if_else(recurrence == "prim",
                         true = "Initial",
                         false = "Recurrence"),
         `Summed Abundances` = Sum_All)
```

### Generate plot % of abundance of semi-specific peptides

``` r
plotSumEM <- ggplot(pept_sum_summary, 
                    aes(x = Stage, 
                        y = Percentage, fill = Stage, 
                        cex.axis = 1.5)) +
  geom_boxplot() +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5) +
  # Box plot with jittered points
  # 0.2 : degree of jitter in x direction
  # geom_jitter(shape=16, position=position_jitter(0.2))
  ylab("Proportional intensity of semi-tryptic peptides [%]") +
  geom_signif(
    comparisons = list(c("Initial", "Recurrence")),
    map_signif_level = TRUE
  ) + 
  stat_compare_means(method="t.test") +
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0, size = 6, angle = 360),
        axis.text.y = element_text(hjust = 0.95, vjust = 0.2, size = 8),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 8),
        legend.key.height= unit(3, 'mm'),
        legend.key.width= unit(3, 'mm'),
        legend.position="bottom") 
```

``` r
ggsave(plot = plotSumEM, 
       filename = here::here("figures/box_perc_semispec_pepts.tiff"), 
       device = "tiff",
       units = "mm",
       width = 40,
       height = 50)

ggsave(plot = plotSumEM, 
       filename = here::here("figures/box_perc_semispec_pepts.eps"), 
       device = "eps",
       units = "mm",
       width = 40,
       height = 50)
```

### Generate plot of sum of abundances semi-specific peptides

``` r
sum_semi_abunds <- ggplot(pept_summary_semi_1, 
                    aes(x = Stage, 
                        y = `Summed Abundances`, fill = Stage, 
                        cex.axis = 1.5)) +
  geom_boxplot() +
  #geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5,
  #             ) +
  geom_jitter(position=position_jitter(0.2)) + 
  # Box plot with jittered points
  # 0.2 : degree of jitter in x direction
  # geom_jitter(shape=16, position=position_jitter(0.2))
  ylab("Sum of log2-intensities of semi-tryptic peptides") +
  geom_signif(
    comparisons = list(c("Initial", "Recurrence")),
    map_signif_level = TRUE
  ) + 
  stat_compare_means(method="t.test") +
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0, size = 6, angle = 360),
        axis.text.y = element_text(hjust = 0.95, vjust = 0.2, size = 8),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 8),
        legend.key.height= unit(3, 'mm'),
        legend.key.width= unit(3, 'mm'),
        legend.position="bottom")
```

``` r
ggsave(plot = sum_semi_abunds, 
       filename = here::here("figures/box_sum_semispec_pepts.tiff"), 
       device = "tiff",
       units = "mm",
       width = 40,
       height = 50)

ggsave(plot = sum_semi_abunds, 
       filename = here::here("figures/box_sum_semispec_pepts.eps"), 
       device = "eps",
       units = "mm",
       width = 40,
       height = 50)
```

## Comparative analysis of semi-specific peptides vs protein abundance

### Prepare data

#### Get log2FC Summary of abundances of semi-specific peptides

``` r
#taking output from limma semi-specific
log2_semi_tryptic <- dplyr::select(compar_tab_interestingwo6_pept,
                                   logFC, Protein, peptide, specificity, 
                                   is_terminal,
                                   Gene) %>%
  filter(specificity == "semi_specific",
         is_terminal == "not_terminal") %>%
  dplyr::rename(logFC_semi_peptide = logFC) %>%
  dplyr::select(-c(specificity, is_terminal))
```

#### Get log2 Summary of protein abundances based only on fully-tryptic peptides

``` r
med_prot_quant_specific <- pep_quant_presummary %>%
  group_by(Protein, paired_id) %>%
  summarise(`median(Abundance)` = median(Abundance, na.rm = TRUE)) %>%
  ungroup()
```

##### Prep wide matrix of median protein abundances from only fully-tryptic peptides

``` r
# prep wide matrix for median proteins abundances (arousing only from fully-tryptic peptides)
wide_premat_prot_only_tryp <- med_prot_quant_specific %>%
  pivot_wider(id_cols = c("Protein"),
              values_from = `median(Abundance)`, 
              names_from = paired_id) %>%
  # keep only proteins present in 2 / 3 of the samples
  filter(!Protein %in% proteins2exclude) 

mat_prot_only_tryp <- wide_premat_prot_only_tryp %>%
  column_to_rownames("Protein") %>%
  as.matrix()
```

**missForest-based imputation**

``` r
t_mat_prot_only_tryp <- t(mat_prot_only_tryp)

if(!file.exists(here("rds/missforest_imp_mat_prot_only_tryp.Rds"))){
  
 imp_prot_only_tryp <- missForest::missForest(t_mat_prot_only_tryp)
 t_mat_imp_prot_only_tryp <- imp_prot_only_tryp$ximp
 
 mat_prot_only_tryp <- t(t_mat_prot_only_tryp)
 
 saveRDS(imp_prot_only_tryp,
         file = here("rds/missforest_imp_mat_prot_only_tryp.Rds"))
 
} else {
  
 imp_prot_only_tryp <- readRDS(here("rds/missforest_imp_mat_prot_only_tryp.Rds"))
 t_mat_imp_prot_only_tryp <- imp_prot_only_tryp$ximp
 
 mat_imp_prot_only_tryp <- t(t_mat_imp_prot_only_tryp)
  
}
```

**Run limma on this new data to generate log2 FCs of proteins from only
fully-tryptic peptides**

``` r
split <- str_split_fixed(colnames(mat_imp_prot_only_tryp), pattern = "\\_", n = 2)

patient_tryp <- split[,1]
recurrence <- split[,2]

design_tryp <- model.matrix(~patient_tryp+recurrence)

rownames(design_tryp) <- colnames(mat_imp_prot_only_tryp)
```

``` r
rownames(design_tryp) 
```

    ##  [1] "x1_prim"  "x1_rec"   "x10_prim" "x10_rec"  "x11_prim" "x11_rec" 
    ##  [7] "x2_prim"  "x2_rec"   "x3_prim"  "x3_rec"   "x4_prim"  "x4_rec"  
    ## [13] "x5_prim"  "x5_rec"   "x6_prim"  "x6_rec"   "x7_prim"  "x7_rec"  
    ## [19] "x8_prim"  "x8_rec"   "x9_prim"  "x9_rec"

``` r
colnames(mat_prot_only_tryp)
```

    ##  [1] "x1_prim"  "x1_rec"   "x10_prim" "x10_rec"  "x11_prim" "x11_rec" 
    ##  [7] "x2_prim"  "x2_rec"   "x3_prim"  "x3_rec"   "x4_prim"  "x4_rec"  
    ## [13] "x5_prim"  "x5_rec"   "x6_prim"  "x6_rec"   "x7_prim"  "x7_rec"  
    ## [19] "x8_prim"  "x8_rec"   "x9_prim"  "x9_rec"

``` r
fit_limmawo6_tryp <- function(mat, design, method,Limma, prot2gene){
  limmafit <- lmFit(mat, design, method = method)
  limmafit <- eBayes(limmafit)
  limma_tab <- topTable(limmafit, coef = "recurrencerec", number = Inf, adjust.method = "BH") %>%
  mutate(Protein = rownames(.),
         Limma = Limma)
  
  limma_tab <- left_join(limma_tab, prot2gene)
  
  return(limma_tab)
}
```

``` r
limma_tab_wo6tryp <- fit_limmawo6_tryp(mat_imp_prot_only_tryp, 
                         design_tryp, 
                         method = 'robust', 
                         Limma = "Robust - Proteins w only tryptic peptides",
                         prot2gene = prot2gene)
```

``` r
log2_fully_tryptic <- dplyr::select(limma_tab_wo6tryp,
                                   logFC, Protein, Gene) %>%
  dplyr::rename(logFC_fully_tryp_protein = logFC) 
```

#### Merge log2FCs of semi-tryptic peptides vs proteins with fully tryptic peptides

``` r
log2semipept2_log2protein_spec <- left_join(log2_fully_tryptic, log2_semi_tryptic) %>%
  na.omit() %>%
  mutate(DA_Protein = if_else(Protein %in% tab_sig_prots_tryptic$protein,
                              true = TRUE,
                              false = FALSE),
         DA_peptide = if_else(peptide %in% recwo6_pept22$peptide,
                              true = TRUE,
                              false = FALSE))
```

#### log2FC Proteins vs Semi-tryp peptides (scatter plot)

``` r
inter_scater <- log2semipept2_log2protein_spec %>%
  filter(Protein %in% check_increased_pept$protein_id)

inter_scater$Gene[duplicated(inter_scater$Gene)] <- NA

#geom_point(data = dataarg %>% filter(logFC < -FC_cutoff,
#                                           adj.P.Val < pval_cutoff),
#                 mapping = aes(x = logFC, y = -log10(adj.P.Val)), color = "red") +
#      geom_point(data = dataarg %>% filter(logFC > FC_cutoff,
#                                           adj.P.Val > pval_cutoff),
#                 mapping = aes(x = logFC, y = -log10(adj.P.Val)), color = "#2a9d8f") +

scatter_proteolysis_plot <- ggplot(log2semipept2_log2protein_spec, 
                                   aes(x = logFC_fully_tryp_protein, 
                                       y = logFC_semi_peptide)) + 
  geom_smooth(method=lm, se = FALSE, linetype="dashed", size = 0.5, 
              color = "black") + 
  geom_point(aes(color = DA_Protein, shape = DA_peptide), size = 0.6) +
  scale_color_manual(values = c("#2a9d8f", "red")) +
  xlab("log2(FC) - Protein abundances") + 
  ylab("log2(FC) - Semi-specific peptides") + 
  ggrepel::geom_text_repel(data = inter_scater,
                               aes(label = Gene), 
                           size = 1,
                           box.padding = 0.5,
                           max.overlaps = 25) +
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0, size = 8, angle = 360),
        axis.text.y = element_text(hjust = 0.95, vjust = 0.2, size = 8),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 8),
        legend.key.height= unit(3, 'mm'),
        legend.key.width= unit(3, 'mm'),
        legend.position="bottom") 

cor.test(log2semipept2_log2protein_spec$logFC_fully_tryp_protein,
    log2semipept2_log2protein_spec$logFC_semi_peptide, method = "pearson")
```

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  log2semipept2_log2protein_spec$logFC_fully_tryp_protein and log2semipept2_log2protein_spec$logFC_semi_peptide
    ## t = 16.098, df = 773, p-value < 2.2e-16
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  0.4464034 0.5520192
    ## sample estimates:
    ##       cor 
    ## 0.5010747

``` r
ggsave(plot = scatter_proteolysis_plot, 
       filename = here::here("figures/scatter_proteolysis_plot_gbm.tiff"), 
       device = "tiff",
       units = "mm",
       width = 63,
       height = 63)

ggsave(plot = scatter_proteolysis_plot, 
       filename = here::here("figures/scatter_proteolysis_plot_gbm.eps"), 
       device = "eps",
       units = "mm",
       width = 63,
       height = 63)
```

## Analysis of differential amino acid usage (iceLogo)

``` r
source(here("scr/get_cleavage_area.R"))
```

``` r
cleave_areas <- get_cleave_area(cleavage_annoated_peptides)
```

``` r
library(dagLogo)
```

**Load background proteome**

Used the identified protein sequences as background proteome.

``` r
proteome_ided <- prepareProteome(fasta = here("data/semi_specific_search_fragpipe17/protein_combined.fas"), 
                                 species = "Homo sapiens")
```

``` r
increased_Rec <- filter(cleave_areas$cleave_area20,
                            Peptide %in% increased_recwo6_pept) 

decreased_Rec <- filter(cleave_areas$cleave_area20,
                            Peptide %in% decreased_recwo6_pept)

write_tsv(x = increased_Rec,
          file = here("results/semi_tryptic/proteolytic_products_increased_in_Recurrent.tsv"))

write_tsv(x = decreased_Rec,
          file = here("results/semi_tryptic/proteolytic_products_decreased_in_Recurrent.tsv"))


increased_Rec_4ice <- filter(cleave_areas$cleave_area20,
                            Peptide %in% increased_recwo6_pept) %>%
                    pull(cleave_area20)

decreased_Rec_4ice <- filter(cleave_areas$cleave_area20,
                            Peptide %in% decreased_recwo6_pept) %>%
                    pull(cleave_area20)
```

**Format peptide sequences**

``` r
if(!file.exists(here("results/semi_tryptic/formated_pept_daglogo_increasedREC.rds"))){
  
  form_peptidesincreased_4ice <- formatSequence(increased_Rec_4ice, 
                                  proteome = proteome_ided)

  write_rds(form_peptidesincreased_4ice, file = here("results/semi_tryptic/formated_pept_daglogo_increasedREC.rds"))
} else {
  
  form_peptidesincreased_4ice <- read_rds(here("results/semi_tryptic/formated_pept_daglogo_increasedREC.rds"))
  
}
```

``` r
if(!file.exists(here("results/semi_tryptic/formated_pept_daglogo_decreasedREC.rds"))){
  
  form_peptidesdecreased_4ice <- formatSequence(decreased_Rec_4ice, 
                                  proteome = proteome_ided)

  write_rds(form_peptidesdecreased_4ice, file = here("results/semi_tryptic/formated_pept_daglogo_decreasedREC.rds"))
} else {
  
  form_peptidesdecreased_4ice <- read_rds(here("results/semi_tryptic/formated_pept_daglogo_decreasedREC.rds"))
  
}
```

#### DAU increased

``` r
bg_mod_ztest_increased <- buildBackgroundModel(form_peptidesincreased_4ice,
                                           proteome = proteome_ided,
                           background = "wholeProteome",
                           testType = "ztest")
```

``` r
dau_nogroup_increased <- testDAU(form_peptidesincreased_4ice, 
                           dagBackground = bg_mod_ztest_increased)
```

``` r
dagHeatmap(dau_nogroup_increased) 
```

![](gbm_proteomics_data_analysis_proteolytic_products_prep_figures_v3_gh_format_files/figure-gfm/unnamed-chunk-136-1.png)<!-- -->

#### DAU decreased

``` r
bg_mod_ztest_decreased <- buildBackgroundModel(form_peptidesdecreased_4ice,
                                           proteome = proteome_ided,
                           background = "wholeProteome",
                           testType = "ztest")
```

``` r
dau_nogroup_decreased <- testDAU(form_peptidesdecreased_4ice, 
                           dagBackground = bg_mod_ztest_decreased)
```

``` r
dagHeatmap(dau_nogroup_decreased) 
```

![](gbm_proteomics_data_analysis_proteolytic_products_prep_figures_v3_gh_format_files/figure-gfm/unnamed-chunk-139-1.png)<!-- -->

# Proteogenomics

## Initial data loading and wrangling

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
  mutate(recurrence = str_extract(paired_id_semi, "[^0-9]+$")) %>%
  mutate(patient = str_remove_all(paired_id_semi, "[^0-9]+$"))

long_mat_progen_2 <- long_mat_progen %>%
  mutate(paired_id = paste(patient, recurrence, sep = "_")) 

paired_annot_uniq_progen <- long_mat_progen_2 %>%
  dplyr::select(patient, recurrence, paired_id) %>%
  distinct()

# wide peptide abundace matrix with NAs
expr_matrix_pept_progen <- dplyr::select(tmt_peptdata_frag_progen,
                             index, peptide, starts_with("x")) %>%
  dplyr::select(!starts_with("x6"))

expr_matrixnona_pept_progen <- expr_matrix_pept_progen %>%
  na.omit()

long_mat_pept_progen <- pivot_longer(expr_matrix_pept_progen,
                         cols = starts_with("x"),
                         names_to = c("paired_id_semi"),
                         values_to = "Abundance") %>%
  mutate(recurrence = str_extract(paired_id_semi, "[^0-9]+$")) %>%
  mutate(patient = str_remove_all(paired_id_semi, "[^0-9]+$"))

long_mat_pept_progen_2 <- long_mat_pept_progen %>%
  mutate(paired_id = paste(patient, recurrence, sep = "_")) 

paired_annot_uniq_pept_progen <- long_mat_pept_progen_2 %>%
  dplyr::select(patient, recurrence, paired_id) %>%
  distinct()

## Annotated abundance data in long format ----

quant_annotated_progen <- left_join(long_mat_progen_2,
                             sample_annotation2) %>%
  mutate(Channel_mix = paste(mixture,channel,  sep = "_"))

quant_annotated_pept_progen <- left_join(long_mat_pept_progen_2,
                             sample_annotation2) %>%
  mutate(Channel_mix = paste(mixture,channel,  sep = "_"))


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
                   into = c("sp", "Genes", "ensembl"),
                   sep = "\\|") %>%
          dplyr::select(Genes, Peptide) %>%
          distinct() %>%
          mutate(Peptide = str_remove_all(Peptide, "\\[.+?\\]")) %>%
          mutate(Peptide = str_remove_all(Peptide, "^n"))

# peptide to protein mappong from tmt abundance report (quantified and normalized)
protein2peptide_tmtrept <- tmt_peptdata_frag_progen %>%
  dplyr::select(index, protein_id, peptide)
```

## Localization of identified peptides and mapping to variant regions

``` r
# list of identified peptides associated with 'variant' proteins
var_prots <- protein2peptideall %>%
                    mutate(variant = if_else(str_detect(Genes, "\\_"),
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

### Localized peptides within protein sequence and confirm identification of SAAVs

``` r
# tabulate location of the variant
annotated_w_variant_location <- cleavage_annoated_peptides %>%
          as_tibble() %>%
          mutate(protein_id_var = protein_id,
                 variant_type = if_else(str_detect(protein_id, "\\:"), 
                                        true = "other",
                                        false = "saav")) %>%
          mutate(protein_id = str_replace_all(protein_id, "\\.", "\\_")) %>%
          mutate(protein_id = str_replace_all(protein_id, "\\:", "\\_")) %>%
          separate(protein_id, 
                   into = c("ense_id", "variant_location"), 
                   sep = "\\_",
                   extra = "merge") %>%
                    na.omit() %>%
          mutate(variant_type = if_else(condition = str_detect(variant_location, "[0-9]$"),
                                        true = "other",
                                        false = variant_type))

# generate a long list of variant locations 

locat2prot <- str_split(annotated_w_variant_location$variant_location,
                        pattern = "\\_")

names(locat2prot) <- annotated_w_variant_location$ense_id %>%
                    na.omit()

locat2protdf <- list2DF(locat2prot)

long2prot <- locat2protdf %>%
          pivot_longer(everything(),
                       values_to = "variant_location",
                       names_to = "ense_id") %>% 
          distinct() %>%
          mutate(variant_num_location = readr::parse_number(variant_location))

long2prot_worivar <- left_join(long2prot, 
                               dplyr::select(annotated_w_variant_location,
                                             protein_id_var, ense_id)) %>%
  distinct()


annotated_wo_variant_location <- annotated_w_variant_location %>%
          dplyr::select(-variant_location)

# merge long table of variant locations with the annotation of peptide position within 
# their associated protein sequences

protein2varlocat <- left_join(annotated_wo_variant_location, long2prot_worivar) %>%
          distinct() %>%
          mutate(start_num = as.numeric(start_position),
                 end_num = as.numeric(end_position)) %>%
          mutate(variant_peptide = variant_num_location >= as.numeric(start_position) & variant_num_location <= as.numeric(end_position)) %>% # evaluate if the falls on the position of a variant sequence
          dplyr::select(protein_id_var, ense_id, protein_description, Peptide, start_position, 
                        start_num,end_position,end_num, variant_location, 
                        variant_num_location, variant_peptide, variant_type)

pept2variant <- protein2varlocat %>%
  dplyr::select(Peptide, variant_peptide, variant_type) %>%
  distinct() %>%
  na.omit() %>%
  filter(variant_type == "saav")
```

``` r
table(pept2variant$variant_peptide)
```

    ## 
    ## FALSE  TRUE 
    ## 16699   232

232 peptides were identified as mapping to single amino acid variant
regions of their corresponding annotated proteins sequence.

## Visualize quantitative features of identified SAAVs

### Prep abundance matrices

``` r
pep_prematwo6 <- dplyr::select(expr_matrix_pept_progen,
                               index, peptide, rownames(design_limmawo6))  

pep_matwo6 <- dplyr::select(expr_matrix_pept_progen,
                             index, rownames(design_limmawo6))  %>%
  column_to_rownames("index") %>%
  as.matrix()

pep_matnona_wo6 <- pep_matwo6 %>%
  na.omit()
```

### Filter, merge and summarize SAAV abundances

``` r
# extract indexing, peptide and sequence information from normalized data
index2genepeptide <- tmt_peptdata_frag_progen %>% 
  dplyr::select(index, peptide, gene, protein_id)

# table of interesting features (saavs)
features1 <- pept2variant %>%
  dplyr::rename(peptide = Peptide) %>%
  left_join(., index2genepeptide) %>%
  na.omit()

interesting_features1 <- features1 %>%
                    filter(variant_peptide == TRUE) 

# data frame of identified SAAVs and their annotation
protein2varlocat_saav <- protein2varlocat %>%
  dplyr::rename(peptide = Peptide) %>%
  filter(peptide %in% interesting_features1$peptide,
         variant_peptide == TRUE) %>%
  left_join(., index2genepeptide) %>%
  mutate(gene_variant = paste(gene, variant_location, sep = "_")) 

# mapping variant to index id
genevariant2index <- protein2varlocat_saav %>%
  dplyr::select(index, gene_variant)

# abundance matrix of peptides mapping to variant locations
pep_premat_variantwo6 <- filter(pep_prematwo6,
                                peptide %in% protein2varlocat_saav$peptide) %>%
  left_join(., genevariant2index) %>%
  relocate(gene_variant)

# summarize SAAVs (gene+variant location)
# i.e. take the median abundance of all the peptides mapping to the same SAAV

pep_premat_sumgenarwo6 <- pep_premat_variantwo6 %>%
  group_by(gene_variant) %>%
  summarise_at(vars(starts_with("x")), .funs = function(x){median(x, na.rm = TRUE)}) %>%
  ungroup() 

full_info_sumgenvar <- pep_premat_sumgenarwo6 %>%
  left_join(., protein2varlocat_saav) %>%
  filter(variant_type == "saav")

write_tsv(protein2varlocat_saav,
          here("results/proteogenomics/cellreps_2021/all_peptides_mapping_to_saavs_gbm.tsv"))

write_tsv(pep_premat_sumgenarwo6,
          here("results/proteogenomics/cellreps_2021/summarized_median_abundance_of_saavs.tsv"))

write_tsv(full_info_sumgenvar,
          here("results/proteogenomics/cellreps_2021/summarized_median_abundance_of_saavs_full.tsv"))

pep_matvariantswo6or <- pep_premat_sumgenarwo6 %>%
  column_to_rownames("gene_variant") %>%
  as.matrix()

pep_matvariantswo6 <- pep_premat_sumgenarwo6 %>%
  filter(str_detect(gene_variant, "CAP1_", 
                    negate = TRUE)) %>%
  column_to_rownames("gene_variant") %>%
  as.matrix()

pep_matvariantswo6nona <- pep_matvariantswo6 %>%
  na.omit()
```

The variant-containing peptides identified are able to identify 214
single amino acid mutations (SAAVs).

Among these identified SAAVs, 30 were identified in all mixtures and
samples.

### Heatmap of all SAAVs observed

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
  str_replace(., pattern = "prim", 
              replacement = "Init") %>%
  str_replace(., pattern = "rec", 
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
col_fun_or = colorRamp2(c(min(pep_matvariantswo6or_hm, na.rm = TRUE), 
                       median(pep_matvariantswo6or_hm, na.rm = TRUE), 
                       max(pep_matvariantswo6or_hm, na.rm = TRUE)), 
                     c("#2a9d8f", "white", "red"))

Heatmap(pep_matvariantswo6or_hm, 
        col = col_fun_or,
        column_split = stage, 
        row_names_gp = gpar(fontsize = 5.5),
        column_names_gp = gpar(fontsize = 7),
        cluster_rows = FALSE, cluster_columns = FALSE,
        heatmap_legend_param = list(
        title = "log2-\nNormalized\nAbundance"
    ))
```

![](gbm_proteomics_data_analysis_proteolytic_products_prep_figures_v3_gh_format_files/figure-gfm/unnamed-chunk-152-1.png)<!-- -->

### Heatmap of SAAVs observed on every sample

``` r
col_fun_sel = colorRamp2(c(min(pep_matvariantswo6or_hm, na.rm = TRUE), 
                       median(pep_matvariantswo6or_hm, na.rm = TRUE), 
                       max(pep_matvariantswo6or_hm, na.rm = TRUE)), 
                     c("#2a9d8f", "white", "red"))


Heatmap(pep_matvariantswo6nona_hm, 
        col = col_fun_sel,
        column_split = stage, 
        row_names_gp = gpar(fontsize = 7),
        column_names_gp = gpar(fontsize = 7),
        heatmap_legend_param = list(
        title = "log2-\nNormalized\nAbundance"
    ))
```

![](gbm_proteomics_data_analysis_proteolytic_products_prep_figures_v3_gh_format_files/figure-gfm/unnamed-chunk-153-1.png)<!-- -->

### Limma on SAAVs

### Prep design matrix

``` r
design_limmawo6_saav <- model.matrix(~patient+stage)

sample_annotationwo6 <- sample_annotationwo6 %>%
  mutate(paired_id_semi = paste(patient, recurrence, sep = ""))

rownames(design_limmawo6_saav) <- sample_annotationwo6$paired_id_semi
```

### Prep abundance matrices

``` r
colnames(pep_matvariantswo6nona)
```

    ##  [1] "x1prim"  "x1rec"   "x2prim"  "x2rec"   "x3prim"  "x3rec"   "x4prim" 
    ##  [8] "x4rec"   "x5prim"  "x5rec"   "x7prim"  "x7rec"   "x8prim"  "x8rec"  
    ## [15] "x9prim"  "x9rec"   "x10prim" "x10rec"  "x11prim" "x11rec"

``` r
rownames(design_limmawo6_saav)
```

    ##  [1] "x1prim"  "x1rec"   "x2prim"  "x2rec"   "x3prim"  "x3rec"   "x4prim" 
    ##  [8] "x4rec"   "x5prim"  "x5rec"   "x7prim"  "x7rec"   "x8prim"  "x8rec"  
    ## [15] "x9prim"  "x9rec"   "x10prim" "x10rec"  "x11prim" "x11rec"

#### SAAV-level model fit

``` r
fit_limmapeptwo6 <- function(mat, design, method,Limma){
  limmafit <- lmFit(mat, design, method = method)
  limmafit <- eBayes(limmafit)
  limma_tab <- topTable(limmafit, coef = "stageRecurrent", number = Inf, adjust.method = "BH") %>%
  mutate(Protein = rownames(.),
         Limma = Limma)
  
  return(limma_tab)
}
```

``` r
fit_mod_saavs <- fit_limmapeptwo6(pep_matvariantswo6nona, 
                                  design_limmawo6_saav, 
                                  method = 'robust', 
                                  Limma = "Robust - w Patient effect") 

compar_tab_saavs <-  fit_mod_saavs %>% 
  rownames_to_column("gene_variant") 
```

``` r
increased_recwo6_saavs <- compar_tab_saavs %>%
                    filter(logFC>0,
                           adj.P.Val < 0.05) %>% 
  pull(gene_variant)

decreased_recwo6_saavs <- compar_tab_saavs %>%
                    filter(logFC<0,
                           adj.P.Val < 0.05) %>% 
  pull(gene_variant)
```

``` r
length(increased_recwo6_saavs)
```

    ## [1] 1

``` r
length(decreased_recwo6_saavs)
```

    ## [1] 0

#### Volcano SAAVs

``` r
diff_saavs <- compar_tab_saavs %>%
                    filter(adj.P.Val < 0.05)

size <- 1.5

volcano_saavs <- ggplot(data = compar_tab_saavs,
                      mapping = aes(x = logFC, y = -log10(adj.P.Val))) +
      geom_point(data = compar_tab_saavs %>% filter(logFC > 0,
                                           adj.P.Val < 0.05),
                 mapping = aes(x = logFC, y = -log10(adj.P.Val)), 
                 color = "red",
                 size = size)+
      geom_point(data = compar_tab_saavs %>% filter(logFC < -0,
                                           adj.P.Val < 0.05),
                 mapping = aes(x = logFC, y = -log10(adj.P.Val)), color = "red",
                 size = size) +
      geom_point(data = compar_tab_saavs %>% filter(logFC > 0,
                                           adj.P.Val > 0.05),
                 mapping = aes(x = logFC, y = -log10(adj.P.Val)), color = "#2a9d8f",
                 size = size) +
      geom_point(data = compar_tab_saavs %>% filter(logFC < -0,
                                           adj.P.Val > 0.05),
                 mapping = aes(x = logFC, y = -log10(adj.P.Val)), color = "#2a9d8f",
                 size = size) +
      geom_hline(yintercept = -log10(0.05),
                 color = "red", linetype = "dashed") +
      ggrepel::geom_text_repel(data = diff_saavs,
                               aes(label = gene_variant), size = 4) +
      xlab("logFC - Recurrent / Primary") + 
      theme(axis.text.x = element_text(hjust = 0.5, vjust = 0, size = 6, angle = 360),
        axis.text.y = element_text(hjust = 0.95, vjust = 0.2, size = 8),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 8),
        legend.key.height= unit(3, 'mm'),
        legend.key.width= unit(3, 'mm'),
        legend.position="bottom") 
```

### Boxplot of differentially abundant SAAV

``` r
prelongsaav <- pep_matvariantswo6nona %>%
  as.data.frame() %>%
  rownames_to_column("gene_variant") %>%
   pivot_longer(cols = starts_with("x"),
                         names_to = c("patient"),
                         values_to = "Abundance") %>%
  mutate(Stage = str_extract(patient, "[^0-9]+$")) %>%
  mutate(Patient = str_remove_all(patient, "[^0-9]+$")) %>%
  filter(gene_variant == increased_recwo6_saavs)
```

``` r
diff_abund_saav <- ggplot(prelongsaav, 
                    aes(x = Stage, 
                        y = Abundance, fill = Stage, 
                        cex.axis = 1.5)) +
  geom_boxplot() +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.2) +
  # Box plot with jittered points
  # 0.2 : degree of jitter in x direction
  # geom_jitter(shape=16, position=position_jitter(0.2))
  ylab("log2-Normalized Abundance") +
  geom_signif(
    comparisons = list(c("Initial", "Recurrence")),
    map_signif_level = TRUE
  ) + 
  stat_compare_means(method="t.test") +
  ggtitle(increased_recwo6_saavs) +
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0, size = 6, angle = 360),
        axis.text.y = element_text(hjust = 0.95, vjust = 0.2, size = 8),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 8),
        legend.key.height= unit(3, 'mm'),
        legend.key.width= unit(3, 'mm'),
        legend.position="bottom") 
```

### Boxplot DA SAAV

``` r
  slim_data2saav <- pivot_wider(prelongsaav %>% dplyr::select(-patient), 
                          values_from = Abundance, 
                          names_from = Stage) #%>%
  #dplyr::select(-patient)

  slim_data3saav <- slim_data2saav %>% 
    mutate(Abs_diff = rec-prim) %>%
    mutate(rank = row_number(Abs_diff)) %>% 
    dplyr::rename(Initial = prim, Recurrent = rec) %>%
    ungroup() 
```

``` r
ggpaired(slim_data3saav,
         cond1 = "Initial",
         cond2 = "Recurrent", 
         fill = "condition",
         ylab = "log2-Normalized Abundance",
         label = NULL,
         repel = TRUE, point.size = 0.2, line.size = 0.1) +
      theme(axis.text.x = element_text(hjust = 0.5, vjust = 0, size = 5, angle = 360),
            axis.text.y = element_text(hjust = 0.95, vjust = 0.2, size = 5),
            panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=0.5),
            axis.title = element_text(size = 8),
            legend.text = element_text(size = 6),
            legend.title = element_text(size = 8),
            legend.key.height= unit(3, 'mm'),
            legend.key.width= unit(3, 'mm'),
            legend.position="bottom",
            strip.text = element_text(hjust = 0.5, vjust = 0, size = 8, angle = 360)) + 
  ggtitle(increased_recwo6_saavs)
```

![](gbm_proteomics_data_analysis_proteolytic_products_prep_figures_v3_gh_format_files/figure-gfm/unnamed-chunk-167-1.png)<!-- -->

``` r
paired_box_saavs <- ggpaired(slim_data3saav,
         cond1 = "Initial",
         cond2 = "Recurrent", 
         fill = "condition",
         ylab = "log2-Normalized Abundance",
         label = NULL,
         repel = TRUE, point.size = 0.2, line.size = 0.1) +
      theme(axis.text.x = element_text(hjust = 0.5, vjust = 0, size = 5, angle = 360),
            axis.text.y = element_text(hjust = 0.95, vjust = 0.2, size = 5),
            panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=0.5),
            axis.title = element_text(size = 8),
            legend.text = element_text(size = 6),
            legend.title = element_text(size = 8),
            legend.key.height= unit(3, 'mm'),
            legend.key.width= unit(3, 'mm'),
            legend.position="bottom",
            strip.text = element_text(hjust = 0.5, vjust = 0, size = 8, angle = 360)) 
```

``` r
sessionInfo()
```

    ## R version 4.1.0 (2021-05-18)
    ## Platform: x86_64-w64-mingw32/x64 (64-bit)
    ## Running under: Windows 10 x64 (build 19043)
    ## 
    ## Matrix products: default
    ## 
    ## locale:
    ## [1] LC_COLLATE=English_United States.1252 
    ## [2] LC_CTYPE=English_United States.1252   
    ## [3] LC_MONETARY=English_United States.1252
    ## [4] LC_NUMERIC=C                          
    ## [5] LC_TIME=English_United States.1252    
    ## 
    ## attached base packages:
    ##  [1] grid      parallel  stats4    stats     graphics  grDevices utils    
    ##  [8] datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] ComplexHeatmap_2.8.0  circlize_0.4.13       dagLogo_1.30.0       
    ##  [4] ggpubr_0.4.0          seqinr_4.2-8          drawProteins_1.12.0  
    ##  [7] janitor_2.1.0         here_1.0.1            DT_0.20              
    ## [10] org.Hs.eg.db_3.13.0   AnnotationDbi_1.54.1  IRanges_2.26.0       
    ## [13] S4Vectors_0.30.0      Biobase_2.52.0        BiocGenerics_0.38.0  
    ## [16] ReactomePA_1.36.0     clusterProfiler_4.0.5 missForest_1.4       
    ## [19] itertools_0.1-3       iterators_1.0.13      foreach_1.5.1        
    ## [22] randomForest_4.6-14   naniar_0.6.1          limma_3.48.3         
    ## [25] sva_3.40.0            BiocParallel_1.26.2   genefilter_1.74.0    
    ## [28] mgcv_1.8-36           nlme_3.1-152          kableExtra_1.3.4     
    ## [31] fs_1.5.0              mixOmics_6.16.3       lattice_0.20-44      
    ## [34] MASS_7.3-54           forcats_0.5.1         stringr_1.4.0        
    ## [37] dplyr_1.0.7           purrr_0.3.4           readr_2.0.1          
    ## [40] tidyr_1.1.3           tibble_3.1.4          ggplot2_3.3.5        
    ## [43] tidyverse_1.3.1      
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] rappdirs_0.3.3         visdat_0.5.3           bit64_4.0.5           
    ##   [4] knitr_1.36             UniProt.ws_2.32.0      data.table_1.14.0     
    ##   [7] KEGGREST_1.32.0        RCurl_1.98-1.4         doParallel_1.0.16     
    ##  [10] generics_0.1.1         cowplot_1.1.1          RSQLite_2.2.8         
    ##  [13] shadowtext_0.0.9       bit_4.0.4              tzdb_0.1.2            
    ##  [16] enrichplot_1.12.3      webshot_0.5.2          xml2_1.3.2            
    ##  [19] lubridate_1.7.10       assertthat_0.2.1       viridis_0.6.2         
    ##  [22] xfun_0.25              hms_1.1.1              evaluate_0.14         
    ##  [25] fansi_0.5.0            progress_1.2.2         dbplyr_2.1.1          
    ##  [28] readxl_1.3.1           igraph_1.2.6           DBI_1.1.1             
    ##  [31] htmlwidgets_1.5.4      rARPACK_0.11-0         ellipsis_0.3.2        
    ##  [34] RSpectra_0.16-0        backports_1.2.1        annotate_1.70.0       
    ##  [37] biomaRt_2.48.3         vctrs_0.3.8            Cairo_1.5-12.2        
    ##  [40] abind_1.4-5            cachem_1.0.6           withr_2.4.2           
    ##  [43] ggforce_0.3.3          checkmate_2.0.0        vroom_1.5.4           
    ##  [46] treeio_1.16.2          prettyunits_1.1.1      cluster_2.1.2         
    ##  [49] svglite_2.0.0          DOSE_3.18.3            ape_5.5               
    ##  [52] lazyeval_0.2.2         crayon_1.4.2           ellipse_0.4.2         
    ##  [55] edgeR_3.34.0           pkgconfig_2.0.3        labeling_0.4.2        
    ##  [58] tweenr_1.0.2           GenomeInfoDb_1.28.4    rlang_0.4.11          
    ##  [61] lifecycle_1.0.1        downloader_0.4         filelock_1.0.2        
    ##  [64] BiocFileCache_2.0.0    modelr_0.1.8           cellranger_1.1.0      
    ##  [67] rprojroot_2.0.2        polyclip_1.10-0        matrixStats_0.60.1    
    ##  [70] graph_1.70.0           Matrix_1.3-4           aplot_0.1.1           
    ##  [73] carData_3.0-4          reprex_2.0.1           GlobalOptions_0.1.2   
    ##  [76] pheatmap_1.0.12        png_0.1-7              viridisLite_0.4.0     
    ##  [79] rjson_0.2.20           bitops_1.0-7           Biostrings_2.60.2     
    ##  [82] blob_1.2.2             shape_1.4.6            qvalue_2.24.0         
    ##  [85] rstatix_0.7.0          gridGraphics_0.5-1     ggsignif_0.6.3        
    ##  [88] reactome.db_1.76.0     scales_1.1.1           memoise_2.0.0         
    ##  [91] graphite_1.38.0        magrittr_2.0.1         plyr_1.8.6            
    ##  [94] zlibbioc_1.38.0        compiler_4.1.0         scatterpie_0.1.7      
    ##  [97] RColorBrewer_1.1-2     clue_0.3-60            snakecase_0.11.0      
    ## [100] cli_3.1.0              ade4_1.7-17            XVector_0.32.0        
    ## [103] patchwork_1.1.1        tidyselect_1.1.1       stringi_1.7.4         
    ## [106] highr_0.9              yaml_2.2.1             GOSemSim_2.18.1       
    ## [109] locfit_1.5-9.4         ggrepel_0.9.1          fastmatch_1.1-3       
    ## [112] tools_4.1.0            rstudioapi_0.13        gridExtra_2.3         
    ## [115] farver_2.1.0           ggraph_2.0.5           digest_0.6.27         
    ## [118] Rcpp_1.0.7             car_3.0-12             broom_0.7.10          
    ## [121] httr_1.4.2             motifStack_1.36.1      colorspace_2.0-2      
    ## [124] rvest_1.0.2            XML_3.99-0.7           splines_4.1.0         
    ## [127] yulab.utils_0.0.4      tidytree_0.3.6         graphlayouts_0.7.1    
    ## [130] ggplotify_0.1.0        systemfonts_1.0.2      xtable_1.8-4          
    ## [133] jsonlite_1.7.2         ggtree_3.0.4           tidygraph_1.2.0       
    ## [136] corpcor_1.6.10         ggfun_0.0.4            R6_2.5.1              
    ## [139] pillar_1.6.4           htmltools_0.5.2        glue_1.4.2            
    ## [142] fastmap_1.1.0          codetools_0.2-18       fgsea_1.18.0          
    ## [145] utf8_1.2.2             curl_4.3.2             magick_2.7.3          
    ## [148] GO.db_3.13.0           survival_3.2-13        rmarkdown_2.11        
    ## [151] munsell_0.5.0          DO.db_2.9              GetoptLong_1.0.5      
    ## [154] GenomeInfoDbData_1.2.6 haven_2.4.3            reshape2_1.4.4        
    ## [157] gtable_0.3.0
