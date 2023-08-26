# Reproducible Data Analysis for the Manuscript 'Quantitative Proteomics Highlights Enrichment of key proteins in Recurrent Glioblastoma'

[![DOI](https://zenodo.org/badge/427062720.svg)](https://zenodo.org/badge/latestdoi/427062720) 

## Reproducible reports:

### General Proteomics analysis

The general data analysis and processing of the search results after peptide/protein identification and quantitation can be accessed via the [general proteomics reproducible report](https://github.com/MiguelCos/gbm_manuscript_data_analysis/blob/main/reports/gbm_recurrence_general_proteomics_manuscript_report.md) in this repo. 

### Large scale analyses of proteolytic processing

The reproducible report, containing code for data preprocessing, statisical analysis and intermediary plots for the large-scale differential analysis of proteolytic processing in recurrent glioblastoma can be found accessed via: [large-scalge proteolytic processing analysis reproducible report](https://github.com/MiguelCos/gbm_manuscript_data_analysis/blob/main/reports/gbm_recurrence_large_scale_proteolysis_manuscript_report.md)

### Integrative lipidomics and proteomics  

The reproducible report, containing code for data preprocessing, statisical analysis and intermediary plots for the large-scale differential analysis of proteolytic processing in recurrent glioblastoma can be found accessed via: [integrative lipidomics+proteomics reproducible report](https://github.com/MiguelCos/gbm_manuscript_data_analysis/blob/main/reports/gbm_recurrence_lipidomics_analysis.md)

### Proteogenomics analyses

The reproducible report, containing code for data preprocessing, statisical analysis and intermediary plots for the large-scale differential analysis of proteolytic processing in recurrent glioblastoma can be found accessed via: [proteogenomics reproducible report](https://github.com/MiguelCos/gbm_manuscript_data_analysis/blob/main/reports/gbm_recurrence_proteogenomics_manuscript_report.md)

### [Analysis of Plasma ELISA Proteomics](https://github.com/MiguelCos/gbm_manuscript_data_analysis/blob/main/results/elisa_plasma_proteomics/manuscript_elisa_figure_final.pdf)

To regenerate analysis for Plasma ELISA Proteomics, knit the [manuscript_elisa_figure_final.Rmd](https://github.com/MiguelCos/gbm_manuscript_data_analysis/blob/main/results/elisa_plasma_proteomics/manuscript_elisa_figure_final.Rmd) r notebook. Re-analysis of the Cox proportional hazards model including Age and Sex and covariates can be found [here](https://github.com/MiguelCos/gbm_manuscript_data_analysis/blob/main/reports/gbm_recurrence_cox_phm_incl_age_sex_revision.md)

### Single-cell RNAseq data mining

We mined the single-cell RNAseq dataset on recurrent GBM [GSM4972210](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4972210), to explore the expression of ASAH1 by annotated cell type. Reproducible report available [here](https://github.com/MiguelCos/gbm_manuscript_data_analysis/blob/main/reports/gbm_recurrence_mine_single_cell_data_visualization.md)

## Specific R functions for data processing and analysis

### Functions for peptide annotation and analysis of proteolytic processing  

We have written a series of functions for the annotation of peptides according to their location in the protein sequence and proteolytic specificity. These would support the quantitative analysis of products of proteolytic processing, specifically in the context of TMT-based isobaric quantitation.

#### [Annotate peptides](https://github.com/MiguelCos/gbm_manuscript_data_analysis/blob/main/scr/annotate_peptides.R): 

Used to map peptides to their corresponding protein in a fasta file and annotate them in terms of their position within the protein sequence, amino acids before and after, and specificity type. This function is also used to evaluate if an identified peptide contains a single amino acid variant, based on its called position.

#### [Annotate N-term](https://github.com/MiguelCos/gbm_manuscript_data_analysis/blob/main/scr/annotate_nterm.R):

Evaluates if the peptide N-term contains a TMT-tag, acetylation or if it is free. 

#### [Get cleavage area](https://github.com/MiguelCos/gbm_manuscript_data_analysis/blob/main/scr/get_cleavage_area.R):

After identifying a set of interesting semi-specific peptides/proteolytic products, this function was used to generate peptide sequences that catch the residues in the vicinity (i.e. 10 amino acids after and before) of the non-tryptic cleavage area. This, as a preparation for the analysis of sequence motifs. 

