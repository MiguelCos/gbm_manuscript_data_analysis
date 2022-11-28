# Reproducible Data Analysis for the Manuscript 'Quantitative Proteomics Highlights Enrichment of key proteins in Recurrent Glioblastoma'

## Reproducible report:

The general data analysis and processing of the search results after peptide/protein identification and quantitation can be accessed via the [reproducible report](https://github.com/MiguelCos/gbm_manuscript_data_analysis/blob/main/scr/gbm_proteomics_data_analysis_proteolytic_products_prep_figures_v3_gh_format.md) in this repo. 

### Functions for peptide annotation and analysis of proteolytic processing  

We have written a series of functions for the annotation of peptides according to their location in the protein sequence and proteolytic specificity. These would support the quantitative analysis of products of proteolytic processing, specifically in the context of TMT-based isobaric quantitation.

#### [Annotate peptides](https://github.com/MiguelCos/gbm_manuscript_data_analysis/blob/main/scr/annotate_peptides.R): 

Used to map peptides to their corresponding protein in a fasta file and annotate them in terms of their position within the protein sequence, amino acids before and after, and specificity type. This function is also used to evaluate if an identified peptide contains a single amino acid variant, based on its called position.

#### [Annotate N-term](https://github.com/MiguelCos/gbm_manuscript_data_analysis/blob/main/scr/annotate_nterm.R):

Evaluates if the peptide N-term contains a TMT-tag, acetylation or if it is free. 

#### [Feature-specific FDR correction](https://github.com/MiguelCos/gbm_manuscript_data_analysis/blob/main/scr/feature_specific_fdr_correction.R):

After fitting the peptide-level linear model with limma, this function applies FDR correction on specific features of interest. In our case, we focus on semi-specific peptides. We based this approach on the notions of Bourgon et al. 2010 [10.1073/pnas.0914005107](https://doi.org/10.1073/pnas.0914005107).

#### [Get cleavage area](https://github.com/MiguelCos/gbm_manuscript_data_analysis/blob/main/scr/get_cleavage_area.R):

After identifying a set of interesting semi-specific peptides/proteolytic products, this function was used to generate peptide sequences that catch the residues in the vicinity (i.e. 10 amino acids after and before) of the non-tryptic cleavage area. This, as a preparation for the analysis of sequence motifs. 

### [Analysis of Plasma ELISA Proteomics](https://github.com/MiguelCos/gbm_manuscript_data_analysis/blob/main/results/elisa_plasma_proteomics/manuscript_elisa_figure_final.pdf)

To regenerate analysis for Plasma ELISA Proteomics, knit the [manuscript_elisa_figure_final.Rmd](https://github.com/MiguelCos/gbm_manuscript_data_analysis/blob/main/results/elisa_plasma_proteomics/manuscript_elisa_figure_final.Rmd) r notebook.