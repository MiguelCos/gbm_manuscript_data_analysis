# function to summarize peptide counts per annotated feature in terms of
# N-terminal modification and proteolytic specificity

summarize_peptide_counts <- function(tmt_report_annotated){
  
  require(dplyr)
  require(tibble)
  
  to_count_info <- tmt_report_annotated %>% 
    dplyr::select(peptide, specificity, nterm, semi_type, tmt_tag)
  
  to_count_info_acetyl <- to_count_info %>% 
    dplyr::filter(nterm == "acetylated")
  
  n_semi <- dplyr::count(to_count_info, specificity) %>% 
    dplyr::rename(feature_type = specificity) %>%
    dplyr::mutate(category = "specificity")
  
  n_term <- dplyr::count(to_count_info, nterm) %>% 
    dplyr::rename(feature_type = nterm) %>%
    dplyr::mutate(category = "N-term")
  
  n_semi_type <- dplyr::count(to_count_info, semi_type) %>% 
    dplyr::rename(feature_type = semi_type) %>%
    dplyr::mutate(category = "Semi type")
  
  n_tmt_tag <- dplyr::count(to_count_info, tmt_tag) %>% 
    dplyr::rename(feature_type = tmt_tag)  %>%
    dplyr::mutate(category = "TMT location")
  
  n_total <- tibble(feature_type = "Total",
                    n = nrow(to_count_info),
                    category = "Total")
  
  summary_count <- bind_rows(n_semi,
                             n_term,
                             n_semi_type,
                             n_tmt_tag,
                             n_total)
  
  return(summary_count)
  
}