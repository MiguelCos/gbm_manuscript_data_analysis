## function for feature-specific FDR correction
## used specifically in this case to apply FDR correction only on semi-specific peptides

feature_fdr_correction <- function(toptable,
                                   features_table,
                                   interesting_features_table){
  
  # merge limma output with feature annotation  ----
  
  tab_limma_feature_annot <- left_join(toptable, 
                                       features_table,
                                       by = c("peptide","index"))
  
  # filter limma table based on interesting feature ----
  
  tab_limma_subsetin <- filter(tab_limma_feature_annot,
                               peptide %in% interesting_features_table$peptide) %>%
    # apply FDR correction only on the subset
    mutate(adj.P.Val = p.adjust(p = P.Value, 
                                method = "BH")) %>%
    mutate(fdr_correction = 'feature-specific')
  
  # create a table of non-interesting features containing the globally adjusted p-values
  tablimma_subsetout <- filter(tab_limma_feature_annot,
                               !peptide %in% interesting_features_table$peptide) %>%
    mutate(fdr_correction = 'global')
  
  output_limma3 <- bind_rows(tab_limma_subsetin,
                             tablimma_subsetout) 
  
  return(output_limma3)
  
}