## peptide coverage plots

get_coverage <- function(annotated_peptides, id){
  
  mappedid <- dplyr::filter(annotated_peptides,
                            protein_id == id) %>% 
    dplyr::mutate(type = "PEPTIDE",
                  begin = as.numeric(start_position),
                  end = as.numeric(end_position), 
                  length = as.numeric(peptide_length)) %>%
    dplyr::rename(description = Peptide,
                  accession = protein_id) %>% 
    dplyr::select(type, description, begin, end, length, accession) %>%
    dplyr::arrange(end) %>%
    dplyr::mutate(order = 1)
  
}

add_peptides <- function(feature_df,
                         peptide_coverage_data){
  
  taxid <- feature_df$taxid[1]
  entryName <- feature_df$entryName[1]
  
  pept_feat <- dplyr::mutate(peptide_coverage_data,
                             taxid = taxid,
                             entryName = entryName) %>%
    dplyr::select(names(feature_df)) # added cleave window information
  
  feature_df <- mutate(feature_df)
  
  wpept_features <- bind_rows(feature_df, pept_feat)
  
  return(wpept_features)
  
}

draw_peptides <- function(p, data){
  data$order <- 1.4
  p + ggplot2::geom_rect(data = data,
                         mapping=ggplot2::aes(xmin=begin,
                                              xmax=end,
                                              ymin=order-0.05,
                                              ymax=order+0.05),
                         colour = "black") +
    annotate("text", x = -50, y = 1.4, label = "Peptides") +
    theme_bw(base_size = 12) + # white background
    theme(panel.grid.minor=element_blank(),
          panel.grid.major=element_blank()) +
    theme(axis.ticks = element_blank(),
          axis.text.y = element_blank()) +
    theme(panel.border = element_blank())
}