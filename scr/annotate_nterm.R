# function for annotation of N-terminal modifications

annotate_nterm <- function(peptidestsv, # peptide.tsv table
                           tmtmass = 304.2072, # either 304.2072 for 16plex or 229.1629 for 10/11plex
                           protease_specificity = "R|K") # or R for argc 
{
  
  require(dplyr)
  require(stringr)
  
  if (tmtmass == 304.2072){
    
    nterm_tmt <- "N-term\\(304.2072\\)"
    ktmt <- "K\\(304.2072\\)"
    
  } else if (tmtmass == 229.1629){
    
    nterm_tmt <- "N-term\\(229.1629\\)"
    ktmt <- "K\\(229.1629\\)"                    
    
  } else {
    
    error("Please check your tmtmass argument input.")
    
  }
  
  nterm_annotated <- peptidestsv %>% 
    mutate(nterm = case_when(str_detect(assigned_modifications, nterm_tmt) ~ "TMT-labelled",
                             str_detect(assigned_modifications, "N-term\\(42.0106\\)") ~ "acetylated",
                             TRUE ~ "free")) %>%
    mutate(tmt_tag = case_when(str_detect(assigned_modifications, nterm_tmt) ~ "nterm",
                               str_detect(assigned_modifications, ktmt) ~ "lysine",
                               str_detect(assigned_modifications, ktmt,
                                          negate = TRUE) & str_detect(assigned_modifications, "N-term\\(42.0106\\)") ~ "untagged_acetylated",
                               str_detect(assigned_modifications, ktmt,
                                          negate = TRUE) & nterm == "acetylated" ~ "untagged_acetylated",
                               str_detect(assigned_modifications, ktmt,
                                          negate = TRUE) & nterm == "free" ~ "untagged_free",
                               TRUE ~ "untagged"))
  
  return(nterm_annotated)
  
}