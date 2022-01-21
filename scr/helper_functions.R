### helper functions

## sample proteins for profile plots and QC plots ----
sample_proteins <- function(x, # char vector of protein IDs
                            size, # number of proteins to sample
                            seed = 363 # random seed number 
) {
  
  # select/sample N proteins to observe their abundance distribution  
  set.seed(seed)
  
  proteins <- x %>% 
    unique() %>%
    str_remove_all(pattern = "Biognosys") # get unique protein IDs
  
  which_prots <- sample(proteins, 
                        size = size, 
                        replace = FALSE)
  
}


## eliminate proteins based on missing values ----

na_filter <- function(data, max_na_per_group, logtrans, minInt_inpute, fraction_Inpute = 5){
  
  data_long <- pivot_longer(data = data,
                            cols = matches("^S"),
                            values_to = "Intensity",
                            names_to = c("Group-Sample")) %>%
    separate(col = "Group-Sample", 
             into = c("Group", "Sample"),
             sep = "\\-")
  
  if (logtrans){
    data_long <- mutate(data_long,
                        Intensity = log2(Intensity))
  }
  
  na_count <- group_by(data_long,
                       ID, Group) %>%
    summarise(na_count = sum(is.na(Intensity)),
              total = n()) %>% 
    ungroup() %>% 
    mutate(NA_fraction = na_count/total)
  
  nafraction <- group_by(na_count,
                         ID) %>%
    summarise(max_nafraction = max(NA_fraction),
              min_nafraction = min(NA_fraction)) %>%
    ungroup()
  
  included_prots <- dplyr::filter(nafraction,
                                  max_nafraction <= max_na_per_group)
  
  
  
  if (minInt_inpute){
    
    min_intensityprot <- data_long %>% 
      group_by(ID) %>% 
      mutate(min_int = min(Intensity, na.rm = TRUE)) %>%
      mutate(th_of_min_int = min_int/fraction_Inpute)
    
    data_long_nonas1 <- min_intensityprot %>% 
      mutate(Intensity = ifelse(is.na(Intensity),
                                yes = th_of_min_int,
                                no = Intensity),
             `Group-Sample` = paste(Group,"-",Sample, sep = "")) %>%
      dplyr::select(ID, `Group-Sample`, Intensity) %>%
      dplyr::filter(ID %in% included_prots$ID)
    
    data_maxperc_nas <- pivot_wider(data_long_nonas1,
                                    names_from = `Group-Sample`,
                                    values_from = Intensity)
  } else {
    data_long_nonas1 <- filter(data_long,
                               ID %in% included_prots$ID) %>% 
      mutate(`Group-Sample` = paste(Group,"-",Sample, sep = "")) %>%
      dplyr::select(ID, `Group-Sample`, Intensity)
    
    data_maxperc_nas <- pivot_wider(data_long_nonas1,
                                    names_from = `Group-Sample`,
                                    values_from = Intensity)
  }
  
  
}

### volcano plots ----  

our_volcano <- function(dataarg, 
                        FC_cutoff = 1, 
                        pval_cutoff = 0.05, 
                        color_diffex = "red", 
                        color_nondifex = "#2a9d8f", 
                        interesting_proteins = NULL, 
                        vert_line_col = "red",
                        hline_col = "red", 
                        hline_pos = 0.05, 
                        linetype = "dashed",
                        #vprot_ylim = c(0,2),
                        #vprot_xlim = c(-2,2),
                        increased_in, 
                        comparison_title) {
  
  signi_hits <- filter(dataarg, 
                       adj.P.Val <= pval_cutoff) %>% 
    pull(Protein)
  
  if(is_null(interesting_proteins)){
  
  protlabels <- dataarg %>%
    filter(adj.P.Val <= pval_cutoff)
  
  } else {
    
  protlabels <- dataarg %>%
      filter(Gene %in% interesting_proteins)
    
  }
  

    volcano <- ggplot(data = dataarg,
                      mapping = aes(x = logFC, y = -log10(adj.P.Val))) +
      geom_point()+
      geom_point(data = dataarg %>% filter(logFC > FC_cutoff,
                                           adj.P.Val < pval_cutoff),
                 mapping = aes(x = logFC, y = -log10(adj.P.Val)), color = "red")+
      geom_point(data = dataarg %>% filter(logFC < -FC_cutoff,
                                           adj.P.Val < pval_cutoff),
                 mapping = aes(x = logFC, y = -log10(adj.P.Val)), color = "red") +
      geom_point(data = dataarg %>% filter(logFC > FC_cutoff,
                                           adj.P.Val > pval_cutoff),
                 mapping = aes(x = logFC, y = -log10(adj.P.Val)), color = "#2a9d8f") +
      geom_point(data = dataarg %>% filter(logFC < -FC_cutoff,
                                           adj.P.Val > pval_cutoff),
                 mapping = aes(x = logFC, y = -log10(adj.P.Val)), color = "#2a9d8f") +
      ggrepel::geom_text_repel(data = protlabels,
                               aes(label = Gene)) +
      geom_hline(yintercept = -log10(hline_pos),
                 color = "red", linetype = "dashed") +
      ggtitle(paste0("Differentially abundant proteins:\n",comparison_title),
              paste0("Positive values = Increased in ", increased_in)) + 
      labs(caption = paste0("Number of significant hits: ",length(signi_hits)))+
      theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.1, size = 14),
            panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=1.5),
            axis.title=element_text(size=12,face="bold")) 
      
  return(volcano)
  
}
