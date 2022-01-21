# function to create volcano plots for visualization in the HTML report

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
