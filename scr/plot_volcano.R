## function to generate a volcano plot

plot_volcano <- function(output_limma3){
  tovolc_all <- output_limma3 %>% 
    mutate(Differentially_expressed = case_when(adj.P.Val <= 0.05 ~ TRUE,
                                                TRUE ~ FALSE))
  
  tovolc_featu <- output_limma3 %>% 
    filter(fdr_correction == "feature-specific") %>%
    mutate(Differentially_expressed = case_when(adj.P.Val <= 0.05 ~ TRUE,
                                                TRUE ~ FALSE))
  volcanoes <- ggplot(data = tovolc_all, 
                      mapping = aes(x = logFC,
                                    y = -log10(adj.P.Val),
                                    color = Differentially_expressed,
                                    label = ID)) + 
    geom_point(alpha = 0.5) + 
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    #ggtitle("Volcano plot of all limma-tested features") + 
    theme(legend.position = "none")
  
  volcano_features <- ggplot(data = tovolc_featu, 
                             mapping = aes(x = logFC,
                                           y = -log10(adj.P.Val),
                                           color = Differentially_expressed,
                                           label = ID)) + 
    geom_point(alpha = 0.5) + 
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    #ggtitle("Volcano plot of selected/interesting limma-tested features") + 
    theme(legend.position = "none")
  
  print(volcano_features)
}