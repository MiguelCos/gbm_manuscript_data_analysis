# Function to fit the linear model with limma, generate the summary table and merge 
# with gene annotation per protein

fit_limmawo6 <- function(mat, design, method,Limma, prot2gene){
  limmafit <- lmFit(mat, design, method = method)
  limmafit <- eBayes(limmafit)
  limma_tab <- topTable(limmafit, coef = "recurrencewo6rec", number = Inf, adjust.method = "BH") %>%
    mutate(Protein = rownames(.),
           Limma = Limma)
  
  limma_tab <- left_join(limma_tab, prot2gene)
  
  return(limma_tab)
}

