## function for feature-specific FDR correction
## used specifically in this case to apply FDR correction only on semi-specific peptides

feature_fdr_correction_mst <- function(toptable,
                                       features_table,
                                       interesting_features_table){

                    # merge limma output with feature annotation  ----

                    tab_limma_feature_annot <- left_join(toptable,
                                                         features_table,
                                                         by = c("peptide"))

                    # filter limma table based on interesting feature ----

                    tab_limma_subsetin <- filter(tab_limma_feature_annot,
                                                 peptide %in% interesting_features_table$peptide) %>%
                                        # apply FDR correction only on the subset
                                        mutate(adj.pvalue = p.adjust(p = pvalue,
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

get_da_peptides <- function(tabular_dea_result,
                            type) {

                    if (type == "limma"){

                                        increased <- tabular_dea_result %>%
                                                            filter(logFC > 0,
                                                                   adj.P.Val < 0.05) %>%
                                                            pull(peptide)

                                        decreased <- tabular_dea_result %>%
                                                            filter(logFC < 0,
                                                                   adj.P.Val < 0.05) %>%
                                                            pull(peptide)

                                        peptides_da <- list(increased = increased,
                                                            decreased = decreased)

                    } else if (type == "msstats"){

                                        increased <- tabular_dea_result %>%
                                                            filter(log2FC > 0,
                                                                   adj.pvalue < 0.05) %>%
                                                            pull(peptide)

                                        decreased <- tabular_dea_result %>%
                                                            filter(log2FC < 0,
                                                                   adj.pvalue < 0.05) %>%
                                                            pull(peptide)

                                        peptides_da <- list(increased = increased,
                                                            decreased = decreased)

                    } else if (type == "msstats_limma") {

                                        increased <- tabular_dea_result %>%
                                                            filter(adj.P.Val < 0.05,
                                                                   logFC > 0) %>%
                                                            pull(index)

                                        decreased <- tabular_dea_result %>%
                                                            filter(adj.P.Val < 0.05,
                                                                   logFC < 0) %>%
                                                            pull(index)

                                        peptides_da <- list(increased = increased,
                                                            decreased = decreased)
                    } else {

                                        errorCondition("Error: Check your type input")

                    }

                    return(peptides_da)
}
