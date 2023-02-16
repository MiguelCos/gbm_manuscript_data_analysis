# merge fasta of identified sequences from semi-specific search into a single fasta.

fasta1 <- seqinr::read.fasta(here("data/semi_specific_search_fragpipe17/mix_1/protein.fas"),
                             seqtype = "AA", 
                             as.string = TRUE)

fasta2 <- seqinr::read.fasta(here("data/semi_specific_search_fragpipe17/mix_2/protein.fas"),
                             seqtype = "AA", 
                             as.string = TRUE)

fasta3 <- seqinr::read.fasta(here("data/semi_specific_search_fragpipe17/mix_3/protein.fas"),
                             seqtype = "AA", 
                             as.string = TRUE)

fasta22 <- append(fasta1, 
                  fasta2)

fasta <- append(fasta22, 
                fasta3)

fasta <- discard(fasta,
                 duplicated(names(fasta)))

attr(fasta[[1]], "Annot")

extr_annt <- function(x){
  ext <- attr(x, "Annot")
}

names <- purrr::map_chr(fasta, 
                        .f = extr_annt)

seqinr::write.fasta(fasta,
                    file.out = here("data/semi_specific_search_fragpipe17/protein_combined.fas"),
                    names = names)