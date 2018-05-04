#!/usr/bin/env Rscript
  
extract_psg_snp <- function(path_aln, path_results){
  # Check whether necessary libraries are available
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("dplyr package needed for this function to work. Please install it.", 
         call. = FALSE)
  }
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("Biostrings package needed for this function to work. Please install it.", 
         call. = FALSE)
  }
  
  # Import posigene result file
  result_file <- list.files(path_results, pattern = "results.tsv",
                            full.names = TRUE)
  results_posi <- read.table(result_file, 
                             header = TRUE, fill = TRUE, sep = "\t")
  
  # Extract positions in alignment file for each gene
  df_posi <- cbind(results_posi[c("Transcript")], 
                      apply(results_posi[, 27:ncol(results_posi)], 2, FUN = function(x) extract_aa(x, pos = c(10:11)))
                      )
  
  # Loop through all gene alignment files to extract SNPs for each PSG
  # Read prank best alignment file
  align_set <- Biostrings::readDNAStringSet(filepath = pathx)
  
  # Extract alignment positions under positive selection
  nucleotides <- lapply(align_set, FUN =  function(x) paste(x[pos_df$positions]))
  nucleotides <- data.frame(sequence = do.call(rbind, nucleotides),
                            gene = names(nucleotides))
  #
  return(nucleotides)
}

pos_df <- data.frame(positions = c(41,42,43))
path_aln <- "/Users/rprops/Desktop/test.aln"


