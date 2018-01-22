#!/usr/bin/env Rscript

# Author: Ruben Props
# Date: 01-22-2018

userprefs <- commandArgs(trailingOnly = TRUE)
file_name <- userprefs[1]

library("phangorn")
library("seqLogo")
library("dplyr")
library("tidyr")

################################################################################
### Required functions
################################################################################
# Small function to extract sites under positive selection in PRANK alignment
# from result file of PosiGene
extract_nt.pos <- function(x, pos = 8, gene_names){
  x_c <- paste(gene_names, x, sep = "; ")
  tmp <- matrix(do.call(rbind, strsplit(as.character(x_c), ";"))[, c(1, pos)], 
                ncol = length(c(1, pos)))
  tmp[,2] <-do.call(rbind, strsplit(tmp[,2], " "))[, 2]
  return(tmp[,2])
}

################################################################################
### Step 0: Read posigene file and extract PSG sites in best PRANK alignment
################################################################################
# Import and format posigene data
results_posi <- read.table("./posigene_analysis/result_tables/Ramlibacter_MAG_results.tsv", 
                        header = TRUE, fill = TRUE, sep = "\t")
gene_names <- as.character(results_posi$Transcript)
sites_results_posi <- results_posi %>% 
  dplyr::select(contains("Site.under.positve.Selection"))
colnames(sites_results_posi)[1] <- "Site.under.positve.Selection.1"

# Extract site positions under positive selection for each gene
df_posi <- apply(sites_results_posi, 2, 
                       FUN = function(x) 
                         extract_nt.pos(x, pos = 8, gene_names = gene_names))
rownames(df_posi) <- gene_names

# Convert to long format 
df_posi_long <- melt(df_posi)
colnames(df_posi_long) <- c("Gene", "PSG_number", "PSG_position")
df_posi_long$PSG_position <- as.numeric(as.character(df_posi_long$PSG_position))
df_posi_long$Gene <- factor(df_posi_long$Gene)
df_posi_long <- df_posi_long[!is.na(df_posi_long$PSG_position), ]

################################################################################
### Step 1: Import alignment file and construct parsimony tree
################################################################################

# Import bests Prank alignments
phydat_aln <- read.phyDat("/Users/rprops/desktop/prank.best.fas", format = "fasta")
df_posi_long_sb <- df_posi_long %>% filter(Gene == "2727790668")

# Import tree? // optional for now
# MyTree <- read.tree("/Users/rprops/desktop/2727790668_tree.newick")
# MyTree$tip.label <- c("gene144", "gene3558", "gene1596", "2727790668")

# Create parsiomny tree with Fitch's algorithm
tree <- pratchet(phydat_aln, trace = 0)

# Return maximum parsimony score
parsimony(tree, phydat_aln)

# Marginal reconstruction of the ancestral character states.
anc.mpr = ancestral.pars(tree, phydat_aln, "MPR")

# Plot character states of LCA/genes at designated site
# plotAnc(tree, anc.mpr, attr(anc.mpr, "index")[3])

# Rename LCA nodes according to the genes they encompass
# This necesary to later select the most appropriate LCA for the target species
# gene

index <- (1:length(desc_list))[unlist(lapply(desc_list, FUN = function(x) length(x)>1))]
desc_list <- phangorn::allDescendants(tree)
desc_df <- matrix(ncol = 2, nrow = sum(unlist(lapply(desc_list, FUN = function(x) length(x)>1))))
for(i in 1:nrow(desc_df)){
    desc_df[i, 2] <- print(paste("LCA", 
                                     paste(desc_list[[index[i]]], collapse="_", sep = "_"), sep = "_"))
    desc_df[i, 1] <- index[i]
}
desc_df <- data.frame(desc_df); colnames(desc_df) <- c("node_name", "LCA_name")

################################################################################
### Step 2: Extract site patterns from PSG sites
################################################################################
# Extract PSG nucleotide variations across reference and target gene
sp_index <- attr(anc.mpr, "index") # index indicating site patterns
sp_index_sb <- sp_index[df_posi_long_sb$PSG_position]

for(sequence in 1:length(anc.mpr)){
  psg_df <- anc.mpr[[sequence]][sp_index_sb, ]
  for(i_nt in 1:4){
    if(i_nt == 1) psg_df[psg_df[, i_nt] > 0, i_nt] <- "A"
    if(i_nt == 2) psg_df[psg_df[, i_nt] > 0, i_nt] <- "C"
    if(i_nt == 3) psg_df[psg_df[, i_nt] > 0, i_nt] <- "G"
    if(i_nt == 4) psg_df[psg_df[, i_nt] > 0, i_nt] <- "T"
  }
  psg_df[psg_df=="0"] <- " "
  psg_df <- data.frame(psg_df); colnames(psg_df) <- c("A","C","G","T")
  # psg_df <- psg_df %>% tidyr::unite(col = c("A","C","G","T"), sep = ";")
  # print(names(anc.mpr)[sequence])
  psg_df <- psg_df %>% tidyr::unite('all_cols', A:T, sep = ";")
  psg_df$all_cols <- gsub(" ;", "", psg_df$all_cols, fixed = TRUE)
  psg_df <- data.frame(Gene = names(anc.mpr)[sequence], site = df_posi_long_sb$PSG_position, 
                       nucleotide = psg_df$all_cols)
  if(sequence == 1) psg_final <- psg_df else psg_final <- rbind(psg_final, psg_df)
}

################################################################################
### Step 4: Format and export results
################################################################################
