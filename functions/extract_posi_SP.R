#!/usr/bin/env Rscript

# Author: Ruben Props
# Date: 01-22-2018

userprefs <- commandArgs(trailingOnly = TRUE)
file_name <- userprefs[1]

library("phangorn")
library("seqLogo")

################################################################################
### Step 0: Read posigene file and extract PSG sites in PRANK alignment
################################################################################
results_posi <- read.table("./posigene_analysis/result_tables/Ramlibacter_MAG_results.tsv", 
                        header = TRUE, fill = TRUE, sep = "\t")
# Extract positions in alignment file for each gene
df_posi <- cbind(results_posi[c("Transcript")], 
                 apply(results_posi[, 27:ncol(results_posi)], 2, FUN = function(x) extract_nt.pos(x, pos = 8))
)


tmp2 <- matrix(do.call(rbind, strsplit(as.character(x), ";"))[, 8], 
              ncol = length(pos))
tmp2 <- as.numeric(matrix(apply(tmp2, 2, FUN = function(x) do.call(rbind, strsplit(as.character(x), " "))[,2]), 
                         ncol = length(pos))
)
tmp <- paste(tmp[,1], tmp[,2], sep = "-")

################################################################################
### Step 1: Import alignment file and construct parsimony tree
################################################################################

# Import bests Prank alignments
phydat_aln <- read.phyDat("/Users/rprops/desktop/prank.best.fas", format = "fasta")

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

# 

################################################################################
### Step 2: Extract site patterns from PSG sites
################################################################################

for(sequence in 1:length(anc.mpr)){
  print(sequence)
}

################################################################################
### Step 3: Extract site patterns from PSG sites in LCAs
################################################################################


################################################################################
### Step 4: Format and export results
################################################################################
