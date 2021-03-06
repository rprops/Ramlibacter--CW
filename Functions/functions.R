# Function for creating a common legend for 2 ggplot2 figures.
grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {
  
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  grid.newpage()
  grid.draw(combined)
  
}


# Function to extract %GC and geneID from gff file

extract_gc_from_gff <- function(inputgff, outputFolder){
  
  for(i in inputgff){
    d1 <- read.table(i, as.is = TRUE, sep = "\t")
    d2 <- t(sapply(strsplit(d1[[9]],";"), `[`, c(1:4)))
    d <- cbind(d1[,1], d2[,c(1,3)])
    colnames(d) <- c("contig", "contig_geneID", "GC")
    d <- data.frame(d)
    
    # Remove unwanted characters
    d$GC <- gsub(d$GC, pattern = "gc_cont=", replacement = "")
    d$contig_geneID <- gsub(d$contig_geneID, pattern = "ID=", replacement = "")
    d$contig_geneID <- gsub(d$contig_geneID, pattern = ".", replacement = "",
                            fixed = TRUE)
    d$GC <- as.numeric(d$GC); d <- d[!is.na(d$GC),]
    d <- data.frame(d, Genome = rep(gsub(".*/","",i), nrow(d)))
    
    if(i!=inputgff[1]) d.tot <- rbind(d.tot, d) else d.tot <- d
    # Export data
    write.table(d.tot, paste(outputFolder,"/","seqid_GC_",gsub(".*/","",i), ".tsv", sep = ""), row.names = FALSE,
                quote = FALSE)
  }
  
  # Plot
  # p <- easyGgplot2::ggplot2.histogram(data = d.tot, xName = 'GC',
  #                   groupName = 'Genome', alpha = 0.5,
  #                   legendPosition = "top", binwidth = 0.01)
  # print(p)
}

# Function to link gc content with gene functions

gc2function <- function(seq_id_gc, gene_id_seq_id, functions, gc_thresh = 0.75,
                        output = FALSE){
  # Import data
  seqid_GC <- read.table(seq_id_gc,
                         header = TRUE, stringsAsFactors = FALSE, sep = " ")
  gene_oid_2_seq_id <- read.table(gene_id_seq_id,
                                  stringsAsFactors = TRUE, col.names = c("gene_oid", "seq_id"),
                                  colClasses = "character", sep = " ")
  function_file <- read.table(functions,
                              stringsAsFactors = FALSE,  sep = "\t", header = TRUE,
                              quote=NULL, comment='')
  function_file$gene_oid <- as.character(function_file$gene_oid)
  
  # Merge files
  merged_df <- dplyr::inner_join(seqid_GC, gene_oid_2_seq_id, by = c("contig_geneID" = "seq_id"))
  merged_df <- dplyr::inner_join(merged_df, function_file, by = "gene_oid")
  tot_genes <- nrow(merged_df)
  
  # Filter out based on threshold
  merged_df <- merged_df[merged_df$GC>gc_thresh,]
  
  # Rank from large to small
  merged_df <- merged_df[order(merged_df$GC), ]
  
  # Print brief summary
  cat(date(), ' --- There are', nrow(merged_df), "genes with >", gc_thresh, "%\n" ,sep = " ")
  cat(date(), ' --- This is', round(100*nrow(merged_df)/tot_genes,2), "% of all genes\n",sep = " ")
  cat(date(), ' --- The 10 genes with the highest GC% are:\n', sep = " ")
  output <- tail(data.frame(function_id = merged_df[, ncol(merged_df)-2], function_name = merged_df[, ncol(merged_df)-1], 
                            GC = 100*merged_df$GC), n = 10)
  print(output)
  # Write output table
  if(output == TRUE){
    write.table(merged_df, "GC_function.tsv",
                row.names = FALSE, quote = FALSE)
  }
  return(merged_df)
}

# Modified igraph functions

plot_network_custom <- function (g, physeq = NULL, type = "samples", color = NULL, shape = NULL, 
                                 point_size = 4, alpha = 1, label = "value", hjust = 1.35, 
                                 line_weight = 0.5, line_color = color, line_alpha = 0.4, 
                                 layout.method = layout.fruchterman.reingold, title = NULL, label_size = 3) 
{
  if (vcount(g) < 2) {
    stop("The graph you provided, `g`, has too few vertices. \\n         Check your graph, or the output of `make_network` and try again.")
  }
  if (type %in% c("taxa", "species", "OTUs", "otus", "otu")) {
    type <- "taxa"
  }
  edgeDF <- data.frame(get.edgelist(g))
  edgeDF$id <- 1:length(edgeDF[, 1])
  vertDF <- layout.method(g)
  colnames(vertDF) <- c("x", "y")
  vertDF <- data.frame(value = get.vertex.attribute(g, "name"), 
                       vertDF)
  if (!is.null(physeq)) {
    extraData <- NULL
    if (type == "samples" & !is.null(sample_data(physeq, 
                                                 FALSE))) {
      extraData = data.frame(sample_data(physeq))[as.character(vertDF$value), 
                                                  , drop = FALSE]
    }
    else if (type == "taxa" & !is.null(tax_table(physeq, 
                                                 FALSE))) {
      extraData = data.frame(tax_table(physeq))[as.character(vertDF$value), 
                                                , drop = FALSE]
    }
    if (!is.null(extraData)) {
      vertDF <- data.frame(vertDF, extraData)
    }
  }
  graphDF <- merge(reshape2::melt(edgeDF, id = "id"), vertDF, 
                   by = "value")
  p <- ggplot(vertDF, aes(x, y))
  p <- p + theme_bw() + theme(panel.grid.major = element_blank(), 
                              panel.grid.minor = element_blank(), axis.text.x = element_blank(), 
                              axis.text.y = element_blank(), axis.title.x = element_blank(), 
                              axis.title.y = element_blank(), axis.ticks = element_blank(), 
                              panel.border = element_blank())
  p <- p + geom_point(aes_string(color = color, shape = shape), 
                      size = point_size, na.rm = TRUE)
  if (!is.null(label)) {
    p <- p + geom_text(aes_string(label = label), size = label_size, 
                       hjust = hjust, na.rm = TRUE)
  }
  p <- p + geom_line(aes_string(group = "id", color = line_color), 
                     graphDF, size = line_weight, alpha = line_alpha, na.rm = TRUE)
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  return(p)
}


# This function will read in codonO output for multiple genome CDS fasta files
# and merge them in to a single dataframe.

codonO_2_df <- function(pathx, patternx = "codonO", patterngene = "gene|Ga"){
  result_list <- list()
  files <- list.files(pathx, pattern = patternx)
  for(i in 1:length(files)){
    input_file <- read.table(paste(pathx, files[i], sep = ""), sep = ",", blank.lines.skip = TRUE, 
                             allowEscapes = FALSE, skipNul = TRUE)
    input_file$V1 <- gsub(input_file$V1, pattern = "\t", replacement = "")
    Genes <- do.call(rbind, strsplit(input_file$V1[grep(patterngene, input_file$V1)], " "))[,1]
    input_file$V1 <- gsub(input_file$V1, pattern = " ", replacement = "")
    input_file <- data.frame(GC = input_file$V1[grep(x = input_file$V1, pattern = "GC*.*=")], 
                             SCUO = rep(input_file$V1[grep(x = input_file$V1, pattern = "SCUO*")], each = 4),
                             Gene = rep(Genes, each = 4),
                             GCx = rep(c("GC_mean", "GC1", "GC2", "GC3"), length(Genes)),
                             Genome = rep(files[i], (length(Genes)*4))
    )
    input_file$Gene <- gsub(".*>", "",input_file$Gene)
    input_file$GC <- as.numeric(gsub(input_file$GC, pattern = ".*=", replacement = "")) 
    input_file$SCUO <- as.numeric(gsub(input_file$SCUO, pattern = ".*=", replacement = ""))
    result_list[[i]] <- input_file
  }
  
  # Reformat list into long format dataframe
  result_df <- do.call("rbind", result_list)
  
  # Return dataframe
  return(result_df)
}

# Get lower triangle of the correlation matrix
get_lower_tri <- function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}


# This function formats the KEGG orthology mapping file (ko_id to hierarchical 
# functional annotation)
format_ko <- function(path = "./Mapping_files/ko00000.keg"){
  # Import data
  ko_path_df <- read.table(path, header = FALSE, sep = ";",
                           skip = 3, quote = "", fill = TRUE, 
                           col.names = c("Level", "KO", "Function_abbrev", "Function_spec"))
  ko_path_df <- ko_path_df[1:(nrow(ko_path_df)-1), ] # remove tailing "!" at the end of file
  
  # Remove empty rows
  ko_path_df$KO <- as.character(ko_path_df$KO);  ko_path_df$Level <- as.character( ko_path_df$Level)
  ko_path_df$KO[grep("A",ko_path_df$Level)] <- ko_path_df$Level[grep("A",ko_path_df$Level)]
  ko_path_df$Level[grep("A",ko_path_df$Level)] <- "A"
  ko_path_df <- ko_path_df[!ko_path_df$KO == "", ]
  
  ko_path_df <- data.frame(ko_path_df, level_A = "A", level_B = "B", level_C = "C",
                           stringsAsFactors = FALSE)
  ko_path_df$KO <- as.character(ko_path_df$KO)
  
  # Get positions where to replicate higher hierarichcal level
  pos_A <- c(c(1:nrow(ko_path_df))[ko_path_df$Level %in% "A"], nrow(ko_path_df)+1)
  pos_B <- c(c(1:nrow(ko_path_df))[ko_path_df$Level %in% "B"], nrow(ko_path_df)+1)
  pos_C <- c(c(1:nrow(ko_path_df))[ko_path_df$Level %in% "C"], nrow(ko_path_df)+1)
  
  
  for(i in 1:(length(pos_A)-1)){
    ko_path_df$level_A[pos_A[i]:(pos_A[i+1]-1)] <- ko_path_df$KO[pos_A[i]]
  }
  
  for(i in 1:(length(pos_B)-1)){
    ko_path_df$level_B[pos_B[i]:(pos_B[i+1]-1)] <- ko_path_df$KO[pos_B[i]]
  }
  
  for(i in 1:(length(pos_C)-1)){
    ko_path_df$level_C[pos_C[i]:(pos_C[i+1]-1)] <- ko_path_df$KO[pos_C[i]]
  }
  
  # Remove all rows with level A, B, C - and level column
  ko_path_df <- ko_path_df[!ko_path_df$Level %in% c("A", "B", "C"), ]
  ko_path_df$level_A <- gsub(ko_path_df$level_A, pattern = "<b>|</b>", replacement = "")
  ko_path_df$level_B <- gsub(ko_path_df$level_B, pattern = "<b>|</b>", replacement = "")
  ko_path_df <- ko_path_df[, -1]
  
  # Remove redundant ID before first space
  ko_path_df$level_A <- gsub("^[^ ]* ", "", ko_path_df$level_A)
  ko_path_df$level_B <- gsub("^[^ ]* ", "", ko_path_df$level_B)
  ko_path_df$level_C <- gsub("^[^ ]* ", "", ko_path_df$level_C)
  # remove [PATH**] pattern
  ko_path_df$level_C <- gsub("\\[[^\\]]*\\]", "", ko_path_df$level_C, perl = TRUE)
  
  colnames(ko_path_df) <- c("ko_id", "ko_function_abbrev", "ko_function_spec",
                            "ko_level_A", "ko_level_B", "ko_level_C")
  return(ko_path_df)
  
}

# Small function to extract amino acid code and codon (anchor species) 
# from result file from PosiGene
# Usage on imported posigene dataframe:
# data_posi2 <- apply(data_posi2[, 27:105], 2, FUN = function(x) extract_aa(x))
extract_aa <- function(x, pos = c(10:11)){
  tmp <- matrix(do.call(rbind, strsplit(as.character(x), ";"))[, pos], 
                ncol = length(pos))
  tmp <- matrix(apply(tmp, 2, FUN = function(x) do.call(rbind, strsplit(as.character(x), " "))[,2]), 
                ncol = length(pos))
  tmp <- paste(tmp[,1], tmp[,2], sep = "-")
  return(tmp)
}


extract_nt.pos <- function(x, pos = 8){
  tmp <- matrix(do.call(rbind, strsplit(as.character(x), ";"))[, pos], 
                ncol = length(pos))
  tmp <- as.numeric(matrix(apply(tmp, 2, FUN = function(x) do.call(rbind, strsplit(as.character(x), " "))[,2]), 
                ncol = length(pos))
  )
  return(tmp)
}

# results_posi$Site.under.positve.Selection.1.1..Probability.to.be.under.positive.selection.2..Position.in.amino.acid.sequence.of.anchor.species.3..Position.in.nucleotide.sequence.of.anchor.species.4.Position.in.main.protein.alignment.fastp.clustalw.aln.5..Position.in.protein.clustal.subalignment...aln.6..Position.in.to.codon.backtranslated.clustal.subalignment..Codon_aln.fasta..7..Position.in.prank.alignment..protein...8..Position.in.prank.alignment..codon...prank.best.fas..9..Position.in.PAML.ready.prank.alignment..prank.best.fas.prepared.gb..10..Amino.acid.in.anchor.species.11..Codon.in.anchor.species.
# 
# 
# x <- results_posi[,27]
# extract_aa(x)
# 
# for(i in 27:215){
#   x <- results_posi[,i]
#   tmpx <- extract_nt.pos(x, pos = 11)
#   str(tmpx)
#   print(i)
# }

# Function to read in multiple posigene result files and concatenate them
# into a single dataframe
read_posi <- function(pathx, patt = "_results.tsv"){
  files_short <- list.files(pathx, pattern = patt)
  files_long <- list.files(pathx, full.names = TRUE)
  for(i in 1:length(files_short)){
    print(files_short[i])
    print(files_long[i])
    data_tmp <- read.table(files_long[i],
                            header = TRUE, fill = TRUE, sep = "\t")
    data_tmp$Transcript <- as.factor(data_tmp$Transcript)
    data_tmp <- data.frame(sample_file = files_short[i], data_tmp)
    if(i == 1){
      data_posi <- data_tmp
    } else {
      data_posi <- bind_rows(data_posi, data_tmp)
    }
  }
  data_posi$sample_file <- gsub(patt, "", data_posi$sample_file)
  return(data_posi)
}

# Function to format IMG compare-genome files to STAMPS compatible format

format_STAMPS <- function(pathx, patho){
  return_first <- function(x){
    y <- x[1]
    return(y)
  }
  
  df <- read.delim(pathx)
  dim(df)
  
  # temp files
  tmp1 <- do.call(rbind,by(df[, 3:ncol(df)], INDICES = factor(df$Func_name), 
                           FUN = colSums))
  
  tmp2 <- aggregate(Func_id~Func_name, data = df, FUN = return_first)
  
  # merge
  output <- cbind(tmp2, tmp1)
  output <- output[,c(2,1,3:ncol(output))]
  
  # write
  write.table(output, file = patho, sep = "\t", row.names = FALSE, quote=FALSE)
}

# format_STAMPS(pathx = "./IMG_annotation/STAMP_profiles/abundance_cog_118624.tab.xls",
              # patho = "./IMG_annotation/STAMP_profiles/STAMPS_abundance_cog_118624.tab.xls")
# 
# format_STAMPS(pathx = "./IMG_annotation/STAMP_profiles/abundance_ko_71619.tab.xls",
#               patho = "./IMG_annotation/STAMP_profiles/STAMPS_abundance_ko_71619.tab.xls")

# Computing the p-value of correlations
# http://www.sthda.com/english/wiki/visualize-correlation-matrix-using-correlogram
# To compute the matrix of p-value, a custom R function is used :
## mat : is a matrix of data
## ... : further arguments to pass to the native R cor.test function
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

# Querry a gene list for PSGs
# querry_PSG <- function(ref, querry, plot = FALSE){
#   
# }


# Querry dataframe for module completness
get_mq <- function(x, genome_label, ko_label){
  for(genome in unique(x[, genome_label])){
    print(genome)
    KO_metqy <- x[x[, genome_label] %in% genome, ]
    KO_metqy <- data.frame(genome_id = genome,
                                 KOs = paste(KO_metqy[, ko_label], collapse = ";"),
                                 stringsAsFactors = FALSE
    )
    KO_metqy$KOs <- gsub("KO:","",KO_metqy$KOs)
    
    # Run MetQy (takes a few seconds)
    query_output <- MetQy::query_genomes_to_modules(KO_metqy, splitBy='[;]',
                                                    GENOME_ID_COL = "genome_id", GENES_COL = "KOs", 
                                                    META_OUT = TRUE)
    
    modules_table = cbind(t(query_output$MATRIX), query_output$METADATA)
    modules_table <- modules_table[, -c(8, 9)]
    colnames(modules_table)[1] <- "module_completeness"
    modules_table <- data.frame(modules_table, Genome = genome)
    if(genome == unique(x[, genome_label])[1]){
      results <- modules_table
    } else{
      results <- rbind(results, modules_table)
    }
  }
  return(results)
}

