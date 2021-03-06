---
title: "Metagenomic analysis of secondary cooling water microbial communities"
author: "Ruben Props"
date: "19 december, 2018"
output:
  html_document:
    code_folding: show
    highlight: haddock
    keep_md: yes
    theme: united
    toc: yes
    toc_float:
      collapsed: no
      smooth_scroll: yes
      toc_depth: 2
    css: report_styles.css
    df_print: paged
editor_options: 
  chunk_output_type: console
---




# A. 16S analysis
## Network analysis on relative abundances  

In the first analysis we utilize the raw counts. As such this network will be built using relative abundance data.


## Network analysis on absolute abundances  

The second network will be built using absolute abundance data by multiplying the relative taxon abundances by the total cell density. The final obtained counts will expressed as nr. of cells measured in 50 µL samples. We also only consider the OTUs that were left after the prevalence filtering conducted in the network construction with relative abundance data.  


```r
# Import cell count data
cell_counts <- read.csv("16S_data/cell_counts.csv")
cell_counts$sample_title <- as.factor(cell_counts$sample_title)

# Import metadata
meta_16S <- read.csv("16s_data/Metadata.csv")[1:77,]; rownames(meta_16S) <- meta_16S$sample_title

# Calculate proportions
phy_df_rel <- transform_sample_counts(phy_df, function(x) x/sum(x))

# Add metadata
sample_data(phy_df_rel) <- sample_data(meta_16S)

# Select samples for which corresponding counts are available
cell_counts <- cell_counts[cell_counts$sample_title %in% sample_names(phy_df_rel), ]
cell_counts <- droplevels(cell_counts)
phy_df_abs <- prune_samples(sample_names(phy_df_rel) %in% cell_counts$sample_title, phy_df_rel)

# Multiply with cell counts in 50 µL of sample
otu_table(phy_df_abs) <- otu_table(phy_df_abs) * cell_counts$Number_of_cells

# Select taxa that were selected based on prevalence in previous chunk
phy_df_abs <- prune_taxa(taxa_names(phy_df_filtered), phy_df_abs)

# Round absolute abundances to integers
otu_table(phy_df_abs) <- round(otu_table(phy_df_abs), 0)

# Construct network
sp_easi_abs <- spiec.easi(phy_df_abs, method='mb', lambda.min.ratio=1e-2,
                           nlambda=20, icov.select.params=list(rep.num=50))
```

```
## Normalizing/clr transformation of data with pseudocount ...
```

```
## Inverse Covariance Estimation with mb ...
```

```
## Model selection with stars ...
```

```
## Done!
```

```r
ig.mb_abs <- adj2igraph(sp_easi_abs$refit,  vertex.attr = list(name=taxa_names(phy_df_abs)))
vsize_abs <- Biobase::rowMedians(clr(otu_table(phy_df_abs), 1))+15
Lineage_abs <- tax_table(phy_df_abs)[,"Lineage"]
Lineage_abs <- factor(Lineage_abs, levels = unique(Lineage_abs))
vweights_abs <- summary(symBeta(getOptBeta(sp_easi_abs), mode='maxabs'))
MAGs <- c(); MAGs[taxa_names(phy_df_abs)=="Otu00001"]  <- "Ramlibacter sp. MAG"
MAGs[taxa_names(phy_df_abs)=="Otu00002"]  <- "Bacteroidetes sp. MAG1"
MAGs[taxa_names(phy_df_abs)=="Otu00003"]  <- "Bacteroidetes sp. MAG2"
MAGs[is.na(MAGs)] <-""
```



```r
# Plot network inferred from absolute abundances
# png(file = "./Figures/Figures_network/NETWORK-ABS-CX-C30-A25.png", width = 9, height = 9, res = 500, units = "in")
plot_network_custom(ig.mb_abs, phy_df_abs, type='taxa',
             line_weight = 2, hjust = 0.5,
             point_size = 0.1, alpha = 0.01, label_size = 3.95)+
  # scale_fill_brewer(palette = "Paired")+
  # scale_color_brewer(palette = "Paired")+
  # scale_fill_manual(values = c("#e2a2fd", brewer.pal(n = 12, "Paired")[c(3:8,1:2,11,12)]) )+
    scale_fill_manual(values = c(col_RAMLI, brewer.pal(n = 12, "Paired")[c(6,3,4,4,7,8,1:2,11,12)]) )+
  geom_point(aes(size = vsize_abs, fill = Lineage_abs), alpha = 0.5,
             colour="black", shape=21)+
  guides(size = FALSE,
    fill  = guide_legend(title = "Lineage", override.aes = list(size = 5),
                         nrow = 4),
    color = FALSE)+
  theme(legend.position="bottom", legend.text=element_text(size=12),
        text = element_text(size = 12),
        plot.margin = unit(c(1,1,1,1), "cm"))+
  scale_size(range = c(5, 15))+
  geom_label_repel(aes(label = MAGs), fontface = 'bold', color = 'black',
                   box.padding = 0.35, point.padding = 0.5,
                   segment.color = 'black',
                   size = 4,
                       # Width of the line segments.
                   segment.size = 1.5,
                   # Draw an arrow from the label to the data point.
                   arrow = arrow(length = unit(0.015, 'npc')),
                   nudge_x = -0.1,
                   nudge_y = 0.6
  )
```

<img src="Figures/cached/network-analysis-absolute-plot-1.png" style="display: block; margin: auto;" />

```r
# dev.off()
```

## Plots  


```r
# Plot absolute OTU dynamics of OTU1
df_abs <- psmelt(phy_df_abs)
col_RAMLI <- "#887CAF"

# Need to account for dilution factor of 2 and 50 µL volume measured
p_abs_otu1 <- df_abs %>% dplyr::filter(OTU == "Otu00001") %>% 
  ggplot(aes(x = Timepoint, y = 2*Abundance/50))+
  facet_grid(.~Reactor.cycle, scales = "free")+
  scale_shape_manual(values = c(21,24))+
  geom_line(size = 1.5, linetype = 2, color = adjustcolor("#000000", 0.5))+
  geom_point(size = 4, fill = col_RAMLI, aes(shape = Reactor_status, 
                                            alpha = Reactor_status),
             color = "black")+
  scale_alpha_manual(values = c(0.5,1))+
  theme_bw()+
  ylab(expression("Otu00001 abundance - cells µL"^"-1"))+
  xlab("Time relative to reactor start - days")+
  scale_y_continuous(breaks = seq(0,50e3,5e3)*2/50, limits = c(0,30e3)*2/50)+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=14),
        legend.background = element_rect(fill="transparent"),
        strip.text.x=element_text(size=18),
        legend.position = "top")+
  guides(shape = guide_legend(title="Reactor status", ncol =1),
         alpha = FALSE)
print(p_abs_otu1)
```

<img src="Figures/cached/OTU1-dynamics-1.png" style="display: block; margin: auto;" />


```r
# Boxplot of absolute abundances
col_RAMLI <- "#887CAF"

p_abs_box <- df_abs %>% 
  dplyr::filter(OTU %in% c("Otu00001","Otu00002","Otu00003")) %>% 
  ggplot(aes(x = OTU, y = 2*Abundance/50, fill = OTU))+
  geom_jitter(size = 2,
             color = "black", shape = 21, width = 0.1, alpha = 0.5)+
  geom_violin(alpha = 0.4, adjust = 1, draw_quantiles = TRUE,
              trim = TRUE)+
  stat_summary(fun.data=mean_sdl, fun.args = list(mult = 1), 
                 geom="pointrange", color="black")+
  # geom_boxplot(width = 0.5, alpha = 0.4, outlier.shape = NA)+
  scale_fill_manual(values = c(col_RAMLI, col_bac1, col_bac2))+
  theme_bw()+
  ylab(expression("Cells µL"^"-1"))+
  xlab("")+
  # ylab("")+
  scale_y_continuous(breaks = seq(0,50e3,10e3)*2/50, limits = c(-5e3,40e3)*2/50)+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        axis.text.x=element_text(size=14),
        title=element_text(size=16), legend.text=element_text(size=14),
        legend.background = element_rect(fill="transparent"),
        strip.text.x=element_text(size=18),
        legend.position = "sqrt")+
  guides(shape = guide_legend(title="Reactor status", ncol =1),
         fill = FALSE)+
    theme(axis.line = element_line(size = 1, colour = "grey80"),
        panel.border = element_blank())
  # ggtitle(expression("Cells mL"^"-1"))
  # ggtitle(expression("Absolute abundance (cells mL"^"-1)"))


print(p_abs_box)
```

<img src="Figures/cached/OTU1-box-abs-1.png" style="display: block; margin: auto;" />


```r
# psmelt relative abundance data
df_rel <- psmelt(phy_df_rel)
col_RAMLI <- "#887CAF"

# Boxplot of relative abundances
p_rel_box <- df_rel %>% 
  dplyr::filter(OTU %in% c("Otu00001","Otu00002","Otu00003")) %>% 
  ggplot(aes(x = OTU, y = 100*Abundance, fill = OTU))+
  geom_jitter(size = 2,
             color = "black", shape = 21, width = 0.1, alpha = 0.5)+
  geom_violin(alpha = 0.4, adjust = 1, draw_quantiles = TRUE)+
  stat_summary(fun.data=mean_sdl, fun.args = list(mult = 1), 
                 geom="pointrange", color="black")+
  # geom_boxplot(width = 0.5, alpha = 0.4, outlier.shape = NA)+
  scale_fill_manual(values = c(col_RAMLI, col_bac1, col_bac2))+
  theme_bw()+
  # ylab(expression("Relative abundance - %"))+
  xlab("")+
  ylab("")+
  scale_y_continuous(breaks = 100*seq(0,1,0.2), limits = 100*c(-0.05,1))+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        axis.text.x=element_text(size=14, angle = 0),
        title=element_text(size=16), legend.text=element_text(size=14),
        legend.background = element_rect(fill="transparent"),
        strip.text.x=element_text(size=18),
        legend.position = "top")+
  guides(shape = guide_legend(title="Reactor status", ncol =1),
         fill = FALSE)+
    theme(axis.line = element_line(size = 1, colour = "grey80"),
        panel.border = element_blank())+
  ggtitle("Relative abundance (%)")


print(p_rel_box)
```

<img src="Figures/cached/OTU1-box-rel-1.png" style="display: block; margin: auto;" />

## Physicochemistry


```r
  # Assess correlations
  Physico_df_trim_abs_OTU1 <- df_abs[, c("OTU", "Timepoint","Dilution_factor", 
                              "Formate", "Acetate", "Sulfate", 
                              "Nitrate", "Chloride", "Temperature", 
                              "Conductivity", "pH", "Sample", "Abundance")] %>%
  dplyr::filter(!is.na(Formate), OTU == "Otu00001")
  
  colnames(Physico_df_trim_abs_OTU1)[13] <- "OTU1 - abundance"
  
  Physico_df_trim_abs_OTU2 <- df_abs[,  c("OTU", "Sample", "Abundance")] %>%
  dplyr::filter(OTU == "Otu00002")
  colnames(Physico_df_trim_abs_OTU2)[3] <- "OTU2 - abundance"
  
  Physico_df_trim_abs_OTU3 <- df_abs[, c("OTU", "Sample", "Abundance")]%>%
  dplyr::filter(OTU == "Otu00003")
  colnames(Physico_df_trim_abs_OTU3)[3] <- "OTU3 - abundance"
  
  Physico_df_trim_abs_merged <- left_join(Physico_df_trim_abs_OTU1, 
                                            Physico_df_trim_abs_OTU2,
                                            by = "Sample")
  
  Physico_df_trim_abs_merged <- left_join(Physico_df_trim_abs_merged, 
                                            Physico_df_trim_abs_OTU3,
                                            by = "Sample")
  # Same for relative abundances
  Physico_df_trim_rel_OTU1 <- df_rel[, c("OTU", "Timepoint","Dilution_factor", 
                              "Formate", "Acetate", "Sulfate", 
                              "Nitrate", "Chloride", "Temperature", 
                              "Conductivity", "pH", "Sample", "Abundance")] %>%
  dplyr::filter(!is.na(Formate), OTU == "Otu00001")
  
  colnames(Physico_df_trim_rel_OTU1)[13] <- "OTU1 - abundance"
  
  Physico_df_trim_rel_OTU2 <- df_rel[,  c("OTU", "Sample", "Abundance")] %>%
  dplyr::filter(OTU == "Otu00002")
  colnames(Physico_df_trim_rel_OTU2)[3] <- "OTU2 - abundance"
  
  Physico_df_trim_rel_OTU3 <- df_rel[, c("OTU", "Sample", "Abundance")]%>%
  dplyr::filter(OTU == "Otu00003")
  colnames(Physico_df_trim_rel_OTU3)[3] <- "OTU3 - abundance"
  
  Physico_df_trim_rel_merged <- left_join(Physico_df_trim_rel_OTU1, 
                                            Physico_df_trim_rel_OTU2,
                                            by = "Sample")
  
  Physico_df_trim_rel_merged <- left_join(Physico_df_trim_rel_merged, 
                                            Physico_df_trim_rel_OTU3,
                                            by = "Sample")
  
  
  # Physico_df_trim_rel <- df_rel[, c("OTU", "Timepoint","Dilution_factor", 
  #                             "Formate", "Acetate", "Sulfate", 
  #                             "Nitrate", "Chloride", "Temperature", 
  #                             "Conductivity", "pH", "Abundance")] %>%
  # dplyr::filter(!is.na(Formate), OTU == otu)

  # Get P-values of correlations between OTUs and physicochemistry
  p.mat_abs <- cor.mtest(Physico_df_trim_abs_merged[, -c(1:3, 12,14,16)],
                                            method = "kendall")
  
  p.mat_rel <- cor.mtest(Physico_df_trim_rel_merged[, -c(1:3, 12,14,16)],
                                            method = "kendall")

  # Making correlation plots
  corrplot::corrplot(cor(Physico_df_trim_abs_merged[, -c(1:3, 12,14,16)],
                         method = "kendall"),
                     type="upper",
                     tl.col="black", tl.srt=45,
                     p.mat = p.mat_abs, sig.level = 0.05,
                     diag=FALSE, 
                     title = paste("Absolute abundances"))
```

<img src="Figures/cached/physico-data-1.png" style="display: block; margin: auto;" />

```r
  corrplot::corrplot(cor(Physico_df_trim_rel_merged[, -c(1:3, 12,14,16)],
                         method = "pearson"),
                     type="upper",
                     tl.col="black", tl.srt=45,
                     p.mat = p.mat_rel, sig.level = 0.05,
                     diag=FALSE, 
                     title = paste("Relative abundances"))
```

<img src="Figures/cached/physico-data-2.png" style="display: block; margin: auto;" />

```r
# Calculate element-wise N/P ratio assuming 1 µg/L of PO4
molP <- (1/94.9714)/1000 # in mmol/L of upper detection limit
min((Physico_df_trim_rel_OTU1$Nitrate/62.0049)/molP) # in mmol/L
```

```
## [1] 96.49557
```

```r
sd((Physico_df_trim_rel_OTU1$Nitrate/62.0049)/molP)
```

```
## [1] 391.1701
```

```r
# mean((14*Physico_df_trim_rel_OTU1$Nitrate/62.0049))
# sd((14*Physico_df_trim_rel_OTU1$Nitrate/62.0049))
```

# B. MetaG analysis

```r
# Read data
mean_coverage <- read.table("./SAMPLES-SUMMARY/bins_across_samples/mean_coverage.txt", header = TRUE)
std_coverage <- read.table("./SAMPLES-SUMMARY/bins_across_samples/std_coverage.txt", header = TRUE)
bin_size <- read.table("./SAMPLES-SUMMARY/general_bins_summary.txt", header = TRUE)[, c(1, 3, 6, 9)]
total_reads <- read.table("./Mapping_files/sample_reads.tsv", header = TRUE)
read_length <- 300

# From wide to long format
mean_coverage_long <- gather(mean_coverage, Sample_ID, coverage, 
                             SAMPLE_16:SAMPLE_65, factor_key=TRUE)

std_coverage_long <- gather(std_coverage, Sample_ID, std_coverage, 
                            SAMPLE_16:SAMPLE_65, 
                            factor_key=TRUE)

coverage_data <- data.frame(mean_coverage_long, 
                            std_coverage = std_coverage_long[,3])

# Read and add metadata
# meta <- read.csv2("metadata.csv")
# meta$Sample_ID <- gsub(meta$Sample_ID, pattern = ".", replacement = "_", fixed = TRUE)
data_total <- left_join(coverage_data, total_reads, by = "Sample_ID")
data_total <- left_join(data_total, bin_size, by = "bins")
# data_total <- left_join(data_total, meta, by =  "Sample_ID")
data_total$bins <- plyr::revalue(data_total$bins, c("BetIa_bin"="Ramlibacter sp. MAG",
                                                    "bacIa_vizbin1"="Bacteroidetes sp. MAG1",
                                                    "bacIa_vizbin2"="Bacteroidetes sp. MAG2"))
# Calculate relative abundance of the bins
data_total$mean_rel_abundance <- 100*(data_total$coverage*data_total$bin_size)/(read_length*data_total$Total_reads)
data_total$upper_rel_abundance <- 100*((data_total$coverage+data_total$std_coverage)*data_total$bin_size)/(read_length*data_total$Total_reads)
data_total$lower_rel_abundance <- 100*((data_total$coverage-data_total$std_coverage)*data_total$bin_size)/(read_length*data_total$Total_reads)

data_total$mean_rel_abundance_map <- 100*(data_total$coverage*data_total$bin_size)/(read_length*data_total$Mapped_reads)
data_total$upper_rel_abundance_map <- 100*((data_total$coverage+data_total$std_coverage)*data_total$bin_size)/(read_length*data_total$Mapped_reads)
data_total$lower_rel_abundance_map <- 100*((data_total$coverage-data_total$std_coverage)*data_total$bin_size)/(read_length*data_total$Mapped_reads)

# Add additional column that assigns PNCs to correct MAG
data_total$Genome_id <- factor(rep(c("Ramlibacter sp. MAG", "Ramlibacter sp. MAG", "Ramlibacter sp. MAG", "Bacteroidetes sp. MAG2", "Bacteroidetes sp. MAG1", "Bacteroidetes sp. MAG2"), 4))

# Plot genome size for all three genomes
data_total[data_total$bins %in% c("Bacteroidetes sp. MAG1",
                                "Bacteroidetes sp. MAG2","Ramlibacter sp. MAG"), ] %>% 
  ggplot(aes(x = bins, y = bin_size/(4*1e6), fill = Genome_id))+
  theme_bw()+
  geom_bar(width = 0.5, stat="identity", alpha = 0.7)+
  coord_flip()+
  scale_fill_manual(values = c(col_bac1, col_bac2, col_RAMLI))+
  theme(axis.text.x = element_text(angle = 0, size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 18),
        legend.title = element_text("Genome bin"), legend.position = "top")+
  xlab("")+
  ylab("Genome size (Mbp)")+
  guides(fill = FALSE)+
  geom_hline(yintercept = 1, size = 2, linetype="dotted")
```

<img src="Figures/cached/read-format data-1.png" style="display: block; margin: auto;" />

```r
# Plot irep for all 3 genomes
data_total[data_total$bins %in% c("Bacteroidetes sp. MAG1",
                                "Bacteroidetes sp. MAG2","Ramlibacter sp. MAG"), ] %>% 
  ggplot(aes(x = bins, y = mean_irep/4, fill = Genome_id))+
  theme_bw()+
  geom_bar(width = 0.5, stat="identity", alpha = 0.7)+
  coord_flip()+
  scale_fill_manual(values = c(col_bac1, col_bac2, col_RAMLI))+
  theme(axis.text.x = element_text(angle = 0, size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 18),
        legend.title = element_text(""), legend.position = "top")+
  xlab("")+
  ylab("Index of Replication (iRep)")+
  guides(fill = FALSE)+
  geom_hline(yintercept = 1.34, size = 2, linetype="dotted")
```

<img src="Figures/cached/read-format data-2.png" style="display: block; margin: auto;" />

# 1. Phylogenetic tree
## Ramlibacter sp.
<!-- ![RAxML tree for Ramlibacter sp. MAG](./Tree/ANI_tree_concat-lowres2.png) -->

# *2. Investigate MAG- and 16S-based abundances*
It is clear that there is significant %GC coverage bias present. The estimated relative abundances
from metagenomics do not quantitatively match with the V3-V4 16S rRNA gene amplicon data. This is probably due to the significant %GC bias that is associate with the MAG-based assessment.  

$$Relative\ abundance =100*(\frac{mean\ coverage * bin\ size}{read\ length*total\ sample\ reads })$$
Another option is to calculate relative to mapped number of reads:
$$Relative\ abundance =100*(\frac{mean\ coverage * bin\ size}{read\ length*total\ sample\ reads * \%mapped\ reads})$$


Import reference relative abundances from 16S data set in order to directly compare with metagenomic data set.

```r
df_16S <- read.delim("./Mapping_files/relative_abundance_16S.tsv",
                     header = TRUE, sep = "\t")
df_16S_long <- gather(df_16S, Sample_ID, relative_abundance_16S, 
                             SAMPLE_16:SAMPLE_65, factor_key=TRUE)
```



```r
# Subset for only the three complete genomes (not PNCs).
data_total_sb <- data_total[data_total$bins %in% c("Ramlibacter sp. MAG", "Bacteroidetes sp. MAG1", "Bacteroidetes sp. MAG2"),]

p_meta <- ggplot(data = data_total_sb, aes(x = bins, y = mean_rel_abundance, fill = bins))+
  geom_point(size = 4, shape = 21, alpha = 0.7)+
  scale_fill_manual(values = c(col_bac1, col_bac2, col_RAMLI))+
  theme_bw()+
  geom_errorbar(aes(ymin=lower_rel_abundance, 
                    ymax=upper_rel_abundance, 
                    width=.1))+
  facet_grid(.~Sample_ID)+
  # ylim(0,1)+ 
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=14),
        legend.background = element_rect(fill="transparent"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.x=element_text(size=18))+
  ylab("Mean relative abundance (%)")+
  ylim(-1,100)+
  ggtitle("Metagenomic - total reads")

# Corrected for mapped N° of reads
p_meta_mapped <- ggplot(data = data_total_sb, aes(x = bins, y = mean_rel_abundance_map, fill = bins))+
  geom_point(size = 4, shape = 21, alpha = 0.7)+
  scale_fill_manual(values = c(col_bac1, col_bac2, col_RAMLI))+
  theme_bw()+
  geom_errorbar(aes(ymin=lower_rel_abundance_map, 
                    ymax=upper_rel_abundance_map, 
                    width=.1))+
  facet_grid(.~Sample_ID)+
  # ylim(0,1)+ 
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=14),
        legend.background = element_rect(fill="transparent"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.x=element_text(size=18))+
  ylab("Mean relative abundance (%)")+
  ylim(-1,100)+
  ggtitle("Metagenomic - mapped reads")

p_16S <- ggplot(data = df_16S_long, aes(x = bins, y = relative_abundance_16S, fill = bins))+
  geom_point(size = 4, shape = 21, alpha = 0.7)+
  scale_fill_manual(values = c(col_bac1, col_bac2, col_RAMLI))+
  theme_bw()+
  facet_grid(.~Sample_ID)+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=14),
        legend.background = element_rect(fill="transparent"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.x=element_text(size=18))+
  ylab("Mean relative abundance (%)")+
  ylim(-1,100)+
  ggtitle("V3-V4 16S")

grid_arrange_shared_legend(p_meta, p_meta_mapped, p_16S, ncol = 1, nrow = 3)
```

<img src="Figures/cached/plot-data-1.png" style="display: block; margin: auto;" />

# *3. Investigate sequence characteristics within coding DNA sequences (CDS)*




# *4. Analysis of gene length distribution*
Here we use the dataframe made in the previous section to see if there is a significant difference in the gene length of the COGs within these three consensus genomes.  

Observation: They have very small genes: on average < 500bp.



We can do the same for the Pfams.


# 5. Identify unique functional genes (COG/Pfams)





# 6. COG functional categories
Get COG ID to COG functional category mapping file here: ftp://ftp.ncbi.nih.gov/pub/wolf/COGs/COG0303/cogs.csv    

The exact statistical analysis to compare genomes based on these profiles should be performed in STAMP.



# 7. KO pathways

* Get reference file that maps KO ids to pathways here: http://www.genome.jp/kegg-bin/get_htext?ko0000.keg (download htext).



# 8.  Synonymous Codon Usage Bias analysis using CodonO

* Get CodonO for linux here: http://sysbio.cvm.msstate.edu/software/CodonO/download

* Run `CU.linux input.genes.fna`  
* Output from CodonO:  
  inputfile.ok  ---- the SCUO units for each input sequence based on its order  
	inputfile.fik ---- the composition ratio of the i-th amino acid in the k-th sequence  
	inputfile.hijk ---- the frequency of the j-th degenerate codon for amino acid i in each sequence  

* This output was not sufficient to perform any analysis. Therefore I used the web browser to calculate the synonymous codon bias in the genes of the three genomes.

** Codon bias is proportional to mRNA production (Wan et al 2004)  

** However, strong correlation between SCUO and %GC has been reported, thereby
confounding possible "biological" effects.. Beware..



```r
# Import and format data for Ramlibacter sp. MAG
SCUO_RAMLI <- read.table("./IMG_annotation/IMG_2724679690_Ramlibacter_bin/Annotation/2724679690.genes.fna.codonO.output",
                         sep = ",", blank.lines.skip = TRUE, allowEscapes = FALSE, skipNul = TRUE
                         )
SCUO_RAMLI$V1 <- gsub(SCUO_RAMLI$V1, pattern = "\t", replacement = "")
Gene_RAMLI <- do.call(rbind, strsplit(SCUO_RAMLI$V1[grep("Ga0*", SCUO_RAMLI$V1)], " "))[,2]
SCUO_RAMLI$V1 <- gsub(SCUO_RAMLI$V1, pattern = " ", replacement = "")
SCUO_RAMLI <- data.frame(GC = SCUO_RAMLI$V1[grep(x = SCUO_RAMLI$V1, pattern = "GC*.*=")], 
                         SCUO = rep(SCUO_RAMLI$V1[grep(x = SCUO_RAMLI$V1, pattern = "SCUO*")], each = 4),
                         Gene = rep(Gene_RAMLI, each = 4)
)
SCUO_RAMLI$GC <- as.numeric(gsub(SCUO_RAMLI$GC, pattern = ".*=", replacement = ""))
SCUO_RAMLI$SCUO <- as.numeric(gsub(SCUO_RAMLI$SCUO, pattern = ".*=", replacement = ""))

# Import and format data for Bacteroidetes MAG1
SCUO_BAC1 <- read.table("./IMG_annotation/IMG_2724679691_Bacteroidetes_bin1/Annotation/2724679691.genes.fna.codonO.output",
                         sep = ",", blank.lines.skip = TRUE, allowEscapes = FALSE, skipNul = TRUE)
SCUO_BAC1$V1 <- gsub(SCUO_BAC1$V1, pattern = "\t", replacement = "")
Gene_BAC1 <- do.call(rbind, strsplit(SCUO_BAC1$V1[grep("Ga0*", SCUO_BAC1$V1)], " "))[,2]
SCUO_BAC1$V1 <- gsub(SCUO_BAC1$V1, pattern = " ", replacement = "")
SCUO_BAC1 <- data.frame(GC = SCUO_BAC1$V1[grep(x = SCUO_BAC1$V1, pattern = "GC*.*=")], 
                         SCUO = rep(SCUO_BAC1$V1[grep(x = SCUO_BAC1$V1, pattern = "SCUO*")], each = 4),
                         Gene = rep(Gene_BAC1, each = 4)
)
SCUO_BAC1$GC <- as.numeric(gsub(SCUO_BAC1$GC, pattern = ".*=", replacement = ""))
SCUO_BAC1$SCUO <- as.numeric(gsub(SCUO_BAC1$SCUO, pattern = ".*=", replacement = ""))


# Import and format data for Bacteroidetes MAG2
SCUO_BAC2 <- read.table("./IMG_annotation/IMG_2724679698_Bacteroidetes_bin2/Annotation/2724679698.genes.fna.codonO.output",
                         sep = ",", blank.lines.skip = TRUE, allowEscapes = FALSE, skipNul = TRUE)
SCUO_BAC2$V1 <- gsub(SCUO_BAC2$V1, pattern = "\t", replacement = "")
Gene_BAC2 <- do.call(rbind, strsplit(SCUO_BAC2$V1[grep("Ga0*", SCUO_BAC2$V1)], " "))[,2]
SCUO_BAC2$V1 <- gsub(SCUO_BAC2$V1, pattern = " ", replacement = "")
SCUO_BAC2 <- data.frame(GC = SCUO_BAC2$V1[grep(x = SCUO_BAC2$V1, pattern = "GC*.*=")], 
                         SCUO = rep(SCUO_BAC2$V1[grep(x = SCUO_BAC2$V1, pattern = "SCUO*")], each = 4),
                         Gene = rep(Gene_BAC2, each = 4)
)
SCUO_BAC2$GC <- as.numeric(gsub(SCUO_BAC2$GC, pattern = ".*=", replacement = ""))
SCUO_BAC2$SCUO <- as.numeric(gsub(SCUO_BAC2$SCUO, pattern = ".*=", replacement = ""))

# Merge data to one dataframe
SCUO_merged_gen <- data.frame(rbind(SCUO_RAMLI, SCUO_BAC1, SCUO_BAC2),
           Genome_ID = c(rep("Ramlibacter sp. MAG", nrow(SCUO_RAMLI)), rep("Bacteroidetes MAG1", nrow(SCUO_BAC1)), 
           rep("Bacteroidetes MAG2", nrow(SCUO_BAC2))),
            GCx = rep(c("GC_mean", "GC1", "GC2", "GC3"), 
                     (nrow(SCUO_RAMLI) + nrow(SCUO_BAC1) + nrow(SCUO_BAC2))/4
            )
)

# Merge codon bias data with KO pathway annotation
SCUO_merged <- dplyr::left_join(SCUO_merged_gen, merged_gc_ko[, c(1:2,4 ,14:21)], by = c("Gene" = "contig_geneID"))
```

```
## Error in `[.data.frame`(merged_gc_ko, , c(1:2, 4, 14:21)): undefined columns selected
```

```r
# Visualize differences in codon bias
p_SCUO.1 <- ggplot(data = SCUO_merged, aes (x = 100*GC, y = SCUO, fill = Genome_ID))+
  geom_point(size = 4, shape = 21, alpha = 0.7)+
  scale_fill_manual(values = c(col_bac1, col_bac2, col_RAMLI))+
  theme_bw()+
  facet_wrap(~GCx, ncol = 2)+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=14),
        legend.background = element_rect(fill="transparent"),
        axis.text.x = element_text(angle = 0, hjust = 1),
        strip.text.x=element_text(size=18),
        legend.position = "bottom")+
  ylab("SCUO")+
  xlab("%GC")+
  ylim(0,1)
```

```
## Error in ggplot(data = SCUO_merged, aes(x = 100 * GC, y = SCUO, fill = Genome_ID)): object 'SCUO_merged' not found
```

```r
print(p_SCUO.1)
```

```
## Error in print(p_SCUO.1): object 'p_SCUO.1' not found
```

```r
# Visualize differences in codon bias per codon position
p_SCUO.2 <- SCUO_merged %>% filter(GCx != "GC_mean") %>% 
  ggplot(aes (x = GCx, y = 100*GC, fill = Genome_ID))+
  geom_jitter(size = 4, shape = 21, alpha = 0.1, width = 0.2)+
  geom_boxplot(alpha = 0.2, size = 1.2, color = "darkorange")+
  scale_fill_manual("", values = c(col_bac1, col_bac2, col_RAMLI))+
  theme_bw()+
  facet_wrap(~Genome_ID, ncol = 2)+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=14),
        legend.background = element_rect(fill="transparent"),
        strip.text.x=element_text(size=18),
        legend.position = "bottom")+
  ylab("%GC")+
  xlab("Codon position")+
  guides(fill = FALSE)+
  ylim(0,100)
```

```
## Error in eval(lhs, parent, parent): object 'SCUO_merged' not found
```

```r
print(p_SCUO.2)
```

```
## Error in print(p_SCUO.2): object 'p_SCUO.2' not found
```

```r
# Subset to genes for which ko annotation is available
SCUO_merged_sb <- SCUO_merged[!is.na(SCUO_merged$ko_level_A), ]
```

```
## Error in eval(expr, envir, enclos): object 'SCUO_merged' not found
```

```r
SCUO_merged_sb <- SCUO_merged_sb[SCUO_merged_sb$GCx == "GC_mean", ]
```

```
## Error in eval(expr, envir, enclos): object 'SCUO_merged_sb' not found
```

```r
# Look at pathways enriched in high %GC
# SCUO_merged_sb[]


SCUO_merged_gen_gcmean <- SCUO_merged_gen %>% dplyr::filter(GCx == "GC_mean")

p_SCUO.3 <- ggplot(data = SCUO_merged_gen_gcmean, aes (x = Genome_ID, y = SCUO))+
  geom_jitter(size = 3, alpha = 0.3, shape = 21, aes(fill = Genome_ID))+
  geom_boxplot(alpha=0, size =1.5, color = "darkorange")+
  # scale_fill_brewer(palette = "Accent")+
  scale_fill_manual(values = c(col_bac1, col_bac2, col_RAMLI))+
  theme_bw()+
  # facet_wrap(Genome_ID~GCx)+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=14),
        legend.background = element_rect(fill="transparent"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.x=element_text(size=18))+
  ylab("Codon bias - SCUO")+
  xlab("")+
  ylim(0,1)+ 
  guides(fill=FALSE)+ 
  scale_x_discrete(labels=c("Bacteroidetes MAG1" = paste("Bacteroidetes MAG1 (n=",table(SCUO_merged_gen_gcmean$Genome_ID)[1],")", sep = ""),
                            "Bacteroidetes MAG2" = paste("Bacteroidetes MAG2 (n=",table(SCUO_merged_gen_gcmean$Genome_ID)[2],")", sep = ""),
                            "Ramlibacter sp. MAG" = paste("Ramlibacter sp. MAG (n=",table(SCUO_merged_gen_gcmean$Genome_ID)[3],")", sep = ""))
  )
# 
print(p_SCUO.3)
```

<img src="Figures/cached/Codon bias-1.png" style="display: block; margin: auto;" />

```r
tmp <- SCUO_merged_sb$genome_id
```

```
## Error in eval(expr, envir, enclos): object 'SCUO_merged_sb' not found
```

```r
tmp2 <- cbind(SCUO_merged_sb$ko_id, 
      c(rep(col_RAMLI, table(tmp)[3]), rep(col_bac1, table(tmp)[1]), rep(col_bac2, table(tmp)[2]))
)
```

```
## Error in cbind(SCUO_merged_sb$ko_id, c(rep(col_RAMLI, table(tmp)[3]), : object 'SCUO_merged_sb' not found
```

```r
write.table(tmp2, file = "All_KO.tsv", quote = FALSE,
            col.names = FALSE, row.names = FALSE)
```

```
## Error in is.data.frame(x): object 'tmp2' not found
```

```r
# merge with annotation
tmp_SCUO <- dplyr::left_join(SCUO_merged_gen, 
                             merged_gc_ko, by = c("Gene" = "contig_geneID"))

tmp_result_scuo <- tmp_SCUO %>% dplyr::filter(GCx == "GC_mean") %>% 
  dplyr::filter(grepl("ribosomal", ko_name) & genome_id == "Ramlibacter sp. MAG") %>% 
  dplyr::select(SCUO, Gene, ko_name, genome_id) %>% 
  distinct %>% 
  ggplot(., aes(x = genome_id, y = SCUO))+
  geom_boxplot()


tmp_SCUO %>% dplyr::filter(GCx == "GC_mean") %>% 
  dplyr::filter(grepl("ribosomal", ko_name) & genome_id == "Ramlibacter sp. MAG") %>% 
  dplyr::select(SCUO, Gene, ko_name, genome_id) %>% 
  distinct %>% 
  summarise(mean(SCUO), sd(SCUO))
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["mean(SCUO)"],"name":[1],"type":["dbl"],"align":["right"]},{"label":["sd(SCUO)"],"name":[2],"type":["dbl"],"align":["right"]}],"data":[{"1":"0.6010922","2":"0.08877455"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>


```r
# Import codonO results of reference genomes
ref_SCUO <- codonO_2_df(pathx = "./IMG_annotation/References/codonO_output/")

# reimport codonO results of RAMLI genome
Ramli_SCUO <- codonO_2_df(pathx = "./IMG_annotation/IMG_2724679690_Ramlibacter_bin/Annotation/",
                          patternx = "codonO")

# Merge dataframes
ref_RAMLI_SCUO <- rbind(ref_SCUO, Ramli_SCUO)

# Replace genome names by better annotated names
map_scuo <- read.delim("./Mapping_files/codonO_ref_names.tsv", stringsAsFactors = FALSE)
ref_RAMLI_SCUO$Genome <- as.character(ref_RAMLI_SCUO$Genome)
for(i in 1:nrow(map_scuo)){
  ref_RAMLI_SCUO$Genome[ref_RAMLI_SCUO$Genome %in% map_scuo$codon_file[i]] <- map_scuo$ref_name[i]
}
ref_RAMLI_SCUO$Genome[ref_RAMLI_SCUO$Genome %in% "2724679690.genes.fna.codonO.output"] <- "Ramli. sp. MAG"

# order data according to mean SCUO
# sum_scuo <- ref_RAMLI_SCUO %>% group_by(Genome) %>% summarize(mean_scuo = mean(SCUO))
# ref_RAMLI_SCUO$Genome <- factor(ref_RAMLI_SCUO$Genome, levels = sum_scuo$Genome[order(sum_scuo$mean_scuo)])

ord_list_bin <- c("Lim. sp. Rim11", "Lim. sp. 103DPR2",
                  "Lim. sp. 2KL-27", "Lim. sp. Rim47",
                  "Lim. sp. II-D5", "Lim. sp. 2KL-3",
                  "Rhodo. sp. ED16","Rhodo. sp. T118",
                  "Curvi. sp. ATCC", "Curvi. sp. PAE-UM",
                  "Vario. sp. 110B", "Vario. sp. EPS",
                  "Ramli. sp. Leaf400", "Ramli. sp. TTB310",
                  "Ramli. sp. MAG", 
                  "Ramli. sp. 5-10"
                  )

# order rows and columns
ref_RAMLI_SCUO$Genome <- factor(ref_RAMLI_SCUO$Genome, levels = ord_list_bin)


# ref_RAMLI_SCUO$new <- ref_RAMLI_SCUO$Genome=="Lim. MAG (2724679690)"

# Plot SCUO profiles
p_ramli_SCUO <- ref_RAMLI_SCUO %>% dplyr::filter(GCx == "GC_mean") %>% 
  ggplot(aes(x = Genome, y = SCUO, fill = Genome))+
  theme_bw()+
    # geom_rect(data = tp, aes(fill = new), xmin = -Inf, xmax = Inf,
            # ymin = -Inf,ymax = Inf, alpha = 0.005, show.legend =FALSE, inherit.aes = FALSE)+
  geom_violin(alpha = 0.4, adjust = 1, draw_quantiles = TRUE)+
  scale_fill_manual(values = c(rep(adjustcolor("#c8c8ff",0.8),6), rep("#f8cf94",2), 
                               rep("#adf7ad",2), rep(adjustcolor("#000000",0.21),2),
                               rep("#e2a2fd",4)))+
  theme(axis.text=element_text(size=16), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=14),
        legend.background = element_rect(fill="transparent"),
        # axis.text.x = element_text(angle = 65, hjust = 1),
        strip.text.x=element_text(size=18),
        legend.position="bottom",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  guides(fill=FALSE)+
  stat_summary(fun.data=mean_sdl, fun.args = list(mult = 1), 
                 geom="pointrange", color="black")+
  xlab("")+
  ylab("")+
  # geom_vline(xintercept = 9.5, col = 'black', lwd = 2, linetype = 1, alpha = 0.6)+
  ylim(0,1)

p_ramli_SCUO
```

<img src="Figures/cached/compare-ramli-CB-1.png" style="display: block; margin: auto;" />

```r
p_ramli_GC <- ref_RAMLI_SCUO  %>% dplyr::filter(GCx == "GC_mean") %>% 
  ggplot(aes(x = Genome, y = GC, fill = Genome))+
  theme_bw()+
  # geom_hline(yintercept = 0.7, col = 'black', lwd = 1, linetype = 2, alpha = 0.6)+
    # geom_rect(data = tp, aes(fill = new), xmin = -Inf, xmax = Inf,
            # ymin = -Inf,ymax = Inf, alpha = 0.005, show.legend =FALSE, inherit.aes = FALSE)+
  geom_violin(alpha = 0.4, adjust = 1, draw_quantiles = TRUE)+
  scale_fill_manual(values = c(rep(adjustcolor("#c8c8ff",0.8),6), rep("#f8cf94",2), 
                               rep("#adf7ad",2), rep(adjustcolor("#000000",0.21),2),
                               rep("#e2a2fd",4)))+
  theme(axis.text=element_text(size=16), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=14),
        legend.background = element_rect(fill="transparent"),
        # axis.text.x = element_text(angle = 65, hjust = 1),
        strip.text.x=element_text(size=18),
        legend.position="bottom",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  guides(fill=FALSE)+
  stat_summary(fun.data=mean_sdl, fun.args = list(mult = 1), 
                 geom="pointrange", color="black")+
  xlab("")+
  ylab("")+
  scale_y_continuous(labels = function(x) sprintf("%.2f", x), breaks = seq(0.20,0.90,0.10),
                     limits = c(0.2,0.9))
  # scale_y_continuous(breaks = seq(0.20,0.90,0.10), limits = c(0.2,0.9))

p_ramli_GC
```

<img src="Figures/cached/compare-ramli-CB-2.png" style="display: block; margin: auto;" />

```r
# Plot number of gene distributions
ref_RAMLI_genes <- ref_RAMLI_SCUO  %>% dplyr::filter(GCx == "GC_mean") %>% group_by(Genome) %>% 
  summarise(Genes = n())

p_ramli_genes <- ref_RAMLI_genes  %>% 
  ggplot(aes(x = Genome, y = Genes, fill = Genome))+
  theme_bw()+
  # geom_hline(yintercept = 0.7, col = 'black', lwd = 1, linetype = 2, alpha = 0.6)+
    # geom_rect(data = tp, aes(fill = new), xmin = -Inf, xmax = Inf,
            # ymin = -Inf,ymax = Inf, alpha = 0.005, show.legend =FALSE, inherit.aes = FALSE)+
  geom_bar(stat = "identity", alpha = 0.4, color = "black", size = 1.1)+
  scale_fill_manual(values = c(rep(adjustcolor("#c8c8ff",0.8),6), rep("#f8cf94",2), 
                               rep("#adf7ad",2), rep(adjustcolor("#000000",0.21),2),
                               rep("#e2a2fd",4)))+
  theme(axis.text=element_text(size=16), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=14),
        legend.background = element_rect(fill="transparent"),
        # axis.text.x = element_text(angle = 65, hjust = 1),
        strip.text.x=element_text(size=18),
        legend.position="bottom",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  guides(fill=FALSE)+
  stat_summary(fun.data=mean_sdl, fun.args = list(mult = 1), 
                 geom="pointrange", color="black")+
  xlab("")+
  ylab("")
  # scale_y_continuous(labels = function(x) sprintf("%.2f", x), breaks = seq(0.20,0.90,0.10),
                     # limits = c(0.2,0.9))
  # scale_y_continuous(breaks = seq(0.20,0.90,0.10), limits = c(0.2,0.9))

p_ramli_genes
```

<img src="Figures/cached/compare-ramli-CB-3.png" style="display: block; margin: auto;" />


# 9.  PosiGene analysis for identifying genes under positive selection in the Ramlibacter sp. MAG

Unrooted phylogenomic tree used in codeML analysis. Green branch indicates the genome for which PSG were tested.  

![Phylogenetic tree used in PosiGene analysis](./posigene_analysis/CodeML_tree_anchor_species=Ramlibacter_MAG.png)
Alternatively, the entire Ramlibacter clade can be tested for PSGs:
Green branch indicates the clade for which the LCA was tested for PSGs.

PosiGene also tries to minimize false positives/negatives through additional filtering:  

**This filtering step can also be seen as an instrument to reduce false negatives. Few badly conserved sequences can force the first mentioned filter to delete large parts of the MSA reducing the power of the test and potentially removing positively selected sites. Third, entire MSAs can be discarded if they are considered unreliable for the following reasons, if: (i) a small absolute number or a small percentage of alignment columns or anchor species codons remain after the first filtering step, (ii) few sequences remain after the second filtering step, (iii) disproportional dN/dS ratios (e.g. ≥100 in foreground branch) were calculated by CODEML or (iv) an implausibly high fraction of positively selected sites was inferred. Additionally, MSAs will only be considered if at least one species from the sister taxon (i.e. the most closely-related species/clade) of the examined branch is represented in it. Without this condition it is not possible to say whether potentially observed selective pressure worked on the branch of interest or before in evolution .**

## 9.1. Analysis  


```r
# Import results file of genome-specific PSGs only
data_posi <- read.table("./posigene_analysis/result_tables/Ramlibacter_MAG_results.tsv", header = TRUE, fill = TRUE, sep = "\t")[, c("Transcript", "P.Value","FDR", "HA.foreground.omega", "HA.kappa",                                                                            "Number.of.Sites.under.positive.Selection")]

# Import results file of Ramlibacter clade-specific PSGs
data_posi_clade <- read.table("./posigene_analysis/result_tables_clade/Ramlibacter_MAG_results.tsv",
                              header = TRUE, fill = TRUE, sep = "\t")[, c("Transcript",
                                                                          "P.Value","FDR",
                                                                          "HA.foreground.omega", 
                                                                          "HA.kappa",                                                                        "Number.of.Sites.under.positive.Selection")]

colnames(data_posi)[1] <- "Gene"; colnames(data_posi_clade)[1] <- "Gene"
data_posi$Gene <- as.character(data_posi$Gene)
data_posi_clade$Gene <- as.character(data_posi_clade$Gene)

# Import mapping file to link gene IDs from PosiGene to 
# those used by the IMG annotation (both are in headers from .genes.fna)
map_posi <- read.table("./Mapping_files/posiG_mapping.tsv")
colnames(map_posi) <- c("posi_geneID", "IMG_geneID")
map_posi$posi_geneID <- as.character(map_posi$posi_geneID)

# Filter out the genes that are under positive selection
# Taking threshold of adjusted p.value of 0.05 and FDR < 0.05
# Also remove dN/dS ratios of less than 15.
data_posi <- data_posi %>% filter(P.Value < 0.05 & 
                                    FDR < 0.05 & HA.foreground.omega < 30)
data_posi_clade <- data_posi_clade %>% filter(P.Value < 0.05 & 
                                    FDR < 0.05 & HA.foreground.omega < 30)

# Report summaries
cat(paste("There are ", nrow(data_posi), " genes under positive selection in this MAG (P<0.05). This is ",round(100*nrow(data_posi)/nrow(map_posi),1), "% of all genes",  sep = "")
)
```

```
## There are 485 genes under positive selection in this MAG (P<0.05). This is 12.6% of all genes
```

```r
cat(paste("There are ", nrow(data_posi_clade), 
          " genes under positive selection in the Ramlibacter clade (P<0.05). This is ",
          round(100*nrow(data_posi_clade)/nrow(map_posi),1), 
          "% of all genes in the Ramlibacter sp. MAG which was used as reference/anchor species",  
          sep = "")
    )
```

```
## There are 291 genes under positive selection in the Ramlibacter clade (P<0.05). This is 7.6% of all genes in the Ramlibacter sp. MAG which was used as reference/anchor species
```

```r
cat(paste("There are ", sum(!unique(data_posi$Gene)%in% unique(data_posi_clade$Gene)), 
          " genes unique to this MAG under positive selection (P<0.05). This is ",
          round(100*sum(!unique(data_posi$Gene)%in% unique(data_posi_clade$Gene))/nrow(map_posi),1),
          "% of all genes in this MAG",  sep = "")
    )
```

```
## There are 335 genes unique to this MAG under positive selection (P<0.05). This is 8.7% of all genes in this MAG
```

```r
# Merge this data with the functional annotation (i.e. KO and COG) of these genes
data_posi <- left_join(data_posi, map_posi, by = c("Gene" = "posi_geneID"))
data_posi_clade <- left_join(data_posi_clade, map_posi, by = c("Gene" = "posi_geneID"))
data_posi_KO <- left_join(data_posi, merged_gc_ko, 
                          by = c("IMG_geneID" = "contig_geneID")) %>% distinct()
data_posi_clade_KO <- left_join(data_posi_clade, merged_gc_ko, 
                                by = c("IMG_geneID" = "contig_geneID")) %>% distinct()

# Also add COG annotation to both data sets
data_posi_COG <- left_join(data_posi, merged_gc_cog, 
                           by = c("IMG_geneID" = "contig_geneID"))

data_posi_clade_COG <- left_join(data_posi_clade, merged_gc_cog, 
                                 by = c("IMG_geneID" = "contig_geneID"))

# Also add Pfam annotation to both data sets

# Retain clade or MAG only genes in the respective dataframes
pos <- !data_posi_KO$Gene %in% data_posi_clade_KO$Gene
pos2 <- !data_posi_clade_KO$Gene %in% data_posi_KO$Gene
pos3 <- data_posi_clade_KO$Gene %in% data_posi_KO$Gene
data_posi_KO <- data_posi_KO[pos, ]
data_posi_clade_MAG_KO <- data_posi_clade_KO[pos3, ]
data_posi_clade_KO <- data_posi_clade_KO[pos2, ]

# Merge dataframes to plot
# Remove levels without "n_level" number of genes
n_level = 15

# make contigency tables for functional level C
posiG_p_df_clade  <- table(data_posi_clade_KO$ko_level_C)
posiG_p_df_clade <- data.frame(posiG_p_df_clade) 
posiG_p_df_MAG  <- table(data_posi_KO$ko_level_C)
posiG_p_df_MAG <- data.frame(posiG_p_df_MAG) 
posiG_p_df_MAG_clade  <- table(data_posi_clade_MAG_KO$ko_level_C)
posiG_p_df_MAG_clade <- data.frame(posiG_p_df_MAG_clade)

# Merge tables and extract labels of interest
test <- left_join(posiG_p_df_MAG, posiG_p_df_MAG_clade, 
                            by = "Var1")
test <- left_join(test, posiG_p_df_clade, 
                            by = "Var1")
test[is.na(test)] <- 0
test <- data.frame(labels = test$Var1, totalSums = rowSums(test[, 2:4]))
labels_2_vis <- test$labels[test$totalSums > n_level]

posiG_p_df_MAG <- posiG_p_df_MAG %>% dplyr::filter(Var1 %in% labels_2_vis)
posiG_p_df_MAG_clade <- posiG_p_df_MAG_clade %>% dplyr::filter(Var1 %in% labels_2_vis)
posiG_p_df_clade <- posiG_p_df_clade %>% dplyr::filter(Var1 %in% labels_2_vis)

# Get total of ko annotated genes per category
total_ko_genes <- table(merged_gc_ko$ko_level_C)
total_ko_genes <- data.frame(total_ko_genes)
colnames(total_ko_genes)[2] <- "TotalGenes"

# Merge dataframes
posiG_p_df_merged <- data.frame(rbind(posiG_p_df_clade, posiG_p_df_MAG, 
                                      posiG_p_df_MAG_clade),
                                branch = factor(c(rep("LCA", nrow(posiG_p_df_clade)), 
                                           rep("MAG", nrow(posiG_p_df_MAG)),
                                           rep( "LCA+MAG", nrow(posiG_p_df_MAG_clade))),
                                levels = c("MAG", "LCA",  "LCA+MAG"))
)
data_posi_KO_merge <- rbind(data_posi_KO, data_posi_clade_KO, data_posi_clade_MAG_KO)

# Add total ko annotated genes per category
posiG_p_df_merged <- dplyr::left_join(posiG_p_df_merged, total_ko_genes,
                                       by = "Var1")

# Merge with level B annotation
posiG_p_df_merged <- left_join(posiG_p_df_merged, data_posi_KO_merge[, c("ko_level_A","ko_level_B","ko_level_C")],
                       by = c("Var1" = "ko_level_C")) %>% distinct()
posiG_p_df_merged$ko_level_B[posiG_p_df_merged$ko_level_B == "Cellular community - prokaryotes"] <- "Biofilm formation & quorum sensing"

# Sort according to frequency
posiG_p_df_merged$Var1 <- factor(posiG_p_df_merged$Var1, levels = unique(posiG_p_df_merged$Var1[rev(order(posiG_p_df_merged$Freq))]))

posiG_p_df_merged$ko_level_B <- factor(posiG_p_df_merged$ko_level_B, levels = unique(posiG_p_df_merged$ko_level_B[rev(order(posiG_p_df_merged$Freq))]))

# Add extra column to shorten number of labels in legend
# Only show 12 most frequent categories
posiG_p_df_merged <- posiG_p_df_merged %>% mutate(ko_level_C_short = 
                                                    Var1)
tmp <- levels(posiG_p_df_merged$Var1)[1:11]
posiG_p_df_merged$ko_level_C_short <- as.character(posiG_p_df_merged$ko_level_C_short)
posiG_p_df_merged$ko_level_C_short[!posiG_p_df_merged$ko_level_C_short %in% tmp] <- "Other"
posiG_p_df_merged$ko_level_C_short <- factor(posiG_p_df_merged$ko_level_C_short,
                                          levels = c(tmp,"Other"))

# selected_KO <- c("Membrane transport", "Amino acid metabolism",
#                  "Amino acid metabolism", "Translation",
#                  "Lipid metabolism", "Energy metabolism",
#                  "Metabolism of cofactors and vitamins")
# Make plot
posiG_p_df_merged$branch <- factor(as.character(posiG_p_df_merged$branch),
                                   levels = c("MAG", "LCA+MAG","LCA"))
# Get order of ko_level_b correct based on summed frequencies
posiG_p_df_merged <- posiG_p_df_merged %>% 
  # dplyr::filter(branch == "MAG") %>% 
  # dplyr::filter(ko_level_B %in% selected_KO) %>%
  arrange(ko_level_B, Freq) %>%               # sort your dataframe
  mutate(Var1 = factor(Var1, unique(Var1))) # reset your factor-column based on that order

# order_var1 <- posiG_p_df_merged %>% group_by(Var1, ko_level_B) %>% summarise(sum(Freq))
# order_var1 <- as.character(order_var1[rev(order(order_var1$`sum(Freq)`)), ]$Var1)
# posiG_p_df_merged$Var1 <- factor(as.character(posiG_p_df_merged$Var1),
#                                  levels = order_var1)

# Only select the functions for which there has been selection in the MAG
relevant_Var1 <- posiG_p_df_merged %>% dplyr::filter(branch == "MAG")

# Order withing ko_level_b based on summed frequency (for visualization)
within_ko_order <- posiG_p_df_merged %>% 
  # dplyr::filter(Var1 %in% unique(relevant_Var1$Var1)) %>%
  group_by(Var1) %>% 
  mutate(sum_freq = sum(Freq))

# Make plot
p_KO_posi <- within_ko_order %>% 
  # dplyr::filter(Var1 %in% unique(relevant_Var1$Var1)) %>%
  droplevels() %>% 
  group_by(ko_level_B) %>% 
  dplyr::arrange(desc(sum_freq), .by_group = TRUE) %>% 
  ungroup() %>% 
  mutate(Var1 = as.character(Var1)) %>% 
  mutate(Var1 = factor(Var1, unique(Var1))) %>% 
  ggplot(aes(x = Var1, y = Freq, fill = branch, group = branch))+
  geom_bar(stat="identity", color = "black")+
  theme_bw()+
  scale_fill_manual(values=c(col_RAMLI, brewer.pal(11, "BrBG")[c(7,10)]))+
  ylab("PSG frequency") + xlab("")+
  theme(axis.text.y=element_text(size=12.5), axis.title=element_text(size = 18),
        legend.text=element_text(size=18),
        legend.background = element_rect(fill="transparent"),
        strip.text=element_text(size=18),
        plot.margin = unit(c(1,1,1,1), "cm"), legend.title = element_blank()
        ,legend.position = c(0.87, 0.9)
        )+
  scale_y_continuous(breaks = seq(0,80,20), limits = c(0,80))+
  # facet_grid(.~ko_level_B)+
  coord_flip()+
  theme(axis.text.x = element_text(size=18, angle = 0))+
  theme(axis.line = element_line(size = 1, colour = "grey80"),
        panel.border = element_blank())

print(p_KO_posi)
```

<img src="Figures/cached/posigene-selection-1.png" style="display: block; margin: auto;" />

```r
# Make plot
p_KO_posi_relative <- within_ko_order %>% 
  # dplyr::filter(Var1 %in% unique(relevant_Var1$Var1)) %>%
  droplevels() %>% 
  group_by(ko_level_B) %>% 
  dplyr::arrange(desc(sum_freq), .by_group = TRUE) %>% 
  ungroup() %>% 
  mutate(Var1 = as.character(Var1)) %>% 
  mutate(Var1 = factor(Var1, unique(Var1))) %>% 
  ggplot(aes(x = Var1, y = 100*Freq/TotalGenes, fill = branch, group = branch))+
  geom_bar(stat="identity", color = "black")+
  theme_bw()+
  scale_fill_manual(values=c(col_RAMLI, brewer.pal(11, "BrBG")[c(7,10)]))+
  ylab("PSG frequency (%)") + xlab("")+
  theme(axis.text.y=element_text(size=12.5), axis.title=element_text(size = 18),
        legend.text=element_text(size=18),
        legend.background = element_rect(fill="transparent"),
        strip.text=element_text(size=18),
        plot.margin = unit(c(1,1,1,1), "cm"), legend.title = element_blank()
        ,legend.position = c(0.87, 0.9)
        )+
  scale_y_continuous(breaks = seq(0,80,20), limits = c(0,80))+
  # facet_grid(.~ko_level_B)+
  coord_flip()+
  theme(axis.text.x = element_text(size=18, angle = 0))+
  theme(axis.line = element_line(size = 1, colour = "grey80"),
        panel.border = element_blank())

print(p_KO_posi_relative)
```

<img src="Figures/cached/posigene-selection-2.png" style="display: block; margin: auto;" />

```r
within_ko_order_sb <- posiG_p_df_merged %>% 
  dplyr::filter(Var1 %in% unique(relevant_Var1$Var1)) %>%
  group_by(Var1) %>% 
  mutate(sum_freq = sum(Freq))

# Make plot
p_KO_posi_subset <- within_ko_order_sb %>% 
  dplyr::filter(Var1 %in% unique(relevant_Var1$Var1)) %>%
  droplevels() %>% 
  group_by(ko_level_B) %>% 
  dplyr::arrange(desc(sum_freq), .by_group = TRUE) %>% 
  ungroup() %>% 
  mutate(Var1 = as.character(Var1)) %>% 
  mutate(Var1 = factor(Var1, unique(Var1))) %>% 
  ggplot(aes(x = Var1, y = Freq, fill = branch, group = branch))+
  geom_bar(stat="identity", color = "black")+
  theme_bw()+
  scale_fill_manual(values=c(col_RAMLI, brewer.pal(11, "BrBG")[c(7,10)]))+
  ylab("PSG frequency") + xlab("")+
  theme(axis.text.y=element_text(size=15), axis.title=element_text(size = 20),
        legend.text=element_text(size=18),
        legend.background = element_rect(fill="transparent"),
        strip.text=element_text(size=18),
        plot.margin = unit(c(1,1,1,1), "cm"), legend.title = element_blank()
        ,legend.position = c(0.87, 0.9)
        )+
  scale_y_continuous(breaks = seq(0,80,20), limits = c(0,80))+
  # facet_grid(.~ko_level_B)+
  coord_flip()+
  theme(axis.text.x = element_text(size=18, angle = 0))+
  theme(axis.line = element_line(size = 1, colour = "grey80"),
        panel.border = element_blank())


print(p_KO_posi_subset)
```

<img src="Figures/cached/posigene-selection-3.png" style="display: block; margin: auto;" />

```r
p_KO_posi_subset_relative <- within_ko_order_sb %>% 
  dplyr::filter(Var1 %in% unique(relevant_Var1$Var1)) %>%
  droplevels() %>% 
  group_by(ko_level_B) %>% 
  dplyr::arrange(desc(sum_freq), .by_group = TRUE) %>% 
  ungroup() %>% 
  mutate(Var1 = as.character(Var1)) %>% 
  mutate(Var1 = factor(Var1, unique(Var1))) %>% 
  ggplot(aes(x = Var1, y = 100*Freq/TotalGenes, fill = branch, group = branch))+
  geom_bar(stat="identity", color = "black")+
  theme_bw()+
  scale_fill_manual(values=c(col_RAMLI, brewer.pal(11, "BrBG")[c(7,10)]))+
  ylab("PSG frequency (%)") + xlab("")+
  theme(axis.text.y=element_text(size=15), axis.title=element_text(size = 20),
        legend.text=element_text(size=18),
        legend.background = element_rect(fill="transparent"),
        strip.text=element_text(size=18),
        plot.margin = unit(c(1,1,1,1), "cm"), legend.title = element_blank()
        ,legend.position = c(0.87, 0.9)
        )+
  scale_y_continuous(breaks = seq(0,80,20), limits = c(0,80))+
  # facet_grid(.~ko_level_B)+
  coord_flip()+
  theme(axis.text.x = element_text(size=18, angle = 0))+
  theme(axis.line = element_line(size = 1, colour = "grey80"),
        panel.border = element_blank())

print(p_KO_posi_subset_relative)
```

<img src="Figures/cached/posigene-selection-4.png" style="display: block; margin: auto;" />



```r
# Also add COG annotation to both data sets
data_posi_COG <- left_join(data_posi, merged_gc_cog, 
                           by = c("IMG_geneID" = "contig_geneID"))

data_posi_clade_COG <- left_join(data_posi_clade, merged_gc_cog, 
                                 by = c("IMG_geneID" = "contig_geneID"))


# Same analysis as above but for the COG annotated PSG dataset
# Retain clade or MAG only genes in the respective dataframes
pos <- !data_posi_COG$Gene %in% data_posi_clade_COG$Gene
pos2 <- !data_posi_clade_COG$Gene %in% data_posi_COG$Gene
pos3 <- data_posi_clade_COG$Gene %in% data_posi_COG$Gene
data_posi_COG <- data_posi_COG[pos, ]
data_posi_clade_COG <- data_posi_clade_COG[pos2, ]
data_posi_clade_MAG_COG <- data_posi_clade_COG[pos3, ]

# Optional: write table for quick view in iPath v2
# write.table(file = "KO_posiG.tsv", unique(data_posi_KO$ko_id), quote = FALSE,
#             row.names = FALSE, col.names = FALSE)
# write.table(file = "KO_posiG_clade.tsv", unique(data_posi_clade_KO$ko_id), quote = FALSE,
#             row.names = FALSE, col.names = FALSE)

# Merge dataframes to plot
# Remove levels without "n_level" number of genes
n_level <- round(0.01*sum(table(data_posi_clade_COG$COG_functional_category)),0)
posiG_p_df_clade  <- table(data_posi_clade_COG$COG_functional_category)[table(data_posi_clade_COG$COG_functional_category)>n_level]
posiG_p_df_clade <- data.frame(posiG_p_df_clade); posiG_p_df_clade$Var1 <- as.character(posiG_p_df_clade$Var1)

n_level <- round(0.01*sum(table(data_posi_COG$COG_functional_category)),0)
posiG_p_df_MAG  <- table(data_posi_COG$COG_functional_category)[table(data_posi_COG$COG_functional_category)>n_level]
posiG_p_df_MAG <- data.frame(posiG_p_df_MAG); posiG_p_df_MAG$Var1 <- as.character(posiG_p_df_MAG$Var1)

n_level <- round(0.01*sum(table(data_posi_clade_MAG_COG$COG_functional_category)),0)
posiG_p_df_MAG_clade  <- table(data_posi_clade_MAG_COG$COG_functional_category)[table(data_posi_clade_MAG_COG$COG_functional_category)>n_level]
posiG_p_df_MAG_clade <- data.frame(posiG_p_df_MAG_clade); posiG_p_df_MAG_clade$Var1 <- as.character(posiG_p_df_MAG_clade$Var1)

# Merge dataframes
posiG_p_df_merged <- data.frame(rbind(posiG_p_df_clade, posiG_p_df_MAG, 
                                      posiG_p_df_MAG_clade),
                                branch = factor(c(rep("LCA", nrow(posiG_p_df_clade)), 
                                           rep("MAG", nrow(posiG_p_df_MAG)),
                                           rep( "LCA+MAG", nrow(posiG_p_df_MAG_clade))),
                                levels = c("MAG", "LCA",  "LCA+MAG"))
)
data_posi_cog_merge <- rbind(data_posi_COG, data_posi_clade_COG, data_posi_clade_MAG_COG)

# Merge with level B annotation
posiG_p_df_merged <- left_join(posiG_p_df_merged, data_posi_cog_merge[, 
                                                                      c("COG_functional_category",
                                                                        "COG_functional_cluster",
                                                                        "COG_class")],
                       by = c("Var1" = "COG_functional_category")) %>% distinct()

# posiG_p_df_merged$ko_level_B[posiG_p_df_merged$ko_level_B == "Cellular community - prokaryotes"] <- "Biofilm formation & quorum sensing"

# Sort according to frequency
posiG_p_df_merged$Var1 <- factor(posiG_p_df_merged$Var1, levels = unique(posiG_p_df_merged$Var1[rev(order(posiG_p_df_merged$Freq))]))

# posiG_p_df_merged$ko_level_B <- factor(posiG_p_df_merged$COG_functional_cluster, levels = unique(posiG_p_df_merged$ko_level_B[rev(order(posiG_p_df_merged$Freq))]))

# Add extra column to shorten number of labels in legend
# Only show 12 most frequent categories
# posiG_p_df_merged <- posiG_p_df_merged %>% mutate(ko_level_C_short = 
#                                                     Var1)
# tmp <- levels(posiG_p_df_merged$Var1)[1:11]
# posiG_p_df_merged$ko_level_C_short <- as.character(posiG_p_df_merged$ko_level_C_short)
# posiG_p_df_merged$ko_level_C_short[!posiG_p_df_merged$ko_level_C_short %in% tmp] <- "Other"
# posiG_p_df_merged$ko_level_C_short <- factor(posiG_p_df_merged$ko_level_C_short,
#                                           levels = c(tmp,"Other"))

# selected_KO <- c("Membrane transport", "Amino acid metabolism",
#                  "Amino acid metabolism", "Translation",
#                  "Lipid metabolism", "Energy metabolism",
#                  "Metabolism of cofactors and vitamins")
# Make plot
p_cog_posi <- posiG_p_df_merged %>% 
  # dplyr::filter(ko_level_B %in% selected_KO) %>% 
  ggplot(aes(x = Var1, y = Freq, fill = COG_functional_cluster))+
  geom_bar(stat="identity", color = "black")+
  theme_bw()+
  scale_fill_brewer(palette="Paired")+
  ggtitle("Number of genes")+
  ylab("") + xlab("")+
  facet_grid(branch~.)+
  theme(axis.text=element_text(size=12.5), axis.title=element_text(18),
        title=element_text(size=18), legend.text=element_text(size=14),
        legend.background = element_rect(fill="transparent"),
        axis.text.x = element_text(angle = 65, hjust = 1),
        strip.text=element_text(size=18),
        plot.margin = unit(c(1,1,1,1), "cm"), legend.title = element_blank()
        # ,legend.position = c(0.87, 0.85)
        )

print(p_cog_posi)
```

<img src="Figures/cached/posigene-selection-2-1.png" style="display: block; margin: auto;" />

As we had significant concerns related that there was a %GC dependency of dN/dS (`HA.foreground.omega`) and the associated `P.value` for the branch-site codon model we also checked this specifically for the branch- and clade-tests.   
**Conclusions: There was no strong pattern/correlation to be observed between** 
**%GC and dN/dS and the P-value** 


```r
data_GC_posi <- SCUO_merged_gen
data_GC_posi$posi_select <- data_GC_posi$Gene %in% data_posi$IMG_geneID
data_GC_posi$posi_select <- factor(data_GC_posi$posi_select, levels = c("FALSE", "TRUE"))

data_GC_posi <- left_join(data_GC_posi,
                            data_posi, by = c("Gene" = "IMG_geneID"))

# Visualize differences in codon bias per codon position
p_SCUO_GC <- data_GC_posi %>% filter(GCx != "GC_mean", Genome_ID == "Ramlibacter sp. MAG") %>% 
  ggplot(aes(x = posi_select, y = 100*GC))+
  geom_jitter(size = 4, shape = 21, alpha = 0.1, width = 0.2, fill = col_RAMLI)+
  geom_boxplot(alpha = 0.2, size = 1.05, color = "darkorange")+
  scale_fill_manual("", values = c(col_RAMLI))+
  theme_bw()+
  facet_wrap(~GCx, ncol = 3)+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=14),
        legend.background = element_rect(fill="transparent"),
        strip.text.x=element_text(size=18),
        legend.position = "bottom")+
  ylab("%GC")+
  xlab("Gene under positive selection (PSG)")+
  guides(fill = FALSE)+
  ylim(0,100)+
  ggtitle("Ramlibacter sp. MAG")

print(p_SCUO_GC)
```

<img src="Figures/cached/posigene-test-GC-1.png" style="display: block; margin: auto;" />


```r
# Correlate genes under positive selection with their GC content
data_SCUO_posi <- SCUO_merged_gen_gcmean %>% filter(Genome_ID == "Ramlibacter sp. MAG")
data_SCUO_posi <- droplevels(data_SCUO_posi)
data_SCUO_posi$posi_select <- data_SCUO_posi$Gene %in% data_posi$IMG_geneID
data_SCUO_posi$posi_select <- factor(data_SCUO_posi$posi_select, levels = c("FALSE", "TRUE"))

# merge with dN/dS ratio data
data_SCUO_posi <- left_join(data_SCUO_posi,
                            data_posi, by = c("Gene" = "IMG_geneID"))
data_SCUO_posi_only <- data_SCUO_posi %>% filter(posi_select == "TRUE")

# Selected annotation level (KO)
selected_KOlevels <- c("Signal transduction",
"Amino acid metabolism",
"Metabolism of other amino acids",
"Metabolism of cofactors and vitamins",
"Carbohydrate metabolism",
"Energy metabolism",
"Translation",
"Membrane transport",
"Biosynthesis of other secondary metabolites",
"Lipid metabolism",
"Membrane transport",
"Unknown"
)

# Extract ko_level_C labels of genes in top 5 dN/ds ratio
# First filter out the genes with multiple annotations by
# selecting the annotation with the highest identity.
# data_posi_KO_unique <- data_posi_KO %>% group_by(Gene) %>%
  # filter(percent_identity == max(percent_identity))


# Remove everything before first space
# data_posi_KO[, 24:26] <- apply(data_posi_KO[, 24:26], 2 ,  
#                                function(x) sub(".*? (.+)", "\\1", x))
# Remove everything between brackets ([*])
# data_posi_KO[, 24:26] <- apply(data_posi_KO[, 24:26], 2 ,  
#                                function(x) gsub(pattern = " *\\[.*?\\] *", 
#                                  replacement = "", x))

top <- 10
thresh <- sort(data_posi_KO[data_posi_KO$ko_level_B %in% selected_KOlevels, ]$HA.foreground.omega, decreasing = TRUE)[top+1]
selected_KOlevels_labels <- data_posi_KO[data_posi_KO$ko_level_B %in% selected_KOlevels, ]$ko_level_C
selected_KOlevels_labels[data_posi_KO[data_posi_KO$ko_level_B %in% selected_KOlevels, ]$HA.foreground.omega < thresh] <- ""

p_SCUO_posi1 <- data_posi_KO[data_posi_KO$ko_level_B %in% selected_KOlevels, ] %>% 
  ggplot(aes(x = ko_level_B, y = HA.foreground.omega,
                                                   fill = ko_level_B))+
  # geom_jitter(shape = 21, aes(fill = posi_select, size = posi_select), width = 0.2)+
  geom_jitter(shape = 21, size = 3, width = 0.2, alpha = 0.5)+
  geom_boxplot(alpha=0.5, size =1, color = "black", width = 0.35,
               outlier.shape = NA)+
  scale_fill_brewer(palette = "Set3")+
  # scale_fill_manual(values = c(adjustcolor(col_RAMLI, 0.1), adjustcolor("red", 0.5)))+
  # scale_size_manual(values = c(2.5,4))+
  theme_bw()+
  # facet_wrap(Genome_ID~GCx)+
  theme(axis.text=element_text(size=12.5), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=14),
        legend.background = element_rect(fill="transparent"),
        axis.text.x = element_text(angle = 60, hjust = 1),
        strip.text.x=element_text(size=18),
        legend.position = "bottom")+
  xlab("")+ 
  ylab(expression(dN/dS))+
  ylim(0,20)+
  guides(size = FALSE, fill = FALSE)+
  geom_label_repel(aes(label = selected_KOlevels_labels), fontface = 'bold', color = 'black',
                   box.padding = 0.35, point.padding = 0.5,
                   segment.color = "#3F4656",
                   size = 3,
                       # Width of the line segments.
                   segment.size = 1,
                   # Draw an arrow from the label to the data point.
                   arrow = arrow(length = unit(0.02, 'npc'), type = "closed"),
                   nudge_x = -0.1,
                   nudge_y = 0.6,
                   force = 10
                   
  )

p_SCUO_posi1
```

<img src="Figures/cached/posigene-scuo-1.png" style="display: block; margin: auto;" />

```r
# Extract ko_level_C labels of genes in top 5 dN/ds ratio
top <- 10
thresh <- sort(data_posi_KO[data_posi_KO$ko_level_B %in% selected_KOlevels, ]$HA.foreground.omega, decreasing = TRUE)[top+1]
selected_KOlevels_labels <- data_posi_KO[data_posi_KO$ko_level_B %in% selected_KOlevels, ]$ko_name
selected_KOlevels_labels[data_posi_KO[data_posi_KO$ko_level_B %in% selected_KOlevels, ]$HA.foreground.omega < thresh] <- ""

p_SCUO_posi2 <- data_posi_KO[data_posi_KO$ko_level_B %in% selected_KOlevels, ] %>% 
  ggplot(aes(x = ko_level_B, y = HA.foreground.omega,
                                                   fill = ko_level_B))+
  # geom_jitter(shape = 21, aes(fill = posi_select, size = posi_select), width = 0.2)+
  geom_jitter(shape = 21, size = 3, width = 0.2, alpha = 0.5)+
  geom_boxplot(alpha=0.5, size =1, color = "black", width = 0.35,
               outlier.shape = NA)+
  scale_fill_brewer(palette = "Set3")+
  # scale_fill_manual(values = c(adjustcolor(col_RAMLI, 0.1), adjustcolor("red", 0.5)))+
  # scale_size_manual(values = c(2.5,4))+
  theme_bw()+
  # facet_wrap(Genome_ID~GCx)+
  theme(axis.text=element_text(size=12.5), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=14),
        legend.background = element_rect(fill="transparent"),
        axis.text.x = element_text(angle = 60, hjust = 1),
        strip.text.x=element_text(size=18),
        legend.position = "bottom")+
  xlab("")+ 
  ylab(expression(dN/dS))+
  ylim(0,20)+
  guides(size = FALSE, fill = FALSE)+
  geom_label_repel(aes(label = selected_KOlevels_labels), fontface = 'bold', color = 'black',
                   box.padding = 0.35, point.padding = 0.5,
                   segment.color = "#3F4656",
                   size = 3,
                       # Width of the line segments.
                   segment.size = 1,
                   # Draw an arrow from the label to the data point.
                   arrow = arrow(length = unit(0.02, 'npc'), type = "closed"),
                   nudge_x = -0.1,
                   nudge_y = 0.6,
                   force = 10
                   
  )

p_SCUO_posi2
```

<img src="Figures/cached/posigene-scuo-2.png" style="display: block; margin: auto;" />

```r
# Test to see if there is correlation %GC and dN/dS
p_posi_test_GC <- ggplot(data_posi_KO, aes(x = GC, y = HA.foreground.omega))+
  # geom_jitter(shape = 21, aes(fill = posi_select, size = posi_select), width = 0.2)+
  geom_point(shape = 21, size = 3, alpha = 1, fill = "black")+
  scale_fill_brewer(palette = "Set3")+
  # scale_fill_manual(values = c(adjustcolor(col_RAMLI, 0.1), adjustcolor("red", 0.5)))+
  # scale_size_manual(values = c(2.5,4))+
  theme_bw()+
  # facet_wrap(Genome_ID~GCx)+
  theme(axis.text=element_text(size=12.5), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=14),
        legend.background = element_rect(fill="transparent"),
        axis.text.x = element_text(angle = 60, hjust = 1),
        strip.text.x=element_text(size=18),
        legend.position = "bottom")+
  xlab("%GC")+ 
  ylab(expression(dN/dS))+
  # ylim(0,20)+
  guides(size = FALSE, fill = FALSE)+
  geom_smooth(color = "black", size = 2)

p_posi_test_GC
```

```
## `geom_smooth()` using method = 'loess'
```

<img src="Figures/cached/posigene-scuo-3.png" style="display: block; margin: auto;" />

```r
# Do the same thing for the Ramlibacter clade-level PSGs
# Remove everything before first space
data_posi_clade_KO[, 24:26] <- apply(data_posi_clade_KO[, 24:26], 2 ,  
                               function(x) sub(".*? (.+)", "\\1", x))
# Remove everything between brackets ([*])
data_posi_clade_KO[, 24:26] <- apply(data_posi_clade_KO[, 24:26], 2 ,  
                               function(x) gsub(pattern = " *\\[.*?\\] *", 
                                 replacement = "", x))

# Test to see if there is correlation %GC and dN/dS
p_posi_test_GC_clade <- ggplot(data_posi_clade_KO, aes(x = GC, y = HA.foreground.omega))+
  # geom_jitter(shape = 21, aes(fill = posi_select, size = posi_select), width = 0.2)+
  geom_point(shape = 21, size = 3, alpha = 1, fill = "black")+
  scale_fill_brewer(palette = "Set3")+
  # scale_fill_manual(values = c(adjustcolor(col_RAMLI, 0.1), adjustcolor("red", 0.5)))+
  # scale_size_manual(values = c(2.5,4))+
  theme_bw()+
  # facet_wrap(Genome_ID~GCx)+
  theme(axis.text=element_text(size=12.5), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=14),
        legend.background = element_rect(fill="transparent"),
        axis.text.x = element_text(angle = 60, hjust = 1),
        strip.text.x=element_text(size=18),
        legend.position = "bottom")+
  xlab("%GC")+ 
  ylab(expression(dN/dS))+
  # ylim(0,20)+
  guides(size = FALSE, fill = FALSE)+
  geom_smooth(color = "black", size = 2)

print(p_posi_test_GC_clade)
```

```
## `geom_smooth()` using method = 'loess'
```

<img src="Figures/cached/posigene-scuo-4.png" style="display: block; margin: auto;" />

```r
top <- 10
thresh <- sort(data_posi_clade_KO[data_posi_clade_KO$ko_level_B %in% selected_KOlevels, ]$HA.foreground.omega, decreasing = TRUE)[top+1]
selected_KOlevels_labels <- data_posi_clade_KO[data_posi_clade_KO$ko_level_B %in% selected_KOlevels, ]$ko_level_C
selected_KOlevels_labels[data_posi_clade_KO[data_posi_clade_KO$ko_level_B %in% selected_KOlevels, ]$HA.foreground.omega < thresh] <- ""

p_SCUO_posi3 <- data_posi_clade_KO[data_posi_clade_KO$ko_level_B %in% selected_KOlevels, ] %>% 
  ggplot(aes(x = ko_level_B, y = HA.foreground.omega,
                                                   fill = ko_level_B))+
  # geom_jitter(shape = 21, aes(fill = posi_select, size = posi_select), width = 0.2)+
  geom_jitter(shape = 21, size = 3, width = 0.2, alpha = 0.5)+
  geom_boxplot(alpha=0.5, size =1, color = "black", width = 0.35,
               outlier.shape = NA)+
  scale_fill_brewer(palette = "Set3")+
  # scale_fill_manual(values = c(adjustcolor(col_RAMLI, 0.1), adjustcolor("red", 0.5)))+
  # scale_size_manual(values = c(2.5,4))+
  theme_bw()+
  # facet_wrap(Genome_ID~GCx)+
  theme(axis.text=element_text(size=12.5), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=14),
        legend.background = element_rect(fill="transparent"),
        axis.text.x = element_text(angle = 60, hjust = 1),
        strip.text.x=element_text(size=18),
        legend.position = "bottom")+
  xlab("")+ 
  ylab(expression(dN/dS))+
  ylim(0,20)+
  guides(size = FALSE, fill = FALSE)+
  geom_label_repel(aes(label = selected_KOlevels_labels), fontface = 'bold', color = 'black',
                   box.padding = 0.35, point.padding = 0.5,
                   segment.color = "#3F4656",
                   size = 3,
                       # Width of the line segments.
                   segment.size = 1,
                   # Draw an arrow from the label to the data point.
                   arrow = arrow(length = unit(0.02, 'npc'), type = "closed"),
                   nudge_x = -0.1,
                   nudge_y = 0.6,
                   force = 10
                   
  )

p_SCUO_posi3
```

<img src="Figures/cached/posigene-scuo-5.png" style="display: block; margin: auto;" />

```r
# Extract ko_level_C labels of genes in top 5 dN/ds ratio
top <- 10
thresh <- sort(data_posi_clade_KO[data_posi_clade_KO$ko_level_B %in% selected_KOlevels, ]$HA.foreground.omega, decreasing = TRUE)[top+1]
selected_KOlevels_labels <- data_posi_clade_KO[data_posi_clade_KO$ko_level_B %in% selected_KOlevels, ]$ko_name
selected_KOlevels_labels[data_posi_clade_KO[data_posi_clade_KO$ko_level_B %in% selected_KOlevels, ]$HA.foreground.omega < thresh] <- ""

p_SCUO_posi4 <- data_posi_clade_KO[data_posi_clade_KO$ko_level_B %in% selected_KOlevels, ] %>% 
  ggplot(aes(x = ko_level_B, y = HA.foreground.omega,
                                                   fill = ko_level_B))+
  # geom_jitter(shape = 21, aes(fill = posi_select, size = posi_select), width = 0.2)+
  geom_jitter(shape = 21, size = 3, width = 0.2, alpha = 0.5)+
  geom_boxplot(alpha=0.5, size =1, color = "black", width = 0.35,
               outlier.shape = NA)+
  scale_fill_brewer(palette = "Set3")+
  # scale_fill_manual(values = c(adjustcolor(col_RAMLI, 0.1), adjustcolor("red", 0.5)))+
  # scale_size_manual(values = c(2.5,4))+
  theme_bw()+
  # facet_wrap(Genome_ID~GCx)+
  theme(axis.text=element_text(size=12.5), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=14),
        legend.background = element_rect(fill="transparent"),
        axis.text.x = element_text(angle = 60, hjust = 1),
        strip.text.x=element_text(size=18),
        legend.position = "bottom")+
  xlab("")+ 
  ylab(expression(dN/dS))+
  ylim(0,20)+
  guides(size = FALSE, fill = FALSE)+
  geom_label_repel(aes(label = selected_KOlevels_labels), fontface = 'bold', color = 'black',
                   box.padding = 0.35, point.padding = 0.5,
                   segment.color = "#3F4656",
                   size = 3,
                       # Width of the line segments.
                   segment.size = 1,
                   # Draw an arrow from the label to the data point.
                   arrow = arrow(length = unit(0.02, 'npc'), type = "closed"),
                   nudge_x = -0.1,
                   nudge_y = 0.6,
                   force = 10
                   
  )

p_SCUO_posi4
```

<img src="Figures/cached/posigene-scuo-6.png" style="display: block; margin: auto;" />


```r
# Import results file of genome-specific PSGs
data_posi2 <- read.table("./posigene_analysis/result_tables/Ramlibacter_MAG_results.tsv", header = TRUE, fill = TRUE, sep = "\t")

data_posi2 <- data_posi2 %>% filter(P.Value < 0.05 & 
                                    FDR < 0.05,
                                    HA.foreground.omega < 10)

# Remove empty columns
data_posi2 <- data_posi2 %>%
    mutate_all(funs(na_if(., ""))) %>%
    remove_empty_cols()
  
# Reformat sites under positive selection
data_posi2 <- cbind(data_posi2[c("Transcript", "HA.foreground.omega", "HA.kappa")], 
                    apply(data_posi2[, 27:105], 2, FUN = function(x) extract_aa(x)))
colnames(data_posi2)[4] <- "Site.under.positve.Selection.1"

# Get this into long dataframe format
data_posi2_long <- gather(data_posi2, Site, AA-codon, 
                          Site.under.positve.Selection.1:Site.under.positve.Selection.79, factor_key=TRUE)
# Split codon and AA
data_posi2_long <- cbind(data_posi2_long[-4], do.call(rbind, strsplit(data_posi2_long$`AA - codon`, "-"))) 
colnames(data_posi2_long)[c(5,6)] <- c("AA", "codon")

# remove NA rows
data_posi2_long <- data_posi2_long[!data_posi2_long$AA %in% "NA", ]
data_posi2_long <- droplevels(data_posi2_long)

# order AA according to abundance
order_AA <- levels(data_posi2_long$AA)[order(table(data_posi2_long$AA), decreasing = TRUE)]
data_posi2_long$AA <- as.character(data_posi2_long$AA)
data_posi2_long$AA <- factor(data_posi2_long$AA, levels = order_AA)

# Ready to make plots
p_aa_summary <- ggplot(data_posi2_long, aes(x = AA))+
  geom_bar()+
    theme(axis.text=element_text(size=13), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=14),
        legend.background = element_rect(fill="transparent"),
        # axis.text.x = element_text(angle = 65, hjust = 1),
        strip.text.x=element_text(size=18),
        legend.position="bottom")+
  xlab("Amino acid")
print(p_aa_summary)
```

<img src="Figures/cached/posigene-aa-1.png" style="display: block; margin: auto;" />

```r
p_omega <- ggplot(data_posi2_long, aes(x = AA, y = HA.foreground.omega))+
  geom_violin(alpha = 0.4, adjust = 1, draw_quantiles = TRUE)+
  stat_summary(fun.data=mean_sdl, fun.args = list(mult = 1), 
                 geom="pointrange", color="black")+
    theme(axis.text=element_text(size=13), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=14),
        legend.background = element_rect(fill="transparent"),
        # axis.text.x = element_text(angle = 65, hjust = 1),
        strip.text.x=element_text(size=18),
        legend.position="bottom")+
  xlab("Amino acid")
print(p_omega)
```

<img src="Figures/cached/posigene-aa-2.png" style="display: block; margin: auto;" />

```r
p_aa_kappa <- ggplot(data_posi2_long, aes(x = AA, y = HA.kappa))+
  geom_violin(alpha = 0.4, adjust = 1, draw_quantiles = TRUE)+
  stat_summary(fun.data=mean_sdl, fun.args = list(mult = 1), 
                 geom="pointrange", color="black")+
    theme(axis.text=element_text(size=13), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=14),
        legend.background = element_rect(fill="transparent"),
        # axis.text.x = element_text(angle = 65, hjust = 1),
        strip.text.x=element_text(size=18),
        legend.position="bottom")+
  xlab("Amino acid")
print(p_aa_kappa)
```

<img src="Figures/cached/posigene-aa-3.png" style="display: block; margin: auto;" />


## 9.2. Evaluate nucleotide transitions


```r
# Read in global nucleotide alignment file

# Read in gene/PSG site information

# Search for nucleotide composition and location of PSGs

# Store in long dataframe

# Investigate specific transitions for the whole genome as well as each gene
```

# Pangenome analysis  

* Genomes were annotated with COG ids through `anvi-run-ncbi-cogs` by `blast` searches.  
* Pangenome was made using the `--use-ncbi-blast` flag as recommended by the developers.  
* Summary of panG protein clusters were manually selected in the interactive interface.  
* The amino acid sequences in the protein cluster of the *Ramlibacter* MAG were further annotated with the KEGG orthology.  

``` 
for gene in `cat gene_list_pan.tsv`; do
	anvi-get-dna-sequences-for-gene-calls -c 121950_assembled.db --gene-caller-ids $gene -o genes_tmp.fa
	cat genes_tmp.fa >> genes_pan.fa
done
```

![Pangenome visualized in anvio](./panG/viewport.png)


```r
# Import and format panG data
panG <- read.table("./panG/SUMMARY_Ramli_PCs/panG-ramli_protein_clusters_summary.txt", header = TRUE, fill = TRUE, sep = "\t")[ , c("bin_name",	"genome_name",	"gene_callers_id",	"COG_CATEGORY_ACC",	"COG_CATEGORY",	"COG_FUNCTION_ACC", "COG_FUNCTION", "aa_sequence")]
panG$aa_sequence <- gsub("-", "", panG$aa_sequence)
panG_MAG <- panG %>% filter(bin_name %in% c("MAG_PC", "Ramli_5-10_PC",
                                          "Ramli_Leaf400_PC", "Ramli_TTB310_PC",
                                          "CORE_PC", "Mixed_PCs"))
panG_MAG$gene_callers_id <- as.character(panG_MAG$gene_callers_id)

# Shorten names
panG_MAG$genome_name <- plyr::revalue(panG_MAG$genome_name,
                            replace = c("Ramlibacter_MAG" = "Ramli-MAG",
                                        "Ramlibacter_sp_Leaf400" = "Ramli-Leaf400",
                                        "Ramlibacter_sp_TTB310" = "Ramli-TTB310",
                                        "Ramlibacter_sp_5_10" = "Ramli-5-10"))

# Add unique gene identifier 
panG_MAG$unique_gene_callers_id <- interaction(panG_MAG$genome_name, 
                                               panG_MAG$gene_callers_id,
                                               sep = ".")

# Export all AA sequences for annotation with KAAS or for blast
# for(i in 1:length(unique(panG_MAG$genome_name))){
#   tmp <- panG_MAG %>% dplyr::filter(genome_name == unique(panG_MAG$genome_name)[i]) %>% 
#     droplevels()
#     write.table(file = paste0("./panG/", unique(panG_MAG$genome_name)[i],
#                               "_aa_export.fa"), 
#            paste0(">", tmp$unique_gene_callers_id, 
#       "\n", tmp$aa_sequence, sep = ""),
#       quote = FALSE, row.names = FALSE, col.names = FALSE
#       )
# }

# Export Ramlibacter sp. auxillary genome COG annotation 
# for visualization in ipath2
# write.table(paste(panG_MAG$COG_FUNCTION_ACC, "W10", "#e2a2fd", panG_MAG$genome_name,
#                   sep=" "), "./panG/cog_id_panG.txt", row.names = FALSE,
#             quote = FALSE, col.names = FALSE)

# Import KEGG annotation through KAAS (http://www.genome.jp/tools/kaas/) of amino
# acid sequences
ko_files <- list.files(".", pattern = "KO-annotation.tsv",
                       recursive = TRUE)
panG_ko <- data.frame()
for(ko_file in ko_files){
  tmp <- read.delim(ko_file, fill = TRUE)
  tmp <- data.frame(tmp, Genome = do.call(rbind, strsplit(ko_file, "_"))[,2])
  if(ko_file == ko_files[1]) panG_ko <- tmp else{
    panG_ko <- rbind(panG_ko, tmp)
  }
}
panG_ko <- panG_ko[panG_ko$ko_id != "",]
panG_ko$ko_id <- gsub(" ","", panG_ko$ko_id)
  
# Annotate KO_IDs with hierarchy
panG_ko <- dplyr::left_join(panG_ko, ko_path_df, by = "ko_id")

# join with corresponding COG ids
panG_ko_cog <- dplyr::left_join(panG_MAG, panG_ko, 
                                by = "unique_gene_callers_id")

# Shorten/change ko_level_B annotation a bit
panG_ko_cog$ko_level_B[panG_ko_cog$ko_level_B == "Cellular community - prokaryotes"] <- "Biofilm formation & quorum sensing"
panG_ko_cog$ko_level_B[panG_ko_cog$ko_level_B == "Xenobiotics biodegradation and metabolism"] <- "Xenobiotics degradation"
panG_ko_cog$ko_level_C[panG_ko_cog$ko_level_C == "Biofilm formation - Escherichia coli "] <- "Biofilm formation"
panG_ko_cog$ko_level_C[panG_ko_cog$ko_level_C == "Biofilm formation - Pseudomonas aeruginosa "] <- "Biofilm formation"

# Add column denoting whether it is core/mixed or accessory
panG_ko_cog$bin_core <- factor(panG_ko_cog$bin_name == "CORE_PC" | panG_ko_cog$bin_name == "Mixed_PCs")
panG_ko_cog$bin_core <- plyr::revalue(panG_ko_cog$bin_core, replace = c("TRUE" = "CORE/Mixed", "FALSE" = "Accessory"))
```



```r
# Plot upset plot flagellar assembly genes
panG_ko_flagel <- panG_ko_cog %>% 
  dplyr::filter(grepl("flagel", ko_function_spec)) %>% 
  dplyr::select(ko_id, bin_name, ko_function_abbrev, ko_function_spec,
                unique_gene_callers_id) %>%
  distinct()

# Make presence column
panG_ko_flagel$Presence <- 1
# panG_ko_flagel$bin_name <- gsub("-", "_", panG_ko_flagel$bin_name)

# Revalue
panG_ko_flagel$bin_name <- plyr::revalue(panG_ko_flagel$bin_name,
                                  replace = c("MAG_PC" = "R.aquaticus", 
                                              "Ramli_5-10_PC" = 'R.5.10', 
                                              "Ramli_TTB310_PC" = "R.TTB310",
                                              "Mixed_PCs" = "Mixed_PCs"))

# From long to wide format
panG_ko_flagel <- panG_ko_flagel %>% 
  group_by(ko_id, bin_name) %>% 
  mutate(sum_presence = sum(Presence)) %>% 
  select(-unique_gene_callers_id, -Presence) %>%
  distinct() %>% 
  droplevels()

panG_ko_flagel$ko_id_function <- paste0(panG_ko_flagel$ko_function_abbrev, 
                                        " (", 
                                        panG_ko_flagel$ko_id,
                                        ")")

# Make heatmap plot of flagel assembly genes
hm_flagel <- panG_ko_flagel %>%
  select(ko_id_function,bin_name,sum_presence) %>% 
  complete(., bin_name, ko_id_function, fill = list(sum_presence=0)) %>% 
  ungroup() %>% 
  mutate(bin_name = factor(bin_name, 
                           levels = c("R.aquaticus","R.5.10",
                                      "R.TTB310","Mixed_PCs"))) %>% 
  distinct() %>% 
  ggplot(aes(y = bin_name, x = ko_id_function)) + # x and y axes => Var1 and Var2
  geom_tile(aes(fill = sum_presence), col = "lightgrey") + 
  geom_text(aes(label = round(sum_presence, 0)), size = 3) + # write the values
  scale_fill_distiller(palette="YlOrRd", na.value="lightgrey", trans = "sqrt",
                       direction = 1, limits = c(0,4)) +
  theme(panel.grid.major.x=element_blank(), #no gridlines
        panel.grid.minor.x=element_blank(), 
        panel.grid.major.y=element_blank(), 
        panel.grid.minor.y=element_blank(),
        panel.background=element_rect(fill="white"), # background=white
        axis.text.x = element_text(angle=45, hjust = 1, vjust=1, size = 12,
                                   face = "bold"),
        plot.title = element_text(size=20,face="bold"),
        axis.text.y = element_text(size = 12,face = "bold"))+
  theme(legend.title=element_text(face="bold", size=14)) + 
  scale_x_discrete(name="") +
  scale_y_discrete(name="") +
  labs(fill="Gene\ncount")
```

```
## Adding missing grouping variables: `ko_id`
```

```r
# Plot upset plot chemotaxis genes
panG_ko_chemo <- panG_ko_cog %>% 
  dplyr::filter(grepl("chemotaxis", ko_level_C)) %>% 
  dplyr::select(ko_id, bin_name, ko_function_abbrev, 
                ko_function_spec, unique_gene_callers_id) %>%
  distinct()

panG_ko_chemo$ko_id_function <- paste0(panG_ko_chemo$ko_function_abbrev, 
                                        " (", 
                                        panG_ko_chemo$ko_id,
                                        ")")


# Make presence column
panG_ko_chemo$Presence <- 1
# panG_ko_chemo$bin_name <- gsub("-", "_", panG_ko_chemo$bin_name)

# Revalue
panG_ko_chemo$bin_name <- plyr::revalue(panG_ko_chemo$bin_name,
                                  replace = c("MAG_PC" = "R.aquaticus", 
                                              "Ramli_5-10_PC" = 'R.5.10', 
                                              "Ramli_TTB310_PC" = "R.TTB310",
                                              "Mixed_PCs" = "Mixed_PCs"))
# From long to wide format
panG_ko_chemo <- panG_ko_chemo %>% 
  group_by(ko_function_abbrev, bin_name) %>% 
  mutate(sum_presence = sum(Presence)) %>% 
  select(-unique_gene_callers_id, -Presence) %>%
  distinct() %>% 
  droplevels()

# Make heatmap plot of flagel assembly genes
hm_chemo <- panG_ko_chemo %>%
  select(ko_id_function,bin_name,sum_presence) %>% 
  complete(., bin_name, ko_id_function, fill = list(sum_presence=0)) %>% 
  ungroup() %>% 
  mutate(bin_name = factor(bin_name, 
                           levels = c("R.aquaticus","R.5.10",
                                      "R.TTB310","Mixed_PCs"))) %>% 
  distinct() %>% 
  ggplot(aes(y = bin_name, x = ko_id_function)) + # x and y axes => Var1 and Var2
  geom_tile(aes(fill = sum_presence), col = "lightgrey") + # background colours are mapped according to the value column
  geom_text(aes(label = round(sum_presence, 0)), size = 3) + # write the values
  scale_fill_distiller(palette="YlOrRd", na.value="lightgrey", trans = "sqrt",
                       direction = 1, limits = c(0,4)) +
  theme(panel.grid.major.x=element_blank(), #no gridlines
        panel.grid.minor.x=element_blank(), 
        panel.grid.major.y=element_blank(), 
        panel.grid.minor.y=element_blank(),
        panel.background=element_rect(fill="white"), # background=white
        axis.text.x = element_text(angle=45, hjust = 1, vjust=1, size = 12,
                                   face = "bold"),
        plot.title = element_text(size=20,face="bold"),
        axis.text.y = element_text(size = 12,face = "bold"))+
  theme(legend.title=element_text(face="bold", size=14)) + 
  scale_x_discrete(name="") +
  scale_y_discrete(name="") +
  labs(fill="Gene\ncount")
```

```
## Adding missing grouping variables: `ko_function_abbrev`
```

```r
cowplot::plot_grid(hm_flagel, hm_chemo, nrow = 2, align = "v",
                   labels = c("A.", "B."))
```

<img src="Figures/cached/panG-analysis-2-1.png" style="display: block; margin: auto;" />

### Enrichment in panG


```r
# Perform hypergeometric test to see if there is functional enrichment in the pangenome
for(i_genome in 1:length(unique(panG_ko$Genome))){
  bg_gsea <- panG_ko %>% filter(Genome == unique(panG_ko$Genome)[i_genome]) %>% 
  select(ko_level_C,unique_gene_callers_id) %>% 
  distinct()
  
  panG_gseq <- panG_ko_cog %>% 
    filter(Genome == unique(panG_ko$Genome)[i_genome] & bin_core == "Accessory") %>%
    select(ko_level_C,unique_gene_callers_id) %>% 
    distinct()
  
  panG_gsea <- enricher(gene = panG_gseq$unique_gene_callers_id,
         universe = bg_gsea$unique_gene_callers_id, 
         TERM2GENE = bg_gsea,
         pvalueCutoff = 0.05,
         qvalueCutoff = 0.2)

  tmp_gseq <- data.frame(panG_gsea@result, Genome = unique(panG_ko$Genome)[i_genome])
  if(i_genome == 1) results_gseq <- tmp_gseq else {
    results_gseq <- rbind(results_gseq, tmp_gseq)
  }
  genes_per_gseq <- panG_gseq %>% 
     filter(ko_level_C %in% panG_gsea@result$Description) %>%
     group_by(ko_level_C) %>% 
     summarize(Counts = n())
  genes_per_gseq <- data.frame(genes_per_gseq, 
                               Genome = unique(panG_ko$Genome)[i_genome])
  if(i_genome == 1) results_pergene <- genes_per_gseq else {
    results_pergene <- rbind(results_pergene, genes_per_gseq)
  }
}

# Calculate fraction of enriched genes vs. annotated panG genes
results_pergene_sum <- results_pergene %>% 
  group_by(Genome) %>% 
  summarize(enrich_counts = sum(Counts))

annotated_fraction_panG <- panG_ko_cog %>% 
  filter(bin_core == "Accessory" & !is.na(ko_level_C)) %>% 
  group_by(Genome) %>% 
  select(unique_gene_callers_id) %>% 
  distinct() %>% 
  summarize(panG_counts = n())
```

```
## Adding missing grouping variables: `Genome`
```

```r
results_pergene_sum <- left_join(results_pergene_sum, annotated_fraction_panG,
                                 by = "Genome")

# Merge actual unique gene counts per gsea category with the accessory genome sizes
results_pergene_sum <- results_pergene_sum %>% 
  mutate(frac_enrich = enrich_counts/panG_counts)

print(results_pergene_sum)
```

```
## # A tibble: 4 x 4
##   Genome        enrich_counts panG_counts frac_enrich
##   <fct>                 <int>       <int>       <dbl>
## 1 Ramli-5-10               42         236      0.178 
## 2 Ramli-Leaf400             6         148      0.0405
## 3 Ramli-MAG               154         174      0.885 
## 4 Ramli-TTB310             66         136      0.485
```

```r
# Revalue gsea results
results_gseq$Genome <- plyr::revalue(results_gseq$Genome,
                                  replace = c("Ramli-MAG" = "R.aquaticus", 
                                              "Ramli-5-10" = 'R. solisilvae 5-10', 
                                              "Ramli-Leaf400" = "R. Leaf400",
                                              "Ramli-TTB310" = "R. tatouinensis TTB310"))

# Make plot
p_panG5 <- results_gseq %>% 
  ggplot(aes(x = Description, fill = Genome, y = Count))+
  geom_bar(color = "black", stat = "identity")+
  theme_bw()+
  scale_fill_brewer(palette="Paired")+
  ylab("Number of genes") + xlab("")+
  facet_grid(.~Genome, scales = "free")+
  theme(axis.text.y=element_text(size=12.5), axis.title.y =element_text(size = 18),
        legend.background = element_rect(fill="transparent"),
        axis.text.x = element_text(size=12.5, angle = 45, hjust = 1),
        plot.margin = unit(c(1.1,1.1,1.1,1.1), "cm"),
        strip.text = element_text(size = 12))+
  guides(fill = FALSE)

print(p_panG5)
```

<img src="Figures/cached/panG-analysis-3-1.png" style="display: block; margin: auto;" />

## Enrichment of PSGs  

* No functional groups were significantly enriched in the positively selected gene pool


```r
posi_df_gsea <- data_posi_KO %>% 
  select(ko_level_C, Gene) %>% 
  filter(!is.na(ko_level_C)) %>% 
  distinct()

bg_posi_gsea <- merged_gc_ko %>% filter(genome_id == "Ramlibacter sp. MAG") %>% 
  select(ko_level_C, gene_oid) %>% 
  distinct()

posiG_gsea <- enricher(gene = posi_df_gsea$Gene,
         universe = unique(bg_posi_gsea$gene_oid), 
         TERM2GENE = bg_posi_gsea,
         pvalueCutoff = 0.05,
         qvalueCutoff = 0.2)
```

## C/N ratio of AA


```r
# Import amino acid information
meta_aa <- read.csv("./mapping_files/aa_comp.csv")

# Get amino acid names individually
aa_seq_split <- strsplit(panG_ko_cog$aa_sequence, split = "")
ncol <- max(sapply(aa_seq_split,length))
aa_seq_split <- as.data.table(lapply(1:ncol, function(i) sapply(aa_seq_split, "[", i)))

df_aa_seq_split <- data.table::data.table(bin_name = panG_ko_cog$bin_name, 
                 genome_name = panG_ko_cog$genome_name,
                 unique_gene_callers_id = panG_ko_cog$unique_gene_callers_id,
                 aa_seq_split)

remove(aa_seq_split)

# Wide to long format
df_aa_seq_split <- gather(df_aa_seq_split, codon_position,
                          aa_ID, V1:V2467, factor_key=TRUE)
# Annotate AA sequences
df_aa_seq_split <- left_join(df_aa_seq_split, meta_aa, by = c("aa_ID" = "AA_abbrev2"))
df_aa_seq_split$unique_gene_callers_id <- factor(df_aa_seq_split$unique_gene_callers_id)

# Calculate C and N content of AA residuals
data_C_N <- df_aa_seq_split %>%
  filter(!is.na(C_elem)) %>% 
  distinct() %>% 
  group_by(unique_gene_callers_id) %>% 
  summarize(N_sum = sum(N_elem),
            C_sum = sum(C_elem),
            aa_length = n())
remove(df_aa_seq_split)

# Merge this information with initial dataframe
final_df_arsc <- left_join(panG_ko_cog, data_C_N, by = "unique_gene_callers_id")

# Calculate N_ARSC and C_ARSC
final_df_arsc <- final_df_arsc %>% 
  select(unique_gene_callers_id, bin_name, genome_name, aa_length, N_sum, C_sum) %>% 
  distinct() %>% 
  mutate(N_ARSC = N_sum/aa_length, C_ARSC = C_sum/aa_length)
```



```r
p_aa_C <- final_df_arsc %>% 
  ggplot(aes(x = bin_name, y = C_ARSC))+
   geom_violin(alpha = 0.2, fill = col_RAMLI, draw_quantiles = TRUE)+
  stat_summary(fun.data=mean_sdl, fun.args = list(mult = 1), 
                 geom="pointrange", color="black", size = 1.5, alpha = 0.75)+
  xlab("")+ ylab("C-ARSC")+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
        plot.title = element_text(size = 20, hjust = 0.5))+
  labs(title = "")

p_aa_N <- final_df_arsc %>% 
  ggplot(aes(x = bin_name, y = N_ARSC))+
   geom_violin(alpha = 0.2, fill = col_RAMLI)+
  stat_summary(fun.data=mean_sdl, fun.args = list(mult = 1), 
                 geom="pointrange", color="black", size = 1.5, alpha = 0.75)+
  xlab("")+ ylab("N-ARSC")+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
        plot.title = element_text(size = 20, hjust = 0.5))+
  labs(title = "")

cowplot::plot_grid(p_aa_C, p_aa_N, nrow = 2,
                   labels = c("A","B"))
```

<img src="Figures/cached/panG-analysis-6-1.png" style="display: block; margin: auto;" />

```r
# Summary statistics
final_df_arsc %>% group_by(bin_name) %>% 
  summarize(mean_N_ARSC = mean(N_ARSC),
            mean_C_ARSC = mean(C_ARSC),
            sd_N_ARSC = sd(N_ARSC),
            sd_C_ARSC = sd(C_ARSC))
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["bin_name"],"name":[1],"type":["fctr"],"align":["left"]},{"label":["mean_N_ARSC"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["mean_C_ARSC"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["sd_N_ARSC"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["sd_C_ARSC"],"name":[5],"type":["dbl"],"align":["right"]}],"data":[{"1":"CORE_PC","2":"0.3831969","3":"2.842217","4":"0.07753881","5":"0.1783426"},{"1":"MAG_PC","2":"0.3808549","3":"2.766043","4":"0.09238682","5":"0.2231605"},{"1":"Mixed_PCs","2":"0.3764602","3":"2.820810","4":"0.07799701","5":"0.1896899"},{"1":"Ramli_5-10_PC","2":"0.3731752","3":"2.816771","4":"0.08403378","5":"0.2047465"},{"1":"Ramli_Leaf400_PC","2":"0.3892738","3":"2.760116","4":"0.09388833","5":"0.2313926"},{"1":"Ramli_TTB310_PC","2":"0.3906899","3":"2.784405","4":"0.09000061","5":"0.2242059"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>



```r
# Lets evaluate if there is a difference between the fraction under positive selection and the rest of the genome
Ramli_aa <- data.frame(readDNAStringSet("./IMG_annotation/IMG_2724679690_Ramlibacter_bin/Annotation/2724679690.genes.faa"))
```

```
## Error in readDNAStringSet("./IMG_annotation/IMG_2724679690_Ramlibacter_bin/Annotation/2724679690.genes.faa"): could not find function "readDNAStringSet"
```

```r
colnames(Ramli_aa)[1] <- "AA_sequence"
```

```
## Error in colnames(Ramli_aa)[1] <- "AA_sequence": object 'Ramli_aa' not found
```

```r
Ramli_aa <- data.frame(Gene = do.call(rbind, strsplit(rownames(Ramli_aa), " "))[,1],
                       Ramli_aa)
```

```
## Error in rownames(Ramli_aa): object 'Ramli_aa' not found
```

```r
# Get amino acid names individually
aa_seq_split <- strsplit(Ramli_aa[,2], split = "")
```

```
## Error in strsplit(Ramli_aa[, 2], split = ""): object 'Ramli_aa' not found
```

```r
ncol <- max(sapply(aa_seq_split, length))
```

```
## Error in lapply(X = X, FUN = FUN, ...): object 'aa_seq_split' not found
```

```r
aa_seq_split <- as.data.table(lapply(1:ncol, function(i) sapply(aa_seq_split, "[", i)))
```

```
## Error in lapply(X = X, FUN = FUN, ...): object 'aa_seq_split' not found
```

```r
df_aa_seq_split <- data.table::data.table(unique_gene_callers_id = Ramli_aa$Gene,
                 aa_seq_split)
```

```
## Error in data.table::data.table(unique_gene_callers_id = Ramli_aa$Gene, : object 'Ramli_aa' not found
```

```r
# Wide to long format
df_aa_seq_split <- gather(df_aa_seq_split, codon_position,
                          aa_ID, V1:V1993, factor_key=TRUE)
```

```
## Error in gather(df_aa_seq_split, codon_position, aa_ID, V1:V1993, factor_key = TRUE): object 'df_aa_seq_split' not found
```

```r
# Annotate AA sequences
df_aa_seq_split <- left_join(df_aa_seq_split, meta_aa, by = c("aa_ID" = "AA_abbrev2"))
```

```
## Error in left_join(df_aa_seq_split, meta_aa, by = c(aa_ID = "AA_abbrev2")): object 'df_aa_seq_split' not found
```

```r
# Calculate C-ARSC and N-ARSC
data_C_N <- df_aa_seq_split %>%
  filter(!is.na(C_elem)) %>%
  distinct() %>%
  group_by(unique_gene_callers_id) %>%
  summarize(N_sum = sum(N_elem),
            C_sum = sum(C_elem),
            aa_length = n())
```

```
## Error in eval(lhs, parent, parent): object 'df_aa_seq_split' not found
```

```r
remove(df_aa_seq_split)

# Calculate N_ARSC and C_ARSC
final_df_arsc_ramli <- data_C_N %>% 
  select(unique_gene_callers_id, aa_length, N_sum, C_sum) %>% 
  distinct() %>% 
  mutate(N_ARSC = N_sum/aa_length, C_ARSC = C_sum/aa_length)

# Add posigene information
final_df_arsc_ramli$PSG <- final_df_arsc_ramli$unique_gene_callers_id %in% data_posi2_long$Transcript

# Plot N-ARSC
p_posi_N_ARSC <- final_df_arsc_ramli %>% 
  ggplot(aes(x = PSG, y = N_ARSC))+
   geom_violin(alpha = 0.2, fill = col_RAMLI)+
  stat_summary(fun.data=mean_sdl, fun.args = list(mult = 1), 
                 geom="pointrange", color="black", size = 1.5, alpha = 0.75)+
  xlab("")+ ylab("N-ARSC")+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
        plot.title = element_text(size = 20, hjust = 0.5))+
  labs(title = "")


# Plot C-ARSC
p_posi_C_ARSC <- final_df_arsc_ramli %>% 
  ggplot(aes(x = PSG, y = C_ARSC))+
   geom_violin(alpha = 0.2, fill = col_RAMLI)+
  stat_summary(fun.data=mean_sdl, fun.args = list(mult = 1), 
                 geom="pointrange", color="black", size = 1.5, alpha = 0.75)+
  xlab("")+ ylab("C-ARSC")+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
        plot.title = element_text(size = 20, hjust = 0.5))+
  labs(title = "")

# bac1_aa <- data.frame(readDNAStringSet("./IMG_annotation/IMG_2724679691_Bacteroidetes_bin1/Annotation/2724679691.genes.faa"))
# bac2_aa <- data.frame(readDNAStringSet("./IMG_annotation/IMG_2724679698_Bacteroidetes_bin2/Annotation/2724679698.genes.faa"))
# 
# # Get amino acid names individually
# aa_seq_split <- strsplit(bac1_aa[,1], split = "")
# ncol <- max(sapply(aa_seq_split,length))
# aa_seq_split <- as.data.table(lapply(1:ncol, function(i) sapply(aa_seq_split, "[", i)))
# 
# # Wide to long format
# df_aa_seq_split <- gather(aa_seq_split, codon_position,
#                           aa_ID, V1:V2145, factor_key=TRUE)
# # Annotate AA sequences
# df_aa_seq_split <- left_join(df_aa_seq_split, meta_aa, by = c("aa_ID" = "AA_abbrev2"))
# 
# # Calculate C-ARSC and N-ARSC
# data_C_N <- df_aa_seq_split %>%
#   filter(!is.na(C_elem)) %>% 
#   distinct() %>% 
#   group_by(unique_gene_callers_id) %>% 
#   summarize(N_sum = sum(N_elem),
#             C_sum = sum(C_elem),
#             aa_length = n())
# remove(df_aa_seq_split)
```


# Module completeness analysis

**Using the R-package [MetQy](https://github.com/OSS-Lab/MetQy) we will assess**
**the Kegg Orthology module completeness across all Ramlibacter genomes**  

* This will require the KO annotation of the reference Ramlibacter sp. genomes
* We will do this based on the gene calling conducted in `anvi'o` and annotation of those genes using the `KAAS` webserver

***  

Citation:  
Andrea Martinez-Vernon, Fred Farrell, Orkun Soyer. _MetQy: an R package to query metabolic functions of genes and genomes._ bioRxiv 215525; doi:https://doi.org/10.1101/215525  

***  



```r
# Get list of all KO annotated genes from IMG
KO_profiles <- read.table("./IMG_annotation/STAMP_profiles/abundance_ko_38065.tab.xls", header = TRUE, sep = "\t", quote = "")

# Remove zero rows
KO_profiles <- KO_profiles[rowSums(KO_profiles[, - c(1:2)])>0, ]

# Wide to long format
KO_profiles_long <- tidyr::gather(KO_profiles, Genome, 
                                       Counts, 
                          Curvibacter_sp._ATCC:Variovorax_sp._EPS,
                          factor_key = TRUE)

KO_profiles_long$Genome <- gsub("_", " ", KO_profiles_long$Genome)
KO_profiles_long$Genome <- gsub("5.10", "5-10", KO_profiles_long$Genome)

# Remove rows with 0 counts
KO_profiles_long <- KO_profiles_long %>% dplyr::filter(Counts > 0)

# Remove frequencies (not important for module completeness estimate)
KO_profiles_long <- KO_profiles_long[, c(1:3)] %>% distinct()

# Get module completeness
modules_table <- get_mq(KO_profiles_long, genome_label = "Genome", ko_label = "Func_id")
```

```
## [1] "Curvibacter sp. ATCC"
## [1] "Curvibacter sp. PAE"
## [1] "Limnohabitans sp. 63ED37"
## [1] "Limnohabitans sp. Rim28"
## [1] "Limnohabitans sp. Rim47"
## [1] "Ramlibacter sp. MAG"
## [1] "Ramlibacter sp. Leaf400"
## [1] "Ramlibacter sp. 5-10"
## [1] "Ramlibacter sp. TTB310"
## [1] "Rhodoferax sp. T118"
## [1] "Rhodoferax sp. ED16"
## [1] "Bacteroidetes sp. MAG1"
## [1] "Bacteroidetes sp. MAG2"
## [1] "Variovorax sp. 110B"
## [1] "Variovorax sp. EPS"
```

```r
# Remove absent modules
# modules_table <- modules_table %>% dplyr::filter(module_completeness > 0)


# Select subset of genomes for phosphate scavenging
selected_genomes <- c("Ramlibacter sp. TTB310",
                      "Ramlibacter sp. Leaf400",
                      "Ramlibacter sp. 5-10",
                      "Ramlibacter sp. MAG",
                      "Bacteroidetes sp. MAG1", 
                      "Bacteroidetes sp. MAG2")

# Visualize output statistics for Carbon fixation
p_mod1 <- modules_table %>% dplyr::filter(module_completeness > 0 & Genome %in% selected_genomes, CLASS_III == "Carbon fixation") %>% 
  ggplot(aes(x = MODULE_ID, y = module_completeness, fill = Genome))+
  geom_point(shape = 21, size = 4, 
             position = position_dodge(width = 1))+ theme_bw()+
  scale_fill_brewer(palette = "Accent")+
  xlab("")+
  ylab("Module completeness")+
  ylim(0,1)+
  # facet_grid(Genome~.)+
    theme(axis.title=element_text(size=16), strip.text.x=element_text(size=16),
        legend.title=element_text(size=15), legend.text=element_text(size=14),
        axis.text.y = element_text(size=16),
        axis.text.x = element_text(size=16, angle = 300, hjust = 0),
        title=element_text(size=20),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA)
  )

print(p_mod1)
```

<img src="Figures/cached/m-genome-completeness-1.png" style="display: block; margin: auto;" />

```r
# Specific modules
mod_spec1 <- c("Cell signaling",
               "Quorum sensing 2-CRS (Qse)",
               "Quorum sensing 2-CRS (Lux)")

# Visualize output statistics for 
p_mod2 <- modules_table %>% 
  dplyr::filter(Genome %in% selected_genomes & 
                  NAME_SHORT %in% mod_spec1) %>% 
  ggplot(aes(x = NAME_SHORT, y = module_completeness, fill = Genome))+
  geom_point(shape = 21, size = 4, 
             position = position_dodge(width = 1))+ theme_bw()+
  scale_fill_brewer(palette = "Accent")+
  xlab("")+
  ylab("Module completeness")+
  ylim(0,1)+
  # facet_grid(Genome~.)+
    theme(axis.title=element_text(size=16), strip.text.x=element_text(size=16),
        legend.title=element_text(size=15), legend.text=element_text(size=14),
        axis.text.y = element_text(size=16),
        axis.text.x = element_text(size=16, angle = 0),
        title=element_text(size=20),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA)
  )


print(p_mod2)
```

<img src="Figures/cached/m-genome-completeness-2.png" style="display: block; margin: auto;" />

```r
# Specific modules
mod_spec3 <- c("Cell signaling",
               "Quorum sensing 2-CRS (Qse)",
               "Quorum sensing 2-CRS (Lux)")

# # Visualize output statistics for 
# p_mod3 <- modules_table %>% 
#   dplyr::filter(Genome %in% selected_genomes & 
#                   NAME_SHORT %in% mod_spec3) %>% 
#   ggplot(aes(x = NAME_SHORT, y = module_completeness, fill = Genome))+
#   geom_point(shape = 21, size = 4, 
#              position = position_dodge(width = 1))+ theme_bw()+
#   scale_fill_brewer(palette = "Accent")+
#   xlab("")+
#   ylab("Module completeness")+
#   ylim(0,1)+
#   # facet_grid(Genome~.)+
#     theme(axis.title=element_text(size=16), strip.text.x=element_text(size=16),
#         legend.title=element_text(size=15), legend.text=element_text(size=14),
#         axis.text.y = element_text(size=16),
#         axis.text.x = element_text(size=16, angle = 0),
#         title=element_text(size=20),
#         panel.background = element_rect(fill = "transparent",colour = NA),
#         plot.background = element_rect(fill = "transparent",colour = NA)
#   )
# 
# 
# print(p_mod3)


# Heatmap
# make heatmap for all class_III levels
# 
# for(class_i in unique(modules_table$CLASS_III)){
#   hm_mod <- modules_table %>%
#   dplyr::filter(CLASS_III %in%
#                   class_i) %>%
#   ggplot(aes(y = Genome, x= MODULE_ID)) + # x and y axes => Var1 and Var2
#   geom_tile(aes(fill = module_completeness), col = "lightgrey") + # background colours are mapped according to the value column
#   # geom_text(aes(label = round(module_completeness, 1)), size = 3) + # write the values
#   # scale_fill_gradientn(colours = terrain.colors(10), trans = "log1p")+
#   # scale_fill_gradient(low = "lightblue", high = "darkslategray", na.value="white",
#                       # trans = "log1p", limits=c(1, 40)) +
#   scale_fill_distiller(palette="YlOrRd", na.value="lightgrey",
#                        direction = 1, limits = c(0,1)) +
#   theme(panel.grid.major.x=element_blank(), #no gridlines
#         panel.grid.minor.x=element_blank(),
#         panel.grid.major.y=element_blank(),
#         panel.grid.minor.y=element_blank(),
#         panel.background=element_rect(fill="white"), # background=white
#         axis.text.x = element_text(angle=45, hjust = 1, vjust=1, size = 12,face = "bold"),
#         plot.title = element_text(size=20,face="bold"),
#         axis.text.y = element_text(size = 12,face = "bold"))+
#   theme(legend.title=element_text(face="bold", size=14)) +
#   scale_x_discrete(name="") +
#   scale_y_discrete(name="") +
#   # facet_grid(CLASS_III~.)+
#   labs(fill="Module\ncompleteness\n")+
#   ggtitle(class_i)
# class_i <- gsub("/"," or ", class_i)
# png(paste("./Figures/cached/MetQ/metQ_", class_i, ".png", sep = ""), res = 500, width = 12, height = 10, unit = "in")
# print(hm_mod)
# dev.off()
# }
```

<!-- ### panG module completeness -->
<!-- ```{r m-panG-completeness, dpi = 400, fig.width = 13, fig.height = 10, warning = FALSE} -->
<!-- # Import KO annotated accessory genome for all other Ramlibacter genomes -->
<!-- panG_annot_files <- list.files("./panG/", pattern = "_KO-annotation.tsv") -->
<!-- for(i in 1:length(panG_annot_files)){ -->
<!--   tmp <- read.table(paste("./panG/", panG_annot_files[[i]], sep = ""), sep = "\t", -->
<!--              fill = TRUE, header = FALSE, stringsAsFactors = FALSE) -->
<!--   tmp <- cbind(tmp[,2], panG_annot_files[i]) -->
<!--   colnames(tmp) <- c("ko_id", "panG_genome_id") -->

<!--   cat(paste("KO annotated", round(100*sum(tmp[,1] != "")/nrow(tmp),0),"% of accessory genome", panG_annot_files[i],"\n", sep = " ")) -->

<!--   if(i == 1) {panG_annot <- tmp } else{ -->
<!--                                           panG_annot <- rbind(panG_annot, tmp) -->
<!--   } -->
<!-- } -->
<!-- panG_annot <- data.frame(panG_annot) -->
<!-- panG_annot$panG_genome_id <- gsub("_KO-annotation.tsv", "", -->
<!--                                   panG_annot$panG_genome_id) -->
<!-- panG_annot <- panG_annot[panG_annot$ko_id !="" , ] -->
<!-- panG_annot <- panG_annot %>% dplyr::filter(ko_id != "ko_id") -->

<!-- # Format imported data for MetQy -->
<!-- for(i in unique(panG_annot$panG_genome_id)){ -->
<!--   KO_tmp <- panG_annot %>% dplyr::filter(panG_genome_id == i) -->
<!--   KO_tmp <- droplevels(KO_tmp) -->
<!--   KO_tmp_metqy <- data.frame(genome_id = unique(KO_tmp$panG_genome_id), -->
<!--                              KOs = paste(KO_tmp$ko_id, collapse = ";"), -->
<!--                              stringsAsFactors = FALSE -->
<!--   ) -->
<!--   if(i == unique(panG_annot$panG_genome_id)[1]){ -->
<!--     KO_panG_metqy <- KO_tmp_metqy -->
<!--   } else{ KO_panG_metqy <- rbind(KO_panG_metqy, KO_tmp_metqy) -->
<!--   } -->
<!--   KO_panG_metqy$KOs <- gsub(" ", "", KO_panG_metqy$KOs) -->
<!-- } -->
<!-- remove(KO_tmp) -->

<!-- # Run query -->
<!-- query_output_panG <- MetQy::query_genomes_to_modules(KO_panG_metqy, splitBy='[;]', -->
<!--                                 GENOME_ID_COL = 1, GENES_COL = 2,  -->
<!--                                 META_OUT = TRUE) -->

<!-- # Visualize differences in module completeness across accessory genomes -->
<!-- panG_modules_table <- data.frame(genome_id = rownames(query_output_panG$MATRIX), -->
<!--                                  query_output_panG$MATRIX) -->

<!-- # now reshape so that panG id is also a column -->
<!-- panG_modules_table <- tidyr::gather(panG_modules_table, module,  -->
<!--                                     module_completeness, M00001:M00822,  -->
<!--                                     factor_key=TRUE) -->
<!-- panG_modules_table <- left_join(panG_modules_table, query_output_panG$METADATA, -->
<!--                                 by = c("module" = "MODULE_ID")) -->
<!-- panG_modules_table <- panG_modules_table[, -c(9, 10)] -->

<!-- # Visualize..  -->
<!-- p_mod3 <- panG_modules_table %>% dplyr::filter(CLASS_III %in% mod_spec1,  -->
<!--                                           module_completeness > 0) %>%  -->
<!--   ggplot(aes(x = CLASS_III, y = module_completeness, fill = CLASS_II))+ -->
<!--   geom_jitter(shape = 21, size = 4, width = 0.3)+ -->
<!--   theme_bw()+ -->
<!--   # geom_boxplot(alpha = 0.4)+ -->
<!--   scale_fill_brewer(palette = "Accent")+ -->
<!--   xlab("")+ -->
<!--   ylab("Module completeness")+ -->
<!--   ylim(0,1)+ -->
<!--   facet_grid(~genome_id)+ -->
<!--     theme(axis.title=element_text(size=16), strip.text.x=element_text(size=16), -->
<!--         legend.title=element_text(size=15), legend.text=element_text(size=14), -->
<!--         axis.text.y = element_text(size=16), -->
<!--         axis.text.x = element_text(size=16, angle = 300, hjust = 0), -->
<!--         title=element_text(size=20), -->
<!--         panel.background = element_rect(fill = "transparent", colour = NA), -->
<!--         plot.background = element_rect(fill = "transparent", colour = NA) -->
<!--   ) -->

<!-- print(p_mod3) -->

<!-- # Lets focus on the two-component regulatory systems (which are highly abundant -->
<!-- # in all Ramlibacter genome annotations) -->

<!-- p_mod4 <- panG_modules_table %>% dplyr::filter(CLASS_III == "Two-component regulatory system" & -->
<!--                                                  module_completeness >0) %>%  -->
<!--   ggplot(aes(x = NAME_SHORT, y = module_completeness))+ -->
<!--   # geom_jitter(shape = 21, size = 4, width = 0.3)+ -->
<!--   theme_bw()+ -->
<!--   geom_bar(alpha = 1, stat = "identity", color = "black", fill = brewer.pal(n = 8, "Accent")[1])+ -->
<!--   xlab("")+ -->
<!--   ylab("Module completeness")+ -->
<!--   ylim(0,1)+ -->
<!--   facet_grid(genome_id~.)+ -->
<!--     theme(axis.title=element_text(size=16), strip.text.x=element_text(size=16), -->
<!--         legend.title=element_text(size=15), -->
<!--         axis.text.y = element_text(size=16), -->
<!--         axis.text.x = element_text(size=16, angle = 65, hjust = 1), -->
<!--         title=element_text(size=20), -->
<!--         panel.background = element_rect(fill = "transparent",colour = NA), -->
<!--         plot.background = element_rect(fill = "transparent",colour = NA), -->
<!--         strip.text = element_text(size = 14) -->
<!--   )+ -->
<!--   ggtitle("Two-component regulatory system")+ -->
<!--   guides(fill = FALSE) -->

<!-- print(p_mod4) -->

<!-- ``` -->

# 10. ANI analysis using `pyani`


```r
# read file
# data.ANI <- read.table("./ANI/ANIb_percentage_identity.tab")
data.ANI <- read.table("./ANI/ANIb_percentage_identity.tab")


# Replace genome names by better annotated names
map_ani <- read.delim("./Mapping_files/pyani_ref_names.tsv", stringsAsFactors = FALSE)
for(i in 1:nrow(map_ani)){
 colnames(data.ANI)[colnames(data.ANI) %in% map_ani$ani_file[i]] <- map_ani$ref_name[i]
 rownames(data.ANI)[rownames(data.ANI) %in% map_ani$ani_file[i]] <- map_ani$ref_name[i]
}

# Order y-axis according to phylogenetic tree order
ord_list_bin <- c("Lim. sp. Rim11", "Lim. sp. 103DPR2",
                  "Lim. sp. 2KL-27", "Lim. sp. Rim47",
                  "Lim. sp. II-D5", "Lim. sp. 2KL-3",
                  "Rhodo. sp. ED16","Rhodo. sp. T118",
                  "Curvi. sp. ATCC", "Curvi. sp. PAE-UM",
                  "Vario. sp. 110B", "Vario. sp. EPS",
                  "Ramli. sp. Leaf400", "Ramli. sp. TTB310",
                  "Ramli. sp. MAG", 
                  "Ramli. sp. 5-10"
                  )

# order rows and columns
data.ANI <- data.ANI[, order(match(colnames(data.ANI), ord_list_bin))]
data.ANI <- data.ANI[order(match(rownames(data.ANI), ord_list_bin)), ]

# Get upper triangle 
data.ANI <- get_upper_tri(data.ANI)

# melt to dataframe
df_pyani <- melt(as.matrix(data.ANI), na.rm = TRUE) # reshape into dataframe
names(df_pyani)[c(1:3)] <- c("bin1", "bin2", "ANI")
df_pyani[, 1] <- as.character(df_pyani[, 1]); df_pyani[, 2] <- as.character(df_pyani[, 2])
df_pyani$bin2 <- factor(df_pyani$bin2, levels = ord_list_bin)
df_pyani$bin1 <- factor(df_pyani$bin1, levels = rev(ord_list_bin))

# make plot
p_ani <- ggplot(data = df_pyani,
       aes(x = bin2, y = bin1, fill = ANI))+
  theme(axis.title=element_text(size=16), strip.text.x=element_text(size=16),
        legend.title=element_text(size=15), legend.text=element_text(size=14),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size=16, angle = 55, hjust = 0),
        title=element_text(size=20),
        plot.margin = unit(c(1.1,1.1,1.1,1.1), "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA)
  )+ 
  geom_raster()+
  scale_fill_distiller(palette = "RdBu", name = "Average\nNucleotide\nIdentity\n",
                       limits = c(0.75,1.0), oob=squish) +
  geom_text(aes(label = round(ANI, 2)), size = 4.5)+
  xlab("")+
  ylab("")+
  scale_x_discrete(position = "top") 

print(p_ani)
```

<img src="Figures/cached/ani-analysis-1.png" style="display: block; margin: auto;" />


# Plasmid reconstruction using `Recycler`


```r
files <- list.files("./recycler_plasmid/", pattern="*_length_cov")
files <- paste0("./recycler_plasmid/", files)
df_plasm_16 <- read.table(files[1], header = FALSE)
df_plasm_17 <- read.table(files[2], header = FALSE)
df_plasm_64 <- read.table(files[3], header = FALSE)
df_plasm_65 <- read.table(files[4], header = FALSE)
df_plasm <- rbind(df_plasm_16, df_plasm_17, df_plasm_64, df_plasm_65)

df_plasm <- as.character(df_plasm$V1); tmp <- df_plasm
df_plasm <- do.call(rbind, strsplit(df_plasm, "_"))
df_plasm <- df_plasm[, -c(1,3,5,7,8)]
colnames(df_plasm) <- c("NODE_ID", "Length", "avg_coverage")
df_plasm <- apply(df_plasm, 2, function(x) as.integer(x))
df_plasm <- data.frame(df_plasm, node_IDs = tmp)

# Filter out plasmid scaffolds lower than 1000 bp
p_plasm1 <- df_plasm %>% filter(Length > 1000) %>% 
  ggplot(aes(x = Length/1000, y = sqrt(avg_coverage), fill = Length))+
  geom_point(shape = 21, size = 4)+ theme_bw()+
  # scale_x_log10()+
  # scale_y_log10()+
  scale_fill_continuous()+
  ylab(expression(sqrt(Average_coverage)))+
  xlab("Length (kbp)")+
  ggtitle("Coverage-Length distribution of reconstructed plasmid fragments")

p_plasm1
```

<img src="Figures/cached/plasmid-reconstruction-1.png" style="display: block; margin: auto;" />

# Sequence discrete population analysis 

This analysis was run in order to evaluate if there were different populations in the different operational cycles. A peak at <95 % nucleotide identity may indicate the presence of a second Ramlibacter population that is concealed in the current MAG.  


```r
# Import contig_ids of each bin
contig_ids <- rbind(
  cbind(read.table("./IMG_annotation/IMG_2724679690_Ramlibacter_bin/121950.assembled.names_map", stringsAsFactors = FALSE)[,1], "Ramlibacter sp. MAG"),
  cbind(read.table("./IMG_annotation/IMG_2724679691_Bacteroidetes_bin1/121951.assembled.names_map", stringsAsFactors = FALSE)[,1], "Bacteroidetes sp. MAG1"),
  cbind(read.table("./IMG_annotation/IMG_2724679698_Bacteroidetes_bin2/121960.assembled.names_map", stringsAsFactors = FALSE)[,1], "Bacteroidetes sp. MAG2")
)
contig_ids <- data.frame(contig_ids)
colnames(contig_ids) <- c("contig_id", "bin")

# Import blast results
blast_df <- read.table("./SEQ_discrete/merged_blast.tsv", header = FALSE)
colnames(blast_df) <- c("Sample", "Contig", "Identity")

# Filter out blast hits not to MAGs of interest
blast_df <- blast_df %>% dplyr::filter(Contig %in% contig_ids$contig_id)
blast_df <- dplyr::left_join(blast_df, contig_ids, by = c("Contig" = "contig_id"))

# Bin into proportions
blast_df$Identity <- round(blast_df$Identity, 0)
blast_df_sum <- blast_df %>% group_by(Sample, bin) %>% dplyr::count(Identity)
blast_df_sum$Sample <- gsub("_blast.tsv", "", blast_df_sum$Sample)

# Normalize mapped reads per sample based on genome size
bin_sizes <- read.table("./SAMPLES-SUMMARY/general_bins_summary.txt", header = TRUE)[, c(1, 3)]
bin_sizes$bins <- plyr::revalue(bin_sizes$bins , 
                                c("BetIa_bin"="Ramlibacter sp. MAG",
                                  "bacIa_vizbin1"="Bacteroidetes sp. MAG1",
                                  "bacIa_vizbin2"="Bacteroidetes sp. MAG2")
                                )
blast_df_sum <- left_join(blast_df_sum, bin_sizes, by = c("bin" = "bins"))
blast_df_sum <- blast_df_sum %>% group_by(bin) %>%
  mutate(n_norm = 1e6*n/bin_size)

# Rename sample IDs
blast_df_sum$Sample <- plyr::revalue(blast_df_sum$Sample,
                                     c("DNA-Pieter-16"="Sample-16",
                                       "DNA-Pieter-17"="Sample-17",
                                       "DNA-Pieter-64"="Sample-64",
                                       "DNA-Pieter-65"="Sample-65")
)

# Plot sequence discrete populations
p_blast_sdisc <- blast_df_sum %>% 
     ggplot(aes(x = Identity, y = n_norm, color = bin))+
      theme_bw()+
      scale_color_manual(values = c(col_bac1, col_bac2, col_RAMLI))+
      facet_grid(bin~Sample)+
      geom_line(size = 1.5)+
      guides(color = FALSE)+
      theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=14),
        legend.background = element_rect(fill="transparent"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text=element_text(size=12))+
      ylab("Reads per Mbp")+
      xlab("Nucleotide identity (%)")+
    xlim(75,100)+
  coord_trans(y = "sqrt")

print(p_blast_sdisc)
```

<img src="Figures/cached/seq-discrete-1-1.png" style="display: block; margin: auto;" />


# Predict MGT  

Compare the growth rate with the minimum generation time estimated from the MAG
using Growthpred. Dotted line in the graph of the predicted optimal growth temperature
indicates the mean +/- st.dev of the operational temperature in the cooling water system.  


```r
MGT_df <- read.table("./Growthpred/GP_results.tsv", header = TRUE,
                     stringsAsFactors = FALSE, sep = "\t")

MGT_df <- MGT_df[-c(4:5), ] %>% droplevels()

# Order genome_ids according to the phylogenetic clustering
ord_full_list_bin <- c("Limnohabitans sp. Rim11", "Limnohabitans sp. 103DPR2",
                  "Limnohabitans sp. 2KL-27", "Limnohabitans sp. Rim47",
                  "Limnohabitans sp. II-D5", "Limnohabitans sp. 2KL-3",
                  "Rhodoferax sp. ED16","Rhodoferax sp. T118",
                  "Curvibacter sp. ATCC", "Curvibacter sp. PAE-UM",
                  "Variovorax sp. 110B", "Variovorax sp. EPS",
                  "Ramlibacter sp. Leaf400", "Ramlibacter sp. TTB310",
                  "Ramlibacter sp. MAG", 
                  "Ramlibacter sp. 5-10",
                  "Bacteroidetes sp. MAG1",
                  "Bacteroidetes sp. MAG2"
                  )

MGT_df$Genome_ID <- factor(MGT_df$Genome_ID, levels = ord_full_list_bin)

# Make barplot with st.dev to visualize MGT and optimal temperature
selected_points <- data.frame(Genome_ID = MGT_df$Genome_ID, 
                              ypos = c(rep(7.5, 3), rep(NA,3)))
p_MGT_1 <- MGT_df[-c(2:3),] %>% 
  ggplot(aes(x = Genome_ID, y = log(2)/MGT, fill = Genome_ID))+
  theme_bw()+
  geom_bar(alpha = 0.4, stat = "identity", color = "black",
           position = position_dodge(width = 1), width = 0.7)+
  scale_fill_manual(values = c(rep(adjustcolor("#c8c8ff",0.8),6), rep("#f8cf94",2), 
                               rep("#adf7ad",2), rep(adjustcolor("#000000",0.21),2),
                               rep("#e2a2fd",4), col_bac1, col_bac2))+
  theme(axis.text=element_text(size=13), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=14),
        legend.background = element_rect(fill="transparent"),
        # axis.text.x = element_text(angle = 65, hjust = 1),
        strip.text.x=element_text(size = 18),
        legend.position="bottom",
        axis.text.x=element_text(size = 9, angle =55, hjust= 1),
        axis.title.x=element_blank(),
        plot.title = element_text(hjust = 0, size=18))+
  guides(fill=FALSE)+
  geom_errorbar(aes(ymin = log(2)/(MGT - sd.MGT),
                    ymax = log(2)/(MGT + sd.MGT)), width = 0.15,
                position = position_dodge(width = 1))+
  ylab("")+
  xlab("")+
  ggtitle(expression("Specific growth rate - h"^-1))+
  ylim(0,0.75)

p_MGT_2 <- MGT_df[c(2:3),] %>% 
  ggplot(aes(x = Genome_ID, y = log(2)/MGT, fill = Genome_ID))+
  theme_bw()+
  geom_bar(alpha = 0.4, stat = "identity", color = "black",
           position = position_dodge(width = 1), width = 0.7)+
  scale_fill_manual(values = c(col_bac1, col_bac2))+
  theme(axis.text=element_text(size=11), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=14),
        legend.background = element_rect(fill="transparent"),
        # axis.text.x = element_text(angle = 65, hjust = 1),
        strip.text.x=element_text(size = 18),
        legend.position="bottom",
        axis.text.x=element_text(size = 9, angle =55, hjust= 1),
        axis.title.x=element_blank(),
        plot.title = element_text(hjust = 0, size=18))+
  guides(fill=FALSE)+
  geom_errorbar(aes(ymin = log(2)/(MGT - sd.MGT),
                    ymax = log(2)/(MGT + sd.MGT)), width = 0.15,
                position = position_dodge(width = 1))+
  ylab("")+
  xlab("")+
  ggtitle("")+
  ylim(0,0.75)

cowplot::plot_grid(p_MGT_1, p_MGT_2, align = "hv", rel_widths = c(4,1))
```

<img src="Figures/cached/MGT-1-1.png" style="display: block; margin: auto;" />

```r
selected_points <- data.frame(Genome_ID = MGT_df$Genome_ID, 
                              ypos = c(rep(35, 3), rep(NA,3)))
selected_band <- data.frame(Genome_ID = MGT_df$Genome_ID, 
                            Topt = MGT_df$Topt,
                            lwr = rep(26.89 - 1.26, 6),
                            upr = rep(26.89 + 1.26, 6),
                            y_val = 26.89)

p_Topt_1 <- MGT_df[-c(2:3),] %>% 
  ggplot(aes(x = Genome_ID, y = Topt, fill = Genome_ID))+
  theme_bw()+
  geom_bar(alpha = 0.4, stat = "identity", color = "black",
           position = position_dodge(width = 1), width = 0.7)+
  scale_fill_manual(values = c(rep(adjustcolor("#c8c8ff",0.8),6), rep("#f8cf94",2), 
                               rep("#adf7ad",2), rep(adjustcolor("#000000",0.21),2),
                               rep("#e2a2fd",4), col_bac1, col_bac2))+
  theme(axis.text=element_text(size=11), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=14),
        legend.background = element_rect(fill="transparent"),
        # axis.text.x = element_text(angle = 65, hjust = 1),
        strip.text.x=element_text(size = 18),
        legend.position="bottom",
        axis.text.x=element_text(size = 9, angle =55, hjust= 1),
        axis.title.x=element_blank(),
        plot.title = element_text(hjust = 0, size=18))+
  guides(fill=FALSE)+
  ylab("")+
  xlab("")+
  ggtitle("Optimal growth temperature (°C)")+
  ylim(0,35)
  # geom_point(data = selected_points, aes(x = Genome_ID, y = ypos),
             # shape = 25, fill = "black", col = "black", size = 3)+
  # geom_hline(yintercept = 26.89, linetype = 2, col = "black")+
  # geom_rect(data = selected_band, aes(ymin = lwr,
  #                                     ymax = upr,
  #                                     xmin = 0.5,
  #                                     xmax = 18.5), alpha=0.01, fill = "gray",
  #           col = NA)

p_Topt_2 <- MGT_df[c(2:3),] %>% 
  ggplot(aes(x = Genome_ID, y = Topt, fill = Genome_ID))+
  theme_bw()+
  geom_bar(alpha = 0.4, stat = "identity", color = "black",
           position = position_dodge(width = 1), width = 0.7)+
  scale_fill_manual(values = c(col_bac1, col_bac2))+
  theme(axis.text=element_text(size=11), axis.title=element_text(size=11),
        title=element_text(size=20), legend.text=element_text(size=14),
        legend.background = element_rect(fill="transparent"),
        # axis.text.x = element_text(angle = 65, hjust = 1),
        strip.text.x=element_text(size = 18),
        legend.position="bottom",
        axis.text.x=element_text(size = 9, angle =55, hjust= 1),
        axis.title.x=element_blank(),
        plot.title = element_text(hjust = 0, size=18))+
  guides(fill=FALSE)+
  ylab("")+
  xlab("")+
  ggtitle("")+
  ylim(0,35)

# Print plots
cowplot::plot_grid(p_Topt_1, p_Topt_2, align = "h", rel_widths = c(4,1))
```

<img src="Figures/cached/MGT-1-2.png" style="display: block; margin: auto;" />

# iRep

Estimate the population replication rate using the iRep tool.


```r
df_irep <- read.table("./iRep/results_irep.tsv", header = TRUE, stringsAsFactors = FALSE,
                      sep = "\t") %>% 
  dplyr::filter(iRep_type == "iRep" & !is.na(iRep_value))

df_irep$Genome_bin <- factor(df_irep$Genome_bin, levels = unique(df_irep$Genome_bin)[c(3,1,2)])

p_iRep_1 <- df_irep %>% 
  ggplot(aes(x = Reactor_cycle, y = iRep_value, shape = Manual_exception))+
  geom_point(size = 5, color = "black", aes(fill = Genome_bin),
             alpha = 0.7)+
  scale_shape_manual(values = c(21, 25))+
  theme_bw()+
  scale_fill_manual(values = c(col_RAMLI, col_bac1, col_bac2))+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=14),
        legend.background = element_rect(fill="transparent"),
        # axis.text.x = element_text(angle = 65, hjust = 1),
        strip.text.x=element_text(size = 18),
        legend.position=c(0.75,0.9),
        axis.text.x=element_text(size = 14, angle =0, hjust= 0.5),
        axis.title.x=element_blank(),
        plot.title = element_text(hjust = 0, size=18))+
  ylim(1.0,1.70)+
  labs(y = "iRep", scales = "free")+
  guides(shape = FALSE,
    fill  = guide_legend(title = "", override.aes = list(size = 5, shape = 21),
                         nrow = 4)
   )
  
print(p_iRep_1)
```

<img src="Figures/cached/irep-1-1.png" style="display: block; margin: auto;" />

```r
p_iRep_2 <- df_irep %>% 
  ggplot(aes(x = Genome_bin, y = iRep_value,fill = Genome_bin))+
    geom_boxplot(width = 0.3, alpha = 0.2)+
  geom_point(size = 5, color = "black", alpha = 0.7, aes(shape = Manual_exception))+
  scale_shape_manual(values = c(21, 25))+
  theme_bw()+
  scale_fill_manual(values = c(col_RAMLI, col_bac1, col_bac2))+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=14),
        legend.background = element_rect(fill="transparent"),
        # axis.text.x = element_text(angle = 65, hjust = 1),
        strip.text.x=element_text(size = 18),
        legend.position=c(0.75,0.9),
        axis.text.x=element_text(size = 14, angle =45, hjust= 1),
        axis.title.x=element_blank(),
        plot.title = element_text(hjust = 0, size=18))+
  ylim(1.0,1.70)+
  labs(y = "iRep", scales = "free")+
  guides(shape = FALSE,
    fill  = FALSE
   )
  
print(p_iRep_2)
```

<img src="Figures/cached/irep-1-2.png" style="display: block; margin: auto;" />

```r
# Print summary irep statistics
df_irep %>% group_by(Genome_bin) %>% summarize(mean(iRep_value))
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["Genome_bin"],"name":[1],"type":["fctr"],"align":["left"]},{"label":["mean(iRep_value)"],"name":[2],"type":["dbl"],"align":["right"]}],"data":[{"1":"Ramlibacter sp. MAG","2":"1.482997"},{"1":"Bacteroidetes sp. MAG1","2":"1.629398"},{"1":"Bacteroidetes sp. MAG2","2":"1.402130"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

```r
df_irep %>% group_by(Genome_bin) %>% summarize(sd(iRep_value))
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["Genome_bin"],"name":[1],"type":["fctr"],"align":["left"]},{"label":["sd(iRep_value)"],"name":[2],"type":["dbl"],"align":["right"]}],"data":[{"1":"Ramlibacter sp. MAG","2":"0.079643772"},{"1":"Bacteroidetes sp. MAG1","2":"0.037571020"},{"1":"Bacteroidetes sp. MAG2","2":"0.005680642"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

```r
# For figure
```


```r
p_iRep_3 <- df_irep %>% 
  ggplot(aes(x = X16S_growthrate_future, y = iRep_value, fill = Genome_bin))+
  geom_point(size = 5, color = "black", alpha = 0.7, aes(shape = Manual_exception))+
  scale_shape_manual(values = c(21, 25))+
  theme_bw()+
  scale_fill_manual(values = c(col_RAMLI, col_bac1, col_bac2))+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=12),
        legend.background = element_rect(fill="transparent"),
        # axis.text.x = element_text(angle = 65, hjust = 1),
        strip.text.x=element_text(size = 18),
        legend.position=c(0.7,0.1),
        axis.text.x=element_text(size = 14),
        plot.title = element_text(hjust = 0, size=18))+
  # ylim(1.2,1.8)+ 
  labs(y = "iRep", x = expression("Growth rate (cells µL"^-1*"d"^-1*")"))+
  guides(shape = FALSE,
    fill  = guide_legend(title = "", override.aes = list(size = 5, shape = 21),
                         nrow = 4)
   )+
  geom_smooth(method="lm", color="black", fill = "lightgray",
              alpha=0.2)


p_iRep_4 <- df_irep %>% 
  ggplot(aes(x = relative_abundance_16S, y = iRep_value, fill = Genome_bin))+
  geom_point(size = 5, color = "black", alpha = 0.7, aes(shape = Manual_exception))+
  scale_shape_manual(values = c(21, 25))+
  theme_bw()+
  scale_fill_manual(values = c(col_RAMLI, col_bac1, col_bac2))+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=12),
        legend.background = element_rect(fill="transparent"),
        # axis.text.x = element_text(angle = 65, hjust = 1),
        strip.text.x=element_text(size = 18),
        legend.position=c(0.7,0.1),
        axis.text.x=element_text(size = 14),
        plot.title = element_text(hjust = 0, size=18))+
  ylim(1.2,1.8)+ xlim(0,100)+
  labs(y = "iRep", x = "Relative abundance (%)")+
  guides(shape = FALSE,
    fill  = guide_legend(title = "", override.aes = list(size = 5, shape = 21),
                         nrow = 4)
   )+
  geom_smooth(method="lm", color="black", fill = "lightgray",
              alpha=0.2, formula = y ~ splines::ns(x,df=3))


p_iRep_5 <- df_irep %>% 
  ggplot(aes(x = absolute_abundance_16S, y = iRep_value, fill = Genome_bin))+
  geom_point(size = 5, color = "black", alpha = 0.7, aes(shape = Manual_exception))+
  scale_shape_manual(values = c(21, 25))+
  theme_bw()+
  scale_fill_manual(values = c(col_RAMLI, col_bac1, col_bac2))+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=12),
        legend.background = element_rect(fill="transparent"),
        # axis.text.x = element_text(angle = 65, hjust = 1),
        strip.text.x=element_text(size = 18),
        legend.position=c(0.7,0.1),
        axis.text.x=element_text(size = 14),
        plot.title = element_text(hjust = 0, size=18))+
  ylim(1.2,1.8)+ xlim(0,800)+
  labs(y = "iRep", x = expression("Absolute abundance (cells µL"^-1*")"))+
  guides(shape = FALSE,
    fill  = guide_legend(title = "", override.aes = list(size = 5, shape = 21),
                         nrow = 4)
   )+
  geom_smooth(method = "lm",color="black", fill ="lightgray", formula=y~x)
  
cowplot::plot_grid(p_iRep_4, p_iRep_5,
                   nrow = 1, align = "hv")
```

<img src="Figures/cached/irep-2-1.png" style="display: block; margin: auto;" />

```r
# ggpubr::ggarrange(p_iRep_3,                                                 #
#           ggpubr::ggarrange(p_iRep_4, p_iRep_5, ncol = 2), # Second row with box and dot plots
#           nrow = 2,                                       # Labels of the scatter plot
#           align = "h") 
```

# Trophic strategy markers  

Here we are searching the "core" and "accessory" genome for COGs/Pfams that were
shown to be important for Carbon-cycling in freshwater systems


```r
Tstrategy_df <- read.csv("./IMG_annotation/markers.Tstrategy.csv")

merged_gc_cog <- merged_gc_cog %>% dplyr::group_by(Genome) %>% 
  mutate(Noligo = sum(cog_id %in% Tstrategy_df$COG_id[Tstrategy_df$Trophic_strategy == "Oligotroph"]),
         Ncopio = sum(cog_id %in% Tstrategy_df$COG_id[Tstrategy_df$Trophic_strategy == "Copiotroph"]))

merged_gc_cog <- merged_gc_cog %>% mutate(Nratio = Noligo/Ncopio)

toMatch_oligo <- Tstrategy_df$COG_id[Tstrategy_df$Trophic_strategy == "Oligotroph"]
toMatch_copio <- Tstrategy_df$COG_id[Tstrategy_df$Trophic_strategy == "Copiotroph"]

test <- panG %>% dplyr::group_by(genome_name) %>% 
  mutate(Noligo = length(grep(paste(toMatch_oligo, collapse="|"), 
                        COG_FUNCTION_ACC, value=TRUE)),
         Ncopio = length(grep(paste(toMatch_copio, collapse="|"), 
                        COG_FUNCTION_ACC, value=TRUE))) %>% 
  mutate(Nratio = Noligo/Ncopio)



p_markers1 <- ggplot(test, aes(x = genome_name, y = Noligo, fill = genome_name))+
  theme_bw()+
  geom_bar(alpha = 0.4, stat = "identity", color = "black",
           position = position_dodge(width = 1), width = 0.7)+
  scale_fill_brewer(palette = "Paired")+
  theme(axis.text=element_text(size=13), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=14),
        legend.background = element_rect(fill="transparent"),
        # axis.text.x = element_text(angle = 65, hjust = 1),
        strip.text.x=element_text(size = 18),
        legend.position="bottom",
        axis.text.x=element_text(size = 13, angle =45, hjust= 1),
        axis.title.x=element_blank(),
        plot.title = element_text(hjust = 0, size=18))+
  ylab("")+
  xlab("")+
  ggtitle("Noligo")

# print(p_markers1)

p_markers2 <- ggplot(test, aes(x = genome_name, y = Ncopio, fill = genome_name))+
  theme_bw()+
  geom_bar(alpha = 0.4, stat = "identity", color = "black",
           position = position_dodge(width = 1), width = 0.7)+
  scale_fill_brewer(palette = "Paired")+
  theme(axis.text=element_text(size=13), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=14),
        legend.background = element_rect(fill="transparent"),
        # axis.text.x = element_text(angle = 65, hjust = 1),
        strip.text.x=element_text(size = 18),
        legend.position="bottom",
        axis.text.x=element_text(size = 13, angle =45, hjust= 1),
        axis.title.x=element_blank(),
        plot.title = element_text(hjust = 0, size=18))+
  ylab("")+
  xlab("")+
  ggtitle("Ncopio")

# print(p_markers2)

cowplot::plot_grid(p_markers1, p_markers2, nrow = 2, align = "v")
```

<img src="Figures/cached/markers-1-1.png" style="display: block; margin: auto;" />

# CAZymes  

### Exploration  


```r
# Import CAZY reference annotation
cazy_reference <- read.table("./IMG_annotation/CAZY_annotation.tsv", 
                             stringsAsFactors = FALSE, header = TRUE,
                            sep = "\t", fill = TRUE)

# Import CAZY annotation of IMG gene calls
cazy_files <- list.files(".", pattern = "*dbCAN*",
                       recursive = TRUE)
cazy_name <- c("Ramlibacter sp. MAG", "Bacteroidetes sp. MAG2",
               "Bacteroidetes sp. MAG1")
for(i in 1:length(cazy_files)){
  cazy_df <- read.table(cazy_files[i], stringsAsFactors = FALSE,
                         check.names =FALSE, header = TRUE)
  cazy_df$Genome <- cazy_name[i]
  if(i == 1) cazy_results_df <- cazy_df else cazy_results_df <- rbind(cazy_results_df, cazy_df)
}

# Summary of CAZY genes found in each genome
cazy_results_df %>% group_by(Genome) %>% summarise(CAZY_counts = n())
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["Genome"],"name":[1],"type":["chr"],"align":["left"]},{"label":["CAZY_counts"],"name":[2],"type":["int"],"align":["right"]}],"data":[{"1":"Bacteroidetes sp. MAG1","2":"293"},{"1":"Bacteroidetes sp. MAG2","2":"249"},{"1":"Ramlibacter sp. MAG","2":"108"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

```r
# Replace ".hmm" extension from reference
cazy_results_df$Subject <- gsub(".hmm", "", cazy_results_df$Subject)
cazy_results_df$Query <- as.character(cazy_results_df$Query)

# Combine with kegg annotation
# Identify genes that are under positive selection & that have DOM_usage annotation
merged_gc_cog_psg_cazy <- left_join(cazy_results_df, merged_gc_cog, by = c("Query"
                                                                      = "gene_oid"))
# Combine with CAZY reference information
merged_gc_cog_psg_cazy <- left_join(merged_gc_cog_psg_cazy, cazy_reference, 
                                    by = c("Subject"= "Family"))
merged_gc_cog_psg_cazy <- merged_gc_cog_psg_cazy %>% 
  mutate(PSG = Query %in% data_posi$Gene)

# Add AA annotation with was missing from the reference file
merged_gc_cog_psg_cazy$cazy_class[is.na(merged_gc_cog_psg_cazy$cazy_class)] <- "Auxiliary activities"

# Relabel other classes for interpretation
merged_gc_cog_psg_cazy$cazy_class <- plyr::revalue(merged_gc_cog_psg_cazy$cazy_class,
                                             c("CBM" = "Carbohydrate-binding modules",
                                               "CE" = "Carbohydrate esterases",
                                               "GH" = "Glycoside hydrolases",
                                               "GT" = "Glycosyltransferases",
                                               "PL" = "Polysaccharide lyases"))
# Customize factor levels
within_cazy_order <- merged_gc_cog_psg_cazy %>% 
  dplyr::filter(cazy_class != "Cellulosome") %>%
  group_by(cazy_class) %>% 
  summarise(cazy_class_freq = n()) %>% 
  arrange(desc(cazy_class_freq))

merged_gc_cog_psg_cazy$cazy_class <- factor(merged_gc_cog_psg_cazy$cazy_class,
                                            levels = within_cazy_order$cazy_class)

merged_gc_cog_psg_cazy$Genome.x <- factor(merged_gc_cog_psg_cazy$Genome.x ,
                                            levels = c("Ramlibacter sp. MAG",
                                                       "Bacteroidetes sp. MAG1",
                                                       "Bacteroidetes sp. MAG2"))

# Perform tSNE for the full CAZyme profiles
tsne_df_cazym <- merged_gc_cog_psg_cazy %>% 
  group_by(Genome.x, Subject) %>% 
  summarize(n_gene = n()) %>% 
  spread(., key = Genome.x, n_gene)

tsne_df_cazym[is.na(tsne_df_cazym)] <- 0

# perform tsne
set.seed(777)
dist_tsne <- dist(t(tsne_df_cazym[, -c(1)]))
df_tsne <- cmdscale(dist_tsne)

# From wide to long
# df_sigma <- gather(df_sigma, Genome, Counts, Curvibacter_sp._ATCC:Variovorax_sp._EPS)
df_tsne <- data.frame(Genome = colnames(tsne_df_cazym)[-c(1)], 
                       X = df_tsne[,1],
                       Y = df_tsne[,2])

p_ord <- ggplot(df_tsne, aes(x = X, y = Y, fill = Genome))+
  geom_point(size = 4, alpha = 0.7, shape = 21)+
  scale_fill_manual(values = c(col_bac1, col_bac2,"#c8c8ff"))+
  theme_bw()+
  labs(x = "Axis 1", y = "Axis 2")+
    theme(axis.text=element_text(size=13), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=14),
        legend.background = element_rect(fill="transparent"),
        # axis.text.x = element_text(angle = 65, hjust = 1),
        strip.text.x=element_text(size = 18),
        legend.position="bottom",
        plot.title = element_text(hjust = 0, size=18))

print(p_ord)
```

<img src="Figures/cached/CAZy-1-1.png" style="display: block; margin: auto;" />

### Plots

```r
# Make plot of cazy classes
p_cazy1 <- merged_gc_cog_psg_cazy %>% 
  dplyr::filter(cazy_class != "Cellulosome") %>%
  droplevels() %>% 
  ggplot(aes(x = cazy_class, fill = cazy_class))+
  geom_bar(stat="count", color = "black", width = 0.5)+
  theme_bw()+
  scale_fill_brewer(palette = "Accent")+
  ylab("Number of genes") + xlab("")+
  theme(axis.text.y=element_text(size=15), axis.title=element_text(size = 20),
        legend.text=element_text(size=18),
        legend.background = element_rect(fill="transparent"),
        strip.text=element_text(size=14), legend.title = element_blank()
        ,legend.position = c(0.87, 0.9)
        )+
  scale_y_continuous(breaks = seq(0,100,20), limits = c(0,80))+
  facet_grid(.~Genome.x)+
  coord_flip()+
  theme(axis.text.x = element_text(size=16, angle = 0))+
  # theme(axis.line = element_line(size = 1, colour = "grey80"),
  #       panel.border = element_blank())+
guides(fill = FALSE)

print(p_cazy1)
```

<img src="Figures/cached/CAZy-2-1.png" style="display: block; margin: auto;" />


```r
# Make plot of GH alone
# Customize factor levels
within_GH_order <- merged_gc_cog_psg_cazy %>% 
  dplyr::filter(cazy_class == "Glycoside hydrolases") %>%
  group_by(Subject) %>% 
  summarise(Subject_Freq = n()) %>% 
  arrange(desc(Subject_Freq))

merged_gc_cog_psg_cazy_GH <- merged_gc_cog_psg_cazy %>% 
  dplyr::filter(cazy_class == "Glycoside hydrolases") %>% 
  droplevels()
merged_gc_cog_psg_cazy_GH$Subject <- factor(merged_gc_cog_psg_cazy_GH$Subject,
                                            levels = within_GH_order$Subject)
p_cazy2 <- merged_gc_cog_psg_cazy_GH %>% 
  dplyr::filter(cazy_class == "Glycoside hydrolases") %>%
  droplevels() %>% 
  ggplot(aes(x = Subject, fill = Genome.x))+
  geom_bar(stat="count", color = "black")+
  theme_bw()+
  scale_fill_manual(values=c(col_RAMLI, col_bac1, col_bac2))+
  ylab("Number of genes") + xlab("Glycoside hydrolase family")+
  theme(axis.text.y=element_text(size=13), axis.title=element_text(size = 20),
        legend.text=element_text(size=18),
        legend.background = element_rect(fill="transparent"),
        strip.text=element_text(size=14), legend.title = element_blank()
        ,legend.position = c(0.87, 0.9)
        )+
  scale_y_continuous(breaks = seq(0,20,5))+
  facet_grid(.~Genome.x)+
  coord_flip()+
  theme(axis.text.x = element_text(size=16, angle = 0))+
  # theme(axis.line = element_line(size = 1, colour = "grey80"),
  #       panel.border = element_blank())+
guides(fill = FALSE)

print(p_cazy2)
```

<img src="Figures/cached/CAZy-3-1.png" style="display: block; margin: auto;" />

```r
# Make plot of GT alone
# Customize factor levels
within_GT_order <- merged_gc_cog_psg_cazy %>% 
  dplyr::filter(cazy_class == "Glycosyltransferases") %>%
  group_by(Subject) %>% 
  summarise(Subject_Freq = n()) %>% 
  arrange(desc(Subject_Freq))

merged_gc_cog_psg_cazy_GT <-  merged_gc_cog_psg_cazy %>% 
  dplyr::filter(cazy_class == "Glycosyltransferases") %>% 
  droplevels()
merged_gc_cog_psg_cazy_GT$Subject <- factor(merged_gc_cog_psg_cazy_GT$Subject,
                                            levels = within_GT_order$Subject)
p_cazy3 <- merged_gc_cog_psg_cazy_GT %>% 
  dplyr::filter(cazy_class == "Glycosyltransferases") %>%
  droplevels() %>% 
  ggplot(aes(x = Subject, fill = Genome.x))+
  geom_bar(stat="count", color = "black")+
  theme_bw()+
  scale_fill_manual(values=c(col_RAMLI, col_bac1, col_bac2))+
  ylab("Number of genes") + xlab("Glycosyltransferases")+
  theme(axis.text.y=element_text(size=13), axis.title=element_text(size = 20),
        legend.text=element_text(size=18),
        legend.background = element_rect(fill="transparent"),
        strip.text=element_text(size=14), legend.title = element_blank()
        ,legend.position = c(0.87, 0.9)
        )+
  scale_y_continuous(breaks = seq(0,20,5))+
  facet_grid(.~Genome.x)+
  coord_flip()+
  theme(axis.text.x = element_text(size=16, angle = 0))+
  # theme(axis.line = element_line(size = 1, colour = "grey80"),
  #       panel.border = element_blank())+
guides(fill = FALSE)

print(p_cazy3)
```

<img src="Figures/cached/CAZy-3-2.png" style="display: block; margin: auto;" />

### Heatmap  


```r
# make heatmap of GH
hm_GH <- merged_gc_cog_psg_cazy_GH %>%
  group_by(Subject, Genome.x) %>% 
  dplyr::summarise(Counts = n()) %>% 
  complete(., Genome.x, Subject, fill = list(Counts=0)) %>% 
  ggplot(aes(y = Genome.x, x = Subject)) + # x and y axes => Var1 and Var2
  geom_tile(aes(fill = Counts), col = "lightgrey") + # background colours are mapped according to the value column
  geom_text(aes(label = round(Counts, 0)), size = 3) + # write the values
  # scale_fill_gradientn(colours = terrain.colors(10), trans = "log1p")+
  # scale_fill_gradient(low = "lightblue", high = "darkslategray", na.value="white",
                      # trans = "log1p", limits=c(1, 40)) +
  scale_fill_distiller(palette="YlOrRd", na.value="lightgrey", trans = "sqrt",
                       direction = 1) +
  theme(panel.grid.major.x=element_blank(), #no gridlines
        panel.grid.minor.x=element_blank(), 
        panel.grid.major.y=element_blank(), 
        panel.grid.minor.y=element_blank(),
        panel.background=element_rect(fill="white"), # background=white
        axis.text.x = element_text(angle=45, hjust = 1, vjust=1, size = 12,face = "bold"),
        plot.title = element_text(size=20,face="bold"),
        axis.text.y = element_text(size = 12,face = "bold"))+
  theme(legend.title=element_text(face="bold", size=14)) + 
  scale_x_discrete(name="") +
  scale_y_discrete(name="") +
  labs(fill="Gene\ncount")

print(hm_GH)
```

<img src="Figures/cached/Cazy-4-1.png" style="display: block; margin: auto;" />

# DOC-transporters  

### Exploration


```r
# Import DOM usage table
DOM_usage_df <- read.csv("./IMG_annotation/DOM_usage.csv", stringsAsFactors = FALSE,
                         check.names =FALSE)
DOM_usage_df$DOM_type <- gsub("\\n","\n", DOM_usage_df$DOM_type, fixed = TRUE)

# Import IMG exported COG-annotation of all genomes in phylogenomic tree
COG_profiles <- read.table("./IMG_annotation/STAMP_profiles/abundance_cog_37895.tab.xls", header = TRUE, sep = "\t", quote = "")

# Retain COG_ids that are found in DOM_usage_df list
COG_profiles_sub <- COG_profiles %>% dplyr::filter(Func_id %in% DOM_usage_df$COG_ID)

COG_profiles_sub <- dplyr::left_join(COG_profiles_sub, DOM_usage_df, 
                                     by = c("Func_id" = "COG_ID")) %>% 
  dplyr::select(contains("Func"),"DOM_type", contains("Ramli"),contains("Bacteroidetes"))
rownames(COG_profiles_sub) <- COG_profiles_sub$Func_id

# Remove zero rows
COG_profiles_sub <- COG_profiles_sub[rowSums(COG_profiles_sub[, - c(1:3)])>0, ]

# Wide to long format
COG_profiles_sub_long <- tidyr::gather(COG_profiles_sub, Genome, 
                                       Counts, 
                          Ramlibacter_sp._MAG:Bacteroidetes_sp._MAG1,
                          factor_key = TRUE)
COG_profiles_sub_long$Genome <- gsub("_", " ", COG_profiles_sub_long$Genome)
COG_profiles_sub_long$Genome <- gsub("5.10", "5-10", COG_profiles_sub_long$Genome)
```



```r
# Identify genes that are under positive selection & that have DOM_usage annotation
merged_gc_cog_psg <- merged_gc_cog %>% 
  dplyr::filter(Genome == "121950.assembled.gff") %>% 
  dplyr::mutate(PSG = gene_oid %in% data_posi$Gene) 

# Label DOM-usage genes
merged_gc_cog_psg <- merged_gc_cog_psg %>% 
  dplyr::left_join(., DOM_usage_df, by = c("cog_id" = "COG_ID"))

# Overview of DOM-genes
merged_gc_cog_psg_POS <- 
  merged_gc_cog_psg %>% dplyr::filter(PSG  == TRUE & !is.na(DOM_type))

# Order according to COG counts
merged_gc_cog_psg_POS$cog_id <- factor(merged_gc_cog_psg_POS$cog_id,
                                       levels = names(table(merged_gc_cog_psg_POS$cog_id))[rev(order(table(merged_gc_cog_psg_POS$cog_id)))])


p_DOM_psg <- ggplot2::ggplot(merged_gc_cog_psg_POS, aes(x = DOM_type, fill = cog_name))+
  geom_bar(alpha = 0.4, stat = "count", color = "black", width = 0.5)+
  scale_fill_brewer("", palette = "Paired")+
  theme_bw()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=10),
        legend.background = element_rect(fill="transparent"),
        # axis.text.x = element_text(angle = 65, hjust = 1),
        strip.text.x=element_text(size = 18),
        legend.position="top",
        axis.text.x=element_text(size = 14, angle =45, hjust= 1),
        axis.title.x=element_blank(),
        plot.title = element_text(hjust = 0, size=18))+
  guides(fill=guide_legend(ncol=1))+
  ylab("Number of PSGs")

print(p_DOM_psg)
```

<img src="Figures/cached/Transp-2-1.png" style="display: block; margin: auto;" />

```r
# Negative for PSG DOM-uptake genes
merged_gc_cog_psg_NEG <-
  merged_gc_cog_psg %>% dplyr::filter(PSG  == FALSE & !is.na(DOM_type))

# Order according to COG counts
merged_gc_cog_psg_NEG$cog_id <- factor(merged_gc_cog_psg_NEG$cog_id,
                                       levels = names(table(merged_gc_cog_psg_NEG$cog_id))[rev(order(table(merged_gc_cog_psg_NEG$cog_id)))])

# p_DOM_no_psg <- ggplot2::ggplot(merged_gc_cog_psg_NEG, aes(x = cog_id, fill = DOM_type))+
#   geom_bar(alpha = 0.4, stat = "count", color = "black",
#            position = position_dodge(width = 1), width = 0.7)+
#   scale_fill_brewer("", palette = "Paired")+
#   theme_bw()+
#   theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
#         title=element_text(size=20), legend.text=element_text(size=12),
#         legend.background = element_rect(fill="transparent"),
#         # axis.text.x = element_text(angle = 65, hjust = 1),
#         strip.text.x=element_text(size = 18),
#         legend.position="top",
#         axis.text.x=element_text(size = 14, angle =45, hjust= 1),
#         axis.title.x=element_blank(),
#         plot.title = element_text(hjust = 0, size=18))+
#   guides(fill=guide_legend(ncol=1))+
#   ylab("Number of PSGs")
# 
# print(p_DOM_no_psg)
```

### Heatmap

```r
# Select subset of genomes for phosphate scavenging
selected_genomes <- c("Ramlibacter sp. TTB310",
                      "Ramlibacter sp. Leaf400",
                      "Ramlibacter sp. 5-10",
                      "Ramlibacter sp. MAG",
                      "Bacteroidetes sp. MAG1", 
                      "Bacteroidetes sp. MAG2")
# Order fators
COG_profiles_sub_long$Genome <- factor(as.character(COG_profiles_sub_long$Genome),
                                 levels = selected_genomes)

COG_order <- DOM_usage_df %>% 
  dplyr::filter(COG_ID %in% unique(COG_profiles_sub_long$Func_id))
COG_profiles_sub_long$Func_id <- factor(COG_profiles_sub_long$Func_id,
                                        levels = COG_order$COG_ID)
# make heatmap
hm_DOC1 <- COG_profiles_sub_long %>% 
  # dplyr::filter(Genome %in% selected_genomes) %>% 
  ggplot(aes(y= Genome, x= Func_id)) + # x and y axes => Var1 and Var2
  geom_tile(aes(fill = Counts), col = "lightgrey") + # background colours are mapped according to the value column
  geom_text(aes(label = round(Counts, 0)), size = 3) + # write the values
  # scale_fill_gradientn(colours = terrain.colors(10), trans = "log1p")+
  # scale_fill_gradient(low = "lightblue", high = "darkslategray", na.value="white",
                      # trans = "log1p", limits=c(1, 40)) +
  scale_fill_distiller(palette="YlOrRd", na.value="lightgrey", trans = "sqrt",
                       direction = 1) +
  theme(panel.grid.major.x=element_blank(), #no gridlines
        panel.grid.minor.x=element_blank(), 
        panel.grid.major.y=element_blank(), 
        panel.grid.minor.y=element_blank(),
        panel.background=element_rect(fill="white"), # background=white
        axis.text.x = element_text(angle=45, hjust = 1, vjust=1, size = 12,face = "bold"),
        plot.title = element_text(size=20,face="bold"),
        axis.text.y = element_text(size = 12,face = "bold"))+
  theme(legend.title=element_text(face="bold", size=14)) + 
  scale_x_discrete(name="") +
  scale_y_discrete(name="") +
  labs(fill="Gene\ncount")

print(hm_DOC1)
```

<img src="Figures/cached/Transp-3-1.png" style="display: block; margin: auto;" />


# Phosphate scavenging genes  

### All genomes 



```r
df_pho <- read.csv("./IMG_annotation/custom_tables/Annotation_P.csv", header = TRUE,
                        stringsAsFactors = FALSE)
# df_pho[, -1] <- scale(df_pho[, -1])
df_pho_long <- df_pho %>% melt()
```

```
## Using Genome as id variables
```

```r
df_pho_long$variable <- gsub("_", " - ", df_pho_long$variable)
colnames(df_pho) <- gsub("_", " - ", colnames(df_pho))


# Order fators
df_pho_long$Genome <- factor(as.character(df_pho_long$Genome),
                                 levels = df_pho$Genome)

df_pho_long$variable <- factor(as.character(df_pho_long$variable),
                                 levels = colnames(df_pho)[-1])

# make heatmap
hm_pho1 <- df_pho_long %>% 
  ggplot(aes(y= Genome, x= variable)) + # x and y axes => Var1 and Var2
  geom_tile(aes(fill = value), col = "lightgrey") + # background colours are mapped according to the value column
  geom_text(aes(label = round(df_pho_long$value, 0))) + # write the values
  # scale_fill_gradientn(colours = terrain.colors(10), trans = "log1p")+
  # scale_fill_gradient(low = "lightblue", high = "darkslategray", na.value="white",
                      # trans = "log1p", limits=c(1, 40)) +
  scale_fill_distiller(palette="YlOrRd", na.value="lightgrey", trans = "sqrt",
                       direction = 1) +
  theme(panel.grid.major.x=element_blank(), #no gridlines
        panel.grid.minor.x=element_blank(), 
        panel.grid.major.y=element_blank(), 
        panel.grid.minor.y=element_blank(),
        panel.background=element_rect(fill="white"), # background=white
        axis.text.x = element_text(angle=45, hjust = 1, vjust=1,size = 12,face = "bold"),
        plot.title = element_text(size=20,face="bold"),
        axis.text.y = element_text(size = 12,face = "bold"))+
  theme(legend.title=element_text(face="bold", size=14)) + 
  scale_x_discrete(name="") +
  scale_y_discrete(name="") +
  labs(fill="Gene\ncount")

print(hm_pho1)
```

<img src="Figures/cached/Pho-1-1.png" style="display: block; margin: auto;" />

```r
# changed again
```

### Subset of genomes  


```r
df_pho_long_sel <- df_pho %>% melt()
```

```
## Using Genome as id variables
```

```r
df_pho_long_sel$variable <- gsub("_", " - ", df_pho_long_sel$variable)

# Select subset of genomes for phosphate scavenging
selected_genomes <- c("Ramlibacter sp. TTB310",
                      "Ramlibacter sp. Leaf400",
                      "Ramlibacter sp. 5-10",
                      "Ramlibacter sp. MAG",
                      "Unclassified Bacteroidetes sp. MAG1", 
                      "Unclassified Bacteroidetes sp. MAG2")
df_pho_long_sel <- df_pho_long %>% 
  dplyr::filter(Genome %in% selected_genomes)

# Order fators
df_pho_long_sel$Genome <- factor(as.character(df_pho_long_sel$Genome),
                                 levels = selected_genomes)

df_pho_long_sel$variable <- factor(as.character(df_pho_long_sel$variable),
                                 levels = colnames(df_pho)[-1])

# make heatmap
hm_pho2 <- df_pho_long_sel %>% 
  dplyr::filter(Genome %in% selected_genomes) %>% 
  ggplot(aes(y= Genome, x= variable)) + # x and y axes => Var1 and Var2
  geom_tile(aes(fill = value), col = "lightgrey") +
  geom_text(aes(label = round(value, 0))) + # write the values
  # scale_fill_gradientn(colours = terrain.colors(10), trans = "log1p")+
  # scale_fill_gradient(low = "lightblue", high = "darkslategray", na.value="white",
                      # trans = "log1p", limits=c(1, 40)) +
  scale_fill_distiller(palette="YlOrRd", na.value="lightgrey", trans = "sqrt",
                       direction = 1) +
  theme(panel.grid.major.x=element_blank(), #no gridlines
        panel.grid.minor.x=element_blank(), 
        panel.grid.major.y=element_blank(), 
        panel.grid.minor.y=element_blank(),
        panel.background=element_rect(fill="white"), # background=white
        axis.text.x = element_text(angle=45, hjust = 1,vjust=1,size = 12,face = "bold"),
        plot.title = element_text(size=20,face="bold"),
        axis.text.y = element_text(size = 12,face = "bold"))+
  theme(legend.title=element_text(face="bold", size=14)) + 
  scale_x_discrete(name="") +
  scale_y_discrete(name="") +
  labs(fill="Gene\ncount")

print(hm_pho2)
```

<img src="Figures/cached/Pho-2-1.png" style="display: block; margin: auto;" />

```r
# changed again
```


### Check for PSGss

```r
# Get list of COGs/KO
p_cog_ko_list <- do.call(rbind, base::strsplit(split = " - ", colnames(df_pho)[-1]))[,2]

# Get P COGs that are PSG
print(data_posi_COG %>% dplyr::filter(cog_id %in% p_cog_ko_list) %>% distinct())
```

```
##         Gene   P.Value         FDR HA.foreground.omega HA.kappa
## 1 2727790921 0.0006156 0.002841803            1.585493  2.07378
## 2 2727790928 0.0042980 0.014352524           12.180261  1.94745
## 3 2727790930 0.0107200 0.031259847            1.075120  1.75049
##   Number.of.Sites.under.positive.Selection       IMG_geneID         contig
## 1                                       11 Ga0182879_101228 Ga0182879_1012
## 2                                       11 Ga0182879_101235 Ga0182879_1012
## 3                                        2 Ga0182879_101237 Ga0182879_1012
##      GC               Genome   gene_oid gene_length percent_identity
## 1 0.702 121950.assembled.gff 2727790921         697            48.76
## 2 0.689 121950.assembled.gff 2727790928         299            47.20
## 3 0.655 121950.assembled.gff 2727790930         234            37.23
##   query_start query_end subj_start subj_end evalue bit_score  cog_id
## 1           6       685          2      683  0e+00       741 COG0855
## 2          21       299          4      289  4e-67       248 COG0581
## 3           4       234          7      235  3e-36       145 COG0704
##                                                  cog_name cog_length
## 1                                    Polyphosphate kinase        696
## 2 ABC-type phosphate transport system, permease component        292
## 3                              Phosphate uptake regulator        240
##             genome_id COG_class                COG_functional_category
## 1 Ramlibacter sp. MAG         P Inorganic ion transport and metabolism
## 2 Ramlibacter sp. MAG         P Inorganic ion transport and metabolism
## 3 Ramlibacter sp. MAG         P Inorganic ion transport and metabolism
##   COG_functional_cluster
## 1             METABOLISM
## 2             METABOLISM
## 3             METABOLISM
```

```r
# Get P KOs that are PSG
print(data_posi_KO %>% dplyr::filter(ko_id %in% p_cog_ko_list))
```

```
##  [1] Gene                                    
##  [2] P.Value                                 
##  [3] FDR                                     
##  [4] HA.foreground.omega                     
##  [5] HA.kappa                                
##  [6] Number.of.Sites.under.positive.Selection
##  [7] IMG_geneID                              
##  [8] contig                                  
##  [9] GC                                      
## [10] Genome                                  
## [11] gene_oid                                
## [12] gene_length                             
## [13] percent_identity                        
## [14] query_start                             
## [15] query_end                               
## [16] subj_start                              
## [17] subj_end                                
## [18] evalue                                  
## [19] bit_score                               
## [20] ko_id                                   
## [21] ko_name                                 
## [22] img_ko_flag                             
## [23] genome_id                               
## [24] ko_level_A                              
## [25] ko_level_B                              
## [26] ko_level_C                              
## <0 rows> (or 0-length row.names)
```

# Sigma factors  


```r
# We use KEGG annotation because of its more defined sigma factor annotations
df_sigma <- read.table("./IMG_annotation/STAMP_profiles/abundance_ko_38065.tab.xls",
                       header = TRUE, sep = "\t", quote = "") %>% 
  dplyr::filter(grepl("sigma", Func_name))

vegan::specnumber(df_sigma[, 3: ncol(df_sigma)])
```

```
##  [1] 10 15 10 15 13 15 10 10 15  1  3 12  4
```

```r
# From wide to long
# df_sigma <- gather(df_sigma, Genome, Counts, Curvibacter_sp._ATCC:Variovorax_sp._EPS)
df_sigma <- data.frame(Genome = colnames(df_sigma)[-c(1,2)], 
                       Counts = colSums(df_sigma[,-c(1,2)]))

# Order genome_ids according to the phylogenetic clustering
df_sigma$Genome <- gsub("_", " ", df_sigma$Genome)
ord_full_list_bin_sigma <- c("Limnohabitans sp. Rim47", "Limnohabitans sp. Rim28",
                  "Limnohabitans sp. 63ED37", 
                  "Rhodoferax sp. ED16","Rhodoferax sp. T118",
                  "Curvibacter sp. ATCC", "Curvibacter sp. PAE",
                  "Variovorax sp. 110B", "Variovorax sp. EPS",
                  "Ramlibacter sp. Leaf400", "Ramlibacter sp. TTB310",
                  "Ramlibacter sp. MAG", 
                  "Ramlibacter sp. 5.10",
                  "Bacteroidetes sp. MAG1",
                  "Bacteroidetes sp. MAG2"
                  )
df_sigma$Genome <- factor(df_sigma$Genome, levels = ord_full_list_bin_sigma)

# Make plot
p_sigma_1 <- df_sigma %>% dplyr::filter(Genome != "Bacteroidetes sp. MAG1" &
                                        Genome != "Bacteroidetes sp. MAG2") %>% 
  ggplot(aes(x = Genome, y = Counts, fill = Genome))+
  theme_bw()+
  geom_bar(alpha = 0.4, stat = "identity", color = "black",
           position = position_dodge(width = 1), width = 0.7)+
  scale_fill_manual(values = c(rep(adjustcolor("#c8c8ff",0.8),3), rep("#f8cf94",2), 
                               rep("#adf7ad",2), rep(adjustcolor("#000000",0.21),2),
                               rep("#e2a2fd",4), col_bac1, col_bac2))+
  theme(axis.text=element_text(size=13), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=14),
        legend.background = element_rect(fill="transparent"),
        # axis.text.x = element_text(angle = 65, hjust = 1),
        strip.text.x=element_text(size = 18),
        legend.position="bottom",
        axis.text.x=element_text(size = 9, angle =55, hjust= 1),
        axis.title.x=element_blank(),
        plot.title = element_text(hjust = 0, size=18))+
  guides(fill=FALSE)+
  ylab("")+
  xlab("")+
  ggtitle(expression("Number of"~sigma~"factor homologs"))+
  ylim(c(0,30))

p_sigma_2 <- df_sigma %>% dplyr::filter(Genome == "Bacteroidetes sp. MAG1" |
                                        Genome == "Bacteroidetes sp. MAG2") %>% 
  ggplot(aes(x = Genome, y = Counts, fill = Genome))+
  theme_bw()+
  geom_bar(alpha = 0.4, stat = "identity", color = "black",
           position = position_dodge(width = 1), width = 0.7)+
  scale_fill_manual(values = c(col_bac1, col_bac2))+
  theme(axis.text=element_text(size=13), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=14),
        legend.background = element_rect(fill="transparent"),
        # axis.text.x = element_text(angle = 65, hjust = 1),
        strip.text.x=element_text(size = 18),
        legend.position="bottom",
        axis.text.x=element_text(size = 9, angle =55, hjust= 1),
        axis.title.x=element_blank(),
        plot.title = element_text(hjust = 0, size=18))+
  guides(fill=FALSE)+
  ylab("")+
  xlab("")+
  ggtitle("")+
  ylim(c(0,30))

cowplot::plot_grid(p_sigma_1, p_sigma_2, align = "hv", rel_widths = c(4,1))
```

<img src="Figures/cached/sigma-1-1.png" style="display: block; margin: auto;" />

# Ordination


```r
# We use KEGG annotation because of its more defined sigma factor annotations
df_ord <- read.table("./IMG_annotation/STAMP_profiles/abundance_ko_38065.tab.xls",
                       header = TRUE, sep = "\t", quote = "")
# perform tsne
set.seed(777)
dist_tsne <- dist(t(df_ord[, -c(1,2)]))
df_tsne <- Rtsne::Rtsne(dist_tsne, perplexity = 1)

# From wide to long
# df_sigma <- gather(df_sigma, Genome, Counts, Curvibacter_sp._ATCC:Variovorax_sp._EPS)
df_tsne <- data.frame(Genome = colnames(df_ord)[-c(1,2)], 
                       X = df_tsne$Y[,1],
                       Y = df_tsne$Y[,2])

# Order genome_ids according to the phylogenetic clustering
df_tsne$Genome <- gsub("_", " ", df_tsne$Genome)
df_tsne$Genome <- factor(df_tsne$Genome, levels = ord_full_list_bin_sigma)

p_ord <- ggplot(df_tsne, aes(x = X, y = Y, fill = Genome))+
  geom_point(size = 4, alpha = 0.7, shape = 21)+
  scale_fill_manual(values = c(rep(adjustcolor("#c8c8ff",0.8),3), rep("#f8cf94",2), 
                               rep("#adf7ad",2), rep(adjustcolor("#000000",0.21),2),
                               rep("#e2a2fd",4), col_bac1, col_bac2))+
  theme_bw()+
  labs(x = "Axis 1", y = "Axis 2")+
    theme(axis.text=element_text(size=13), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=14),
        legend.background = element_rect(fill="transparent"),
        # axis.text.x = element_text(angle = 65, hjust = 1),
        strip.text.x=element_text(size = 18),
        legend.position="bottom",
        plot.title = element_text(hjust = 0, size=18))

print(p_ord)
```

<img src="Figures/cached/ord-1-1.png" style="display: block; margin: auto;" />

# Size distribution R. aquaticus


```r
size_files <- list.files(".", pattern = "RP_",
                       recursive = TRUE)
Raquat_size <- data.frame()
for(size_file in size_files){
  tmp <- read.delim(size_file, fill = TRUE)
  tmp <- data.frame(tmp, sample = do.call(rbind, strsplit(size_file, "_"))[,2])
  if(size_file == size_files[1]) Raquat_size <- tmp else{
    Raquat_size <- rbind(Raquat_size, tmp)
  }
}

# Make ggplot plot

p_size <- ggplot(Raquat_size, aes(x= "R. aquaticus", y = Length))+
  # geom_jitter(shape = 21, width = 0.025, size = 2.5, fill = col_RAMLI, alpha = 0.5)+
  geom_violin(width = 0.5, alpha = 0.2, fill = col_RAMLI, bw = 0.1, 
              adjust = 1, draw_quantiles = TRUE,
              trim = TRUE, scale = "count")+
  stat_summary(fun.data=mean_sdl, fun.args = list(mult = 1), 
                 geom="pointrange", color="black", size = 1.5, alpha = 0.75)+
  xlab("")+ ylab("")+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        plot.title = element_text(size = 20, hjust = 0.5))+
  scale_y_sqrt(breaks = c(0.25,0.5,1,2,4,6,8), labels = c(0.25,0.5,1,2,4,6,8), 
               lim = c(0.25,8.1),
               position = "right")+
  labs(title = "Length (µm)")

# Calculate mean & st dev
Raquat_size %>% dplyr::filter(Length<=1) %>% summarize(mean(Length)) %>% round(.,1)
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["mean(Length)"],"name":[1],"type":["dbl"],"align":["right"]}],"data":[{"1":"0.7"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

```r
Raquat_size %>% dplyr::filter(Length>1) %>% summarize(mean(Length))%>% round(.,1)
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["mean(Length)"],"name":[1],"type":["dbl"],"align":["right"]}],"data":[{"1":"2.5"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

```r
Raquat_size %>% dplyr::filter(Length<=1) %>% summarize(sd(Length))%>% round(.,1)
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["sd(Length)"],"name":[1],"type":["dbl"],"align":["right"]}],"data":[{"1":"0.2"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

```r
Raquat_size %>% dplyr::filter(Length>1) %>% summarize(sd(Length))%>% round(.,1)
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["sd(Length)"],"name":[1],"type":["dbl"],"align":["right"]}],"data":[{"1":"1.6"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

```r
# Partition around medoids clustering
Raquat_size$cluster_alloc <- factor(cluster::pam(Raquat_size$Length, 3)$cluster)
tmp.si <- c()
for(i in 2:10){
    tmp.si[i] <- cluster::pam(Raquat_size$Length, k = i)$silinfo$avg.width
}
nr_clusters_bacteria <- which(tmp.si == max(tmp.si, na.rm = TRUE))
Raquat_size$cluster_alloc <- factor(cluster::pam(Raquat_size$Length,
                                        nr_clusters_bacteria)$cluster)

p_size2 <- ggplot(Raquat_size, aes(x= "R. aquaticus", y = Length))+
  geom_jitter(shape = 21, width = 0.1, size = 4, aes(fill = cluster_alloc)
              , alpha = 0.7)+
  scale_fill_brewer(palette = "Paired")+
  xlab("")+ 
  ylab("")+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        plot.title = element_text(size = 20, hjust = 0.5))+
  scale_y_sqrt(breaks = c(0.25,0.5,1,2,4,6,8), labels = c(0.25,0.5,1,2,4,6,8), 
               lim = c(0.25,8.1),
               position = "right")+
  labs(title = "Length (µm)")+
  guides(fill = FALSE)

cowplot::plot_grid(p_size, p_size2, ncol =2, align = "h")
```

<img src="Figures/cached/size-1-1.png" style="display: block; margin: auto;" />


# External metagenomes


```r
df_BlastExt <- data.table::fread("./ExternalData/merged_blastfiles.tsv", 
                         header = TRUE)[, -4]
# Import metadata
meta_Ext <- xlsx::read.xlsx("./ExternalData/ismej2017156x1.xlsx", sheetIndex = 1)
meta_Ext <- meta_Ext[!is.na(meta_Ext$metagenomic.sample),]
meta_Ext <- meta_Ext[, c(1, 2, 4, 5)]

# Import contig_ids of each bin
contig_ids <- rbind(
  cbind(read.table("./IMG_annotation/IMG_2724679690_Ramlibacter_bin/121950.assembled.names_map", stringsAsFactors = FALSE)[,1], "Ramlibacter sp. MAG"),
  cbind(read.table("./IMG_annotation/IMG_2724679691_Bacteroidetes_bin1/121951.assembled.names_map", stringsAsFactors = FALSE)[,1], "Bacteroidetes sp. MAG1"),
  cbind(read.table("./IMG_annotation/IMG_2724679698_Bacteroidetes_bin2/121960.assembled.names_map", stringsAsFactors = FALSE)[,1], "Bacteroidetes sp. MAG2")
)
contig_ids <- data.frame(contig_ids)
colnames(contig_ids) <- c("contig_id", "bin")

# Clean up SRA identifier label
df_BlastExt$SRA <- gsub("_blast.tsv", "", df_BlastExt$SRA)

# Remove contigs not in contig_ids
df_BlastExt <- df_BlastExt %>% dplyr::filter(contig_id %in% contig_ids$contig_id)

# Add contig labels of MAGs
df_BlastExt <- dplyr::left_join(df_BlastExt, contig_ids, by = c("contig_id"))

# Bin the %Identity in intervals of 0.5%
df_BlastExt_sum <- transform(df_BlastExt, bin_group = cut(identity,  breaks=seq(0, 100, 0.5)))

# Add extra column that converts the binning range to a numeric x-coordinate that
# is positioned in the middle of the binning interval
df_BlastExt_sum$bin_group <- gsub("\\(|]", "", 
                               as.character(df_BlastExt_sum$bin_group))
df_BlastExt_sum$bin_xcoord <- as.numeric(do.call(rbind,
                                              strsplit(df_BlastExt_sum$bin_group,
                                                       ","))[,1])+0.25

# Add metadata
df_BlastExt_sum <- left_join(df_BlastExt_sum, meta_Ext, by = c("SRA" = "accession.number"))

SRA_ranking <- df_BlastExt_sum %>% dplyr::filter(bin == "Ramlibacter sp. MAG" & identity > 95) %>%
  dplyr::select(SRA) %>% table()
SRA_ranking <- SRA_ranking[rev(order(SRA_ranking))]

# Plot sequence discrete populations
p_blast_sdisc <- df_BlastExt_sum %>% 
     ggplot(aes(x = bin_xcoord, ..scaled.., fill = bin))+
      theme_bw()+
      scale_fill_manual("", values = c(col_bac1,col_bac2,col_RAMLI))+
      facet_grid(bin~Water_type)+
      geom_density(color = "black")+
      guides(color = FALSE)+
      theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=14),
        legend.background = element_rect(fill="transparent"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text=element_text(size=14))+
      ylab("Density")+
      xlab("Nucleotide identity (%)")+
    xlim(75,100)

print(p_blast_sdisc)
```

<img src="Figures/cached/SDP_external-1.png" style="display: block; margin: auto;" />

```r
# Plot for all bins density plots
p_blast_all_dens <- df_BlastExt_sum %>% 
  ggplot(aes(x = bin_xcoord, shape = SRA))+
  theme_bw()+
  geom_density(alpha = 0.4, size = 0.4, color = "#333333",
               bw = "nrd0")+
  scale_color_brewer("", palette = "Paired")+
  guides(fill = FALSE)+
  facet_grid(Water_type~bin)+
  theme(axis.text.y=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=14),
        legend.background = element_rect(fill="transparent"),
        axis.text.x = element_text(size = 14),
        strip.text=element_text(size=14),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")+
  ylab("")+
  xlab("% Identity")+
  scale_x_continuous(limits = c(80,100))+
  guides(shape = FALSE)

p_blast_all_dens
```

<img src="Figures/cached/SDP_external-2.png" style="display: block; margin: auto;" />


```r
# Plot % reads corrected for genome size over threshold of 0.95
blast_df_sum_comp <- df_BlastExt_sum %>% group_by(SRA, bin) %>% dplyr::count(bin_xcoord)

id_thresh <- 95-0.25
map_disc_cum <- blast_df_sum_comp  %>% 
  dplyr::filter(bin_xcoord > id_thresh) %>% group_by(SRA, bin) %>% 
  mutate(cum_rel_reads_mapped = cumsum(n))%>% 
  dplyr::filter(bin_xcoord == 100-0.25)

map_disc_cum <- left_join(map_disc_cum, meta_Ext, by = c("SRA" = "accession.number"))

p_sdisc_cum3 <- ggplot(map_disc_cum, aes(x = bin, 
                                         y = 100*cum_rel_reads_mapped/1e6, 
                                        fill = bin))+
  theme_bw()+
  scale_fill_manual("", values = c(col_bac1,col_bac2,col_RAMLI))+
  geom_jitter(size = 4, shape = 21, color = "black", alpha = 0.7, width = 0.15)+
  geom_boxplot(fill = NA, outlier.shape = NA)+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
      title=element_text(size=20), legend.text=element_text(size=12),
      legend.background = element_rect(fill="transparent"),
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
      strip.text=element_text(size=14), legend.position = "bottom",
      strip.background = element_rect(fill = adjustcolor("gray", 0.15)))+
  ylab(paste0("Relative abundance\n ( > ", id_thresh-0.25, "% NI)"))+
  xlab("")+
  guides(fill=FALSE)+
  facet_grid(.~Water_type, scales ="free")

print(p_sdisc_cum3)
```

<img src="Figures/cached/SDP_external-2-1.png" style="display: block; margin: auto;" />


```r
df_BlastExt2 <- data.table::fread("./ExternalData/merged_blastfiles.tsv", 
                         header = TRUE)

# Import contig lengths
contigs_bac1 <- data.table::fread("./ExternalData/bacI_contig_lengths.txt", 
                         header = FALSE) %>% 
  top_n(3, V2) %>% 
  dplyr::select(V1)
contigs_bac2 <- data.table::fread("./ExternalData/bacII_contig_lengths.txt", 
                         header = FALSE) %>% 
  top_n(3, V2) %>% 
  dplyr::select(V1)
contigs_Ramli <- data.table::fread("./ExternalData/Ramli_contig_lengths.txt", 
                         header = FALSE) %>% 
  top_n(3, V2) %>% 
  dplyr::select(V1)

# Filter out metagenomes with at least 1,000 reads
metaG_5k <- df_BlastExt2 %>% 
  dplyr::select(SRA) %>% 
  table() %>% 
  data.frame() %>% 
  dplyr::filter(Freq>5e3) %>% 
  dplyr::select(".")
colnames(metaG_5k) <- "Filtered_metaG"
df_BlastExt2 <- df_BlastExt2 %>% dplyr::filter(SRA %in% metaG_5k$Filtered_metaG)

# Filter out largest contigs (top N)
df_BlastExt2 <- df_BlastExt2 %>% dplyr::filter(contig_id %in% rbind(contigs_Ramli,
                                                                contigs_bac1,
                                                                contigs_bac2)$V1)


# add metadata
df_BlastExt2 <- dplyr::left_join(df_BlastExt2, contig_ids, by = "contig_id")

# Plot first alignment point
p_sdp_aln1 <- df_BlastExt2 %>% dplyr::filter(bin == "Ramlibacter sp. MAG") %>% 
  ggplot(., aes(x = sstart, y = identity, fill = bin))+
  theme_bw()+
  scale_fill_manual("", values = c(col_RAMLI))+
  geom_point(shape = 21, size = 2)+
  geom_smooth(color = "black", size = 3)+
 # geom_density(aes(x=sstart, y=..scaled.., fill=bin),
 #               alpha = 0.4, size = 0.4, color = "#333333",
 #               bw = "nrd0", trim = TRUE)+
  facet_grid(bin~contig_id, scales = "free_x")+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
      title=element_text(size=20), legend.text=element_text(size=12),
      legend.background = element_rect(fill="transparent"),
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
      strip.text=element_text(size=12), legend.position = "bottom",
      strip.background = element_rect(fill = adjustcolor("gray", 0.15)))+
  ylab(paste0("% identity"))+
  xlab("Position in alignment")+
  guides(fill=FALSE)

p_sdp_aln2 <- df_BlastExt2 %>% dplyr::filter(bin == "Bacteroidetes sp. MAG1") %>% 
  ggplot(., aes(x = sstart, y = identity, fill = bin))+
  theme_bw()+
  scale_fill_manual("", values = c(col_bac1))+
  geom_point(shape = 21, size = 2)+
  geom_smooth(color = "black", size = 3)+
  # geom_density(aes(x=sstart, y=..scaled.., fill=bin),
  #              alpha = 0.4, size = 0.4, color = "#333333",
  #              bw = "nrd0", trim = TRUE)+
  facet_grid(bin~contig_id, scales = "free_x")+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
      title=element_text(size=20), legend.text=element_text(size=12),
      legend.background = element_rect(fill="transparent"),
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
      strip.text=element_text(size=12), legend.position = "bottom",
      strip.background = element_rect(fill = adjustcolor("gray", 0.15)))+
  ylab(paste0("% identity"))+
  xlab("Position in alignment")+
  guides(fill=FALSE)


p_sdp_aln3 <- df_BlastExt2 %>% dplyr::filter(bin == "Bacteroidetes sp. MAG2") %>% 
  ggplot(., aes(x = sstart, y = identity, fill = bin))+
  theme_bw()+
  scale_fill_manual("", values = c(col_bac2))+
  geom_point(shape = 21, size = 2)+
  geom_smooth(color = "black", size = 3)+
  # geom_density(aes(x=sstart, y=..scaled.., fill=bin),
  #              alpha = 0.4, size = 0.4, color = "#333333",
  #              bw = "nrd0", trim = TRUE)+
  facet_grid(bin~contig_id, scales = "free_x")+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
      title=element_text(size=20), legend.text=element_text(size=12),
      legend.background = element_rect(fill="transparent"),
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
      strip.text=element_text(size=12), legend.position = "bottom",
      strip.background = element_rect(fill = adjustcolor("gray", 0.15)))+
  ylab(paste0("% identity"))+
  xlab("Position in alignment")+
  guides(fill=FALSE)

cowplot::plot_grid(p_sdp_aln1, p_sdp_aln2 ,p_sdp_aln3, align = "hv", nrow = 3)
```

```
## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'
```

```
## `geom_smooth()` using method = 'loess' and formula 'y ~ x'
```

```
## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'
```

<img src="Figures/cached/SDP_alignment-check-1-1.png" style="display: block; margin: auto;" />

```r
# Plot for each SRA data set seperately
for(dataset in 1:length(metaG_5k$Filtered_metaG)){
  # Plot first alignment point
p_sdp_aln1 <- df_BlastExt2 %>% 
  dplyr::filter(SRA == metaG_5k$Filtered_metaG[dataset]) %>% 
  dplyr::filter(bin == "Ramlibacter sp. MAG") %>% 
  ggplot(., aes(x = sstart, y = identity, fill = bin))+
  theme_bw()+
  scale_fill_manual("", values = c(col_RAMLI))+
  geom_point(shape = 21, size = 2)+
  geom_smooth(color = "black", size = 3)+
  facet_grid(bin~contig_id, scales = "free_x")+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
      title=element_text(size=14), legend.text=element_text(size=12),
      legend.background = element_rect(fill="transparent"),
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
      strip.text=element_text(size=12), legend.position = "bottom",
      strip.background = element_rect(fill = adjustcolor("gray", 0.15)))+
  ylab(paste0("% identity"))+
  xlab("Position in alignment")+
  guides(fill=FALSE)+
  ggtitle(metaG_5k$Filtered_metaG[dataset])+
  ylim(75,100)

p_sdp_aln2 <- df_BlastExt2 %>% 
  dplyr::filter(SRA == metaG_5k$Filtered_metaG[dataset]) %>% 
  dplyr::filter(bin == "Bacteroidetes sp. MAG1") %>% 
  ggplot(., aes(x = sstart, y = identity, fill = bin))+
  theme_bw()+
  scale_fill_manual("", values = c(col_bac1))+
  geom_point(shape = 21, size = 2)+
  geom_smooth(color = "black", size = 3)+
  facet_grid(bin~contig_id, scales = "free_x")+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
      title=element_text(size=14), legend.text=element_text(size=12),
      legend.background = element_rect(fill="transparent"),
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
      strip.text=element_text(size=12), legend.position = "bottom",
      strip.background = element_rect(fill = adjustcolor("gray", 0.15)))+
  ylab(paste0("% identity"))+
  xlab("Position in alignment")+
  guides(fill=FALSE)+
  ggtitle(metaG_5k$Filtered_metaG[dataset])+
  ylim(75,100)

p_sdp_aln3 <- df_BlastExt2 %>% 
  dplyr::filter(SRA == metaG_5k$Filtered_metaG[dataset]) %>% 
  dplyr::filter(bin == "Bacteroidetes sp. MAG2") %>% 
  ggplot(., aes(x = sstart, y = identity, fill = bin))+
  theme_bw()+
  scale_fill_manual("", values = c(col_bac2))+
  geom_point(shape = 21, size = 2)+
  geom_smooth(color = "black", size = 3)+
  facet_grid(bin~contig_id, scales = "free_x")+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
      title=element_text(size=14), legend.text=element_text(size=12),
      legend.background = element_rect(fill="transparent"),
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
      strip.text=element_text(size=12), legend.position = "bottom",
      strip.background = element_rect(fill = adjustcolor("gray", 0.15)))+
  ylab(paste0("% identity"))+
  xlab("Position in alignment")+
  guides(fill=FALSE)+
  ggtitle(metaG_5k$Filtered_metaG[dataset])+
  ylim(75,100)

# now add the title
# title <- ggdraw() + draw_label(metaG_5k$Filtered_metaG[dataset], fontface='bold')
# plot_grid(title, p, ncol=1, rel_heights=c(0.1, 1)) 
png(paste0("./ExternalData/Figures/",gsub(".tsv","",metaG_5k$Filtered_metaG[dataset]),
           ".png"), height = 5, width = 10, res = 500, units = "in")
print(p_sdp_aln1)
# print(cowplot::plot_grid(p_sdp_aln1, p_sdp_aln2, p_sdp_aln3, align = "hv", nrow = 3))
dev.off()
}
```

```
## `geom_smooth()` using method = 'loess' and formula 'y ~ x'
```

```
## `geom_smooth()` using method = 'loess' and formula 'y ~ x'
## `geom_smooth()` using method = 'loess' and formula 'y ~ x'
## `geom_smooth()` using method = 'loess' and formula 'y ~ x'
## `geom_smooth()` using method = 'loess' and formula 'y ~ x'
## `geom_smooth()` using method = 'loess' and formula 'y ~ x'
## `geom_smooth()` using method = 'loess' and formula 'y ~ x'
## `geom_smooth()` using method = 'loess' and formula 'y ~ x'
## `geom_smooth()` using method = 'loess' and formula 'y ~ x'
## `geom_smooth()` using method = 'loess' and formula 'y ~ x'
## `geom_smooth()` using method = 'loess' and formula 'y ~ x'
## `geom_smooth()` using method = 'loess' and formula 'y ~ x'
## `geom_smooth()` using method = 'loess' and formula 'y ~ x'
## `geom_smooth()` using method = 'loess' and formula 'y ~ x'
## `geom_smooth()` using method = 'loess' and formula 'y ~ x'
## `geom_smooth()` using method = 'loess' and formula 'y ~ x'
## `geom_smooth()` using method = 'loess' and formula 'y ~ x'
## `geom_smooth()` using method = 'loess' and formula 'y ~ x'
## `geom_smooth()` using method = 'loess' and formula 'y ~ x'
## `geom_smooth()` using method = 'loess' and formula 'y ~ x'
```

```
## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'
## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'
## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'
## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'
## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'
## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'
## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'
## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'
## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'
## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'
## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'
```

```
## `geom_smooth()` using method = 'loess' and formula 'y ~ x'
```

```
## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'
```

```
## `geom_smooth()` using method = 'loess' and formula 'y ~ x'
```

```
## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'
```

```
## `geom_smooth()` using method = 'loess' and formula 'y ~ x'
```

```
## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'
```

```
## `geom_smooth()` using method = 'loess' and formula 'y ~ x'
## `geom_smooth()` using method = 'loess' and formula 'y ~ x'
```

```
## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'
```

```
## `geom_smooth()` using method = 'loess' and formula 'y ~ x'
```

```
## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'
## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'
## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'
## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'
```

```
## `geom_smooth()` using method = 'loess' and formula 'y ~ x'
## `geom_smooth()` using method = 'loess' and formula 'y ~ x'
## `geom_smooth()` using method = 'loess' and formula 'y ~ x'
## `geom_smooth()` using method = 'loess' and formula 'y ~ x'
## `geom_smooth()` using method = 'loess' and formula 'y ~ x'
## `geom_smooth()` using method = 'loess' and formula 'y ~ x'
## `geom_smooth()` using method = 'loess' and formula 'y ~ x'
```

```
## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'
```

```
## `geom_smooth()` using method = 'loess' and formula 'y ~ x'
```

```
## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'
```


```r
# Only Ramlibacter results now
p_sdp_aln1
```

<img src="Figures/cached/SDP_alignment-check-2-1.png" style="display: block; margin: auto;" />
