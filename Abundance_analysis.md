# Metagenomic analysis of secondary cooling water microbial communities
Ruben Props  
August 17 2017  





```r
# Read data
mean_coverage <- read.table("./SAMPLES-SUMMARY/bins_across_samples/mean_coverage.txt", header = TRUE)
std_coverage <- read.table("./SAMPLES-SUMMARY/bins_across_samples/std_coverage.txt", header = TRUE)
bin_size <- read.table("./SAMPLES-SUMMARY/general_bins_summary.txt", header = TRUE)[, c(1,3)]
total_reads <- read.table("./sample_reads.tsv", header = TRUE)
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
data_total$bins <- plyr::revalue(data_total$bins, c("BetIa_bin"="Limnohabitans_bin", "bacIa_vizbin1"="Bacteroidetes_bin1",
                                              "bacIa_vizbin2"="Bacteroidetes_bin2"))
# Calculate relative abundance of the bins
data_total$mean_rel_abundance <- 100*(data_total$coverage*data_total$bin_size)/(read_length*data_total$Total_reads)
data_total$upper_rel_abundance <- 100*((data_total$coverage+data_total$std_coverage)*data_total$bin_size)/(read_length*data_total$Total_reads)
data_total$lower_rel_abundance <- 100*((data_total$coverage-data_total$std_coverage)*data_total$bin_size)/(read_length*data_total$Total_reads)

data_total$mean_rel_abundance_map <- 100*(data_total$coverage*data_total$bin_size)/(read_length*data_total$Mapped_reads)
data_total$upper_rel_abundance_map <- 100*((data_total$coverage+data_total$std_coverage)*data_total$bin_size)/(read_length*data_total$Mapped_reads)
data_total$lower_rel_abundance_map <- 100*((data_total$coverage-data_total$std_coverage)*data_total$bin_size)/(read_length*data_total$Mapped_reads)
```

# *1. Investigate MAG- and 16S-based abundances*
It is clear that there is significant %GC coverage bias present. The estimated relative abundances
from metagenomics do not quantitatively match with the V3-V4 16S rRNA gene amplicon data.
$$Relative\ abundance =100*(\frac{mean\ coverage * bin\ size}{read\ length*total\ sample\ reads })$$
Another option is to calculate relative to mapped number of reads:
$$Relative\ abundance =100*(\frac{mean\ coverage * bin\ size}{read\ length*total\ sample\ reads * \%mapped\ reads})$$


Import reference relative abundances from 16S data set in order to directly compare with metagenomic data set.

```r
df_16S <- read.table("relative_abundance_16S.tsv", header = TRUE)
df_16S_long <- gather(df_16S, Sample_ID, relative_abundance_16S, 
                             SAMPLE_16:SAMPLE_65, factor_key=TRUE)
```



```r
# Subset for only the three complete genomes (not PNCs).
data_total_sb <- data_total[data_total$bins %in% c("Limnohabitans_bin", "Bacteroidetes_bin1", "Bacteroidetes_bin2"),]

p_meta <- ggplot(data = data_total_sb, aes(x = bins, y = mean_rel_abundance, fill = bins))+
  geom_point(size = 4, shape = 21, alpha = 0.7)+
  scale_fill_manual(values = c(col_bac1, col_bac2, col_limno))+
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
  scale_fill_manual(values = c(col_bac1, col_bac2, col_limno))+
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
  scale_fill_manual(values = c(col_bac1, col_bac2, col_limno))+
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

# *2. Investigate sequence characteristics within coding DNA sequences (CDS)*


```r
# First we need the files that map the gene ID to the sequence ID (linux cmd: https://github.com/rprops/MetaG_lakeMI/wiki/11.-Genome-annotation)
# These are stored in the IMG_annotation data for each genome bin

# Next, extract the %GC of each gene from the gff file
extract_gc_from_gff("./IMG_annotation/IMG_2724679690_Limnohabitans/121950.assembled.gff", outputFolder = "GC_analysis")
extract_gc_from_gff("./IMG_annotation/IMG_2724679691_Bacteroidetes_bin1/121951.assembled.gff", outputFolder = "GC_analysis")
extract_gc_from_gff("./IMG_annotation/IMG_2724679698_Bacteroidetes_bin2/121960.assembled.gff", outputFolder = "GC_analysis")

# Use these files to make dataframes mapping function (COGs) and %GC
LIMNO_gc_cog <- gc2function(seq_id_gc = "GC_analysis/seqid_GC_121950.assembled.gff.tsv", gene_id_seq_id ="./IMG_annotation/IMG_2724679690_Limnohabitans/Annotation/2724679690_gene_oid_2_seq_id.txt", 
                             functions = "./IMG_annotation/IMG_2724679690_Limnohabitans/Annotation/2724679690.cog.tab.txt", gc_thresh = 0.1, output = FALSE)
```

```
## Fri Aug 25 18:52:48 2017  --- There are 2830 genes with > 0.1 %
## Fri Aug 25 18:52:48 2017  --- This is 100 % of all genes
## Fri Aug 25 18:52:48 2017  --- The 10 genes with the highest GC% are:
##      function_id                                             function_name
## 2821     COG0405                              Gamma-glutamyltranspeptidase
## 2822     COG2755                  Lysophospholipase L1 or related esterase
## 2823     COG0642                      Signal transduction histidine kinase
## 2824     COG1514                                          2'-5' RNA ligase
## 2825     COG1767                Triphosphoribosyl-dephospho-CoA synthetase
## 2826     COG1261         Flagella basal body P-ring formation protein FlgA
## 2827     COG0665                Glycine/D-amino acid oxidase (deaminating)
## 2828     COG0810 Periplasmic protein TonB, links inner and outer membranes
## 2829     COG1040                 Predicted amidophosphoribosyltransferases
## 2830     COG1240                                 Mg-chelatase subunit ChlD
##        GC
## 2821 78.5
## 2822 78.5
## 2823 78.5
## 2824 78.5
## 2825 78.6
## 2826 78.6
## 2827 79.1
## 2828 79.5
## 2829 79.7
## 2830 80.6
```

```r
BAC1_gc_cog <- gc2function(seq_id_gc = "GC_analysis/seqid_GC_121951.assembled.gff.tsv", gene_id_seq_id ="./IMG_annotation/IMG_2724679691_Bacteroidetes_bin1/Annotation/2724679691_gene_oid_2_seq_id.txt", 
                             functions = "./IMG_annotation/IMG_2724679691_Bacteroidetes_bin1//Annotation/2724679691.cog.tab.txt", gc_thresh = 0.1, output = FALSE)
```

```
## Fri Aug 25 18:52:48 2017  --- There are 1889 genes with > 0.1 %
## Fri Aug 25 18:52:48 2017  --- This is 100 % of all genes
## Fri Aug 25 18:52:48 2017  --- The 10 genes with the highest GC% are:
##      function_id
## 1880     COG0052
## 1881     COG0183
## 1882     COG0629
## 1883     COG0509
## 1884     COG1734
## 1885     COG0636
## 1886     COG3502
## 1887     COG1765
## 1888     COG0377
## 1889     COG4104
##                                                                               function_name
## 1880                                                                   Ribosomal protein S2
## 1881                                                           Acetyl-CoA acetyltransferase
## 1882                                                    Single-stranded DNA-binding protein
## 1883                                    Glycine cleavage system H protein (lipoate-binding)
## 1884                                       RNA polymerase-binding transcription factor DksA
## 1885 FoF1-type ATP synthase, membrane subunit c/Archaeal/vacuolar-type H+-ATPase, subunit K
## 1886                                       Uncharacterized conserved protein, DUF952 family
## 1887                                                   Uncharacterized OsmC-related protein
## 1888 NADH:ubiquinone oxidoreductase 20 kD subunit (chhain B) or related Fe-S oxidoreductase
## 1889                 Zn-binding Pro-Ala-Ala-Arg (PAAR) domain, incolved in TypeVI secretion
##        GC
## 1880 51.6
## 1881 51.6
## 1882 52.0
## 1883 52.6
## 1884 52.7
## 1885 52.8
## 1886 52.9
## 1887 53.6
## 1888 53.6
## 1889 59.3
```

```r
BAC2_gc_cog <- gc2function(seq_id_gc = "GC_analysis/seqid_GC_121960.assembled.gff.tsv", gene_id_seq_id ="./IMG_annotation/IMG_2724679698_Bacteroidetes_bin2/Annotation/2724679698_gene_oid_2_seq_id.txt", 
                             functions = "./IMG_annotation/IMG_2724679698_Bacteroidetes_bin2/Annotation/2724679698.cog.tab.txt", gc_thresh = 0.1, output = FALSE)
```

```
## Fri Aug 25 18:52:48 2017  --- There are 1797 genes with > 0.1 %
## Fri Aug 25 18:52:48 2017  --- This is 100 % of all genes
## Fri Aug 25 18:52:48 2017  --- The 10 genes with the highest GC% are:
##      function_id
## 1788     COG4675
## 1789     COG0636
## 1790     COG1501
## 1791     COG0725
## 1792     COG5434
## 1793     COG2115
## 1794     COG4225
## 1795     COG3669
## 1796     COG4588
## 1797     COG3258
##                                                                               function_name
## 1788                                      Microcystin-dependent protein  (function unknown)
## 1789 FoF1-type ATP synthase, membrane subunit c/Archaeal/vacuolar-type H+-ATPase, subunit K
## 1790                                      Alpha-glucosidase, glycosyl hydrolase family GH31
## 1791                             ABC-type molybdate transport system, periplasmic component
## 1792                                                                      Polygalacturonase
## 1793                                                                       Xylose isomerase
## 1794                                                      Rhamnogalacturonyl hydrolase YesR
## 1795                                                                     Alpha-L-fucosidase
## 1796               Accessory colonization factor AcfC, contains ABC-type periplasmic domain
## 1797                                                                           Cytochrome c
##        GC
## 1788 44.5
## 1789 44.6
## 1790 44.9
## 1791 45.0
## 1792 45.1
## 1793 45.1
## 1794 45.4
## 1795 45.5
## 1796 46.6
## 1797 46.8
```

```r
LIMNO_gc_pfam <- gc2function(seq_id_gc = "GC_analysis/seqid_GC_121950.assembled.gff.tsv", gene_id_seq_id ="./IMG_annotation/IMG_2724679690_Limnohabitans/Annotation/2724679690_gene_oid_2_seq_id.txt", 
                             functions = "./IMG_annotation/IMG_2724679690_Limnohabitans/Annotation/2724679690.pfam.tab.txt", gc_thresh = 0.1, output = FALSE)
```

```
## Fri Aug 25 18:52:48 2017  --- There are 4954 genes with > 0.1 %
## Fri Aug 25 18:52:48 2017  --- This is 100 % of all genes
## Fri Aug 25 18:52:48 2017  --- The 10 genes with the highest GC% are:
##      function_id function_name   GC
## 4945   pfam13202     EF-hand_5 79.0
## 4946   pfam16537         T2SSB 79.0
## 4947   pfam01266           DAO 79.1
## 4948   pfam03544        TonB_C 79.5
## 4949   pfam00156  Pribosyltran 79.7
## 4950   pfam11142       DUF2917 79.8
## 4951   pfam13318       DUF4089 79.8
## 4952   pfam13519         VWA_2 80.6
## 4953   pfam02120      Flg_hook 80.9
## 4954   pfam03023          MVIN 81.7
```

```r
BAC1_gc_pfam <- gc2function(seq_id_gc = "GC_analysis/seqid_GC_121951.assembled.gff.tsv", gene_id_seq_id ="./IMG_annotation/IMG_2724679691_Bacteroidetes_bin1/Annotation/2724679691_gene_oid_2_seq_id.txt", 
                             functions = "./IMG_annotation/IMG_2724679691_Bacteroidetes_bin1//Annotation/2724679691.pfam.tab.txt", gc_thresh = 0.1, output = FALSE)
```

```
## Fri Aug 25 18:52:49 2017  --- There are 3929 genes with > 0.1 %
## Fri Aug 25 18:52:49 2017  --- This is 100 % of all genes
## Fri Aug 25 18:52:49 2017  --- The 10 genes with the highest GC% are:
##      function_id function_name   GC
## 3920   pfam02803    Thiolase_C 51.6
## 3921   pfam00436           SSB 52.0
## 3922   pfam00171        Aldedh 52.2
## 3923   pfam01597         GCV_H 52.6
## 3924   pfam01258  zf-dskA_traR 52.7
## 3925   pfam00137    ATP-synt_C 52.8
## 3926   pfam06108        DUF952 52.9
## 3927   pfam02566          OsmC 53.6
## 3928   pfam01058   Oxidored_q6 53.6
## 3929   pfam05488    PAAR_motif 59.3
```

```r
BAC2_gc_pfam <- gc2function(seq_id_gc = "GC_analysis/seqid_GC_121960.assembled.gff.tsv", gene_id_seq_id ="./IMG_annotation/IMG_2724679698_Bacteroidetes_bin2/Annotation/2724679698_gene_oid_2_seq_id.txt", 
                             functions = "./IMG_annotation/IMG_2724679698_Bacteroidetes_bin2/Annotation/2724679698.pfam.tab.txt", gc_thresh = 0.1, output = FALSE)
```

```
## Fri Aug 25 18:52:49 2017  --- There are 3573 genes with > 0.1 %
## Fri Aug 25 18:52:49 2017  --- This is 100 % of all genes
## Fri Aug 25 18:52:49 2017  --- The 10 genes with the highest GC% are:
##      function_id   function_name   GC
## 3564   pfam13531      SBP_bac_11 46.6
## 3565   pfam13442 Cytochrome_CBB3 46.8
## 3566   pfam12779          YXWGXW 50.0
## 3567   pfam12779          YXWGXW 50.0
## 3568   pfam01391        Collagen 50.4
## 3569   pfam01391        Collagen 50.4
## 3570   pfam01391        Collagen 50.4
## 3571   pfam01391        Collagen 50.4
## 3572   pfam01391        Collagen 50.4
## 3573   pfam01391        Collagen 50.4
```

Motivation: For COGs there exists a hierarchy allowing us to investigate whether there is a conservation of high/low %GC in certain functional gene groups. In order to do this we need to incorporate this hierarchy into the genome dataframes we have now.
  




# *3. Analysis of gene length distribution*
Here we use the dataframe made in the previous section to see if there is a significant difference in the gene length of the COGs within these three consensus genomes.  

Observation: They have very small genes: on average < 500bp.


```r
# Limnohabitans MAG gene length distribution 
p_LIMNO_length <- easyGgplot2::ggplot2.histogram(data = LIMNO_gc_cog, xName = 'gene_length',
                  groupName = 'Genome', alpha = 0.5,
                  legendPosition = "top", binwidth = 0.15,
                  groupColors = col_limno,addMeanLine=TRUE, meanLineColor="black",
                  meanLineType="dashed")+ theme_bw()+ ylim(0,15)+
  labs(x = "Gene length (bp)", y = "Count")+ theme(legend.position="none")+
  ggtitle("Limnohabitans MAG")+ xlim(0,2000)

# Bacteroidetes MAG1 gene length distribution 
p_BAC1_length <- easyGgplot2::ggplot2.histogram(data = BAC1_gc_cog, xName = 'gene_length',
                  groupName = 'Genome', alpha = 0.5,
                  legendPosition = "top", binwidth = 0.15, 
                  groupColors = col_bac1,addMeanLine=TRUE, meanLineColor="black",
                  meanLineType="dashed")+ theme_bw()+ ylim(0,15)+
  labs(x = "Gene length (bp)", y = "Count")+ theme(legend.position="none")+
  ggtitle("Bacteroidetes MAG1")+ xlim(0,2000)

# Bacteroidetes MAG2 gene length distribution 
p_BAC2_length <- easyGgplot2::ggplot2.histogram(data = BAC2_gc_cog, xName = 'gene_length',
                  groupName = 'Genome', alpha = 0.5,
                  legendPosition = "top", binwidth = 0.15,
                  groupColors = col_bac2,addMeanLine=TRUE, meanLineColor="black",
                  meanLineType="dashed")+ theme_bw()+ ylim(0,15)+
  labs(x = "Gene length (bp)", y = "Count")+ theme(legend.position="none")+
  ggtitle("Bacteroidetes MAG2")+ xlim(0,2000)

grid.arrange(p_LIMNO_length, p_BAC1_length, p_BAC2_length, ncol = 3)
```

<img src="Figures/cached/gene-length-analysis_COG-1.png" style="display: block; margin: auto;" />

We can do the same for the Pfams.

```r
# Limnohabitans MAG gene length distribution 
p_LIMNO_length <- easyGgplot2::ggplot2.histogram(data = LIMNO_gc_pfam, xName = 'gene_length',
                  groupName = 'Genome', alpha = 0.5,
                  legendPosition = "top", binwidth = 0.15,
                  groupColors = col_limno,addMeanLine=TRUE, meanLineColor="black",
                  meanLineType="dashed")+ theme_bw()+ ylim(0,25)+
  labs(x = "Gene length (bp)", y = "Count")+ theme(legend.position="none")+
  ggtitle("Limnohabitans MAG")+ xlim(0,3000)

# Bacteroidetes MAG1 gene length distribution 
p_BAC1_length <- easyGgplot2::ggplot2.histogram(data = BAC1_gc_pfam, xName = 'gene_length',
                  groupName = 'Genome', alpha = 0.5,
                  legendPosition = "top", binwidth = 0.15, 
                  groupColors = col_bac1,addMeanLine=TRUE, meanLineColor="black",
                  meanLineType="dashed")+ theme_bw()+ ylim(0,25)+
  labs(x = "Gene length (bp)", y = "Count")+ theme(legend.position="none")+
  ggtitle("Bacteroidetes MAG1")+ xlim(0,3000)

# Bacteroidetes MAG2 gene length distribution 
p_BAC2_length <- easyGgplot2::ggplot2.histogram(data = BAC2_gc_pfam, xName = 'gene_length',
                  groupName = 'Genome', alpha = 0.5,
                  legendPosition = "top", binwidth = 0.15,
                  groupColors = col_bac2,addMeanLine=TRUE, meanLineColor="black",
                  meanLineType="dashed")+ theme_bw()+ ylim(0,25)+
  labs(x = "Gene length (bp)", y = "Count")+ theme(legend.position="none")+
  ggtitle("Bacteroidetes MAG2")+ xlim(0,3000)

grid.arrange(p_LIMNO_length, p_BAC1_length, p_BAC2_length, ncol = 3)
```

<img src="Figures/cached/gene-length-analysis_Pfam-1.png" style="display: block; margin: auto;" />

# 5. Identify unique functional genes (COG/Pfams)

```r
# Find unique functions in Limnohabitans MAG vs Bacteroidetes MAG1
unique_LIMNO_BAC1 <- dplyr::anti_join(LIMNO_gc_cog, BAC1_gc_cog, by = "cog_id")
cat("There are", paste(nrow(unique_LIMNO_BAC1)), "unique COGs in Limnohabitans MAG vs Bacteroidetes MAG1")
```

```
## There are 1178 unique COGs in Limnohabitans MAG vs Bacteroidetes MAG1
```

```r
# Find unique functions in Limnohabitans MAG vs Bacteroidetes MAG2
unique_LIMNO_BAC2 <- dplyr::anti_join(LIMNO_gc_cog, BAC2_gc_cog, by = "cog_id")
cat("There are", paste(nrow(unique_LIMNO_BAC2)), "unique COGs in Limnohabitans MAG vs Bacteroidetes MAG2")
```

```
## There are 1143 unique COGs in Limnohabitans MAG vs Bacteroidetes MAG2
```

```r
# Find unique functions in Bacteroidetes MAG1 vs Bacteroidetes MAG2
unique_BAC1_BAC2 <- dplyr::anti_join(BAC1_gc_cog, BAC2_gc_cog, by = "cog_id")
cat("There are", paste(nrow(unique_BAC1_BAC2)), "unique COGs in Bacteroidetes MAG1 vs Bacteroidetes MAG2")
```

```
## There are 153 unique COGs in Bacteroidetes MAG1 vs Bacteroidetes MAG2
```

```r
# Find unique functions in Bacteroidetes MAG1 vs Bacteroidetes MAG2
unique_BAC2_BAC1 <- dplyr::anti_join(BAC2_gc_cog, BAC1_gc_cog, by = "cog_id")
cat("There are", paste(nrow(unique_BAC2_BAC1)), "unique COGs in Bacteroidetes MAG2 vs Bacteroidetes MAG1")
```

```
## There are 143 unique COGs in Bacteroidetes MAG2 vs Bacteroidetes MAG1
```


```r
# Find unique functions in Limnohabitans MAG vs Bacteroidetes MAG1
unique_pfam_LIMNO_BAC1 <- dplyr::anti_join(LIMNO_gc_pfam, BAC1_gc_pfam, by = "pfam_id")
cat("There are", paste(nrow(unique_pfam_LIMNO_BAC1)), "unique Pfams in Limnohabitans MAG vs Bacteroidetes MAG1")
```

```
## There are 1474 unique Pfams in Limnohabitans MAG vs Bacteroidetes MAG1
```

```r
# Find unique functions in Limnohabitans MAG vs Bacteroidetes MAG2
unique_pfam_LIMNO_BAC2 <- dplyr::anti_join(LIMNO_gc_pfam, BAC2_gc_pfam, by = "pfam_id")
cat("There are", paste(nrow(unique_pfam_LIMNO_BAC2)), "unique Pfams in Limnohabitans MAG vs Bacteroidetes MAG2")
```

```
## There are 1424 unique Pfams in Limnohabitans MAG vs Bacteroidetes MAG2
```

```r
# Find unique functions in Bacteroidetes MAG1 vs Bacteroidetes MAG2
unique_pfam_BAC1_BAC2 <- dplyr::anti_join(BAC1_gc_pfam, BAC2_gc_pfam, by = "pfam_id")
cat("There are", paste(nrow(unique_pfam_BAC1_BAC2)), "unique Pfams in Bacteroidetes MAG1 vs Bacteroidetes MAG2")
```

```
## There are 285 unique Pfams in Bacteroidetes MAG1 vs Bacteroidetes MAG2
```

```r
# Find unique functions in Bacteroidetes MAG1 vs Bacteroidetes MAG2
unique_pfam_BAC2_BAC1 <- dplyr::anti_join(BAC2_gc_pfam, BAC1_gc_pfam, by = "pfam_id")
cat("There are", paste(nrow(unique_pfam_BAC2_BAC1)), "unique Pfams in Bacteroidetes MAG2 vs Bacteroidetes MAG1")
```

```
## There are 233 unique Pfams in Bacteroidetes MAG2 vs Bacteroidetes MAG1
```


# 6. COG functional categories
Get COG ID to COG functional category mapping file here: ftp://ftp.ncbi.nih.gov/pub/wolf/COGs/COG0303/cogs.csv    

The exact statistical analysis to compare genomes based on these profiles should be performed in STAMP.


```r
# Import COG mapping file
cogid_2_cogcat <- read.csv("./Mapping_files/cogid_2_cogcat.csv", sep = ",", header = FALSE, fill = TRUE,col.names = c("COG_ID", "COG_class", "function"))[, 1:2]
cogid_2_cogcat <- cogid_2_cogcat[(cogid_2_cogcat$COG_class)!="", ]
cogid_2_cogcat <- droplevels(cogid_2_cogcat)

# Read COG category file
cog_categories <- read.table("./Mapping_files/cog_categories.tsv", header = TRUE, sep = "\t")

# Merge COG metadata
cog_meta <- dplyr::left_join(cog_categories, cogid_2_cogcat, by = c("COG_class" = "COG_class"))
cog_meta <- droplevels(cog_meta)

# Merge genome information of all genome bins
merged_gc_cog <- rbind(LIMNO_gc_cog, BAC1_gc_cog, BAC2_gc_cog)
merged_gc_cog <- data.frame(merged_gc_cog, genome_id = c(rep("Limnohabitans MAG", nrow(LIMNO_gc_cog)), 
                                             rep("Bacteroidetes MAG1", nrow(BAC1_gc_cog)),
                                             rep("Bacteroidetes MAG2", nrow(BAC2_gc_cog)))
                            )

# Merge this metadata with the genome data from before
# COGs with multiple classifications are currently still NA - work on this.
merged_gc_cog <- dplyr::left_join(merged_gc_cog, cog_meta, by = c("cog_id" = "COG_ID"))
merged_gc_cog <- merged_gc_cog[!is.na(merged_gc_cog$COG_functional_category),]

# Visualize distribution across major metabolism functional COG groups per genome.
p_cog_func_group <- ggplot(data = merged_gc_cog, aes(x=COG_functional_category, fill = COG_functional_cluster))+
  geom_bar(stat="count", width=0.7, color = "black", size = 0.75)+
  theme_bw()+
  facet_grid(genome_id~.)+
  scale_fill_brewer(palette = "Accent")+
  labs(x = "Gene length (bp)", y = "Count")+ theme(legend.position="bottom", axis.text.x = element_text(angle = 90, hjust = 1),
                                                   legend.text = element_text(size = 7))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

print(p_cog_func_group)
```

<img src="Figures/cached/COG functional categories-1.png" style="display: block; margin: auto;" />

```r
p_cog_func_clust <- ggplot(data = merged_gc_cog, aes(x=COG_functional_cluster, fill = COG_functional_cluster))+
  geom_bar(stat="count", width=0.7, color = "black", size = 0.75)+
  theme_bw()+
  facet_grid(genome_id~.)+
  scale_fill_brewer(palette = "Accent")+
  labs(x = "Gene length (bp)", y = "Count")+ theme(legend.position="bottom",axis.text.x = element_text(angle = 90, hjust = 1),
                                                   legend.text = element_text(size = 7))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

print(p_cog_func_clust)
```

<img src="Figures/cached/COG functional categories-2.png" style="display: block; margin: auto;" />

