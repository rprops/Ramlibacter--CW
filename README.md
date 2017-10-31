# Limno_CoolingWater
Repository for the analysis of reconstructed genomes from secondary cooling water of the BR2 nuclear test reactor (Mol, Antwerp, Belgium).

# Analysis steps

## 1. QC and assembly

### Coassembly with IDBA-ID with on interleaved fasta with following parameters


### Map reads to coassembly
```
#/bin/bash

set -e

for file in `cat SAMPLE.list`
do
bwa mem -t 20 ./anvio/contigs-fixed.fa ../${file}/derep_scythe_sickle_fwd_${file}.fastq_pairs_R1.fastq ../${file}/derep_scythe_sickle_rev_${file}.fastq_pairs_R2.fastq > ./anvio/map/${file}.sam
samtools view -h -b -S ./anvio/map/${file}.sam > ./anvio/map/${file}-raw.bam
anvi-init-bam ./anvio/map/${file}-raw.bam -o ./anvio/map/${file}.bam
name=${file//-/_}
anvi-profile -i ./anvio/map/${file}.bam -c ./anvio/contigs.db --min-contig-length 1000 --output-dir ./anvio/$file --sample-name $name --min-coverage-for-variability 5 
done
```


## 2. Binning strategy

## 3. Binning QC and refinement

## 4. Calculation of relative abundances

# Reconstruction of full-length 16S sequences (`EMIRGE`)
`EMIRGE` was run on the merged QC'd fastq files with insert size 500 and st.dev 500 in order to achieve maximum mapping of reads. Futher, `EMIRGE` was run both with the full NR database clustered at 97% using `UCLUST` and with the small manually curated freshwater database (FWDB) available from [here](https://github.com/McMahonLab/TaxAss/tree/master/FreshTrain-files). In addition emirge was run separately with `-j` on 0.97 and 1.0 allowing both 97% consensus sequences and unique sequences to be reconstructed. From all these runs a merged fasta file was constructed. Finally, all the reconstructed sequences were clustered at 97% identity using `UCLUST`. `EMIRGE` reconstructed sequences with a normalized prior abundance of less than 5% were removed, and sequences were ordered from high to low abundance before clustering because `UCLUST` is dependent on the order of the sequences (see [here](https://www.drive5.com/usearch/manual/uclust_algo.html)). Classification of the sequences was performed using the `TaxAss` pipeline as described [here](https://github.com/rprops/Mothur_oligo_batch). Full-length sequences were classified using the FWDB database if their top matches to the FWDB were higher than 95%.  

**Emirge reconstructed sequences**  

```
NR_1.0
15|HQ222271.1.1558_Prior=0.791329_Length=1536_NormPrior=0.784664        Bacteria(100);Proteobacteria(100);Betaproteobacteria(100);Burkholderiales(100);Comamonadaceae(100);Comamonadaceae_unclassified(100);
16|EF516083.1.1452_Prior=0.020773_Length=1440_NormPrior=0.021971        Bacteria(100);Proteobacteria(100);Betaproteobacteria(100);Burkholderiales(100);Comamonadaceae(100);Comamonadaceae_unclassified(100);
0|FR853751.1.1492_Prior=0.187898_Length=1480_NormPrior=0.193365 Bacteria(100);Bacteroidetes(100);Sphingobacteriia(100);Sphingobacteriales(100);Chitinophagaceae(100);Sediminibacterium(100);

NR_0.97
15|HQ222271.1.1558_Prior=0.812102_Length=1536_NormPrior=0.806369        Bacteria(100);Proteobacteria(100);Betaproteobacteria(100);Burkholderiales(100);Comamonadaceae(100);Comamonadaceae_unclassified(100);
0|FR853751.1.1492_Prior=0.187898_Length=1480_NormPrior=0.193631 Bacteria(100);Bacteroidetes(100);Sphingobacteriia(100);Sphingobacteriales(100);Chitinophagaceae(100);Sediminibacterium(100);

FWDB_0.97
5_Bctrm474_Prior=1.000000_Length=1481_NormPrior=1.000000        Bacteria(100);Proteobacteria(100);Betaproteobacteria(100);Burkholderiales(100);Comamonadaceae(100);Ramlibacter(91);

FWDB_1.0
5|Bctrm474_Prior=0.097627_Length=1481_NormPrior=0.097627        Bacteria(100);Proteobacteria(100);Betaproteobacteria(100);Burkholderiales(100);Comamonadaceae(100);Comamonadaceae_unclassified(100);
22|LimCurv2_Prior=0.015439_Length=1481_NormPrior=0.015439       Bacteria(100);Proteobacteria(100);Betaproteobacteria(100);Burkholderiales(100);Comamonadaceae(100);Comamonadaceae_unclassified(100);
18|LimSpe14_Prior=0.886934_Length=1481_NormPrior=0.886934       Bacteria(100);Proteobacteria(100);Betaproteobacteria(100);Burkholderiales(100);Comamonadaceae(100);Comamonadaceae_unclassified(100);
```

**Clustering of merged sequence file**  

```
usearch -cluster_fast merged.emirged.cons.fasta -id 0.97 -centroids merged.emirged.cons.0.97.fasta -uc merged.emirged.cons.0.97.clusters.uc
```

**Final dereplicated consensus EMIRGE sequences**
```
5_Bctrm474_Prior=1.000000_Length=1481_NormPrior=1.000000        Bacteria(100);Proteobacteria(100);Betaproteobacteria(100);Burkholderiales(100);Comamonadaceae(100);Ramlibacter(91);
0|FR853751.1.1492_Prior=0.187898_Length=1480_NormPrior=0.193365 Bacteria(100);Bacteroidetes(100);Sphingobacteriia(100);Sphingobacteriales(100);Chitinophagaceae(100);Sediminibacterium(100);
```
