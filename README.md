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
