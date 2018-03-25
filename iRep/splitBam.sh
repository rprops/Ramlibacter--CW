#!/bin/bash
set -e
for bin in `ls rebinned_bins/final-refined_bins/B*.tsv`
do
	bin_stub1=${bin##*/}
	bin_stub=${bin_stub1%.fa.tsv}
	echo $bin_stub
	if [ ! -d "./splitBam/${bin_stub}" ]; then
  mkdir ./splitBam/${bin_stub}
  fi
	for sample_bam in `ls /scratch/vdenef_fluxm/rprops/DESMAN/metaG/vizbin_rebin_anvio_v230/map/*.mapped.bam`
		do
			stub1=${sample_bam##*/}
			stub=${stub1%.mapped.bam}
			echo $stub
			echo $bin
			echo $sample_bam
			samtools view -hL ${bin}.tsv $sample_bam > ./splitBam/${bin_stub}/${stub}.filtered.sam
	done
done

#!/bin/bash
set -e
for bin in `ls rebinned_bins/final-refined_bins/B*.fa`
do
	python /home/rprops/DESMAN/scripts/Lengths.py -i ${bin} > ${bin}.len
	grep ">" ${bin} | sed 's/>//g' > ${bin}.txt
	/home/rprops/StrainMetaSim/scripts/AddLengths.pl ${bin}.len < ${bin}.txt > ${bin}.tsv
do
done

python /home/rprops/DESMAN/scripts/Lengths.py -i B52_Sp13.BD.MM15.SD_rebin1.fa > B52_Sp13.BD.MM15.SD_rebin1.len
grep ">" B52_Sp13.BD.MM15.SD_rebin1.fa | sed 's/>//g' > B52_Sp13.BD.MM15.SD_rebin1.txt
/home/rprops/StrainMetaSim/scripts/AddLengths.pl B52_Sp13.BD.MM15.SD_rebin1.len < B52_Sp13.BD.MM15.SD_rebin1.txt > B52_Sp13.BD.MM15.SD_rebin1.tsv

rebinned_bins/final-refined_bins/B21_Fa13.BD.MLB.SN_rebin4.fa.tsv
/scratch/vdenef_fluxm/rprops/DESMAN/metaG/vizbin_rebin_anvio_v230/map/Fa13.BD.MM110.SD.mapped.bam
samtools view -hL rebinned_bins/final-refined_bins/B21_Fa13.BD.MLB.SN_rebin4.fa.tsv /scratch/vdenef_fluxm/rprops/DESMAN/metaG/vizbin_rebin_anvio_v230/map/Fa13.BD.MM110.SD.mapped.bam > test.sam


#!/bin/bash
set -e
for results_bin in `ls ./results_B*/`
do
	for sample in `ls ${results_bin}/*.tsv`
	do
		sed -n '3p;7p' < ./${sample} >> irep_${results_bin}.tsv
	done
cd -
done
