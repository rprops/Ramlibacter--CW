#!/bin/bash
set -e
for bin in `ls ./bins/*.tsv`
do
	bin_stub1=${bin##*/}
	bin_stub=${bin_stub1%-contigs.tsv}
	echo $bin_stub
	mkdir ./splitBam/${bin_stub}
	for sample_sam in `ls /scratch/vdenef_fluxm/rprops/metaG_SCK/data/coassembly/anvio_v2.2.2/map/*.sam`
		do
			stub1=${sample_sam##*/}
			stub=${stub1%.sam}
			echo $stub
			echo $bin
			echo $sample_sam
			samtools view -hL $bin $sample_sam > ./splitBam/${bin_stub}/${stub}.sam
	done
done

