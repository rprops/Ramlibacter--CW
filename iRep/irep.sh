#!/bin/bash
set -e
for bin in `ls ./bins/*.fa`
do
	bin_stub1=${bin##*/}
	bin_stub=${bin_stub1%-contigs.fa}
	mkdir ./results_${bin_stub}
	for sample_sam in `ls ./splitBam/${bin_stub}/*.sam`
		do
			stub1=${sample_sam##*/}
			stub=${stub1%.sam}
			echo $stub
			xvfb-run -e ./xvfb.${PBS_JOBID}.err -f ./xvfb.${PBS_JOBID}.auth -w 10 -s "-screen 0 1600x1200x24" iRep -f ${bin} -s ${sample_sam} -o ./results_${bin_stub}/iRep_${bin_stub}_${stub} -t 20 -ff --sort
	done
done

