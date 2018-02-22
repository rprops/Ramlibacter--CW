#!/bin/bash
set -e
for results_file in `ls results*/iRep*tsv`
do
	echo ${results_file}
	sed -n '3p;7p' < ${results_file} >> results_irep.tsv
done

