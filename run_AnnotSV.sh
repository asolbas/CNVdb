#!/bin/bash

help() {
	echo -e "The following programs must be installed: AnnotSV, bedtools, bcftools and tcl\n"
	echo -e	"The following environment variable must be created: \$SANNOTSV: AnnotSV installation directory \n"
	echo -e "Usage: run_AnnotSV input_file ouput_name\n"
	echo -e "Input file must have a header and at least three columns (tab-separated): chrm start end \n"
        }

while getopts "help:h" option; do
	case $option in
		h | help) #display help
			help
			exit;;
	esac
done


input=$1 #input file: must have at least three columns (chr start end) and header
output=$2 #name of output files

#header
#chr start end
awk '{ if (NR==1) print "#"$1"\t"$2"\t"$3"\t"$4 }' $input > $output.bed

#convert to bed format

awk '{if (NR!=1) print $1"\t"$2"\t"$3"\t"$4 }' ${input}  >> $output.bed

#run AnnotSV

$ANNOTSV/bin/AnnotSV -SVinputFile $output.bed -outputFile $output.annotated.tsv -svtBEDcol 4 -genomeBuild

#CREATE GENE AND ACMG CLASS LIST

awk '{ if (NR==1) { print $2"\t"$3"\t"$4"\t"$9"\t"$NF }}' $output.annotated.tsv >  ${output}_gene_cnv.tsv

awk '{if ($7=="split") {print $2"\t"$3"\t"$4"\t"$9"\t"$NF }} ' $output.annotated.tsv | sort -nk 1,1 -nk 2,2 > aux.temp 

awk '{ gsub(/full=/,"", $5); print $1"\t"$2"\t"$3"\t"$4"\t"$5 } ' aux.temp >> ${output}_gene_cnv.tsv


#delete temporal files

rm $output.bed 
rm aux.temp 

