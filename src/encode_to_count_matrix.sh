#!/bin/bash
# 
# Author: Nicole Gay, Anna Scherbina 
# Updated: 7 May 2020 
# Script : encode_to_count_matrix.sh 
# Purpose: generate peak x sample read counts matrix from tagAlign and hammock files 
# Run pass_extract_from_gcp.sh before this 

#get a merged peak file       

module load miniconda/3
module load bedtools

indir=/projects/motrpac/PASS1A/ATAC/NOVASEQ_BATCH2/outputs
outdir=${indir}/merged_peaks
srcdir=/oak/stanford/groups/smontgom/nicolerg/src/MOTRPAC/ATAC/PASS
cores=10 # 25G per core. 10G was too low 

mkdir -p ${outdir}
cd ${indir}

# concatenate peaks (narrowpeak.gz)
cat $(find -path "./peak/*narrowPeak.gz") > ${outdir}/overlap.optimal_peak.narrowPeak.bed.gz

#truncate peaks to 200 bp around summit
python ${srcdir}/truncate_narrowpeak_200bp_summit.py --infile ${outdir}/overlap.optimal_peak.narrowPeak.bed.gz --outfile ${outdir}/overlap.optimal_peak.narrowPeak.200.bed.gz

# sort and merge peaks --> master peak file 
zcat ${outdir}/overlap.optimal_peak.narrowPeak.200.bed.gz | bedtools sort | bedtools merge > ${outdir}/overlap.optimal_peak.narrowPeak.200.sorted.merged.bed

# intersect with tagalign files 
mkdir -p counts_matrix
intersect_tag () {
	local tag=$1
	local viallabel=$(basename $tag | sed "s/_.*//")
	echo ${viallabel} > counts_matrix/counts.${viallabel}.txt
	bedtools coverage -nonamecheck -counts -a merged_peaks/overlap.optimal_peak.narrowPeak.200.sorted.merged.bed -b ${tag} | cut -f4 >> counts_matrix/counts.${viallabel}.txt
}
export -f intersect_tag

tag_align=$(find -path "./tagalign/*tagAlign.gz")
parallel --verbose --jobs ${cores} intersect_tag ::: $(echo ${tag_align}) 

echo -e $'chrom\tstart\tend' > index
cat merged_peaks/overlap.optimal_peak.narrowPeak.200.sorted.merged.bed >> index
paste index counts_matrix/counts* > batch_20200318.atac.counts.txt
gzip batch_20200318.atac.counts.txt
rm index
