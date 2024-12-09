#!/bin/bash
set -Eeuxo pipefail
trap "echo ERR trap fired!" ERR 
# Author: Nicole Gay, Anna Scherbina 
# Updated: 7 May 2020 
# Script : encode_to_count_matrix.sh ${indir} ${srcdir} ${batch_file} ${cores} ${final_results_dir}
# Purpose: generate peak x sample read counts matrix from tagAlign and hammock files 
# Run pass_extract_from_gcp.sh before this 

#get a merged peak file       

#module load miniconda/3 # for python3
#module load bedtools 

############################################################
## USER-DEFINED VARIABLES 
#Examples for user defined input arguments
#indir=~/test_mnt/PASS/atac-seq/stanford #base path to the location of outputs from all batches 
#srcdir=~/motrpac-atac-seq-pipeline/src # directory with truncate_narrowpeak_200bp_summit.py
#batch_file=/home/araja7/motrpac-atac-seq-pipeline_code_dev/src/test_batch.txt #file containing the list of batches to merge the peak files
#example contents of batch file
#batch5_2020092
#final_results_dir=pass1b_atac_final
#3 worked on gcp , 8 crashes the system current gcp vm with 60gb ram
#cores=3 # number of cores allocated for parallelization 
# need ~25G per core. 10G was too low
indir=$1
srcdir=$2
batch_file=$3
cores=$4
final_results_dir=$5
mode=$6
############################################################


#make the same code usable for generating counts from single or multiple batches

outdir=${indir}/${final_results_dir}/merged_peaks
echo ${outdir}
mkdir -p ${outdir}


# intersect with tagalign files 
mkdir -p ${indir}/${final_results_dir}/counts_matrix
intersect_tag () {
	local tag=$1
	#change below
	local results_dir="/home/araja7/test_mnt/PASS/atac-seq/rn7/pass1ac-06/Output/test_split_by_tissue/lung" #assumes the results folder has the word final
	#below loction was used for pass1ac
	#local results_dir=$(ls|grep "pass1c_merged_counts_v2")
	echo "results dir is" ${results_dir}
	echo "tag is" ${tag}		
	local viallabel=$(basename $tag | sed "s/_.*//")
	echo ${viallabel} > ${results_dir}/counts_matrix/counts.${viallabel}.txt
	bedtools coverage -nonamecheck -counts -a ${results_dir}/merged_peaks/overlap.optimal_peak.narrowPeak.200.sorted.merged.bed -b ${tag} | cut -f4 >> ${results_dir}/counts_matrix/counts.${viallabel}.txt
}
export -f intersect_tag

for i in `cat ${batch_file}`;do
	#tag_align=$(ls ${indir}/$i/Output/final/tagalign/*68??_R1*tagAlign.gz)
	#modify the below to match with the tissue
	tag_align=$(ls ${indir}/$i/Output/final/tagalign/*66??_R1*tagAlign.gz)
	echo ${tag_align}
	echo ${final_results_dir}
	parallel --verbose --progress --bar --jobs ${cores} intersect_tag ::: $(echo ${tag_align})
done

echo "Success generating sample level counts matrix"

echo -e $'chrom\tstart\tend' > ${indir}/${final_results_dir}/index
cat ${outdir}/overlap.optimal_peak.narrowPeak.200.sorted.merged.bed >> ${indir}/${final_results_dir}/index

#split the results counts matrix by tissue

cd ${indir}/${final_results_dir}/counts_matrix
ls *|awk -F "." '{print $2}'|awk '{print substr($1,8,2)}'|cut -f1|sort|uniq >>${indir}/${final_results_dir}/tmp_tids.txt
for i in `cat ${indir}/${final_results_dir}/tmp_tids.txt`;do
	paste ${indir}/${final_results_dir}/index counts.*$i??.txt >${indir}/${final_results_dir}/t$i.atac.counts.txt
	gzip ${indir}/${final_results_dir}/t$i.atac.counts.txt
done
rm ${indir}/${final_results_dir}/tmp_tids.txt
rm ${indir}/${final_results_dir}/index

echo "Success generating merged tissue level counts matrix"
