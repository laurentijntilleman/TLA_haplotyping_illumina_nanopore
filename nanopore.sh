#!/bin/bash

# parameters

threads=40

# references

minimap2_index= # path to minimap2 index
ref= # path to reference genome

# path

nanopore_basecall= # path to raw nanopore fastq files
nanopore_trim= # path to trimed fastq files
nanopore_mapped= # path to mapped data
nanopore_variant_calling= # path to files with the called variants
nanopore_phase= # path to phased variants
variant_path= # path to the variants file
regions= # path to the bed file with the target region

# variables

depth= # minimal sequencing depth
chromosome= # chromosome where the gene is located
gene= # gene name
start= # gene start position
stop= # gene stop position
min_variants_phase_block= # minimal variants per phasing block


# programs

hap_py= # path to hap.py


source venv/bin/activate

# concat data

cat ${nanopore_basecall}minion/fastq_pass/pass/* > ${nanopore_basecall}nanopore_all_pass.fastq

cat ${nanopore_basecall}nanopore_all_pass.fastq | paste - - - - | awk '{$1=$1"_"$2; print $1" "$3" "$4" "$5" "$6"\n"$7"\n"$8"\n"$9}' > ${nanopore_trim}nanopore_all_pass_name.fastq

# mapping data

round=1

map_output=${nanopore_mapped}round_${round}

minimap2 -t ${threads} -a -o ${map_output}.sam ${minimap2_index} ${nanopore_trim}nanopore_all_pass_name.fastq 2>> ${map_output}.out
samtools view -Sb ${map_output}.sam > ${map_output}.bam
samtools sort -@ ${threads} -n -o ${map_output}_sorted.bam ${map_output}.bam
samtools view ${map_output}_sorted.bam > ${map_output}_sorted.sam

python split_reads.py -n round -f ${nanopore_mapped} -r ${round}

grep ^@ ${map_output}.sam > ${nanopore_mapped}header.sam

cat_command="cat ${nanopore_mapped}header.sam ${nanopore_mapped}round_${round}_split.sam"

while [ -s ${nanopore_mapped}round_${round}.fastq ]
do

fastq_input=${nanopore_mapped}round_${round}
((round++))
map_output=${nanopore_mapped}round_${round}

minimap2 -t ${threads} -a -o ${map_output}.sam ${minimap2_index} ${fastq_input}.fastq 2>> ${map_output}.out
samtools view -Sb ${map_output}.sam > ${map_output}.bam
samtools sort -@ ${threads} -n -o ${map_output}_sorted.bam ${map_output}.bam
samtools view ${map_output}_sorted.bam > ${map_output}_sorted.sam

python split_reads.py -n round -f ${nanopore_mapped} -r ${round}

cat_command+=" ${nanopore_mapped}round_${round}_split.sam"

done

$cat_command > ${nanopore_mapped}all_reads.sam

samtools view -Sb ${nanopore_mapped}all_reads.sam > ${nanopore_mapped}all_reads.bam
samtools sort -@ ${threads} -o ${nanopore_mapped}all_reads_sorted.bam ${nanopore_mapped}all_reads.bam
samtools index ${nanopore_mapped}all_reads_sorted.bam

python rename_reads.py -n all_reads_sorted -f ${nanopore_mapped}

samtools index ${nanopore_mapped}all_reads_sorted_name.bam

deactivate

source ~/miniconda3/etc/profile.d/conda.sh
conda activate clair3

run_clair3.sh --bam_fn=${nanopore_mapped}all_reads_sorted.bam --ref_fn=${ref} --threads=${threads} --platform="ont" --model_path=~/miniconda3/envs/clair3/bin/models/ont_guppy5 --output=${nanopore_variant_calling}

conda deactivate

bedtools coverage -a ${regions}target_${chromosome}.bed -b ${nanopore_mapped}all_reads_sorted.bam -d > ${nanopore_mapped}${gene}_coverage_${chromosome}.txt

source venv/bin/activate

python filter_vcf.py -m ${nanopore_mapped}${gene}_ -v ${variant_path} -c ${chromosome} -d ${depth}

whatshap phase -o ${nanopore_phase}${gene}_${chromosome}_filter_${depth}_phased.vcf --include-homozygous --chromosome ${chromosome} --indels --distrust-genotypes --ignore-read-groups --reference=${ref} ${variant_path}_${chromosome}_filter_${depth}.vcf ${nanopore_mapped}all_reads_sorted.bam

~/tools/bin/bgzip -f ${nanopore_phase}${gene}_${chromosome}_filter_${depth}_phased.vcf
~/tools/bin/tabix -f ${nanopore_phase}${gene}_${chromosome}_filter_${depth}_phased.vcf.gz
~/tools/bin/bgzip -f ${variant_path}_${chromosome}_ref.vcf
~/tools/bin/tabix -f ${variant_path}_${chromosome}_ref.vcf.gz

deactivate

python2 ${hap_py} ${variant_path}_${chromosome}_ref.vcf.gz ${nanopore_phase}${gene}_${chromosome}_filter_${depth}_phased.vcf.gz --threads 2 -r ${ref} -R ${regions}target_${gene}.bed  -o ${variant_path}_${chromosome}_hap.out

whatshap haplotag -o ${nanopore_phase}${gene}_${chromosome}_filter_${depth}_phased.bam --ignore-read-groups --reference=${ref} ${nanopore_phase}${gene}_${chromosome}_filter_${depth}_phased.vcf.gz ${nanopore_mapped}all_reads_sorted.bam

samtools index ${nanopore_phase}${gene}_${chromosome}_filter_${depth}_phased.bam

python haplotyper3.0.py ${nanopore_phase}${gene}_${chromosome}_filter_${depth}_phased.bam ${nanopore_phase}${gene}_${chromosome}_filter_${depth}_phased.vcf.gz  ${nanopore_phase}${gene}_ ${chromosome} ${gene_start} ${gene_stop} ${min_variants_phase_block} 0 True

~/tools/bin/bgzip -f ${nanopore_phase}${gene}_${chromosome}_phased_haplo.vcf
~/tools/bin/tabix -f ${nanopore_phase}${gene}_${chromosome}_phased_haplo.vcf.gz

deactivate

python2 ${hap_py} ${ref} ${nanopore_phase}${gene}_${chromosome}_phased_haplo.vcf.gz --threads 2 -r ${ref} -R ${regions}target_${gene}.bed  -o ${nanopore_phase}${gene}_${chromosome}_phased_haplo_hap.out


source venv/bin/activate

python Visualise.py --vcf ${nanopore_phase}${gene}_${chromosome}_phased_haplo.vcf.gz -r ref/hg38.hybrid.vcf.gz --covered_ref ${variant_path}_${chromosome}_ref.vcf.gz --coverage ${nanopore_mapped}${gene}_coverage_${chromosome}_depth_${depth}.bed --start ${gene_start} --stop ${gene_stop} -c ${chromosome} -f ${nanopore_phase}${gene}_variant_plot -q False --haplo

deactivate
