#!/bin/bash

# references

bwa_index= # path to bwa index
ref= # path to reference genome
hybrid_ref= # path to the reference variants

# program paths

freebayes= # path to freebayes
hap_py= # path to hap.py

# path

illumina_raw= # path to illumina raw reads
illumina_trim= # path to illumina trimmed reads 
illumina_mapped= # path to illumina mapped reads
illumina_variant_calling= # path to the variant calling files
illumina_phasing= # path to the phase variant files
regions= # path to the bed file with the target region

sample= # sample name

chromosome= # chromosome where the gene is located
chrom_len= # length of the chromosome where the gene is located 
gene= # gene name
start= # gene start position
stop= # gene stop position

depth= # minimum sequencing depth

# activate python3 environment
source venv/bin/activate

# concat data
cat ${illumina_raw}${sample}*.fastq.gz > ${illumina_raw}${sample}.fastq.gz
 
cutadapt --max-n 0 -q 20,20 -g CTGTCTCTTATACACATCT -e 0.1 -O 5 -m 15 -o ${illumina_trim}${sample}_trim.fastq.gz ${illumina_raw}${sample}.fastq.gz

zcat ${illumina_trim}${sample}_trim.fastq.gz | paste - - - - | awk '{$1=$1"_"$2; print $1"\n"$3"\n"$4"\n"$5}' | gzip > ${illumina_trim}${sample}_trim_name.fastq.gz


# mapping data

round=1

map_output=${illumina_mapped}${sample}_round_${round}

bwa bwasw -b 7 -t 40 ${bwa_index} ${illumina_trim}${sample}_trim_name.fastq.gz > ${map_output}.sam 2> ${map_output}.out
samtools view -Sb ${map_output}.sam > ${map_output}.bam
samtools sort -@ 40 -n -o ${map_output}_sorted.bam ${map_output}.bam
samtools view ${map_output}_sorted.bam > ${map_output}_sorted.sam

python split_reads.py -n ${sample}_round -f ${illumina_mapped} -r ${round}

grep ^@ ${map_output}.sam > ${illumina_mapped}${sample}_header.sam

cat_command="cat ${illumina_mapped}${sample}_header.sam ${illumina_mapped}${sample}_round_${round}_split.sam"

while [ -s ${illumina_mapped}${sample}_round_${round}.fastq ]
do

fastq_input=${illumina_mapped}${sample}_round_${round}
((round++))
map_output=${illumina_mapped}${sample}_round_${round}

bwa bwasw -b 7 -t 40 ${bwa_index} ${fastq_input}.fastq > ${map_output}.sam 2> ${map_output}.out
samtools view -Sb ${map_output}.sam > ${map_output}.bam
samtools sort -@ 40 -n -o ${map_output}_sorted.bam ${map_output}.bam
samtools view ${map_output}_sorted.bam > ${map_output}_sorted.sam

python split_reads.py -n ${sample}_round -f ${illumina_mapped} -r ${round}

cat_command+=" ${illumina_mapped}${sample}_round_${round}_split.sam"

done

$cat_command > ${illumina_mapped}${sample}_all_reads.sam

samtools view -Sb ${illumina_mapped}${sample}_all_reads.sam > ${illumina_mapped}${sample}_all_reads.bam
samtools sort -@ 40 -o ${illumina_mapped}${sample}_all_reads_sorted.bam ${illumina_mapped}${sample}_all_reads.bam
samtools index ${illumina_mapped}${sample}_all_reads_sorted.bam

python primery_reads.py -n ${illumina_mapped}${sample}_all_reads_sorted

java -jar /opt/tools/picard.jar AddOrReplaceReadGroups I=${illumina_mapped}${sample}_all_reads_sorted_prim.bam O=${illumina_mapped}${sample}_all_reads_sorted_prim_RG.bam RGID=${sample} RGLB=${sample} RGPL=ILLUMINA RGPU=unit1 RGSM=${sample}
samtools index ${illumina_mapped}${sample}_all_reads_sorted_prim_RG.bam
gatk_4.1.4.0 HaplotypeCaller -R ${ref} -L ${chromosome} -I ${illumina_mapped}${sample}_all_reads_sorted_prim_RG.bam -O ${illumina_variant_calling}${sample}_prim_gatk.vcf.gz  --dont-use-soft-clipped-bases

deactivate

python2 ${hap_py} ${hybrid_ref} ${illumina_variant_calling}${sample}_prim_gatk.vcf.gz --threads 2 -r ${ref} -R ${regions}target_${gene}.bed  -o ${illumina_variant_calling}${sample}_prim_gatk_hap.out

source venv/bin/activate

bedtools coverage -a ${regions}target_${chromosome}.bed -b ${illumina_mapped}${sample}_all_reads_sorted_prim.bam -d > ${illumina_mapped}${sample}_coverage_${chromosome}.txt

python filter_vcf.py -m ${illumina_mapped}${gene}_ -v ${illumina_variant_calling}${sample}_prim_gatk -c ${chromosome} -d ${depth}

~/tools/bin/bgzip -f ${illumina_variant_calling}${sample}_prim_gatk_${chromosome}_filter_${depth}.vcf
~/tools/bin/tabix -f ${illumina_variant_calling}${sample}_prim_gatk_${chromosome}_filter_${depth}.vcf.gz

~/tools/bin/bgzip -f ${illumina_variant_calling}${sample}_prim_gatk_${chromosome}_ref.vcf
~/tools/bin/tabix -f ${illumina_variant_calling}${sample}_prim_gatk_${chromosome}_ref.vcf.gz

deactivate

python2 ${hap_py} ${illumina_variant_calling}${sample}_prim_gatk_${chromosome}_ref.vcf.gz ${illumina_variant_calling}${sample}_prim_gatk_${chromosome}_filter_${depth}.vcf.gz --threads 2 -r ${ref} -R regions/target_${gene}.bed  -o ${illumina_variant_calling}${sample}_prim_gatk_hap.out

source venv/bin/activate

python rename_reads.py -n ${sample}_all_reads_sorted_prim -f ${illumina_mapped}
samtools index ${illumina_mapped}${sample}_all_reads_sorted_prim_name.bam

zcat ${illumina_variant_calling}${sample}_prim_gatk_${chromosome}_filter_${depth}.vcf.gz | awk '! /\#/' | awk '{if(length($4) > length($5)) print $1"\t"($2-1)"\t"($2+length($4)-1); else print $1"\t"($2-1)"\t"($2+length($5)-1)}' > ${variant_path}${sample}_phasing.bed
samtools mpileup -B ${illumina_mapped}${sample}_all_reads_sorted_name.bam -l ${variant_path}${sample}_phasing.bed -q 1 -Q 20 -f ${ref} > phasing/illumina/${sample}_temp_CV.pileup
perl create_pgSnp.pl ${illumina_phasing}${sample}_temp_CV.pileup 0 > ${illumina_phasing}${sample}.pgsnp
samtools view -bq 1 ${illumina_variant_calling}${sample}_all_reads_sorted_name.bam ${chromosome}:${start}-${stop} > ${illumina_phasing}${sample}_seq_CV.bam
samtools view  ${illumina_phasing}${sample}_seq_CV.bam | cut -f 1,3,4,6,10 | sort -k3n >> ${illumina_phasing}${sample}_seq_CV.sam
python selectOverlappingReads.py ${illumina_phasing}${sample}_seq_CV.sam phasing/illumina/${sample}.pgsnp | sort > ${illumina_phasing}${sample}_total_snp_reads_CV.txt
perl linkSNP.pl ${illumina_phasing}${sample}_total_snp_reads_CV.txt > ${illumina_phasing}${sample}.hap
python haplotyper2.0.py ${illumina_phasing}${sample}.hap ${illumina_phasing}${sample}.pgsnp ${illumina_variant_calling}${sample}_prim_gatk_${chromosome}_filter_${depth}.vcf.gz ${chromosome} 0
Rscript -e "source('plot_linked_alleles_short_1.2.R'); png('phasing/illumina/${sample}.png', width=800, height=700); plot_linked_alleles(networkfile = 'phasing/illumina/${sample}_network.txt', locus = '${chromosome}:${start}-${stop}'); dev.off();"

deactivate

~/tools/bin/bgzip -f ${illumina_phasing}${sample}_variant.vcf
~/tools/bin/tabix -f ${illumina_phasing}${sample}_variant.vcf.gz

python2 ${hap_py} ${hybrid_ref} phasing/illumina/${sample}_variant.vcf.gz --threads 2 -r ${ref} -R ${regions}target_${gene}.bed  -o ${illumina_phasing}${sample}_variant_hap.out 

source venv/bin/activate

python Visualise.py --vcf ${illumina_phasing}${sample}_variant.vcf.gz -r ${hybrid_ref} --covered_ref ${illumina_variant_calling}${sample}_prim_gatk_${chromosome}_ref.vcf.gz --coverage ${illumina_mapped}${sample}_coverage_${chromosome}_depth_${depth}.bed --start ${start} --stop ${stop}  -c ${chromosome} -f ${illumina_phasing}${sample}_variant_plot -q False

deactivate

