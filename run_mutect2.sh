#!/bin/bash
#SBATCH --job-name=Mutect2
#SBATCH --output=jobName.out
#SBATCH --error=jobName.err
#SBATCH --time=16:59:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=40gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L


module load GATK

input_tumor=$1

fastq_name_path=$(echo $1 | cut -d"." -f1)
wes_name_path=$(echo $2 | cut -d"." -f1)
fastq_name=$(echo ${fastq_name_path##*/} | cut -d'_' -f2-)
WES_name=$(echo ${wes_name_path##*/} | cut -d'_' -f2-)


input_WES=$2

reference=/groups/umcg-weersma/tmp01/Iwan/datastorage/reference/hg38/Homo_sapiens_assembly38.fasta


uniq_run_id=$3

echo $input_tumor
echo $input_WES

gatk Mutect2 \
-R $reference \
-I $input_tumor \
-I $input_WES \
--normal $WES_name \
-O /groups/umcg-weersma/tmp01/Iwan/pipeline/output_mutect/${uniq_run_id}_${WES_name}_${fastq_name}.vcf.gz \
--dont-use-soft-clipped-bases true \
--germline-resource /groups/umcg-weersma/tmp01/Iwan/datastorage/reference/hg38/somatic-hg38_af-only-gnomad.hg38.vcf \
--panel-of-normals /groups/umcg-weersma/tmp01/Iwan/datastorage/reference/hg38/somatic-hg38_1000g_pon.hg38.vcf

echo "cleaning up"
rm $1
rm $2

#take a vcf input, annotate using funcotator default annotations, then output both vcf and maf files

#module load GATK

reference=/groups/umcg-weersma/tmp01/Iwan/datastorage/reference/hg38/Homo_sapiens_assembly38.fasta
gunzip /groups/umcg-weersma/tmp01/Iwan/pipeline/output_mutect/${uniq_run_id}_${WES_name}_${fastq_name}.vcf.gz

input_file=/groups/umcg-weersma/tmp01/Iwan/pipeline/output_mutect/${uniq_run_id}_${WES_name}_${fastq_name}.vcf
#hard coded paths only
input_name_path=$(echo ${input_file} | cut -d"." -f1)
input_name=$(echo ${input_name_path##*/})

gatk Funcotator \
     --variant ${input_file} \
     --reference ${reference} \
     --ref-version hg38 \
     --data-sources-path /groups/umcg-weersma/tmp01/Iwan/pipeline/funcotator/funcotator_dataSources.v1.6.20190124s \
     --output /groups/umcg-weersma/tmp01/Iwan/pipeline/output_funcotator/${input_name}_funcotated.vcf \
     --output-file-format VCF


if [ -f "/groups/umcg-weersma/tmp01/Iwan/pipeline/output_funcotator/${input_name}_funcotated.vcf" ]; then
        echo "$FILE exists."
else
        echo "$FILE run failed output not detected" >> /groups/umcg-weersma/tmp01/Iwan/pipeline/log.txt
fi
