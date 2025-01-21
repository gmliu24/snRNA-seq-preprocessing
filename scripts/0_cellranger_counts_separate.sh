#!/bin/bash
#$ -N cellrangercount
#$ -cwd
#$ -pe omp 28
#$ -m bea
#$ -M gmliu@bu.edu

module load bcl2fastq/2.20                                                                                                                                   
module load cellranger/7.2.0

cd /projectnb/mccall/guangmeiliu/snrnaseq_gmliu/analysis

cellranger count --id=$1 \
   --fastqs=/projectnb/mccall/guangmeiliu/snrnaseq_gmliu/data \
   --sample=$2 \
   --transcriptome=/projectnb/mccall/guangmeiliu/snrnaseq_gmliu/reference/BDGP6.32.109 \
   --include-introns true \
   --expect-cells 10000
