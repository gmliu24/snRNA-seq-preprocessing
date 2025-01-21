#!/bin/bash
#$ -N cellranger
#$ -cwd
#$ -pe omp 28
#$ -m bea
#$ -M gmliu@bu.edu

module load bcl2fastq/2.20                                                                                                                                   
module load cellranger/6.1.2

cd /projectnb/mccall/guangmeiliu/snrnaseq_gmliu/reference/

#Remove previous runs
#rm -rf /projectnb/mccall/guangmeiliu/snrnaseq_gmliu/reference/BDGP6.32.109


cellranger mkgtf /projectnb/mccall/guangmeiliu/snrnaseq_gmliu/reference/Drosophila_melanogaster.BDGP6.32.109.gtf /projectnb/mccall/guangmeiliu/snrnaseq_gmliu/reference/Drosophila_melanogaster.BDGP6.32.109.filtered.gtf \
   --attribute=gene_biotype:protein_coding \
   --attribute=gene_biotype:ncRNA 

cellranger mkref --genome=BDGP6.32.109 \
    --fasta=/projectnb/mccall/guangmeiliu/snrnaseq_gmliu/reference/Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa \
    --genes=/projectnb/mccall/guangmeiliu/snrnaseq_gmliu/reference/Drosophila_melanogaster.BDGP6.32.109.filtered.gtf \
    --nthreads=28
