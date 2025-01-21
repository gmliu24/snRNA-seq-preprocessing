#!/bin/bash
#$ -N fastqc
#$ -cwd
#$ -pe omp 10
#$ -m bea
#$ -M gmliu@bu.edu

module load fastqc

fastqfiledir=$1
outputdir=$2

cd $fastqfiledir
for i in `ls $fastqfiledir`;
    do fastqc $fastqfiledir/$i --outdir $outputdir;
done 


#echo fastqc $fastqfile --outdir $outputdir
#fastqc $fastqfile --outdir $outputdir
