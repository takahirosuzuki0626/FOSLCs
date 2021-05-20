#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 4

## $1: id
## $2: fastqs

cellranger count --id=$1 \
                 --transcriptome=$HOME/ProgramFiles/Cell_Ranger/refdata-cellranger-mm10-3.0.0 \
                 --fastqs=$2 \
                 --expect-cells=8000 \
                 --localcores=1
                 
