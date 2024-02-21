#!/usr/bin/env bash
RENEE_DIR=/data/sevillas2/PIPELINES/RENEE
FASTQDIR=/data/CCBR/projects/ccbr1324/rawdata/rna
OUTPUTDIR=/data/CCBR/projects/ccbr1324/renee_240214

if [[ ! -d $OUTPUTDIR ]]; then mkdir -p $OUTPUTDIR; fi

module load ccbrpipeliner
$RENEE_DIR/bin/renee run \
    --input $FASTQDIR/*fastq.gz \
    --output $OUTPUTDIR \
    --genome mm10_M25 \
    --mode slurm \
    --sif-cache /data/CCBR_Pipeliner/SIFS $1
