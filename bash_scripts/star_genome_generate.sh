#!/bin/bash
#PBS -q hotel
#PBS -N star_genome_generate
#PBS -l nodes=1:ppn=16
#PBS -l walltime=4:00:00
#PBS -o star_genome_generate.out
#PBS -e star_genome_generate.err
#PBS -M agfoote@health.ucsd.edu
#PBS -m ae

STAR --runThreadN 16 --runMode genomeGenerate --genomeDir ~/scratch/mouse_genome_mm10/star --genomeFastaFiles ~/scratch/mouse_genome_mm10/allchrom.fa --sjdbGTFfile ~/scratch/mouse_genome_mm10/gencode.vM23.annotation.gtf
