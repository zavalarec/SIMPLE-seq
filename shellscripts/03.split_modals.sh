#!/bin/sh
s=$1

module load samtools
module load perl

perl perlscripts/02.split_modality.pl ${s}_mapped_sorted.bam
