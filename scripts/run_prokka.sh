#!/bin/bash

prokka_dir=/Users/williamrshoemaker/GitHub/SporeCosts/data/genome_analysis/prokka_annotations/

for fasta in /Users/williamrshoemaker/GitHub/SporeCosts/data/genome_analysis/genomes/*.fna; do
    genome_id=$(echo $fasta | rev | cut -f1 -d/ | rev | cut -f1 -d.)
    echo $genome_id
    
    if test -d ${prokka_dir}${genome_id}; then
        continue
    fi

    prokka --outdir ${prokka_dir}${genome_id} --force --compliant --rfam --prefix $genome_id $fasta
    # --gram pos
done





#prokka --outdir mydir --compliant --rfam --gram pos --prefix mygenome contigs.fa
