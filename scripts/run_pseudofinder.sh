#!/bin/bash

#conda activate pseudofinder


# loop through genomes



for d in /Users/williamrshoemaker/GitHub/SporeCosts/data/genomes/prokka_annotations/*/ ; do
    echo "$d"
    genome_id=$(echo $d | rev | cut -f2 -d/ | rev)

    echo "$genome_id"

    pseudofinder.py annotate -g test_data/Mycobacterium_leprae_TN.gbff -db test_data/combined_mycobacteria.faa -op testAnnotate --diamond -ref test_data/Mycobacterium_tuberculosis_H37Rv.gbff


    #gbff_path="${d}"
done