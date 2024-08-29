#!/usr/bin/env python3

import os
import pandas as pd
from Bio import SeqIO
import logging

# Setup logging
logging.basicConfig(filename='sporulation_gene_analysis.log', level=logging.INFO, format='%(asctime)s - %(message)s')
logging.info("Script started")

# Define paths
sporulation_genes_file = 'sporulation_genes.fasta'  # Sporulation genes file
annotation_dir = 'prokka_annotations'  # Directory containing annotation files
map_file = 'map_file.csv'  # Map file with locus tags and synonyms
output_dir = 'output'

# Check directories
logging.info(f"Sporulation genes file: {sporulation_genes_file}")
logging.info(f"Annotation directory: {annotation_dir}")
logging.info(f"Map file: {map_file}")
logging.info(f"Output directory: {output_dir}")

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Load the map file
map_df = pd.read_csv(map_file)
gene_to_synonyms = {}
for _, row in map_df.iterrows():
    gene_name = row['geneP']
    synonyms = row[['locus_tag', 'locus_tag2', 'gene1', 'gene3', 'gene4', 'gene5']].dropna().tolist()
    gene_to_synonyms[gene_name] = synonyms

# Load sporulation genes and calculate total nucleotides
sporulation_genes = list(SeqIO.parse(sporulation_genes_file, "fasta"))
total_sporulation_nucleotides = sum(len(gene.seq) for gene in sporulation_genes)
logging.info(f"Total nucleotides in sporulation genes: {total_sporulation_nucleotides}")

# Prepare a list of all sporulation gene names and their sizes
sporulation_gene_names = [gene.id for gene in sporulation_genes]
gene_sizes = {gene.id: len(gene.seq) for gene in sporulation_genes}

# Initialize results
results = []
presence_absence_table = []

# Iterate over each annotation file
for genome_id in os.listdir(annotation_dir):
    annotation_folder = os.path.join(annotation_dir, genome_id)
    annotation_file = os.path.join(annotation_folder, f"{genome_id}.gbk")
    
    if not os.path.exists(annotation_file):
        logging.info(f"Annotation file not found for genome {genome_id}. Skipping...")
        continue
    
    # Load the annotation file
    present_nucleotides = 0
    gene_presence = {gene_name: 0 for gene_name in sporulation_gene_names}
    try:
        for record in SeqIO.parse(annotation_file, "genbank"):
            for feature in record.features:
                if feature.type == "CDS":
                    gene_name = feature.qualifiers.get('gene', [None])[0]
                    locus_tag = feature.qualifiers.get('locus_tag', [None])[0]
                    if locus_tag and locus_tag.startswith('BSU') and '_' not in locus_tag:
                        locus_tag = f'BSU_{locus_tag[3:]}'
                    all_names = [gene_name, locus_tag] if gene_name and locus_tag else [gene_name or locus_tag]
                    if gene_name:
                        all_names.extend(gene_to_synonyms.get(gene_name, []))
                    for sporulation_gene in sporulation_genes:
                        if sporulation_gene.id in all_names:
                            gene_length = len(feature.location)
                            present_nucleotides += gene_length
                            gene_presence[sporulation_gene.id] = 1
    except Exception as e:
        logging.error(f"Error processing annotation file {annotation_file}: {e}")
        continue
    
    lost_nucleotides = total_sporulation_nucleotides - present_nucleotides
    gene_presence["present_nucleotides"] = present_nucleotides
    gene_presence["lost_nucleotides"] = lost_nucleotides
    gene_presence["genome"] = genome_id
    presence_absence_table.append(gene_presence)

# Add a row for gene sizes
gene_sizes["genome"] = "Gene Sizes"
gene_sizes["present_nucleotides"] = total_sporulation_nucleotides
gene_sizes["lost_nucleotides"] = 0
presence_absence_table.insert(0, gene_sizes)

# Convert results to DataFrame
presence_absence_df = pd.DataFrame(presence_absence_table)

# Save results to CSV
presence_absence_df.to_csv(os.path.join(output_dir, "sporulation_gene_presence_absence.csv"), index=False)
logging.info("Analysis complete. Results saved to 'sporulation_gene_presence_absence.csv'.")
