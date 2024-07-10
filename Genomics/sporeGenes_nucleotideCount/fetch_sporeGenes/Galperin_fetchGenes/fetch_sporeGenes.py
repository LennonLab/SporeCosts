import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import os

# Load the gene list
gene_list_file = 'Galperin_spore_genes.csv'
with open(gene_list_file, 'r') as file:
    gene_list = [line.strip().lower() for line in file]

# Directory containing Prokka annotation files
annotation_dir = 'prokka_annotations'
output_dir = 'output'
os.makedirs(output_dir, exist_ok=True)

# Function to normalize identifiers
def normalize_identifier(identifier):
    return identifier.strip().lower() if identifier else None

# Function to extract gene sequences from a GenBank file
def extract_gene_sequences(genbank_file, gene_list):
    gene_sequences = []
    for record in SeqIO.parse(genbank_file, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                gene_name = normalize_identifier(feature.qualifiers.get('gene', [None])[0])
                locus_tag = normalize_identifier(feature.qualifiers.get('locus_tag', [None])[0])
                
                all_names = [gene_name, locus_tag] if gene_name and locus_tag else [gene_name or locus_tag]
                all_names = [normalize_identifier(name) for name in all_names if name]
                
                for name in all_names:
                    if name in gene_list:
                        seq_record = SeqRecord(
                            feature.extract(record.seq),
                            id=name,
                            description=f"Extracted from {genbank_file}"
                        )
                        gene_sequences.append(seq_record)
    return gene_sequences

# Extract sequences for Bacillus subtilis
bacillus_genbank_file = 'genomic.gbff'
bacillus_sequences = extract_gene_sequences(bacillus_genbank_file, gene_list)

# Save sequences to a FASTA file
fasta_file = os.path.join(output_dir, 'sporulation_genes.fasta')
SeqIO.write(bacillus_sequences, fasta_file, 'fasta')

# Print the summary
print(f"Total genes in the list: {len(gene_list)}")
print(f"Total gene sequences found: {len(bacillus_sequences)}")
print(f"FASTA file saved to {fasta_file}")
