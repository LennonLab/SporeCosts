import os
from Bio import SeqIO
import pandas as pd

# Directory containing Prokka annotation files
annotation_dir = 'prokka_annotations'
phylo_names_file = 'phylo_names.csv'
output_file = 'genome_sizes_and_noncoding_nucleotides.csv'

# Load the mapping file
phylo_names_df = pd.read_csv(phylo_names_file)

# Initialize lists to store results
results = []

# Function to normalize identifiers
def normalize_identifier(identifier):
    return identifier.strip().lower() if identifier else None

# Iterate over each annotation folder
for index, row in phylo_names_df.iterrows():
    taxid = str(row['Taxid'])
    genome_folder = f'genome_{taxid}'
    annotation_file = os.path.join(annotation_dir, genome_folder, f"{genome_folder}.gbk")
    
    if not os.path.exists(annotation_file):
        print(f"Annotation file not found for genome {genome_folder}. Skipping...")
        continue

    genome_size = 0
    coding_nucleotides = 0
    
    try:
        for record in SeqIO.parse(annotation_file, "genbank"):
            genome_size += len(record.seq)
            for feature in record.features:
                if feature.type == "CDS":
                    coding_nucleotides += len(feature.location)
    except Exception as e:
        print(f"Error processing annotation file {annotation_file} for genome {genome_folder}: {e}")
    
    noncoding_nucleotides = genome_size - coding_nucleotides
    results.append({'Taxid': taxid, 'Genome_Size': genome_size, 'Noncoding_Nucleotides': noncoding_nucleotides})

# Convert results to DataFrame
results_df = pd.DataFrame(results)

# Save results to CSV
results_df.to_csv(output_file, index=False)

print(f"Genome sizes and noncoding nucleotides saved to {output_file}")
