import pandas as pd

# Load the original Galperin presence-absence table
presence_absence_file = 'Galperin2022_jb.00079-22-s0003.csv'
presence_absence_df = pd.read_csv(presence_absence_file)

# Load the updated gene lengths
gene_lengths_file = 'updated_galperin_gene_lengths.csv'
gene_lengths_df = pd.read_csv(gene_lengths_file)

# Create a dictionary of gene lengths
gene_lengths_dict = dict(zip(gene_lengths_df['gene'], gene_lengths_df['length']))

# Initialize lists to store results
genomes = []
present_nucleotides = []
lost_nucleotides = []

# Iterate over each genome in the presence-absence table
for idx, row in presence_absence_df.iterrows():
    genome = row['name']
    total_present = 0
    total_length = 0
    
    # Calculate the total length of all genes and the total length of present genes
    for gene, presence in row.items():
        if gene in ['name', 'spore_former']:
            continue
        if gene in gene_lengths_dict:
            gene_length = gene_lengths_dict[gene]
            total_length += gene_length if pd.notna(gene_length) else 0
            if presence == '1':
                total_present += gene_length if pd.notna(gene_length) else 0
    
    genomes.append(genome)
    present_nucleotides.append(total_present)
    lost_nucleotides.append(total_length - total_present)

# Create a DataFrame to store the results
nucleotide_counts_df = pd.DataFrame({
    'genome': genomes,
    'present_nucleotides': present_nucleotides,
    'lost_nucleotides': lost_nucleotides
})

# Save the results to a CSV file
output_file = 'nucleotide_presence_absence.csv'
nucleotide_counts_df.to_csv(output_file, index=False)

print(f"Nucleotide presence-absence table saved to {output_file}")
