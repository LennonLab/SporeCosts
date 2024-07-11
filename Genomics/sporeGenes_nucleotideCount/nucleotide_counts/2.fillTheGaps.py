import pandas as pd

# Load the gene list and gene lengths from the previous step
gene_list_file = 'Galperin_spore_genes.csv'
gene_lengths_file = 'galperin_gene_lengths.csv'
gene_list_df = pd.read_csv(gene_list_file, header=None, names=['gene'])
gene_lengths_df = pd.read_csv(gene_lengths_file)

# Find genes with length 0 (not found in the previous step)
remaining_genes = gene_lengths_df[gene_lengths_df['length'] == 0]['gene'].tolist()

# Load the map file
map_file = 'map_file.csv'
map_df = pd.read_csv(map_file)

# Function to normalize identifiers
def normalize_identifier(identifier):
    return str(identifier).strip().lower() if pd.notna(identifier) else None

# Normalize identifiers in the map file
map_df['geneP'] = map_df['geneP'].apply(normalize_identifier)
map_df['locus_tag'] = map_df['locus_tag'].apply(normalize_identifier)
map_df['locus_tag2'] = map_df['locus_tag2'].apply(normalize_identifier)
map_df['gene1'] = map_df['gene1'].apply(normalize_identifier)
map_df['gene3'] = map_df['gene3'].apply(normalize_identifier)
map_df['gene4'] = map_df['gene4'].apply(normalize_identifier)
map_df['gene5'] = map_df['gene5'].apply(normalize_identifier)

# Create a dictionary to map synonyms and locus tags to gene lengths
map_dict = {}
for idx, row in map_df.iterrows():
    for col in ['geneP', 'locus_tag', 'locus_tag2', 'gene1', 'gene3', 'gene4', 'gene5']:
        if pd.notna(row[col]):
            map_dict[normalize_identifier(row[col])] = row['length'] if 'length' in row else 0

# Match remaining genes with the map file
for gene in remaining_genes:
    gene_synonyms = [normalize_identifier(name) for name in gene.replace('(', '/').replace(')', '').replace('/', ',').split(',')]
    found_length = None
    for synonym in gene_synonyms:
        if synonym in map_dict:
            found_length = map_dict[synonym]
            break
    gene_lengths_df.loc[gene_lengths_df['gene'] == gene, 'length'] = found_length if found_length else 0

# Save the updated gene lengths to CSV
updated_output_file = 'updated_galperin_gene_lengths.csv'
gene_lengths_df.to_csv(updated_output_file, index=False)

print(f"Updated gene lengths table saved to {updated_output_file}")
