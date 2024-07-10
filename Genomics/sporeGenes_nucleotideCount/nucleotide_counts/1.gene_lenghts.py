import pandas as pd

# Load the gene list from Galperin_spore_genes.csv
gene_list_file = 'Galperin_spore_genes.csv'
gene_list_df = pd.read_csv(gene_list_file, header=None, names=['gene'])

# Load the local gene lengths database from SubtiWiki
subtiwiki_db_file = 'subtiwiki.gene.export.2022-05-11.csv'
subtiwiki_db_df = pd.read_csv(subtiwiki_db_file)

# Function to normalize identifiers
def normalize_identifier(identifier):
    return identifier.strip().lower() if identifier else None

# Normalize gene names in the database
subtiwiki_db_df['gene'] = subtiwiki_db_df['gene'].apply(normalize_identifier)
subtiwiki_db_df['locus_tag'] = subtiwiki_db_df['locus_tag'].apply(normalize_identifier)
gene_lengths_dict = dict(zip(subtiwiki_db_df['gene'], subtiwiki_db_df['gene_length']))
locus_tag_lengths_dict = dict(zip(subtiwiki_db_df['locus_tag'], subtiwiki_db_df['gene_length']))

# Match gene names from the Galperin gene list with the lengths from the database
gene_lengths = []
for gene in gene_list_df['gene']:
    gene_synonyms = [normalize_identifier(name) for name in gene.replace('(', '/').replace(')', '').split('/')]
    found_length = None
    for synonym in gene_synonyms:
        if synonym in gene_lengths_dict:
            found_length = gene_lengths_dict[synonym]
            break
        elif synonym in locus_tag_lengths_dict:
            found_length = locus_tag_lengths_dict[synonym]
            break
    gene_lengths.append({'gene': gene, 'length': found_length if found_length else 0})

# Convert the results to a DataFrame and save to CSV
gene_lengths_df = pd.DataFrame(gene_lengths)
output_file = 'galperin_gene_lengths.csv'
gene_lengths_df.to_csv(output_file, index=False)

print(f"Gene lengths table saved to {output_file}")
