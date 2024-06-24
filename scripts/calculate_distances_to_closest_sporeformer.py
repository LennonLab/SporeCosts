import pandas as pd
from ete3 import Tree

# Load the Newick tree from the file
tree_file_path = 'RAxML_parsimonyTree.output_tree'
with open(tree_file_path, 'r') as file:
    newick_tree = file.read()

# Parse the Newick tree
tree = Tree(newick_tree)

# Read the taxon labels
taxon_labels_file = 'taxon_labels.csv'
taxon_labels_df = pd.read_csv(taxon_labels_file)

# Verify the content of taxon labels
print("Taxon Labels DataFrame:")
print(taxon_labels_df.head())

# Read the identifier mapping
mapping_df = pd.read_csv('identifier_mapping.txt', sep='\t', header=None, names=['Gene_ID', 'Genome_RibosomalGene'])

# Extract the Genome_ID from the Genome_RibosomalGene
mapping_df['Genome_ID'] = mapping_df['Genome_RibosomalGene'].apply(lambda x: "_".join(x.split("_")[:2]))

# Verify the content of the identifier mapping
print("Identifier Mapping DataFrame:")
print(mapping_df.head())

# Merge the taxon labels with the identifier mapping
taxon_labels_df['Genome_ID'] = taxon_labels_df['taxid'].apply(lambda x: f"genome_{x}")
merged_df = pd.merge(taxon_labels_df, mapping_df, left_on='Genome_ID', right_on='Genome_ID', how='inner')

# Verify the merged DataFrame
print("Merged DataFrame:")
print(merged_df.head())

# Separate spore-formers and lost taxa
spore_formers = merged_df[merged_df['spore_former'].str.strip() == 'spore-former']['Genome_ID'].unique().tolist()
lost_taxa = merged_df[merged_df['spore_former'].str.strip() == 'lost']['Genome_ID'].unique().tolist()

# Debugging statements
print("Spore Formers:", spore_formers)
print("Lost Taxa:", lost_taxa)

# Function to calculate the average distance for each genome
def calculate_average_distances(tree, taxa):
    genome_distances = {}
    for taxon in taxa:
        ribosomal_genes = mapping_df[mapping_df['Genome_ID'] == taxon]['Gene_ID'].tolist()
        total_distance = 0
        count = 0
        for gene in ribosomal_genes:
            try:
                distance = tree.get_distance(gene)
                total_distance += distance
                count += 1
            except:
                print(f"Error calculating distance for gene {gene} in taxon {taxon}")
        if count > 0:
            average_distance = total_distance / count
            genome_distances[taxon] = average_distance
        else:
            print(f"No valid distances for taxon {taxon}")
    return genome_distances

# Calculate the average distances for spore-formers and lost taxa
spore_former_distances = calculate_average_distances(tree, spore_formers)
lost_taxa_distances = calculate_average_distances(tree, lost_taxa)

# Function to calculate the distance from a lost taxon to the closest spore-former
def calculate_distance_to_closest_sporeformer(lost_taxa_distances, spore_former_distances):
    distances = []
    for lost_taxon, lost_distance in lost_taxa_distances.items():
        closest_distance = float('inf')
        closest_spore_former = None
        for spore_former, spore_distance in spore_former_distances.items():
            distance = abs(lost_distance - spore_distance)
            if distance < closest_distance:
                closest_distance = distance
                closest_spore_former = spore_former
        distances.append((lost_taxon, closest_spore_former, closest_distance))
    return distances

# Calculate distances to closest spore-former for each lost taxon
distances_to_closest_sporeformer = calculate_distance_to_closest_sporeformer(lost_taxa_distances, spore_former_distances)

# Convert distances to a DataFrame
distances_df = pd.DataFrame(distances_to_closest_sporeformer, columns=['Lost_Taxon', 'Closest_Spore_Former', 'Distance'])

# Save the result to a CSV file
distances_df.to_csv('distances_to_closest_sporeformer.csv', index=False)

# Display the DataFrame
print("Distances to the closest spore-former have been calculated and saved to 'distances_to_closest_sporeformer.csv'.")
print(distances_df.head())

