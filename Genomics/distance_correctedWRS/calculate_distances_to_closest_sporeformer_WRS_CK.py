import pandas as pd
from ete3 import Tree
from Bio import Phylo

# Load the bootstrap trees from the .bootstraps file
bootstrap_file_path = 'ribosomal_genes_aligned_concat.fasta.raxml.bestModel'
bootstrap_trees = []

with open(bootstrap_file_path, 'r') as file:
    for line in file:
        newick_tree = line.strip()
        bootstrap_trees.append(Tree(newick_tree))

# Print all leaf names in one of the bootstrap trees (for debugging)
print("Leaf names in one of the bootstrap trees:")
for leaf in bootstrap_trees[0].iter_leaves():
    print(leaf.name)

# Load the taxon labels
taxon_labels_file = 'taxon_labels.csv'
taxon_labels_df = pd.read_csv(taxon_labels_file)

# Separate spore formers and lost taxa
spore_formers = taxon_labels_df[taxon_labels_df['spore_former'] == 'spore-former']['taxid'].astype(str).tolist()
lost_taxa = taxon_labels_df[taxon_labels_df['spore_former'] == 'lost ']['taxid'].astype(str).tolist()

# Add prefix 'genome_' to taxon IDs for matching with the tree leaves
spore_formers = ['genome_' + taxon for taxon in spore_formers]
lost_taxa = ['genome_' + taxon for taxon in lost_taxa]

# Function to find the closest spore former for each lost taxon in all bootstrap trees
def find_closest_spore_former_in_bootstrap_trees(lost_taxa, spore_formers, bootstrap_trees):
    closest_spore_former_distances = []

    for lost in lost_taxa:
        min_distance = float('inf')
        closest_spore_former = None

        for tree in bootstrap_trees:
            for spore_former in spore_formers:
                try:
                    distance = tree.get_distance(lost, spore_former)
                    if distance < min_distance:
                        min_distance = distance
                        closest_spore_former = spore_former
                except:
                    pass

        if closest_spore_former is not None:
            closest_spore_former_distances.append((lost, closest_spore_former, min_distance))
        else:
            print(f"Lost taxon {lost} not found in the tree.")

    return closest_spore_former_distances

# Find the closest spore former for each lost taxon in bootstrap trees
closest_spore_former_distances = find_closest_spore_former_in_bootstrap_trees(lost_taxa, spore_formers, bootstrap_trees)

# Convert to DataFrame
closest_spore_former_df = pd.DataFrame(closest_spore_former_distances, columns=['Lost_Taxon', 'Closest_Spore_Former', 'Distance'])

# Save the result to a CSV file
closest_spore_former_df.to_csv('distances_to_closest_sporeformer_bootstrap.csv', index=False)

# Display the DataFrame
print("Distances to the closest spore-former have been calculated and saved to 'distances_to_closest_sporeformer_bootstrap.csv'.")
print(closest_spore_former_df.head())
