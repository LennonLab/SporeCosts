import os
import pandas as pd
from ete3 import Tree


data_directory = os.path.expanduser("~/GitHub/SporeCosts/data/")
scripts_directory = os.path.expanduser("~/GitHub/SporeCosts/scripts/")

tree_file_path = '%sribosomal_genes_aligned_concat.fasta.raxml.bestTree' % data_directory
# Load the ML tree
tree = Tree(tree_file_path)

# Read taxon labels
taxon_labels_file = '%staxon_labels.csv' % data_directory
taxon_labels_df = pd.read_csv(taxon_labels_file)


# Read the identifier mapping
mapping_df = pd.read_csv('%sidentifier_mapping.txt' % data_directory, sep='\t', header=None, names=['Gene_ID', 'Genome_RibosomalGene'])

# Extract the Genome_ID from the Genome_RibosomalGene
mapping_df['Genome_ID'] = mapping_df['Genome_RibosomalGene'].apply(lambda x: "_".join(x.split("_")[:2]))

# Merge the taxon labels with the identifier mapping
taxon_labels_df['Genome_ID'] = taxon_labels_df['taxid'].apply(lambda x: f"genome_{x}")
merged_df = pd.merge(taxon_labels_df, mapping_df, left_on='Genome_ID', right_on='Genome_ID', how='inner')

# Separate spore-formers and lost taxa
spore_formers = merged_df[merged_df['spore_former'].str.strip() == 'spore-former']['Genome_ID'].unique().tolist()
lost_taxa = merged_df[merged_df['spore_former'].str.strip() == 'lost']['Genome_ID'].unique().tolist()


# get intersection of lost_taxa and taxa that are in the tree
tree_leaf_names = list(tree.get_leaf_names())
tree_leaf_names = [str(x) for x in tree_leaf_names]
spore_formers_tree = list(set(tree_leaf_names) & set(spore_formers))
lost_taxa_tree = list(set(tree_leaf_names) & set(lost_taxa))


# loop through each taxa that lost spore formation
# identify the spore-forming taxon it lost spore formation from as the spore-forming taxon with the smallest phylogeneic distance
# Identify the internal node for those two taxa (root)
# Calculate the distance between non-spore-forming and internal node.


output_file = open('%sdistances_to_closest_sporeformer_bootstrap.csv' % data_directory, 'w')
output_file.write(','.join(['Lost_taxon', 'Closest_Spore_Former', 'Distance']))
output_file.write('\n')

for lost_taxa_tree_i in lost_taxa_tree:

    dist_to_spore_i = [tree.get_distance(lost_taxa_tree_i, spore_formers_tree_j) for spore_formers_tree_j in spore_formers_tree ]
    spore_formers_tree_min_i = spore_formers_tree[dist_to_spore_i.index(min(dist_to_spore_i))]
    spore_formers_mrca_i = tree.get_common_ancestor(lost_taxa_tree_i, spore_formers_tree_min_i)
    # spore_formers_mrca_i.get_topology_id()
    dist_to_mrca_i = tree.get_distance(lost_taxa_tree_i, spore_formers_mrca_i)

    output_file.write(','.join([lost_taxa_tree_i, spore_formers_tree_min_i, str(dist_to_mrca_i)]))
    output_file.write('\n')





output_file.close()


