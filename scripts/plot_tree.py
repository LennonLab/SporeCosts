from __future__ import division

import numpy
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import config
import os

import math
import sys

import ete3



data_directory = os.path.expanduser("~/GitHub/SporeCosts/data/")
scripts_directory = os.path.expanduser("~/GitHub/SporeCosts/scripts/")

#consumption_blocks = utils.identify_consumption_blocks()

#block_color_dict = {0: 'red', 1:'blue', 2:'green'}
#asv_color_dict = {}

# genome_562

#n_asv = 0
#for consumption_block_idx, consumption_block in enumerate(consumption_blocks):
#    for asv in consumption_block:
#        asv_color_dict[str(asv.name)] = block_color_dict[consumption_block_idx]
#        n_asv+=1


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

#print(lost_taxa)


tree_file_path = '%sribosomal_genes_aligned_concat.fasta.raxml.bestTree' % data_directory
#tree_file_path = '%scogs_tree/cogs_concat_alignment_outgroup.fasta.raxml.bestTree' % data_directory
#tree_file_path = '%scogs_tree/cogs_concat_alignment_outgroup.fasta.raxml.support' % data_directory

# Load the ML tree
tree = ete3.Tree(tree_file_path)

tree.show_branch_support = True

#networkx_graph = to_networkx(tree)
ts = ete3.TreeStyle()
ts.show_leaf_name = True

ts.show_branch_support = True

asvs = numpy.asarray(tree.get_leaf_names())


node2rootdist = {tree:0}
for node in tree.iter_descendants('preorder'):
    node2rootdist[node] = node.dist + node2rootdist[node.up]




spore_dist = []
non_spore_dist = []

for n in tree.traverse():
    nstyle = ete3.NodeStyle()

    #nstyle["fgcolor"] = "red"
    #nstyle["size"] = 15
    #nstyle['name'] = ''
    #print(vars(nstyle))

    if n.is_leaf() == True:
        
        # outgroup
        if n.name == 'genome_562':
            nstyle["fgcolor"] = matplotlib.colors.to_hex('k')

        elif str(n.name) in lost_taxa:
            # red
            non_spore_dist.append(node2rootdist[n])
            nstyle["fgcolor"] = '#d62728'

        else:
            # blue
            spore_dist.append(node2rootdist[n])
            nstyle["fgcolor"] = '#1f77b4'
        
        #if str(n.name) in asv_color_dict:
        #    nstyle["fgcolor"] = asv_color_dict[str(n.name)]
    else:
        nstyle['size'] = 0

    #else:
    #    nstyle["fgcolor"] = 'transparent'

    n.name = ''
    n.set_style(nstyle)


print(numpy.mean(spore_dist), numpy.mean(non_spore_dist))
print(numpy.std(spore_dist)/numpy.sqrt(len(spore_dist)), numpy.std(non_spore_dist)/numpy.sqrt(len(spore_dist)))



#fig_name = "%stree_cog.png" % config.analysis_directory
#tree.render(fig_name, w=100, units="mm", dpi=600, tree_style=ts)