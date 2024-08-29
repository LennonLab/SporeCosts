import os
import pandas as pd
#from ete3 import Tree
from Bio import SeqIO
import config

import numpy

import matplotlib.pyplot as plt


import sys
import pickle
import utils

# conda activate pseudofinder

n_fna_characters = 80

data_directory = os.path.expanduser("~/GitHub/SporeCosts/data/")
scripts_directory = os.path.expanduser("~/GitHub/SporeCosts/scripts/")
directory = os.path.expanduser("~/GitHub/SporeCosts/")


spore_genome_statistics_dict_path = '%sspore_genome_statistics_dict.pickle' % data_directory


def get_spore_forming_cogs():

    # list of spore-forming COGs
    cogs_path = '%sCOGs_Galperin/speciesMatrix.csv' % directory
    cogs_file = open(cogs_path, 'r')
    cogs_header = cogs_file.readline()
    cogs_all = [line.strip().split(',')[1] for line in cogs_file]
    cogs_file.close()

    return cogs_all


def build_pseudofinder_database_per_genome():

    cogs_all = get_spore_forming_cogs()

    mrca_file = open('%sdistances_to_closest_sporeformer_bootstrap.csv' % data_directory,  'r')
    header = mrca_file.readline()

    for line in mrca_file:
        line = line.strip().split(',')

        non_spore_genome, mrca_spore_genome = line[0], line[1]

        #if non_spore_genome != 'genome_176280':
        #    continue

        # load gbk file
        non_spore_gbk_path = '%sgenome_analysis/prokka_annotations/%s/%s.gbk' % (data_directory, non_spore_genome, non_spore_genome)
        spore_gbk_path = '%sgenome_analysis/prokka_annotations/%s/%s.gbk' % (data_directory, mrca_spore_genome, mrca_spore_genome)

        spore_cog_faa_path = '%sgenome_analysis/pseudofinder_database/%s.faa' % (data_directory, mrca_spore_genome)
        spore_cog_faa_file = open(spore_cog_faa_path, 'w')

        # build .faa of MRCA
        for record in SeqIO.parse(spore_gbk_path, "genbank"):
            for feature in record.features:
                if feature.type == "CDS":
                    cog_name = feature.qualifiers.get('db_xref', [None])[0]
                    
                    if cog_name == None:
                        continue

                    cog_name = cog_name.split(':')[1]

                    if cog_name not in cogs_all:
                        continue
                    
                    aa_seq = feature.qualifiers['translation'][0]

                    spore_cog_faa_file.write('>%s\n' % cog_name)

                    for i in range(0, len(aa_seq), n_fna_characters):
                        sequence_i = aa_seq[i : i + n_fna_characters]
                        spore_cog_faa_file.write('%s\n' % sequence_i)
    
                    spore_cog_faa_file.write('\n')


        spore_cog_faa_file.close()

        # make database
        os.system('makeblastdb -dbtype prot -in %s -out %s' % (spore_cog_faa_path, spore_cog_faa_path))
        #makeblastdb -dbtype prot -in /Users/williamrshoemaker/GitHub/SporeCosts/data/genome_analysis/pseudofinder_database/genome_241244.faa -out /Users/williamrshoemaker/GitHub/SporeCosts/data/genome_analysis/pseudofinder_database/genome_241244.faa
        #makeblastdb -dbtype prot -in ofinder_database/genome_241244.faa -out /Users/williamrshoemaker/GitHub/SporeCosts/data/genome_analysis/pseudofinder_database/genome_241244.faa

        # run pseudofinder
        pseudofinder_output_path = '%sgenome_analysis/pseudofinder_output/%s' % (data_directory, non_spore_genome)
        os.system('pseudofinder.py annotate -g %s -db %s -op %s' % (non_spore_gbk_path, spore_cog_faa_path, pseudofinder_output_path))

        # load mrca_spore_genome genome




def build_spore_genome_statistics_dict():

    spore_forming_cogs = get_spore_forming_cogs()

    # list non-spore-forming lineages
    non_spore_former_genomes = []
    mrca_spore_former_genomes = []
    mrca_distances = []
    mrca_file = open('%sdistances_to_closest_sporeformer_bootstrap.csv' % data_directory,  'r')
    header = mrca_file.readline()
    
    for line in mrca_file:
        line = line.strip().split(',')
        non_spore_former_genomes.append(line[0])
        mrca_spore_former_genomes.append(line[1])
        mrca_distances.append(float(line[2]))
    mrca_file.close()

    # get number of pseudogenized spore-forming genes.
    non_spore_former_genomes_i = non_spore_former_genomes[0]
    non_spore_genome_stats_dict = {}

    for non_spore_former_genomes_i_idx, non_spore_former_genomes_i in enumerate(non_spore_former_genomes):

        pseudofiner_fasta_path = '%sgenome_analysis/pseudofinder_output/%s_pseudos.fasta' % (data_directory, non_spore_former_genomes_i)

        pseudofiner_fasta = utils.classFASTA(pseudofiner_fasta_path).readFASTA()
        n_pseudogenes = len(pseudofiner_fasta)

        mrca_spore_former_fasta = utils.classFASTA('%sgenome_analysis/pseudofinder_database/%s.faa' % (config.data_directory, mrca_spore_former_genomes[non_spore_former_genomes_i_idx])).readFASTA()
        mrca_spore_former_cogs = [s[0] for s in mrca_spore_former_fasta]

        n_spore_forming_cogs_nucleotides = 0
        n_cds_nucleotides = 0
        #n_genome_nucleotides = 0
        spore_forming_cogs_i = []
        # spore-forming COGs

        non_spore_gbk_path = '%sgenome_analysis/prokka_annotations/%s/%s.gbk' % (data_directory, non_spore_former_genomes_i, non_spore_former_genomes_i)

        for record in SeqIO.parse(non_spore_gbk_path, "genbank"):
            for feature in record.features:
                if feature.type == "CDS":
                    cog_name = feature.qualifiers.get('db_xref', [None])[0]

                    n_cds_nucleotides += len(str(feature.location.extract(record).seq))
                    
                    if cog_name == None:
                        continue

                    cog_name = cog_name.split(':')[1]

                    if cog_name not in spore_forming_cogs:
                        continue

                    #n_spore_forming_cogs += 1
                    spore_forming_cogs_i.append(cog_name)
                    n_spore_forming_cogs_nucleotides += len(str(feature.location.extract(record).seq))


        # number of spore-forming genes lost from MRCA
        n_lost_spore_forming = len(set(mrca_spore_former_cogs) - set(spore_forming_cogs_i))

        non_spore_genome_stats_dict[non_spore_former_genomes_i] = {}
        non_spore_genome_stats_dict[non_spore_former_genomes_i]['n_pseudogenes'] = n_pseudogenes
        non_spore_genome_stats_dict[non_spore_former_genomes_i]['n_spore_forming_cogs'] = len(spore_forming_cogs_i)
        non_spore_genome_stats_dict[non_spore_former_genomes_i]['n_spore_forming_cogs_nucleotides'] = n_spore_forming_cogs_nucleotides
        non_spore_genome_stats_dict[non_spore_former_genomes_i]['n_lost_spore_forming'] = n_lost_spore_forming
        non_spore_genome_stats_dict[non_spore_former_genomes_i]['mrca_distance'] = mrca_distances[non_spore_former_genomes_i_idx]
        non_spore_genome_stats_dict[non_spore_former_genomes_i]['mrca_spore_former_genome'] = mrca_spore_former_genomes[non_spore_former_genomes_i_idx]

    #return non_spore_genome_stats_dict
    #n_genome_nucleotides += len(str(record.seq))

    sys.stderr.write("Saving parameter dictionary...\n")
    with open(spore_genome_statistics_dict_path, 'wb') as outfile:
        pickle.dump(non_spore_genome_stats_dict, outfile, protocol=pickle.HIGHEST_PROTOCOL)
    sys.stderr.write("Done!\n")




def plot_spore_genome_statistics(mean_over_descendants=False):

    spore_genome_statistics_dict = pickle.load(open(spore_genome_statistics_dict_path, "rb"))

    genomes = list(spore_genome_statistics_dict.keys())
    n_pseudogenes = [spore_genome_statistics_dict[g]['n_pseudogenes'] for g in genomes]
    n_spore_forming_cogs = [spore_genome_statistics_dict[g]['n_spore_forming_cogs'] for g in genomes]
    n_spore_forming_cogs_nucleotides = [spore_genome_statistics_dict[g]['n_spore_forming_cogs_nucleotides'] for g in genomes]
    mrca_distance = [spore_genome_statistics_dict[g]['mrca_distance'] for g in genomes]
    mrca_spore_former_genome = [spore_genome_statistics_dict[g]['mrca_spore_former_genome'] for g in genomes]

    mrca_spore_former_genome_set = list(set(mrca_spore_former_genome))

    genomes = numpy.asarray(genomes)
    n_pseudogenes = numpy.asarray(n_pseudogenes)
    n_spore_forming_cogs = numpy.asarray(n_spore_forming_cogs)
    n_spore_forming_cogs_nucleotides = numpy.asarray(n_spore_forming_cogs_nucleotides)
    mrca_distance = numpy.asarray(mrca_distance)
    mrca_spore_former_genome = numpy.asarray(mrca_spore_former_genome)


    n_pseudogenes_to_plot = []
    n_spore_forming_cogs_to_plot = []
    n_spore_forming_cogs_nucleotides_to_plot = []
    mrca_distance_to_plot = []

    if mean_over_descendants == True:

        #for m in mrca_spore_former_genome_set:

        #    m_idx = (mrca_spore_former_genome == m)
        #    if sum(m_idx) == 1:

        n_pseudogenes = [numpy.mean(n_pseudogenes[mrca_spore_former_genome==g]) for g in mrca_spore_former_genome_set]
        n_spore_forming_cogs = [numpy.mean(n_spore_forming_cogs[mrca_spore_former_genome==g]) for g in mrca_spore_former_genome_set]
        n_spore_forming_cogs_nucleotides = [numpy.mean(n_spore_forming_cogs_nucleotides[mrca_spore_former_genome==g]) for g in mrca_spore_former_genome_set]
        mrca_distance = [numpy.mean(mrca_distance[mrca_spore_former_genome==g]) for g in mrca_spore_former_genome_set]

        print(n_pseudogenes)



    # get mean for each ancestor.


    fig = plt.figure(figsize = (12, 4))
    fig.subplots_adjust(bottom= 0.15)

    ax_pseudo = plt.subplot2grid((1, 3), (0, 0), colspan=1)
    ax_n_cogs = plt.subplot2grid((1, 3), (0, 1), colspan=1)
    ax_n_bases = plt.subplot2grid((1, 3), (0, 2), colspan=1)

    ax_pseudo.scatter(mrca_distance, n_pseudogenes)
    ax_n_cogs.scatter(mrca_distance, n_spore_forming_cogs)
    ax_n_bases.scatter(mrca_distance, n_spore_forming_cogs_nucleotides)

    ax_pseudo.set_xlabel("Root-to-spore-forming-ancestor distance", fontsize = 10)
    ax_pseudo.set_ylabel("Number of pseudogenes", fontsize = 10)

    ax_n_cogs.set_xlabel("Root-to-spore-forming-ancestor distance", fontsize = 10)
    ax_n_cogs.set_ylabel("Number of spore COGs", fontsize = 10)

    ax_n_bases.set_xlabel("Root-to-spore-forming-ancestor distance", fontsize = 10)
    ax_n_bases.set_ylabel("Number of spore COG nucleotides", fontsize = 10)


    ax_n_bases.set_yscale('log', basey=10)

    fig.subplots_adjust(hspace=0.25, wspace=0.25)
    fig_name = "%sspore_genome_statistic.png" % (config.analysis_directory)
    fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()





#build_spore_genome_statistics_dict()

plot_spore_genome_statistics()

