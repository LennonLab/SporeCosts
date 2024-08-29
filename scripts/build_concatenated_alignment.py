import config
import pickle
import itertools
import collections
import os
import sys
from collections import Counter
#import matplotlib.pyplot as plt

import numpy
import matplotlib.pyplot as plt
from Bio import SeqIO

import utils



data_directory = os.path.expanduser("~/GitHub/SporeCosts/data/")
scripts_directory = os.path.expanduser("~/GitHub/SporeCosts/scripts/")
figures_directory = os.path.expanduser("~/GitHub/SporeCosts/figures/")


prokka_annotations_path = '%sgenome_analysis/prokka_annotations/' % data_directory
cog_dict_path = data_directory + 'cog_dict.pickle'
cog_fasta_dict_path = data_directory + 'cog_fasta_dict%s.pickle'

outgroup_label_dict = {True: '_outgroup', False: ''}


n_fna_characters = 80




def clean_alignment_old(muscle_path, muscle_clean_path, min_n_sites=100, max_fraction_empty=0.8):

    # removes all sites where the fraction of empty bases across ASVs is greater than max_fraction_empty (putatively uninformative)
    # removes a sequencies with fewer than min_n_sites informative sites

    frn_aligned = utils.classFASTA(muscle_path).readFASTA()

    n = len(frn_aligned)

    frn_aligned_seqs = [x[1] for x in frn_aligned]
    frn_aligned_seqs_names = [x[0] for x in frn_aligned]

    frns = []
    for site in zip(*frn_aligned_seqs):

        fraction_empty = site.count('-')/n

        if fraction_empty > max_fraction_empty:
            continue

        # skip site if it is uninformative
        if len(set([s for s in site if s != '-'])) == 1:
            continue

        frns.append(site)

    if len(frns) < min_n_sites:
        exit()

    clean_sites_list = zip(*frns)

    # skip if there are too few sites
    #if len(list(clean_sites_list)[0]) < min_n_sites:
    #    continue

    frn_aligned_clean = open(muscle_clean_path, 'w')

    for clean_sites_idx, clean_sites in enumerate(clean_sites_list):
        clean_sites_species = frn_aligned_seqs_names[clean_sites_idx]
        clean_sites_seq = "".join(clean_sites)

        frn_aligned_clean.write('>%s\n' % clean_sites_species)

        clean_sites_seq_split = [clean_sites_seq[i:i+n_fna_characters] for i in range(0, len(clean_sites_seq), n_fna_characters)]

        for seq in clean_sites_seq_split:

            frn_aligned_clean.write('%s\n' % seq)

        frn_aligned_clean.write('\n')


    frn_aligned_clean.close()




def clean_alignment(muscle_path, muscle_clean_path, max_fraction_empty=0.9):
    
    #  min_n_sites=100,

    # removes all sites where the fraction of empty bases across ASVs is greater than max_fraction_empty (putatively uninformative)
    # removes a sequencies with fewer than min_n_sites informative sites

    frn_aligned = utils.classFASTA(muscle_path).readFASTA()

    n = len(frn_aligned)

    frn_aligned_seqs = [x[1] for x in frn_aligned]
    frn_aligned_seqs_names = [x[0] for x in frn_aligned]

    frns = []
    for site in zip(*frn_aligned_seqs):

        fraction_empty = site.count('-')/n

        if fraction_empty > max_fraction_empty:
            continue

        # skip site if it is uninformative
        if len(set([s for s in site if s != '-'])) == 1:
            continue

        frns.append(site)

    #if len(frns) < min_n_sites:
    #    exit()

    clean_sites_list = zip(*frns)

    # skip if there are too few sites
    #if len(list(clean_sites_list)[0]) < min_n_sites:
    #    continue

    frn_aligned_clean = open(muscle_clean_path, 'w')

    for clean_sites_idx, clean_sites in enumerate(clean_sites_list):
        clean_sites_species = frn_aligned_seqs_names[clean_sites_idx]
        clean_sites_seq = "".join(clean_sites)

        frn_aligned_clean.write('>%s\n' % clean_sites_species)

        clean_sites_seq_split = [clean_sites_seq[i:i+n_fna_characters] for i in range(0, len(clean_sites_seq), n_fna_characters)]

        for seq in clean_sites_seq_split:

            frn_aligned_clean.write('%s\n' % seq)

        frn_aligned_clean.write('\n')


    frn_aligned_clean.close()




def split_fasta_file_and_run_muscle(min_n_genomes=160):

    fasta_path = '%sribosomal_genes_nucleotide.fasta' % data_directory

    fasta_object = utils.classFASTA(fasta_path).readFASTA()
    fasta_dict = {x[0]: x[1] for x in fasta_object}

    genes_all = [x[0].split('_')[-1] for x in fasta_object]
    #genes_set = list(set(genes_all))
    genomes_set = list(set(['_'.join(x[0].split('_')[:2]) for x in fasta_object]))

    # count genes
    genes_count_dict = dict(Counter(genes_all))
    genes_set = [x[0] for x in genes_count_dict.items() if x[1] >= min_n_genomes] 


    # plot counter
    #fig, ax = plt.subplots(figsize=(4,4))
    #n_gene_counts = genes_count_dict.values()
    #ax.hist(n_gene_counts)
    #ax.set_xlabel('Number of ribosomal protein coding genes')
    #ax.set_ylabel('Number of genomes')
    #fig.savefig('%sribosomal_gene_counts.png' % figures_directory, format='png', bbox_inches = "tight", pad_inches = 0.5, dpi = 600)

    # get set of genomes that have all the genes.
    genomes_to_keep = set(genomes_set.copy())
    for gene in genes_set:

        genomes_with_gene = []
        for g in list(genomes_to_keep):

            if '%s_%s' % (g, gene) in fasta_dict:

                genomes_with_gene.append(g)


        genomes_to_keep = genomes_to_keep & set(genomes_with_gene)


    genomes_to_keep = list(genomes_to_keep)


    for gene in genes_set:

        fasta_file_gene_path = '%sribosomal_genes_fasta/%s.fasta' % (data_directory, gene)
        fasta_file_gene = open(fasta_file_gene_path, 'w')

        for genome in genomes_to_keep:

            fasta_file_gene.write('>%s\n' % genome)
            sequence = fasta_dict['%s_%s' % (genome, gene)]

            for i in range(0, len(sequence), n_fna_characters):
                sequence_i = sequence[i : i + n_fna_characters]
                fasta_file_gene.write('%s\n' % sequence_i)
        
            fasta_file_gene.write('\n')

        fasta_file_gene.close()


        # run muscle
        fasta_file_gene_muscle_path = '%sribosomal_genes_fasta/%s_muscle.fasta' % (data_directory, gene)

        #os.system('muscle -in %s -out %s' % (fasta_file_gene_path, fasta_file_gene_muscle_path))

        #fasta_file_gene_muscle_clean_path = '%sribosomal_genes_fasta/%s_muscle_clean.fasta' % (data_directory, gene)

        #clean_alignment(fasta_file_gene_muscle_path, fasta_file_gene_muscle_clean_path)



    # concatenate alignment
    concat_dict = {}
    for gene in genes_set:

        fasta_file_gene_muscle_clean_path = '%sribosomal_genes_fasta/%s_muscle_clean.fasta' % (data_directory, gene)

        fasta_object = utils.classFASTA(fasta_file_gene_muscle_clean_path).readFASTA()

        for x in fasta_object:

            if x[0] not in concat_dict:
                concat_dict[x[0]] = x[1]

            else:
                concat_dict[x[0]] += x[1]


    fasta_file_concat_path = '%sribosomal_genes_aligned_concat.fasta' % (data_directory)
    fasta_file_concat = open(fasta_file_concat_path, 'w')

    for genome in concat_dict.keys():

        fasta_file_concat.write('>%s\n' % genome)
        sequence = concat_dict[genome]

        for i in range(0, len(sequence), n_fna_characters):
            sequence_i = sequence[i : i + n_fna_characters]
            fasta_file_concat.write('%s\n' % sequence_i)
    
        fasta_file_concat.write('\n')

    fasta_file_concat.close()





def build_cog_dict(require_outgroup=True):

    cog_dict = {}

    for genome_id in os.listdir(prokka_annotations_path):
        
        if genome_id == '.DS_Store':
            continue

        annotation_folder = os.path.join(prokka_annotations_path, genome_id)
        annotation_file = os.path.join(annotation_folder, f"{genome_id}.gbk")

        cog_dict[genome_id] = []
        for record in SeqIO.parse(annotation_file, "genbank"):
            for feature in record.features:
                if feature.type == "CDS":
                    cog_name = feature.qualifiers.get('db_xref', [None])[0]
                    
                    if cog_name == None:
                        continue

                    cog_name = cog_name.split(':')[1]

                    if cog_name not in cog_dict[genome_id]:
                        cog_dict[genome_id].append(cog_name)


    sys.stderr.write("Saving parameter dictionary...\n")
    with open(cog_dict_path, 'wb') as outfile:
        pickle.dump(cog_dict, outfile, protocol=pickle.HIGHEST_PROTOCOL)
    sys.stderr.write("Done!\n")





def plot_cog_count_hist(require_outgroup=True):
    
    cog_dict = pickle.load(open(cog_dict_path, "rb"))

    cog_all_flat = list(itertools.chain(*list(cog_dict.values())))

    cog_counts_dict = dict(collections.Counter(cog_all_flat))
    cog_counts_list = numpy.asarray(list(cog_counts_dict.values()))

    n_genomes_range = list(range(1, len(cog_dict)+1))

    n_genes_with_count = [sum(cog_counts_list>=i) for i in n_genomes_range]

    fig, ax = plt.subplots(figsize=(4,4))

    ax.plot(n_genomes_range, n_genes_with_count, lw=2, ls='-', c='k', color='dodgerblue')

    ax.set_xlabel('Number of genomes, ' + r'$n$', fontsize=12)
    ax.set_ylabel('Number of COGs present in ' + r'$\geq n$' +  ' genomes', fontsize=12)

    fig.subplots_adjust(hspace=0.15,wspace=0.2)
    fig_name = "%scog_count_survival.png" % config.analysis_directory
    fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()



def identify_cogs_in_all_genomes(min_n_genomes=170, require_outgroup=True):

    genome_key_cog_value_dict = pickle.load(open(cog_dict_path, "rb"))

    cog_all_flat = list(itertools.chain(*list(genome_key_cog_value_dict.values())))

    # only consider COGs present in outgroup
    if require_outgroup == True:
        cog_all_flat = [g for g in cog_all_flat if g in genome_key_cog_value_dict['genome_562']]

    cog_counts_dict = dict(collections.Counter(cog_all_flat))
    #cog_counts_list = numpy.asarray(list(cog_counts_dict.values()))

    cogs_to_keep = [key for key in cog_counts_dict.keys() if cog_counts_dict[key] >= min_n_genomes]

    # find genomes that have these COGs
    cog_key_genome_value = {}
    for c in cogs_to_keep:
        cog_key_genome_value[c] = []
        for g in genome_key_cog_value_dict.keys():

            if c in genome_key_cog_value_dict[g]:
                cog_key_genome_value[c].append(g)
    
    # intersection of all 
    genomes_all_flat = list(itertools.chain(*list(cog_key_genome_value.values())))
    genomes_counts_dict = dict(collections.Counter(genomes_all_flat))

    # find genomes with all the COGs
    genomes_to_keep = [key for key in genomes_counts_dict.keys() if genomes_counts_dict[key] == len(cog_key_genome_value)]

    return list(cog_key_genome_value.keys()), genomes_to_keep



def plot_cog_genomes_threshold(min_n_genomes=170, require_outgroup=True):

    genome_key_cog_value_dict = pickle.load(open(cog_dict_path, "rb"))

    min_n_genomes_range = list(range(100, len(genome_key_cog_value_dict)))

    n_cogs_all = []
    n_genomes_all = []
    for i in min_n_genomes_range:

        cog_list, genome_list = identify_cogs_in_all_genomes(min_n_genomes=i, require_outgroup=require_outgroup)
        n_cogs_all.append(len(cog_list))
        n_genomes_all.append(len(genome_list))


    fig, ax = plt.subplots(figsize=(4,4))
    ax2 = ax.twinx()


    ax.plot(min_n_genomes_range, n_cogs_all, lw=2, ls='-', c='b')
    ax2.plot(min_n_genomes_range, n_genomes_all, lw=2, ls='-', c='r')

    ax.set_xlabel('Min. number of genomes where COG is present')
    ax.set_ylabel('Number of COGs', color='b')
    ax2.set_ylabel('Number of genomes harboring all COGs', color='r')

    fig.subplots_adjust(hspace=0.15,wspace=0.2)
    fig_name = "%scog_genomes_threshold%s.png" % (config.analysis_directory, outgroup_label_dict[require_outgroup])
    fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()




def generate_cog_fasta_and_run_muscle(min_n_genomes=170, make_dict=False, require_outgroup=True):

    cog_list, genome_list = identify_cogs_in_all_genomes(min_n_genomes=min_n_genomes, require_outgroup=require_outgroup)
    sys.stderr.write("# COGs = %d, # genomes = %d\n" % (len(cog_list), len(genome_list)))

    if make_dict == True:

        seq_dict = {}
        for c in cog_list:
            seq_dict[c] = {}

        for genome_id in genome_list:

            annotation_folder = os.path.join(prokka_annotations_path, genome_id)
            annotation_file = os.path.join(annotation_folder, f"{genome_id}.gbk")

            for record in SeqIO.parse(annotation_file, "genbank"):
                for feature in record.features:
                    if feature.type == "CDS":
                        cog_name = feature.qualifiers.get('db_xref', [None])[0]

                        if cog_name != None:
                            cog_name = cog_name.split(':')[1]

                            if cog_name in cog_list:
                                seq_dict[cog_name][genome_id] = str(feature.location.extract(record).seq)


        sys.stderr.write("Saving parameter dictionary...\n")
        with open(cog_fasta_dict_path % outgroup_label_dict[require_outgroup], 'wb') as outfile:
            pickle.dump(seq_dict, outfile, protocol=pickle.HIGHEST_PROTOCOL)
        sys.stderr.write("Done!\n")


    seq_dict = pickle.load(open(cog_fasta_dict_path % outgroup_label_dict[require_outgroup], "rb"))

    # make individual fasta files
    for c in cog_list:

        print(c)

        fasta_file_gene_path = '%scogs_fasta%s/%s.fasta' % (data_directory, outgroup_label_dict[require_outgroup], c)
        fasta_file_gene = open(fasta_file_gene_path, 'w')

        for genome_id in genome_list:

            fasta_file_gene.write('>%s\n' % genome_id)
            sequence = seq_dict[c][genome_id]

            for i in range(0, len(sequence), n_fna_characters):
                sequence_i = sequence[i : i + n_fna_characters]
                fasta_file_gene.write('%s\n' % sequence_i)

            fasta_file_gene.write('\n')

        fasta_file_gene.close()

        # run muscle
        fasta_file_gene_muscle_path = '%scogs_fasta%s/%s_muscle.fasta' % (data_directory,  outgroup_label_dict[require_outgroup], c)
        os.system('muscle -in %s -out %s' % (fasta_file_gene_path, fasta_file_gene_muscle_path))
        fasta_file_gene_muscle_clean_path = '%scogs_fasta%s/%s_muscle_clean.fasta' % (data_directory,  outgroup_label_dict[require_outgroup], c)
        clean_alignment(fasta_file_gene_muscle_path, fasta_file_gene_muscle_clean_path)


    # concatenate alignment
    concat_dict = {}
    for c in cog_list:

        fasta_file_gene_muscle_clean_path = '%scogs_fasta%s/%s_muscle_clean.fasta' % (data_directory,  outgroup_label_dict[require_outgroup], c)

        fasta_object = utils.classFASTA(fasta_file_gene_muscle_clean_path).readFASTA()

        for x in fasta_object:

            if x[0] not in concat_dict:
                concat_dict[x[0]] = x[1]
            else:
                concat_dict[x[0]] += x[1]


    fasta_file_concat_path = '%scogs_concat_alignment%s.fasta' % (data_directory, outgroup_label_dict[require_outgroup])
    fasta_file_concat = open(fasta_file_concat_path, 'w')

    for genome in concat_dict.keys():

        fasta_file_concat.write('>%s\n' % genome)
        sequence = concat_dict[genome]

        for i in range(0, len(sequence), n_fna_characters):
            sequence_i = sequence[i : i + n_fna_characters]
            fasta_file_concat.write('%s\n' % sequence_i)
    
        fasta_file_concat.write('\n')

    fasta_file_concat.close()



#def run_muscle_cogs():






if __name__ == "__main__":

    #plot_cog_genomes_threshold()


    #split_fasta_file_and_run_muscle()

    generate_cog_fasta_and_run_muscle(make_dict=False)

    #generate_cog_fasta()

