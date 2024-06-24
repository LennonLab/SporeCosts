import config
import os
from collections import Counter
#import matplotlib.pyplot as plt


data_directory = os.path.expanduser("~/GitHub/SporeCosts/data/")
scripts_directory = os.path.expanduser("~/GitHub/SporeCosts/scripts/")
figures_directory = os.path.expanduser("~/GitHub/SporeCosts/figures/")



n_fna_characters = 80


class classFASTA:

    # class to load FASTA file

    def __init__(self, fileFASTA):
        self.fileFASTA = fileFASTA

    def readFASTA(self):
        '''Checks for fasta by file extension'''
        file_lower = self.fileFASTA.lower()
        '''Check for three most common fasta file extensions'''
        if file_lower.endswith('.txt') or file_lower.endswith('.fa') or \
        file_lower.endswith('.fasta') or file_lower.endswith('.fna') or \
        file_lower.endswith('.fasta') or file_lower.endswith('.frn') or \
        file_lower.endswith('.faa') or file_lower.endswith('.ffn'):
            with open(self.fileFASTA, "r") as f:
                return self.ParseFASTA(f)
        else:
            print("Not in FASTA format.")

    def ParseFASTA(self, fileFASTA):
        '''Gets the sequence name and sequence from a FASTA formatted file'''
        fasta_list=[]
        for line in fileFASTA:
            if line[0] == '>':
                try:
                    fasta_list.append(current_dna)
            	#pass if an error comes up
                except UnboundLocalError:
                    #print "Inproper file format."
                    pass
                current_dna = [line.lstrip('>').rstrip('\n'),'']
            else:
                current_dna[1] += "".join(line.split())
        fasta_list.append(current_dna)
        '''Returns fasa as nested list, containing line identifier \
            and sequence'''
        return fasta_list



def clean_alignment(muscle_path, muscle_clean_path, min_n_sites=100, max_fraction_empty=0.8):

    # removes all sites where the fraction of empty bases across ASVs is greater than max_fraction_empty (putatively uninformative)
    # removes a sequencies with fewer than min_n_sites informative sites

    frn_aligned = classFASTA(muscle_path).readFASTA()

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

    frn_aligned = classFASTA(muscle_path).readFASTA()

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

    fasta_object = classFASTA(fasta_path).readFASTA()
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

        fasta_object = classFASTA(fasta_file_gene_muscle_clean_path).readFASTA()

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








split_fasta_file_and_run_muscle()

