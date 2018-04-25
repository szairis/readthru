from itertools import combinations
from glob import glob
import numpy as np
import pandas as pd
from scipy.spatial.distance import hamming
from sklearn import metrics, tree, ensemble, svm, model_selection
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.stats import fisher_exact
from itertools import combinations, product
from Bio import pairwise2, Seq, SeqUtils
import networkx as nx

def sense_codons(seq_nt):
    ### returns the number of sense amino acids after stop codon and the peptide itself
    seq_aa = str(Seq.Seq(seq_nt).translate())
    stop_codon_pos = seq_aa.find('*')
    if stop_codon_pos == -1:
        return len(seq_aa), seq_aa
    else:
        return stop_codon_pos, seq_aa[0 : stop_codon_pos + 1]

def get_sizes_and_intersection(motif_1, motif_2, seq_list):
    ### notation will be '1_CAA' for specifying codon 1
    start_1 = 3 * (int(motif_1[0]) - 1)
    codon_1 = motif_1[2:].replace('U','T')
    stop_1 = start_1 + 3
    start_2 = 3 * (int(motif_2[0]) - 1)
    codon_2 = motif_2[2:].replace('U','T')
    stop_2 = start_2 + 3
    
    size_1 = len([seq for seq in seq_list if seq[start_1 : stop_1] == codon_1])
    size_2 = len([seq for seq in seq_list if seq[start_2 : stop_2] == codon_2])
    size_intersection = len([seq for seq in seq_list if (seq[start_1 : stop_1] == codon_1 and \
                                                         seq[start_2 : stop_2] == codon_2)])
    return size_1, size_2, size_intersection

def gen_psdb_fs(struct_list, start_pos=1, end_pos=20):
    ps_dotbrackets = ['{}_('.format(x) for x in map(str, range(start_pos, end_pos + 1))] + \
                     ['{}_.'.format(x) for x in map(str, range(start_pos, end_pos + 1))] + \
                     ['{}_)'.format(x) for x in map(str, range(start_pos, end_pos + 1))]
    feature_space = pd.DataFrame(np.zeros((len(struct_list),len(ps_dotbrackets))), columns=ps_dotbrackets)
    for n in range(len(struct_list)):
        if (n+1) % 1000 == 0:
            print('{0} / {1}'.format(n + 1, len(struct_list)))
        for pos in range(start_pos, end_pos + 1):
            db = '{0}_{1}'.format(pos, struct_list[n][pos - 1])
            feature_space[db][n] = 1
    return feature_space

def ps_nucleotides_present(seq_list, first_N=10):
    if first_N == 6:
        return ['1_A', '1_C', '1_G', '1_T',
                '2_A', '2_C', '2_G', '2_T',
                '3_A', '3_C', '3_G', '3_T',
                '4_A', '4_C', '4_G', '4_T',
                '5_A', '5_C', '5_G', '5_T',
                '6_A', '6_C', '6_G', '6_T'
               ]
    start = 0
    stop = first_N
    nucleotides = []
    for seq in seq_list:
        for pos in range(start, stop):
            nucleotides.append('{0}_{1}'.format(pos+1, seq[pos]))
    return list(set(nucleotides))

def gen_psnt_fs(seq_list, ps_nucleotides, first_N=10):
    start = 0
    stop = first_N
    feature_space = pd.DataFrame(np.zeros((len(seq_list),len(ps_nucleotides))), columns=ps_nucleotides)
    for n in range(len(seq_list)):
        if (n+1) % 1000 == 0:
            print('{0} / {1}'.format(n + 1, len(seq_list)))
        for pos in range(start, stop):
            nt = '{0}_{1}'.format(pos + 1, seq_list[n][pos])
            feature_space[nt][n] = 1
    return feature_space

def kmers_present(seq_list, k, first_N=30):
    start = 0
    stop = first_N - k
    kmers = []
    for seq in seq_list:
        N = len(seq)
        for pos in range(start, stop):
            kmers.append('{0}'.format(seq[pos:pos+k]))
    return list(set(kmers))

def gen_mk_fs(seq_list, kmers, m, first_N=30):
    k = len(kmers[0])
    start = 0
    stop = first_N - k
    feature_space = pd.DataFrame(np.zeros((len(seq_list),len(kmers))), columns=kmers)
    for n in range(len(seq_list)):
        if (n+1) % 100 == 0:
            print('{0} / {1}'.format(n + 1, len(seq_list)))
        N = len(seq_list[n])
        for pos in range(start, stop):
            kmer = seq_list[n][pos:pos+k]
            for feature in kmers:
                if hamming(kmer, feature) <= m:
                    feature_space[feature][n] = 1
    return feature_space

def base_pairings(dotbracket):
    pairing_list = []
    dotbracket = list(dotbracket)
    for j in range(len(dotbracket)):
        if dotbracket[j] == ')':
            for i in range(j - 1, -1, -1):
                if dotbracket[i] == '(':
                    pairing = (i, j)
                    pairing_list.append(pairing)
                    dotbracket[i] = '.'
                    break
    return pairing_list

def characterize_first_stem(struct, seq, matrix_string, skip_N=0, first_N=80, max_length=40, min_bp=8, min_bp_percent=0.9):
    pairing = base_pairings(struct)
    matrix = np.reshape(np.array(map(float, matrix_string.split(','))), (135,136))
    
    ###############################################
    
    found_a_real_stem = False
    best_position = 0
    start_search_at = max(struct.find('('), skip_N)
    for position in range(start_search_at, first_N):
        best_length = 1
        stem_length_list = []
        if struct[position] != '(':
            continue
        for length in range(1, max_length + 1):
            percent_left_paren = struct[position : position + length].count('(') / float(length)
            if percent_left_paren >= min_bp_percent and length > min_bp and struct[position + length - 1] != '.' \
            and struct[position + length - 1] != ')':
                stem_length_list.append(length)
        if len(stem_length_list): best_length = max(stem_length_list)
        if best_length > min_bp:
            best_position = position
            found_a_real_stem = True
            break
    
    if not found_a_real_stem: return [-1, 0, 0.0, 0.0, 0.0, 0]
    GC = round(SeqUtils.GC(Seq.Seq(seq[best_position : best_position + best_length])), 1)
    
    ###############################################
    
    startpos_partner = [pair[1] for pair in pairing if pair[0] == best_position][0]
    endpos_partner = [pair[1] for pair in pairing if pair[0] == best_position + best_length - 1][0]
    num_bulge = struct[best_position : best_position + best_length].count('.')
    num_bulge += struct[endpos_partner : startpos_partner].count('.')
    
    ################################################
    
    bpp_list = []
    for i in range(best_position, best_position + best_length):
        if struct[i] == '.':
            continue
        j = [pair[1] for pair in pairing if pair[0] == i][0]
        bpp_list.append(matrix[i,j])
    bpp_array = np.array(bpp_list)
    bpp_mean = round(bpp_array.mean(), 3)
    bpp_std = round(bpp_array.std(), 3)
    
    ################################################
    
    return [best_position + 1, best_length, GC, bpp_mean, bpp_std, num_bulge]

def gen_stem_fs(struct_list, seq_list, matrix_list):
    feature_space = pd.DataFrame(np.zeros((len(struct_list), 6)),
                                 columns=['position','length','GC','bpp_mean','bpp_std','bulges'])
    for n in range(len(struct_list)):
        if (n + 1) % 1000 == 0:
            print('{0} / {1}'.format(n + 1, len(struct_list)))
        feature_space.ix[n, :] = characterize_first_stem(struct_list[n], seq_list[n], matrix_list[n])
    return feature_space

def get_prob_matrix(base_path, file_number):
    prob_matrix = np.zeros((135,136))
    with open('{0}/seq_{1}.txt'.format(base_path, file_number)) as fh:
        for line in fh:
            i, j, value = line.strip().split()
            i = int(i) - 1
            j = int(j) - 1
            value = float(value) ** 2
            prob_matrix[i, j] = value
            prob_matrix[j, i] = value
    largest_sum = 0
    for row in range(135):
        if sum(prob_matrix[row, :]) > largest_sum:
            largest_sum = sum(prob_matrix[row, :])
    for row in range(135):
            prob_matrix[row, :] = prob_matrix[row, :] / largest_sum
            prob_matrix[row, 135] = 1.0 - sum(prob_matrix[row, :])
    return prob_matrix

def positional_entropy(prob_matrix):
    epsilon = 10**-10
    H_i = -np.sum(prob_matrix * np.log(prob_matrix + epsilon), axis=1)
    return H_i