#!/usr/bin/env python

from timeit import default_timer as timer
from scipy.spatial.distance import hamming
import numpy as np
import pandas as pd
import re
import sys
import os

## Creates sets of unique sequences and abundances for length 75 (N-free) and length 72 (N-free) sequences ##

os.mkdir('tmp')
filename  = sys.argv[1].strip()
read_set = open(filename, 'r')
read_list = []
N_list = []
no_N_list = []
list_72 = []

for line in read_set:
    read = str(line)
    read = read.rstrip()
    read_list.append(read)
    
for sequence in read_list:
    if 'N' in sequence:
        N_list.append(sequence)
    else:
        no_N_list.append(sequence)
        
for sequence in read_list:
    seq_72 = sequence[:72]
    list_72.append(seq_72)

no_N_read_array = np.array(no_N_list)
array_72 = np.array(list_72)

no_N_uniq_seqs, no_N_abunds = np.unique(no_N_read_array, return_counts=True)
uniq_seqs_72, abunds_72 = np.unique(array_72, return_counts=True)

no_N_seqs_df = pd.DataFrame({'seq' : no_N_uniq_seqs, 'abundance' : no_N_abunds})
seqs_72_df = pd.DataFrame({'seq' : uniq_seqs_72, 'abundance' : abunds_72})

no_N_seqs_df = no_N_seqs_df.sort_values(by=['abundance'], ascending=False)
seqs_72_df = seqs_72_df.sort_values(by=['abundance'], ascending=False)

np.savetxt('tmp/{0}_no_N.txt'.format(filename[:-4]), no_N_seqs_df.values, fmt='%s', delimiter="\t")
np.savetxt('tmp/{0}_72.txt'.format(filename[:-4]), seqs_72_df.values, fmt='%s', delimiter="\t")

print("unique sequence datasets generated")

## processes list of n unique sequences of length 72 (N-free) to remove those with point mutations from parent seq ##

# set the number of sequences to be processed
n = 50000

df_all = pd.read_table('tmp/{0}_72.txt'.format(filename[:-4]), delim_whitespace=True, names=('abundance', 'sequence'))
df = df_all[:n]

# splits sequences into 3 portions for sequential sorting and comparison
split_df = df.assign(seq_1=df.sequence.str[0:24], seq_2=df.sequence.str[24:48], seq_3=df.sequence.str[48:72])

temp_df = pd.DataFrame({'seq_3' : [],'seq_2' : [],'seq_1' : [],'sequence' : [], 'abundance' : []})  
clean_df = pd.DataFrame({'seq_3' : [],'seq_2' : [],'seq_1' : [],'sequence' : [], 'abundance' : []}) 
new_df = pd.DataFrame({'seq_3' : [],'seq_2' : [],'seq_1' : [],'sequence' : [], 'abundance' : []}) 
sort_list = ('seq_1', 'seq_2', 'seq_3')

new_df = new_df.append(split_df)

counter = 0
hits = 0

print(len(new_df), "sequences to process")
start = timer()

for i in range(3):
    label = sort_list[i]
    sorted_df = new_df.sort_values(by=[label], ascending=True)
    counter = 0
    for row in sorted_df.itertuples():
        counter += 1
        if counter % 10000 == 0:
            current = timer()
            print("round", i+1, counter, "sequences searched", "in", current-start, "seconds")
        row_df = sorted_df[counter-1:counter]
        seq1 = row[i+2]
        if len(temp_df) == 0:
            temp_df = temp_df.append(row_df)
            founder_seq = seq1
        else:
            if len(founder_seq) * hamming(founder_seq, seq1) < 3:
                hits += 1
                temp_df = temp_df.append(row_df)
            else:
                temp_df = temp_df.sort_values(by=['abundance'], ascending=False)
                winner = temp_df[0:1]
                clean_df = clean_df.append(winner)
                temp_df = temp_df.iloc[0:0]
                temp_df = temp_df.append(row_df)
                founder_seq = seq1
            if counter == len(sorted_df):
                temp_df = temp_df.sort_values(by=['abundance'], ascending=False)
                winner = temp_df[0:1]
                clean_df = clean_df.append(winner)
                temp_df = temp_df.iloc[0:0]
    new_df = clean_df
    clean_df = clean_df.iloc[0:0]
    current = timer()
    print("round", i+1, "complete", len(new_df), "sequences remaining")
    
new_df = new_df.sort_values(by=['abundance'], ascending=False)
np.savetxt('tmp/cleaned_{0}_72_50k.txt'.format(filename[:-4]), new_df.values, fmt='%s', delimiter="\t")

print("redundant sequences removed from length 72 dataset")

## processes list of n unique sequences of length 75 (N-free) to remove those with point mutations from parent seq ##

start = timer()

# set the number of sequences to be processed
n = 50000

df_all_75 = pd.read_table('tmp/{0}_no_N.txt'.format(filename[:-4]), delim_whitespace=True, names=('abundance', 'sequence'))
df_75 = df_all_75[:n]

# splits sequences into 3 portions for sequential sorting and comparison
split_df_75 = df_75.assign(seq_1=df_75.sequence.str[0:25], seq_2=df_75.sequence.str[25:50], seq_3=df_75.sequence.str[50:75])
sorted_1_df_75 = split_df_75[:n].sort_values(by=['seq_1'], ascending=True)

temp_df_75 = pd.DataFrame({'seq_3' : [],'seq_2' : [],'seq_1' : [],'sequence' : [], 'abundance' : []})  
clean_df_75 = pd.DataFrame({'seq_3' : [],'seq_2' : [],'seq_1' : [],'sequence' : [], 'abundance' : []}) 
new_df_75 = pd.DataFrame({'seq_3' : [],'seq_2' : [],'seq_1' : [],'sequence' : [], 'abundance' : []}) 
sort_list = ('seq_1', 'seq_2', 'seq_3')

new_df_75 = new_df_75.append(sorted_1_df_75)

counter = 0
hits = 0
#print new_df_75.shape
current = timer()
print("ready to start", "time", current-start)

for i in range(3):
    label = sort_list[i]
    sorted_df_75 = new_df_75.sort_values(by=[label], ascending=True)
    counter = 0
    for row in sorted_df_75.itertuples():
        counter += 1
        if counter % 10000 == 0:
            current = timer()
            print("round", i+1, counter, "sequences searched", "time", current-start)
        row_df_75 = sorted_df_75[counter-1:counter]
        seq1 = row[i+2]
        if len(temp_df_75) == 0:
            temp_df_75 = temp_df_75.append(row_df_75)
            founder_seq = seq1
        else:
            if len(founder_seq) * hamming(founder_seq, seq1) < 3:
                hits += 1
                temp_df_75 = temp_df_75.append(row_df_75)
            else:
                temp_df_75 = temp_df_75.sort_values(by=['abundance'], ascending=False)
                winner = temp_df_75[0:1]
                clean_df_75 = clean_df_75.append(winner)
                temp_df_75 = temp_df_75.iloc[0:0]
                temp_df_75 = temp_df_75.append(row_df_75)
                founder_seq = seq1
            if counter == len(sorted_df_75):
                temp_df_75 = temp_df_75.sort_values(by=['abundance'], ascending=False)
                winner = temp_df_75[0:1]
                clean_df_75 = clean_df_75.append(winner)
                temp_df_75 = temp_df_75.iloc[0:0]
    new_df_75 = clean_df_75
    clean_df_75 = clean_df_75.iloc[0:0]
    current = timer()
    print("round", i+1, "complete", len(new_df_75), "sequences remaining")

new_df_75 = new_df_75.sort_values(by=['abundance'], ascending=False)
np.savetxt('tmp/cleaned_{0}_no_N_50k.txt'.format(filename[:-4]), new_df_75.values, fmt='%s', delimiter="\t")

print("redundant sequences removed from length 75 dataset")

## assigns abundances from length 72 sequences to corresponding length 75 (N-free) sequences ##

new_df_75 = pd.read_table('tmp/cleaned_{0}_no_N_50k.txt'.format(filename[:-4]), delim_whitespace=True, names=('abundance', 'seq1', 'seq2', 'seq3', 'sequence'))
new_df = pd.read_table('tmp/cleaned_{0}_72_50k.txt'.format(filename[:-4]), delim_whitespace=True, names=('abundance', 'seq1', 'seq2', 'seq3', 'sequence'))

print(len(new_df), "sequences and", len(new_df_75), "sequences")

final_df = pd.DataFrame({'abundance' : [], 'sequence' : []})

start = timer()
total = 0
for row_a in new_df.itertuples():
    total += 1
    found = 0
    if total % 1000 == 0:
        current = timer()
        print(total, "sequences searched in", current-start, "seconds")
    counter = 0
    ab_a = row_a[1]
    if ab_a > 44:
        seq_a = row_a[5]
        for row_b in new_df_75.itertuples():
            counter += 1
            seq_b = row_b[5]
            if seq_a in seq_b:
                new_row_df = pd.DataFrame({'abundance' : [ab_a], 'sequence' : [seq_b]})
                final_df = final_df.append(new_row_df)
                new_df_75 = new_df_75.drop(new_df_75.index[counter-1])
                found = 1
                break 
    if found == 0:
        counter = 0
        for row_b in new_df_75.itertuples():
            counter += 1
            seq_b = row_b[5]
            if len(seq_a) * hamming(seq_a, seq_b[:72]) < 3:
                new_row_df = pd.DataFrame({'abundance' : [ab_a], 'sequence' : [seq_b]})
                final_df = final_df.append(new_row_df)
                new_df_75 = new_df_75.drop(new_df_75.index[counter-1])
                found = 1
                break 

print(len(final_df), "sequences in final dataset", len(new_df_75), "sequences discarded")

final_df = final_df.sort_values(by=['abundance'], ascending=False)
# np.savetxt('unique_seq_{0}_75.txt'.format(filename[:-4]), final_df.values, fmt='%s', delimiter="\t")


## generates length 135 sequences for structural analysis ##

final_df_long = pd.DataFrame({'abundance' : [], 'sequence' : []})

for row in final_df.itertuples():
    ab = row[1]
    seq = row[2] + 'GGCAGCGGCCATCATCACCATCACCACGGCGGTTCTATGGGAATGTCTGGATCTGGCTAT'
    long_row = pd.DataFrame({'abundance' : [ab], 'sequence' : [seq]})
    final_df_long = final_df_long.append(long_row)

final_df_long = final_df_long.sort_values(by=['abundance'], ascending=False)
np.savetxt('unique_seq_{0}.txt'.format(filename[:-4]), final_df_long.values, fmt='%s', delimiter="\t")

os.remove('tmp/{0}_no_N.txt'.format(filename[:-4]))
os.remove('tmp/{0}_72.txt'.format(filename[:-4]))
os.remove('tmp/cleaned_{0}_72_50k.txt'.format(filename[:-4]))
os.remove('tmp/cleaned_{0}_no_N_50k.txt'.format(filename[:-4]))
os.rmdir('tmp')

print("final datasets generated")
print("finished")
