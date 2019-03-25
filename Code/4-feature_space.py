#!/usr/bin/env python

import argparse
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser(description='feature space construction')
parser.add_argument('-f1', '--filename1', required=True, help='sequence data')
parser.add_argument('-f2', '--filename2', required=True, help='structural data')
parser.add_argument('-f3', '--filename3', required=True, help='bpp matrix')
parser.add_argument('-nr', '--numrows', required=True, help='number of rows to featurize')
parser.add_argument('-o', '--outfile', required=True, help='csv file name for the output feature space')
args = parser.parse_args()

num_rows = int(args.numrows)

data_seq = pd.read_csv(args.filename1, nrows=num_rows, sep='\t', header=None, names=['copy_number','sequence'])
data_struct = pd.read_csv(args.filename2, nrows=num_rows, sep='\t', header=None, names=['structure','mfe'])
data_bpp = pd.read_csv(args.filename3, header=None, sep='\t', names=['matrix'])
total = pd.concat([noflag_seq, noflag_struct, noflag_bpp], axis=1, join='inner')
print(total.shape)

seq_list = list(total['sequence'])
struct_list = list(total['structure'])
matrix_list = list(total['matrix'])

print('Constructing PSNT feature space')
n = 6
featspace_psnt = gen_psnt_fs(seq_list, ps_nucleotides_present(seq_list, first_N=n), n)

print('Constructing Stem feature space')
featspace_stem = gen_stem_fs(struct_list, seq_list, matrix_list)

featspace_combined = pd.merge(featspace_stem, featspace_psnt, left_index=True, right_index=True)
print(featspace_combined_utr.shape)

featspace_combined.to_csv(args.outfile, index=False, header=False)
