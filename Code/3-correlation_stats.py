#!/usr/bin/env python

from scipy import stats as stats
from scipy.stats import chi2_contingency
from scipy.stats.contingency import expected_freq
from scipy.stats.contingency import margins
from statsmodels.sandbox.stats.multicomp import fdrcorrection0
from statsmodels.compat.python import range
from statsmodels.compat.collections import OrderedDict
import pandas as pd
import numpy as np
import itertools
import argparse

parser = argparse.ArgumentParser(description='calculate codon co-occurrence statistics')
parser.add_argument('-f', '--filename', required=True, help='unique sequence abundance text file, two column tab delimited')
parser.add_argument('-s1s', '--seq1start', required=True, help='subsequence 1 starting position')
parser.add_argument('-s1e', '--seq1end', required=True, help='subsequence 1 ending position')
parser.add_argument('-s2s', '--seq2start', required=True, help='subsequence 2 starting position')
parser.add_argument('-s2e', '--seq2end', required=True, help='subsequence 2 ending position')
parser.add_argument('-n', '--numseq', required=False, default=1000, help='number of most abundant sequences to analyze')
args = parser.parse_args()

q = open(args.filename, 'r')

seq_list = []
str1_list = []
str2_list = []
rows = []
cols = []

for line in q:
	sequence = line.split()[1].strip()
	seq_list.append(sequence)
seq_list = seq_list[:int(args.numseq)]

for seq in seq_list:
	str1 = seq[int(args.seq1start):int(args.seq1end)]
	str1_list.append(str1)
	str2 = seq[int(args.seq2start):int(args.seq2end)]
	str2_list.append(str2)

for substring in str1_list:
	if substring not in rows:
		rows.append(substring)
for substring in str2_list:
	if substring not in cols:
		cols.append(substring)

# rows is a list of the unique sequence substrings for string 1
# cols is a list of the unique sequence substrings for string 2

paired_str = zip(str1_list, str2_list)

# paired_str is a list of all tuples of string 1, string 2; has len equal to int(args.numseq)

obs = np.zeros((len(rows), len(cols)))
odds_table = np.zeros((len(rows), len(cols)))
pvalue_table = np.zeros((len(rows), len(cols)))
adj_pvalue_table = np.zeros((len(rows), len(cols)))
phi_table = np.zeros((len(rows), len(cols)))
chi2_table = np.zeros((len(rows), len(cols)))
residual_table = np.zeros((len(rows), len(cols)))

for pair in paired_str:
	str1 = pair[0]
	str2 = pair[1]
	r1 = rows.index(str1)
	c1 = cols.index(str2)
	obs[r1,c1] += 1
exp = expected_freq(obs)

# obs is a numpy array containing the number of observations for each unique pairing of string 1 and string 2
# exp is a numpy array containing the expected number of observation 

m0_obs, m1_obs = margins(obs)
m0_exp, m1_exp = margins(exp)

for row_idx in range(len(rows)):
	for col_idx in range(len(cols)):
		contingency_table = np.zeros((2,2))
		row_col_obs = obs[row_idx,col_idx]
		row_else_obs = m0_obs[row_idx,0] - row_col_obs
		else_col_obs = m1_obs[0,col_idx] - row_col_obs
		else_else_obs = int(args.numseq) - m0_obs[row_idx,0] - else_col_obs
		contingency_table[0,0] = row_col_obs
		contingency_table[0,1] = row_else_obs
		contingency_table[1,0] = else_col_obs
		contingency_table[1,1] = else_else_obs
		row_col_exp = exp[row_idx,col_idx]
		chi2, p, dof, ex = chi2_contingency(contingency_table)
		phi = (chi2/int(args.numseq))**0.5
		oddsratio, pvalue = stats.fisher_exact(contingency_table)
		odds_table[row_idx, col_idx] = oddsratio
		pvalue_table[row_idx, col_idx] = pvalue
		adj_pvalue = 1-((1-pvalue)**(len(rows)*len(cols)))
		adj_pvalue_table[row_idx, col_idx] = adj_pvalue
		row_col_diff = row_col_obs - row_col_exp
		residual = row_col_diff/(row_col_exp**0.5)
		if oddsratio < 1:
			phi = -phi
			chi2 = -chi2
		phi_table[row_idx, col_idx] = phi
		chi2_table[row_idx, col_idx] = chi2
		residual_table[row_idx, col_idx] = residual

df = pd.DataFrame({'pvalue' : [], 'adj_pvalue' : [], 'chi2' : [],'phi' : [],'residual' : [],'string_1' : [],'string_2' : [],'sequence' : [],'odds_ratio' : [],'expected' : [],'observed' : []})

for row in range(len(rows)):
	for col in range(len(cols)):
		pvalue = pvalue_table[row,col]
		adj_pvalue = adj_pvalue_table[row,col]
		odds = odds_table[row,col]
		phi = phi_table[row,col]
		chi2 = chi2_table[row,col]
		residual = residual_table[row,col]
		str1 = rows[row]
		str2 = cols[col]
		sequence = str1+str2
		expected = exp[row,col]
		observed = obs[row,col]
		df_row = pd.DataFrame({'pvalue' : [pvalue], 'adj_pvalue' : [adj_pvalue], 'chi2' : [chi2],'phi' : [phi],'residual' : [residual],'string_1' : [str1],'string_2' : [str2],'sequence' : [sequence],'odds_ratio' : [odds],'expected' : [expected],'observed' : [observed]})
		df = df.append(df_row)

df.to_csv('../Output/{0}_{1}_{2}_{3}_{4}_{5}.csv'.format(args.filename[:-4],
	int(args.seq1start), int(args.seq1end), int(args.seq2start), int(args.seq2end), int(args.numseq)), index=False)
