## Interrogation of eukaryotic stop codon readthrough signals by in vitro RNA selection ##

This computational workflow accompanies the [published manuscript](https://pubs.acs.org/doi/abs/10.1021/acs.biochem.8b01280) with the aforementioned title.
These scripts recapitulate the analysis presented for detecting stop codon readthrough contexts from deep sequencing of a randomized RNA library after 3 rounds of in vitro selection.

### Requirements: ###
---

- UNIX like operating system
- wget
- Vienna RNAfold
- python 2.7.x or 3.x
    - numpy
    - pandas
    - matplotlib
    - biopython
    - scipy
    - scikit-learn

### Setup: ###
---

To begin the workflow at Step 1, the following compressed fastq file will be required:

<https://www.dropbox.com/s/94z0o9c9hklhouz/RTseq.fastq.gz?dl=0>

To begin the workflow at Step 2, with the reads already divided into RT and decoy sets, the following two compressed text files will be required:

<https://www.dropbox.com/s/kma4qcg508cjod9/reads_noflag.txt.gz?dl=0>
<https://www.dropbox.com/s/zu4gr2zv495a1b8/reads_flag.txt.gz?dl=0>

To begin the workflow at Step 3, with the unique sequence abundances already computed, the "unique_seq_reads_*" text files contained within the Data directory are provided. 

### Usage: ###
---

Step 1. Navigate to the Data directory, download the RTseq.fastq.gz file, and execute the 1-preprocess.sh script.

```bash
$ cd Data
$ wget -O RTseq.fastq.gz https://www.dropbox.com/s/94z0o9c9hklhouz/RTseq.fastq.gz?dl=0
$ ../Code/1-preprocess.sh
```

Step 2. If proceeding from Step 1, there should now be reads_noflag.txt and reads_flag.txt files in the Data directory. If beginning at Step 2, download the required files and uncompress before executing 2-unique_abundance.py.

```bash
$ wget -O reads_noflag.txt.gz https://www.dropbox.com/s/kma4qcg508cjod9/reads_noflag.txt.gz?dl=0
$ wget -O reads_flag.txt.gz https://www.dropbox.com/s/zu4gr2zv495a1b8/reads_flag.txt.gz?dl=0
$ gunzip reads_noflag.txt.gz
$ gunzip reads_flag.txt.gz
$ ../Code/2-unique_abundance.py reads_noflag.txt
$ ../Code/2-unique_abundance.py reads_flag.txt
```

Step 3. If proceeding from Step 2, there should now be unique_seq_reads_noflag.txt and unique_seq_reads_flag.txt files in the Data directory. If beginning at Step 3, use the provided precomputed_unique_seq_reads_noflag.txt and precomputed_unique_seq_reads_flag.txt files as input to 3-correlation_stats.py.

```bash
$ ../Code/3-correlation_stats.py -f unique_seq_reads_noflag.txt -n 38000 -s1s 0 -s1e 3 -s2s 3 -s2e 6
$ ../Code/3-correlation_stats.py -f unique_seq_reads_flag.txt -n 3000 -s1s 0 -s1e 3 -s2s 3 -s2e 6
```

Step 4. 

Step 5. provide a motif sequence "m" to be used for nucleotide variant analysis.


### Final Repo Structure After Steps 1-6: ###
---

```bash
$ tree readthru
    Code/
        1-preprocess.sh
        2-unique_abundance.py
        3-correlation_stats.py
        4-feature_space.py
        5-classifier_training.py
        6-classifier_humanUTR.py
    Data/
        precomputed_unique_seq_reads_flag.txt
        precomputed_unique_seq_reads_noflag.txt
        reads_flag.txt
        reads_noflag.txt
        RTseq.fastq.gz
        unique_seq_reads_flag.txt
        unique_seq_reads_noflag.txt
    Output/
        precomputed_unique_seq_reads_flag_0_3_3_6_3000.csv
        precomputed_unique_seq_reads_noflag_0_3_3_6_38000.csv
        unique_seq_reads_flag_0_3_3_6_3000.csv
        unique_seq_reads_noflag_0_3_3_6_38000.csv
    README.md
```
