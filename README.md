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

Step 1: Navigate to the Data directory, download the RTseq.fastq.gz file, and execute the 1-preprocess.sh script.

```bash
$ cd Data
$ wget -O RTseq.fastq.gz https://www.dropbox.com/s/94z0o9c9hklhouz/RTseq.fastq.gz?dl=0
$ ../Code/1-preprocess.sh
```

Step 2: If proceeding from Step 1, there should now be reads_noflag.txt and reads_flag.txt files in the Data directory. If beginning at Step 2, download the required files and uncompress before executing 2-unique_abundance.py.

```bash
$ wget -O reads_noflag.txt.gz https://www.dropbox.com/s/kma4qcg508cjod9/reads_noflag.txt.gz?dl=0
$ wget -O reads_flag.txt.gz https://www.dropbox.com/s/zu4gr2zv495a1b8/reads_flag.txt.gz?dl=0
$ gunzip reads_noflag.txt.gz
$ gunzip reads_flag.txt.gz
$ ../Code/2-unique_abundance.py reads_noflag.txt
$ ../Code/2-unique_abundance.py reads_flag.txt
```

Step 3: Any unique sequence count files generated from Step 2 can now be explored for subsequence correlations, eg pairwise codon contingency tables.

```bash
$ ../Code/3-subsequence_stats.py -f unique_seq_reads_noflag.txt -n 38000 -s1s 0 -s1e 3 -s2s 3 -s2e 6
```

Step 4: This step takes in unique sequence counts and produces a feature space array. Each input sequence is represented as a row vector in a column space of position specific nucleotide identity as well as structural features of stem-loop formation. NOTE this step assumes the user has already prepared certain computationally expensive outputs from the ViennaRNA program RNAfold. Bypassing this step, feature spaces for positive and negative training examples from the in vitro selection can be downloaded below, as well as featurized 3' UTRs from the human transcriptome:

[positive examples / no FLAG](https://www.dropbox.com/s/8abnh7upkru5yd3/featspace_top10k_training_noflag.csv?dl=0)

[negative examples / FLAG](https://www.dropbox.com/s/ipqarqnj0mbbdsu/featspace_top10k_training_flag.csv?dl=0)

[unlabeled human 3' UTRs](https://www.dropbox.com/s/mhebiy6jszya7xd/featspace_utr_human.csv?dl=0)


Step 5: Assuming the existence of a feature space representation of the sequence data (position specific nucleotide identity, structural features), this step can be run in either training or prediction mode. If the featurized data comes in two files corresponding to positive / negative labels, then a classifier can be trained. If a trained classifier already exists from a prior run and the featurized data is unlabeled, then the classifier can be applied to predict labels.

```bash
$ ../Code/5-classifier.py -m train -ft1 featspace_train_pos.txt -ft2 featspace_train_neg -n1 1000 -n2 5000
$ ../Code/5-classifier.py -m predict -fp featspace_pred.txt -clf ../Out/classifier.pkl
```

### Sample repository structure after running steps 1-5: ###
---

```bash
$ tree readthru
    Code/
        1-preprocess.sh
        2-unique_abundance.py
        3-subsequence_stats.py
        4-feature_space.py
        5-classifier.py
    Data/
        featspace_pred.csv
        featspace_train_neg.csv
        featspace_train_pos.csv
        reads_flag.txt
        reads_noflag.txt
        RTseq.fastq.gz
        unique_seq_pred.txt
        unique_seq_train_pos.txt
        unique_seq_train_neg.txt
        unique_struct_pred.txt
        unique_struct_train_pos.txt
        unique_struct_train_neg.txt
        unique_BPP_pred.csv
        unique_BPP_train_pos.csv
        unique_BPP_train_neg.csv
    Output/
        boosting_importances.pdf
        boosting_loss_func.pdf
        classifier.pkl
        subsequence_stats_0_3_3_6_38000.csv
        predicted_labels.txt
    README.md
```
