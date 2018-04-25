## STEP 1: Separate fastq data into two files, those reads lacking vs. flexibly matching the flag epitope sequence (DYKDDDDK) downstream of an initiation codon ##

# (uncompress data stream) -> (select 75 nucleotide library) -> (remove sequences with 3 consecutive N) -> (remove flag epitope) -> (remove stop codon containing sequences) -> (write to text file)
gzcat RTseq.fastq.gz | grep -o -E "GACGACAAGTAG(.)+" | cut -c 13-87 | grep -x '.\{75\}' | grep -v NNN | grep -v -E "[AC]TG(...)*GA[TC]TA[TC]AA[AG]|[AC]TG(...)*TA[TC]AA[AG]GA[TC]|[AC]TG(...)*GA[TC]...AA[AG]GA[TC]" | grep -v -E "^(...)*T[AG]A|^(...)*TAG" > reads_noflag.txt

# (uncompress data stream) -> (select 75 nucleotide library) -> (remove sequences with 3 consecutive N) -> (select flag epitope containing sequences) -> (write to text file)
gzcat RTseq.fastq | grep -o -E "GACGACAAGTAG(.)+" | cut -c 13-87 | grep -x '.\{75\}' | grep -v NNN | grep -E "[AC]TG(...)*GA[TC]TA[TC]AA[AG]|[AC]TG(...)*TA[TC]AA[AG]GA[TC]" > reads_flag.txt
