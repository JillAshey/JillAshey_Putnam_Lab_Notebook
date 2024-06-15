On my local computer, make fasta file with the given miRNA names. 

```
cd /Users/jillashey/Desktop/PutnamLab/Repositories/deep-dive/DEF-cross-species/output/10-shortRNA-ShortStack-comparison

awk -F, 'NR>1 {print ">" $22 "\n" $11}' Apul_results_mature_named.csv > Apul_results_mature_named.fasta

awk -F, 'NR>1 {print ">" $22 "\n" $11}' Peve_results_mature_named.csv > Peve_results_mature_named.fasta

awk -F, 'NR>1 {print ">" $22 "\n" $11}' Ptuh_results_mature_named.csv > Ptuh_results_mature_named.fasta
```

Copy the fasta files to the server 

```
cd /data/putnamlab/jillashey/e5
mkdir mirna_seqs
cd mirna_seqs
```

To run miranda, I need to identify the 3' UTR ends. 

First, identify counts of each feature from gff file 

```
GFF_FILE="/data/putnamlab/jillashey/genome/Amil_v2.01/Amil.all.maker.noseq.gff"
genome="/data/putnamlab/jillashey/genome/Amil_v2.01/Amil.v2.01.chrs.fasta"

grep -v '^#' ${GFF_FILE} | cut -s -f 3 | sort | uniq -c | sort -rn > all_features.txt

1572090 match_part
 558048 match
 201232 exon
 176923 CDS
  88262 expressed_sequence_match
  84724 protein_match
  38247 gene
  28188 mRNA
  19507 three_prime_UTR
  15776 five_prime_UTR
  10059 tRNA
    854 contig
```

There are a decent amount of three prime UTRs annotated in the Amil gff. Calculate the length of the 3' UTRs. 

```
awk '$3 == "three_prime_UTR" {print $0, $5 - $4 + 1}' OFS="\t" Amil.all.maker.noseq.gff > three_prime_UTR_lengths.txt

awk '{sum += $NF} END {if (NR > 0) print sum / NR}' three_prime_UTR_lengths.txt
```

Average length of 3' UTRs is 424.277. 