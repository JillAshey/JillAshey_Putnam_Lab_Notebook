
Options for 3' UTR prediction 

- Gene Ext (see [github](https://github.com/sebepedroslab/GeneExt))
- Calling 1000 bases ahead of gene the UTR and making sure no gene overlap occurs 
- Looking at other cnidarian species for 3' UTR lengths 


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

### Apul

Since the Apul genome is assembled and the annotation is finishing up, the rest of e5 data analysis will proceed with the Apul genome instead of the Amil genome. Therefore, I need to identify the 3'UTR ends because that is what miranda requires to run. 

Identify the counts of each feature from the gff file 

```
cd /data/putnamlab/jillashey/e5/refs

GFF_FILE="/data/putnamlab/tconn/predict_results/Acropora_pulchra.gff3"
genome="/data/putnamlab/REFS/Apul/apul.hifiasm.s55_pa.p_ctg.fa.k32.w100.z1000.ntLink.5rounds.masked.fa"

grep -v '^#' ${GFF_FILE} | cut -s -f 3 | sort | uniq -c | sort -rn > all_features_apul.txt

 209537 exon
 201613 CDS
  44371 gene
  36447 mRNA
   7924 tRNA
```

No UTRs annotated. This is expected, as Trinity is still finishing up the Apul genome annotation. Extract gene as a feature types and generate individual gff for the gene feature. The mRNA feature will be used to as a spatial reference to create 3'UTRs. Gene will not be used because the genes also code for tRNAs. 

```
grep $'\tmRNA\t' ${GFF_FILE} > apul_GFFannotation.mRNA.gff
```

Extract scaffold lengths
```
cat is ${genome} | awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > apul.Chromosome_lenghts.txt

wc -l apul.Chromosome_lenghts.txt 
174 apul.Chromosome_lenghts.txt
```

Extract scaffold names 

```
awk -F" " '{print $1}' apul.Chromosome_lenghts.txt > apul.Chromosome_names.txt
```

The following script will sort the mRNA gff, extract 1kb down the 3' end of mRNA, subtract portions of the 1kb flank (representing the 3'UTR) from any overlapping mRNA, and make fasta file of the 3'UTRs. In the e5 scripts folder: `nano bed_apul.sh`

```
#!/bin/bash 
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=5
#SBATCH --export=NONE
#SBATCH --mem=250GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/e5/scripts
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

module load BEDTools/2.30.0-GCC-11.3.0

cd /data/putnamlab/jillashey/e5/refs

echo "Sorting gffs by chromosome" $(date)

sortBed -faidx apul.Chromosome_names.txt -i apul_GFFannotation.mRNA.gff > apul_GFFannotation.mRNA_sorted.gff

echo "Sorting complete!" $(date)

echo "Extracting 1kb 3' UTRs" $(date)

bedtools flank -i apul_GFFannotation.mRNA_sorted.gff -g apul.Chromosome_lenghts.txt -l 0 -r 1000 -s | awk '{gsub("gene","3prime_UTR",$3); print $0 }' | awk '{if($5-$4 > 3)print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' | tr ' ' '\t' > apul.GFFannotation.3UTR_1kb.gff

echo "Subtract portions of UTRs that overlap nearby genes" $(date)

bedtools subtract -a apul.GFFannotation.3UTR_1kb.gff -b apul_GFFannotation.mRNA_sorted.gff > apul.GFFannotation.3UTR_1kb_corrected.gff 
echo "3' UTRs identified!" $(date)

echo "Extracting 3' UTR sequences" $(date)

bedtools getfasta -fi /data/putnamlab/REFS/Apul/apul.hifiasm.s55_pa.p_ctg.fa.k32.w100.z1000.ntLink.5rounds.masked.fa -bed apul.GFFannotation.3UTR_1kb_corrected.gff -fo apul_3UTR_1kb.fasta -fullHeader

echo "Sequence extraction complete!" $(date)
```

Submitted batch job 341153. Ran very fast. 

```
zgrep -c ">" apul_3UTR_1kb.fasta 
37359
```

Time to run miranda for Apul. In the e5 scripts folder: `nano miranda_strict_all_1kb_apul.sh`

```
#!/bin/bash -i
#SBATCH -t 48:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --export=NONE
#SBATCH --mem=500GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/e5/scripts
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

echo "Apul starting miranda run with all genes and miRNAs with score cutoff >100, energy cutoff <-10, and strict binding invoked"$(date)

module load Miniconda3/4.9.2
conda activate /data/putnamlab/conda/miranda 

miranda /data/putnamlab/jillashey/e5/mirna_seqs/Apul_results_mature_named.fasta /data/putnamlab/jillashey/e5/refs/apul_3UTR_1kb.fasta -sc 100 -en -10 -strict -out /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_1kb_apul.tab

conda deactivate

echo "miranda run finished!"$(date)
echo "counting number of putative interactions predicted" $(date)

zgrep -c "Performing Scan" /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_1kb_apul.tab

echo "Parsing output" $(date)
grep -A 1 "Scores for this hit:" /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_1kb_apul.tab | sort | grep '>' > /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_1kb_parsed_apul.txt

echo "counting number of putative interactions predicted" $(date)
wc -l /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_1kb_parsed_apul.txt

echo "Apul miranda script complete" $(date)
```

Submitted batch job 341156. Started running immediately. 

#### using the Amil genome ?

Identify counts of each feature from gff 

```
cd /data/putnamlab/jillashey/genome/Amil_v2.01

grep -v '^#' GCF_013753865.1_Amil_v2.1_genomic.gff | cut -s -f 3 | sort | uniq -c | sort -rn > all_features_apul_amil.txt

390533 exon
 317969 CDS
  41860 mRNA
  36904 gene
  22387 cDNA_match
   6128 lnc_RNA
   5871 pseudogene
   2066 transcript
   1413 tRNA
    854 region
    283 rRNA
    170 snRNA
     62 snoRNA
      1 guide_RNA
```

Generate gff for genes

```
grep $'\tgene\t' GCF_013753865.1_Amil_v2.1_genomic.gff > Amil_gene.gtf
```

Extract scaffold lengths. 

```
sed 's/^\(>[^ ]*\).*/\1/' GCF_013753865.1_Amil_v2.1_genomic.fna > GCF_013753865.1_Amil_v2.1_genomic_modified.fna


cat is GCF_013753865.1_Amil_v2.1_genomic_modified.fna | awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > Amil_Chromosome_lengths.txt
```

Extract scaffold names

```
awk -F" " '{print $1}' Amil_Chromosome_lengths.txt > Amil_Chromosome_names.txt
```

Create 1kb 3'UTR in gff

```
interactive

module load BEDTools/2.30.0-GCC-11.3.0

bedtools flank -i Amil_gene.gtf -g Amil_Chromosome_lengths.txt -l 0 -r 1000 -s | awk '{gsub("gene","3prime_UTR",$3); print $0 }' | awk '{if($5-$4 > 3)print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' | tr ' ' '\t' > Amil_3UTR_1kb.gtf
```

Subtract portions of 3'UTR that overlap nearby genes

```
bedtools subtract -a Amil_3UTR_1kb.gtf -b Amil_gene.gtf > Amil_3UTR_1kb_corrected.gtf
```

Extract 3'UTR sequences from genome 

```
awk '{print $1 "\t" $4-1 "\t" $5 "\t" $9 "\t" "." "\t" $7}' Amil_3UTR_1kb_corrected.gtf | sed 's/"//g' > Amil_3UTR_1kb_corrected.bed

bedtools getfasta -fi GCF_013753865.1_Amil_v2.1_genomic_modified.fna -bed Amil_3UTR_1kb_corrected.bed -fo Amil_3UTR_1kb.fasta -name
```

Run miranda for Apul data using Amil genome. In the scripts folder: `nano miranda_strict_all_1kb_amil_apul.sh`


```
#!/bin/bash -i
#SBATCH -t 48:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --export=NONE
#SBATCH --mem=500GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/e5/scripts
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

echo "Apul using Amil genome starting miranda run with all genes and miRNAs with energy cutoff <-20 and strict binding invoked"$(date)

module load Miniconda3/4.9.2
conda activate /data/putnamlab/conda/miranda 

miranda /data/putnamlab/jillashey/e5/mirna_seqs/Apul_results_mature_named.fasta /data/putnamlab/jillashey/genome/Amil_v2.01/Amil_3UTR_1kb.fasta -en -20 -strict -out /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_1kb_amil_apul.tab

conda deactivate

echo "miranda run finished!"$(date)
echo "counting number of putative interactions predicted" $(date)

zgrep -c "Performing Scan" /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_1kb_amil_apul.tab

echo "Parsing output" $(date)
grep -A 1 "Scores for this hit:" /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_1kb_amil_apul.tab | sort | grep '>' > /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_1kb_parsed_amil_apul.txt

echo "counting number of putative interactions predicted" $(date)
wc -l /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_1kb_parsed_amil_apul.txt

echo "Apul using Amil genome miranda script complete" $(date)
```

Submitted batch job 346077. Results: 

```
counting number of putative interactions predicted Wed Oct 30 16:54:24 EDT 2024
1227134
Parsing output Wed Oct 30 16:54:29 EDT 2024
counting number of putative interactions predicted Wed Oct 30 16:54:30 EDT 2024
4144 /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_1kb_parsed_amil_apul.txt
```

### Peve

To run miranda, I need to identify the 3' UTR ends. 

First, identify counts of each feature from gff file 

```
GFF_FILE="/data/putnamlab/jillashey/genome/Peve/Porites_evermanni_v1.annot.gff"
genome="/data/putnamlab/jillashey/genome/Peve/Porites_evermanni_v1.fa"

grep -v '^#' ${GFF_FILE} | cut -s -f 3 | sort | uniq -c | sort -rn > all_features_peve.txt

 231320 CDS
  40389 mRNA
  15098 UTR
```

There are a decent amount of UTRs annotated in the Peve gff but they aren't differentiated between 5' and 3'. Calculate the length of the 3' UTRs. 

```
awk '$3 == "UTR" {print $0, $5 - $4 + 1}' OFS="\t" Porites_evermanni_v1.annot.gff > UTR_lengths_Peve.txt

awk '{sum += $NF} END {if (NR > 0) print sum / NR}' UTR_lengths_Peve.txt
```

Average length of UTRs is 296.414. Extract feature types and generate individual gffs for each feature. 

```
grep $'\tCDS\t' ${GFF_FILE} > peve_GFFannotation.CDS.gff
grep $'\tmRNA\t' ${GFF_FILE} > peve_GFFannotation.mRNA.gff
grep $'\tUTR\t' ${GFF_FILE} > peve_GFFannotation.UTR.gff
```

Extract scaffold lengths 

```
cat is ${genome} | awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > peve.Chromosome_lenghts.txt

wc -l peve.Chromosome_lenghts.txt 
8186 peve.Chromosome_lenghts.txt
```

Extract scaffold names 

```
awk -F" " '{print $1}' peve.Chromosome_lenghts.txt > peve.Chromosome_names.txt
```

I'm going to combine all of the bed stuff into one script. In the e5 scripts folder: `nano bed_peve.sh`

```
#!/bin/bash 
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=5
#SBATCH --export=NONE
#SBATCH --mem=250GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/e5/scripts
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

module load BEDTools/2.30.0-GCC-11.3.0

cd /data/putnamlab/jillashey/genome/Peve

echo "Sorting gffs by chromosome" $(date)

sortBed -faidx peve.Chromosome_names.txt -i peve_GFFannotation.CDS.gff > peve_GFFannotation.CDS_sorted.gff
sortBed -faidx peve.Chromosome_names.txt -i peve_GFFannotation.mRNA.gff > peve_GFFannotation.mRNA_sorted.gff
sortBed -faidx peve.Chromosome_names.txt -i peve_GFFannotation.UTR.gff > peve_GFFannotation.UTR_sorted.gff

echo "Sorting complete!" $(date)

echo "Extracting 1kb 3' UTRs" $(date)

bedtools flank -i peve_GFFannotation.mRNA_sorted.gff -g peve.Chromosome_lenghts.txt -l 0 -r 1000 -s | awk '{gsub("gene","3prime_UTR",$3); print $0 }' | awk '{if($5-$4 > 3)print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' | tr ' ' '\t' > peve.GFFannotation.3UTR_1kb.gff

echo "Subtract portions of UTRs that overlap nearby genes" $(date)

bedtools subtract -a peve.GFFannotation.3UTR_1kb.gff -b peve_GFFannotation.mRNA_sorted.gff > peve.GFFannotation.3UTR_1kb_corrected.gff 

echo "3' UTRs identified!" $(date)

echo "Extracting 3' UTR sequences" $(date)

bedtools getfasta -fi Porites_evermanni_v1.fa -bed peve.GFFannotation.3UTR_1kb_corrected.gff -fo peve_3UTR_1kb.fasta -fullHeader

echo "Sequence extraction complete!" $(date)
```

Submitted batch job 340834. Done! Now run miranda for the Peve data. In the scripts folder: `nano miranda_strict_all_1kb_peve.sh`

```
#!/bin/bash -i
#SBATCH -t 200:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --export=NONE
#SBATCH --mem=500GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/e5/scripts
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

echo "PEVE starting miranda run with all genes and miRNAs with score cutoff >100, energy cutoff <-10, and strict binding invoked"$(date)

module load Miniconda3/4.9.2
conda activate /data/putnamlab/conda/miranda 

miranda /data/putnamlab/jillashey/e5/mirna_seqs/Peve_results_mature_named.fasta /data/putnamlab/jillashey/genome/Peve/peve_3UTR_1kb.fasta -sc 100 -en -10 -strict -out /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_1kb_peve.tab

conda deactivate

echo "miranda run finished!"$(date)
echo "counting number of putative interactions predicted" $(date)

zgrep -c "Performing Scan" /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_1kb_peve.tab

echo "Parsing output" $(date)
grep -A 1 "Scores for this hit:" /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_1kb_peve.tab | sort | grep '>' > /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_1kb_parsed_peve.txt

echo "counting number of putative interactions predicted" $(date)
wc -l /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_1kb_parsed_peve.txt

echo "PEVE miranda script complete" $(date)
```

Submitted batch job 340835. Ran in about 1.5 hours. 

```
head miranda_strict_all_1kb_parsed_peve.txt
>peve-mir-100	Porites_evermani_scaffold_1005:62085-63085	144.00	-16.87	2 17	1 16	15	66.67%	73.33%
>peve-mir-100	Porites_evermani_scaffold_100:97079-98079	140.00	-16.46	2 9	376 395	7	100.00%	100.00%
>peve-mir-100	Porites_evermani_scaffold_101:185725-186725	142.00	-14.24	2 11	850 869	9	88.89%	88.89%
>peve-mir-100	Porites_evermani_scaffold_1027:48934-49934	140.00	-17.00	2 9	543 562	7	100.00%	100.00%
>peve-mir-100	Porites_evermani_scaffold_104:288813-289813	146.00	-19.73	2 19	931 947	17	70.59%	76.47%
>peve-mir-100	Porites_evermani_scaffold_1047:141562-142562	142.00	-16.36	2 11	148 167	9	88.89%	88.89%
>peve-mir-100	Porites_evermani_scaffold_1050:86266-86716	145.00	-16.55	2 10	333 352	8	100.00%	100.00%
>peve-mir-100	Porites_evermani_scaffold_106:362578-363578	144.00	-17.12	2 15	221 239	13	76.92%	84.62%
>peve-mir-100	Porites_evermani_scaffold_1066:63731-64731	150.00	-17.45	2 18	600 617	16	75.00%	81.25%
>peve-mir-100	Porites_evermani_scaffold_1082:85736-86736	155.00	-21.77	2 16	964 983	14	78.57%	85.71%

wc -l miranda_strict_all_1kb_parsed_peve.txt
97782 miranda_strict_all_1kb_parsed_peve.txt
``` 

Redoing Peve as of 10/30/24 - adding gene ids to 3'UTR fasta and more stringent cutoffs for miranda 

```
cd /data/putnamlab/jillashey/genome/Peve
module load BEDTools/2.30.0-GCC-11.3.0

awk '{print $1 "\t" $4-1 "\t" $5 "\t" $9 "\t" "." "\t" $7}' peve.GFFannotation.3UTR_1kb_corrected.gff | sed 's/"//g' > peve.GFFannotation.3UTR_1kb_corrected.bed

bedtools getfasta -fi Porites_evermanni_v1.fa -bed peve.GFFannotation.3UTR_1kb_corrected.bed -fo peve_3UTR_1kb.fasta -name
```

Rerun miranda with `-en -20 -strict` and the rest as defaults. Submitted batch job 346029

```
counting number of putative interactions predicted Wed Oct 30 12:23:23 EDT 2024
1789400
Parsing output Wed Oct 30 12:23:42 EDT 2024
counting number of putative interactions predicted Wed Oct 30 12:23:42 EDT 2024
7187 /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_1kb_parsed_peve.txt
```

### Pmea 

aka Ptuahensis 

To run miranda, I need to identify the 3' UTR ends.

First, identify counts of each feature from gff file

```
cd /data/putnamlab/jillashey/genome/Pmea

grep -v '^#' Pocillopora_meandrina_HIv1.genes_sorted.gff3 | cut -s -f 3 | sort | uniq -c | sort -rn > all_features_pmea.txt

 208535 exon
 208535 CDS
  31840 transcript
```

Generate individual gff for transcripts. 

```
grep $'\ttranscript\t' Pocillopora_meandrina_HIv1.genes_sorted.gff3 > Pmea_transcript.gtf
```

Extract scaffold lengths 

```
cat is Pocillopora_meandrina_HIv1.assembly.fasta | awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > Pmea.Chromosome_lengths.txt

wc -l Pmea.Chromosome_lengths.txt 
212 Pmea.Chromosome_lengths.txt
```

Extract scaffold names 

```
awk -F" " '{print $1}' Pmea.Chromosome_lengths.txt > Pmea.Chromosome_names.txt
```

Create 1kb 3'UTR in gff 

```
interactive

module load BEDTools/2.30.0-GCC-11.3.0

bedtools flank -i Pmea_transcript.gtf -g Pmea.Chromosome_lengths.txt -l 0 -r 1000 -s | awk '{gsub("transcript","3prime_UTR",$3); print $0 }' | awk '{if($5-$4 > 3)print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' | tr ' ' '\t' > Pmea_3UTR_1kb.gtf
```

Subtract portions of 3'UTRs that overlap nearby genes 

```
bedtools subtract -a Pmea_3UTR_1kb.gtf -b Pmea_transcript.gtf > Pmea_3UTR_1kb_corrected.gtf
```

Extract 3'UTR sequences from genome 

```
awk '{print $1 "\t" $4-1 "\t" $5 "\t" $9 "\t" "." "\t" $7}' Pmea_3UTR_1kb_corrected.gtf | sed 's/"//g' > Pmea_3UTR_1kb_corrected.bed

bedtools getfasta -fi Pocillopora_meandrina_HIv1.assembly.fasta -bed Pmea_3UTR_1kb_corrected.bed -fo Pmea_3UTR_1kb.fasta -name
```

Run miranda for Pmea data. In the scripts folder: `nano miranda_strict_all_1kb_pmea.sh`

```
#!/bin/bash -i
#SBATCH -t 48:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --export=NONE
#SBATCH --mem=500GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/e5/scripts
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

echo "PMEA aka PTUH starting miranda run with all genes and miRNAs with energy cutoff <-20 and strict binding invoked"$(date)

module load Miniconda3/4.9.2
conda activate /data/putnamlab/conda/miranda 

miranda /data/putnamlab/jillashey/e5/mirna_seqs/Ptuh_results_mature_named.fasta /data/putnamlab/jillashey/genome/Pmea/Pmea_3UTR_1kb.fasta -en -20 -strict -out /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_1kb_ptuh.tab

conda deactivate

echo "miranda run finished!"$(date)
echo "counting number of putative interactions predicted" $(date)

zgrep -c "Performing Scan" /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_1kb_ptuh.tab

echo "Parsing output" $(date)
grep -A 1 "Scores for this hit:" /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_1kb_ptuh.tab | sort | grep '>' > /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_1kb_parsed_ptuh.txt

echo "counting number of putative interactions predicted" $(date)
wc -l /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_1kb_parsed_ptuh.txt

echo "PMEA aka PTUH miranda script complete" $(date)
```

Submitted batch job 346032. Results:

```
counting number of putative interactions predicted Wed Oct 30 12:42:23 EDT 2024
1204979
Parsing output Wed Oct 30 12:42:34 EDT 2024
counting number of putative interactions predicted Wed Oct 30 12:42:35 EDT 2024
3863 /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_1kb_parsed_ptuh.txt
```
