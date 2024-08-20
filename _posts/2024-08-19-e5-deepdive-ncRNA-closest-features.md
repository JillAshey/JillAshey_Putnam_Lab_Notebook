---
layout: post
title: e5 closest
date: '2024-08-19'
categories: Analysis
tags: [Bioinformatics]
projects: e5 deep dive
---

## e5 deep dive - investigating ncRNA protein machinery in 3 species of corals 

The e5 deep dive project is examining ncRNA dynamics in 3 species of coral from Moorea, French Polynesia. The github for the project is [here](https://github.com/urol-e5/deep-dive). 

In this post, I will be assessing what genomic features are closest to the ncRNAs that we have identified through the deep dive. I will need gffs and fasta files for this task. 

##### Important miRNA files

- [Apul miRNA fasta](https://github.com/urol-e5/deep-dive/blob/main/D-Apul/output/13.2.1-Apul-sRNAseq-ShortStack-31bp-fastp-merged-cnidarian_miRBase/ShortStack_out/mir.fasta); [Apul miRNA gff](https://github.com/urol-e5/deep-dive/blob/main/D-Apul/output/13.2.1-Apul-sRNAseq-ShortStack-31bp-fastp-merged-cnidarian_miRBase/ShortStack_out/Results.gff3)
- [Ptuh miRNA fasta](https://github.com/urol-e5/deep-dive/blob/main/F-Pmea/output/13.2.1-Pmea-sRNAseq-ShortStack-31bp-fastp-merged-cnidarian_miRBase/ShortStack_out/mir.fasta); [Ptuh miRNA gff](https://github.com/urol-e5/deep-dive/blob/main/F-Pmea/output/13.2.1-Pmea-sRNAseq-ShortStack-31bp-fastp-merged-cnidarian_miRBase/ShortStack_out/Results.gff3)
- [Peve miRNA fasta](https://github.com/urol-e5/deep-dive/blob/main/E-Peve/output/08.2-Peve-sRNAseq-ShortStack-31bp-fastp-merged/ShortStack_out/mir.fasta); [Peve miRNA gff](https://github.com/urol-e5/deep-dive/blob/main/E-Peve/output/08.2-Peve-sRNAseq-ShortStack-31bp-fastp-merged/ShortStack_out/Results.gff3)

Because I don't own the deep dive repo, I copied all needed files onto my local computer and then put them on Andromeda. I also added the species prefix to each file. 

```
cd /data/putnamlab/jillashey/e5
mkdir ncRNA_gff
cd ncRNA_gff
```

### Apul miRNA

Sort the genome gff file for Apul. We are using the Amil genome for reference. 

```
cd /data/putnamlab/jillashey/genome/Amil_v2.01/
wget http://gannet.fish.washington.edu/seashell/snaps/GCF_013753865.1_Amil_v2.1_genomic.gff
sort -k1,1 -k4,4n GCF_013753865.1_Amil_v2.1_genomic.gff > GCF_013753865.1_Amil_v2.1_genomic_sorted.gff 
```

Sort miRNA gff file and run bed closest.

```
cd /data/putnamlab/jillashey/e5/ncRNA_gff

sort -k1,1 -k4,4n Apul_Results.gff3 > Apul_Results_sorted.gff3

interactive 
module load BEDTools/2.30.0-GCC-11.3.0
bedtools closest -a Apul_Results_sorted.gff3 -b /data/putnamlab/jillashey/genome/Amil_v2.01/GCF_013753865.1_Amil_v2.1_genomic_sorted.gff > Apul_output.bed

wc -l Apul_output.bed 
61599 Apul_output.bed

head Apul_output.bed 
NC_058066.1	ShortStack	Unknown_sRNA_locus	152483	152910	140	-	.	ID=Cluster_1;DicerCall=N;MIRNA=N	NC_058066.1	RefSeq	region	1	39361238	.	+	.	ID=NC_058066.1:1..39361238;Dbxref=taxon:45264;Name=1;chromosome=1;collection-date=2017;country=Indonesia;gbkey=Src;genome=chromosome;isolate=JS-1;isolation-source=Whole tissue;mol_type=genomic DNA;tissue-type=Adult tissue
NC_058066.1	ShortStack	Unknown_sRNA_locus	152483	152910	140	-	.	ID=Cluster_1;DicerCall=N;MIRNA=N	NC_058066.1	Gnomon	gene	92732	195229	.	+	.	ID=gene-LOC114963509;Dbxref=GeneID:114963509;Name=LOC114963509;gbkey=Gene;gene=LOC114963509;gene_biotype=protein_coding
NC_058066.1	ShortStack	Unknown_sRNA_locus	152483	152910	140	-	.	ID=Cluster_1;DicerCall=N;MIRNA=N	NC_058066.1	Gnomon	mRNA	92732	195229	.	+	.	ID=rna-XM_044317725.1;Parent=gene-LOC114963509;Dbxref=GeneID:114963509,Genbank:XM_044317725.1;Name=XM_044317725.1;Note=The sequence of the model RefSeq transcript was modified relative to this genomic sequence to represent the inferred CDS: deleted 1 base in 1 codon;exception=unclassified transcription discrepancy;experiment=COORDINATES: polyA evidence [ECO:0006239];gbkey=mRNA;gene=LOC114963509;model_evidence=Supporting evidence includes similarity to: 2 mRNAs%2C 3 ESTs%2C 12 Proteins%2C and 100%25 coverage of the annotated genomic feature by RNAseq alignments%2C including 1 sample with support for all annotated introns;product=uncharacterized LOC114963509;transcript_id=XM_044317725.1
NC_058066.1	ShortStack	Unknown_sRNA_locus	152483	152910	140	-	.	ID=Cluster_1;DicerCall=N;MIRNA=N	NC_058066.1	Gnomon	gene	145277	165521	.	+	.	ID=gene-LOC122957574;Dbxref=GeneID:122957574;Name=LOC122957574;gbkey=Gene;gene=LOC122957574;gene_biotype=protein_coding
NC_058066.1	ShortStack	Unknown_sRNA_locus	152483	152910	140	-	.	ID=Cluster_1;DicerCall=N;MIRNA=N	NC_058066.1	Gnomon	mRNA	145277	165521	.	+	.	ID=rna-XM_044317744.1;Parent=gene-LOC122957574;Dbxref=GeneID:122957574,Genbank:XM_044317744.1;Name=XM_044317744.1;gbkey=mRNA;gene=LOC122957574;model_evidence=Supporting evidence includes similarity to: 5 Proteins%2C and 99%25 coverage of the annotated genomic feature by RNAseq alignments%2C including 2 samples with support for all annotated introns;product=uncharacterized LOC122957574;transcript_id=XM_044317744.1
NC_058066.1	ShortStack	Unknown_sRNA_locus	152483	152910	140	-	.	ID=Cluster_1;DicerCall=N;MIRNA=N	NC_058066.1	Gnomon	CDS	152268	152501	.	+	2	ID=cds-XP_044173679.1;Parent=rna-XM_044317744.1;Dbxref=GeneID:122957574,Genbank:XP_044173679.1;Name=XP_044173679.1;gbkey=CDS;gene=LOC122957574;product=uncharacterized protein LOC122957574;protein_id=XP_044173679.1
NC_058066.1	ShortStack	Unknown_sRNA_locus	152483	152910	140	-	.	ID=Cluster_1;DicerCall=N;MIRNA=N	NC_058066.1	Gnomon	exon	152268	152501	.	+	.	ID=exon-XM_044317744.1-4;Parent=rna-XM_044317744.1;Dbxref=GeneID:122957574,Genbank:XM_044317744.1;gbkey=mRNA;gene=LOC122957574;product=uncharacterized LOC122957574;transcript_id=XM_044317744.1
NC_058066.1	ShortStack	Unknown_sRNA_locus	161064	161674	549	.	.	ID=Cluster_2;DicerCall=N;MIRNA=N	NC_058066.1	RefSeq	region	1	39361238	.	+	.	ID=NC_058066.1:1..39361238;Dbxref=taxon:45264;Name=1;chromosome=1;collection-date=2017;country=Indonesia;gbkey=Src;genome=chromosome;isolate=JS-1;isolation-source=Whole tissue;mol_type=genomic DNA;tissue-type=Adult tissue
NC_058066.1	ShortStack	Unknown_sRNA_locus	161064	161674	549	.	.	ID=Cluster_2;DicerCall=N;MIRNA=N	NC_058066.1	Gnomon	gene	92732	195229	.	+	.	ID=gene-LOC114963509;Dbxref=GeneID:114963509;Name=LOC114963509;gbkey=Gene;gene=LOC114963509;gene_biotype=protein_coding
NC_058066.1	ShortStack	Unknown_sRNA_locus	161064	161674	549	.	.	ID=Cluster_2;DicerCall=N;MIRNA=N	NC_058066.1	Gnomon	mRNA	92732	195229	.	+	.	ID=rna-XM_044317725.1;Parent=gene-LOC114963509;Dbxref=GeneID:114963509,Genbank:XM_044317725.1;Name=XM_044317725.1;Note=The sequence of the model RefSeq transcript was modified relative to this genomic sequence to represent the inferred CDS: deleted 1 base in 1 codon;exception=unclassified transcription discrepancy;experiment=COORDINATES: polyA evidence [ECO:0006239];gbkey=mRNA;gene=LOC114963509;model_evidence=Supporting evidence includes similarity to: 2 mRNAs%2C 3 ESTs%2C 12 Proteins%2C and 100%25 coverage of the annotated genomic feature by RNAseq alignments%2C including 1 sample with support for all annotated introns;product=uncharacterized LOC114963509;transcript_id=XM_044317725.1
```

Remove the unknown siRNA loci

```
awk 'BEGIN {OFS="\t"} $3 != "Unknown_sRNA_locus"' Apul_output.bed > filtered_Apul_output.bed

wc -l filtered_Apul_output.bed
755 filtered_Apul_output.bed

head filtered_Apul_output.bed 
NC_058066.1	ShortStack	siRNA24_locus	5224792	5225215	1264	+	.	ID=Cluster_184;DicerCall=24;MIRNA=N	NC_058066.1	RefSeq	region	1	39361238	.	+	.	ID=NC_058066.1:1..39361238;Dbxref=taxon:45264;Name=1;chromosome=1;collection-date=2017;country=Indonesia;gbkey=Src;genome=chromosome;isolate=JS-1;isolation-source=Whole tissue;mol_type=genomic DNA;tissue-type=Adult tissue
NC_058066.1	ShortStack	siRNA24_locus	5224792	5225215	1264	+	.	ID=Cluster_184;DicerCall=24;MIRNA=N	NC_058066.1	Gnomon	gene	5153290	5231353	.	-	.	ID=gene-LOC114950433;Dbxref=GeneID:114950433;Name=LOC114950433;gbkey=Gene;gene=LOC114950433;gene_biotype=protein_coding
NC_058066.1	ShortStack	siRNA24_locus	5224792	5225215	1264	+	.	ID=Cluster_184;DicerCall=24;MIRNA=N	NC_058066.1	Gnomon	mRNA	5153290	5231353	.	-	.	ID=rna-XM_044310280.1;Parent=gene-LOC114950433;Dbxref=GeneID:114950433,Genbank:XM_044310280.1;Name=XM_044310280.1;experiment=COORDINATES: polyA evidence [ECO:0006239];gbkey=mRNA;gene=LOC114950433;model_evidence=Supporting evidence includes similarity to: 15 mRNAs%2C 31 Proteins%2C and 99%25 coverage of the annotated genomic feature by RNAseq alignments;product=uncharacterized LOC114950433;transcript_id=XM_044310280.1
NC_058066.1	ShortStack	siRNA21_locus	7563377	7563797	118	+	.	ID=Cluster_225;DicerCall=21;MIRNA=N	NC_058066.1	RefSeq	region	1	39361238	.	+	.	ID=NC_058066.1:1..39361238;Dbxref=taxon:45264;Name=1;chromosome=1;collection-date=2017;country=Indonesia;gbkey=Src;genome=chromosome;isolate=JS-1;isolation-source=Whole tissue;mol_type=genomic DNA;tissue-type=Adult tissue
NC_058066.1	ShortStack	siRNA21_locus	7563377	7563797	118	+	.	ID=Cluster_225;DicerCall=21;MIRNA=N	NC_058066.1	Gnomon	gene	7523354	7569602	.	-	.	ID=gene-LOC114970982;Dbxref=GeneID:114970982;Name=LOC114970982;gbkey=Gene;gene=LOC114970982;gene_biotype=protein_coding
NC_058066.1	ShortStack	siRNA21_locus	7563377	7563797	118	+	.	ID=Cluster_225;DicerCall=21;MIRNA=N	NC_058066.1	Gnomon	mRNA	7523354	7569602	.	-	.	ID=rna-XM_044320669.1;Parent=gene-LOC114970982;Dbxref=GeneID:114970982,Genbank:XM_044320669.1;Name=XM_044320669.1;Note=The sequence of the model RefSeq transcript was modified relative to this genomic sequence to represent the inferred CDS: inserted 1 base in 1 codon;exception=unclassified transcription discrepancy;gbkey=mRNA;gene=LOC114970982;model_evidence=Supporting evidence includes similarity to: 6 Proteins%2C and 98%25 coverage of the annotated genomic feature by RNAseq alignments%2C including 1 sample with support for all annotated introns;product=uncharacterized LOC114970982;transcript_id=XM_044320669.1
NC_058066.1	ShortStack	siRNA22_locus	8905068	8905484	102	-	.	ID=Cluster_251;DicerCall=22;MIRNA=N	NC_058066.1	RefSeq	region	1	39361238	.	+	.	ID=NC_058066.1:1..39361238;Dbxref=taxon:45264;Name=1;chromosome=1;collection-date=2017;country=Indonesia;gbkey=Src;genome=chromosome;isolate=JS-1;isolation-source=Whole tissue;mol_type=genomic DNA;tissue-type=Adult tissue
NC_058066.1	ShortStack	MIRNA_hairpin	12757125	12757218	8293	-	.	ID=Cluster_316;DicerCall=23;MIRNA=Y	NC_058066.1	RefSeq	region	1	39361238	.	+	.	ID=NC_058066.1:1..39361238;Dbxref=taxon:45264;Name=1;chromosome=1;collection-date=2017;country=Indonesia;gbkey=Src;genome=chromosome;isolate=JS-1;isolation-source=Whole tissue;mol_type=genomic DNA;tissue-type=Adult tissue
NC_058066.1	ShortStack	MIRNA_hairpin	12757125	12757218	8293	-	.	ID=Cluster_316;DicerCall=23;MIRNA=Y	NC_058066.1	Gnomon	gene	12755159	12764546	.	-	.	ID=gene-LOC114961148;Dbxref=GeneID:114961148;Name=LOC114961148;gbkey=Gene;gene=LOC114961148;gene_biotype=protein_coding
NC_058066.1	ShortStack	MIRNA_hairpin	12757125	12757218	8293	-	.	ID=Cluster_316;DicerCall=23;MIRNA=Y	NC_058066.1	Gnomon	mRNA	12755159	12764546	.	-	.	ID=rna-XM_029339755.2;Parent=gene-LOC114961148;Dbxref=GeneID:114961148,Genbank:XM_029339755.2;Name=XM_029339755.2;experiment=COORDINATES: polyA evidence [ECO:0006239];gbkey=mRNA;gene=LOC114961148;model_evidence=Supporting evidence includes similarity to: 7 mRNAs%2C 8 Proteins%2C and 100%25 coverage of the annotated genomic feature by RNAseq alignments%2C including 87 samples with support for all annotated introns;product=zinc finger MYND domain-containing protein 19-like;transcript_id=XM_029339755.2
```

If I subtract the 5th column (end of miRNA feature) by the 13th column (beginning of genomic feature), this will tell me if there is overlap in the two features. If the number is positive, it means there is overlap. If it is negative, it means there is no overlap. I did not run bed intersect because I attempted that with Ptuh in code below and it did not give me any meaningful output. 

```
awk '{ $(NF+1) = $13 - $5 }1' OFS='\t' filtered_Apul_output.bed > filtered_Apul_output_overlap.bed

head filtered_Apul_output_overlap.bed 
NC_058066.1	ShortStack	siRNA24_locus	5224792	5225215	1264	+	.	ID=Cluster_184;DicerCall=24;MIRNA=N	NC_058066.1	RefSeq	region	1	39361238	.	+	.	ID=NC_058066.1:1..39361238;Dbxref=taxon:45264;Name=1;chromosome=1;collection-date=2017;country=Indonesia;gbkey=Src;genome=chromosome;isolate=JS-1;isolation-source=Whole	tissue;mol_type=genomic	DNA;tissue-type=Adult	tissue	-5225214
NC_058066.1	ShortStack	siRNA24_locus	5224792	5225215	1264	+	.	ID=Cluster_184;DicerCall=24;MIRNA=N	NC_058066.1	Gnomon	gene	5153290	5231353	.	-	.	ID=gene-LOC114950433;Dbxref=GeneID:114950433;Name=LOC114950433;gbkey=Gene;gene=LOC114950433;gene_biotype=protein_coding	-71925
NC_058066.1	ShortStack	siRNA24_locus	5224792	5225215	1264	+	.	ID=Cluster_184;DicerCall=24;MIRNA=N	NC_058066.1	Gnomon	mRNA	5153290	5231353	.	-	.	ID=rna-XM_044310280.1;Parent=gene-LOC114950433;Dbxref=GeneID:114950433,Genbank:XM_044310280.1;Name=XM_044310280.1;experiment=COORDINATES:	polyA	evidence	[ECO:0006239];gbkey=mRNA;gene=LOC114950433;model_evidence=Supporting	evidence	includes	similarity	to:	15	mRNAs%2C	31	Proteins%2C	and	99%25	coverage	of	the	annotated	genomic	feature	by	RNAseq	alignments;product=uncharacterized	LOC114950433;transcript_id=XM_044310280.1	-71925
NC_058066.1	ShortStack	siRNA21_locus	7563377	7563797	118	+	.	ID=Cluster_225;DicerCall=21;MIRNA=N	NC_058066.1	RefSeq	region	1	39361238	.	+	.	ID=NC_058066.1:1..39361238;Dbxref=taxon:45264;Name=1;chromosome=1;collection-date=2017;country=Indonesia;gbkey=Src;genome=chromosome;isolate=JS-1;isolation-source=Whole	tissue;mol_type=genomic	DNA;tissue-type=Adult	tissue	-7563796
NC_058066.1	ShortStack	siRNA21_locus	7563377	7563797	118	+	.	ID=Cluster_225;DicerCall=21;MIRNA=N	NC_058066.1	Gnomon	gene	7523354	7569602	.	-	.	ID=gene-LOC114970982;Dbxref=GeneID:114970982;Name=LOC114970982;gbkey=Gene;gene=LOC114970982;gene_biotype=protein_coding	-40443
NC_058066.1	ShortStack	siRNA21_locus	7563377	7563797	118	+	.	ID=Cluster_225;DicerCall=21;MIRNA=N	NC_058066.1	Gnomon	mRNA	7523354	7569602	.	-	.	ID=rna-XM_044320669.1;Parent=gene-LOC114970982;Dbxref=GeneID:114970982,Genbank:XM_044320669.1;Name=XM_044320669.1;Note=The	sequence	of	the	model	RefSeq	transcript	was	modified	relative	to	this	genomic	sequence	to	represent	the	inferred	CDS:	inserted	base	in	1	codon;exception=unclassified	transcription	discrepancy;gbkey=mRNA;gene=LOC114970982;model_evidence=Supporting	evidence	includes	similarity	to:	6	Proteins%2C	and	98%25	coverage	of	theannotated	genomic	feature	by	RNAseq	alignments%2C	including	1	sample	with	support	for	all	annotated	introns;product=uncharacterized	LOC114970982;transcript_id=XM_044320669.1	-40443
NC_058066.1	ShortStack	siRNA22_locus	8905068	8905484	102	-	.	ID=Cluster_251;DicerCall=22;MIRNA=N	NC_058066.1	RefSeq	region	1	39361238	.	+	.	ID=NC_058066.1:1..39361238;Dbxref=taxon:45264;Name=1;chromosome=1;collection-date=2017;country=Indonesia;gbkey=Src;genome=chromosome;isolate=JS-1;isolation-source=Whole	tissue;mol_type=genomic	DNA;tissue-type=Adult	tissue	-8905483
NC_058066.1	ShortStack	MIRNA_hairpin	12757125	12757218	8293	-	.	ID=Cluster_316;DicerCall=23;MIRNA=Y	NC_058066.1	RefSeq	region	1	39361238	.	+	.	ID=NC_058066.1:1..39361238;Dbxref=taxon:45264;Name=1;chromosome=1;collection-date=2017;country=Indonesia;gbkey=Src;genome=chromosome;isolate=JS-1;isolation-source=Whole	tissue;mol_type=genomic	DNA;tissue-type=Adult	tissue	-12757217
NC_058066.1	ShortStack	MIRNA_hairpin	12757125	12757218	8293	-	.	ID=Cluster_316;DicerCall=23;MIRNA=Y	NC_058066.1	Gnomon	gene	12755159	12764546	.	-	.	ID=gene-LOC114961148;Dbxref=GeneID:114961148;Name=LOC114961148;gbkey=Gene;gene=LOC114961148;gene_biotype=protein_coding	-2059
NC_058066.1	ShortStack	MIRNA_hairpin	12757125	12757218	8293	-	.	ID=Cluster_316;DicerCall=23;MIRNA=Y	NC_058066.1	Gnomon	mRNA	12755159	12764546	.	-	.	ID=rna-XM_029339755.2;Parent=gene-LOC114961148;Dbxref=GeneID:114961148,Genbank:XM_029339755.2;Name=XM_029339755.2;experiment=COORDINATES:	polyA	evidence	[ECO:0006239];gbkey=mRNA;gene=LOC114961148;model_evidence=Supporting	evidence	includes	similarity	to:	7	mRNAs%2C	8	Proteins%2C	and	100%25	coverage	of	the	annotated	genomic	feature	by	RNAseq	alignments%2C	including	87	samples	with	support	for	all	annotated	introns;product=zinc	finger	MYNdomain-containing	protein	19-like;transcript_id=XM_029339755.2	-2059
```

Count number of positive and negative numbers

```
awk '{if ($NF < 0) neg++; else if ($NF > 0) pos++} END {print "Positive count:", pos; print "Negative count:", neg}' filtered_Apul_output_overlap.bed

Positive count: 
Negative count: 755
```

Negative = no overlap though I should check my math with someone to make sure I'm not dumb. But looks like no overlap between the miRNA features and the genomic features. 

### Peve miRNA

Sort Peve genome gff file

```
cd /data/putnamlab/jillashey/genome/Peve/
sort -k1,1 -k4,4n Porites_evermanni_v1.annot.gff > Porites_evermanni_v1.annot_sorted.gff
``` 

Sort miRNA gff file and run `bed closest`. 

```
cd /data/putnamlab/jillashey/e5/ncRNA_gff

sort -k1,1 -k4,4n Peve_Results.gff3 > Peve_Results_sorted.gff3

interactive 
module load BEDTools/2.30.0-GCC-11.3.0
bedtools closest -a Peve_Results_sorted.gff3 -b /data/putnamlab/jillashey/genome/Peve/Porites_evermanni_v1.annot_sorted.gff > Peve_output.bed

wc -l Peve_output.bed
32806 Peve_output.bed

head Peve_output.bed
Porites_evermani_scaffold_1	ShortStack	Unknown_sRNA_locus	45711	46131	88	+	.	ID=Cluster_1;DicerCall=N;MIRNA=N	Porites_evermani_scaffold_1	Gmove	mRNA	32616	67628	399	-	.	ID=Peve_00000122;Name=Peve_00000122;start=1;stop=1;cds_size=399
Porites_evermani_scaffold_1	ShortStack	Unknown_sRNA_locus	201507	201931	58	-	.	ID=Cluster_2;DicerCall=N;MIRNA=N	Porites_evermani_scaffold_1	Gmove	CDS	205241	208000	.	+	.	Parent=Peve_00000106
Porites_evermani_scaffold_1	ShortStack	Unknown_sRNA_locus	201507	201931	58	-	.	ID=Cluster_2;DicerCall=N;MIRNA=N	Porites_evermani_scaffold_1	Gmove	mRNA	205241	208000	276	+	.	ID=Peve_00000106;Name=Peve_00000106;start=1;stop=1;cds_size=2760
Porites_evermani_scaffold_1	ShortStack	Unknown_sRNA_locus	313446	313846	50	-	.	ID=Cluster_3;DicerCall=N;MIRNA=N	Porites_evermani_scaffold_1	Gmove	mRNA	307343	313927	924.6	-	.	ID=Peve_00000114;Name=Peve_00000114;start=1;stop=1;cds_size=1206
Porites_evermani_scaffold_1	ShortStack	Unknown_sRNA_locus	313446	313846	50	-	.	ID=Cluster_3;DicerCall=N;MIRNA=N	Porites_evermani_scaffold_1	Gmove	UTR	313287	313927	.	-	.	Parent=Peve_00000114
Porites_evermani_scaffold_1	ShortStack	Unknown_sRNA_locus	406146	406734	175	-	.	ID=Cluster_4;DicerCall=N;MIRNA=N	Porites_evermani_scaffold_1	Gmove	mRNA	384175	413351	1590.6	-	.	ID=Peve_00000121;Name=Peve_00000121;start=1;stop=1;cds_size=1860
Porites_evermani_scaffold_1	ShortStack	Unknown_sRNA_locus	409839	410269	169	-	.	ID=Cluster_5;DicerCall=N;MIRNA=N	Porites_evermani_scaffold_1	Gmove	mRNA	384175	413351	1590.6	-	.	ID=Peve_00000121;Name=Peve_00000121;start=1;stop=1;cds_size=1860
Porites_evermani_scaffold_1	ShortStack	Unknown_sRNA_locus	465244	465668	169	-	.	ID=Cluster_6;DicerCall=N;MIRNA=N	Porites_evermani_scaffold_1	Gmove	mRNA	462457	477071	1669.32	-	.	ID=Peve_00000006;Name=Peve_00000006;start=1;stop=1;cds_size=2034
Porites_evermani_scaffold_1	ShortStack	Unknown_sRNA_locus	465244	465668	169	-	.	ID=Cluster_6;DicerCall=N;MIRNA=N	Porites_evermani_scaffold_1	Gmove	CDS	465272	465508	.	-	.	Parent=Peve_00000006
Porites_evermani_scaffold_1	ShortStack	Unknown_sRNA_locus	468473	468950	91900	-	.	ID=Cluster_7;DicerCall=N;MIRNA=N	Porites_evermani_scaffold_1	Gmove	mRNA	462457	477071	1669.32	-	.	ID=Peve_00000006;Name=Peve_00000006;start=1;stop=1;cds_size=2034
```

Remove anything that is unknown 

```
awk 'BEGIN {OFS="\t"} $3 != "Unknown_sRNA_locus"' Peve_output.bed > filtered_Peve_output.bed

wc -l filtered_Peve_output.bed
449 filtered_Peve_output.bed

head filtered_Peve_output.bed
Porites_evermani_scaffold_1	ShortStack	MIRNA_hairpin	1404250	1404342	9574	-	.	ID=Cluster_29;DicerCall=N;MIRNA=Y	Porites_evermani_scaffold_1	Gmove	mRNA	1380413	1416448	1270.53	-	.	ID=Peve_00000077;Name=Peve_00000077;start=1;stop=1;cds_size=1851
Porites_evermani_scaffold_1	ShortStack	mature_miRNA	1404272	1404293	3403	-	.	ID=Cluster_29.mature;Parent=Cluster_29	Porites_evermani_scaffold_1	Gmove	mRNA	1380413	1416448	1270.53	-	.	ID=Peve_00000077;Name=Peve_00000077;start=1;stop=1;cds_size=1851
Porites_evermani_scaffold_1	ShortStack	miRNA-star	1404301	1404322	23	-	.	ID=Cluster_29.star;Parent=Cluster_29	Porites_evermani_scaffold_1	Gmove	mRNA	1380413	1416448	1270.53	-	.	ID=Peve_00000077;Name=Peve_00000077;start=1;stop=1;cds_size=1851
Porites_evermani_scaffold_10	ShortStack	siRNA21_locus	565492	565912	76	+	.	ID=Cluster_353;DicerCall=21;MIRNA=N	Porites_evermani_scaffold_10	Gmove	mRNA	562195	564790	948	+	.	ID=Peve_00000127;Name=Peve_00000127;start=1;stop=1;cds_size=474
Porites_evermani_scaffold_10	ShortStack	siRNA21_locus	565492	565912	76	+	.	ID=Cluster_353;DicerCall=21;MIRNA=N	Porites_evermani_scaffold_10	Gmove	UTR	564704	564790	.	+	.	Parent=Peve_00000127
Porites_evermani_scaffold_1005	ShortStack	siRNA23_locus	126975	127397	168	+	.	ID=Cluster_9183;DicerCall=23;MIRNA=N	Porites_evermani_scaffold_1005	Gmove	mRNA	114473	133413	606	+	.	ID=Peve_00000326;Name=Peve_00000326;start=0;stop=1;cds_size=606
Porites_evermani_scaffold_1060	ShortStack	siRNA22_locus	77378	77799	99	+	.	ID=Cluster_9583;DicerCall=22;MIRNA=N	Porites_evermani_scaffold_1060	Gmove	mRNA	64377	75861	817.5	+	.	ID=Peve_00001166;Name=Peve_00001166;start=1;stop=1;cds_size=1275
Porites_evermani_scaffold_1060	ShortStack	siRNA22_locus	77378	77799	99	+	.	ID=Cluster_9583;DicerCall=22;MIRNA=N	Porites_evermani_scaffold_1060	Gmove	CDS	75822	75861	.	+	.	Parent=Peve_00001166
Porites_evermani_scaffold_108	ShortStack	siRNA24_locus	306073	306496	53	-	.	ID=Cluster_2199;DicerCall=24;MIRNA=N	Porites_evermani_scaffold_108	Gmove	mRNA	278980	311091	2217.57	-	.	ID=Peve_00001444;Name=Peve_00001444;start=1;stop=1;cds_size=4026
Porites_evermani_scaffold_108	ShortStack	siRNA24_locus	306073	306496	53	-	.	ID=Cluster_2199;DicerCall=24;MIRNA=N	Porites_evermani_scaffold_108	Gmove	CDS	306249	306300	.	-	.	Parent=Peve_00001444
```

If I subtract the 5th column (end of miRNA feature) by the 13th column (beginning of genomic feature), this will tell me if there is overlap in the two features. If the number is positive, it means there is overlap. If it is negative, it means there is no overlap. I did not run bed intersect because I attempted that with Ptuh in code below and it did not give me any meaningful output.

```
awk '{ $(NF+1) = $13 - $5 }1' OFS='\t' filtered_Peve_output.bed > filtered_Peve_output_overlap.bed

head filtered_Peve_output_overlap.bed
Porites_evermani_scaffold_1	ShortStack	MIRNA_hairpin	1404250	1404342	9574	-	.	ID=Cluster_29;DicerCall=N;MIRNA=Y	Porites_evermani_scaffold_1	Gmove	mRNA	1380413	1416448	1270.53	-	.	ID=Peve_00000077;Name=Peve_00000077;start=1;stop=1;cds_size=1851	-23929
Porites_evermani_scaffold_1	ShortStack	mature_miRNA	1404272	1404293	3403	-	.	ID=Cluster_29.mature;Parent=Cluster_29	Porites_evermani_scaffold_1	Gmove	mRNA	1380413	1416448	1270.53	-	.	ID=Peve_00000077;Name=Peve_00000077;start=1;stop=1;cds_size=1851	-23880
Porites_evermani_scaffold_1	ShortStack	miRNA-star	1404301	1404322	23	-	.	ID=Cluster_29.star;Parent=Cluster_29	Porites_evermani_scaffold_1	Gmove	mRNA	1380413	1416448	1270.53	-	.	ID=Peve_00000077;Name=Peve_00000077;start=1;stop=1;cds_size=1851	-23909
Porites_evermani_scaffold_10	ShortStack	siRNA21_locus	565492	565912	76	+	.	ID=Cluster_353;DicerCall=21;MIRNA=N	Porites_evermani_scaffold_10	Gmove	mRNA	562195	564790	948	+	.	ID=Peve_00000127;Name=Peve_00000127;start=1;stop=1;cds_size=474	-3717
Porites_evermani_scaffold_10	ShortStack	siRNA21_locus	565492	565912	76	+	.	ID=Cluster_353;DicerCall=21;MIRNA=N	Porites_evermani_scaffold_10	Gmove	UTR	564704	564790	.	+	.	Parent=Peve_0000012-1208
Porites_evermani_scaffold_1005	ShortStack	siRNA23_locus	126975	127397	168	+	.	ID=Cluster_9183;DicerCall=23;MIRNA=N	Porites_evermani_scaffold_1005	Gmove	mRNA	114473	133413	606	+	.	ID=Peve_00000326;Name=Peve_00000326;start=0;stop=1;cds_size=606	-12924
Porites_evermani_scaffold_1060	ShortStack	siRNA22_locus	77378	77799	99	+	.	ID=Cluster_9583;DicerCall=22;MIRNA=N	Porites_evermani_scaffold_1060	Gmove	mRNA	64377	75861	817.5	+	.	ID=Peve_00001166;Name=Peve_00001166;start=1;stop=1;cds_size=1275	-13422
Porites_evermani_scaffold_1060	ShortStack	siRNA22_locus	77378	77799	99	+	.	ID=Cluster_9583;DicerCall=22;MIRNA=N	Porites_evermani_scaffold_1060	Gmove	CDS	75822	75861	.	+	.	Parent=Peve_00001166	-1977
Porites_evermani_scaffold_108	ShortStack	siRNA24_locus	306073	306496	53	-	.	ID=Cluster_2199;DicerCall=24;MIRNA=N	Porites_evermani_scaffold_108	Gmove	mRNA	278980	311091	2217.57	-	.	ID=Peve_00001444;Name=Peve_00001444;start=1;stop=1;cds_size=4026	-27516
Porites_evermani_scaffold_108	ShortStack	siRNA24_locus	306073	306496	53	-	.	ID=Cluster_2199;DicerCall=24;MIRNA=N	Porites_evermani_scaffold_108	Gmove	CDS	306249	306300	.	-	.	Parent=Peve_00001444	-247
```

Count number of positive and negative numbers

```
awk '{if ($NF < 0) neg++; else if ($NF > 0) pos++} END {print "Positive count:", pos; print "Negative count:", neg}' filtered_Peve_output_overlap.bed

Positive count: 100
Negative count: 349
```

### Ptuh miRNA

Sort Ptuh genome gff file; in this project, we used the Pmea genome. 

```
cd /data/putnamlab/jillashey/genome/Pmea
sort -k1,1 -k4,4n Pocillopora_meandrina_HIv1.genes.gff3 > Pocillopora_meandrina_HIv1.genes_sorted.gff3
```

Sort miRNA gff file and run `bed closest`. 

```
cd /data/putnamlab/jillashey/e5/ncRNA_gff

sort -k1,1 -k4,4n Ptuh_Results.gff3 > Ptuh_Results_sorted.gff3

interactive 
module load BEDTools/2.30.0-GCC-11.3.0
bedtools closest -a Ptuh_Results_sorted.gff3 -b /data/putnamlab/jillashey/genome/Pmea/Pocillopora_meandrina_HIv1.genes_sorted.gff3 > Ptuh_output.bed
```

Success!!!

```
wc -l Ptuh_output.bed 
21604 Ptuh_output.bed

head Ptuh_output.bed 
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	Unknown_sRNA_locus	9092	9521	10813	+	.	ID=Cluster_1;DicerCall=N;MIRNA=N	Pocillopora_meandrina_HIv1___Sc0000000	AUGUSTUS	CDS	10771	11117	.	+	0	Parent=Pocillopora_meandrina_HIv1___RNAseq.g20902.t1
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	Unknown_sRNA_locus	9092	9521	10813	+	.	ID=Cluster_1;DicerCall=N;MIRNA=N	Pocillopora_meandrina_HIv1___Sc0000000	AUGUSTUS	exon	10771	11117	.	+	0	Parent=Pocillopora_meandrina_HIv1___RNAseq.g20902.t1
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	Unknown_sRNA_locus	9092	9521	10813	+	.	ID=Cluster_1;DicerCall=N;MIRNA=N	Pocillopora_meandrina_HIv1___Sc0000000	AUGUSTUS	transcript	10771	23652	.	+	.	ID=Pocillopora_meandrina_HIv1___RNAseq.g20902.t1
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	Unknown_sRNA_locus	53578	53997	287	+	.	ID=Cluster_2;DicerCall=N;MIRNA=N	Pocillopora_meandrina_HIv1___Sc0000000	AUGUSTUS	transcript	52050	53624	.	-	.	ID=Pocillopora_meandrina_HIv1___RNAseq.g20906.t1
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	Unknown_sRNA_locus	53578	53997	287	+	.	ID=Cluster_2;DicerCall=N;MIRNA=N	Pocillopora_meandrina_HIv1___Sc0000000	AUGUSTUS	CDS	53573	53624	.	-	0	Parent=Pocillopora_meandrina_HIv1___RNAseq.g20906.t1
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	Unknown_sRNA_locus	53578	53997	287	+	.	ID=Cluster_2;DicerCall=N;MIRNA=N	Pocillopora_meandrina_HIv1___Sc0000000	AUGUSTUS	exon	53573	53624	.	-	0	Parent=Pocillopora_meandrina_HIv1___RNAseq.g20906.t1
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	Unknown_sRNA_locus	150243	150718	2549	-	.	ID=Cluster_3;DicerCall=N;MIRNA=N	Pocillopora_meandrina_HIv1___Sc0000000	AUGUSTUS	transcript	143552	155669	.	-	.	ID=Pocillopora_meandrina_HIv1___RNAseq.g20914.t1
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	Unknown_sRNA_locus	150243	150718	2549	-	.	ID=Cluster_3;DicerCall=N;MIRNA=N	Pocillopora_meandrina_HIv1___Sc0000000	AUGUSTUS	CDS	150290	150371	.	-	1	Parent=Pocillopora_meandrina_HIv1___RNAseq.g20914.t1
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	Unknown_sRNA_locus	150243	150718	2549	-	.	ID=Cluster_3;DicerCall=N;MIRNA=N	Pocillopora_meandrina_HIv1___Sc0000000	AUGUSTUS	exon	150290	150371	.	-	1	Parent=Pocillopora_meandrina_HIv1___RNAseq.g20914.t1
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	Unknown_sRNA_locus	150243	150718	2549	-	.	ID=Cluster_3;DicerCall=N;MIRNA=N	Pocillopora_meandrina_HIv1___Sc0000000	AUGUSTUS	CDS	150573	150661	.	-	0	Parent=Pocillopora_meandrina_HIv1___RNAseq.g20914.t1
```

Remove everything that is unknown. 

```
awk 'BEGIN {OFS="\t"} $3 != "Unknown_sRNA_locus"' Ptuh_output.bed > filtered_Ptuh_output.bed
wc -l filtered_Ptuh_output.bed 
405 filtered_Ptuh_output.bed

head filtered_Ptuh_output.bed
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	siRNA22_locus	173728	174150	1257	+	.	ID=Cluster_4;DicerCall=22;MIRNA=N	Pocillopora_meandrina_HIv1___Sc0000000	AUGUSTUS	CDS	174509	175333	.	-	0	Parent=Pocillopora_meandrina_HIv1___RNAseq.g20918.t1
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	siRNA22_locus	173728	174150	1257	+	.	ID=Cluster_4;DicerCall=22;MIRNA=N	Pocillopora_meandrina_HIv1___Sc0000000	AUGUSTUS	exon	174509	175333	.	-	0	Parent=Pocillopora_meandrina_HIv1___RNAseq.g20918.t1
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	siRNA22_locus	173728	174150	1257	+	.	ID=Cluster_4;DicerCall=22;MIRNA=N	Pocillopora_meandrina_HIv1___Sc0000000	AUGUSTUS	transcript	174509	176444	.	-	.	ID=Pocillopora_meandrina_HIv1___RNAseq.g20918.t1
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	MIRNA_hairpin	818027	818120	12096	+	.	ID=Cluster_19;DicerCall=23;MIRNA=Y	Pocillopora_meandrina_HIv1___Sc0000000	AUGUSTUS	transcript	816355	820160	.	+	.	ID=Pocillopora_meandrina_HIv1___RNAseq.g21001.t1
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	mature_miRNA	818049	818070	3240	+	.	ID=Cluster_19.mature;Parent=Cluster_19	Pocillopora_meandrina_HIv1___Sc0000000	AUGUSTUS	transcript	816355	820160	.	+	.	ID=Pocillopora_meandrina_HIv1___RNAseq.g21001.t1
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	miRNA-star	818079	818100	9	+	.	ID=Cluster_19.star;Parent=Cluster_19	Pocillopora_meandrina_HIv1___Sc0000000	AUGUSTUS	transcript	816355	820160	.	+	.	ID=Pocillopora_meandrina_HIv1___RNAseq.g21001.t1
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	MIRNA_hairpin	2872019	2872110	177	+	.	ID=Cluster_34;DicerCall=21;MIRNA=Y	Pocillopora_meandrina_HIv1___Sc0000000	AUGUSTUS	transcript	2868586	2871318	.	+	.	ID=Pocillopora_meandrina_HIv1___TS.g25957.t1
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	MIRNA_hairpin	2872019	2872110	177	+	.	ID=Cluster_34;DicerCall=21;MIRNA=Y	Pocillopora_meandrina_HIv1___Sc0000000	AUGUSTUS	CDS	2870522	2871318	.	+	2	Parent=Pocillopora_meandrina_HIv1___TS.g25957.t1
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	MIRNA_hairpin	2872019	2872110	177	+	.	ID=Cluster_34;DicerCall=21;MIRNA=Y	Pocillopora_meandrina_HIv1___Sc0000000	AUGUSTUS	exon	2870522	2871318	.	+	2	Parent=Pocillopora_meandrina_HIv1___TS.g25957.t1
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	mature_miRNA	2872041	2872061	110	+	.	ID=Cluster_34.mature;Parent=Cluster_34	Pocillopora_meandrina_HIv1___Sc0000000	AUGUSTUS	transcript	2868586	2871318	.	+	.	ID=Pocillopora_meandrina_HIv1___TS.g25957.t1
```

Do any features overlap?

```
bedtools intersect -a Ptuh_Results_sorted.gff3 -b /data/putnamlab/jillashey/genome/Pmea/Pocillopora_meandrina_HIv1.genes_sorted.gff3 > Ptuh_overlaps.bed

awk 'BEGIN {OFS="\t"} $3 != "Unknown_sRNA_locus"' Ptuh_overlaps.bed > filtered_Ptuh_overlaps.bed

wc -l filtered_Ptuh_overlaps.bed
174 filtered_Ptuh_overlaps.bed

head filtered_Ptuh_overlaps.bed
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	MIRNA_hairpin	818027	818120	12096	+	.	ID=Cluster_19;DicerCall=23;MIRNA=Y
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	mature_miRNA	818049	818070	3240	+	.	ID=Cluster_19.mature;Parent=Cluster_19
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	miRNA-star	818079	818100	9	+	.	ID=Cluster_19.star;Parent=Cluster_19
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	siRNA24_locus	5624826	5625249	125	+	.	ID=Cluster_73;DicerCall=24;MIRNA=N
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	siRNA24_locus	7205683	7206107	1021	+	.	ID=Cluster_100;DicerCall=24;MIRNA=N
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	siRNA22_locus	13015507	13015928	82	+	.	ID=Cluster_275;DicerCall=22;MIRNA=N
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	siRNA22_locus	16197870	16198683	17179	+	.	ID=Cluster_306;DicerCall=22;MIRNA=N
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	MIRNA_hairpin	20372416	20372510	522228	-	.	ID=Cluster_356;DicerCall=22;MIRNA=Y
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	miRNA-star	20372436	20372456	19166	-	.	ID=Cluster_356.star;Parent=Cluster_356
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	mature_miRNA	20372469	20372490	769	-	.	ID=Cluster_356.mature;Parent=Cluster_356
```

This doesn't tell me if the files are overlapping...it looks exactly like the results gff file. Does this indicate that most of the features have some kind of overlap??

If I subtract the 5th column (end of miRNA feature) by the 13th column (beginning of genomic feature), this will tell me if there is overlap in the two features. If the number is positive, it means there is overlap. If it is negative, it means there is no overlap. 

```
awk '{ $(NF+1) = $13 - $5 }1' OFS='\t' filtered_Ptuh_output.bed > filtered_Ptuh_output_overlap.bed

head filtered_Ptuh_output_overlap.bed
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	siRNA22_locus	173728	174150	1257	+	.	ID=Cluster_4;DicerCall=22;MIRNA=Pocillopora_meandrina_HIv1___Sc0000000	AUGUSTUS	CDS	174509	175333	.	-	0	Parent=Pocillopora_meandrina_HIv1___RNAseq.g20918.t1	359
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	siRNA22_locus	173728	174150	1257	+	.	ID=Cluster_4;DicerCall=22;MIRNA=Pocillopora_meandrina_HIv1___Sc0000000	AUGUSTUS	exon	174509	175333	.	-	0	Parent=Pocillopora_meandrina_HIv1___RNAseq.g20918.t1	359
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	siRNA22_locus	173728	174150	1257	+	.	ID=Cluster_4;DicerCall=22;MIRNA=Pocillopora_meandrina_HIv1___Sc0000000	AUGUSTUS	transcript	174509	176444	.	-	.	ID=Pocillopora_meandrina_HIv1___RNAseq.g20918.t1	359
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	MIRNA_hairpin	818027	818120	12096	+	.	ID=Cluster_19;DicerCall=23;MIRNA=Y	Pocillopora_meandrina_HIv1___Sc0000000	AUGUSTUS	transcript	816355	820160	.	+	.	ID=Pocillopora_meandrina_HIv1___RNAseq.g21001.t1	-1765
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	mature_miRNA	818049	818070	3240	+	.	ID=Cluster_19.mature;Parent=Cluster_19	Pocillopora_meandrina_HIv1___Sc0000000	AUGUSTUS	transcript	816355	820160	.	+	.	ID=Pocillopora_meandrina_HIv1___RNAseq.g21001.t1	-1715
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	miRNA-star	818079	818100	9	+	.	ID=Cluster_19.star;Parent=Cluster_19	Pocillopora_meandrina_HIv1___Sc0000000	AUGUSTUS	transcript	816355	820160	.	+	.	ID=Pocillopora_meandrina_HIv1___RNAseq.g21001.t1	-1745
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	MIRNA_hairpin	2872019	2872110	177	+	.	ID=Cluster_34;DicerCall=21;MIRNA=Y	Pocillopora_meandrina_HIv1___Sc0000000	AUGUSTUS	transcript	2868586	2871318	.	+	.	ID=Pocillopora_meandrina_HIv1___TS.g25957.t1	-3524
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	MIRNA_hairpin	2872019	2872110	177	+	.	ID=Cluster_34;DicerCall=21;MIRNA=Y	Pocillopora_meandrina_HIv1___Sc0000000	AUGUSTUS	CDS	2870522	2871318	.	+	2	Parent=Pocillopora_meandrina_HIv1___TS.g25957.t1	-1588
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	MIRNA_hairpin	2872019	2872110	177	+	.	ID=Cluster_34;DicerCall=21;MIRNA=Y	Pocillopora_meandrina_HIv1___Sc0000000	AUGUSTUS	exon	2870522	2871318	.	+	2	Parent=Pocillopora_meandrina_HIv1___TS.g25957.t1	-1588
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	mature_miRNA	2872041	2872061	110	+	.	ID=Cluster_34.mature;Parent=Cluster_34	Pocillopora_meandrina_HIv1___Sc0000000	AUGUSTUS	transcript	2868586	2871318	.	+	.	ID=Pocillopora_meandrina_HIv1___TS.g25957.t1	-3475
```

Count number of positive and negative numbers

```
awk '{if ($NF < 0) neg++; else if ($NF > 0) pos++} END {print "Positive count:", pos; print "Negative count:", neg}' filtered_Ptuh_output_overlap.bed

Positive count: 99
Negative count: 306
```

Negative = no overlap though I should check my math with someone to make sure I'm not dumb. 

For the miRNA data for all species, make a new miRNA directory and move everything there. 

```
cd /data/putnamlab/jillashey/e5/ncRNA_gff
mkdir miRNA
mv * miRNA/
```

Also renaming the output files with miRNA at the beginning of the file name. 




Questions to answer now that I have this data: 

- Which miRNA features are overlapping with their closest genomic feature? 
- On average, how far away are miRNAs from their closest genomic feature? 