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

##### Important piRNA files

- [Apul piRNA bed](https://github.com/urol-e5/deep-dive/blob/main/D-Apul/output/18-Apul-piRNA/0_piRNA_pipeline_proTRAC/proTRAC_results_APUL/APUL.merged.clusters.bed)
- [Peve piRNA bed](https://github.com/urol-e5/deep-dive/blob/main/E-Peve/output/18-Peve-piRNA/0_piRNA_pipeline_proTRAC/proTRAC_results_PEVE/PEVE.merged.clusters.bed)
- [Ptuh piRNA bed](https://github.com/urol-e5/deep-dive/blob/main/F-Pmea/output/18-Pmea-piRNA/0_piRNA_pipeline_proTRAC/proTRAC_results_PMEA/PMEA.merged.clusters.bed)

##### Important lncRNA files

- [Apul lncRNA bed](https://github.com/urol-e5/deep-dive/blob/main/D-Apul/output/05.33-lncRNA-discovery/Apul_lncRNA.bed); [Apul lncRNA fasta](https://github.com/urol-e5/deep-dive/blob/main/D-Apul/output/05.33-lncRNA-discovery/Apul_lncRNA.fasta)
- [Peve lncRNA bed](https://github.com/urol-e5/deep-dive/blob/main/E-Peve/output/Peve_lncRNA.bed); [Peve lncRNA fasta](https://github.com/urol-e5/deep-dive/blob/main/E-Peve/output/Peve_lncRNA.fasta)
- [Ptuh lncRNA bed](https://github.com/urol-e5/deep-dive/blob/main/F-Pmea/output/02-lncRNA-discovery/Pmea_lncRNA.bed); [Ptuh lncRNA fasta](https://github.com/urol-e5/deep-dive/blob/main/F-Pmea/output/02-lncRNA-discovery/Pmea_lncRNA.fasta)

Because I don't own the deep dive repo, I copied all needed files onto my local computer and then put them on Andromeda. I also added the species prefix to each file. 

```
cd /data/putnamlab/jillashey/e5
mkdir ncRNA_gff
cd ncRNA_gff
mkdir miRNA piRNA lncRNA
```

## Acropora pulchra 

### Apul miRNA

Sort the genome gff file for Apul. We are using the Amil genome for reference. 

```
cd /data/putnamlab/jillashey/genome/Amil_v2.01/
wget http://gannet.fish.washington.edu/seashell/snaps/GCF_013753865.1_Amil_v2.1_genomic.gff
sort -k1,1 -k4,4n GCF_013753865.1_Amil_v2.1_genomic.gff > GCF_013753865.1_Amil_v2.1_genomic_sorted.gff 
```

Sort miRNA gff file and run bed closest.

```
cd /data/putnamlab/jillashey/e5/ncRNA_gff/miRNA

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

Some of the closest features to the miRNAs are "regions", which seems to be large sections of chromosomes or chromosomes themselves. Because I only care about which mRNAs are closest to the miRNA features, I'm going to subset the gff by mRNA. 

```
cd /data/putnamlab/jillashey/genome/Amil_v2.01

awk '$3 == "mRNA"' GCF_013753865.1_Amil_v2.1_genomic_sorted.gff > GCF_013753865.1_Amil_v2.1_genomic_sorted_mRNA.gff
```

Run bed closest 

```
cd /data/putnamlab/jillashey/e5/ncRNA_gff/miRNA

interactive 
module load BEDTools/2.30.0-GCC-11.3.0
bedtools closest -a Apul_Results_sorted.gff3 -b /data/putnamlab/jillashey/genome/Amil_v2.01/GCF_013753865.1_Amil_v2.1_genomic_sorted_mRNA.gff > Apul_output_mRNA_only.bed

wc -l Apul_output.bed 
24039 Apul_output_mRNA_only.bed
```

Remove the unknown siRNA loci

```
awk 'BEGIN {OFS="\t"} $3 != "Unknown_sRNA_locus"' Apul_output_mRNA_only.bed > filtered_Apul_output_mRNA_only.bed
```

### Apul piRNA

The gff file is already sorted above and the piRNA bed file is also already sorted. Run bed closest 

```
cd /data/putnamlab/jillashey/e5/ncRNA_gff/piRNA

interactive 
module load BEDTools/2.30.0-GCC-11.3.0
bedtools closest -a APUL.merged.clusters.bed -b /data/putnamlab/jillashey/genome/Amil_v2.01/GCF_013753865.1_Amil_v2.1_genomic_sorted_mRNA.gff > Apul_piRNA_output_mRNA_only.bed

wc -l Apul_piRNA_output_mRNA_only.bed
143 Apul_piRNA_output_mRNA_only.bed

head Apul_piRNA_output_mRNA_only.bed
NC_058066.1	17726050	17734960	NC_058066.1	Gnomon	mRNA	17729627	17732829	.	+	.	ID=rna-XM_029347425.2;Parent=gene-LOC114967414;Dbxref=GeneID:114967414,Genbank:XM_029347425.2;Name=XM_029347425.2;gbkey=mRNA;gene=LOC114967414;model_evidence=Supporting evidence includes similarity to: 8 mRNAs%2C and 100%25 coverage of the annotated genomic feature by RNAseq alignments;product=zinc finger protein 862-like;transcript_id=XM_029347425.2
NC_058066.1	27441463	27447983	NC_058066.1	Gnomon	mRNA	27462364	27464454	.	+	.	ID=rna-XM_029327730.2;Parent=gene-LOC114951564;Dbxref=GeneID:114951564,Genbank:XM_029327730.2;Name=XM_029327730.2;gbkey=mRNA;gene=LOC114951564;model_evidence=Supporting evidence includes similarity to: 3 Proteins;product=zinc finger protein 862-like;transcript_id=XM_029327730.2
NC_058066.1	28121256	28125982	NC_058066.1	Gnomon	mRNA	28126954	28127481	.	+	.	ID=rna-XM_044328014.1;Parent=gene-LOC122964458;Dbxref=GeneID:122964458,Genbank:XM_044328014.1;Name=XM_044328014.1;gbkey=mRNA;gene=LOC122964458;model_evidence=Supporting evidence includes similarity to: 100%25 coverage of the annotated genomic feature by RNAseq alignments;product=uncharacterized protein K02A2.6-like;transcript_id=XM_044328014.1
NC_058066.1	28290198	28297001	NC_058066.1	Gnomon	mRNA	28275708	28276278	.	+	.	ID=rna-XM_029347200.2;Parent=gene-LOC114967202;Dbxref=GeneID:114967202,Genbank:XM_029347200.2;Name=XM_029347200.2;gbkey=mRNA;gene=LOC114967202;model_evidence=Supporting evidence includes similarity to: 1 mRNA%2C and 16%25 coverage of the annotated genomic feature by RNAseq alignments;product=piggyBac transposable element-derived protein 4-like;transcript_id=XM_029347200.2
NC_058066.1	28445323	28452636	NC_058066.1	Gnomon	mRNA	28444999	28446236	.	+	.	ID=rna-XM_029353270.2;Parent=gene-LOC114972849;Dbxref=GeneID:114972849,Genbank:XM_029353270.2;Name=XM_029353270.2;gbkey=mRNA;gene=LOC114972849;model_evidence=Supporting evidence includes similarity to: 67%25 coverage of the annotated genomic feature by RNAseq alignments;product=uncharacterized LOC114972849;transcript_id=XM_029353270.2
NC_058066.1	28445323	28452636	NC_058066.1	Gnomon	mRNA	28446503	28451426	.	-	.	ID=rna-XM_044328060.1;Parent=gene-LOC114972848;Dbxref=GeneID:114972848,Genbank:XM_044328060.1;Name=XM_044328060.1;gbkey=mRNA;gene=LOC114972848;model_evidence=Supporting evidence includes similarity to: 2 mRNAs%2C 1 Protein%2C and 100%25 coverage of the annotated genomic feature by RNAseq alignments%2C including 2 samples with support for all annotated introns;product=uncharacterized LOC114972848;transcript_id=XM_044328060.1
NC_058066.1	29297022	29310880	NC_058066.1	Gnomon	mRNA	29313159	29314285	.	+	.	ID=rna-XM_044317430.1;Parent=gene-LOC114972063;Dbxref=GeneID:114972063,Genbank:XM_044317430.1;Name=XM_044317430.1;experiment=COORDINATES: polyA evidence [ECO:0006239];gbkey=mRNA;gene=LOC114972063;model_evidence=Supporting evidence includes similarity to: 1 Protein%2C and 84%25 coverage of the annotated genomic feature by RNAseq alignments;product=uncharacterized LOC114972063;transcript_id=XM_044317430.1
NC_058067.1	14128407	14136847	NC_058067.1	Gnomon	mRNA	14118443	14130434	.	+	.	ID=rna-XM_044316573.1;Parent=gene-LOC122956872;Dbxref=GeneID:122956872,Genbank:XM_044316573.1;Name=XM_044316573.1;gbkey=mRNA;gene=LOC122956872;model_evidence=Supporting evidence includes similarity to: 3 mRNAs%2C 3 ESTs%2C 1 Protein%2C and 99%25 coverage of the annotated genomic feature by RNAseq alignments%2C including 35 samples with support for all annotated introns;product=uncharacterized LOC122956872;transcript_id=XM_044316573.1
NC_058067.1	14155111	14163931	NC_058067.1	Gnomon	mRNA	14157497	14159293	.	+	.	ID=rna-XM_044316312.1;Parent=gene-LOC122956629;Dbxref=GeneID:122956629,Genbank:XM_044316312.1;Name=XM_044316312.1;gbkey=mRNA;gene=LOC122956629;model_evidence=Supporting evidence includes similarity to: 4 Proteins;product=uncharacterized protein K02A2.6-like;transcript_id=XM_044316312.1
NC_058067.1	30302082	30308566	NC_058067.1	Gnomon	mRNA	30314470	30317762	.	+	.	ID=rna-XM_044317247.1;Parent=gene-LOC114948565;Dbxref=GeneID:114948565,Genbank:XM_044317247.1;Name=XM_044317247.1;gbkey=mRNA;gene=LOC114948565;model_evidence=Supporting evidence includes similarity to: 7 mRNAs%2C and 100%25 coverage of the annotated genomic feature by RNAseq alignments;product=uncharacterized LOC114948565;transcript_id=XM_044317247.1
```

### Apul lncRNA

The gff file is already sorted above. Sort lncRNa bed file. 

```
cd /data/putnamlab/jillashey/e5/ncRNA_gff/lncRNA
sort -k1,1 -k2,2n -k3,3n Apul_lncRNA.bed > Apul_lncRNA_sorted.bed
```

There were some issues with negative numbers being present in the starting coordinate position in some of the bed files. When I look at the fasta files for these specific lncRNAs, it says that the start position is 0. I'm going to change any instances of negative numbers to a 0

```
awk '{if ($2 < 0) $2 = 0; print $1 "\t" $2 "\t" $3}' Apul_lncRNA_sorted.bed > Apul_lncRNA_sorted_fixed.bed
```

Run bed closest

```
interactive 
module load BEDTools/2.30.0-GCC-11.3.0

bedtools closest -a Apul_lncRNA_sorted_fixed.bed -b /data/putnamlab/jillashey/genome/Amil_v2.01/GCF_013753865.1_Amil_v2.1_genomic_sorted_mRNA.gff > Apul_lncRNA_output_mRNA_only.bed

wc -l Apul_lncRNA_output_mRNA_only.bed
18475 Apul_lncRNA_output_mRNA_only.bed

head Apul_lncRNA_output_mRNA_only.bed
head Apul_lncRNA_output_mRNA_only.bed
NC_058066.1	393116	393357	NC_058066.1	Gnomon	mRNA	393884	399722	.	+	.	ID=rna-XM_029329017.2;Parent=gene-LOC114952935;Dbxref=GeneID:114952935,Genbank:XM_029329017.2;Name=XM_029329017.2;experiment=COORDINATES: polyA evidence [ECO:0006239];gbkey=mRNA;gene=LOC114952935;model_evidence=Supporting evidence includes similarity to: 2 mRNAs%2C 3 Proteins%2C and 100%25 coverage of the annotated genomic feature by RNAseq alignments%2C including 61 samples with support for all annotated introns;product=uncharacterized LOC114952935;transcript_id=XM_029329017.2
NC_058066.1	468617	469943	NC_058066.1	Gnomon	mRNA	470084	472938	.	+	.	ID=rna-XM_029329050.2;Parent=gene-LOC114952957;Dbxref=GeneID:114952957,Genbank:XM_029329050.2;Name=XM_029329050.2;experiment=COORDINATES: polyA evidence [ECO:0006239];gbkey=mRNA;gene=LOC114952957;model_evidence=Supporting evidence includes similarity to: 1 mRNA%2C 23 Proteins%2C and 100%25 coverage of the annotated genomic feature by RNAseq alignments;product=trace amine-associated receptor 1-like;transcript_id=XM_029329050.2
NC_058066.1	574074	574816	NC_058066.1	Gnomon	mRNA	566968	573222	.	-	.	ID=rna-XM_029328984.2;Parent=gene-LOC114952908;Dbxref=GeneID:114952908,Genbank:XM_029328984.2;Name=XM_029328984.2;experiment=COORDINATES: polyA evidence [ECO:0006239];gbkey=mRNA;gene=LOC114952908;model_evidence=Supporting evidence includes similarity to: 3 mRNAs%2C 2 ESTs%2C 8 Proteins%2C and 100%25 coverage of the annotated genomic feature by RNAseq alignments%2C including 88 samples with support for all annotated introns;product=titin-like%2C transcript variant X1;transcript_id=XM_029328984.2
NC_058066.1	574074	574816	NC_058066.1	Gnomon	mRNA	566968	573222	.	-	.	ID=rna-XM_044317851.1;Parent=gene-LOC114952908;Dbxref=GeneID:114952908,Genbank:XM_044317851.1;Name=XM_044317851.1;experiment=COORDINATES: polyA evidence [ECO:0006239];gbkey=mRNA;gene=LOC114952908;model_evidence=Supporting evidence includes similarity to: 1 EST%2C 7 Proteins%2C and 100%25 coverage of the annotated genomic feature by RNAseq alignments%2C including 51 samples with support for all annotated introns;product=titin-like%2C transcript variant X2;transcript_id=XM_044317851.1
NC_058066.1	852086	852315	NC_058066.1	Gnomon	mRNA	814602	850828	.	-	.	ID=rna-XM_029328850.2;Parent=gene-LOC114952824;Dbxref=GeneID:114952824,Genbank:XM_029328850.2;Name=XM_029328850.2;experiment=COORDINATES: polyA evidence [ECO:0006239];gbkey=mRNA;gene=LOC114952824;model_evidence=Supporting evidence includes similarity to: 16 mRNAs%2C 1 EST%2C 9 Proteins%2C and 100%25 coverage of the annotated genomic feature by RNAseq alignments%2C including 48 samples with support for all annotated introns;product=ubiquitin carboxyl-terminal hydrolase 24-like;transcript_id=XM_029328850.2
NC_058066.1	853114	853820	NC_058066.1	Gnomon	mRNA	854617	868846	.	+	.	ID=rna-XM_029328890.2;Parent=gene-LOC114952850;Dbxref=GeneID:114952850,Genbank:XM_029328890.2;Name=XM_029328890.2;Note=The sequence of the model RefSeq transcript was modified relative to this genomic sequence to represent the inferred CDS: deleted 1 base in 1 codon;exception=unclassified transcription discrepancy;experiment=COORDINATES: polyA evidence [ECO:0006239];gbkey=mRNA;gene=LOC114952850;model_evidence=Supporting evidence includes similarity to: 1 mRNA%2C 2 ESTs%2C 20 Proteins%2C and 99%25 coverage of the annotated genomic feature by RNAseq alignments%2C including 84 samples with support for all annotated introns;product=endoplasmic reticulum metallopeptidase 1-like;transcript_id=XM_029328890.2
NC_058066.1	946276	946580	NC_058066.1	Gnomon	mRNA	947261	949367	.	+	.	ID=rna-XM_029329051.2;Parent=gene-LOC114952958;Dbxref=GeneID:114952958,Genbank:XM_029329051.2;Name=XM_029329051.2;experiment=COORDINATES: polyA evidence [ECO:0006239];gbkey=mRNA;gene=LOC114952958;model_evidence=Supporting evidence includes similarity to: 8 mRNAs%2C 5 Proteins%2C and 100%25 coverage of the annotated genomic feature by RNAseq alignments%2C including 16 samples with support for all annotated introns;product=dynein light chain Tctex-type 5-B-like;transcript_id=XM_029329051.2
NC_058066.1	1132235	1134678	NC_058066.1	Gnomon	mRNA	1088762	1114844	.	+	.	ID=rna-XM_044318173.1;Parent=gene-LOC114952875;Dbxref=GeneID:114952875,Genbank:XM_044318173.1;Name=XM_044318173.1;experiment=COORDINATES: polyA evidence [ECO:0006239];gbkey=mRNA;gene=LOC114952875;model_evidence=Supporting evidence includes similarity to: 6 mRNAs%2C 2 ESTs%2C 10 Proteins%2C and 99%25 coverage of the annotated genomic feature by RNAseq alignments;product=transient receptor potential cation channel subfamily A member 1-like;transcript_id=XM_044318173.1
NC_058066.1	1135314	1144814	NC_058066.1	Gnomon	mRNA	1088762	1114844	.	+	.	ID=rna-XM_044318173.1;Parent=gene-LOC114952875;Dbxref=GeneID:114952875,Genbank:XM_044318173.1;Name=XM_044318173.1;experiment=COORDINATES: polyA evidence [ECO:0006239];gbkey=mRNA;gene=LOC114952875;model_evidence=Supporting evidence includes similarity to: 6 mRNAs%2C 2 ESTs%2C 10 Proteins%2C and 99%25 coverage of the annotated genomic feature by RNAseq alignments;product=transient receptor potential cation channel subfamily A member 1-like;transcript_id=XM_044318173.1
NC_058066.1	1144882	1148491	NC_058066.1	Gnomon	mRNA	1088762	1114844	.	+	.	ID=rna-XM_044318173.1;Parent=gene-LOC114952875;Dbxref=GeneID:114952875,Genbank:XM_044318173.1;Name=XM_044318173.1;experiment=COORDINATES: polyA evidence [ECO:0006239];gbkey=mRNA;gene=LOC114952875;model_evidence=Supporting evidence includes similarity to: 6 mRNAs%2C 2 ESTs%2C 10 Proteins%2C and 99%25 coverage of the annotated genomic feature by RNAseq alignments;product=transient receptor potential cation channel subfamily A member 1-like;transcript_id=XM_044318173.1
```

### Peve miRNA

Sort Peve genome gff file

```
cd /data/putnamlab/jillashey/genome/Peve/
sort -k1,1 -k4,4n Porites_evermanni_v1.annot.gff > Porites_evermanni_v1.annot_sorted.gff
``` 

Sort miRNA gff file and run `bed closest`. 

```
cd /data/putnamlab/jillashey/e5/ncRNA_gff/miRNA

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

### Peve piRNA

The gff file is already sorted above. Sort piRNA bed file and run bed closest 

```
cd /data/putnamlab/jillashey/e5/ncRNA_gff/piRNA

interactive 
module load BEDTools/2.30.0-GCC-11.3.0

bedtools sort -i PEVE.merged.clusters.bed > PEVE.merged.clusters_sorted.bed

bedtools closest -a PEVE.merged.clusters_sorted.bed -b /data/putnamlab/jillashey/genome/Peve/Porites_evermanni_v1.annot_sorted.gff > Peve_piRNA_output.bed

wc -l Peve_piRNA_output.bed
475 Peve_piRNA_output.bed

head Peve_piRNA_output.bed
Porites_evermani_scaffold_100	452587	464012	Porites_evermani_scaffold_100	Gmove	CDS	467740	467796	.	+	.	Parent=Peve_00000202
Porites_evermani_scaffold_100	452587	464012	Porites_evermani_scaffold_100	Gmove	mRNA	467740	468875	273	+	.	ID=Peve_00000202;Name=Peve_00000202;start=1;stop=1;cds_size=273
Porites_evermani_scaffold_1011	75258	81284	Porites_evermani_scaffold_1011	Gmove	mRNA	42046	82245	666.253	-	.	ID=Peve_00000418;Name=Peve_00000418;start=1;stop=1;cds_size=924
Porites_evermani_scaffold_1011	75258	81284	Porites_evermani_scaffold_1011	Gmove	CDS	78877	78941	.	+	.	Parent=Peve_00000424
Porites_evermani_scaffold_1011	75258	81284	Porites_evermani_scaffold_1011	Gmove	mRNA	78877	80087	1125	+	.	ID=Peve_00000424;Name=Peve_00000424;start=0;stop=1;cds_size=1125
Porites_evermani_scaffold_1011	75258	81284	Porites_evermani_scaffold_1011	Gmove	CDS	79012	79034	.	+	.	Parent=Peve_00000424
Porites_evermani_scaffold_1011	75258	81284	Porites_evermani_scaffold_1011	Gmove	CDS	79051	80087	.	+	.	Parent=Peve_00000424
Porites_evermani_scaffold_1024	29005	37949	Porites_evermani_scaffold_1024	Gmove	mRNA	28469	29122	195	-	.	ID=Peve_00000604;Name=Peve_00000604;start=0;stop=1;cds_size=195
Porites_evermani_scaffold_1024	29005	37949	Porites_evermani_scaffold_1024	Gmove	CDS	29089	29122	.	-	.	Parent=Peve_00000604
Porites_evermani_scaffold_1024	29005	37949	Porites_evermani_scaffold_1024	Gmove	CDS	30636	30770	.	-	.	Parent=Peve_00000611
```

### Ptuh miRNA

Sort Ptuh genome gff file; in this project, we used the Pmea genome. 

```
cd /data/putnamlab/jillashey/genome/Pmea
sort -k1,1 -k4,4n Pocillopora_meandrina_HIv1.genes.gff3 > Pocillopora_meandrina_HIv1.genes_sorted.gff3
```

Sort miRNA gff file and run `bed closest`. 

```
cd /data/putnamlab/jillashey/e5/ncRNA_gff/miRNA

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

Also renaming the output files with miRNA at the beginning of the file name. 

### Ptuh piRNA

The gff file is already sorted above and the piRNA bed file is also sorted. Run bed closest 

```
cd /data/putnamlab/jillashey/e5/ncRNA_gff/piRNA

interactive 
module load BEDTools/2.30.0-GCC-11.3.0

bedtools closest -a PMEA.merged.clusters.bed -b /data/putnamlab/jillashey/genome/Pmea/Pocillopora_meandrina_HIv1.genes_sorted.gff3 > Ptuh_piRNA_output.bed

wc -l Ptuh_piRNA_output.bed
647 Ptuh_piRNA_output.bed

head Ptuh_piRNA_output.bed
head Ptuh_piRNA_output.bed
Pocillopora_meandrina_HIv1___Sc0000000	10376955	10381586	Pocillopora_meandrina_HIv1___Sc0000000	AUGUSTUS	CDS	10377338	10377849	Parent=Pocillopora_meandrina_HIv1___RNAseq.g21904.t1
Pocillopora_meandrina_HIv1___Sc0000000	10376955	10381586	Pocillopora_meandrina_HIv1___Sc0000000	AUGUSTUS	exon	10377338	10377849	Parent=Pocillopora_meandrina_HIv1___RNAseq.g21904.t1
Pocillopora_meandrina_HIv1___Sc0000000	10376955	10381586	Pocillopora_meandrina_HIv1___Sc0000000	AUGUSTUS	transcript	10377338	10381414	.	-	.	ID=Pocillopora_meandrina_HIv1___RNAseq.g21904.t1
Pocillopora_meandrina_HIv1___Sc0000000	10376955	10381586	Pocillopora_meandrina_HIv1___Sc0000000	AUGUSTUS	CDS	10380139	10381414	Parent=Pocillopora_meandrina_HIv1___RNAseq.g21904.t1
Pocillopora_meandrina_HIv1___Sc0000000	10376955	10381586	Pocillopora_meandrina_HIv1___Sc0000000	AUGUSTUS	exon	10380139	10381414	Parent=Pocillopora_meandrina_HIv1___RNAseq.g21904.t1
Pocillopora_meandrina_HIv1___Sc0000001	7491180	7496937	Pocillopora_meandrina_HIv1___Sc0000001	AUGUSTUS	CDS	7491795	7493045	.	+	0	Parent=Pocillopora_meandrina_HIv1___RNAseq.g19232.t1
Pocillopora_meandrina_HIv1___Sc0000001	7491180	7496937	Pocillopora_meandrina_HIv1___Sc0000001	AUGUSTUS	exon	7491795	7493045	.	+	0	Parent=Pocillopora_meandrina_HIv1___RNAseq.g19232.t1
Pocillopora_meandrina_HIv1___Sc0000001	7491180	7496937	Pocillopora_meandrina_HIv1___Sc0000001	AUGUSTUS	transcript	7491795	7493045	.	+	.	ID=Pocillopora_meandrina_HIv1___RNAseq.g19232.t1
Pocillopora_meandrina_HIv1___Sc0000001	7491180	7496937	Pocillopora_meandrina_HIv1___Sc0000001	AUGUSTUS	CDS	7494081	7495580	.	-	0	Parent=Pocillopora_meandrina_HIv1___RNAseq.g19233.t1
Pocillopora_meandrina_HIv1___Sc0000001	7491180	7496937	Pocillopora_meandrina_HIv1___Sc0000001	AUGUSTUS	exon	7494081	7495580	.	-	0	Parent=Pocillopora_meandrina_HIv1___RNAseq.g19233.t1
```

### Ptuh lncRNA

The gff file is already sorted above. Sort the bed file 

