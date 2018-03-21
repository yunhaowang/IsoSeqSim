# IsoSeqSim (version 0.2)
Time-stamp: <2017-10-09 Yunhao Wang, Email: yunhaowang@126.com>


## Introduction

Full-length isoform sequencing (Iso-Seq) technology originally developed by Pacific Biosciences (PacBio) has been widely applied to transcriptome study. In 2014, Oxford Nanopore Technologies (ONT) released its first sequencing platform, MinION, that has been also used to sequence full-length isoforms. IsoSeqSim is an Iso-Seq reads simulator for evaluating the performance of Iso-Seq bioinformatics analysis tools. IsoSeqSim has five modes: 1) "normal" for isoform construction and quantification analysis; 2) "fusion" for gene fusion analysis; 3) "apa" for alternative cleavage and polyadenylation (APA) analysis; 4) "ats" for alternative transcriptional start site (ATS) analysis; and 5) "ase" for allele-specific expression (ASE) analysis.


## Prerequisite

- Linux system

- python 2.7

- Numpy (tested with version 1.9.3)

- Scipy (tested with version 0.17.0)


## Install and Run

1. Download the package (e.g., `IsoSeqSim-0.2.tar.gz`) to a directory (e.g., `/home/`)

2. Unpack it using the command `tar -zxvf /home/IsoSeqSim-0.2.tar.gz`

3. Now, you can run IsoSeqSim by the executable file `/home/IsoSeqSim-0.2/bin/isoseqsim`. Optional, you can add IsoSeqSim into your PATH so that you can run IsoSeqSim without having to specify the entire path. For example, you can add one line `export PATH=/home/IsoSeqSim-0.2/bin:$PATH` to your `~/.bashrc`.


## Input

1. genome sequence file (FASTA format)

2. gene annotation file (GTF format). 

Note: considering gene duplication, wherein multiple genic loci (from same/different chromosome) use same isoform ID (e.g., RefSeq annotation library), make sure the isoform ID is unique for different genic locus. For human, suggest to use Ensembl or GENCODE but NOT RefSeq annotation library.

3. genotype file (VCF format, specifically for ASE analysis)


## Output

### 1. simulated read file (FASTA format)

See the naming rule of simulated reads for different modes:

- 1.1 "normal" mode

'>read_ID isoform_ID'. Space-split

- 1.2 "fusion" mode 

'>read_ID Isoform1+fusionSiteInIsoform1+Isoform2+fusionSiteInIsoform2'. Space-split

- 1.3 "apa" mode

'>read_ID isoform_ID+polyA_site'. Space-split

- 1.4 "ats" mode

'>read_ID isoform_ID+TSS'. Space-split

- 1.5 "ase" mode

'>read_ID isoform_ID allele_ID'. Space-split


### 2. transcript annotation file with read count (GenePred table format)

See format of modified GenePred table file for different modes:

- 2.1 "normal", "apa" and "ats" modes

Tab-split (12 columns)

(1) gene ID

(2) isoform ID

(3) chromosome ID

(4) strand: "+/plus" or "-/minus"

(5) transcription start site (TSS) for plus strand-derived isoform, while transcription terminal site (TTS) for minus strand-derived isoform

(6) TTS for plus strand-derived isoform, while TSS for minus strand-derived isoform

(7) .

(8) .

(9) number of exons

(10) exon start positions

(11) exon end positions

(12) simulated read count

- 2.2 "fusion" mode

Tab-split (23 columns)

(1-11) gpd for upstream isoform of fusion transcript

(12-22) gpd for downstream isoform of fusion transcript

(23) simulated read count

- 2.3 "ase" mode

Tab-split (14 columns)

(1-11) gpd for upstream isoform of fusion transcript

(12) simulated total read count

(13) simulated allele1-specific read count

(14) simulated allele2-specific read count


## Usage and Example

1. "normal" mode

`./bin/isoseqsim -g example/input/genome.fasta -a example/input/gene_annotation.gtf --c5 utilities/5_end_completeness.PacBio-P6-C4.tab --c3 utilities/3_end_completeness.PacBio-P6-C4.tab -o example/simulated_reads_normal.fa -t example/simulated_transcipt_normal.gpd --tempdir example/temp_normal`

2. "fusion" mode

`./bin/isoseqsim -g example/input/genome.fasta -a example/input/gene_annotation.gtf --c5 utilities/5_end_completeness.PacBio-P6-C4.tab --c3 utilities/3_end_completeness.PacBio-P6-C4.tab -o example/simulated_reads_fusion.fa -t example/simulated_transcipt_fusion.gpd --tempdir example/temp_fusion -m fusion --fc 20`

3. "apa" mode

`./bin/isoseqsim -g example/input/genome.fasta -a example/input/gene_annotation.gtf --c5 utilities/5_end_completeness.PacBio-P6-C4.tab --c3 utilities/3_end_completeness.PacBio-P6-C4.tab -o example/simulated_reads_apa.fa -t example/simulated_transcipt_apa.gpd --tempdir example/temp_apa -m apa --dis 50`

4. "ats" mode

`./bin/isoseqsim -g example/input/genome.fasta -a example/input/gene_annotation.gtf --c5 utilities/5_end_completeness.PacBio-P6-C4.tab --c3 utilities/3_end_completeness.PacBio-P6-C4.tab -o example/simulated_reads_ats.fa -t example/simulated_transcipt_ats.gpd --tempdir example/temp_ats -m ats --dis 50`

5. "ase" mode

`./bin/isoseqsim -g example/input/genome.fasta -a example/input/gene_annotation.gtf --c5 utilities/5_end_completeness.PacBio-P6-C4.tab --c3 utilities/3_end_completeness.PacBio-P6-C4.tab -o example/simulated_reads_ase.fa -t example/simulated_transcipt_ase.gpd --tempdir example/temp_ase -m ase --vcf example/input/genotype.vcf --id SIRV`


## Error Rate

- For PacBio Iso-Seq data, first extract ROI (reads of insert) with the full pass >= 0 and the accuracy >=0.7; then do alignment using GMAP aligner; last calculate the error rate and error pattern using AlignQC software

1. PacBio Sequel (Data: our own unpublished Iso-Seq data generated in 2017, 4 Sequel SMRT cells)

  (1) Substitution (mismatch) 1.731%; (2) Deletion 1.090%; (3) Insertion 2.204%

2. PacBio RS II with P6-C4 chemistry (Data: Alzheimer's disease brain Iso-Seq data released by PacBio in 2016, 40 SMRT cells)

  (1) Substitution (mismatch) 2.626%; (2) Deletion 2.156%; (3) Insertion 3.550%

3. PacBio RS II with P5-C3 chemistry (Data: MCF-7 Iso-Seq data released by PacBio in 2015, 28 SMRT cells)

  (1) Substitution (mismatch) 3.266%; (2) Deletion 1.618%; (3) Insertion 5.506%

4. PacBio RS II with P4-C2 chemistry (Data: MCF-7 Iso-Seq data released by PacBio in 2013, 119 SMRT cells)

  (1) Substitution (mismatch) 2.677%; (2) Deletion 1.787%; (3) Insertion 3.971%

- For ONT RNA-Seq data, first download data from a published study (Byrne, A. Nature Communications 2017); then do alignment using GMAP aligner; last calculate the error rate and error pattern using AlignQC software

5. ONT MinION with R7.3 chemistry (2D reads)

  (1) Substitution (mismatch) 5.391%; (2) Deletion 5.555%; (3) Insertion 2.324%

6. ONT MinION with R9.4 chemistry (2D reads)

  (1) Substitution (mismatch) 2.390%; (2) Deletion 5.355%; (3) Insertion 0.884%
