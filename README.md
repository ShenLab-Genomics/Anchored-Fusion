# Anchored-Fusion
The Anchored-fusion is a bioinformatic tool which detects fusion genes bases on paired-end RNA-sequence with ultrahigh sensitivity. By anchoring one fusion partner of interest, anchored-fusion can recover neglected reads in the de novo scenario, allowing for highly sensitive detection of gene fusion events from shallow sequencing data. Anchored-fusion use a deep learning model to predict if the candidate fusion gene is derived from natural chromsome fusion or artifical cDNA fusion during PCR.

Anchored-fusion needs the fasta file of target gene transcript reference. You can download the fasta file you need through NCBI GENE.

## Requirement

- samtools >= 1.16
- bwa >= 0.7.17
- BLAT >= v .35
- bedtools >= v2.30.0

### Environment

- Python 3.8 or higher
- numpy
- pytorch>=1.10.2
