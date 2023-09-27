# Anchored-Fusion
The Anchored-fusion is a bioinformatic tool which detects fusion genes bases on paired-end RNA-sequence with ultrahigh sensitivity. By anchoring one fusion partner of interest, anchored-fusion can recover neglected reads in the de novo scenario, allowing for highly sensitive detection of gene fusion events from shallow sequencing data. Anchored-fusion use a deep learning model to predict if the candidate fusion gene is derived from natural chromsome fusion or artifical cDNA fusion during PCR.

Anchored-fusion target on the gene you choose, it will give the fusion genes including the target gene. Other fusion genes will not be represented in the output file.

If you are using the Conda package manager, execute `conda install git-lfs` to install the tool.

After installing Git-LFS, use the following command:
```bash
git lfs install
git lfs clone https://github.com/ShenLab-Genomics/Anchored-Fusion
```
to clone this repository.

## Requirement

- samtools >= 1.16
- bwa >= 0.7.17
- BLAT >= v .35
- bedtools >= v2.30.0

## Environment

- Python 3.8 or higher
- numpy
- pytorch>=1.10.2

## Genome reference and annotation

(1) Anchored-fusion needs the fasta file of target gene transcript reference. You can download the fasta file you need through NCBI GENE.
(2) Anchored-fusion needs the genome fa file, it can be downloaded from genecode or ensemble.
(3) Anchored-fusion also needs the gene annotation gtf file, hg38 and hg19 are both welcome. Please make sure that the versions of your reference file（fa） and annotation file (gtf) are compatible.

## Filter false fusion gene with deep learning model
（1）You can choose not to use the filter, which will save a lot of time. Please select the '--not_filter_false_positive' option.
（2）You can use our pre-trained model to filter out false fusion genes for your model, which can also save a lot of time. However, since it has not been trained on your input data, the results may have some deviations. Please select the '--not_train_filter_model' option. Please note that in order to do this, you need to provide the trained model. Please specify its location in '--model_file'. By default, we assume the model file is located in the 'data' folder of our program.
（3）For Anchored-fusion, by default, we expect you to retrain a filter using the data we provide, for which you need to provide positive training samples and the homologous gene file. These files are included in the 'data' folder of our program, please do not delete them.

## run
python Anchored_Fusion.py --file_anchored_cds=target_gene.fasta --gene_name=target_gene --fastq1=test_sample_1.fastq.gz  --fastq2=test_sample_1.fastq.gz --out_folder=output_dir --file_ref_seq=reference.genome.fa --file_ref_ann=genome.annotation.gtf
