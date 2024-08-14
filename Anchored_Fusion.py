import argparse
import sys
import os
import re
from Model import Train_model, Test_model
from functions import Gene_co, Find_homo_genes, Find_blocks,del_too_many_reads,contact_reads, Build_candidate_fasta, Find_Anchored_split,  Find_candidate_genes, Final_fusion,Find_fine_block, make_negative_file, get_test_reads
import numpy as np

lastd = sys.argv[0].rfind('/')
if lastd == -1:
    programdir = './'
else:
    programdir = sys.argv[0][:lastd] + '/'

parser = argparse.ArgumentParser(description='Anchor Gene Fusion Detection (c)')
parser.add_argument('--file_anchored_cds', type=str, required=True, default='', help='Target gene fasta file of anchored transcript')
parser.add_argument('--gene_names', type=str, default='', help='The file of target gene names')
parser.add_argument('--fastq1', type=str, default='fastq_1.fastq', help='The fastq1 file to scan')
parser.add_argument('--fastq2', type=str, default='fastq_2.fastq', help='The fastq2 file to scan')
parser.add_argument('--out_folder', type=str, default='./', help='The folder of the output file')
parser.add_argument('--file_ref_seq', type=str, required=True, default='', help='The reference sequence file')
parser.add_argument('--file_ref_ann', type=str, required=True, default='', help='The reference annotation file')
parser.add_argument('--not_filter_false_positive', action='store_true', help='When this parameter is added, possible artificial segments are not filtered out.')
parser.add_argument('--not_train_filter_model', action='store_true', help='When this parameter is added, the filtering model is not trained based on input data, and a pre-trained model file must be provided.')
parser.add_argument('--model_file', type=str, default='./data/model.pt', help='The model file for training model as negative samples')
parser.add_argument('--positive_samples', type=str, default='./data/positive_samples.txt', help='The fusion genes file for training model as positive samples')
parser.add_argument('--homo_gene_file', type=str, default='./data/homo_gene.npy', help='The homo genes file for training model. If you want to train a filtering model from the input data, you must include that file.')
parser.add_argument('--negative_samples', type=str, default='./Model/negative_samples.txt', help='The fusion genes file for training model as negative samples')
parser.add_argument('--thread', type=str, default='1', help='The threads number you want use.')
parser.add_argument('--gpu_number', type=str, default='-1', help='The gpu number you want use.')

args = parser.parse_args()

script_folder = os.path.dirname(os.path.abspath(__file__))
if not os.path.exists(args.positive_samples):
    args.positive_samples = script_folder+'/data/positive_samples.txt'
if not args.not_filter_false_positive:
    if not os.path.exists(args.homo_gene_file):
        args.homo_gene_file = script_folder+'/data/homo_gene.npy'
        if not os.path.exists(args.homo_gene_file):
            if os.path.exists(script_folder+'/data/homo_gene_1.txt.gz') and os.path.exists(script_folder+'/data/homo_gene_1.txt.gz'):
                os.system('gzip -d '+ script_folder+'/data/homo_gene_1.txt.gz')
                os.system('gzip -d '+ script_folder+'/data/homo_gene_2.txt.gz')
            if os.path.exists(script_folder+'/data/homo_gene_1.txt') and os.path.exists(script_folder+'/data/homo_gene_1.txt'):
                os.system('cat '+script_folder+'/data/homo_gene_1.txt '+script_folder+'/data/homo_gene_2.txt > '+script_folder+'/data/homo_gene.txt')
            if os.path.exists(script_folder+'/data/homo_gene.txt'):
                homo_genes = {}
                F = open(script_folder+'/data/homo_gene.txt','r')
                gene_names_now = F.readlines()[0].split(';')
                for line in F.readlines()[1:]:
                    key,values = line.rstrip().split('\t')
                    values = values.split(';')
                    values = [int(v) for v in values]
                    homo_genes[key] = values
                np.save(script_folder+'/data/homo_gene.npy',homo_genes)
                del homo_genes

gene_names = []
if not os.path.exists(args.gene_names) or args.gene_names == '':
    F_target_fasta = open(args.file_anchored_cds,'r')
    for line in F_target_fasta.readlines():
        if line.startswith('>'):
            arr = line.rstrip()[1:].split(' ')
            i = 0
            while i < len(arr):
                if re.match(r'[a-zA-Z]+_\d+\.\d+',arr[i]):
                    del arr[i]
                elif re.search(r'gene', arr[i], re.IGNORECASE) or re.search(r'specie', arr[i], re.IGNORECASE) or re.search(r'trans', arr[i], re.IGNORECASE) or re.search(r'for', arr[i], re.IGNORECASE) or re.search(r'homo', arr[i], re.IGNORECASE) or re.search(r'sapiens', arr[i], re.IGNORECASE):
                    del arr[i]
                else:
                    i += 1
            gene_names.append(arr[0])
    F_target_fasta.close()
else:
    F_target_genes = open(args.gene_names,'r')
    for line in F_target_genes.readlines():
        gene_name_now = line.rstrip()
        if gene_name_now != '':
            gene_names.append(gene_name_now)
    F_target_genes.close()

if not os.path.exists(args.out_folder):
   os.system('mkdir ' + args.out_folder)
model_out_name = args.out_folder+'/model_dir/'
if not os.path.exists(model_out_name):
   os.system('mkdir ' + model_out_name)

gene_co = Gene_co()
gene_co.Build_dic(args.file_ref_ann)
gpu_number =args.gpu_number
if not args.not_filter_false_positive and not args.not_train_filter_model:
    try:
        with open(args.positive_samples, 'r') as file:
            negative_samples = args.negative_samples
            if not os.path.exists(negative_samples):
                negative_samples = model_out_name + '/negative_samples.txt'
                file_unmapped_reads = model_out_name + '/unmapped_reads.bam'
            if not os.path.exists(negative_samples):
                try:
                    with open(args.homo_gene_file,'r') as file:
                        if not os.path.exists(file_unmapped_reads):
                            os.system("bwa mem -M -t " + args.thread + " " + args.file_ref_seq + " " + args.fastq1 + " " + args.fastq2 + " | samtools view -bSu - -o " + file_unmapped_reads)
                        make_negative_file(model_out_name, args.file_ref_seq, file_unmapped_reads ,negative_samples, gene_co, args.homo_gene_file, args.thread, gene_names)
                except FileNotFoundError:
                    print("Error: homo genes file not found!, not performing filter false positives.")
                    args.not_filter_false_positive = True
                    pass
            if not os.path.exists(args.model_file):
                model_file = model_out_name+'/model.pt'
            else:
                model_file = args.model_file
            Train_model(model_out_name, args.positive_samples,negative_samples,model_file,gpu_number)
    except FileNotFoundError:
        args.not_filter_false_positive = True
        print("Error: positive samples file not found!, not performing filter false positives.")
        pass
elif not args.not_filter_false_positive and args.not_train_filter_model:
    if not os.path.exists(args.model_file):
        model_file = model_out_name+'/model.pt'
    else:
        model_file = args.model_file

FAC = open(args.file_anchored_cds, 'r')
fasta_lines = FAC.readlines()
line_num = 1
for gene_name in gene_names:
    out_name = gene_name + '_fusion'
    temp_folder = args.out_folder + '/' + out_name
    temp_work_folder = temp_folder+ "/work_dir/"
    temp_model_folder = temp_folder + "/model_dir/"
    out_dir_name = temp_work_folder + out_name
    out_dir_name_out = temp_folder +'/' + out_name
    gpu_number =args.gpu_number

    if not os.path.exists(temp_folder):
        os.system('mkdir ' + temp_folder)
    if not os.path.exists(temp_work_folder):
        os.system('mkdir ' + temp_work_folder)
    if not os.path.exists(temp_model_folder):
        os.system('mkdir ' + temp_model_folder)

    file_anchored_seq = out_dir_name + "_anchored_gene_sequence.fa"
    file_candidate_seq = out_dir_name + "_candidate_gene_sequence.fa"
    file_realigned = out_dir_name + "_realign_reads.bam"
    file_spanning_reads = out_dir_name + "_spanning_reads.bam"
    file_anchored_reads = out_dir_name + "_anchored_reads.bam"
    file_anchored_reads_filter = out_dir_name + "_anchored_reads.sam"
    negative_samples = out_dir_name + '_negative_samples.txt'
    file_bad_genes = out_dir_name + '_homo_genes.bed'
    file_candidate_seq = out_dir_name + "_candidate_gene_sequence.fa"
    file_predictions_abridged = out_dir_name_out + "_predictions_abridged.txt"
    file_predictions = out_dir_name_out + "_predictions.txt"

    FAS = open(file_anchored_seq, 'w')
    FAS.write(">" + gene_name + "\n")
    while line_num < len(fasta_lines):
        if fasta_lines[line_num].startswith('>'):
            break
        else:
            FAS.write(fasta_lines[line_num])
        line_num += 1
    FAS.close()
    line_num += 1
    if line_num >= len(fasta_lines):
        FAC.close()

    if os.path.exists(file_anchored_seq + ".bwt") and os.path.exists(file_anchored_seq + ".pac") and os.path.exists(
            file_anchored_seq + ".amb") and os.path.exists(file_anchored_seq + ".ann") and os.path.exists(
            file_anchored_seq + ".sa"):
        pass
    else:
        os.system("bwa index " + file_anchored_seq)
    if os.path.exists(args.file_ref_seq + ".bwt") and os.path.exists(args.file_ref_seq + ".pac") and os.path.exists(
            args.file_ref_seq + ".amb") and os.path.exists(args.file_ref_seq + ".ann") and os.path.exists(
            args.file_ref_seq + ".sa"):
        pass
    else:
        os.system("bwa index " + args.file_ref_seq)


    if not os.path.exists(file_realigned):
        os.system("bwa mem -M -t " + args.thread + " " + file_anchored_seq + " " + args.fastq1 + " " + args.fastq2 + " | samtools view -bSu - | samtools sort - -o " + file_realigned)
    if not os.path.exists(file_spanning_reads):
        tmp_fastq1 = out_dir_name + '_tmp_1.fastq'
        tmp_fastq2 = out_dir_name + '_tmp_2.fastq'
        os.system("samtools view -u -f 8 -F 260 -@ " + args.thread + " " + file_realigned + " | samtools fastq - -o " + tmp_fastq1)
        os.system("samtools view -u -f 4 -F 264 -@ " + args.thread + " " + file_realigned + " | samtools fastq - -o " + tmp_fastq2)
        os.system("bwa mem -M -t " + args.thread + " " + args.file_ref_seq + " " + tmp_fastq1 + " " + tmp_fastq2 + " | samtools view -bSu - -o " + file_spanning_reads)
        if os.path.exists(tmp_fastq1):
            os.remove(tmp_fastq1)
        if os.path.exists(tmp_fastq2):
            os.remove(tmp_fastq2)
    if not os.path.exists(file_anchored_reads):
        os.system("samtools view -u -F 772 -h " + file_realigned + " |  samtools sort -@ " + args.thread + " - -o " + file_anchored_reads)

    if not os.path.exists(file_bad_genes):
        Find_homo_genes(args.file_ref_seq, file_anchored_seq, args.file_ref_ann, out_dir_name, file_bad_genes)
    homo_genes = []
    F = open(file_bad_genes,'r')
    for line in F.readlines():
        homo_genes.append(line.split('\t')[3])
    if not os.path.exists(file_anchored_reads_filter):
        del_too_many_reads(file_anchored_reads,file_anchored_reads_filter, out_dir_name, args.file_ref_seq, args.thread)
    blocks_chr = Find_blocks(file_spanning_reads, gene_co, homo_genes)
    blocks_chr = Find_fine_block(file_anchored_reads_filter,args.file_ref_seq,out_dir_name,gene_co, homo_genes,blocks_chr)
    Build_candidate_fasta(file_candidate_seq, out_dir_name, args.file_ref_seq, file_anchored_seq, blocks_chr)
    breakpoints = contact_reads(file_anchored_reads_filter, out_dir_name, args.file_ref_seq,args.thread)
    breakpoint_good_ids = Find_Anchored_split(out_dir_name, file_candidate_seq, blocks_chr, breakpoints,gene_co,file_anchored_seq)
    candidates_new, cnt_max = Find_candidate_genes(out_dir_name,  args.file_ref_seq, breakpoint_good_ids, breakpoints, blocks_chr,  gene_co, args.thread)
    model_out_name = temp_model_folder + out_name
    scores = []
    if not args.not_filter_false_positive and len(candidates_new) != 0:
      try:
        with open(model_file, 'r') as file:
            import torch
            state_dict = torch.load(model_file)
            test_file = model_out_name + '_test_reads.txt'
            get_test_reads(model_out_name,test_file,candidates_new, file_anchored_seq,args.file_ref_seq,gene_co )
            scores = Test_model(test_file, model_file,args.gpu_number)
            for i,candidate_new in enumerate(candidates_new):
                candidate_new.score = scores[i]
      except FileNotFoundError:
        args.not_filter_false_positive = True
        print("Error: model file not found!, not performing filter false positives.")
        pass

    Final_fusion(out_dir_name_out, candidates_new, gene_name, gene_co, scores,cnt_max,args.not_filter_false_positive)

exit()
