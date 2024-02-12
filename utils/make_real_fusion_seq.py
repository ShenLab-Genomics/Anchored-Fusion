import os
import random
import re

class Gene_co():
    def __init__(self):
        self.dic = {}

    def Build_dic(self, ref):
        F = open(ref, "r")
        for line in F.readlines():
            if re.match(r'^##', line):
                continue
            arr = line.split('\t')
            if arr[2] == "exon":
                tmp = re.findall(r'gene_id\s+"(ENSG\d+\S+)";\s+.+gene_type\s+"(\S+)";\s+.+gene_name\s+"(\S+)";\s+transcript_type\s+"protein_coding";\s+',arr[8])
                if len(tmp) != 0:
                    gene_id = tmp[0][0]
                    gene_type = tmp[0][1]
                    gene_name = tmp[0][2]
                    if arr[0] not in self.dic:
                        self.dic[arr[0]] = [[int(arr[3]), int(arr[4])+1, gene_id, gene_name]]
                    else:
                        self.dic[arr[0]].append([int(arr[3]), int(arr[4])+1, gene_id, gene_name])
        for key, value in self.dic.items():
            value.sort()
            i = 0
            while i < len(value) - 1:
                if value[i][1] >= value[i + 1][0]:
                    if value[i][1] >= value[i + 1][1]:
                        del value[i + 1]
                    else:
                        value[i][1] = value[i + 1][1]
                        del value[i + 1]
                    continue
                i += 1
        return

    def Find_gene(self, chrome, start, end):
        if chrome not in self.dic or chrome == 'chrM':
            return ['', '', '', '', '']
        chr_list = self.dic[chrome]
        l, r, m = 0, len(chr_list), 0
        while r - 1 > l:
            m = (l + r) // 2
            if chr_list[m][0] <= start:
                l = m
            else:
                r = m
        m = l
        if chr_list[m][0] <= start and chr_list[m][1] >= end:
            return [chr_list[m][2], chr_list[m][3], chrome, chr_list[m][0], chr_list[m][1]]
        else:
            return ['', '', '', '', '']

    def Find_exon(self, chrome, start, end):
        if chrome not in self.dic or chrome == 'chrM':
            return ['', '', '', '', ''], -1
        chr_list = self.dic[chrome]
        l, r, m = 0, len(chr_list), 0
        while r - 1 > l:
            m = (l + r) // 2
            if chr_list[m][0] <= start:
                l = m
            else:
                r = m
        m = l
        if chr_list[m][0] <= start and chr_list[m][1] >= end:
            return [chr_list[m][2], chr_list[m][3], chrome, chr_list[m][0], chr_list[m][1]],m
        else:
            return ['', '', '', '', ''],-1

def reverse(seq):
    r_dic = {'A': 'T', 'T': 'A', "G": 'C', "C": 'G'}
    seq = seq[::-1]
    seq2 = ''
    for c in seq:
        seq2 += r_dic[c]
    return seq2

def find_positions(gene_co,chrom,pos,length,gene_now,dir):
    gene,exon_num =  gene_co.Find_exon(chrom,pos,pos)
    poses = []
    forward_length = length
    backward_length = length
    if gene[0] == '':
        return poses
    exon_num_now = exon_num
    if dir == 'forward':
        pos_now = pos
        while forward_length > 0 :
           if gene_co.dic[chrom][exon_num_now][3] != gene[1] or gene_co.dic[chrom][exon_num_now][0]>pos_now or gene_co.dic[chrom][exon_num_now][1]<=pos_now:
              exon_num_now -= 1
              if exon_num_now < 0 or exon_num_now >= len(gene_co.dic[chrom]):
                  break
              else:
                  pos_now = gene_co.dic[chrom][exon_num_now][1] - 1
                  continue
           l = pos_now - gene_co.dic[chrom][exon_num_now][0] + 1
           if l >= forward_length:
              poses = [(pos_now - forward_length , pos_now )]+poses
              forward_length -= l
           elif l != 0 :
              forward_length -= l
              poses= [(gene_co.dic[chrom][exon_num_now][0] - 1,pos_now)]+poses
              exon_num_now -= 1
              pos_now = gene_co.dic[chrom][exon_num_now][1] - 1
           else:
              exon_num_now -= 1
              pos_now = gene_co.dic[chrom][exon_num_now][1] - 1

    if dir == 'backward':
        pos_now = pos
        while backward_length > 0 :
          if gene_co.dic[chrom][exon_num_now][3] != gene[1] or gene_co.dic[chrom][exon_num_now][0]>pos_now or gene_co.dic[chrom][exon_num_now][1]<=pos_now:
              exon_num_now += 1
              if  exon_num_now < 0 or exon_num_now >= len(gene_co.dic[chrom]):
                  break
              else:
                  pos_now = gene_co.dic[chrom][exon_num_now][0]
                  continue
          l = gene_co.dic[chrom][exon_num_now][1] - pos_now
          if l >= backward_length:
              poses.append((pos_now-1, pos_now + backward_length-1 ))
              backward_length = 0
          elif l != 0:
              backward_length -= l
              poses.append((pos_now-1 ,gene_co.dic[chrom][exon_num_now][1]-1))
              exon_num_now += 1
              pos_now = gene_co.dic[chrom][exon_num_now][0]
          else:
              exon_num_now += 1
              pos_now = gene_co.dic[chrom][exon_num_now][0]
    return poses

ann_file = './data/gencode.v19.chr_patch_hapl_scaff.annotation.gtf'
reference_fast_file = './data/hg19.fa'
gene_co = Gene_co()
gene_co.Build_dic(ann_file)
input_file = './data/positive/gene.fusions.V1.tsv'
output_file = './data/positive/positive_seq.txt'
tmp_bed = './data/positive/tmp.bed'
tmp_fasta = './data/positive/tmp.fasta'
O = open(output_file,'w')

T_B = open(tmp_bed,'w')
F = open(input_file,'r')
seq_last = ''
c1 = 0
c2 = 0
have_write = set()

for line in F.readlines()[1:]:
    arr = line.split('\t')
    if arr[0]+'$'+arr[7]+'$'+arr[9] not in have_write:
        have_write.add(arr[0]+'$'+arr[7]+'$'+arr[9])
    else:
        continue
    arr[-1] = arr[-1][:-1]
    fusiongene = arr[0]
    gene1 = arr[20].split('^')[1].split(':')
    gene2 = arr[21].split('^')[1].split(':')
    gene1_name,gene2_name = fusiongene.split('->')
    gene1_strand,gene2_strand = gene1[2],gene2[2]
    fusiongene = arr[0]
    if gene1_strand == '+':
        poses1 = find_positions(gene_co, 'chr'+gene1[0] ,int(gene1[1]),100, gene1_name,'forward')
    else:
        poses1 = find_positions(gene_co, 'chr'+gene1[0] ,int(gene1[1]),100, gene1_name,'backward')
    if gene2_strand == '+':
        poses2 = find_positions(gene_co, 'chr'+gene2[0] ,int(gene2[1]),100, gene2_name,'backward')
    else:
        poses2 = find_positions(gene_co, 'chr'+gene2[0] ,int(gene2[1]),100, gene2_name,'forward')
    if len(poses1) == 0 or len(poses2) == 0:
        continue
    flag = 'left'
    for pos in poses1:
        T_B.write('\t'.join(['chr'+gene1[0],str(pos[0]),str(pos[1]),fusiongene+'$'+gene1_name+'$'+gene1[2]+'$'+flag,gene1[2]])+'\n')
    flag = 'right'
    for pos in poses2:
        T_B.write('\t'.join(['chr'+gene2[0],str(pos[0]),str(pos[1]),fusiongene+'$'+gene2_name+'$'+gene2[2]+'$'+flag,gene2[2]])+'\n')
T_B.close()

os.system('bedtools getfasta -s -nameOnly -fi '+reference_fast_file+' -bed '+tmp_bed+' -fo '+tmp_fasta)

O = open(output_file,'w')
T_F = open(tmp_fasta,'r')
flag = 2
last_name = ''
seq_left,seq_right = '',''
last_strand = ''
write_name = ''
fusiongene = ''
last_fusiongene = ''
i = 0
lines = T_F.readlines()
while i < len(lines):
    if re.match(r'^>',lines[i]):
        fusiongene,gene_name,strand,dir = lines[i][1:-3].split('$')
        if gene_name != last_name or fusiongene!=last_fusiongene:
            if flag == 2:
                if last_name != '':
                    if last_strand == '-':
                        seq_right = reverse(seq_right)
                    O.write('N'*(100-len(seq_left))+seq_left+'H'+seq_right+'N'*(100-len(seq_right))+'\t'+write_name+'\n')
                flag = 1
                seq_left,seq_right = '',''
                write_name=gene_name
                last_name = gene_name
                last_strand = strand
                last_fusiongene = fusiongene
            else:
                if last_name != '':
                    if last_strand == '-':
                        seq_left = reverse(seq_left)
                flag = 2
                write_name = write_name + '->' + gene_name
                last_name = gene_name
                last_strand = strand
                last_fusiongene = fusiongene
    else:
          if dir == 'left':
              seq_left += lines[i][:-1].upper()
          if dir == 'right':
              seq_right += lines[i][:-1].upper()
    i += 1
O.close()
T_F.close()
