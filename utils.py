import sys
import os
import re
from Bio import Align
import numpy as np

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
                tmp = re.findall(r'gene_id\s+"(ENSG\d+\S+)";\s+.+gene_type\s+"(\S+)";\s+gene_name\s+"(\S+)";\s+',
                                 arr[8])
                if len(tmp) != 0:
                    gene_id = tmp[0][0]
                    gene_type = tmp[0][1]
                    gene_name = tmp[0][2]
                    if arr[0] not in self.dic:
                        self.dic[arr[0]] = [[int(arr[3]) - 100, int(arr[4]) + 100, gene_id, gene_name]]
                    else:
                        self.dic[arr[0]].append([int(arr[3]) - 100, int(arr[4]) + 100, gene_id, gene_name])
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
        if chrome == 'chr14' and 105586437 - 100 <= start and end <= 106879844 + 100:
            return ['IGH@', 'IGH@', 'chr14', 105586337, 106879944]
        elif chrome == 'chr14' and 21621904 - 100 <= start and end <= 22552232:
            return ['TRA@', 'TRA@', 'chr14', 21621804, 22552332]
        elif chrome == 'KI270846.1' and 0 <= start and end <= 1351393:
            return ['IGH@', 'IGH@', 'KI270846.1', 0, 1351393]
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


class Block():
    def __init__(self, chrom, start, end, gene):
        self.gene = gene
        self.bad = False
        self.chrom = chrom
        self.start = start
        self.end = end
        self.anchored_split_breakpoints = set()
        self.count = 0
        self.reads = []

    def Add_spanning_reads(self, start, end,read):
        if start < self.start:
            self.start = start
        if end > self.end:
            self.end = end
        self.count += 1
        self.reads.append(read)

    def Add_block(self, block):
        if self.start > block.start:
            self.start = block.start
        if self.end < block.end:
            self.end = block.end
        self.count += block.count
        self.reads.extend(block.reads)

class spanning_anchored():
    def __init__(self,type_,left_length,right_length,read_name):
        self.type_ = type_
        self.left_length = left_length
        self.right_length = right_length
        self.read_name = read_name

class Split_reads():
    def __init__(self, chrom, breakpoint, type_, seq_left, seq_right, read):
        self.chrom = chrom
        self.cnt = 1
        self.breakpoint = breakpoint
        self.type_ = type_
        self.seq_left = seq_left
        self.seq_right = seq_right
        self.other_breakpoints = []
        self.reads = [read]

    def add_reads(self, seq_left, seq_right,read):
        if len(seq_left) > len(self.seq_left):
            self.seq_left = seq_left
        if len(seq_right) > len(self.seq_right):
            self.seq_right = seq_right
        self.cnt += 1
        self.reads.append(read)

    def return_all_seq(self):
        return self.seq_left + self.seq_right

    def return_left_seq(self):
        return self.seq_left

    def return_right_seq(self):
        return self.seq_right

    def Add_other_breakpoint(self, chrom, breakpoint, strand, in_breakpoint, cut):
        self.other_breakpoints.append([chrom, breakpoint, strand, in_breakpoint, cut])


class Co_Split_reads():
    def __init__(self, chrom, breakpoint, type_):
        self.turn_dic = {'A': 0, 'T': 1, 'G': 2, 'C': 3}
        self.chrom = chrom
        self.breakpoint = breakpoint
        self.type_ = type_
        self.seq = [[0, 0, 0, 0] for _ in range(400)]
        self.in_breakpoint = 200
        self.l_left = 0
        self.l_right = 0
        self.cnt = 0
        self.reads = []

    def Add_reads(self, seq_left, seq_right, n, read,index):
        if index >= 0:
            seq_left += seq_right[:index]
            seq_right = seq_right[index:]
        else:
            seq_right = seq_left[index:]+seq_right
            seq_left = seq_left[:index]
        i = len(seq_left) - 1
        j = self.in_breakpoint - 1
        while i >= 0:
            if seq_left[i] in self.turn_dic:
                self.seq[j][self.turn_dic[seq_left[i]]] += n
            i -= 1
            j -= 1
        i = 0
        j = self.in_breakpoint
        while i < len(seq_right):
            if seq_right[i] in self.turn_dic:
                self.seq[j][self.turn_dic[seq_right[i]]] += n
            i += 1
            j += 1
        if self.l_left < len(seq_left):
            self.l_left = len(seq_left)
        if self.l_right < len(seq_right):
            self.l_right = len(seq_right)
        self.cnt += n
        self.reads.extend(read)
        return

    def find_max(self, i):
        turn_dic = {0: 'A', 1: 'T', 2: 'G', 3: 'C'}
        l = self.seq[i]
        max_value = max(l)
        max_indices = [j for j, x in enumerate(l) if x == max_value]
        if len(max_indices) > 1:
            return 'N'
        else:
            return turn_dic[max_indices[0]]

    def return_left_seq(self):
        seq_left = ''
        i = self.in_breakpoint - self.l_left
        while i < self.in_breakpoint:
            seq_left += self.find_max(i)
            i += 1
        return seq_left

    def return_right_seq(self):
        seq_right = ''
        i = self.in_breakpoint
        while i < self.in_breakpoint + self.l_right:
            seq_right += self.find_max(i)
            i += 1
        return seq_right



class Candidate_reads():
    def __init__(self, type_):
        self.pos = []
        self.type_ = type_
        self.spanning_reads = []
        self.split_reads = []
        self.l_left = 0
        self.l_right = 0
        self.l_mid = 0
        self.seq_left = [[0, 0, 0, 0] for _ in range(200)]
        self.seq_right = [[0, 0, 0, 0] for _ in range(200)]
        self.seq_mid = [[0, 0, 0, 0] for _ in range(100)]
        self.turn_dic = {'A': 0, 'T': 1, 'G': 2, 'C': 3}
        self.score = 0

    def add_reads(self, target_breakpoint, other_breakpoint, seq_left, seq_right,seq_mid, cnt,spanning_reads,split_reads):
        i = len(seq_left) - 1
        j = 199
        while i >= 0:
            if seq_left[i] in self.turn_dic:
                self.seq_left[j][self.turn_dic[seq_left[i]]] += cnt
            i -= 1
            j -= 1
        i = 0
        j = 0
        while i < len(seq_right):
            if seq_right[i] in self.turn_dic:
                self.seq_right[j][self.turn_dic[seq_right[i]]] += cnt
            i += 1
            j += 1
        i = 0
        j = 0
        while i < len(seq_mid):
            if seq_mid[i] in self.turn_dic:
                self.seq_mid[j][self.turn_dic[seq_mid[i]]] += cnt
            i += 1
            j += 1
        if self.l_left < len(seq_left):
            self.l_left = len(seq_left)
        if self.l_right < len(seq_right):
            self.l_right = len(seq_right)
        if self.l_mid < len(seq_mid):
            self.l_mid = len(seq_mid)

        flag = 0
        for pos in self.pos:
            if pos[0] == target_breakpoint and pos[1] == other_breakpoint[0] and pos[2] == other_breakpoint[1] and pos[
                3] == other_breakpoint[2] and pos[4] == other_breakpoint[3]:
                pos[5] += cnt
                flag = 1
                break
        if flag == 0:
            self.pos.append([target_breakpoint] + other_breakpoint + [cnt])
        self.spanning_reads.extend(spanning_reads)
        self.split_reads.extend(split_reads)
        return

    def find_max(self, i, seq):
        turn_dic = {0: 'A', 1: 'T', 2: 'G', 3: 'C'}
        l = seq[i]
        max_value = max(l)
        max_indices = [j for j, x in enumerate(l) if x == max_value]
        if len(max_indices) > 1:
            return 'N'
        else:
            return turn_dic[max_indices[0]]

    def return_seq_left(self):
        seq_left = ''
        i = 200 - self.l_left
        while i < 200:
            seq_left += self.find_max(i,self.seq_left)
            i += 1
        return seq_left

    def return_seq_right(self):
        seq_right = ''
        i = 0
        while i < self.l_right:
            seq_right += self.find_max(i,self.seq_right)
            i += 1
        return seq_right

    def return_seq_mid(self):
        seq_mid = ''
        i = 0
        while i < self.l_mid:
            seq_mid += self.find_max(i,self.seq_mid)
            i += 1
        return seq_mid

    def find_max_pos(self):
        i = 0
        max_cnt = 0
        max_id = 0
        while i < len(self.pos):
            if self.pos[i][6] > max_cnt:
                max_cnt = self.pos[i][6]
                max_id = i
            i += 1
        seq_left = self.return_seq_left()
        seq_right = self.return_seq_right()
        seq_mid = self.return_seq_mid()
        return self.pos[max_id] + [seq_left, seq_right, self.type_, seq_mid], max_id


def Find_homo_genes(f_ref, f_target, f_ann, out_dir_name, f_bad_genes):
    file_tmp_blat = out_dir_name + "_tmp_blat.psl"
    file_tmp_bed = out_dir_name + "_tmp_bed.bed"
    file_tmp_gtf = out_dir_name + "_tmp_gtf.bed"
    if not os.path.exists(file_tmp_blat):
        os.system("blat -stepSize=3 -repMatch=10000 -minScore=50 -minIdentity=80 " + f_ref + " " + f_target + " " + file_tmp_blat)
    F_blat = open(file_tmp_blat, "r")
    F_bed = open(file_tmp_bed, "w")
    for line in F_blat.readlines():
        if not re.match(r'^\d+', line):
            continue
        arr = line.split('\t')
        F_bed.write(arr[13] + "\t" + arr[15] + "\t" + arr[16] + "\t" + arr[9] + "\t" + arr[8] + "\n")
    F_bed.close()
    F_blat.close()
    F_gtf = open(f_ann, "r")
    F_t_gtf = open(file_tmp_gtf, "w")
    for line in F_gtf.readlines():
        if re.match(r'^##', line):
            continue
        arr = line.split('\t')
        if arr[2] == "gene":
            tmp = re.findall(r'gene_id\s+"(ENSG\d+\S+)";\s+.+gene_name\s+"(\S+)";\s+', arr[8])
            if len(tmp) != 0:
                gene_id = tmp[0][0]
                gene_name = tmp[0][1]
                F_t_gtf.write(
                    arr[0] + "\t" + arr[3] + "\t" + arr[4] + "\t" + gene_id + "\t" + gene_name + "\t" + arr[6] + "\n")
    F_gtf.close()
    F_t_gtf.close()
    os.system('bedtools intersect -a ' + file_tmp_gtf + ' -b ' + file_tmp_bed + ' -wa >' + f_bad_genes)
    if os.path.exists(file_tmp_blat):
        os.remove(file_tmp_blat)
    if os.path.exists(file_tmp_bed):
        os.remove(file_tmp_bed)
    if os.path.exists(file_tmp_gtf):
        os.remove(file_tmp_gtf)
    return


def Find_blocks(f_read, gene_co, homo_genes):
    def find_pos(line):
        arr = line.split('\t')
        letter_id = -1
        read_name, chrom, start, end = arr[0], '', 0, 0
        for index, c in enumerate(arr[5]):
            if c.isalpha():
                if c == 'M':
                    chrom, start, end = arr[2], int(arr[3]), int(arr[3]) + int(arr[5][letter_id + 1:index])
                    break
                letter_id = index
        return (read_name, chrom, start, end)
    blocks_chr = {}
    F = os.popen("samtools view " + f_read).read()
    lines = F.splitlines()
    j = 0
    while j < len(lines):
        k = 0
        poses = []
        chroms = set()
        nums = []
        while j + k < len(lines):
            pos = find_pos(lines[j + k])
            if len(poses) == 0 or pos[0] == poses[0][0]:
                poses.append(pos)
                chroms.add(pos[1])
                nums.extend([pos[2], pos[3]])
            else:
                break
            k += 1
        j += k
        if len(poses) == 1:
            continue
        if len(chroms) == 1 and max(nums) - min(nums) < 2000:
            continue
        gene_in_homo = -1
        gene_not_in_homo = -1
        gene_not_name = ''
        for k, pos in enumerate(poses):
            gene = gene_co.Find_gene(pos[1], pos[2], pos[3])
            if gene[0] in homo_genes:
                gene_in_homo = k
            elif gene[0] != '':
                if gene_not_in_homo == -1:
                    gene_not_in_homo = k
                    gene_not_name = gene
                elif gene[0] != gene_not_name[0]:
                    gene_not_in_homo = -1
                    break
        if gene_in_homo == -1 or gene_not_in_homo == -1:
            continue
        chrom, start, end, gene = poses[gene_not_in_homo][1], poses[gene_not_in_homo][2], poses[gene_not_in_homo][
            3], gene_not_name
        if chrom not in blocks_chr:
            blocks_chr[chrom] = []
        blocks = blocks_chr[chrom]
        i = len(blocks) - 1
        while i >= 0:
            if end < blocks[i].start:
                i -= 1
            else:
                break
        if i >= 0 and blocks[i].end + 300 >= start and blocks[i].gene[0] == gene[0]:
            blocks[i].Add_spanning_reads(start, end, poses[0][0])
        else:
            block = Block(chrom, start, end, gene)
            block.Add_spanning_reads(start, end, poses[0][0])
            if i != -1 and start < blocks[i].start:
                blocks.insert(i, block)
            else:
                blocks.insert(i + 1, block)
                i += 1
        while i < len(blocks) - 1 and blocks[i].gene[0] == blocks[i + 1].gene[0] and blocks[i].end + 300 >= blocks[
            i + 1].start:
            blocks[i].Add_block(blocks[i + 1])
            del blocks[i + 1]
        while i >= 1 and blocks[i].gene[0] == blocks[i - 1].gene[0] and blocks[i - 1].end >= blocks[i].start - 300:
            blocks[i].Add_block(blocks[i - 1])
            del blocks[i - 1]
            i -= 1
    for blocks in blocks_chr.values():
        for block in blocks:
            block.start, block.end = max(block.gene[3], block.start - 300), min(block.gene[4], block.end + 300)
    return blocks_chr

def reverse(seq):
    r_dic = {'A': 'T', 'T': 'A', "G": 'C', "C": 'G','N':'N','H':'H'}
    seq = seq[::-1]
    seq2 = ''
    for c in seq:
        seq2 += r_dic[c]
    return seq2

def Find_fine_block(f_read,file_ref_seq,work_dir,gene_co, homo_genes,block_chr):
    F = open(f_read,'r')
    fa_reads = work_dir + '_prb_spanning.fa'
    f_tmp_blat = work_dir + '_prb_spanning.psl'
    O = open(fa_reads,'w')
    i = 0
    candidates = []
    for line in F.readlines():
        arr = line.split('\t')
        read_name, chrom, start, cigar = arr[0], arr[2], arr[3], arr[5]
        cigar_list, seq = deal_cigar(cigar, arr[9])
        if len(cigar_list) != 2:
            continue
        left_length = cigar_list[0][1]
        right_length = cigar_list[1][1]
        if cigar_list[0][2] == 'S' and cigar_list[1][2] == 'M':
            type_ = 'SM'
        else:
            type_ = 'MS'
        O.write('>'+str(i)+'\n'+seq+'\n')
        candidates.append(spanning_anchored(type_,left_length,right_length,read_name))
        i += 1
    O.close()
    F.close()
    if not os.path.exists(f_tmp_blat):
        os.system("blat -minScore=20 " + file_ref_seq + " " + fa_reads + " " + f_tmp_blat)
    F = open(f_tmp_blat,'r')
    last_id = -1
    bad_flag = 0
    good_flag = 0
    blocks = []
    lines = F.readlines()
    i = 0
    while i <= len(lines):
        if i < len(lines):
            line = lines[i]
            if not re.match(r'^\d+\s', line):
                i+= 1
                continue
            arr = line.split('\t')
            now_id = int(arr[9])
        else:
            now_id = -2
        i += 1
        if now_id != last_id:
            last_id = now_id
            if bad_flag == 1 or good_flag == 0:
                bad_flag = 0
                good_flag = 0
                blocks = []
            else:
                for block in blocks:
                    gene = gene_co.Find_gene(block[0],block[1],block[2])
                    if gene[0] == '' or gene[0] in homo_genes:
                        continue
                    if block[0] not in block_chr:
                        Block_now = Block(block[0],block[1],block[2],gene)
                        Block_now.Add_spanning_reads(block[1], block[2],block[3])
                        block_chr[block[0]] = [Block_now]
                    else:
                        Block_now = block_chr[block[0]]
                        j = len(Block_now) - 1
                        while j >= 0:
                            if block[2] < Block_now[j].start:
                                j -= 1
                            else:
                                break
                        if j >= 0 and Block_now[j].end  >= block[1] and Block_now[j].gene[0] == gene[0]:
                            Block_now[j].Add_spanning_reads(block[1], block[2], block[3])
                        else:
                            block_now = Block(block[0],block[1],block[2],gene)
                            block_now.Add_spanning_reads(block[1], block[2], block[3])
                            if j != -1 and block[1] < Block_now[j].start:
                                Block_now.insert(j, block_now)
                            else:
                                Block_now.insert(j + 1, block_now)
                                j += 1
                        while j < len(Block_now) - 1 and Block_now[j].gene[0] == Block_now[j + 1].gene[0] and Block_now[j].end  >= Block_now[j + 1].start:
                            Block_now[j].Add_block(Block_now[j + 1])
                            del Block_now[j + 1]
                        while j >= 1 and Block_now[j].gene[0] == Block_now[j - 1].gene[0] and Block_now[j - 1].end >= Block_now[j].start :
                            Block_now[j].Add_block(Block_now[j - 1])
                            del Block_now[j - 1]
                            j -= 1
                blocks = []
        if bad_flag != 0:
            continue
        if i > len(lines):
            break
        chrom, start, end, in_start,in_end = arr[13],int(arr[15]),int(arr[16]),int(arr[11]),int(arr[12])
        if end - start > 200:
            continue
        if candidates[now_id].type_ == "MS":
            if in_start <= candidates[now_id].left_length//2 and in_end >= candidates[now_id].left_length + 5:
                bad_flag = 1
            elif candidates[now_id].left_length - 5 <= in_start <= candidates[now_id].left_length + 5 and in_end >= candidates[now_id].left_length + candidates[now_id].right_length  - 5:
                blocks.append([chrom,start,end,candidates[now_id].read_name])
            elif in_start <= 5 and in_end <= candidates[now_id].left_length + 5:
                gene = gene_co.Find_gene(chrom,start,end)
                if gene[0] in homo_genes:
                    good_flag = 1
        elif candidates[now_id].type_ == "SM":
            if candidates[now_id].left_length - 5<=in_end <= candidates[now_id].left_length + 5 and in_start <= 5:
                blocks.append([chrom,start,end,candidates[now_id].read_name])
            elif in_start < candidates[now_id].left_length - 5 and in_end >= candidates[now_id].left_length + candidates[now_id].right_length//2:
                bad_flag = 1
            elif candidates[now_id].left_length - 5 <= in_start <= candidates[now_id].left_length + 5 and in_end >= candidates[now_id].left_length + candidates[now_id].right_length  - 5:
                gene = gene_co.Find_gene(chrom,start,end)
                if gene[0] in homo_genes:
                    good_flag = 1
    return block_chr

def deal_cigar(cigar, seq):
    co = 0
    b_id = -1
    res = []
    for id_, c in enumerate(cigar):
        if c.isalpha():
            num = int(cigar[b_id + 1:id_])
            co += num
            res.append([co, num, c])
            b_id = id_
    i = 0
    while i < len(res):
        if res[i][2] == 'N':
            j = i + 1
            while j < len(res):
                res[j][0] -= res[i][1]
                j += 1
            del res[i]
        elif res[i][2] == 'D':
            if i != len(res) - 1:
                res[i + 1][1] += res[i][1]
            seq = seq[:res[i - 1][0]] + 'N' * res[i][1] + seq[res[i - 1][0]:]
            del res[i]
        elif res[i][2] == 'I':
            j = i + 1
            while j < len(res):
                res[j][0] -= res[i][1]
                j += 1
            seq = seq[:res[i - 1][0]] + seq[res[i][0]:]
            del res[i]
        elif res[i][2] == 'H':
            j = i + 1
            while j < len(res):
                res[j][0] -= res[i][1]
                j += 1
            del res[i]
        else:
            i += 1
    i = 0
    while i < len(res) - 1:
        if res[i][2] == 'M' and res[i + 1][2] == 'M':
            res[i][0] = res[i + 1][0]
            res[i][1] += res[i + 1][1]
            del res[i + 1]
        else:
            i += 1
    return res, seq


def del_too_many_reads(f_read,out_sam, work_dir, f_ref, thread):
    tmp_fa = work_dir + '_del_tmp.fa'
    tmp_sam = work_dir + '_del_tmp.sam'
    F = os.popen("samtools view " + f_read).read()
    FTF = open(tmp_fa, 'w')
    for line in F.splitlines():
        arr = line.split('\t')
        cigar, seq = deal_cigar(arr[5], arr[9])
        if len(cigar) == 2:
            FTF.write('>' + '$'.join([arr[0],arr[2],arr[3],arr[5]]) + '\n' + arr[9] + '\n')
    FTF.close()
    os.system("bwa mem -M -t " + thread + " " + f_ref + " " + tmp_fa + " > " + tmp_sam)
    FS = open(tmp_sam, 'r')
    FO = open(out_sam,'w')
    i = 0
    lines = FS.readlines()
    last_name = ''
    last_seq = ''
    bad_flag = 0
    while i < len(lines):
        if re.match(r'^@', lines[i]):
            i += 1
            continue
        arr = lines[i].split('\t')
        if int(arr[1]) > 15 and bin(int(arr[1]))[-5] == '1':
            arr[9] = reverse(arr[9])
        now_name,now_seq = arr[0],arr[9]
        if now_name.find('SRR8616157.67784184')!=-1:
            print('lalala')
        if now_name != last_name:
            if bad_flag == 0 and last_name != '':
                arr_name = last_name.split('$')
                FO.write(arr_name[0]+'\t0\t'+arr_name[1]+'\t'+arr_name[2]+'\t60\t'+arr_name[3]+'\t=\t1111\t0\t'+last_seq+'\tA\n')
            last_name,bad_flag,last_seq = now_name,0,now_seq
        if bad_flag == 0:
            id_H = arr[5].find('H')
            if id_H > -1:
                arr[5] = arr[5][:id_H]+'S'+arr[5][id_H+1:]
            cigar_now,seq = deal_cigar(arr[5],arr[9])
            cigar_before,seq = deal_cigar(arr[0].split('$')[3],arr[9])
            if int(arr[1]) > 15 and bin(int(arr[1]))[-5] == '1':
                cigar_now = cigar_now[::-1]
                t = 0
                for m in cigar_now:
                    t += m[1]
                    m[0] = t
            if len(cigar_now) == 1:
                bad_flag = 1
            elif len(cigar_now) >= 2:
                for mb in cigar_before:
                    if mb[2] == 'M':
                        for m in cigar_now:
                            if m[2] == 'M' and m[0]-m[1] < mb[0]- mb[1]*0.2 and m[0] > mb[0]+mb[1] * 0.2:
                                bad_flag = 1
        i += 1
    FS.close()
    FO.close()
    #if os.path.exists(tmp_fa):
        #os.remove(tmp_fa)
    #if os.path.exists(tmp_sam):
        #os.remove(tmp_sam)
    return


def combine_split_reads(breakpoints):
    aligner = Align.PairwiseAligner()
    aligner.open_gap_score = -1
    aligner.extend_gap_score = -0.1
    aligner.match_score = 1
    aligner.mismatch_score = -0.5

    def if_similar(seq_left1, seq_right1, seq_left2, seq_right2, threshold, index):
        seq_left1, seq_left2 = seq_left1[::-1], seq_left2[::-1]
        if len(seq_left1)== 0 or len(seq_left2)== 0 or len(seq_right1)== 0 or len(seq_right2)== 0:
            return False
        if index >= 0:
            seq_left2 = seq_left2[index:]
        else:
            seq_left1 = seq_left1[-1*index:]
        c = 0
        if len(seq_left1)== 0 or len(seq_left2)== 0:
            return False
        for i in range(min(len(seq_left1),len(seq_left2))):
            if seq_left1[i] == seq_left2[i]:
                c += 1
        score_left = c/min(len(seq_left1),len(seq_left2))
        if index >= 0:
            seq_right1 = seq_right1[index:]
        else:
            seq_right2 = seq_right2[-1*index:]
        c = 0
        if len(seq_right1)== 0 or len(seq_right2)== 0:
            return False
        for i in range(min(len(seq_right1),len(seq_right2))):
            if seq_right1[i] == seq_right2[i]:
                c += 1
        score_right = c/min(len(seq_right1),len(seq_right2))
        if score_left > threshold and score_right > threshold:
            return True
        return False

    breakpoints_new = []
    for breakpoint_list in breakpoints:
        breakpoints_new.append([breakpoint_list[0], []])
        breakpoint_list_new = breakpoints_new[-1]
        co_Split_reads_now = Co_Split_reads(breakpoint_list[1][0].chrom, breakpoint_list[0],
                                            breakpoint_list[1][0].type_)
        while len(breakpoint_list[1]) > 0:
            type_1 = breakpoint_list[1][0].type_
            seq_left1 = breakpoint_list[1][0].seq_left
            seq_right1 = breakpoint_list[1][0].seq_right
            cnt1 = breakpoint_list[1][0].cnt
            reads1 = breakpoint_list[1][0].reads
            co_Split_reads_now.Add_reads(seq_left1, seq_right1, cnt1,reads1,0)
            j = 1
            while j < len(breakpoint_list[1]):
                type_2 = breakpoint_list[1][j].type_
                seq_left2 = breakpoint_list[1][j].seq_left
                seq_right2 = breakpoint_list[1][j].seq_right
                cnt2 = breakpoint_list[1][j].cnt
                reads2 = breakpoint_list[1][j].reads
                if type_1 == type_2 and if_similar(seq_left1, seq_right1, seq_left2, seq_right2, 0.9, 0):
                    co_Split_reads_now.Add_reads(seq_left2, seq_right2, cnt2, reads2,0)
                    del breakpoint_list[1][j]
                else:
                    j += 1
            breakpoint_list_new[1].append(co_Split_reads_now)
            del breakpoint_list[1][0]
            if len(breakpoint_list[1]) == 0:
                break
            co_Split_reads_now = Co_Split_reads(breakpoint_list[1][0].chrom, breakpoint_list[0],
                                                breakpoint_list[1][0].type_)
    breakpoints = breakpoints_new
    i = 0
    while i < len(breakpoints):
        j = 0
        while j < len(breakpoints[i][1]):
            seq_left1 = breakpoints[i][1][j].return_left_seq()
            seq_right1 = breakpoints[i][1][j].return_right_seq()
            cnt1 = breakpoints[i][1][j].cnt
            reads1 = breakpoints[i][1][j].reads
            type_1 = breakpoints[i][1][j].type_
            z = i + 1
            flag = 0
            while z < len(breakpoints) and breakpoints[z][0] - breakpoints[i][0] <= 3:
                k = 0
                while k < len(breakpoints[z][1]):
                    seq_left2 = breakpoints[z][1][k].return_left_seq()
                    seq_right2 = breakpoints[z][1][k].return_right_seq()
                    cnt2 = breakpoints[z][1][k].cnt
                    reads2 = breakpoints[z][1][k].reads
                    type_2 = breakpoints[z][1][k].type_
                    if type_1 == type_2 and if_similar(seq_left1, seq_right1, seq_left2, seq_right2, 0.9, breakpoints[z][0] - breakpoints[i][0]):
                        if cnt1 > cnt2:
                            seq_right2 = seq_left2[-1 * (breakpoints[z][0] - breakpoints[i][0]):] + seq_right2
                            seq_left2 = seq_left2[:-1 * (breakpoints[z][0] - breakpoints[i][0])]
                            del breakpoints[z][1][k]
                            breakpoints[i][1][j].Add_reads(seq_left2, seq_right2, cnt2, reads2,breakpoints[i][0] - breakpoints[z][0])
                        else:
                            seq_left1 += seq_right1[:breakpoints[z][0] - breakpoints[i][0]]
                            seq_right1 = seq_right1[breakpoints[z][0] - breakpoints[i][0]:]
                            flag = 1
                            del breakpoints[i][1][j]
                            breakpoints[z][1][k].Add_reads(seq_left1, seq_right1, cnt1, reads1,breakpoints[z][0] - breakpoints[i][0])
                            break
                    else:
                        k += 1
                if flag == 1:
                    break
                z += 1
            if flag == 0:
                j += 1
        i += 1

    breakpoints_new = []
    for breakpoint_list in breakpoints:
        for breakpoint_now in breakpoint_list[1]:
            split_reads_now = Split_reads(breakpoint_now.chrom, breakpoint_now.breakpoint, breakpoint_now.type_,
                                          breakpoint_now.return_left_seq(), breakpoint_now.return_right_seq(),'')
            split_reads_now.cnt = breakpoint_now.cnt
            split_reads_now.reads = breakpoint_now.reads
            breakpoints_new.append(split_reads_now)
    return breakpoints_new


def contact_reads(f_read, work_dir, f_ref, thread):
    def biselect(breakpoints, breakpointnow):
        left, right = 0, len(breakpoints)
        if right == 1 and breakpoints[0][0] == breakpointnow:
            return 1, 0
        while left < right - 1:
            mid = (left + right) // 2
            if breakpoints[mid][0] == breakpointnow:
                return 1, mid
            if breakpoints[mid][0] < breakpointnow:
                left = mid
            else:
                right = mid
        return 0, right

    def if_same(seq1_left, seq1_right, seq2_left, seq2_right):
        if seq1_left[max(len(seq1_left) - len(seq2_left), 0):] == seq2_left[max(len(seq2_left) - len(seq1_left), 0):]:
            if seq1_right[:min(len(seq1_right), len(seq2_right))] == seq2_right[:min(len(seq1_right), len(seq2_right))]:
                return True
        return False

    breakpoints = []
    F = open(f_read,'r')
    for line in F.readlines():
        arr = line.split('\t')
        read_name, chrom, start, cigar = arr[0], arr[2], arr[3], arr[5]
        cigar_list, seq = deal_cigar(cigar, arr[9])
        if len(cigar_list) != 2:
            continue
        if cigar_list[0][2] == 'S' and cigar_list[1][2] == 'M':
            type_ = 'SM'
            if cigar_list[0][0] < 15:
                continue
            breakpointnow = int(start)
        else:
            type_ = 'MS'
            if cigar_list[1][0] < 15:
                continue
            breakpointnow = int(start) + cigar_list[0][1]
        flag, index = biselect(breakpoints, breakpointnow)
        if flag == 0:
            breakpoints.insert(index, (
            breakpointnow, [Split_reads(chrom, breakpointnow, type_, seq[:cigar_list[0][0]], seq[cigar_list[0][0]:], read_name)]))
        else:
            breakpointnow_list = breakpoints[index][1]
            seq2_left, seq2_right = seq[:cigar_list[0][0]], seq[cigar_list[0][0]:]
            i = len(breakpointnow_list) - 1
            while i >= 0:
                breakpoint = breakpointnow_list[i]
                if type_ == breakpoint.type_:
                    seq1_left, seq1_right = breakpoint.return_left_seq(), breakpoint.return_right_seq()
                    if if_same(seq1_left, seq1_right, seq2_left, seq2_right):
                        breakpoint.add_reads(seq2_left, seq2_right, read_name)
                        break
                i -= 1
            if i < 0:
                breakpointnow_list.append(
                    Split_reads(chrom, breakpointnow, type_, seq[:cigar_list[0][0]], seq[cigar_list[0][0]:], read_name))
    F.close()
    breakpoints = combine_split_reads(breakpoints)
    return breakpoints


def Build_candidate_fasta(file_candidate_seq, work_dir, f_ref, f_target, blocks_chr):
    tmp_bed = work_dir + "_tmp.bed"
    f_tmp_blat = work_dir + '_tmp_blat.psl'
    TB = open(tmp_bed, 'w')
    for blocks in blocks_chr.values():
        for i, block in enumerate(blocks):
            TB.write(block.chrom + '\t' + str(block.start) + '\t' + str(block.end) + '\t' + str(i) + '\n')
    TB.close()
    os.system("bedtools getfasta -name -fi " + f_ref + " -bed " + tmp_bed + " > " + file_candidate_seq)
    os.system("blat -stepSize=3 -minScore=20 -minMatch=2 -minIdentity=0 " + file_candidate_seq + " " + f_target + " " + f_tmp_blat)
    FB = open(f_tmp_blat, 'r')
    TB = open(tmp_bed, 'w')
    for line in FB.readlines():
        if not re.match(r'^\d+\s', line):
            continue
        arr = line.split('\t')
        block_id, chrom_y = int(arr[13].split('::')[0]), arr[13].split('::')[1].split(':')[0]
        blocks_chr[chrom_y][block_id].bad = True
    for blocks in blocks_chr.values():
        i = 0
        while i < len(blocks):
            if blocks[i].bad == True:
                del blocks[i]
            else:
                TB.write(
                    blocks[i].chrom + '\t' + str(blocks[i].start) + '\t' + str(blocks[i].end) + '\t' + str(i) + '\n')
                i += 1
    TB.close()
    os.system("bedtools getfasta -name -fi " + f_ref + " -bed " + tmp_bed + " > " + file_candidate_seq)
    os.system("bwa index " + file_candidate_seq)
    if os.path.exists(tmp_bed):
        os.remove(tmp_bed)
    if os.path.exists(f_tmp_blat):
        os.remove(f_tmp_blat)
    return


def Find_Anchored_split(work_dir, file_candidate_seq, blocks_chr, breakpoints, gene_co, file_anchored_seq):
    f_tmp_fa = work_dir + "_tmp_anchored_fa.fa"
    f_tmp_blat = work_dir + "_tmp_anchored_blat.psl"
    FA = open(f_tmp_fa, 'w')
    for i, breakpoint in enumerate(breakpoints):
        if breakpoints[i].type_ == 'SM':
            read = breakpoints[i].return_left_seq()
        else:
            read = breakpoints[i].return_right_seq()
        FA.write('>' + str(i) + '\n' + read + '\n')
    FA.close()
    os.system("blat -stepSize=3 -minScore=12 -minMatch=2 -minIdentity=90 " + file_candidate_seq + " " + f_tmp_fa + " " + f_tmp_blat)
    breakpoint_id_good = set()
    FB = open(f_tmp_blat, 'r')
    for line in FB.readlines():
        if not re.match(r'^\d+\s', line):
            continue
        arr = line.split('\t')
        l, s, e = int(arr[10]), int(arr[11]), int(arr[12])
        if s > 5 and e < l - 5:
            continue
        breakpoint_id, block_id, chrom_y, start = int(arr[9]), int(arr[13].split('::')[0]), \
                                                  arr[13].split('::')[1].split(':')[0], int(
            arr[13].split('::')[1].split(':')[1].split('-')[0])
        start_y, end_y, strand = int(arr[15]) + start, int(arr[16]) + start, arr[8]
        gene = gene_co.Find_gene(chrom_y, start_y, end_y)
        #if (s > 5 or e < l - 5) and (gene[0] != 'IGH@' and gene[0]!='TRA@'):
            #continue
        breakpoint = breakpoints[breakpoint_id]
        if breakpoint.type_ == 'SM':
            if strand == '+':
                breakpoint.Add_other_breakpoint(chrom_y, end_y, strand, s, l-e)
            else:
                breakpoint.Add_other_breakpoint(chrom_y, start_y, strand, s, l-e)
        else:
            if strand == '+':
                breakpoint.Add_other_breakpoint(chrom_y, start_y, strand, s, l - e)
            else:
                breakpoint.Add_other_breakpoint(chrom_y, end_y, strand, s, l - e)
        blocks_chr[chrom_y][block_id].anchored_split_breakpoints.add(breakpoint_id)
        breakpoint_id_good.add(breakpoint_id)
    FB.close()
    f_tmp_anchored_fa = work_dir + "_anchored_tmp_fa.fa"
    f_tmp_anchored_psl = work_dir + "_anchored_tmp.psl"
    F = open(f_tmp_anchored_fa, 'w')
    for breakpoint_id in breakpoint_id_good:
        breakpoint = breakpoints[breakpoint_id]
        if breakpoint.type_ == 'SM':
            F.write('>' + str(breakpoint_id) + '\n' + breakpoint.return_right_seq() + '\n')
        elif breakpoint.type_ == 'MS':
            F.write('>' + str(breakpoint_id) + '\n' + breakpoint.return_left_seq() + '\n')
    F.close()
    os.system("blat -stepSize=3 -minScore=12 -minMatch=2 -minIdentity=90 " + file_anchored_seq + " " + f_tmp_anchored_fa + " " + f_tmp_anchored_psl)
    breakpoint_id_good = set()
    F = open(f_tmp_anchored_psl,'r')
    for line in F.readlines():
        if not re.match(r'^\d+\s', line):
            continue
        arr = line.split('\t')
        score, l  = int(arr[0]), int(arr[10])
        if l*0.9<=score:
            breakpoint_id_good.add(int(arr[9]))
    F.close()
    # if os.path.exists(f_tmp_fa):
    # os.remove(f_tmp_fa)
    # if os.path.exists(f_tmp_blat):
    # os.remove(f_tmp_blat)
    return breakpoint_id_good


def Find_candidate_genes(out_dir_name_now, file_ref_seq, breakpoint_id_good, breakpoints, blocks_chr, gene_co, thread):
    aligner = Align.PairwiseAligner()
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5
    aligner.match_score = 1
    aligner.mismatch_score = -0.5

    def if_similar(seq_left1, seq_right1, seq_left2, seq_right2, seq_mid1, seq_mid2, threshold):
        if len(seq_left1)== 0 or len(seq_left2)== 0 or len(seq_right1)== 0 or len(seq_right2)== 0:
            return False
        if len(seq_mid1) != 0 and len(seq_mid2) != 0:
            c = 0
            for i in range(min(len(seq_mid1),len(seq_mid2))):
                if seq_mid1[i] == seq_mid2[i]:
                    c += 1
            if c/min(len(seq_mid1),len(seq_mid2)) < threshold:
                return False
        elif (len(seq_mid1)>3 and len(seq_mid2) == 0) or (len(seq_mid2)>3 and len(seq_mid1) == 0):
            return False
        seq_left1, seq_left2 = seq_left1[::-1], seq_left2[::-1]
        c = 0
        for i in range(min(len(seq_left1),len(seq_left2))):
            if seq_left1[i] == seq_left2[i]:
                c += 1
        score_left = c/min(len(seq_left1),len(seq_left2))
        c = 0
        for i in range(min(len(seq_right1),len(seq_right2))):
            if seq_right1[i] == seq_right2[i]:
                c += 1
        score_right = c/min(len(seq_right1),len(seq_right2))
        if score_left > threshold and score_right > threshold:
            return True
        return False

    candidates = []
    f_tmp_fa = out_dir_name_now + "_tmp_fa.fa"
    f_tmp_psl = out_dir_name_now + "_tmp.psl"
    for blocks in blocks_chr.values():
        for block in blocks:
            for breakpoint_id in block.anchored_split_breakpoints:
                if len(breakpoints[breakpoint_id].other_breakpoints) == 0:
                    continue
                if breakpoint_id not in breakpoint_id_good:
                    continue
                for other_breakpoint in breakpoints[breakpoint_id].other_breakpoints:
                    if other_breakpoint[0] != block.chrom:
                        continue
                    breakpoint = breakpoints[breakpoint_id]
                    type_2 = breakpoint.type_
                    seq_left2 = breakpoint.return_left_seq()
                    seq_right2 = breakpoint.return_right_seq()
                    seq_mid2 = ''
                    cnt2 = breakpoint.cnt
                    split_reads2 = breakpoint.reads
                    if type_2 == 'SM' :
                        seq_left2 = seq_left2[other_breakpoint[-2]:]
                        if other_breakpoint[-1]!=0:
                            seq_mid2 = seq_left2[-1 * other_breakpoint[-1]:]
                            seq_left2 = seq_left2[:-1 * other_breakpoint[-1]]
                    elif type_2 == 'MS':
                        seq_mid2 = seq_right2[:other_breakpoint[-2]]
                        seq_right2 = seq_right2[other_breakpoint[-2]:]
                        if other_breakpoint[-1] !=0 :
                            seq_right2 = seq_right2[:-1 * other_breakpoint[-1]]
                    flag = 0
                    j = len(candidates) - 1
                    while j >= max(0, len(candidates) - 200):
                        candidate = candidates[j]
                        seq_left1, seq_right1, type_1, seq_mid1 = candidate.return_seq_left(), candidate.return_seq_right(), candidate.type_, candidate.return_seq_mid()
                        if type_1 == type_2 and if_similar(seq_left1, seq_right1, seq_left2, seq_right2, seq_mid1,seq_mid2, 0.9):
                            flag = 1
                            candidate.add_reads(breakpoint.breakpoint, other_breakpoint, seq_left2, seq_right2, seq_mid2, cnt2, block.reads, split_reads2 )
                            break
                        j -= 1
                    if flag == 0:
                        candidate_read = Candidate_reads(type_2)
                        candidate_read.add_reads(breakpoint.breakpoint, other_breakpoint, seq_left2, seq_right2,seq_mid2, cnt2, block.reads, split_reads2)
                        candidates.append(candidate_read)
    good_ids, anchored_gene = [], ''
    F = open(f_tmp_fa, 'w')
    for i, candidate in enumerate(candidates):
        F.write('>' + str(i) + '\n' + candidate.return_seq_left() + candidate.return_seq_mid() + candidate.return_seq_right() + '\n')
    F.close()
    if len(candidates) != 0:
        os.system("blat -stepSize=3 -minScore=20 -minMatch=3 -minIdentity=90 " + file_ref_seq + " " + f_tmp_fa + " " + f_tmp_psl)
        F = open(f_tmp_psl, 'r')
        lines = F.readlines()
        i = 0
        name = -1
        out_name = -1
        while i < len(lines):
            if not re.match(r'^\d+\s', lines[i]):
                i += 1
                continue
            j = 0
            bad_flag = 0
            good_flag = 0
            if name != -1:
                left_length = candidates[name].l_left
                mid_length = candidates[name].l_mid
                right_length = candidates[name].l_right
            while i + j <= len(lines):
                if i + j == len(lines):
                    if name != -1:
                        if bad_flag == 0 and good_flag == 3:
                            good_ids.append(name)
                        break
                arr = lines[i+j].split('\t')
                read_name, chrom, start, end = int(arr[9]), arr[13], int(arr[11]),int(arr[12])
                if read_name != name:
                    if name != -1:
                        if bad_flag == 0 and good_flag == 3:
                            good_ids.append(name)
                    name = read_name
                    break
                if (start < left_length * 0.5 and end > left_length * 1.5 + mid_length):
                    bad_flag = 1
                elif (start <= left_length*0.5 and left_length * 0.5 <= end <=left_length * 1.5 ):
                    if good_flag == 0:
                        good_flag = 1
                    elif good_flag == 2:
                        good_flag = 3
                elif ( left_length+mid_length - right_length*0.5<= start <= left_length+mid_length+right_length*0.5 and left_length+mid_length+right_length*0.5<= end):
                    if good_flag == 0:
                        good_flag = 2
                    elif good_flag == 1:
                        good_flag = 3
                j += 1
            i += j
            if i == len(lines) - 1:
                break
        F.close()
    cnt_max = 0
    candidates_new = []
    for i in range(len(candidates)):
        if i in good_ids:
            pos = candidates[i].find_max_pos()[0]
            candidates_new.append([candidates[i], pos])
            if pos[6] > cnt_max:
                cnt_max = pos[6]
    candidates = []
    i = 0
    while i < len(candidates_new):
        pos = candidates_new[i][1]
        if pos[6] * 10 < cnt_max:
            del candidates_new[i]
            continue
        j = 0
        while j < i:
            if pos[0] == candidates_new[j][1][0] and pos[0] == candidates_new[j][1][1] and pos[0] == \
                    candidates_new[j][1][2]:
                del candidates_new[i]
                break
            j += 1
        if j == i:
            candidates.append(candidates_new[i][0])
            i += 1
    # if os.path.exists(f_tmp_fa):
    # os.remove(f_tmp_fa)
    return candidates

def make_negative_file(out_dir_name_now, negative_samples, f_read, gene_co, file_ref_seq, homo_genes, thread, gene_name):
    def if_are_homo_genes(can_gene,homo_gene,homo_gene_list):
        if can_gene not in homo_gene_list:
            return False
        if homo_gene in homo_gene_list[can_gene]:
            return True
        else:
            return False

    def Inspect_name(name1,name2):
        if re.match(r'^IG',name1) or re.match(r'^IG',name2):
            return True
        if re.match(r'^ENSG', name1) or re.match(r'^ENSG', name2):
            return True
        if len(name1) < 3 or len(name2) < 3:
            return False
        if name1[0] == name2[0] and name1[1] == name2[1] and name1[2] == name2[2]:
            return True
        return False

    gene_name = gene_name.upper()
    file1 = out_dir_name_now+'_negative_all_seq2.txt'
    file2 = out_dir_name_now+'_negative_all_seq.psl'
    last_genes = set()
    X_last = ''
    Y_last = []
    O = open(file1, 'w')
    F = os.popen("samtools view -t "+thread+ ' ' + f_read).read()
    homo_gene_list = np.load(homo_genes,allow_pickle=True).item()
    c = 0
    length = {}
    for line in F.splitlines():
        arr = line.strip().split('\t')
        if len(arr) <= 15:
            continue
        X_strand = '+'
        if int(arr[1]) > 15 and bin(int(arr[1]))[-5] == '1':
            X_strand = '-'
        X = (arr[2],arr[3],X_strand,arr[5])
        Y = []
        for tag in arr[15:]:
            if re.match(r'SA:Z:',tag):
                Y = tag[5:].split(';')[0].split(',')
                break
        if len(Y) == 0:
            continue
        seq = arr[9]
        l = len(seq)
        if re.match(r'\d+M\d+S', X[3]) and re.match(r'\d+S\d+M', Y[3]):
            l_X = int(re.findall(r'(\d+)M', X[3])[0])
            l_Y = int(re.findall(r'(\d+)S', Y[3])[0])
            if abs(l_X - l_Y) > 5:
                continue
            if X[2] == '-':
                seq = reverse(seq)
                l_X = len(seq)-l_X
        elif re.match(r'\d+S\d+M', X[3]) and re.match(r'\d+M\d+S', Y[3]):
            l_X = int(re.findall(r'\d+S(\d+)M', X[3])[0])
            l_Y = int(re.findall(r'\d+M(\d+)S', Y[3])[0])
            if abs(l_X - l_Y) > 5:
                continue
            if X[2] == '-':
                seq = reverse(seq)
                l_X = len(seq)-l_X
        else:
            continue
        l = len(seq)
        if l in length:
            length[l] += 1
        else:
            length[l] = 1
        X_gene = gene_co.Find_gene(X[0], int(X[1]), int(X[1]) + 1)
        Y_gene = gene_co.Find_gene(Y[0], int(Y[1]), int(Y[1]) + 1)
        if X_gene[1] == gene_name or Y_gene[1] == gene_name:
            continue
        if X_gene[1] == '' or Y_gene[1] == '' or X_gene[1] == Y_gene[1] or Inspect_name(X_gene[1],Y_gene[1]) or if_are_homo_genes(X_gene[1], Y_gene[1],homo_gene_list):
            continue
        if (X_gene[1] == X_last and Y_gene[1] in Y_last) or (X_gene[1] in Y_last and Y_gene[1] == X_last):
            continue
        fusion_gene = X_gene[1] + '--' + Y_gene[1]
        fusion_gene_reverse = Y_gene[1] + '--' + X_gene[1]
        if fusion_gene in last_genes or fusion_gene_reverse in last_genes:
            continue
        else:
            last_genes.add(fusion_gene)
        if X_gene[1] == X_last:
            Y_last.append(Y_gene[1])
        elif Y_gene[1] == X_last:
            Y_last.append(X_gene[1])
        else:
            X_last, Y_last = X_gene[1], [Y_gene[1]]
        O.write('>'+str(c)+'$'+fusion_gene+'$'+str(l_X)+'\n'+seq+'\n')
        c += 1
    O.close()
    ml,mc = 0,0
    for l,c in length.items():
        if c > mc:
            ml,mc = l,c
    ms = int(ml*0.8)
    os.system("blat -stepSize=3 -minScore="+str(ms)+" -minMatch=3 -minIdentity=90 -maxGap=1 " + file_ref_seq + " " + file1 + " " + file2)
    bad_ids = set()
    F2 = open(file2,'r')
    for line in F2.readlines():
        if not re.match(r'^\d+\s',line):
            continue
        arr = line.split('\t')
        bad_ids.add(int(arr[9].split('$')[0]))
    F2.close()
    F3 = open(file1,'r')
    O2 = open(negative_samples,'a+')
    i = 0
    for line in F3.readlines():
        if re.match(r'^>',line):
            id_now,fusion_gene,h_spot = line.rstrip()[1:].split("$")
            if int(id_now) in bad_ids:
                flag = False
            else:
                flag, h_spot = True, int(h_spot)
        else:
            if flag:
                seq = line.rstrip()
                if len(seq) != ml:
                    continue
                O2.write(seq[:h_spot]+'H'+seq[h_spot:]+'\t'+fusion_gene+'\n')
    F3.close()
    O2.close()
    return


def get_test_reads(test_file,candidates_new,len_seq,state='b'):
    c = 0
    F = open(test_file,'w')
    for candidate_new in candidates_new:
        if state == 'sc':
            l_left,l_mid,l_right = len(candidates_new[1]),len(candidates_new[2]),len(candidates_new[3])
            type_,seq_left,seq_mid,seq_right = candidates_new[0],candidates_new[1],candidates_new[2],candidates_new[3]
        else:
            l_left,l_mid,l_right = candidate_new.l_left , candidate_new.l_mid , candidate_new.l_right
            type_,seq_left,seq_mid,seq_right = candidate_new.type_ , candidate_new.return_seq_left() , candidate_new.return_seq_mid() , candidate_new.return_seq_right()
        if type_ == 'MS':
            if l_left + l_mid + l_right < len_seq:
                need = len_seq - l_left  - l_mid - l_right
                seq = 'N'*(need-need//2) + seq_left+ 'H' + seq_mid + seq_right + 'N'*(need//2)
            elif l_left < (len_seq-l_mid) // 2:
                seq = seq_left + 'H' + seq_mid + seq_right[:len_seq-l_left-l_mid]
            elif l_right < (len_seq-l_mid) // 2:
                seq = seq_left[l_right+l_mid-len_seq:] + 'H' + seq_mid + seq_right
            else:
                l_right = (len_seq-l_mid) // 2
                seq = seq_left[l_right+l_mid-len_seq:] + 'H' + seq_mid + seq_right[:l_right]
        else:
            if l_left + l_mid + l_right < len_seq:
                need = len_seq - l_left  - l_mid - l_right
                seq = 'N'*(need//2) + seq_left + seq_mid+ 'H' + seq_right + 'N'*(need - need//2)
            elif l_left < (len_seq-l_mid) // 2:
                seq = seq_left +  seq_mid +'H' +  seq_right[:len_seq-l_left-l_mid]
            elif l_right < (len_seq-l_mid) // 2:
                seq = seq_left[l_right+l_mid-len_seq:] + seq_mid + 'H' +  seq_right
            else:
                l_right = (len_seq-l_mid) // 2
                seq = seq_left[l_right+l_mid-len_seq:] + seq_mid +'H' +  seq_right[:l_right]
        F.write(seq+'\t'+str(c)+'\n')
        c += 1
    F.close()
    return

def Final_fusion(out_dir_name_now,  candidates_new, gene_name, gene_co,not_fiter=True):
    def deal_pos(pos, gene_name, gene_co):
        target_breakpoint, chrom, other_breakpoint, strand = pos[0], pos[1], pos[2], pos[3]
        mid_seq = pos[10]
        other_gene = gene_co.Find_gene(chrom, other_breakpoint, other_breakpoint + 1)
        if pos[9] == 'SM':
            fusion_gene = other_gene[1] + '--' + gene_name
        else:
            fusion_gene = gene_name + '--' + other_gene[1]
        write_line = fusion_gene + '\t' + gene_name + '\t' + gene_name + ':' + str(target_breakpoint) + '\t' + \
                     other_gene[1] + ':' + other_gene[0] + '\t' + chrom + ':' + str(other_breakpoint)
        return write_line, mid_seq

    def deal_pos2(type_, pos, gene_name, gene_co):
        target_breakpoint, chrom, other_breakpoint, strand = pos[0], pos[1], pos[2], pos[3]
        other_gene = gene_co.Find_gene(chrom, other_breakpoint, other_breakpoint + 1)
        if type_ == 'SM':
            fusion_gene = other_gene[1] + ':' + chrom + ':' + str(other_breakpoint) + '--' + gene_name + ':' + str(
                target_breakpoint)
        else:
            fusion_gene = gene_name + ':' + str(target_breakpoint) + '--' + other_gene[1] + ':' + chrom + ':' + str(
                other_breakpoint)
        return fusion_gene

    file_tmp_fa = out_dir_name_now + "_tmp_fa.fa"
    file_tmp_sam = out_dir_name_now + "_tmp_sam.sam"
    file_tmp_bed = out_dir_name_now + "_tmp.bed"
    file_abridged_output = out_dir_name_now + '_fusion_predictions_abridged.txt'
    file_output = out_dir_name_now + '_fusion_predictions.txt'
    FA = open(file_abridged_output, 'w')
    FO = open(file_output, 'w')
    if not_fiter:
        FA.write("Fusion_gene" + "\t" + "Anchored_gene_X" + "\t" + "X_clip_location" + "\t" + "Partner_gene_Y" + "\t" + "Y_clip_location" + "\t" + "Spanning_read_count" + "\t" + "Breakpoint_read_count" + "\n")
        FO.write("Fusion_gene" + "\t" + "Anchored_gene_X" + "\t" + "X_clip_location" + "\t" + "Partner_gene_Y" + "\t" +  "Y_clip_location" + "\t" + "Spanning_read_count" + "\t" + "Breakpoint_read_count" + "\t" + "Spanning_reads" + "\t" + "Breakpoint_reads" + "\t" + "Breakpoint_site_reads_1" + '\t' + 'Breakpoint_site_reads_2' + '\t' + 'Homo_genes'+"\tLeft_sequence\tMiddle_sequence\tRight_sequence"+'\n')
    else:
        FA.write("Fusion_gene" + "\t" + "Anchored_gene_X" + "\t" + "X_clip_location" + "\t" + "Partner_gene_Y" + "\t" + "Y_clip_location" + "\t"+"Natural_score"+"\t" + "Spanning_read_count" + "\t" + "Breakpoint_read_count" + "\n")
        FO.write("Fusion_gene" + "\t" + "Anchored_gene_X" + "\t" + "X_clip_location" + "\t" + "Partner_gene_Y" + "\t" +  "Y_clip_location" + "\t"+"Natural_score"+"\t" + "Spanning_read_count" + "\t" + "Breakpoint_read_count" + "\t" + "Spanning_reads" + "\t" + "Breakpoint_reads" + "\t" + "Breakpoint_site_reads_1" + '\t' + 'Breakpoint_site_reads_2' + '\t' + 'Homo_genes'+"\tLeft_sequence\tMiddle_sequence\tRight_sequence"+'\n')
    pos_and_mid = []
    for candidate_new in candidates_new:
        pos, max_id = candidate_new.find_max_pos()
        write_line, seq_mid = deal_pos(pos, gene_name, gene_co)
        target_breakpoint, chrom, other_breakpoint, strand = pos[0], pos[1], pos[2], pos[3]
        type_ = candidate_new.type_
        if (target_breakpoint, chrom, other_breakpoint, strand) in pos_and_mid:
            continue
        pos_and_mid.append((target_breakpoint, chrom, other_breakpoint, strand))
        spanning_reads = list(set(candidate_new.spanning_reads))
        split_reads = list(set(candidate_new.split_reads))
        if len(spanning_reads) * 5 < len(split_reads) or len(split_reads) * 5 < len(spanning_reads):
            continue
        if len(spanning_reads) == 0 and len(split_reads) == 0:
            continue
        if not_fiter:
            FA.write(write_line + '\t'+ str(len(spanning_reads)) + '\t' + str(len(split_reads)) + '\n')
        else:
            FA.write(write_line + '\t'+str(candidate_new.score)+'\t'+ str(len(spanning_reads)) + '\t' + str(len(split_reads)) + '\n')
        other_fusion_genes = []
        for i, pos in enumerate(candidate_new.pos):
            if i != max_id:
                other_fusion_genes.append(deal_pos2(type_, pos, gene_name, gene_co))
        seqs =  candidate_new.return_seq_left()+'\t' +  candidate_new.return_seq_mid() + '\t' +candidate_new.return_seq_right()
        if not_fiter:
            FO.write(write_line + '\t' + str(len(spanning_reads)) + '\t' + str(len(split_reads)) + '\t' + ';'.join(spanning_reads) + '\t' + ';'.join(split_reads) + '\t' + ';'.join(other_fusion_genes) +seqs+ '\n')
        else:
            FO.write(write_line + '\t'+str(candidate_new.score)+'\t' + str(len(spanning_reads)) + '\t' + str(len(split_reads)) + '\t' + ';'.join(spanning_reads) + '\t' + ';'.join(split_reads) + '\t' + ';'.join(other_fusion_genes) +seqs+ '\n')
    FA.close()
    FO.close()
    #if os.path.exists(file_tmp_fa):
        #os.remove(file_tmp_fa)
    # if os.path.exists(file_tmp_sam):
    # os.remove(file_tmp_sam)
    # if os.path.exists(file_tmp_bed):
    # os.remove(file_tmp_bed)
    return

