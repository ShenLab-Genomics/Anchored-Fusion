import os
import random

names =  [('k562','BCR'),('TE441T','DUX4'),('nalm6','DUX4'),('NCI_HCC660','TMPRSS2'),('KM12','NTRK1'),('KM12','NTRK1')]
for name in names:
    dir_now = '/data/Anchor_Fusion/downsample/'+name[0]+'/'
    last_size = '0'
    for size in ['2','4','8','16','32','64','128']:
        input_file1 = dir_now+ last_size+'_1.fastq'
        input_file2 = dir_now + last_size + '_2.fastq'
        output_file1 = dir_now + size + '_1.fastq'
        output_file2 = dir_now + size + '_2.fastq'
        if not os.path.exists(output_file1):
            F1 = open(input_file1,'r')
            F2 = open(input_file2,'r')
            O1 = open(output_file1,'w')
            O2 = open(output_file2,'w')
            l_F1 = F1.readlines()
            l_F2 = F2.readlines()
            n_l = len(l_F1)//4
            s_l =  random.sample(list(range(n_l)),n_l//2)
            for i in s_l:
                O1.write(''.join(l_F1[i*4:i*4+4]))
                O2.write(''.join(l_F2[i * 4:i * 4 + 4]))
            del s_l
            F1.close()
            F2.close()
            O1.close()
            O2.close()
        last_size = size

exit()
