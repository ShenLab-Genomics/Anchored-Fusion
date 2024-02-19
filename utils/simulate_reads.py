import os
import random

names = [('k562','BCR','BCR-ABL1'),('TE441T','DUX4','CIC-DUX4'),('nalm6','DUX4','DUX4-IGH'),('NCI_HCC660','TMPRSS2','TMPRSS2-ERG'),('NCIH3122','ALK','EML4-ALK'),('KM12','NTRK1','TPM3-NTRK1')]

for name in names:
    file = '/data/Anchor_Fusion/'+name[2]+'.fasta'
    F = open(file)
    s = ''
    for line in F.readlines()[1:]:
        s += line[:-1]
    l = len(s)
    dir_now = '/data/Anchor_Fusion/sim/'+name[0]+'/'
    if not os.path.exists(dir_now):
        os.system('mkdir '+dir_now)
    sizes = ['128','64','32','16','8','4','2']
    for size in sizes:
        n = int(size)*l//202
        if not os.path.exists(dir_now+size+'/'+size+'_new_1.fastq'):
          os.system('wgsim -d 200 -1 101 -2 101 -N '+str(n)+' '+file + ' '+ dir_now+size+'/'+size+'_new_1.fastq '+dir_now+size+'/'+size+'_new_2.fastq')
exit()
