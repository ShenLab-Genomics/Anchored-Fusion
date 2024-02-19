from keras.models import Sequential,Model
from keras.layers import Embedding,Dropout,Bidirectional,Flatten,Dense,LSTM,TimeDistributed, Activation,Input,concatenate
from keras.callbacks import ModelCheckpoint,CSVLogger
from keras.layers import Conv1D, GlobalAveragePooling1D, MaxPooling1D
from tensorflow.keras.optimizers import Adam
import numpy as np
from tensorflow.keras.utils import to_categorical
import os
import sys
import random
import tensorflow as tf
from sklearn.metrics import  roc_auc_score, accuracy_score

def read_lines_positive(file_name,lacks):
    def turn_tensor(seq,chr_index):
        res = np.zeros((len(seq)))
        for i, c in enumerate(seq):
            if c in chr_index:
                res[i] = chr_index[c]
        return res
    
    char_index = {'A':'0','T':'1','G':'2','C':'3','H':'4'}
    F=open(file_name,"r")
    c = 0
    seqs = list()
    j = 0
    for line in F.readlines():
        arr=line.split("\t")
        seq=arr[0].upper().replace('\n','')
        if j < len(lacks):
            lack = lacks[j]
            left = lack//2
            right = lack - left
            if right != 0:
                seq = 'N'*left + seq[left:-1*right]+ 'N'*right
            j += 1
        c += 1
        s = turn_tensor(seq,char_index)
        seqs.append(s)
    seqs = np.array(seqs)
    List = list(range(c))
    return seqs

def read_lines_negative(file_name,lacks):
    def turn_tensor(seq,chr_index):
        res = np.zeros((len(seq)))
        for i, c in enumerate(seq):
            if c in chr_index:
                res[i] = chr_index[c]
        return res
    
    char_index = {'A':'0','T':'1','G':'2','C':'3','H':'4'}
    F=open(file_name,"r")
    c = 0
    seqs = list()
    for line in F.readlines():
        arr=line.split("\t")
        seq=arr[0].upper().replace('\n','')
        lack = 61 - len(seq)
        lacks.append(lack)
        left = lack//2
        right = lack - left
        seq = 'N'*left + seq +'N'*right
        c += 1
        s = turn_tensor(seq,char_index)
        seqs.append(s)
    seqs = np.array(seqs)
    List = list(range(c))
    return seqs

def read_lines(file_name):
    def turn_tensor(seq,chr_index):
        res = np.zeros((len(seq)))
        for i, c in enumerate(seq):
            if c in chr_index:
                res[i] = chr_index[c]
        return res
    
    char_index = {'A':'0','T':'1','G':'2','C':'3','H':'4'}
    F=open(file_name,"r")
    seqs = list()
    for line in F.readlines():
        arr=line.split("\t")
        seq=arr[0].upper().replace('\n','')
        s = turn_tensor(seq,char_index)
        seqs.append(s)
    seqs = np.array(seqs)
    return seqs

def Cla_LSTM():
# 搭模型
    INPUT1 = Input(shape=(61,))
    #INPUT2 = Input(shape=(60,1))
    INPUT1_enco = Embedding(5,5,input_length=61)(INPUT1)

    # MERGE = merge((INPUT1_enco, INPUT2),mode='concat',concat_axis=-1)
    # MERGE = concatenate([INPUT1,INPUT2])
    LSTM1 = Bidirectional(LSTM(32,return_sequences=True),merge_mode='concat')(INPUT1_enco)
    DROP1 = Dropout(0.5)(LSTM1)
    
    LSTM2 = Bidirectional(LSTM(64,return_sequences=True),merge_mode='concat')(DROP1)
    DROP2 = Dropout(0.5)(LSTM2)
    
    LSTM3 = Bidirectional(LSTM(128,return_sequences=True),merge_mode='concat')(DROP2)
    DROP3 = Dropout(0.5)(LSTM3)
    
    '''
    LSTM4 = Bidirectional(LSTM(64,return_sequences=True),merge_mode='concat')(DROP3)
    DROP4 = Dropout(0.5)(LSTM4)    
    
    LSTM5 = Bidirectional(LSTM(128,return_sequences=True),merge_mode='concat')(DROP4)
    DROP5 = Dropout(0.5)(LSTM5) 

    LSTM6 = Bidirectional(LSTM(128,return_sequences=True),merge_mode='concat')(DROP5)
    DROP6 = Dropout(0.5)(LSTM6) 

    LSTM7 = Bidirectional(LSTM(36,return_sequences=True),merge_mode='concat')(DROP6)
    DROP7 = Dropout(0.5)(LSTM7) 
    '''
    LSTM8 = Bidirectional(LSTM(256, return_sequences=False), merge_mode='concat')(DROP3)
    DENSE1 = Dense(256)(LSTM8)
    DENSE2 = Dense(2)(DENSE1)
    
    ACT1 = Activation('softmax')(DENSE2)
    # model = Model(inputs=[INPUT1,INPUT2],outputs= ACT1)
    model = Model(inputs=INPUT1,outputs= ACT1)
    #model.summary()
    return model

save_dir="/data/Anchor_Fusion/model/"
epochoutdir =  '/data/Anchor_Fusion/k562_output/model_dir/bi-lstm/'
np.random.seed(1122)
random.seed(1122)
tf.random.set_seed(1122)

tes_p_file = '/data/Anchor_Fusion/k562_output/model_dir/positive_tes_seq.txt'
tes_n_file = '/data/Anchor_Fusion/k562_output/model_dir/negative_tes_seq.txt'
tra_p_file = '/data/Anchor_Fusion/k562_output/model_dir/positive_tra_seq.txt'
tra_n_file = '/data/Anchor_Fusion/k562_output/model_dir/negative_tra_seq.txt'
Simu_for_Tra = read_lines(tra_n_file)
Good_for_Tra = read_lines(tra_p_file)
Simu_for_Tst = read_lines(tes_n_file)
Good_for_Tst = read_lines(tes_p_file)
List_n_tra = list(range(Simu_for_Tra.shape[0]))
Simu_for_Tra = Simu_for_Tra[List_n_tra,:]
l_p_tra,l_n_tra = Good_for_Tra.shape[0],Simu_for_Tra.shape[0]
Tra_x = np.squeeze(np.concatenate((Good_for_Tra,Simu_for_Tra),axis=0))
Tra_y = np.concatenate((np.ones((Good_for_Tra.shape[0],1)),np.zeros((Simu_for_Tra.shape[0],1))), axis=0)
Tst_x = np.squeeze(np.concatenate((Good_for_Tst,Simu_for_Tst),axis=0))
Tst_y = np.concatenate((np.ones((Good_for_Tst.shape[0], 1)), np.zeros((Simu_for_Tst.shape[0], 1))), axis=0)
Tra_y = to_categorical(Tra_y)
Tst_y = to_categorical(Tst_y)
LIST = list(range(Tra_x.shape[0]))
np.random.shuffle(LIST)

Tra_x = Tra_x[LIST,:]
Tra_y = Tra_y[LIST,:]

LIST = list(range(Tst_x.shape[0]))
np.random.shuffle(LIST)
Tst_x = Tst_x[LIST,:]
Tst_y = Tst_y[LIST,:]


best_auc = 0.5
batch_size = 512
epoch_num = 200
np.random.seed(1122)
weightfile = '/data/Anchor_Fusion/k562_output/model_dir/new_model.hdf5'
LIST = list(range(Tra_x.shape[0]))
np.random.shuffle(LIST)
Tra_x = Tra_x[LIST,:]
Tra_y = Tra_y[LIST,:]
ADAM = Adam(learning_rate=0.0001)
model_checkpoint = ModelCheckpoint(filepath=epochoutdir + '/RetrainWeight-{epoch:03d}.hdf5', verbose=1, monitor='val_loss', save_best_only=True)
model_checkpoint2 = ModelCheckpoint(filepath=epochoutdir + '/RetrainWeight.hdf5', verbose=1, monitor='val_loss', save_best_only=True)
model.compile(loss='binary_crossentropy', optimizer=ADAM, metrics=['accuracy'])
csv_loger=CSVLogger('log.csv',append=True,separator=';')

# 训练模型
batch_size = 500
epochs = itere
np.random.seed(1122)
model.fit(x=Tra_x, y=Tra_y,batch_size=batch_size,epochs=epochs,validation_data=(Tst_x, Tst_y),verbose=1 ,callbacks=[model_checkpoint,model_checkpoint2, csv_loger], shuffle=False)
