import pickle
from sklearn import svm
from sklearn.datasets import make_classification
from sklearn.model_selection import train_test_split,cross_validate
from sklearn.datasets import make_classification
import tensorflow as tf
from tensorflow.keras.utils import to_categorical
from sklearn.model_selection import cross_val_score
import numpy as np
import random
from sklearn.metrics import  roc_auc_score, accuracy_score

def read_lines_positive(file_name,lacks):
    c = [0,0,0,0]
    l = 0
    F=open(file_name,"r")
    seqs = list()
    seqs_ori = list()
    j = 0
    for line in F.readlines():
        arr=line.split("\t")
        seq =arr[0].upper().replace('\n','')
        if j < len(lacks):
            lack = lacks[j]
            left = lack//2
            right = lack - left
            if right != 0:
                seq = 'N'*left + seq[left:-1*right]+ 'N'*right
            j += 1
        s = np.zeros(len(seq)*5)
        i = 0
        i = seq.find('H')
        s[i*5],s[i*5+1],s[i*5+2],s[i*5+3],s[i*5+4] = 1.0,1.0,1.0,1.0,1.0
        i -= 1
        while i >= 0:
            l += 1
            if seq[i] == 'N':
                i -= 1
                continue
            if seq[i] == 'A':
                c[0] += 1
                s[i*5],s[i*5+4] = 1.0,c[0]/l
            elif seq[i] == 'T':
                c[1] += 1
                s[i*5+1],s[i*5+4] = 1.0,c[1]/l
            elif seq[i] == 'G':
                c[2] += 1
                s[i*5+2],s[i*5+4] = 1.0,c[2]/l
            elif seq[i] == 'C':
                c[3] += 1
                s[i*5+3],s[i*5+4] = 1.0,c[3]/l
            i -= 1
        i= seq.find('H') + 1
        c = [0,0,0,0]
        l = 0
        while i < len(seq):
            l += 1
            if seq[i] == 'N':
                i += 1
                continue
            if seq[i] == 'A':
                c[0] += 1
                s[i*5],s[i*5+4] = 1.0,c[0]/l
            elif seq[i] == 'T':
                c[1] += 1
                s[i*5+1],s[i*5+4] = 1.0,c[1]/l
            elif seq[i] == 'G':
                c[2] += 1
                s[i*5+2],s[i*5+4] = 1.0,c[2]/l
            elif seq[i] == 'C':
                c[3] += 1
                s[i*5+3],s[i*5+4] = 1.0,c[3]/l
            i += 1
        seqs.append(s)
    seqs = np.array(seqs)
    return seqs

def read_lines_negative(file_name,lacks):
    c = [0,0,0,0]
    l = 0
    F=open(file_name,"r")
    seqs = list()
    seqs_ori = list()
    for line in F.readlines():
        arr=line.split("\t")
        seq =arr[0].upper().replace('\n','')
        lack = 102 - len(seq)
        lacks.append(lack)
        left = lack//2
        right = lack - left
        seq = 'N'*left + seq +'N'*right
        s = np.zeros(len(seq)*5)
        i = 0
        i = seq.find('H')
        s[i*5],s[i*5+1],s[i*5+2],s[i*5+3],s[i*5+4] = 1.0,1.0,1.0,1.0,1.0
        i -= 1
        while i >= 0:
            l += 1
            if seq[i] == 'N':
                i -= 1
                continue
            if seq[i] == 'A':
                c[0] += 1
                s[i*5],s[i*5+4] = 1.0,c[0]/l
            elif seq[i] == 'T':
                c[1] += 1
                s[i*5+1],s[i*5+4] = 1.0,c[1]/l
            elif seq[i] == 'G':
                c[2] += 1
                s[i*5+2],s[i*5+4] = 1.0,c[2]/l
            elif seq[i] == 'C':
                c[3] += 1
                s[i*5+3],s[i*5+4] = 1.0,c[3]/l
            i -= 1
        i= seq.find('H') + 1
        c = [0,0,0,0]
        l = 0
        while i < len(seq):
            l += 1
            if seq[i] == 'N':
                i += 1
                continue
            if seq[i] == 'A':
                c[0] += 1
                s[i*5],s[i*5+4] = 1.0,c[0]/l
            elif seq[i] == 'T':
                c[1] += 1
                s[i*5+1],s[i*5+4] = 1.0,c[1]/l
            elif seq[i] == 'G':
                c[2] += 1
                s[i*5+2],s[i*5+4] = 1.0,c[2]/l
            elif seq[i] == 'C':
                c[3] += 1
                s[i*5+3],s[i*5+4] = 1.0,c[3]/l
            i += 1
        seqs.append(s)
    seqs = np.array(seqs)
    return seqs

save_dir="/data/Anchor_Fusion/model/"
np.random.seed(1122)
random.seed(1122)
tes_p_file = '/data/Anchor_Fusion/positive/positive_tes_seq.txt'
tes_n_file = '/data/Anchor_Fusion/negative/negative_tes_seq.txt'
tra_p_file = '/data/Anchor_Fusion/positive/positive_tra_seq.txt'
tra_n_file = '/data/Anchor_Fusion/negative/negative_tra_seq.txt'
lacks= []
Simu_for_Tra = read_lines_negative(tra_n_file,lacks)
Good_for_Tra = read_lines_positive(tra_p_file,lacks)
Simu_for_Tst = read_lines_negative(tes_n_file,lacks)
Good_for_Tst = read_lines_positive(tes_p_file,lacks)
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
Tra_y = Tra_y[LIST,1]

LIST = list(range(Tst_x.shape[0]))
np.random.shuffle(LIST)
Tst_x = Tst_x[LIST,:]
Tst_y = Tst_y[LIST,1]

scores_file = '/data/Anchor_Fusion/model/svm_scores.npy'
labels_file = '/data/Anchor_Fusion/model/svm_labels.npy'
filename = '/data3/yuanxilu/scFusion-main/Anchor_Fusion/model/svm_models.pkl'
LIST = list(range(Tra_x.shape[0]))
np.random.shuffle(LIST)
Tra_x = Tra_x[LIST,:]
Tra_y = Tra_y[LIST]
model = svm.SVC(probability=True)
scoring = ['accuracy', 'precision', 'recall']
results = cross_validate(model, Tra_x, Tra_y, cv=5, scoring=scoring, return_estimator=True)
best_model_index = results['test_accuracy'].argmax()
best_model = results['estimator'][best_model_index]
with open(filename, 'wb') as file:
    pickle.dump(best_model.get_params(), file)
y_pred = best_model.predict(Tst_x)
y_pred_prob = best_model.predict_proba(Tst_x)[:, 1]
accuracy = sum(Tst_y == y_pred) / len(Tst_y)
auc_score = roc_auc_score(Tst_y, y_pred_prob)
np.save(scores_file,y_pred_prob)
np.save(labels_file,Tst_y)
