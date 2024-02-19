import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
from torch.utils import data
import random
import os
from sklearn.metrics import  roc_auc_score, accuracy_score
import torch.distributions as td
import math

os.environ['CUDA_LAUNCH_BLOCKING'] = '1'

def try_gpu(i):
    i = int(i)
    if i != -1 and torch.cuda.device_count() >= i + 1:
        return torch.device(f'cuda:{i}')
    else:
        return torch.device('cpu')

def setup_seed(seed):
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    np.random.seed(seed)
    random.seed(seed)
    torch.backends.cudnn.deterministic = True

setup_seed(1024)

class Loss(nn.Module):
    def __init__(self,a1,a2):
        super(Loss, self).__init__()
        self.a1 = a1
        self.a2 = a2
        self.ce_loss = nn.CrossEntropyLoss(reduction='mean')
        self.kd_loss = F.kl_div
    def forward(self,outputs,output,y):
        ce_loss1 = self.ce_loss(outputs[0],y)
        ce_loss2 = self.ce_loss(outputs[1],y)
        ce_loss = self.ce_loss(output,y)
        kd_loss1 = self.kd_loss(torch.log(outputs[0]),output,reduction = 'batchmean')
        kd_loss2 = self.kd_loss(torch.log(outputs[1]),output,reduction = 'batchmean')
        loss = self.a1*ce_loss1+ (1-self.a1)*kd_loss1 + self.a2*ce_loss2+ (1-self.a2)*kd_loss2 + ce_loss
        return loss

class MLP_2(nn.Module):
    def __init__(self,input_dim,mid_dim,output_dim,dropout = 0.0):
        super(MLP_2, self).__init__()
        self.fc1 = nn.Linear(input_dim,mid_dim)
        self.dropout = nn.Dropout(dropout)
        self.fc2 = nn.Linear(mid_dim,output_dim)
        self.relu = nn.ReLU()
    def forward(self,X):
        X = self.dropout(self.relu(self.fc1(X)))
        X = self.fc2(X)
        return X

class Block(nn.Module):
    def __init__(self,input_dim,embed_dim,output_dim,window,maxpool_dim):
        super(Block, self).__init__()
        self.normal_layer1 = nn.BatchNorm1d(embed_dim)
        padding = window//2
        self.conv1 = nn.Conv1d(input_dim,embed_dim,window,padding=padding)
        self.conv2 = nn.Conv1d(embed_dim,output_dim,window,padding=padding)
        self.relu = nn.ReLU()
        self.avgpool = nn.AvgPool1d(maxpool_dim,stride = maxpool_dim)
    def forward(self,X):
        X = X.transpose(-1,-2)
        X = self.relu(self.normal_layer1(self.conv1(X)))
        X = self.relu(self.conv2(X))
        X = self.avgpool(X)
        X = X.transpose(-1,-2)
        return X

class Classify(nn.Module):
    def __init__(self,input_dim,embed_dim,len_seq,class_shrink_dim,num_class,dropout):
        super(Classify, self).__init__()
        self.prj = nn.Linear(input_dim,input_dim//class_shrink_dim)
        self.num_class = num_class
        self.flatten = nn.Flatten()
        self.classify = MLP_2(len_seq*input_dim//class_shrink_dim,embed_dim,num_class,dropout=0.2)
    def forward(self,X,t1=1):
        X = self.prj(X)
        X = self.flatten(X)
        X = self.classify(X)/t1
        X = F.softmax(X,dim=1)
        return X

class TransformerEncoder(nn.Module):
    def __init__(self, input_dim, len_seq, hidden_dim, num_layers, num_heads):
        super(TransformerEncoder, self).__init__()
        self.len_seq = len_seq
        self.input_embedding = nn.Linear(input_dim,hidden_dim)
        self.position_encoding = nn.Embedding(len_seq, hidden_dim)
        nn.init.normal_(self.position_encoding.weight, std=0.02)
        encoder_layer = nn.TransformerEncoderLayer(hidden_dim, num_heads,dropout=0.0,batch_first=True)
        self.transformer_encoder = nn.TransformerEncoder(encoder_layer, num_layers=num_layers)
        self.relu = nn.ReLU()

    def forward(self, X):
        X = self.input_embedding(X)
        positions = torch.arange(self.len_seq, device=X.device).unsqueeze(0)
        position_encoding = self.position_encoding(positions)
        X = X + position_encoding
        #X = X.transpose(0, 1)
        X = self.relu(self.transformer_encoder(X))
        #X = X.transpose(0, 1)
        return X

class Model(nn.Module):
    def __init__(self,input_dim,block_dim,embed_dim,class_dim,window,maxpool_dim,class_shrink_dim,transformer_dim,num_class,len_seq,dropout):
        super(Model, self).__init__()
        self.embed_dim = embed_dim
        self.num_class = num_class
        self.relu = nn.ReLU()
        self.input_embedding = nn.Linear(input_dim,embed_dim)
        self.block1 = Block(embed_dim,block_dim,embed_dim,window,maxpool_dim)
        self.classify1 = Classify(embed_dim,class_dim,len_seq//maxpool_dim,class_shrink_dim,num_class,dropout)
        self.block2 = Block(embed_dim,block_dim,embed_dim,window,maxpool_dim)
        self.classify2 = Classify(embed_dim,class_dim,len_seq//(maxpool_dim**2),class_shrink_dim,num_class,dropout)
        self.transformer = TransformerEncoder(embed_dim,len_seq//(maxpool_dim**2),transformer_dim,1,2)
        self.classify3 = Classify(transformer_dim,class_dim,len_seq//(maxpool_dim**2),class_shrink_dim,num_class,0.2)
    def forward(self,X):
        X = self.relu(self.input_embedding(X))
        X = self.block1(X)
        output1 = self.classify1(X,0.25)
        X = self.block2(X)
        output2 = self.classify2(X,0.25)
        X = self.transformer(X)
        output3 = self.classify3(X)
        return (output1,output2),output3

def make_train_file(positive_file,negative_file,negative_file_tra,negative_file_tes,positive_file_tra,positive_file_tes):
    Fn = open(negative_file,'r')
    lines_n = Fn.readlines()
    l_negative = len(lines_n)
    Fp = open(positive_file,'r')
    lines_p = Fp.readlines()
    l_positive = len(lines_p)
    l = min(l_negative,l_positive)
    l_seq = len(lines_n[0].split('\t')[0])
    list_negative = list(range(l_negative))
    random.shuffle(list_negative)
    O1 = open(negative_file_tra,'w')
    for i in list_negative[:int(0.7*l)]:
        O1.write(lines_n[i])
    O1.close()
    O2 = open(negative_file_tes,'w')
    for i in list_negative[int(0.7*l):l]:
        O2.write(lines_n[i])
    O2.close()
    list_positive = list(range(l_positive))
    random.shuffle(list_positive)
    l_p_seq = len(lines_p[0].split('\t')[0])//2
    O3 = open(positive_file_tra,'w')
    for i in list_positive[:int(0.7*l)]:
        seq,fusion_gene = lines_p[i].split('\t')
        O3.write(seq+'\t'+fusion_gene)
    O3.close()
    O4 = open(positive_file_tes,'w')
    for i in list_positive[int(0.7*l):l]:
        seq,fusion_gene = lines_p[i].split('\t')
        O4.write(seq+'\t'+fusion_gene)
    O4.close()
    Fn.close()
    Fp.close()
    return


def read_lines(file_name):
    turn_dic = {'A':0,'T':1,'G':2,'C':3,'H':4,'D':5}
    F=open(file_name,"r")
    seqs = list()
    seqs_ori = list()
    for line in F.readlines():
        arr=line.split("\t")
        seq =arr[0].upper().replace('\n','')
        s = np.zeros((len(seq),6))
        seqs_ori.append(seq)
        i = 0
        while i < len(seq):
            if seq[i] in turn_dic:
                s[i][turn_dic[seq[i]]] = 1.0
            i += 1
        seqs.append(s)
    seqs = torch.from_numpy(np.array(seqs))
    return seqs

def data_load(X_tra,X_tes, y_tra,y_tes, batch_size, num_works):
    dataset_tra, dataset_tes = data.TensorDataset(X_tra, y_tra),data.TensorDataset(X_tes, y_tes)
    return data.DataLoader(dataset_tra, batch_size, shuffle=True, num_workers=num_works), \
           data.DataLoader(dataset_tes, batch_size, shuffle=True, num_workers=num_works)

def train_epoch(net, train_iter, loss, updater, device, num):
    y_hat_all = []
    y_true_all = []
    y_hat_all_arg = []
    y_hat_all1 = []
    y_hat_all_arg1 = []
    y_hat_all2 = []
    y_hat_all_arg2 = []
    if isinstance(net, torch.nn.Module):
        net.train()
    for X, y in train_iter:
        X, y = X.to(device),y.to(device)
        outputs,y_hat = net(X)
        l = loss(outputs,y_hat,y)
        updater.zero_grad()
        l.backward()
        updater.step()
        y_hat_arg = y_hat.argmax(axis=1)
        y_hat_all.extend(y_hat[:,1].detach().cpu().numpy())
        y_hat_all_arg.extend(y_hat_arg.detach().cpu().numpy())
        y_true_all.extend(y.detach().cpu().numpy())

    accuracy = accuracy_score(y_true_all, y_hat_all_arg)
    auc = roc_auc_score(y_true_all, y_hat_all)
    print(f'Train {num} epoch: 0 Accuracy:{format(accuracy,".3f")}, Auc:{format(auc,".3f")}')

@torch.no_grad()
def test_epoch(net, test_iter, device, num):
    y_hat_all = []
    y_true_all = []
    y_hat_all_arg = []

    for X, y in test_iter:
        X, y = X.to(device),y.to(device)
        _,y_hat = net(X)
        y_hat_arg = y_hat.argmax(axis=1)
        y_hat_all.extend(y_hat[:,1].detach().cpu().numpy())
        y_hat_all_arg.extend(y_hat_arg.detach().cpu().numpy())
        y_true_all.extend(y.detach().cpu().numpy())

    accuracy = accuracy_score(y_true_all, y_hat_all_arg)
    auc = roc_auc_score(y_true_all, y_hat_all)
    print(f'Test {num} epoch 0 : Accuracy:{format(accuracy,".3f")}, Auc:{format(auc,".3f")}')

    return auc

def train(net, train_iter, test_iter, loss, num_epochs, updater, device, model_file):
    mac_auc = 0
    for epoch in range(num_epochs):
        train_epoch(net, train_iter, loss, updater, device, epoch)
        test_auc = test_epoch(net, test_iter, device, epoch)
        if test_auc > mac_auc:
            torch.save(net.state_dict(),model_file)
            mac_auc = test_auc
    return

@torch.no_grad()
def test(net,X,device):
    X = X.to(device)
    y_hat = net(X)[1][:,1]
    y_hat = y_hat.detach().cpu().numpy().tolist()
    return y_hat

def Train_model(model_out_name, positive_samples,negative_samples, model_file,gpu_number):
    device = try_gpu(gpu_number)
    np.random.seed(0)
    torch.manual_seed(0)
    torch.cuda.manual_seed_all(0)
    negative_file_tra = model_out_name + 'negative_tra_seq.txt'
    negative_file_tes = model_out_name + 'negative_tes_seq.txt'
    positive_file_tra = model_out_name + 'positive_tra_seq.txt'
    positive_file_tes = model_out_name + 'positive_tes_seq.txt'
    make_train_file(positive_samples,negative_samples,negative_file_tra,negative_file_tes,positive_file_tra,positive_file_tes)
    X_n_tra = read_lines(negative_file_tra)
    X_n_tes = read_lines(negative_file_tes)
    X_p_tra = read_lines(positive_file_tra)
    X_p_tes = read_lines(positive_file_tes)
    X_tra = torch.cat((X_p_tra, X_n_tra), 0)
    X_tes = torch.cat((X_p_tes, X_n_tes), 0)
    y_tra = torch.cat((torch.ones(X_p_tra.shape[0],dtype=torch.int64),torch.zeros(X_n_tra.shape[0],dtype=torch.int64)),dim=0)
    y_tes = torch.cat((torch.ones(X_p_tes.shape[0],dtype=torch.int64),torch.zeros(X_n_tes.shape[0],dtype=torch.int64)),dim=0)
    List_tra = list(range(X_tra.shape[0]))
    List_tes = list(range(X_tes.shape[0]))
    random.shuffle(List_tra)
    random.shuffle(List_tes)
    X_tra = X_tra[List_tra,:,:]
    y_tra = y_tra[List_tra]
    X_tes = X_tes[List_tes,:,:]
    y_tes = y_tes[List_tes]
    batch_size = 256
    num_works = 1
    lr = 0.0001
    num_epochs = 30
    loss = Loss(0.5,0.5)
    len_seq = X_n_tra.shape[1]
    input_dim = 6
    block_dim = 256
    embed_dim = 256
    class_dim = 256
    window = 3
    maxpool_dim = 3
    class_shrink_dim = 4
    transformer_dim = 128
    num_class = 2
    dropout = 0.2
    net = Model(input_dim,block_dim,embed_dim,class_dim,window,maxpool_dim,class_shrink_dim,transformer_dim,num_class,len_seq,dropout)
    if os.path.exists(model_file):
        net.load_state_dict(torch.load(model_file))
    net = net.to(device).double()
    updater = torch.optim.Adam(net.parameters(), lr=lr)
    mac_auc = 0
    train_iter, test_iter = data_load(X_tra, X_tes, y_tra, y_tes, batch_size,1)
    for epoch in range(num_epochs):
        train_epoch(net, train_iter, loss, updater, device, epoch)
        test_auc = test_epoch(net, test_iter, device, epoch)
        if test_auc > mac_auc:
            torch.save(net.state_dict(),model_file)
            mac_auc = test_auc
    return  len_seq

def Test_model(test_file, model_file,gpu_number):
    device = try_gpu(gpu_number)
    X = read_lines(test_file)
    len_seq = X.shape[1]
    input_dim = 6
    block_dim = 256
    embed_dim = 256
    class_dim = 256
    window = 3
    maxpool_dim = 3
    class_shrink_dim = 4
    transformer_dim = 128
    num_class = 2
    dropout = 0.2
    net = Model(input_dim,block_dim,embed_dim,class_dim,window,maxpool_dim,class_shrink_dim,transformer_dim,num_class,len_seq,dropout)
    if os.path.exists(model_file):
        net.load_state_dict(torch.load(model_file))
    net = net.to(device).double()
    scores = test(net,X,device)
    return scores

