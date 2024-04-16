import os
import argparse

import numpy as np
import pandas as pd
import scanpy as sc
#from anndata import AnnData
import anndata as ad

import torchmetrics

import pytorch_lightning as pl
from pytorch_lightning.callbacks.progress import RichProgressBar
from pytorch_lightning.loggers import CSVLogger, TensorBoardLogger
from pytorch_lightning.callbacks import ModelCheckpoint, EarlyStopping

import torch
import torch.nn as nn
import torch.optim as optim
import torch.distributions as D
import torch.nn.functional as F
from torch.utils.data import DataLoader, random_split

import scanpy.external as sce

from sklearn.cluster import KMeans
from sklearn.mixture import GaussianMixture

from flowsom import flowsom as flowsom
from sklearn.cluster import AgglomerativeClustering

from sklearn.metrics import auc
from sklearn.metrics import f1_score
from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score
from sklearn.metrics import confusion_matrix
from sklearn.metrics import precision_recall_curve

from sklearn.ensemble import RandomForestClassifier

from utils import core, yaml_config_hook

os.environ["PYTORCH_CUDA_ALLOC_CONF"] = "max_split_size_mb:20000"

parser = argparse.ArgumentParser(description = "Create parameters for starling input")

parser.add_argument(
    "--cohort", type=str, required=True,
    help="which cohort did the data come from?",
)

parser.add_argument(
    "--seg_type", type=str, required=True,
    help="which mask type to use",
)

parser.add_argument(
    "--cofactor", type=int, required=True,
    help="denumuator of arcsinh input data",
)

parser.add_argument(
    "--sample_size", type=int, required=True,
    help="subsample data size (in thousands)",
)

parser.add_argument(
    "--noise_model", type=int, required=True,
    help="normal/0 or student distribution/2",
)

parser.add_argument(
    "--singlet_prop", type=float, required=True,
    help="prop of singlet in synehetic data",
)

parser.add_argument(
    "--lambda_val", type=float, required=True,
    help="regulization term between observed and simulated data",
)

parser.add_argument(
    "--cell_size", type=int, required=True,
    help="include cell size?",
)

parser.add_argument(
    "--relax_rule", type=int, required=True,
    help="relax rule",
)

parser.add_argument(
    "--repetition", type=int, required=True,
    help="number of experiments for each setting",
)

parser.add_argument(
    "--output_theta1", type=str, required=True,
    help="theta1 file name",
)

args = parser.parse_args()

cohort = args.cohort
seg_type = args.seg_type
sample_size = args.sample_size

cofactor = args.cofactor
dist_option = args.noise_model
singlet_prop = args.singlet_prop
model_cell_size = args.cell_size
model_overlap = args.relax_rule
model_regularizer = args.lambda_val
repetition = args.repetition

output_theta1 = args.output_theta1

#cohort = 'basel'
#seg_type = 'cp'
#cofactor = 5
#model_cell_size = 1
#dist_option = 2
#singlet_prop = 0.7
#model_overlap = 1 
#model_regularizer = 1
#sample_size = 10000

OPATH = '/home/campbell/yulee/project/st/{}/analysis/_2bm/{}/'.format(cohort, seg_type)
#channels = yaml_config_hook('/home/campbell/yulee/project/st/{}/channel.yaml'.format(cohort))['pretty_channels']

fn = '/home/campbell/yulee/project/st/{}/exp_mat/{}/exp_mat.h5ad'.format(cohort, seg_type)
if seg_type == 'dc':
    fn = '/home/campbell/yulee/project/st/{}/exp_mat/{}/exp_mat_corrected.h5ad'.format(cohort, seg_type)

df = ad.read_h5ad(fn)

df.var_names = np.array(df.var['marker'])

if cohort == 'basel':    
    exclude = ['Argon', 'RT2', 'RT3', 'RT4', 'RT5', 'RT6', 'RT7', 'RT8', 'RT9', 'DNA1', 'DNA2']
elif cohort == 'meta':
    exclude = ['Argon', 'Xe126', 'I127', 'Xe131', 'Xe134', 'Ce140', 'DNA1', 'DNA2', 'Hg202', 'Pb204', 'Pb206', 'Pb207', 'Pb208']
elif cohort == 'tonsil':
    exclude = ['DNA1', 'DNA2']
else:
    print("cohort is not one of basel/meta/tonsil")

df = df[:, ~df.var_names.isin(exclude)]

if cofactor > 0:
    df.X = np.arcsinh(df.X / cofactor)

sample_size = sample_size * 1000
#idx1 = np.random.choice(df.shape[0], size = sample_size, replace = False) ## for likelihood
idx = np.random.choice(df.shape[0], size = 2 * sample_size, replace = False) ## for train & valid fake data

tr_data = df[idx[:sample_size]]
va_data = df[idx[sample_size:(2*sample_size)]]

Theta1 = []

k = None
for j, initial_clustering_method in enumerate(["PG", 'FS', "KM"]):

    #tr_data = sample_data.tr_h5ad.copy()
    #tr_data = df.copy()

    print(cohort, repetition, initial_clustering_method)
    ONAME = "ia{}_nc{}_cf{}_nm{}_sp{}_lv{}_cs{}_rr{}_r{}".format(initial_clustering_method, k, cofactor, dist_option, 
                                                                singlet_prop, model_regularizer, model_cell_size, model_overlap, repetition)
    
    init_cen, init_var, init_label = core.init_clustering(tr_data.X, initial_clustering_method, k, OPATH + "fs/in_" + ONAME + ".csv")

    if model_cell_size == 1:     
        init_s = []; init_sv = []
        for c in range(init_cen.shape[0]):
            init_s.append(tr_data.obs.iloc[np.where(init_label == str(c))]['area'].mean())
            init_sv.append(tr_data.obs.iloc[np.where(init_label == str(c))]['area'].var())
        #print(init_s)
        #print(init_sv)
    else:
        init_s = None; init_sv = None

    if initial_clustering_method == "PG": ## update k
            k = init_cen.shape[0]
    
    tr_data.obs[initial_clustering_method] = init_label
    st = core.ST(tr_data, va_data, np.array(init_cen), np.array(init_var), np.array(init_s), np.array(init_sv), model_cell_size, dist_option, singlet_prop, model_overlap, model_regularizer)
    
    cb_progress = RichProgressBar()
    cb_early_stopping = EarlyStopping(monitor = 'train_loss', mode = 'min', verbose = False)

    if model_cell_size and model_overlap:
        log_tb = TensorBoardLogger(save_dir = OPATH + 'logs/', name = ONAME)
        trainer = pl.Trainer(max_epochs = 150, accelerator = 'auto', devices = 'auto', callbacks = [cb_progress, cb_early_stopping], logger=[log_tb], default_root_dir = OPATH + 'logs/')
    else:
        trainer = pl.Trainer(max_epochs = 150, accelerator = 'auto', devices = 'auto', callbacks = [cb_progress, cb_early_stopping], default_root_dir = OPATH + 'logs/')
    
    trainer.fit(st)
    st.result()

    Theta1.append(st.model_params)

    ## run random forest
    if model_cell_size == 1:

        train_set = st.train_df.datasets[2]
        valid_set = st.val_df.datasets[2]
    
        train_set1 = torch.hstack((st.train_df.datasets[2], st.train_df.datasets[3].reshape(-1,1)))
        valid_set1 = torch.hstack((st.val_df.datasets[2], st.val_df.datasets[3].reshape(-1,1)))

        train_label = st.train_df.datasets[4]
        valid_label = st.val_df.datasets[4]

        rf = RandomForestClassifier()
        rf.fit(train_set, train_label)
        p_singlet = rf.predict_proba(train_set)[:,1]
        st.rf_cen, st.rf_var, rf_label = core.rf_clustering(train_set, p_singlet, initial_clustering_method, k, OPATH + "fs/rf_" + ONAME + ".csv")
        
        st.ad_tr.obs['rf_{}'.format(initial_clustering_method)] = rf_label
        st.rf_f1 = f1_score(valid_label, np.where(rf.predict_proba(valid_set)[:,1] > 0.5, 1, 0))
        
        rf1 = RandomForestClassifier()
        rf1.fit(train_set1, train_label)
        p_singlet1 = rf1.predict_proba(train_set1)[:,1]
        st.rf1_cen, st.rf1_var, rf1_label = core.rf_clustering(train_set1, p_singlet1, initial_clustering_method, k, OPATH + "fs/rf1_" + ONAME + ".csv")
        
        st.ad_tr.obs['rf1_{}'.format(initial_clustering_method)] = rf1_label
        st.rf1_f1 = f1_score(valid_label, np.where(rf1.predict_proba(valid_set1)[:,1] > 0.5, 1, 0))
        
        fake_pred_loader = torch.utils.data.DataLoader(core.ConcatDataset(st.train_df.datasets[2], st.train_df.datasets[3].reshape(-1,1)), batch_size = 1000, shuffle = False)
        fake_singlet_prob, _, _, _ = core.predict(fake_pred_loader, st.model_params, dist_option, model_cell_size, model_overlap, 0.5)
    else:
        train_set = st.train_df.datasets[1]
        valid_set = st.val_df.datasets[1]
        
        train_label = st.train_df.datasets[2]
        valid_label = st.val_df.datasets[2]

        rf = RandomForestClassifier()
        rf.fit(train_set, train_label)
        p_singlet = rf.predict_proba(train_set)[:,1]
        st.rf_cen, st.rf_var, rf_label = core.rf_clustering(train_set, p_singlet, initial_clustering_method, k, OPATH + "fs/rf_" + ONAME + ".csv")
        
        st.ad_tr.obs['rf_{}'.format(initial_clustering_method)] = rf_label
        st.rf_f1 = f1_score(valid_label, np.where(rf.predict_proba(valid_set)[:,1] > 0.5, 1, 0))

        st.rf1_cen = None
        st.rf1_var = None
        st.rf1_label = None
        st.rf1_f1 = None

        fake_pred_loader = torch.utils.data.DataLoader(train_set, batch_size = 1000, shuffle = False)
        fake_singlet_prob, _, _, _ = core.predict(fake_pred_loader, st.model_params, dist_option, 0, 0, 0.5)

    st.st_f1 = f1_score(valid_label, torch.where(fake_singlet_prob > 0.5, 1, 0).cpu())

    torch.save(st, OPATH + ONAME)

## output model parameters
torch.save(Theta1, output_theta1)


    