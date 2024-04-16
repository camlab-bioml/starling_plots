
import numpy as np
import pandas as pd
import tifffile

#!pip install DeepCell
from deepcell.applications import Mesmer
from deepcell.utils.plot_utils import create_rgb_image, make_outline_overlay

from utils import yaml_config_hook

import ssl
ssl._create_default_https_context = ssl._create_unverified_context

def mask1(raw, channels):
    
    marker_index = {'cyto_channels': [], 'nuc_channels': []}
    for cyto in channels['cyto_channels']: #range(len(cyto_channels)):
        #print(np.where(cyto==np.array(channels))[0])
        marker_index['cyto_channels'].append(np.where(cyto==np.array(channels['all_channels']))[0][0])

    for nuc in channels['nuc_channels']: #range(len(cyto_channels)):
        #print(np.where(nuc==np.array(channels))[0])
        marker_index['nuc_channels'].append(np.where(nuc==np.array(channels['all_channels']))[0][0])
        
    for p in range(raw.shape[0]):
        tmp = raw[p,:,:].copy()
        raw[p,:,:] = (tmp - tmp.min()) / (tmp.max()-tmp.min())
            
    cyto_image = raw[marker_index['cyto_channels'],:,:].mean(0)
    nuc_image = raw[marker_index['nuc_channels'],:,:].mean(0)

    new_img1 = np.stack([nuc_image, cyto_image], axis=2)

    cyto_image = (cyto_image-cyto_image.min())/(cyto_image.max()-cyto_image.min())
    nuc_image = (nuc_image-nuc_image.min())/(nuc_image.max()-nuc_image.min())
    
    new_img2 = np.stack([nuc_image, cyto_image], axis=2)
    
    return new_img1, new_img2

app = Mesmer()

cohort = 'tonsil'

if cohort == 'basel':
    channels = yaml_config_hook('/home/campbell/yulee/project/starling/config/basel_channel.yaml')
elif cohort == 'meta':
    channels = yaml_config_hook('/home/campbell/yulee/project/starling/config/meta_channel.yaml')
elif cohort == 'tonsil':
    channels = yaml_config_hook('/home/campbell/yulee/project/starling/config/tonsil_channel.yaml')

column_names = np.hstack((['id', 'x', 'y', 'area', 'area_convex', 'has_neighbor'], channels['all_channels']))
    
sample_fn = "/home/campbell/yulee/project/st/{}/tiff_mask_sample_file_aligned.tsv".format(cohort)
samples = pd.read_csv(sample_fn, sep="\t")

for i in range(len(samples)):
    
    print(i)
    raw = tifffile.imread(samples.iloc[i]['expression_tiff']).astype(np.float64)
    
    _, m2 = mask1(raw, channels)
    new_img = np.stack([m2])
    
    try:
        segmentation_predictions = app.predict(new_img, image_mpp=1)
        segmentation_predictions_nuc = app.predict(new_img, image_mpp=1, compartment='nuclear')
        np.save(samples.iloc[i]['mask_dc'].split(".tiff")[0] + '_input.npy', new_img)
    except:
        print(i)
        continue

    ## save masks
    tifffile.imwrite(samples.iloc[i]['mask_dc'], segmentation_predictions)
    tifffile.imwrite(samples.iloc[i]['mask_dc'].split(".tiff")[0] + '_nuc.tiff', segmentation_predictions_nuc)