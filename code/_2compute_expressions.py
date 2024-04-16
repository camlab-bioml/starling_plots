import os
import argparse
import tifffile
import numpy as np
import pandas as pd
from skimage.measure import regionprops
from skimage.segmentation import find_boundaries

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
    "--sample_id", type=int, required=True,
    help="sample id",
)

parser.add_argument(
    "--output_fn", type=str, required=True,
    help="output file name",
)

args = parser.parse_args()

cohort = args.cohort
seg_type = args.seg_type
sample_id = args.sample_id
output_fn = args.output_fn


def cell_morphology(view, exp, sample_id, cell_id, has_neighbor, multiple_nuclei, nuclear_larger_than30):
    props = regionprops(np.where(view, 1, 0))
    if len(props) == 0:
        return np.hstack((-1, -1, 0., 0., 0., 0., 0., 0., 0., np.zeros(exp.shape[0])))
    else:
        x, y = props[0].centroid
        return np.hstack((sample_id, cell_id, x, y, props[0].area, props[0].area_convex, has_neighbor, multiple_nuclei, nuclear_larger_than30, exp[:, view == 1].mean(1)))

#i = 4
#cohort = 'tonsil'
#seg_type = 'dc'

sample_fn = "/home/campbell/yulee/project/st/{}/tiff_mask_sample_file_aligned.tsv".format(cohort)
samples = pd.read_csv(sample_fn, sep="\t")

#for sample_id in range(samples.shape[0]):
print(sample_id)

if cohort == 'meta' and sample_id in [109, 112, 129, 131, 151, 168, 187, 188, 235, 238, 317, 350, 369, 420, 424, 426, 488]:
    np.save(output_fn, None)
else:
    exp = tifffile.imread(samples.iloc[sample_id]['expression_tiff']).astype(np.float64)
    if seg_type == 'dc':
        #msk = np.squeeze(tifffile.imread(samples.iloc[sample_id]['mask_dc']).astype(int), axis=(0,3))
        nuc_cell = np.squeeze(tifffile.imread(samples.iloc[sample_id]['mask_dc'].split(".tiff")[0] + '_nuc.tiff').astype(int), axis=(0,3))
        whole_cell = np.squeeze(tifffile.imread(samples.iloc[sample_id]['mask_dc']).astype(int), axis=(0,3))
    else:
        whole_cell = tifffile.imread(samples.iloc[sample_id]['mask_cp']).astype(int)
        
    ii = 0
    cell_ids = np.unique(whole_cell)[1:]
    exp_mat = np.zeros((len(cell_ids), 9 + exp.shape[0])) #cell, channel
    for cell_id in cell_ids:
        tmp_mask = np.zeros_like(whole_cell, dtype=int)
        tmp_mask[whole_cell == cell_id] = 1

        outer_boundary = find_boundaries(tmp_mask, mode='outer', connectivity=1)
        any_id_overlap = whole_cell[outer_boundary]
        has_neighbor = len(np.unique(any_id_overlap[np.where((any_id_overlap != 0) & (any_id_overlap != cell_id))])) > 0
        
        #print(has_neighbor)
        
        multiple_nuclei = -1
        nuclear_larger_than30 = -1
        if seg_type == 'dc':
            #nuc_cell_id = nuc_cell[np.where(whole_cell == cell_id)]
            #remove_zero = len(np.unique(nuc_cell_id[np.where(nuc_cell_id > 0)[0]]))
            
            #nuc_cell_id = nuc_cell[tmp_mask == 1]

            nuc_cell_id = nuc_cell[tmp_mask == 1]
            candidate_id = pd.Series(nuc_cell_id).value_counts(normalize=True)
            remove_background = candidate_id[candidate_id.index > 0]
            multiple_nuclei = len(remove_background)
            nuclear_larger_than30 = len(remove_background[remove_background > 0.3].index)

        exp_mat[ii] = cell_morphology(tmp_mask, exp, sample_id, cell_id, has_neighbor, multiple_nuclei, nuclear_larger_than30) ## original
        ii += 1
        
        #np.save('/home/campbell/yulee/project/st/{}/mat/dc/{}.npy'.format(cohort, sample_id), exp_mat)
        np.save(output_fn, exp_mat)