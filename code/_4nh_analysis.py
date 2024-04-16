import numpy as np
import pandas as pd
import anndata as ad
import scanpy.external as sce
from sklearn.cluster import KMeans
from sklearn.mixture import GaussianMixture

from flowsom import flowsom as flowsom
from sklearn.cluster import AgglomerativeClustering

def plausibility_scores1(mat1, neg_pairs, pos_pairs, thresholds):

    score = []
    for c in range(mat1.shape[0]):
    
        tmp = []
        for pair in neg_pairs:
            x = np.array(mat1[pair[0]][c])
            y = np.array(mat1[pair[1]][c])
            tmp.append((x < thresholds[pair[0]][0]) | (y < thresholds[pair[1]][0]))
        
        for pair in pos_pairs:
            x1 = np.array(mat1[pair[0]][c])
            y1 = np.array(mat1[pair[1]][c])
            tmp.append((x1 > thresholds[pair[0]][1]) | (y1 < thresholds[pair[1]][0]))
    
        #score.append([c, sum(tmp)/len(tmp)])
        score.append(sum(tmp)/len(tmp))
    return score

def plausibility_scores2(mat1, anti, pairs, thresholds):

    out = []
    for pair in pairs:
        x_name = pair[0]; y_name = pair[1]
        x1 = np.array(mat1[x_name]); y1 = np.array(mat1[y_name])

        if anti:
            ## case 1: negative
            score1 = sum((x1 < thresholds[x_name][0]) | (y1 < thresholds[y_name][0]))/len(x1)
            out.append([x_name + '_' + y_name, score1, -1])
        else:
            ## case 2: positive
            score2 = sum((x1 > thresholds[x_name][1]) | (y1 < thresholds[y_name][0]))/len(x1)
            out.append([x_name + '_' + y_name, score2, 1])
    
    return pd.DataFrame(out, columns = ['pair', 'metric', 'type'])  

def ps(cofactor, df, neg_pairs, pos_pairs, thresholds, fn, save_center=0):
    
    out = []
    k = None
    for j, initial_clustering_method in enumerate(['PG', 'FS', 'KM']):
        
        try:
            c, v, l = init_clustering(df, initial_clustering_method, k, fn, cofactor)
        except:
            out.append([cofactor, initial_clustering_method, k, 0, 0])
            continue

        cl = pd.DataFrame(c, columns = df.columns)
        #if save_center:
        #    cl.to_csv(PATH + "/" + initial_clustering_method + "_c" + str(cofactor) + "_m" + str(mask_input) + "_r" + str(iter) + ".csv")

        if initial_clustering_method == "PG": ## update k
            k = c.shape[0]
            print(k)
        
        for i, thres in enumerate(thresholds):
            cluster_score = np.array(plausibility_scores1(cl, neg_pairs, pos_pairs, thres)).mean()
            neg_score = plausibility_scores2(cl, 1, neg_pairs, thres)
            pos_score = plausibility_scores2(cl, 0, pos_pairs, thres)
            score = (neg_score.loc[:,'metric'].sum() + pos_score.loc[:,'metric'].sum()) / (neg_score.shape[0] + pos_score.shape[0])

            out.append([cofactor, initial_clustering_method, k, score, cluster_score])
    
    return out

def init_clustering(X, initial_clustering_method, k, training_file_fn=None, cofactor=None):
    
    if initial_clustering_method == 'KM':
        kms = KMeans(k).fit(X)
        init_l = kms.labels_
        init_label_class = np.unique(init_l)

        init_c = kms.cluster_centers_
        init_v = np.array([np.array(X)[init_l == c,:].var(0) for c in init_label_class])
        
    elif initial_clustering_method == 'GMM':
        gmm = GaussianMixture(n_components = k, covariance_type = 'diag', reg_covar=0.5).fit(X)
        init_l = gmm.predict(X)

        init_c = gmm.means_
        init_v = gmm.covariances_
                 
    elif initial_clustering_method == 'PG':
        init_l, _, _ = sce.tl.phenograph(X)
        
        ## save phenograph centers
        nc = len(np.unique(init_l))
        init_v = np.zeros((nc, X.shape[1]))
        init_c = np.zeros((nc, X.shape[1]))
        for c in range(nc):
            init_v[c,:] = X[init_l==c].var(0)
            init_c[c,:] = X[init_l==c].mean(0)
    elif initial_clustering_method == 'FS':

        ## needs to output to csv first
        #ofn = OPATH + "fs_" + ONAME + ".csv"
        #pd.DataFrame(mat.numpy()).to_csv(training_file_fn)
        #if model_cell_size:
        fsom = flowsom(training_file_fn, if_fcs=False, if_drop=True, drop_col=['Unnamed: 0'])
        #else:
        #    fsom = flowsom(training_file_fn, if_fcs=False, if_drop=True, drop_col=['Unnamed: 0', 'area'])
        
        #if cofactor != 0:
        #    fsom.df = np.arcsinh(fsom.df/cofactor)
            
        fsom.som_mapping(50, # x_n: e.g. 100, the dimension of expected map
            50, # y_n: e.g. 100, the dimension of expected map
            fsom.df.shape[1],
            1, # sigma: e.g 1, the standard deviation of initialized weights
            0.5, # lr: e.g 0.5, learning rate
            1000, # batch_size: 1000, iteration times
            tf_str=None, # string, e.g. hlog', None, etc - the transform algorithm
            if_fcs=False # bool, whethe the imput file is fcs file. If not, it should be a csv file
            # seed = 10, for reproducing
        )
        start = k; fsom_num_cluster = 0
        while fsom_num_cluster < k:
            #print(nc, start, fsom_nc)
            fsom.meta_clustering(AgglomerativeClustering, min_n=start, max_n=start, verbose=False, iter_n=10) # train the meta clustering for cluster in range(40,45)  

            fsom.labeling()
            #fsom.bestk # the best number of clusters within the range of (min_n, max_n)
            fsom_class = np.unique(fsom.df['category'])
            fsom_num_cluster = len(fsom_class)
            start += 1
    
        fsom_labels = np.array(fsom.df['category'])

        i = 0
        init_l = np.zeros(fsom.df.shape[0])
        init_vars = np.zeros((len(fsom_class), fsom.df.shape[1]))
        init_centers = np.zeros((len(fsom_class), fsom.df.shape[1]))
        for row in fsom_class:
            init_l[fsom_labels==row] = i
            init_vars[i,:] = fsom.df[fsom_labels==row].var(0)
            init_centers[i,:] = fsom.df[fsom_labels==row].mean(0)
            i += 1

        init_v = init_vars[:,:-1]
        init_c = init_centers[:,:-1]
        #os.remove(ofn)
        
    return init_c, init_v, init_l.astype(str)

pos_pairs = [['CD45', 'CD3'],
        ['CD45', 'CD20'],
        ['CD45', 'CD68'],
        ['pCK', 'CK5'],
        ['pCK', 'CK7'],
        ['pCK', 'CK8n18'],
        ['pCK', 'CK19']]

cohort = 'basel'
seg_type = 'dc'
sample_size = 10000
PATH = "/home/campbell/yulee/project/st/{}/analysis/_1nh/{}".format(cohort, seg_type)

fn = '/home/campbell/yulee/project/st/{}/exp_mat/{}/exp_mat.h5ad'.format(cohort, seg_type)
if seg_type == 'dc':
    fn = '/home/campbell/yulee/project/st/{}/exp_mat/{}/exp_mat_corrected.h5ad'.format(cohort, seg_type)

df = ad.read_h5ad(fn)

print(df.shape)
print(cohort)
print(seg_type)

thres025 = pd.DataFrame(np.quantile(df.X, [0.25, 0.75], axis=0), columns = df.var['marker'])
thres125 = pd.DataFrame(np.quantile(np.arcsinh(df.X), [0.25, 0.75], axis=0), columns = df.var['marker'])
thres525 = pd.DataFrame(np.quantile(np.arcsinh(df.X/5), [0.25, 0.75], axis=0), columns = df.var['marker'])

#n = df[df.obs['has_neighbor'] == 1]
#if n.shape[0] > sample_size:
#    select_n_cells = np.random.choice(n.shape[0], size = sample_size, replace = False)
#else:
#    select_n_cells = np.random.choice(n.shape[0], size = sample_size, replace = True)

nn = df[df.obs['has_neighbor'] != 1]
if nn.shape[0] > sample_size:
    select_nn_cells = np.random.choice(nn.shape[0], size = sample_size, replace = False)
else:
    select_nn_cells = np.random.choice(nn.shape[0], size = sample_size, replace = True)

#n0_sample = pd.DataFrame(n[select_n_cells].X, columns = np.array(df.var['marker']))
nn0_sample = pd.DataFrame(nn[select_nn_cells].X, columns = np.array(df.var['marker']))

if cohort == 'basel':
    n0_sample = n0_sample.drop(columns=['Argon', 'RT2', 'RT3', 'RT4', 'RT5', 'RT6', 'RT7', 'RT8', 'RT9', 'DNA1', 'DNA2'])
    nn0_sample = nn0_sample.drop(columns=['Argon', 'RT2', 'RT3', 'RT4', 'RT5', 'RT6', 'RT7', 'RT8', 'RT9', 'DNA1', 'DNA2'])

    neg_pairs = [['CD3', 'CD20'], 
        ['CD3', 'CD31'],
        ['CD3', 'CD68'], 
        ['CD3','ECadherin'],
        ['CD3', 'pCK'], 
        ['CD20', 'CD31'],
        ['CD20', 'CD68'], 
        ['CD20','ECadherin'],
        ['CD20', 'pCK'], 
        ['CD31', 'CD68'], 
        ['CD31','ECadherin'],
        ['CD45', 'pCK'], 
        ['CD68','ECadherin'],
        ['pCK', 'Vimentin']]
elif cohort == 'meta':
    #n0_sample = n0_sample.drop(columns=['Argon', 'Xe126', 'I127', 'Xe131', 'Xe134', 'Ce140', 'DNA1', 'DNA2', 'Hg202', 'Pb204', 'Pb206', 'Pb207', 'Pb208'])
    nn0_sample = nn0_sample.drop(columns=['Argon', 'Xe126', 'I127', 'Xe131', 'Xe134', 'Ce140', 'DNA1', 'DNA2', 'Hg202', 'Pb204', 'Pb206', 'Pb207', 'Pb208'])

    neg_pairs = [['CD3', 'CD20'], 
        ['CD3', 'vWFCD31'],
        ['CD3', 'CD68'], 
        ['CD3','ECadherin'],
        ['CD3', 'pCK'], 
        ['CD20', 'vWFCD31'],
        ['CD20', 'CD68'], 
        ['CD20','ECadherin'],
        ['CD20', 'pCK'], 
        ['vWFCD31', 'CD68'], 
        ['vWFCD31','ECadherin'],
        ['CD45', 'pCK'], 
        ['CD68','ECadherin'],
        ['pCK', 'Vimentin']]
else:
    n0_sample = n0_sample.drop(columns=['DNA1', 'DNA2'])
    nn0_sample = nn0_sample.drop(columns=['DNA1', 'DNA2'])

    neg_pairs = [
        ['CD3','CD20'], ['CD3','CD31'], ['CD3','CD68'], ['CD3','E_Cadherin'],
        ['CD4','CD8'], ['CD4','CD20'], ['CD4','CD31'], ['CD4','E_Cadherin'],
        ['CD8','CD20'], ['CD8','CD31'], ['CD8','CD68'], ['CD8','E_Cadherin'],
        ['CD20','CD31'], ['CD20','CD68'], ['CD20','E_Cadherin'],
        ['CD31','CD68'], ['CD31','E_Cadherin'], ['CD68','E_Cadherin']]

    pos_pairs = [
        ['CD3','CD4'], ['CD3','CD8'], ['CD45','CD3'], ['CD45','CD4'],
        ['CD45','CD8'], ['CD45','CD20'], ['CD45','CD68']]
'''
for i in range(50):

    print("n {}".format(i))

    n1_sample = np.arcsinh(n0_sample)
    n5_sample = np.arcsinh(n0_sample/5)

    sfn_n0 = PATH + "/fs/n0_" + str(i) + ".csv"
    sfn_n1 = PATH + "/fs/n1_" + str(i) + ".csv"
    sfn_n5 = PATH + "/fs/n5_" + str(i) + ".csv"
    onfn = PATH + "/n_" + str(i) + ".csv"

    pd.DataFrame(n0_sample).to_csv(sfn_n0)
    pd.DataFrame(n1_sample).to_csv(sfn_n1)
    pd.DataFrame(n5_sample).to_csv(sfn_n5)

    n0 = ps(0, n0_sample, neg_pairs, pos_pairs, [thres025], sfn_n0, 0)
    n1 = ps(1, n1_sample, neg_pairs, pos_pairs, [thres125], sfn_n1, 0)
    n5 = ps(5, n5_sample, neg_pairs, pos_pairs, [thres525], sfn_n5, 0)
    n_res = np.vstack((n0, n1, n5))

    pd.DataFrame(n_res, columns = ['Cofactor', 'Initial', 'Cluster', 'Score', 'Score_c']).to_csv(onfn)

'''
for i in range(50):

    print("nn {}".format(i))

    nn1_sample = np.arcsinh(nn0_sample)
    nn5_sample = np.arcsinh(nn0_sample/5)

    sfn_nn0 = PATH + "/fs/nn0_" + str(i) + ".csv"
    sfn_nn1 = PATH + "/fs/nn1_" + str(i) + ".csv"
    sfn_nn5 = PATH + "/fs/nn5_" + str(i) + ".csv"
    onnfn = PATH + "/nn_" + str(i) + ".csv"

    pd.DataFrame(nn0_sample).to_csv(sfn_nn0)
    pd.DataFrame(nn1_sample).to_csv(sfn_nn1)
    pd.DataFrame(nn5_sample).to_csv(sfn_nn5)

    nn0 = ps(0, nn0_sample, neg_pairs, pos_pairs, [thres025], sfn_nn0, 0)
    nn1 = ps(1, nn1_sample, neg_pairs, pos_pairs, [thres125], sfn_nn1, 0)
    nn5 = ps(5, nn5_sample, neg_pairs, pos_pairs, [thres525], sfn_nn5, 0)
    nn_res = np.vstack((nn0, nn1, nn5))

    pd.DataFrame(nn_res, columns = ['Cofactor', 'Initial', 'Cluster', 'Score', 'Score_c']).to_csv(onnfn)
