import sys
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import seaborn
import anndata as ad
from scipy import stats

import matplotlib
import matplotlib.pyplot as plt

from sklearn import preprocessing
from sklearn import linear_model
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import GridSearchCV
from sklearn.preprocessing import StandardScaler

from sklearn.metrics import homogeneity_completeness_v_measure

font = {# 'family' : 'serif', # Times (source: https://matplotlib.org/tutorials/introductory/customizing.html)
        'family': 'sans-serif', # Helvetica
#         'family': 'monospace',
#         'weight' : 'bold',
        'size'   : 14}
matplotlib.rc('font', **font) 
text = {'usetex': False}
matplotlib.rc('text', **text)
monospace_font = {'fontname':'monospace'}

np.random.seed(3)
sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor="white")

def load_data(load_lsi=True, n_neighbors=20, n_pca=30):
    # motif name
    with open('data/motif_names.txt') as f:
        motif_list = f.readlines()
    motif_list = [item.strip() for item in motif_list]
    
    # matrices
    motif_reg_data_dir = "results/motif_regression_output/"
    motif_dir = f"{motif_reg_data_dir}/matchR_out/"
    adata_atac = sc.read(motif_reg_data_dir + "ArchR_PBMC_peak_cell_matrix.mtx", cache=True)
    adata_atac = adata_atac.T
    adata_atac.layers["raw"] = adata_atac.X.copy()
    adata_motif = sc.read(motif_dir + "MOTIF_PBMC_summit_motif_matrix.mtx", cache=True)
    
    # summit metadata
    matrix_summit_metadata = pd.read_csv(motif_reg_data_dir + 'ArchR_matrix_summit_metadata.csv', index_col=0)
    
    filtered_summit_range_df = pd.read_csv(motif_reg_data_dir + 'SUMMITS-all-atac-summits.bed', sep='\t', 
                                       names=['chrm', 'start', 'end', 'name', '.', '..'])

    chipseeker = pd.read_csv(motif_reg_data_dir + "ChIPseeker_PBMC_annotation.txt", sep = "\t")
    chipseeker = chipseeker.drop(columns=['score'])
    summit_df = pd.merge(matrix_summit_metadata, filtered_summit_range_df,
                         how="left", on=["name"]).merge(chipseeker, how='left', on='name')
    
    # add metadata to ann data
    summit_df.index = summit_df['name']
    adata_atac.var_names = summit_df.index
    adata_motif.obs_names = summit_df.index
    for i in ['name', 'chrm', 'start_x', 'end_x', 'annotation', 'SYMBOL', 'GENENAME']:
        adata_atac.var[i] = summit_df[i]

    adata_motif.var_names = motif_list
    
    peak_cells = pd.read_csv(motif_reg_data_dir + "ArchR_PBMC_cellnames_info.txt", sep = "\t")
    peak_cells.index = [i.split('#')[1] for i in list(peak_cells.index)]
    
    adata_atac.obs_names = peak_cells.index                                                                   
    for i in list(peak_cells):
        adata_atac.obs[i] = peak_cells[i]
    
    # PCA
    if load_lsi:
        lsi = pd.read_csv(motif_reg_data_dir + "ArchR_PBMC_Iterative_LSI.txt", sep = "\t", index_col = 0)
        lsi.index = [i.split('#')[1] for i in list(lsi.index)]
        adata_atac.obsm['X_pca'] = np.array(lsi)
    else:
#         sc.pp.scale(adata_atac, max_value=10)
        sc.tl.pca(adata_atac, n_comps = n_pca, svd_solver='arpack', zero_center=False, 
                  random_state=1, use_highly_variable = False)
    
    # KNN
    if load_lsi:
        sc.pp.neighbors(adata_atac, n_neighbors=n_neighbors, n_pcs=n_pca, metric = 'cosine')
    else:
        sc.pp.neighbors(adata_atac, n_neighbors=n_neighbors, n_pcs=n_pca, metric = 'euclidean')
        
    # UMAP    
    sc.tl.umap(adata_atac)
    
    # UMAP plot of the single cell data
    sc.pl.umap(adata_atac, color = 'Clusters')
#     adata_atac.X = adata_atac.layers["raw"]
    
    return adata_atac, adata_motif

# poissone correction
# https://github.com/theislab/scatac_poisson/tree/abadd31321e985f34cb0e03aa988b9058ccc8fe8

def plot_count_histogram(X):
    plt.hist(X, bins=int(np.max(X)), log=True)
    plt.xlabel('Read (summit values)')
    plt.ylabel('Count (# summits)')
    plt.grid()
    
def poisson_correction(adata):
    adata.X.data = np.ceil(adata.X.data/2)
    adata.layers["poisson_corrected"] = adata.X.copy()
    return adata

def check_poisson_assumption(X):
    C_mean_var = np.zeros((X.shape[0], 2))
    for i in range(X.shape[0]):
        x = X[i, :]
        C_mean_var[i, 0], C_mean_var[i, 1] = np.mean(x), np.var(x)
    return C_mean_var

def plot_poisson_assumpition(X):
    c_mean_var = check_poisson_assumption(X)
    plt.scatter(c_mean_var[:, 0], c_mean_var[:, 1], alpha=0.1)
    plt.plot([0, 10], [0, 10], '--', color='gray')
    plt.xscale("log")
    plt.yscale("log")
    plt.xlim([0.001, 10])
    plt.ylim([0.001, 10])
    
def most_frequent(List):
    return max(set(List), key = List.count)

def group_atac_sneha(adata, n_neighbors=20, max_allowed_overlap=0.9, N_ITER=20):
    mat = adata.obsp['distances'].todense()
    adata.obs['index'] = np.arange(len(adata.obs_names))

    # sample                                                                                                                                                          
    neighbors = {}
    neighbor_clusters = {}
    lib_sizes = {}
    counter = 0
    meta_counter = 0
    d_counter = {}
    # maintain check on which cell has already been selected                                                                                                          
    for obs_name in adata.obs['index']: d_counter[obs_name] = 0

    # max 20 iterations where check is applied for 90%                                                                                                                
    while counter < N_ITER:
        cell = np.random.choice(adata.obs['index'])
        if d_counter[cell] == 1:
            counter += 1
            continue

        # cell explored                                                                                                                                               
        d_counter[cell] = 1

        # get neighbors and append the current index                                                                                                                  
        indices = np.where(mat[cell, :] > 0)[1]
        indices = np.append(indices, cell)

        is_uniq = True
        # check overlap                                                                                                                                               
        for neighbor_k in neighbors:
            overlap = np.intersect1d(neighbors[neighbor_k], indices).shape[0] / n_neighbors
            if overlap >= max_allowed_overlap:
                is_uniq = False
                break

        # if overlap with any existing group is more than allowed overlap                                                                                             
        if not is_uniq:
            counter +=1
            continue

        # add to list                                                                                                                                                 
        neighbors["metacell_" + str(cell)] = indices
        lib_sizes['metacell_' + str(cell)] = np.sum(adata.obs['ReadsInTSS'][indices])
        meta_counter += 1

    # get most common cluster among neighbors and set it as meta cell cluster
    for n in neighbors:
        c = most_frequent(list(adata.obs['Clusters'][neighbors[n]]))
        neighbor_clusters[n] = c
    return neighbors, neighbor_clusters, lib_sizes

def group_atac_sarah(adata, n_neighbors=20, max_allowed_overlap=0.3, N_ITER=1000):
    mat = adata.obsp['distances'].todense()
    adata.obs['index'] = np.arange(len(adata.obs_names))

    # sample                                                                                                                                                          
    neighbors = {}
    neighbor_clusters = {}
    all_neighbors = set([])
    lib_sizes = {}
    counter = 0
    meta_counter = 0
    d_counter = {}
    # maintain check on which cell has already been selected                                                                                                          
    for obs_name in adata.obs['index']: d_counter[obs_name] = 0

    # max 20 iterations where check is applied for 90%                                                                                                                
    while counter < N_ITER:
        cell = np.random.choice(adata.obs['index'])
        if d_counter[cell] == 1:
            counter += 1
            continue

        # cell explored                                                                                                                                               
        d_counter[cell] = 1

        # get neighbors and append the current index                                                                                                                  
        indices = np.where(mat[cell, :] > 0)[1]
        indices = np.append(indices, cell)

        is_uniq = True
        # check overlap                                                                                                                                               
        overlap = len(all_neighbors.intersection(indices)) / n_neighbors
        if overlap >= max_allowed_overlap:
            is_uniq = False

        # IFF fewer than 30% of its neighbors are in set of all metacell neigbhors                                                                                            
        if not is_uniq:
            counter +=1
            continue

        # add to list, count from scratch  
        counter = 0
        all_neighbors.update(indices)
        neighbors["metacell_" + str(cell)] = indices
        lib_sizes['metacell_' + str(cell)] = np.sum(adata.obs['ReadsInTSS'][indices])
        meta_counter += 1

    # get most common cluster among neighbors and set it as meta cell cluster
    for n in neighbors:
        c = most_frequent(list(adata.obs['Clusters'][neighbors[n]]))
        neighbor_clusters[n] = c
    return neighbors, neighbor_clusters, lib_sizes

def meta_atac(groups, group_clusters, adata_atac, merge_mode='sum'):
    meta_cells = sorted(list(groups.keys()))
    n_meta_cells = len(meta_cells)
    X_atac = np.zeros((n_meta_cells, adata_atac.n_vars))
    clusters = [group_clusters[i] for i in meta_cells]
    adata_atac_X = adata_atac.X
    print(np.min(adata_atac_X), np.max(adata_atac_X))
    lib_sizes = np.zeros(X_atac.shape[0])
    for mc in range(n_meta_cells):
        if merge_mode == 'sum':
            X_atac[mc, :] = np.sum(adata_atac_X[groups[meta_cells[mc]], :], axis = 0)
        else:
            X_atac[mc, :] = np.mean(adata_atac_X[groups[meta_cells[mc]], :], axis = 0)
        lib_sizes[mc] = np.sum(adata_atac.obs['ReadsInTSS'][groups[meta_cells[mc]]])

    adata_atac_new = ad.AnnData(X_atac)
    adata_atac_new.obs_names = meta_cells
    adata_atac_new.var_names = adata_atac.var_names
    adata_atac_new.obs['Clusters'] = clusters
    adata_atac_new.obs['ReadsInTSS'] = lib_sizes

#     adata_atac_new.layers["raw"] = adata_atac_new.X.copy()
#     sc.pp.scale(adata_atac_new, max_value=10)
    sc.tl.pca(adata_atac_new, svd_solver='arpack', zero_center=False, use_highly_variable=False)
    sc.pp.neighbors(adata_atac_new, n_neighbors=20, n_pcs=30, metric = "cosine")
    sc.tl.umap(adata_atac_new, min_dist = 0.3)                                                                                                                        
    sc.pl.umap(adata_atac_new, edges = True, color = ['Clusters'])  
#     adata_atac_new.X = adata_atac_new.layers["raw"]

    return adata_atac_new

def ridge_reg(X, Y, if_scale=True, if_grid_search=True):

    # X: peak x TF matrix
    # Y: peak x meta cell (or single cell) matrix
    # W: learned coefficient vectors for separate regression of Y_i = X x W_i

    rhos = []
    mses = []
    W = []

    train_ix, test_ix = train_test_split(np.arange(Y.shape[1]), random_state = 9)
    tuned_parameters = {'alpha': [25000, 40000, 64000, 100000, 160000, 250000]}
    
    if if_scale:
        scaler = StandardScaler() 
        scaler.fit(X[train_ix])
        x_train = scaler.transform(X[train_ix]).astype(np.float32)
        x_test = scaler.transform(X[test_ix]).astype(np.float32)
    else:
        x_train = X[train_ix]
        x_test = X[test_ix]

    STEP = Y.shape[0]//10
    for i in range(Y.shape[0]):
        if i % STEP == 0:
            print(f'{10*i/STEP}%')
        y = Y[i]    
        y_train = y[train_ix]
        y_test = y[test_ix]
        if if_grid_search:
            clf = GridSearchCV(linear_model.Ridge(), tuned_parameters, cv=5, n_jobs=4, scoring="r2")
        else:
            clf = linear_model.Ridge(alpha=1000000)
        clf.fit(x_train, y_train)
        
        if if_grid_search:
            w = clf.best_estimator_.coef_
            print(clf.best_params_)
       	    print(clf.best_score_)
        else:
            w = clf.coef_
        W.append(list(w))
        ypred = clf.predict(x_test)
        rho, _ = stats.spearmanr(ypred, y_test) 
        mse = mean_squared_error(y_test, ypred)
        rhos.append(rho)
        mses.append(mse)

    return np.array(W), np.array(rhos), np.array(mses)


PARAMS = [[20, True, False, False, True, 'chipseeker', False, True, True],
          [90, True, False, False, True, 'chipseeker', False, True, True],
          [20, True, False, False, True, 'no filter', False, True, True],
          [20, False, False, False, True, 'chipseeker', False, True, True],
          [20, True, False, False, True, 'chipseeker', True, True, True],
          [20, True, True, True, True, 'chipseeker', False, True, True],
          [90, True, False, False, False, 'chipseeker', False, True, True],
          [90, True, False, False, False, 'chipseeker', True, True, True],
          [90, True, True, True, False, 'chipseeker', False, True, True]]


JOB_INDEX = int(sys.argv[1])-1
# original sneha 20: 
params_list = PARAMS[JOB_INDEX]
n_neighbors = params_list[0]
if_load_lsi = params_list[1]
if_poisson_correction = params_list[2]
if_pearson_residuals_first = params_list[3]
if_sneha = params_list[4]
if_filter_summits = params_list[5]
if_pearson_residuals_norm = params_list[6]
if_scale=params_list[7]
if_grid_search=params_list[8]

adata_atac, adata_motif = load_data(load_lsi=if_load_lsi, n_neighbors=n_neighbors)
# poisson correction
if if_poisson_correction:
    adata_atac = poisson_correction(adata_atac)
# normal I
if if_poisson_correction and if_pearson_residuals_first:
    sc.experimental.pp.normalize_pearson_residuals(adata_atac, theta=np.Inf)
elif if_pearson_residuals_first:
    sc.experimental.pp.normalize_pearson_residuals(adata_atac)
# meta-cell
if if_sneha:
    groups, group_clusters, lib_sizes = group_atac_sneha(adata_atac, n_neighbors=n_neighbors, 
                                                         max_allowed_overlap=0.9, N_ITER=20)
else:
    groups, group_clusters, lib_sizes = group_atac_sarah(adata_atac, n_neighbors=n_neighbors, 
                                                         max_allowed_overlap=0.3, N_ITER=1000)
if if_pearson_residuals_first:
    adata_atac_meta = meta_atac(groups, group_clusters, adata_atac, merge_mode='mean')
else:
    adata_atac_meta = meta_atac(groups, group_clusters, adata_atac, merge_mode='sum')
# filter
if if_filter_summits == 'chipseeker':
    adata_atac_filtered = adata_atac_meta[:, adata_atac.var['annotation'].isin(["Promoter (1-2kb)", "Promoter (<=1kb)", "Promoter (2-3kb)", "Distal Intergenic"])]
    adata_motif_filtered = adata_motif[adata_atac.var['annotation'].isin(["Promoter (1-2kb)", "Promoter (<=1kb)", "Promoter (2-3kb)", "Distal Intergenic"]), :]
elif if_filter_summits == 'accessibility':
    tmp = adata_atac_meta.X.copy()
    tmp[tmp>0] = 1
    accessibe_count = np.sum(tmp, axis=0)
    INDEX = np.argsort(accessibe_count)[::-1][:accessibe_count.shape[0]//2]
    adata_atac_filtered = adata_atac_meta[:, INDEX]
    adata_motif_filtered = adata_motif[INDEX, :]
else:
    adata_atac_filtered = adata_atac_meta
    adata_motif_filtered = adata_motif
# normal II
if if_pearson_residuals_first is False and if_pearson_residuals_norm:
    sc.experimental.pp.normalize_pearson_residuals(adata_atac_filtered)
    Y = adata_atac_filtered.X
    X = adata_motif_filtered.X.toarray()
elif if_pearson_residuals_first is False:
    Y = adata_atac_filtered.X / adata_atac_meta.obs['ReadsInTSS'].values[:, np.newaxis]
    X = adata_motif_filtered.X.toarray()
else:
    Y = adata_atac_filtered.X
    X = adata_motif_filtered.X.toarray()

train_ix, test_ix = train_test_split(np.arange(Y.shape[1]), random_state = 3)
X_train = X[train_ix]
Y_train = Y[:, train_ix]
X_test = X[test_ix]
Y_test = Y[:, test_ix]

w, r, mse = ridge_reg(X_train, Y_train, if_scale=if_scale, if_grid_search=if_grid_search)

ax = seaborn.displot(r)
plt.axvline(np.median(r), 0, 10000,  linestyle='--', color='r')
plt.savefig(f'model_selection_{JOB_INDEX}_correlation_hist.png', bbox_inches='tight')
plt.show()

# UMAP of learned TF weights for the meta cells
adata_w = sc.AnnData(w)
# sc.pp.scale(adata_w, max_value=10)
adata_w.obs = adata_atac_meta.obs   
adata_w.obs['rho'] = r
sc.tl.pca(adata_w, svd_solver='arpack', use_highly_variable = False)
sc.pp.neighbors(adata_w, n_neighbors=20, n_pcs=30, metric = "euclidean")
sc.tl.umap(adata_w, min_dist = 0.3)    
sc.tl.leiden(adata_w)


# plot using cluster annotations
sc.pl.umap(adata_w, edges = True, color = ['Clusters'])      

# plot Spearman correlation across the meta cells
g = seaborn.scatterplot(adata_w.obsm['X_umap'][:, 0], adata_w.obsm['X_umap'][:, 1],
                    hue=adata_w.obs['rho'], palette='Spectral_r', edgecolor='none',
                    alpha=0.5, legend='brief')
g.legend_.remove()

norm = plt.Normalize(adata_w.obs['rho'].min(), adata_w.obs['rho'].max())
sm = plt.cm.ScalarMappable(cmap="Spectral_r", norm=norm)
sm.set_array([])

# Remove the legend and add a colorbar
g.figure.colorbar(sm)

# sc.tl.leiden(adata_w)
sc.pl.umap(adata_w, edges = True, color = ['leiden'])      
# plt.savefig('foo.png', bbox_inches='tight')
clustering_scores = homogeneity_completeness_v_measure(adata_w.obs['Clusters'], adata_w.obs['leiden'].values)
print(w.shape)
print('{:.4f}-({:.4f})-{:.4f}'.format(np.min(r), np.median(r), np.max(r)))
print(PARAMS[JOB_INDEX])
print('{:.4f} & {:.4f} & {:.4f}'.format(*clustering_scores))
print(np.unique(adata_w.obs['Clusters']))
print(np.unique(adata_w.obs['leiden'].values))
