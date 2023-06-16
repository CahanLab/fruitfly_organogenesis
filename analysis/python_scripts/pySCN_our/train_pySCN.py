
import scanpy as sc
import pandas as pd
import numpy as np 
import anndata
import pySingleCellNet as pySCN

from scipy import sparse, io

import pickle

sparsematrix = io.mmread('raw_train_exp.txt')
m_dense = sparsematrix.toarray()

var_names = np.genfromtxt('raw_train_rownames.txt', dtype=str)
col_names = np.genfromtxt('raw_train_colnames.txt', dtype=str)

# Export to txt:
exp_df = pd.DataFrame(m_dense, columns=col_names, index=var_names)

st_df = pd.read_csv("raw_train_st.csv", index_col=0)
exp_df = exp_df.loc[:, st_df.index]
train_adata = sc.AnnData(exp_df.T)
train_adata.obs = st_df
#train_adata.X = exp_df.T.values

selected_cell_types = train_adata.obs['manual_celltypes'].value_counts().index[train_adata.obs['manual_celltypes'].value_counts() > 2]
selected_cell_types = np.array(selected_cell_types)
selected_cell_types = selected_cell_types[selected_cell_types != 'Unknown']

train_adata = train_adata[train_adata.obs['manual_celltypes'].isin(selected_cell_types), :]
expTrain, expVal = pySCN.splitCommonAnnData(train_adata, ncells=200, dLevel="manual_celltypes")

pickle.dump(expTrain, open("expTrain_adata.pickle", "wb"))
pickle.dump(expVal, open("expVal_adata.pickle", "wb"))

expTrain = pickle.load(open("expTrain_adata.pickle", 'rb'))

[cgenesA, xpairs, tspRF] = pySCN.scn_train(expTrain, nTopGenes = 75, nRand = 100, nTrees = 1500 ,nTopGenePairs = 75, dLevel = "manual_celltypes", stratify=True, limitToHVG=True, normalization = True)
pickle.dump([cgenesA, xpairs, tspRF], open("SCN_classifier_obj.pickle", "wb"))

expVal = pickle.load(open("expVal_adata.pickle", "rb"))
expTrain = pickle.load(open("expTrain_adata.pickle", 'rb'))

[cgenesA, xpairs, tspRF] = pickle.load(open('SCN_classifier_obj.pickle', 'rb'))

adVal = pySCN.scn_classify(expVal, cgenesA, xpairs, tspRF, nrand = 0)

assessment =  pySCN.assess_comm(expTrain, adVal, resolution = 0.005, nRand = 0, dLevelSID = None, classTrain = "manual_celltypes", classQuery = "manual_celltypes")
pickle.dump(assessment, open("assessment.pickle", "wb"))

pySCN.plot_PRs(assessment)
ax = sc.pl.heatmap(adVal, adVal.var_names.values, groupby='SCN_class', cmap='viridis', dendrogram=False, swap_axes=True, figsize = [24, 20], save='heatmap_classification.png')

sub_expVal, other = pySCN.splitCommonAnnData(train_adata, ncells=70, dLevel="manual_celltypes")
adVal = pySCN.scn_classify(sub_expVal, cgenesA, xpairs, tspRF, nrand = 0)

assessment =  pySCN.assess_comm(expTrain, adVal, resolution = 0.005, nRand = 0, dLevelSID = None, classTrain = "manual_celltypes", classQuery = "manual_celltypes")
pickle.dump(assessment, open("sub_assessment.pickle", "wb"))

fig = pySCN.plot_PRs(assessment)
fig.savefig("PR.png")

ax = sc.pl.heatmap(adVal, adVal.var_names.values, groupby='SCN_class', cmap='viridis', dendrogram=False, swap_axes=True, figsize = [15, 15], save='sub_heatmap_classification.png')
