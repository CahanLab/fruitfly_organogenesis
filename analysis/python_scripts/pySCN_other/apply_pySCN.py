import sys 
import os 
from scipy import sparse, io
import pandas as pd 
import numpy as np 
import scanpy as sc
import pickle
import pySingleCellNet as pySCN

[cgenesA, xpairs, tspRF] = pickle.load(open("SCN_classifier_obj.pickle", "rb"))

sparsematrix = io.mmread('raw_query_exp.txt')
m_dense = sparsematrix.toarray()

var_names = np.genfromtxt('raw_query_rownames.txt', dtype=str)
col_names = np.genfromtxt('raw_query_colnames.txt', dtype=str)

exp_df = pd.DataFrame(m_dense, columns=col_names, index=var_names)
query_adata = sc.AnnData(exp_df.T)

adVal = pySCN.scn_classify(query_adata, cgenesA, xpairs, tspRF, nrand = 0)
class_matrix = adVal.to_df()
SCN_label_mat = adVal.obs
result = pd.concat([SCN_label_mat, class_matrix], axis=1)
result.to_csv("SCN_classification.csv")
