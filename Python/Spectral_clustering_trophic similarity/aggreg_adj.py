import numpy as np
import scipy
import pandas as pd
from sklearn.metrics import jaccard_similarity_score
from sklearn.cluster import SpectralClustering

### Read adjacency matrix
clust=pd.read_csv("../../data/classes.csv",sep=";",decimal=".",index_col=0)
adj_mat=pd.read_csv("../../data/adjacence_esp.csv",sep=";",decimal=".",index_col=0)

### Aggregate 
