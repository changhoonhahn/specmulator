'''

code for compressing observables 


'''
import numpy as np
from sklearn.decomposition import PCA

import data as Dat


def PCA_fid(n_comp, mneut=0.0, nreal=1, nzbin=4, obvs='plk', poles=[0], Nmesh=360, rsd=True, 
        HODrange='sinha2017prior_narrow', krange=[0.01, 0.5]): 
    ''' Fiducial PCA that compresses the observable to n components 
    '''
    X = X_fid(nreal, nzbin, obvs=obvs, poles=poles, mneut=mneut, Nmesh=Nmesh, rsd=rsd, 
            HODrange=HODrange, krange=krange)
    mu_X = np.sum(X, axis=0)/np.float(X.shape[0])
    X -= mu_X 
    pca_fid = PCA(n_components=n_comp, whiten=False)
    pca_fid.fit(X)
    return pca_fid


def PCA_decomp(X, n_comp=None, exp_var=False):
    ''' PCA decomposition using sklearn's SVD 
    '''
    n_sample = X.shape[0]
    if n_comp is None: 
        n_comp = X.shape[1]

    mu_X = np.sum(X, axis=0)/np.float(n_sample)
    X -= mu_X 

    pca = PCA(n_components=n_comp, whiten=False)
    X_w = pca.fit_transform(X)
    # X_w = np.dot(X, V)
    V = pca.components_.T / np.sqrt(pca.explained_variance_)
    if not exp_var:
        return X_w, V
    else: 
        return X_w, V, pca.explained_variance_
