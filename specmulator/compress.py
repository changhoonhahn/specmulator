'''

code for compressing observables 


'''
import numpy as np
from sklearn.decomposition import PCA

import data as Dat


class Svd(object): 
    ''' Class for Singular Value Decomposition because the sklearn.decomposition.PCA 
    is gnarly. I just use numpy.linalg.svd
    '''
    def __init__(self, n_comp=None): 
        '''
        '''
        self.n_comp = n_comp

    def fit(self, X): 
        '''
        '''
        self._mean = np.mean(X, axis=0) 
        self._stddev = np.std(X, axis=0) 

        X_w = self._white(X) 

        # X_w  = U * np.diag(sigma) * Vt 
        U, sigma, Vt = np.linalg.svd(X_w) 

        self._U = U 
        self._sigma = sigma
        self._Vt = Vt
        
        if self.n_comp is None: 
            self.n_comp = U.shape[0] 

        exp_var = ((sigma**2) / (X_w.shape[0] - 1))
        exp_var_ratio = exp_var/exp_var.sum()
        
        self.exp_var = exp_var[:self.n_comp]
        self.exp_var_ratio = exp_var_ratio[:self.n_comp]

        return U, sigma, Vt

    def transform(self, X): 
        '''
        '''
        X_w = self._white(X) 

        return np.dot(X_w, self._Vt[:self.n_comp,:].T/np.sqrt(len(self._sigma)))

    def inv_transform(self, pc):  
        '''
        '''
        # check dimensions
        assert pc.shape[1] == self.n_comp

        X_rec = np.dot(pc, self._Vt[:self.n_comp,:]) * np.sqrt(len(self._sigma)) 

        return (X_rec * self._stddev) + self._mean
    
    def _white(self, X): 
        '''
        '''
        return (X - self._mean)/self._stddev


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
