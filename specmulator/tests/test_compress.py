'''

'''
import numpy as np 

import env 
import compress as Comp 


def PCA_components(): 
    '''
    '''
    X = Comp.X_fid(1, 4, krange=[0.01, 0.5])
    _, _, exp_var = Comp.PCA_decomp(X, exp_var=True) 
    f_exp_var = exp_var/exp_var.sum()
    print(f_exp_var)
    n_thresh = np.sum(f_exp_var > 0.001)
    print(str(n_thresh)+' components have greater than 0.1%')
    return None


if __name__=="__main__": 
    PCA_components()
