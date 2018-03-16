'''

'''
import numpy as np 

import env 
import util as UT 
import data as Dat
import compress as Comp 

import matplotlib as mpl 
import matplotlib.pyplot as plt
mpl.rcParams['text.usetex'] = True
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['axes.linewidth'] = 1.5
mpl.rcParams['axes.xmargin'] = 1
mpl.rcParams['xtick.labelsize'] = 'x-large'
mpl.rcParams['xtick.major.size'] = 5
mpl.rcParams['xtick.major.width'] = 1.5
mpl.rcParams['ytick.labelsize'] = 'x-large'
mpl.rcParams['ytick.major.size'] = 5
mpl.rcParams['ytick.major.width'] = 1.5
mpl.rcParams['legend.frameon'] = False


def test_Svd(): 
    '''
    '''
    k_arr, pk_lhd = Dat.X_lhd(0.0, 1, 4, 1, obvs='pk', 
            ell=0, Nmesh=360, rsd=True, krange=[0.01, 0.5], karr=True,
            HODrange='sinha2017prior_narrow', samples=40, method='mdu', silent=True)
    svd = Comp.Svd()
    _ = svd.fit(pk_lhd) 
    pcs = svd.transform(pk_lhd) 
    assert np.allclose(svd.inv_transform(pcs), pk_lhd) 

    fig = plt.figure()
    sub = fig.add_subplot(111)

    sub.plot(range(len(svd.exp_var_ratio)+1), np.concatenate([np.array([0]), np.cumsum(svd.exp_var_ratio)]))
    sub.set_ylim([5e-1, 1.]) 
    sub.set_yscale('log') 

    f = ''.join([UT.fig_dir(), 'tests/plk_svd.exp_var.png']) 
    fig.savefig(f, bbox_inches='tight') 
    return None


def test_PC(): 
    '''
    '''
    k_arr, pk_lhd = Dat.X_lhd(0.0, 1, 4, 1, obvs='pk', 
            ell=0, Nmesh=360, rsd=True, krange=[0.01, 0.5], karr=True,
            HODrange='sinha2017prior_narrow', samples=40, method='mdu', silent=True)
    svd = Comp.Svd()
    _ = svd.fit(pk_lhd) 
    pcs = svd.transform(pk_lhd) 

    fig = plt.figure()
    sub = fig.add_subplot(111)
    for i in range(pcs.shape[0]):  
        sub.plot(range(pcs.shape[1]), pcs[i,:])
    f = ''.join([UT.fig_dir(), 'tests/plk_svd.pc.png']) 
    fig.savefig(f, bbox_inches='tight') 
    return None




if __name__=="__main__": 
    test_PC()
