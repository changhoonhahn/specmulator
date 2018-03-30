'''



'''
import numpy as np 
# -- local --
import env
from specmulator import util as UT 
from specmulator import onlyhod as onlyHOD
# -- plotting --
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


def X_HODLHD():
    ''' Make sure that the powerspectra read in `specmulator.onlyhod.X_HODLHD` 
    are sensible. 
    '''
    # read P(k|theta_i) for theta_i in LHD
    k, X_pk = onlyHOD.X_HODLHD(1, 
            prior='sinha2017prior_narrow', samples=40, karr=True) 

    fig = plt.figure(figsize=(6,6))
    sub = fig.add_subplot(111)
    for i in range(X_pk.shape[0]):
        sub.plot(k, X_pk[i,:]) 
    sub.set_xlabel('$k$', fontsize=20) 
    sub.set_xscale('log') 
    sub.set_xlim([0.01, 0.5]) 
    sub.set_ylabel(r'$P_0(k|\theta_i^{LHD})$', fontsize=20) 
    sub.set_yscale('log') 
    f = ''.join([UT.fig_dir(), 'tests/test.X_HODLHD.png']) 
    fig.savefig(f, bbox_inches='tight') 
    plt.close() 
    return None 


def X_testHODLHD(): 
    ''' Make sure that the powerspectra `specmulator.onlyhod.X_testHODLHD` 
    are sensible. 
    '''
    # read P(k|theta_i) for theta_i in LHD
    k, X_pk = onlyHOD.X_HODLHD(1,  
            prior='sinha2017prior_narrow', samples=40, karr=True) 

    # read P(k|theta_i) for theta_i in LHD test sample 
    X_pk_test = onlyHOD.X_testHODLHD(1, 
            prior='sinha2017prior_narrow', samples=20) 

    fig = plt.figure(figsize=(6,6))
    sub = fig.add_subplot(111)
    for i in range(X_pk.shape[0]):
        sub.plot(k, X_pk[i,:], c='k') 

    for i in range(X_pk_test.shape[0]):
        sub.plot(k, X_pk_test[i,:]) 
    sub.set_xlabel('$k$', fontsize=20) 
    sub.set_xscale('log') 
    sub.set_xlim([0.01, 0.5]) 
    sub.set_ylabel(r'$P_0(k|\theta_i^{LHD})$', fontsize=20) 
    sub.set_yscale('log') 
    f = ''.join([UT.fig_dir(), 'tests/test.X_testHODLHD.png']) 
    fig.savefig(f, bbox_inches='tight') 
    plt.close() 
    return None 


if __name__=="__main__": 
    X_HODLHD()
    X_testHODLHD()
