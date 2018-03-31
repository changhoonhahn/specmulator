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


def HODLHD(test=False, ndim=None): 
    if not test: 
        onlyHOD.HOD_LHD(prior='sinha2017prior_narrow', method='mdu', samples=40, ndim=ndim) 
    else: 
        onlyHOD.testHOD_LHD(prior='sinha2017prior_narrow', samples=20, ndim=ndim) 
    return None 


def check_X_HODLHD(): 
    ''' run specmulator.onlyhod._check_X_HODLHD, which checks which pk files 
    are missing in the LHD
    '''
    for i in range(1,11): 
        onlyHOD._check_X_HODLHD(i, 
                prior='sinha2017prior_narrow', samples=40, karr=True) 
    return None 


def X_HODLHD():
    ''' Make sure that the powerspectra read in `specmulator.onlyhod.X_HODLHD` 
    are sensible. 
    '''
    # read in theta_i's of the LHD 
    lhd = onlyHOD.HOD_LHD(prior='sinha2017prior_narrow', method='mdu', samples=40) 

    # read P(k|theta_i) for theta_i in LHD average over
    # 10 realizations to make less noisy 
    X_pk = []
    for i in range(1,11): 
        k, X_pk_i = onlyHOD.X_HODLHD(i, 
                prior='sinha2017prior_narrow', samples=40, method='mdu', karr=True) 
        X_pk.append(X_pk_i) 
    X_pk = np.mean(X_pk, axis=0) 

    fig = plt.figure(figsize=(12,6))
    sub1 = fig.add_subplot(121) # 2d projections of parameter space 
    sub2 = fig.add_subplot(122)
    for i in range(X_pk.shape[0]):
        sub1.scatter([lhd[i,3]], [lhd[i,4]]) 
        sub2.plot(k, k*X_pk[i,:]) 
    sub1.set_xlabel('log $M_1$', fontsize=20) 
    sub1.set_xlim([12., 14.])
    sub1.set_ylabel(r'$\alpha$', fontsize=20) 
    sub1.set_ylim([0.5, 1.5])

    sub2.set_xlabel('$k$', fontsize=20) 
    sub2.set_xscale('log') 
    sub2.set_xlim([0.01, 0.5]) 
    sub2.set_ylabel(r'$k\, P_0(k|\theta_i^{LHD})$', fontsize=20) 
    sub2.set_yscale('log') 
    f = ''.join([UT.fig_dir(), 'tests/test.X_HODLHD.png']) 
    fig.savefig(f, bbox_inches='tight') 
    plt.close() 
    return None 


def X_testHODLHD(): 
    ''' Make sure that the powerspectra `specmulator.onlyhod.X_testHODLHD` 
    are sensible. 
    '''
    # read in theta_i's of the LHD 
    lhd = onlyHOD.HOD_LHD(prior='sinha2017prior_narrow', method='mdu', samples=40) 
    lhd_test = onlyHOD.testHOD_LHD(prior='sinha2017prior_narrow', samples=20) 
    # read P(k|theta_i) for theta_i in LHD average over
    # 10 realizations to make less noisy 
    X_pk = []
    for i in range(1,11): 
        k, X_pk_i = onlyHOD.X_HODLHD(i, 
                prior='sinha2017prior_narrow', samples=40, karr=True) 
        X_pk.append(X_pk_i) 
    X_pk = np.mean(X_pk, axis=0) 

    # read P(k|theta_i) for theta_i in LHD test sample 
    # average over 10 realization to make less noisy
    X_pk_test = [] 
    for i in range(1, 11): 
        X_pk_test_i = onlyHOD.X_testHODLHD(i, 
                prior='sinha2017prior_narrow', samples=20) 
        X_pk_test.append(X_pk_test_i) 
    X_pk_test = np.mean(X_pk_test, axis=0) 

    fig = plt.figure(figsize=(12,6))
    sub1 = fig.add_subplot(121)
    sub2 = fig.add_subplot(122)
    sub1.scatter(lhd[:,3], lhd[:,4], c='k') 
    for i in range(X_pk.shape[0]):
        sub2.plot(k, k*X_pk[i,:], c='k', ls='--', lw=1) 

    for i in range(X_pk_test.shape[0]):
        sub1.scatter([lhd_test[i,3]], [lhd_test[i,4]])
        sub2.plot(k, k*X_pk_test[i,:]) 
    sub1.set_xlabel('log $M_1$', fontsize=20) 
    sub1.set_xlim([12., 14.])
    sub1.set_ylabel(r'$\alpha$', fontsize=20) 
    sub1.set_ylim([0.5, 1.5])

    sub2.set_xlabel('$k$', fontsize=20) 
    sub2.set_xscale('log') 
    sub2.set_xlim([0.01, 0.5]) 
    sub2.set_ylabel(r'$k P_0(k|\theta_i^{LHD})$', fontsize=20) 
    sub2.set_yscale('log') 
    f = ''.join([UT.fig_dir(), 'tests/test.X_testHODLHD.png']) 
    fig.savefig(f, bbox_inches='tight') 
    plt.close() 
    return None 


if __name__=="__main__": 
    HODLHD(test=False, ndim=1)
