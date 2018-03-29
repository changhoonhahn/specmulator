'''

Test hod.hod.py


'''
import numpy as np 

from hodee import lhd_hod as LHDhod
from specmulator import util as UT 
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


def HOD_LHD(prior='sinha2017prior_narrow', samples=17):
    ''' test parameter values in the LHD and the test samples 
    '''
    hodprior = LHDhod.HODprior(prior) # HOD prior 
    hodrange = hodprior.range() 

    # parameter values of the LHD 
    theta_lhd = LHDhod.HOD_LHD(prior=prior, samples=samples, method='mdu', overwrite=False)
    assert theta_lhd.shape[0] == samples 
    assert theta_lhd.shape[1] == len(hodprior.labels)
    # parameter values of the test sample 
    theta_test = LHDhod.testHOD_LHD(prior=prior, samples=20, overwrite=False)
    assert theta_test.shape[1] == len(hodprior.labels)

    fig = plt.figure(figsize=(4*theta_lhd.shape[1], 4*theta_lhd.shape[1]))
    for i in range(theta_lhd.shape[1]): 
        for j in range(theta_lhd.shape[1]): 
            sub = fig.add_subplot(theta_lhd.shape[1], theta_lhd.shape[1], 1+j+theta_lhd.shape[1]*i) 
            if i != j: 
                sub.scatter(theta_lhd[:,j], theta_lhd[:,i], c='k') 
                sub.scatter(theta_test[:,j], theta_test[:,i], c='C1') 
            sub.set_xlim(np.array(hodrange)[:,j]) 
            sub.set_ylim(np.array(hodrange)[:,i]) 
            if i != theta_lhd.shape[1]-1: sub.set_xticklabels([]) 
            else: sub.set_xlabel(hodprior.labels[j], fontsize=20)
            if j != 0: sub.set_yticklabels([]) 
            else: sub.set_ylabel(hodprior.labels[i], fontsize=20)

    f = ''.join([UT.fig_dir(), 'tests/thetaHOD_lhd.png']) 
    fig.savefig(f, bbox_inches='tight') 
    plt.close() 
    return None 


if __name__=='__main__': 
    HOD_LHD() 
