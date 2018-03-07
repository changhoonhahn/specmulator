'''

test data.py

'''
import numpy as np 

import env 
import util as UT 
import data as Data 

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


def Plk_fiducial(nreal): 
    ''' Test that the methods Data.X_fid and Data.Fiducial_Obvs
    are returning sensible results for obvs = 'plk'
    '''
    karr, Xfid = Data.X_fid(nreal, 4, obvs='plk', karr=True)
    muX = np.sum(Xfid, axis=0)/float(Xfid.shape[0])  
    sigX = (np.diag(np.cov(Xfid.T)))

    fig = plt.figure()
    sub = fig.add_subplot(111)
    sub.errorbar(karr, muX, yerr=sigX, fmt='.k') 
    for i in range(Xfid.shape[0]): 
        sub.plot(karr, Xfid[i,:], c='k', lw=1) 

    sub.set_xlim([1e-2, 0.5]) 
    sub.set_xlabel('$k$', fontsize=20) 
    sub.set_xscale('log') 
    sub.set_ylim([2e3, 2e6]) 
    sub.set_ylabel('$P(k)$', fontsize=20) 
    sub.set_yscale('log') 

    f = ''.join([UT.fig_dir(), 'tests/test_data.Plk_fiducial.png']) 
    fig.savefig(f, bbox_inches='tight') 
    return None


def runNeutParticles(nreal):
    ''' Test that the method NeutrinoHalos works properly 
    '''
    cat = Data.NeutParticles(0.0, nreal, 4, clobber=True)  # mneut = 0.0eV, realization #1, z = 0 
    print type(cat)
    return None 


def runNeutHalos(nreal):
    ''' Test that the method NeutrinoHalos works properly 
    '''
    halocat = Data.NeutHalos(0.0, nreal, 4, clobber=True)  # mneut = 0.0eV, realization #1, z = 0 
    return None 


if __name__=="__main__": 
    Plk_fiducial(1)
    pass
