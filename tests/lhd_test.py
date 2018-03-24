'''
'''
import numpy as np 

from specmulator import lhd
from specmulator import util as UT 

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


def theta_test(HODrange='sinha2017prior_narrow',): 
    ''' Test the theta_test parameters
    '''
    if HODrange == 'sinha2017prior_narrow': # narrow alpha range 
        range_descrip = "Sinha et al (2017) prior with narrow alpha range"
        HOD_range_min = [11., 0.001, 6., 12., 0.5]
        HOD_range_max = [12.2, 1., 14., 14., 1.5]
        theta_lbl = ['log $M_\mathrm{min}$', '$\sigma_{\mathrm{log}\,M}$', 'log $M_0$', 'log $M_1$', r'$\alpha$']
    elif HODrange == 'sinha2017prior':  
        range_descrip = "Sinha et al (2017) prior"
        HOD_range_min = [11., 0.001, 6., 12., 0.001]
        HOD_range_max = [12.2, 1., 14., 14., 2.]
        theta_lbl = ['log $M_\mathrm{min}$', '$\sigma_{\mathrm{log}\,M}$', 'log $M_0$', 'log $M_1$', r'$\alpha$']
    else: 
        raise NotImplementedError
    # LHD thetas
    theta_lhd = lhd.HOD_LHD(HODrange=HODrange, samples=40, method='mdu')
    
    # test thetas
    theta_test = lhd.HOD_LHD_test(HODrange=HODrange, samples=20, seed=1)

    fig = plt.figure(figsize=(4*theta_test.shape[1], 4*theta_test.shape[1]))
    for i in range(theta_test.shape[1]): 
        for j in range(theta_test.shape[1]): 
            sub = fig.add_subplot(theta_test.shape[1], theta_test.shape[1], 1+j+theta_test.shape[1]*i) 
            if i != j: 
                sub.scatter(theta_lhd[:,j], theta_lhd[:,i], c='k') 
                sub.scatter(theta_test[:,j], theta_test[:,i], c='C1') 
            sub.set_xlim([HOD_range_min[j], HOD_range_max[j]]) 
            sub.set_ylim([HOD_range_min[i], HOD_range_max[i]]) 
            if i != theta_test.shape[1]-1: sub.set_xticklabels([]) 
            else: sub.set_xlabel(theta_lbl[j], fontsize=20)
            if j != 0: sub.set_yticklabels([]) 
            else: sub.set_ylabel(theta_lbl[i], fontsize=20)

    f = ''.join([UT.fig_dir(), 'tests/theta_test.png']) 
    fig.savefig(f, bbox_inches='tight') 
    plt.close() 
    return None 


def thetaLHD():  
    ''' ***TESTED***
    Test thetaLHD using the HOD parameter priors 
    '''
    # from Sinha et al (2017) prior 
    theta_lbl = ['log $M_\mathrm{min}$', '$\sigma_{\mathrm{log}\,M}$', 'log $M_0$', 'log $M_1$', r'$\alpha$']
    HOD_range_min = [11., 0.001, 6., 12., 0.001]
    HOD_range_max = [12.2, 1., 14., 14., 2.]
    
    samples = 17
    m_list = ['maximin', 'centermaximin', 'mdu', 'nohl']
    for meth in m_list: 
        lhcube = lhd.thetaLHD([HOD_range_min, HOD_range_max], samples=samples, method=meth)

        fig = plt.figure(figsize=(9,9))
        for i in range(lhcube.shape[1]):
            for j in range(lhcube.shape[1]): 
                if i < j:
                    sub = fig.add_subplot(lhcube.shape[1],lhcube.shape[1],lhcube.shape[1]*j+i+1)
                    sub.scatter(lhcube[:,i], lhcube[:,j])
                    sub.set_xlim([HOD_range_min[i], HOD_range_max[i]])
                    if j == lhcube.shape[1]-1: sub.set_xlabel(theta_lbl[i], fontsize=20)
                    else: sub.set_xticks([])
                    sub.set_ylim([HOD_range_min[j], HOD_range_max[j]])
                    if i == 0: sub.set_ylabel(theta_lbl[j], fontsize=20)
                    else: sub.set_yticks([])            
                elif i == j:
                    sub = fig.add_subplot(lhcube.shape[1],lhcube.shape[1],lhcube.shape[1]*j+i+1)
                    sub.hist(lhcube[:,i], range=[HOD_range_min[i], HOD_range_max[i]], normed=True)
                    sub.set_xlim([HOD_range_min[i], HOD_range_max[i]])
                    if i != 0: sub.set_yticks([])            
                    else: sub.set_ylabel(theta_lbl[j], fontsize=20)
                    if i != lhcube.shape[1]-1: sub.set_xticks([])            
                    else: sub.set_xlabel(theta_lbl[i], fontsize=20)
        f = ''.join([UT.fig_dir(), 'tests/test.thetaLHD.HOD.', meth, '.', str(samples), '.samples.png']) 
        fig.savefig(f, bbox_inches='tight') 
    return None


def LHD(): 
    ''' *** TESTED *** 
    Testing the different methods of LHD 
    '''
    dim = 5 
    samples = 17
    m_list = ['maximin', 'centermaximin', 'mdu', 'nohl']

    for meth in m_list: 
        lhcube = lhd.LHD(dim, samples=samples, method=meth, niter=1000)

        fig = plt.figure(figsize=(9,9))
        for i in range(lhcube.shape[1]):
            for j in range(lhcube.shape[1]): 
                if i < j:
                    sub = fig.add_subplot(lhcube.shape[1],lhcube.shape[1],lhcube.shape[1]*j+i+1)
                    sub.scatter(lhcube[:,i], lhcube[:,j])
                    sub.set_xlim([0., 1.])
                    if j == lhcube.shape[1]-1: 
                        sub.set_xticks([0., 0.5, 1.])
                        sub.set_xlabel(r'$\theta_'+str(i)+'$', fontsize=20)
                    else: 
                        sub.set_xticks([])
                    sub.set_ylim([0., 1.])
                    if i == 0: 
                        sub.set_yticks([0., 0.5, 1.])            
                        sub.set_ylabel(r'$\theta_'+str(j)+'$', fontsize=20)
                    else: 
                        sub.set_yticks([])            
                elif i == j:
                    sub = fig.add_subplot(lhcube.shape[1],lhcube.shape[1],lhcube.shape[1]*j+i+1)
                    sub.hist(lhcube[:,i], range=[0., 1], normed=True)
                    sub.set_xlim([0., 1.])
                    if i != 0: 
                        sub.set_yticks([])            
                    if i != lhcube.shape[1]-1: 
                        sub.set_xticks([])            
        f = ''.join([UT.fig_dir(), 'tests/test.LHD.', meth, '.', str(dim), 'dim.', str(samples), '.samples.png']) 
        fig.savefig(f, bbox_inches='tight') 
    return None


if __name__=="__main__": 
    theta_test()
