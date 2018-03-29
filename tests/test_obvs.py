import numpy as np 

import env
import util as UT
import obvs as Obvs

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


def Plk_halo_mneut_ratio(nzbin=4, zspace=False): 
    ''' Plot the ratio of P_l^mneut(k)/P_l^0.0eV 
    for different neutrino masses 
    '''
    mneuts = [0.0, 0.06, 0.10, 0.15, 0.6] # eV
    p0ks_mneut, p2ks_mneut, p4ks_mneut = [], [], [] 
    for mneut in mneuts: 
        p0ks, p2ks, p4ks = [], [], [] 
        for ireal in range(1, 101): 
            # read all 100 realizations
            plk_i = Obvs.Plk_halo(mneut, ireal, nzbin, zspace=zspace) 
            if ireal == 1: k = plk_i['k']
            p0ks.append(plk_i['p0k']) 
            p2ks.append(plk_i['p2k']) 
            p4ks.append(plk_i['p4k']) 
        # plot the average 
        p0ks_mneut.append(np.average(np.array(p0ks), axis=0)) 
        p2ks_mneut.append(np.average(np.array(p2ks), axis=0)) 
        p4ks_mneut.append(np.average(np.array(p4ks), axis=0)) 

    plks_mneut = [p0ks_mneut, p2ks_mneut, p4ks_mneut]
    fig = plt.figure(figsize=(15, 5)) 
    for i, ell in enumerate([0,2,4]):  
        sub = fig.add_subplot(1,3,i+1)
        for ii in range(len(mneuts)): 
            sub.plot(k, plks_mneut[i][ii]/plks_mneut[i][0], lw=2, label=r'$\sum m_\nu = $ '+str(mneuts[ii])+'eV') 
        if i == 0: 
            sub.legend(loc='lower right', prop={'size': 12}) 
        else: 
            sub.set_yticks([])
        sub.set_xscale('log')
        sub.set_xlim([0.01, 0.5]) 
        sub.set_xlabel('k', fontsize=20) 

        sub.set_ylim([0.9, 1.15]) 
        sub.set_ylabel('$P_{'+str(ell)+'}(k)/P_{'+str(ell)+'}^{0.0\mathrm{eV}}(k)$', fontsize=20) 

    if zspace: str_space = 'z' 
    else: str_space = 'r'
    fig.savefig(''.join([UT.fig_dir(), 'tests/plk_halo.mneuts_ratio.nzbin', str(nzbin), 
        '.', str_space, 'space.png']), bbox_inches='tight') 
    return None 


def Plk_halo_mneut(nzbin=4, zspace=False): 
    ''' Plot P_l(k) for different neutrino masses 
    '''
    mneuts = [0.0, 0.06, 0.10, 0.15, 0.6] # eV
    p0ks_mneut, p2ks_mneut, p4ks_mneut = [], [], [] 
    for mneut in mneuts: 
        p0ks, p2ks, p4ks = [], [], [] 
        for ireal in range(1, 101): 
            # read all 100 realizations
            plk_i = Obvs.Plk_halo(mneut, ireal, nzbin, zspace=zspace) 
            if ireal == 1: k = plk_i['k']
            p0ks.append(plk_i['p0k']) 
            p2ks.append(plk_i['p2k']) 
            p4ks.append(plk_i['p4k']) 
        # plot the average 
        p0ks_mneut.append(np.average(np.array(p0ks), axis=0)) 
        p2ks_mneut.append(np.average(np.array(p2ks), axis=0)) 
        p4ks_mneut.append(np.average(np.array(p4ks), axis=0)) 

    plks_mneut = [p0ks_mneut, p2ks_mneut, p4ks_mneut]
    fig = plt.figure(figsize=(15, 5)) 
    for i, ell in enumerate([0,2,4]):  
        sub = fig.add_subplot(1,3,i+1)
        for mneut, plk in zip(mneuts, plks_mneut[i]):
            sub.plot(k, plk, lw=2, label=r'$\sum m_\nu = $ '+str(mneut)+'eV') 
        if i == 0: 
            sub.legend(loc='lower right', prop={'size': 12}) 
        else: 
            sub.set_yticks([])
        sub.set_xscale('log')
        sub.set_xlim([0.01, 0.15]) 
        sub.set_xlabel('k', fontsize=20) 

        sub.set_ylim([1e3, 1e5]) 
        sub.set_yscale('log') 
        sub.set_ylabel('$k P_{'+str(ell)+'}(k)$', fontsize=20) 

    if zspace: str_space = 'z' 
    else: str_space = 'r'
    fig.savefig(''.join([UT.fig_dir(), 'tests/plk_halo.mneuts.nzbin', str(nzbin), 
        '.', str_space, 'space.png']), bbox_inches='tight') 
    return None 


def Plk_halo(mneut=0.0, nzbin=4, zspace=False): 
    ''' **TESTED --- Nov 7, 2017 ** 
    Test the Plk_halo 
    '''
    p0ks, p2ks, p4ks = [], [], [] 
    for ireal in range(1, 101): 
        # read all 100 realizations
        plk_i = Obvs.Plk_halo(mneut, ireal, nzbin, zspace=zspace) 
        if ireal == 1: k = plk_i['k']
        p0ks.append(plk_i['p0k']) 
        p2ks.append(plk_i['p2k']) 
        p4ks.append(plk_i['p4k']) 
    
    fig = plt.figure() 
    sub = fig.add_subplot(111)
    for p0k, p2k, p4k in zip(p0ks, p2ks, p4ks): 
        sub.plot(k, k * p0k, c='k', lw=0.1) 
        sub.plot(k, k * p2k, c='b', lw=0.1) 
        sub.plot(k, k * p4k, c='r', lw=0.1) 
    # plot the average 
    sub.plot(k, k * np.average(np.array(p0ks), axis=0), c='k', lw=2, ls='--', label='$\ell=0$') 
    sub.plot(k, k * np.average(np.array(p2ks), axis=0), c='b', lw=2, ls='--', label='$\ell=2$') 
    sub.plot(k, k * np.average(np.array(p4ks), axis=0), c='r', lw=2, ls='--', label='$\ell=4$') 

    sub.set_xlim([0.01, 0.15]) 
    sub.set_xlabel('k', fontsize=20) 

    sub.set_ylim([-2000., 2500.]) 
    sub.set_ylabel('$k P(k)$', fontsize=20) 
    sub.legend(loc='lower right', prop={'size': 15}) 
    if zspace: str_space = 'z' 
    else: str_space = 'r'
    fig.savefig(''.join([UT.fig_dir(), 'tests/plk_halo.', str(mneut), 'eV.nzbin', str(nzbin), 
        '.', str_space, 'space.png']), bbox_inches='tight') 
    return None 


if __name__=="__main__": 
    Plk_halo_mneut_ratio(nzbin=4, zspace=False)
    Plk_halo_mneut_ratio(nzbin=4, zspace=True)
    #Plk_halo_mneut(nzbin=4, zspace=False)
    #Plk_halo_mneut(nzbin=4, zspace=True)
    #Plk_halo(mneut=0.6, zspace=False)
    #Plk_halo(mneut=0.6, zspace=True)
