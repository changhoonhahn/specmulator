'''

Tests for forwardmodel.py

'''
import env
import util as UT
import data as Dat
import forwardmodel as FM

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


def Galaxies():
    ''' Test that the Galaxies function correctly populats halo catalogs with an HOD 
    '''
    # HOD parameters (Zheng+2007 Mr < -21)  
    p_hod = {'logM0': 11.92, 'sigma_logM': 0.39, 'logMmin': 12.79, 'alpha': 1.15, 'logM1': 13.94}

    # import Neutrino halo with 0.0eV realization 1 at z=0
    halos = Dat.NeutHalos(0.0, 1, 4) 
    halos['RSDPosition'] = FM.RSD(halos) 
    # measure the P_l(k) of the galaxy catalog
    plk_halo = FM.Observables(halos, observable='plk', rsd=True, Nmesh=256)

    # now run HOD on the halo catalog
    gals = FM.Galaxies(halos, p_hod, seed=10)  
    gals['RSDPosition'] = FM.RSD(gals) 
    
    # measure the P_l(k) of the galaxy catalog
    plk = FM.Observables(gals, observable='plk', rsd=True, Nmesh=256)
    
    fig = plt.figure() 
    sub = fig.add_subplot(121)
    for ell in [0,2,4]:
        sub.plot(plk['k'], plk['k'] * plk['p'+str(ell)+'k']) 
        sub.plot(plk_halo['k'], plk_halo['k'] * plk_halo['p'+str(ell)+'k'], ls='--') 
    sub.set_xlabel('$k$', fontsize=20) 
    sub.set_xlim([0.01, 0.5]) 
    sub.set_xscale('log') 
    sub.set_ylabel('$k P_\ell(k)$', fontsize=20) 
    sub = fig.add_subplot(122)
    sub.plot(plk['k'], plk['p0k']) 
    sub.plot(plk_halo['k'], plk_halo['p0k'], ls='--') 
    sub.set_xlabel('$k$', fontsize=20) 
    sub.set_xlim([0.01, 0.5]) 
    sub.set_xscale('log') 
    sub.set_ylabel('$P_\ell(k)$', fontsize=20) 
    sub.set_yscale('log') 
    fig.savefig(''.join([UT.fig_dir(), 'tests/hodgalaxy_pk.png']), bbox_inches='tight')  
    return None  


if __name__=="__main__": 
    Galaxies() 
