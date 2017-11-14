'''

Code for generating galaxy catalogs from halo catalog  
for the HOD parameters sampled by a LHD 


'''
import os 
import sys as Sys
import numpy as np 

import env 
import lhd
import util as UT 
import data as Dat
import forwardmodel as FM 


def HOD_LHD(HODrange='sinha2017prior_narrow', samples=17): 
    ''' Build HOD_LHD for all the method options 
    ''' 
    m_list = ['maximin', 'centermaximin', 'mdu', 'nohl']
    for meth in m_list: 
        print(meth)
        lhd.HOD_LHD(HODrange=HODrange, samples=samples, method=meth) 
    return None 


def hodlhd_catalogs(mneut, nreal, nzbin, seed_hod, i_p, HODrange='sinha2017prior_narrow', method='mdu', samples=17): 
    ''' Generate HOD catalogs from specified halo catalog 
    based on the LHD sampled by the Sinha M., et al. (2017) HOD parameter priors. 

    parameters
    ----------
    mneut : float, 
        total neutrino mass 

    nreal : int,
        realization number 

    nzbin : int, 
        integer specifying the redshift of the snapshot. 
        nzbin = 0 --> z=3
        nzbin = 1 --> z=2
        nzbin = 2 --> z=1
        nzbin = 3 --> z=0.5
        nzbin = 4 --> z=0
    
    seed_hod : int, 
        random seed for the HOD 

    HODrange : str, optional
        string specifying the HOD range. Default is 'sinha2017prior', which uses
        the prior from Sinha et al. (2017) 

    method : str 
        string specifying the method of LHD. Default is nohl 

    samples : int, optional 
        sample size of the LHD. Default is 17, which is from the fixed sample size of NOHL method.

    LOS : list, optional 
        3 element list of integers 0 or 1 that specify the the line-of-sight direction 
    '''
    _ = lhd.HODLHD_NeutCatalog(mneut, nreal, nzbin, seed_hod, i_p, 
            HODrange=HODrange, method=method, samples=samples)
    return None 


def hodlhd_observables(obvs, mneut, nreal, nzbin, seed_hod, i_p, HODrange='sinha2017prior_narrow', method='nohl', samples=17, 
        Nmesh=360, rsd=True): 
    ''' Calculate and save observables of the HOD LHD catalogs
    '''
    _ = lhd.HODLHD_NeutObvs(obvs, mneut, nreal, nzbin, seed_hod, i_p,
            HODrange=HODrange, method=method, samples=samples, Nmesh=Nmesh, rsd=rsd)
    return None 


if __name__=="__main__": 
    mod = Sys.argv[1]
    if mod == 'lhd':
        # e.g. python run/run_hodlhd_catalog.py lhd 17 
        nsample = Sys.argv[2]
        HOD_LHD(samples=nsample)
    elif mod in ['catalog', 'observable']: 
        # e.g. python run/run_hodlhd_catalog.py catalog 0.0 1 4 1 1 
        # e.g. python run/run_hodlhd_catalog.py observation 0.0 1 4 1 1 plk real 
        mneut = float(Sys.argv[2])
        nreal = int(Sys.argv[3]) 
        nzbin = int(Sys.argv[4]) 
        seed_hod = int(Sys.argv[5]) 
        i_p = int(Sys.argv[6]) 

        print('%f eV'%(mneut))
        print('realization %i'%(nreal))
        print('nzbin %i'%(nzbin))
        print('random seed %i ' % seed_hod)
        print('%i th HOD parameter' % i_p)
    
        if mod == 'catalog': 
            hodlhd_catalogs(mneut, nreal, nzbin, seed_hod, i_p, HODrange='sinha2017prior_narrow', method='mdu', samples=17)
        elif mod == 'observable': 
            obvs = Sys.argv[7] 
            space = Sys.argv[8]
            if space == 'real': 
                rsd_bool = False 
            elif space == 'z': 
                rsd_bool = True 

            hodlhd_observables(obvs, mneut, nreal, nzbin, seed_hod, i_p, 
                    HODrange='sinha2017prior_narrow', method='mdu', samples=17, Nmesh=360, rsd=rsd_bool)
