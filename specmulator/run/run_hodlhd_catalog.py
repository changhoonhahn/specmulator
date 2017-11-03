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


def hodlhd_catalogs(mneut, nreal, nzbin, seed_hod, HODrange='sinha2017prior_narrow', method='mdu', samples=17): 
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
    if HODrange in ['sinha2017prior', 'sinha2017prior_narrow']:  
        keylist = ['logMmin', 'sigma_logM', 'logM0', 'logM1', 'alpha'] 
    lhcube = lhd.HOD_LHD(HODrange=HODrange, samples=samples, method=method)

    # import Neutrino halo with mneut eV, realization # nreal, at z specified by nzbin 
    halos = Dat.NeutHalos(mneut, nreal, nzbin) 

    for i_p in range(samples): 
        print('%i of %i LHD'%(i_p+1,samples))
        p_hod = {} 
        for ik, k in enumerate(keylist): 
            p_hod[k] = lhcube[i_p,ik]
        print(p_hod)
        # populate the halo catalogs using HOD 
        gals = FM.Galaxies(halos, p_hod, seed=seed_hod)  
        
        # RSD position (hardcoded in the z direction) 
        gals['RSDPosition'] = FM.RSD(gals, LOS=[0,0,1]) 
        # save to file 
        folder = ''.join([UT.dat_dir(), 'lhd/', str(mneut), 'eV_', str(nreal), '_z', str(nzbin), '/']) 
        if not os.path.exists(folder): # make directory
            os.mkdir(folder) 
        gals.save('%s/HOD%s_seed%i_%i'%(folder,method,seed_hod,i_p), ('Position', 'Velocity', 'RSDPosition'))
        del gals

    return None 


def hodlhd_observables(obvs, mneut, nreal, nzbin, seed_hod, HODrange='sinha2017prior_narrow', method='nohl', samples=17): 
    ''' Calculate the observables of the HOD LHD catalogs
    '''
    pass


if __name__=="__main__": 
    # e.g. python run/run_hodlhd_catalog.py catalog 0.0 1 4 1
    mod = Sys.argv[1]
    if mod == 'lhd': 
        HOD_LHD()
    elif mod in ['catalog', 'observable']: 
        mneut = float(Sys.argv[2])
        nreal = int(Sys.argv[3]) 
        nzbin = int(Sys.argv[4]) 
        seed_hod = int(Sys.argv[5]) 

        print('%f eV'%(mneut))
        print('realization %i'%(nreal))
        print('nzbin %i'%(nzbin))
        print('random seed',seed_hod)
    
        if mod == 'catalog': 
            hodlhd_catalogs(mneut, nreal, nzbin, seed_hod, HODrange='sinha2017prior_narrow', method='mdu', samples=17)
