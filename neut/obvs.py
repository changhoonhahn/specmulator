'''

code for dealing with observables

'''
import os 
import numpy as np 

import util as UT 


def Plk_halo(mneut, nreal, nzbin, zspace=False): 
    ''' Return the powerspectrum multipoles for halo catalog with
    total neutrino mass `mneut`, realization `nreal`, and redshift 
    bin `nzbin` in either real/redshift space. 

    '''
    if mneut == 0.1: dir = ''.join([UT.dat_dir(), '0.10eV/'])
    else: dir = ''.join([UT.dat_dir(), str(mneut), 'eV/'])
    if zspace: str_space = 'z' # redhsift
    else: str_space = 'r' # real 
    f = ''.join([dir, 'plk.groups.', str(mneut), 'eV.', str(nreal), '.nzbin', str(nzbin), 
        '.', str_space, 'space.dat']) 
    if not os.path.isfile(f): raise ValueError('file does not exist') 

    # read in plk 
    k,p0k,p2k,p4k = np.loadtxt(f, skiprows=3, unpack=True, usecols=[0,1,2,3]) 
    plk = {'k': k, 'p0k': p0k, 'p2k': p2k, 'p4k':p4k} 

    # readin shot-noise from header 
    with open(f) as lines: 
        for i_l, line in enumerate(lines):
            if i_l == 1: 
                str_sn = line
                break 
    plk['shotnoise'] = float(str_sn.strip().split('P_shotnoise')[-1])
    return plk 


def threePCF_halo(mneut, nreal, nzbin, zspace=False): 
    ''' Read in three point data. Need to figure out how I want to do this...
    '''
    pass 
