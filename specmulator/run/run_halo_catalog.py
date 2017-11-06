'''

Code for generating observables for halo catalogs


'''
import os 
import sys as Sys
import numpy as np 

import env 
import lhd
import util as UT 
import data as Dat
import forwardmodel as FM 


def NeutHalo_Plk(mneut, nreal, nzbin): 
    ''' Calculate the powerspectrum multipoles for Paco's Neutrio 
    halo catalogs
    '''
    # import Neutrino halo with mneut eV, realization # nreal, at z specified by nzbin 
    halos = Dat.NeutHalos(mneut, nreal, nzbin) 

    plk = FM.Observables(halos, observable='plk', rsd=False, Nmesh=360)

    # save somehow 
    return None 


def NeutHalo_pre3PCF(mneut, nreal, nzbin, zspace=False): 
    ''' Pre-process halo catalogs for 3PCF run. Read in halo catalog
    and output to input file format for Daniel Eisenstein's code. 

    parameters
    ----------
    zspace : bool, optional 
        if True, calculates redshift space positions. Otherwise real space. 

    '''
    # import Neutrino halo with mneut eV, realization # nreal, at z specified by nzbin 
    halos = Dat.NeutHalos(mneut, nreal, nzbin) 
    if zspace: 
        halos['RSDPosition'] = FM.RSD(halos, LOS=[0,0,1])
    
    # halo positions
    x = halos['Position'][:,0]
    y = halos['Position'][:,1]
    z = halos['Position'][:,2]
    # weights (all ones!) 
    w = np.ones(len(x)) 
    
    # output file 
    if mneut == 0.1: dir = ''.join([UT.dat_dir(), '0.10eV/', str(nreal)])
    else: dir = ''.join([UT.dat_dir(), str(mneut), 'eV/', str(nreal)])
    if zspace: str_space = '.z'
    else: str_space = '.r'
    fout = ''.join([dir, '/groups_', str(nzbin).zfill(3), '/groups', str_space, 'space.dat']) 
    
    # header 
    hdr = ''.join(['m neutrino = ', str(mneut), ' eV, realization ', str(nreal), ', zbin ', str(nzbin)]) 

    # write to file 
    np.savetxt(fout, np.array([x, y, z, w]).T, header=hdr) 
    return None 


if __name__=="__main__": 
    arg1 = Sys.argv[1]    

    if arg1 == '3pcf': 
        mneut = float(Sys.argv[2])
        nreal = int(Sys.argv[3])
        nzbin = int(Sys.argv[4])
        zstr = Sys.argv[5]
        if zstr == 'z': 
            zbool = True
        elif zstr == 'real': 
            zbool = False
        else: 
            raise ValueError

        NeutHalo_pre3PCF(mneut, nreal, nzbin, zspace=zbool)
    else: 
        raise ValueError

