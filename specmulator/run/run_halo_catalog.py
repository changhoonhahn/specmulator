'''

Code for generating observables for halo catalogs


'''
import os 
import time 
import sys as Sys
import numpy as np 

import env 
import lhd
import util as UT 
import data as Dat
import forwardmodel as FM 


def NeutHalo_Plk(mneut, nreal, nzbin, zspace=False): 
    ''' Calculate the powerspectrum multipoles for Paco's Neutrio 
    halo catalogs
    '''
    # import Neutrino halo with mneut eV, realization # nreal, at z specified by nzbin 
    halos = Dat.NeutHalos(mneut, nreal, nzbin) 
    if zspace: 
        halos['RSDPosition'] = FM.RSD(halos, LOS=[0,0,1])
    
    # calculate P_l(k) 
    plk = FM.Observables(halos, observable='plk', rsd=zspace, Nmesh=360)

    # output file 
    if mneut == 0.1: dir = ''.join([UT.dat_dir(), '0.10eV/'])
    else: dir = ''.join([UT.dat_dir(), str(mneut), 'eV/'])
    if zspace: str_space = '.z'
    else: str_space = '.r'
    fout = ''.join([dir, 'plk.groups.', str(mneut), 'eV.', str(nreal), '.nzbin', str(nzbin), str_space, 'space.dat']) 

    # header 
    hdr = ''.join(['P_l(k) measurements for m neutrino = ', str(mneut), ' eV, realization ', str(nreal), ', zbin ', str(nzbin), 
        '\n P_shotnoise ', str(plk['shotnoise']), 
        '\n cols: k, P_0, P_2, P_4']) 

    # write to file 
    np.savetxt(fout, np.array([plk['k'], plk['p0k'], plk['p2k'], plk['p4k']]).T, header=hdr) 
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
    xyz = np.array(halos['Position'])
    x = xyz[:,0]
    y = xyz[:,1]
    z = xyz[:,2]
    # weights (all ones!) 
    w = np.ones(len(x)) 
    
    # output file 
    if mneut == 0.1: dir = ''.join([UT.dat_dir(), '0.10eV/', str(nreal)])
    else: dir = ''.join([UT.dat_dir(), str(mneut), 'eV/', str(nreal)])
    if zspace: str_space = '.z'
    else: str_space = '.r'
    fout = ''.join([dir, '/groups.nzbin', str(nzbin), '.', str_space, 'space.dat']) 
    
    # header 
    hdr = ''.join(['m neutrino = ', str(mneut), ' eV, realization ', str(nreal), ', zbin ', str(nzbin)]) 

    outarr = np.array([x, y, z, w]).T
    # write to file 
    np.savetxt(fout, outarr, header=hdr) 
    print('--- halo written to ---\n %s' % (fout))
    return None 


if __name__=="__main__": 
    arg1 = Sys.argv[1] # 3pcf/plk etc 

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

    if arg1 == '3pcf': 
        NeutHalo_pre3PCF(mneut, nreal, nzbin, zspace=zbool)
    elif arg1 == 'plk': 
        NeutHalo_Plk(mneut, nreal, nzbin, zspace=zbool)
    else: 
        raise ValueError
