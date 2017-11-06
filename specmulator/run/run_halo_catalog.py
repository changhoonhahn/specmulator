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


def NeutHalo_pre3PCF(mneut, nreal, nzbin): 
    ''' Pre-process halo catalogs for 3PCF run. Read in halo catalog
    and output to input file format for Daniel Eisenstein's code. 
    '''
    # import Neutrino halo with mneut eV, realization # nreal, at z specified by nzbin 
    halos = Dat.NeutHalos(mneut, nreal, nzbin) 

    # write to file 
    halos['


