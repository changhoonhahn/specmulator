'''

Code for generating galaxy catalogs from halo catalog  
for the HOD parameters sampled by a LHD 


'''
import os 
import sys
import numpy as np 

from feasibgs import lhd
from feasibgs import util as UT 
from feasibgs import data as Dat
from feasibgs import forwardmodel as FM 


def HODLHD_test(HODrange='sinha2017prior_narrow', samples=20, seed=1, overwrite=True): 
    ''' Build `samples` test parameter points  
    ''' 
    lhd.HOD_LHD_test(HODrange=HODrange, samples=samples, overwrite=overwrite, seed=seed) 
    return None 


def Obvs_HODLHD_test(theta, halo, mneut=0.0, nreal=1, nzbin=4, seed=1, obvs='plk', Nmesh=360): 
    ''' Build the observables `obvs` for {theta_test} given halo catalog
    '''
    folder = ''.join([UT.dat_dir(), 
        'lhd/', str(mneut), 'eV_', str(nreal), '_z', str(nzbin), '_', str(samples), 'samples/', 
        'HOD', method, '_seed', str(seed_hod), '_', str(i_p), '/']) 

    # F(theta) --- i.e. the galaxy catalog generated
    # from the halo catalog
    p_hod = {'logMmin': theta[0], 'sigma_logM': theta[1], 'logM0': theta[2], 'logM1': theta[3], 'alpha': theta[4]}
    g = FM.Galaxies(halos, p_hod, seed=seed)
    g['RSDPosition'] = FM.RSD(g, LOS=[0,0,1]) # impose RSD

    # measure P(k) from F(theta) in real and z space 
    for rsd in [False, True]: 
        if obvs == 'plk': # power spectrum multipole 
            plk = FM.Observables(gals, observable='plk', rsd=rsd, Nmesh=Nmesh) 
            
            # write to file 
            if rsd: str_rsd = '.zspace'
            else: str_rsd = '.rspace'
            fname = ''.join([folder, 
                 'pk.theta_test.menut', str(mneut), '.nreal', str(nreal), '.nzbin', str(nzbin), str_rsd, 
                 '.', str(Nmesh), '.nbkt.dat'])
            
            # save to file 
            f = open(fname, 'w')
            f.write("### header ### \n")
            f.write("# shotnoise %f \n" % plk['shotnoise'])
            f.write("# columns : k , P0, P2, P4 \n")
            f.write('### header ### \n') 

            for ik in range(len(plk['k'])): 
                f.write("%f \t %f \t %f \t %f" % (plk['k'][ik], plk['p0k'][ik], plk['p2k'][ik], plk['p4k'][ik]))
                f.write("\n") 
            f.close() 
            obvs = plk
        else: 
            raise NotImplementedError('only Plk implemented') 
    return None 


if __name__=="__main__": 
    mneut = 0.0 
    nreal = 1  
    nzbin = 4

    mod = sys.argv[1]
    nsample = int(sys.argv[2]) 
    if mod == 'lhd': HODLHD_test(samples=nsample, seed=1)
    elif mod == 'observable': 
        # read in {theta_test} 
        theta_test = lhd.HOD_LHD_test(HODrange='sinha2017prior_narrow', samples=nsample, 
                overwrite=False, seed=1) 
        # Read in halo catalog
        halos = Dat.NeutHalos(mneut, nreal, nzbin)

        for i in range(nsample):  
            Obvs_HODLHD_test(theta_test[i,:], halos, mneut=mneut, nreal=nreal, nzbin=nzbin, 
                    seed=i, obvs='plk', Nmesh=360) 
