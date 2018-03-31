'''

functions for generating all the necessary data for
`specmulator.onlyhod.py`

'''
import os 
import sys 
import numpy as np 
import nbodykit.lab as NBlab
# -- local -- 
import env
from neut import data as Dat
from neut import forwardmodel as FM
from specmulator import util as UT
from specmulator import onlyhod as onlyHOD


def HODLHDpks(seed_hod, i_lhd, 
        prior='sinha2017prior_narrow', method='mdu', samples=17, ndim=None,
        Nmesh=360, rsd=True): 
    ''' Measure and write P(k|theta_i,LHD) for the HOD LHD that spans the `prior` 
    '''
    # RSD flag 
    if rsd: str_rsd = '.zspace'
    else: str_rsd = '.rspace'

    # directory 
    if ndim is None: str_ndim = ''
    else: str_ndim = '_'+str(ndim)+'D'
    parent_dir = ''.join([UT.dat_dir(), 'lhd/onlyHOD/', 
        method, '_', str(samples), '_', prior, str_ndim, 
        '/HOD', method, '_seed', str(seed_hod), '_', str(i_lhd), '/']) 
    fname = ''.join([parent_dir,
        'pk.menut0.0.nreal1.nzbin4', str_rsd, '.', str(Nmesh), '.nbkt.dat'])
     
    if not os.path.exists(parent_dir): # make directory
        raise ValueError("the directory \n %s \n does not exist!" % parent_dir) 

    # read from galaxy catalog  
    gals = NBlab.BigFileCatalog(parent_dir, header='Header')

    plk = FM.Observables(gals, observable='plk', rsd=rsd, Nmesh=Nmesh)

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
    return None 


def HODLHDcatalogs(seed_hod, i_lhd, 
        prior='sinha2017prior_narrow', method='mdu', samples=40, ndim=None): 
    ''' Generate HOD catalog of the i_lhd th parameter of the HOD LHD 
    that spans the Sinha M., et al. (2017) HOD parameter priors. 

    parameters
    ----------
    seed_hod : int, 
        random seed for the HOD 

    prior : str, optional
        string specifying the HOD range. Default is 'sinha2017prior', which uses
        the prior from Sinha et al. (2017) 

    method : str 
        string specifying the method of LHD. Default is nohl 

    samples : int, optional 
        sample size of the LHD. Default is 17, which is from the fixed sample size of NOHL method
    '''
    if i_lhd >= samples: 
        raise ValueError("i_lhd has to be less than samples") 

    # read in  Neutrino halo with m_neut = 0.0 eV, realization #1, at z_bin = 4 
    # this is hardcoded; don't change this 
    halos = Dat.NeutHalos(0.0, 1, 4)

    # read in the HOD LHD 
    lhcube = onlyHOD.HOD_LHD(prior=prior, samples=samples, method=method, ndim=ndim)
        
    print('%i of %i LHD' % (i_lhd+1,samples))
    p_hod = {}
    hodkeys = ['logMmin', 'sigma_logM', 'logM0', 'logM1', 'alpha']
    for ik, k in enumerate(hodkeys):
        p_hod[k] = lhcube[i_lhd,ik]
    print(p_hod)

    # populate the halo catalogs using HOD
    gals = FM.Galaxies(halos, p_hod, seed=seed_hod)

    # RSD position (hardcoded in the z direction)
    gals['RSDPosition'] = FM.RSD(gals, LOS=[0,0,1])

    # directory of the catalogs 
    if ndim is None: str_ndim = ''
    else: str_ndim = '_'+str(ndim)+'D'
    parent_dir = ''.join([UT.dat_dir(), 'lhd/onlyHOD/', 
        method, '_', str(samples), '_', prior, str_ndim,
        '/HOD', method, '_seed', str(seed_hod), '_', str(i_lhd), '/']) 
    print('writing to %s ---------' % parent_dir) 
    if not os.path.exists(parent_dir): # make directory
        os.mkdir(parent_dir)
    # save to file
    gals.save(parent_dir, ('Position', 'Velocity', 'RSDPosition'))
    return None


def testHODLHDpks(seed_hod, i_lhd, 
        prior='sinha2017prior_narrow', samples=17, ndim=None, Nmesh=360, rsd=True): 
    ''' Measure and write P(k|theta_i,LHD) for the HOD LHD that spans the `prior` 
    '''
    # RSD flag 
    if rsd: str_rsd = '.zspace'
    else: str_rsd = '.rspace'

    # directory 
    if ndim is None: str_ndim = ''
    else: str_ndim = '_'+str(ndim)+'D'
    parent_dir = ''.join([UT.dat_dir(), 'lhd/onlyHOD/test_', str(samples), '_', prior, str_ndim,
        '/testHOD_seed', str(seed_hod), '_', str(i_lhd), '/']) 
    fname = ''.join([parent_dir,
        'pk.menut0.0.nreal1.nzbin4', str_rsd, '.', str(Nmesh), '.nbkt.dat'])
    if not os.path.exists(parent_dir): # make directory
        raise ValueError("the directory \n %s \n does not exist!" % parent_dir) 

    # read from galaxy catalog  
    gals = NBlab.BigFileCatalog(parent_dir, header='Header')

    plk = FM.Observables(gals, observable='plk', rsd=rsd, Nmesh=Nmesh)

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
    return None 


def testHODLHDcatalogs(seed_hod, i_lhd, prior='sinha2017prior_narrow', samples=20, ndim=None): 
    ''' Generate HOD catalog of the i_lhd th parameter of the HOD LHD 
    that spans the Sinha M., et al. (2017) HOD parameter priors. 

    parameters
    ----------
    seed_hod : int, 
        random seed for the HOD 

    prior : str, optional
        string specifying the HOD range. Default is 'sinha2017prior', which uses
        the prior from Sinha et al. (2017) 

    samples : int, optional 
        size of test sample 

    LOS : list, optional 
        3 element list of integers 0 or 1 that specify the the line-of-sight direction 
    '''
    if i_lhd >= samples: 
        raise ValueError("i_lhd has to be less than samples") 

    # read in  Neutrino halo with m_neut = 0.0 eV, realization #1, at z_bin = 4 
    # this is hardcoded; don't change this 
    halos = Dat.NeutHalos(0.0, 1, 4)

    # read in test sample for the HOD LHD 
    lhcube = onlyHOD.testHOD_LHD(prior=prior, samples=samples, ndim=ndim)
        
    print('%i of %i LHD' % (i_lhd+1,samples))
    p_hod = {}
    hodkeys = ['logMmin', 'sigma_logM', 'logM0', 'logM1', 'alpha']
    for ik, k in enumerate(hodkeys):
        p_hod[k] = lhcube[i_lhd,ik]
    print(p_hod)

    # populate the halo catalogs using HOD
    gals = FM.Galaxies(halos, p_hod, seed=seed_hod)

    # RSD position (hardcoded in the z direction)
    gals['RSDPosition'] = FM.RSD(gals, LOS=[0,0,1])

    # directory of the catalogs 
    if ndim is None: str_ndim = ''
    else: str_ndim = '_'+str(ndim)+'D'
    parent_dir = ''.join([UT.dat_dir(), 'lhd/onlyHOD/', 
        'test_', str(samples), '_', prior, str_ndim, 
        '/testHOD_seed', str(seed_hod), '_', str(i_lhd), '/']) 
    print('writing to %s ---------' % parent_dir) 
    if not os.path.exists(parent_dir): # make directory
        os.mkdir(parent_dir)
    # save to file
    gals.save(parent_dir, ('Position', 'Velocity', 'RSDPosition'))
    return None


if __name__=="__main__": 
    lhd_or_test = sys.argv[1]
    nsample = int(sys.argv[2])
    seed_hod = int(sys.argv[3]) 
    i_lhd = int(sys.argv[4])
    ndim = int(sys.argv[5]) 
    if ndim == 5: 
        ndim = None 

    if lhd_or_test == 'lhd': 
        print('constructing LHD sample %i' % i_lhd) 
        print('with random seed %i ' % seed_hod)
        # construct catalog
        HODLHDcatalogs(seed_hod, i_lhd, 
                prior='sinha2017prior_narrow', method='mdu', samples=nsample, ndim=ndim)
        # real space P(k)
        HODLHDpks(seed_hod, i_lhd, 
                prior='sinha2017prior_narrow', method='mdu', samples=nsample, ndim=ndim,
                Nmesh=360, rsd=False)
        # redshift space P(k)
        HODLHDpks(seed_hod, i_lhd, 
                prior='sinha2017prior_narrow', method='mdu', samples=nsample, ndim=ndim, 
                Nmesh=360, rsd=True)
    elif lhd_or_test == 'test': 
        print('constructing test sample %i' % i_lhd) 
        print('with random seed %i ' % seed_hod)
        # construct catalog
        testHODLHDcatalogs(seed_hod, i_lhd, 
                prior='sinha2017prior_narrow', samples=nsample, ndim=ndim)
        # real space P(k)
        testHODLHDpks(seed_hod, i_lhd, 
                prior='sinha2017prior_narrow', samples=nsample, ndim=ndim,
                Nmesh=360, rsd=False)
        # redshift space P(k)
        testHODLHDpks(seed_hod, i_lhd, 
                prior='sinha2017prior_narrow', samples=nsample, ndim=ndim,
                Nmesh=360, rsd=True)
