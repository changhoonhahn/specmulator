'''

forward model for HOD parameters

here halo catalog is always fixed to the following: 
    neut.data.NeutHalos(0.0, 1, 4) 


'''
import numpy as np 
# -- local -- 
from . import lhd_hod
from neut import data as NeutDat 
from specmulator import util as UT 


def X_HOD_LHD(seed_hod, obvs='plk', 
        ell=0, Nmesh=360, rsd=True, krange=[0.01, 0.5], karr=False, # kwargs specifying P(k)
        prior='sinha2017prior_narrow', samples=40, method='mdu', # kwargs specifying the LHD 
        silent=True):
    ''' Read in P(k)s for a LHD of HOD parameters (specified using prior) 
    '''
    if ell not in [0,2,4]: 
        raise ValueError("only monopole, quadrupole, and hexadecapole") 

    # read in prior object 
    hod_prior = lhd_hod.HODprior(prior) 
    
    # directory where file is stored 
    f_dir = ''.join(['/Users/chang/projects/specmulator/dat/lhd/onlyHOD/', 
        method, '_', str(samples), '_', prior, '/']) 

    plks = [] 
    for i_p in range(samples): 
        try: 
            plk_i = HODLHD_NeutObvs('plk', mneut, nreal, nzbin, seed_hod, i_p,
                    HODrange=HODrange, method=method, samples=samples, Nmesh=Nmesh, 
                    rsd=rsd, overwrite=False, silent=silent)
        except (ValueError, IOError): 
            print("%i of the LHD not read" % i_p)
            continue 
        if krange is not None: 
            klim = np.where((plk_i['k'] > krange[0]) & (plk_i['k'] < krange[1]))
            plks.append(plk_i['p'+str(ell)+'k'][klim])
        else:
            plks.append(plk_i['p'+str(ell)+'k'])
    if krange is not None: 
        k = plk_i['k'][klim]
    else:
        k = plk_i['k']
    if not karr: 
        return np.array(plks) 
    else: 
        return k, np.array(plks) 


def X_testHOD_LHD(seed_hod, obvs='plk', 
        ell=0, Nmesh=360, rsd=True, krange=[0.01, 0.5], karr=False, # kwargs specifying P(k)
        prior='sinha2017prior_narrow', samples=40, method='mdu', # kwargs specifying the LHD 
        silent=True):
    ''' Read in P(k)s for a LHD of HOD parameters (specified using prior) 
    '''
    if ell not in [0,2,4]: 
        raise ValueError("only monopole, quadrupole, and hexadecapole") 

    # read in prior object 
    hod_prior = lhd_hod.HODprior(prior) 


def HODcatalog(theta, seed_hod, LOS=[0,0,1], write=None, silent=True): 
    ''' Generate HOD catalogs by populating the halo catalog based on 
    given HOD parameter values and random seed. Wrapper for `Neut.forwardmodel.Galaxies`

    Parameters
    ----------
    seed_hod : int, 
        random seed for the HOD 

    LOS : list, optional 
        3 element list of integers 0 or 1 that specify the the line-of-sight direction 
    '''
    if len(theta.shape) != 1: 
        raise ValueError("can only read one theta at a time") 

    # read in hardcoded Neutrino halo (this is hardcoded in, do *not* change) 
    halos = NeutHalos(0.0, 1, 4) 

    # HOD keys
    hodkeys = ['logMmin', 'sigma_logM', 'logM0', 'logM1', 'alpha'] 
    p_hod = {} 
    for i_k, k in enumerate(hodkeys): 
        p_hod[k] = theta[i_k] 
    if not silent: 
        print("HOD parameters") 
        print(p_hod) 
    # populate the halo catalogs using HOD 
    gals = FM.Galaxies(halos, p_hod, seed=seed_hod)  
    
    # RSD position (hardcoded in the z direction) 
    gals['RSDPosition'] = FM.RSD(gals, LOS=[0,0,1]) 
    return gals  
