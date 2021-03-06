'''



'''
import os
import numpy as np 
# --- local --- 
from specmulator import lhd
from specmulator import util as UT 


# --- Data -- 
def X_HODLHD(seed_hod, obvs='plk', 
        ell=0, Nmesh=360, rsd=True, karr=False, # kwargs specifying P(k)
        prior='sinha2017prior_narrow', samples=40, method='mdu', ndim=None, # kwargs for LHD 
        silent=True):
    ''' Read observable (e.g. P(k)) for a LHD of HOD parameters (specified by prior, 
    samples, and method) 
    '''
    # directory of the catalogs 
    if ndim is None: str_ndim = ''
    else: str_ndim = '_'+str(ndim)+'D'
    f_dir = ''.join([UT.dat_dir(), 'lhd/onlyHOD/', 
        method, '_', str(samples), '_', prior, str_ndim, '/']) 

    # rsd flag 
    if rsd: rsd_str = 'z'
    else: rsd_str = 'r'

    # read in observable
    pks = [] 
    for i in range(samples):  
        # 21st LHD HOD theta requires ton of memory 
        # this should be addressed more robustly but for
        # the sake of my sanity at the moment, I will skip 
        # it for now. 
        if i == 21: continue 
        # directory of theta_i,LHD 
        fname = ''.join([f_dir, 'HOD', method, '_seed', str(seed_hod), '_', str(i), '/', 
           'pk.menut0.0.nreal1.nzbin4.', rsd_str, 'space.', str(Nmesh), '.nbkt.dat'])   
        if not silent: 
            print('%i -- reading %s' % (i, fname))  

        # read in plks 
        k, pk = np.loadtxt(fname, skiprows=4, unpack=True, usecols=[0,1+ell/2])
        pks.append(pk)

    if karr: 
        return k, np.array(pks) 
    else: 
        return np.array(pks) 


def X_testHODLHD(seed_hod, obvs='plk', 
        ell=0, Nmesh=360, rsd=True, karr=False,     # kwargs specifying P(k)
        prior='sinha2017prior_narrow', samples=20, ndim=None, # kwargs specifying the LHD 
        silent=True):
    ''' Read observable (e.g. P(k)) for a LHD of HOD parameters (specified by prior, 
    samples, and method) 
    '''
    # directory of the catalogs 
    if ndim is None: str_ndim = ''
    else: str_ndim = '_'+str(ndim)+'D'
    f_dir = ''.join([UT.dat_dir(), 'lhd/onlyHOD/test_', str(samples), '_', prior, str_ndim, '/']) 

    # rsd flag 
    if rsd: rsd_str = 'z'
    else: rsd_str = 'r'

    # read in observable
    pks = [] 
    for i in range(samples):  
        # directory of theta_i,LHD 
        fname = ''.join([f_dir, 'testHOD_seed', str(seed_hod), '_', str(i), '/',
            'pk.menut0.0.nreal1.nzbin4.', rsd_str, 'space.', str(Nmesh), '.nbkt.dat'])   
        if not silent: 
            print('%i -- reading %s' % (i, fname))  

        # read in plks 
        k, pk = np.loadtxt(fname, skiprows=4, unpack=True, usecols=[0,1+ell/2])
        pks.append(pk)

    if karr: 
        return k, np.array(pks) 
    else: 
        return np.array(pks) 


def _check_X_HODLHD(seed_hod, obvs='plk', 
        ell=0, Nmesh=360, rsd=True, karr=False, # kwargs specifying P(k)
        prior='sinha2017prior_narrow', samples=40, method='mdu',    # kwargs specifying the LHD 
        silent=True):
    ''' check whether all the X_HODLHD files are there
    '''
    f_dir = ''.join([UT.dat_dir(), 'lhd/onlyHOD/', method, '_', str(samples), '_', prior, '/']) 
    # rsd flag 
    if rsd: rsd_str = 'z'
    else: rsd_str = 'r'

    # read in observable
    for i in range(samples):  
        if i == 21: continue # 21 requires ton of memory 
        # directory of theta_i,LHD 
        fname = ''.join([f_dir, 'HOD', method, '_seed', str(seed_hod), '_', str(i), '/', 
           'pk.menut0.0.nreal1.nzbin4.', rsd_str, 'space.', str(Nmesh), '.nbkt.dat'])   
        if not os.path.isfile(fname):
            print('%s does not exist' % fname) 
    return None 


# -- HOD LHD -- 
def HOD_LHD(prior=None, samples=None, method=None, ndim=None, overwrite=False):
    ''' Return latin hypercubes with `samples` elements using `method` method 
    that spans the specified prior of the vanilla Zheng et al.(2007) HOD model. 

    References
    -------
    - Zheng Z., et al. (2007) -- arXiv:0512071
    - Sinha M., et al. (2017) -- arXiv:1708.04892 
    '''
    if ndim is None: str_ndim = ''
    else: str_ndim = str(ndim)+'D.'
    fname = ''.join([UT.dat_dir(), 'lhd/onlyHOD/', 
        'HOD_LHD.', str_ndim, prior, '.', method, '.', str(samples), 'samples.dat']) 
    
    if os.path.isfile(fname) and not overwrite: # file exists 
        lhcube = np.loadtxt(fname, skiprows=4)
        if ndim is None: 
            # skip 21st LHD parameter
            indx = np.ones(lhcube.shape[0], dtype=bool)
            indx[21] = False
            lhcube = lhcube[indx,:]
    else: 
        hodprior = HODprior(prior) 
        hod_range = hodprior.range()

        if ndim is None: # Nsample x Ndim latin hypercube
            lhcube = lhd.thetaLHD(hod_range, samples=samples, method=method, ndim=ndim)
        elif ndim == 1: # only vary alpha 
            lhcube = np.zeros((samples, len(hod_range[0])))
            for i in range(4): 
                lhcube[:,i] = 0.5*(hod_range[0]+hod_range[1])[i]
            lhc = lhd.thetaLHD([[hod_range[0][4]], [hod_range[1][4]]],
                    samples=samples, method=method, ndim=ndim)
            lhcube[:,4] = lhc.flatten()

        f = open(fname, 'w') 
        f.write('# '+prior+'\n') 
        f.write('# parameters : '+', '.join(hodprior.labels)+'\n') 
        f.write('# theta min : '+', '.join([str(r) for r in hod_range[0]])+'\n') 
        f.write('# theta max : '+', '.join([str(r) for r in hod_range[1]])+'\n') 

        for i in range(lhcube.shape[0]): 
            f.write('%f' % (lhcube[i,0]))
            for j in range(1, lhcube.shape[1]):
                f.write('\t %f' % (lhcube[i,j]))
            f.write('\n')
        f.close() 
    return lhcube 


def testHOD_LHD(prior=None, samples=None, overwrite=False, ndim=None, seed=1):
    ''' Return parameters of the vanilla Zheng et al.(2007) HOD model 
    for testing the LHD. These values will be significantly more within 
    the boundaries of the parameter space 
    '''
    if ndim is None: str_ndim = ''
    else: str_ndim = str(ndim)+'D.'
    fname = ''.join([UT.dat_dir(), 'lhd/onlyHOD/',
        'HOD_LHDtest.', str_ndim, prior, '.rseed', str(seed), '.', str(samples), 'samples.dat']) 
    
    if os.path.isfile(fname) and not overwrite: # file exists 
        theta_test = np.loadtxt(fname, skiprows=4)
    else: 
        # prior and its range 
        hodprior = HODprior(prior) 
        hod_range = hodprior.range()

        np.random.seed(seed)

        theta_test = np.zeros((samples, len(hod_range[0])))
        if ndim is None: 
            for i in range(theta_test.shape[1]):
                dHOD_range = (hod_range[1] - hod_range[0])[i] 
                # linearly spaced in the inner 50% of the parameter space 
                ti = np.linspace(hod_range[0][i] + 0.25*dHOD_range, 
                        hod_range[1][i]-0.25*dHOD_range, 
                        theta_test.shape[0])
                np.random.shuffle(ti)
                theta_test[:,i] = ti
        elif ndim == 1: 
            # only alpha varies
            dHOD_range = (hod_range[1] - hod_range[0])[4] 
            # linearly spaced in the inner 50% of the parameter space 
            ti = np.linspace(hod_range[0][4] + 0.25*dHOD_range, 
                    hod_range[1][4]-0.25*dHOD_range, 
                    theta_test.shape[0])
            np.random.shuffle(ti)
            theta_test[:,4] = ti
            theta_test[:,0] = 0.5 * (hod_range[0] + hod_range[1])[0]
            theta_test[:,1] = 0.5 * (hod_range[0] + hod_range[1])[1]
            theta_test[:,2] = 0.5 * (hod_range[0] + hod_range[1])[2]
            theta_test[:,3] = 0.5 * (hod_range[0] + hod_range[1])[3]

        f = open(fname, 'w') 
        f.write('# '+prior+'\n') 
        f.write('# parameters : '+', '.join(hodprior.labels)+'\n') 
        f.write('# theta min : '+', '.join([str(r) for r in hod_range[0]])+'\n') 
        f.write('# theta max : '+', '.join([str(r) for r in hod_range[1]])+'\n') 

        for i in range(theta_test.shape[0]): 
            f.write('%f' % (theta_test[i,0]))
            for j in range(1, theta_test.shape[1]):
                f.write('\t %f' % (theta_test[i,j]))
            f.write('\n')
        f.close() 
    return theta_test 


class HODprior(object): 
    ''' class object to handle parameter priors for Zheng et al. (2007)
    HOD model. Parameters are: log M_min, sigma_logM, log M_0, log M_1, alpha
    '''
    def __init__(self, name):  
        if name not in ['sinha2017prior_narrow', 'sinha2017prior']: 
            raise NotImplementedError("%s prior not yet implement" % name) 
        self.name = name 
        self.labels = ['log $M_\mathrm{min}$', '$\sigma_{\mathrm{log}\,M}$', 'log $M_0$', 'log $M_1$', r'$\alpha$']

    def range(self):
        ''' range of the prior distribution 
        '''
        if self.name == 'sinha2017prior_narrow': # narrow alpha range 
            # Sinha et al (2017) prior with narrow alpha range
            range_min = [11., 0.001, 6., 12., 0.5]
            range_max = [12.2, 1., 14., 14., 1.5]
        elif self.name == 'sinha2017prior':  
            # Sinha et al (2017) prior
            range_min = [11., 0.001, 6., 12., 0.001]
            range_max = [12.2, 1., 14., 14., 2.]
        return [np.array(range_min), np.array(range_max)]
