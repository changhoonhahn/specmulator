'''

Code for the Latin Hypercube Design 


'''
import os 
import random
import numpy as np 

import util as UT
import data as Dat

try: 
    import pynolh as pyNOLH
    import pyDOE
    import lhsmdu as lhsMDU
except ModuleNotFoundError: 
    pass 


def HODLHD_NeutCatalog(mneut, nreal, nzbin, seed_hod, i_p, HODrange='sinha2017prior_narrow', method='mdu', samples=17): 
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
    folder = ''.join([UT.dat_dir(), 
        'lhd/', str(mneut), 'eV_', str(nreal), '_z', str(nzbin), '_', str(samples), 'samples/', 
        'HOD', method, '_seed', str(seed_hod), '_', str(i_p), '/']) 
    
    assert i_p < samples
        
    # read in  Neutrino halo with mneut eV, realization # nreal, at z specified by nzbin 
    halos = Dat.NeutHalos(mneut, nreal, nzbin) 

    if not np.all([os.path.exists(folder+subfold+'/') for subfold in ['Position', 'Velocity', 'RSDPosition']]):   
        # generate the LHD HOD catalog  
        if HODrange in ['sinha2017prior', 'sinha2017prior_narrow']:  
            keylist = ['logMmin', 'sigma_logM', 'logM0', 'logM1', 'alpha'] 
        lhcube = HOD_LHD(HODrange=HODrange, samples=samples, method=method)

        #for i_p in range(samples): 
        print('%i of %i LHD'%(i_p+1,samples))
        p_hod = {} 
        for ik, k in enumerate(keylist): 
            p_hod[k] = lhcube[i_p,ik]
        print(p_hod)
        # populate the halo catalogs using HOD 
        gals = FM.Galaxies(halos, p_hod, seed=seed_hod)  
        
        # RSD position (hardcoded in the z direction) 
        gals['RSDPosition'] = FM.RSD(gals, LOS=[0,0,1]) 

        parent_dir = '/'.join(folder[:-1].split('/')[:-1])+'/'
        if not os.path.exists(parent_dir): # make directory
            os.mkdir(parent_dir) 
        # save to file 
        gals.save(folder, ('Position', 'Velocity', 'RSDPosition'))
        del gals
    else:
        # read from file 
        gals = NBlab.BigFileCatalog(folder, header='Header')
        gals.attrs = halos.attrs.copy() 
    return None 


def HOD_LHD(HODrange=None, samples=None, method=None):
    ''' Return latin hypercubes for vanilla Zheng et al.(2007) 
    HOD model. 

    Sources
    -------
    - Zheng Z., et al. (2007) -- arXiv:0512071
    - Sinha M., et al. (2017) -- arXiv:1708.04892 
    '''
    if HODrange == 'sinha2017prior_narrow': # narrow alpha range 
        range_descrip = "Sinha et al (2017) prior with narrow alpha range"
        HOD_range_min = [11., 0.001, 6., 12., 0.5]
        HOD_range_max = [12.2, 1., 14., 14., 1.5]
        theta_lbl = ['log $M_\mathrm{min}$', '$\sigma_{\mathrm{log}\,M}$', 'log $M_0$', 'log $M_1$', r'$\alpha$']
    elif HODrange == 'sinha2017prior':  
        range_descrip = "Sinha et al (2017) prior"
        HOD_range_min = [11., 0.001, 6., 12., 0.001]
        HOD_range_max = [12.2, 1., 14., 14., 2.]
        theta_lbl = ['log $M_\mathrm{min}$', '$\sigma_{\mathrm{log}\,M}$', 'log $M_0$', 'log $M_1$', r'$\alpha$']
    else: 
        raise NotImplementedError
    
    fname = ''.join([UT.dat_dir(), 'lhd/HOD_LHD.', HODrange, '.', method, '.', str(samples), 'samples.dat']) 
    
    if os.path.isfile(fname): # file exists 
        lhcube = np.loadtxt(fname, skiprows=4)
    else: 
        # Nsample x Ndim latin hypercube
        lhcube = thetaLHD([HOD_range_min, HOD_range_max], samples=samples, method=method)

        f = open(fname, 'w') 
        f.write('# '+range_descrip+'\n') 
        f.write('# parameters : '+', '.join(theta_lbl)+'\n') 
        f.write('# theta min : '+', '.join([str(r) for r in HOD_range_min])+'\n') 
        f.write('# theta max : '+', '.join([str(r) for r in HOD_range_max])+'\n') 

        for i in range(lhcube.shape[0]): 
            f.write('%f' % (lhcube[i,0]))
            for j in range(1, lhcube.shape[1]):
                f.write('\t %f' % (lhcube[i,j]))
            f.write('\n')
        f.close() 
    return lhcube 
    

def thetaLHD(theta_range, samples=None, method=None):  
    ''' Given the range of parameters, return a latin hypercube
    of the parameter space based on specified method 

    Parameters 
    ----------
    theta_range : list of two array
        list of two arrays. first specifying the lower bound  
        of the parameter space. The second specifying the upper bound
        of the parameter space. 

    method : str, 
        string specifying the LHD. Current implementation includes
        maximin, centermaximin, correlation, mdu, nohl. maximin and 
        centermaximin maxmizes the minimum distance between the 
        points. correlation minimizes the maximum correlation 
        coefficient. MDU (multidimensional uniformity) 

    samples : int, optional 
        sample size of the LH. NOHL has *fixed* sample size. If not
        specified, we use the fixed sample size of NOHL. 
    '''
    if len(theta_range) != 2: 
        raise ValueError('careful with the dimensions of theta_range input') 
    theta_min, theta_max = theta_range 
    ndim = len(theta_min) # dimensions of the parameter space 
    
    # check that the min and max values make sense 
    for idim in range(ndim): assert theta_min[idim] < theta_max[idim]
    
    # Nsample x Ndim latin hypercube
    lhcube = LHD(ndim, samples=samples, method=method) 

    # now scale to the parameter space
    for idim in range(ndim):
        lhcube[:,idim] *= theta_max[idim] - theta_min[idim]
        lhcube[:,idim] += theta_min[idim]
    return lhcube 


def LHD(dim, samples=None, method=None, niter=1000): 
    ''' Given dimensions and a string specifying the LHD method, 
    return latin hypercube 

    Parameters
    ----------
    dim : int, 
        dimenion of the LH

    method : str, 
        string specifying the LHD. Current implementation includes
        maximin, centermaximin, correlation, mdu, nohl. maximin and 
        centermaximin maxmizes the minimum distance between the 
        points. correlation minimizes the maximum correlation 
        coefficient. MDU (multidimensional uniformity) 

    samples : int, optional 
        sample size of the LH. NOHL has *fixed* sample size. If not
        specified, we use the fixed sample size of NOHL. 

    niter : int, default 1000
        number of iterations to run the pyDOE LHD algorithms. 
        Bigger the better! 

    Sources
    -------
    - https://pythonhosted.org/pyDOE/randomized.html#latin-hypercube 
    - Cioppa T. and Lucas T. (2007) 
    - Deutsch J. and Deutsch C. (2012) 
    '''
    if method is None: 
        raise ValueError("specify the LHD method") 
    m_list = ['maximin', 'centermaximin', 'correlation', 'mdu', 'nohl']
    if method not in m_list: 
        raise NotImplementedError("only methods ("+', '.join(m_list)+") implemented") 
    
    # sample size is not specified or NOHL 
    if samples is None or method == 'nohl': 
        try: 
            m, q, r = pyNOLH.params(dim)
            conf = random.sample(range(q), q)
            remove = random.sample(range(q), r)
            lhcube = pyNOLH.nolh(conf, remove)
        except IndexError: # this is bad
            conf, remove = pyNOLH.CONF[dim]
            lhcube = pyNOLH.nolh(conf, remove)
        samples = lhcube.shape[0]

    if method in ['maximin', 'centermaximin', 'correlation']: # pyDOE implementation 
        lhcube = pyDOE.lhs(dim, samples=samples, criterion=method, iterations=1000)
    elif method == 'mdu': # multidimensional uniformity 
        lhcube = np.array(lhsMDU.sample(dim, samples)).T
    return lhcube 
