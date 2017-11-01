'''

Code for the Latin Hypercube Design 


'''
import os 
import random
import numpy as np 

import util as UT

import pynolh as pyNOLH
import pyDOE
import lhsmdu as lhsMDU


def HOD_LHD(HODrange=None, samples=None, method=None):
    ''' Return latin hypercubes for vanilla Zheng et al.(2007) 
    HOD model. 

    Sources
    -------
    - Zheng Z., et al. (2007) -- arXiv:0512071
    - Sinha M., et al. (2017) -- arXiv:1708.04892 
    '''
    if HODrange is 'sinha2017prior':  
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
