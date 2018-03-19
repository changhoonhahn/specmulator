'''

emulator 

'''
import numpy as np 
import scipy.optimize as op

import george as George
from george import kernels as Kerns

import lhd 


class HODemulator(object): 
    def __init__(self):
        '''
        '''
        self.HODrange = None 
        self.LHDmethod = None
        self.LHDsamples = None
        
        self.mneut = None 
        self.nreal = None 
        self.nzbin = None
    
    def trainGP(self, X, Y, yerr=None):   
        ''' Train Gaussian Process for each component of observable Y 

        Parameters
        ----------
        X : np.ndarray (Nsample x Ndim) 
            parameter values where Y is evaluated  

        Y : np.ndarray (Nsample x Ncomponents)  
            Observable at parameter value X 
        '''
        if X.shape[0] != Y.shape[0]: 
            raise ValueError
        self._x = X
        
        lguess = [(np.max(X[:,i]) - np.min(X[:,i]))/float(X.shape[0]) for i in range(X.shape[1])]

        self.GPs = [] 
        self.kernels = [] 
        for i in range(Y.shape[1]): # optimize a GP for each component 
            self._y = Y[:,i]
            
            kernel = np.var(self._y) * Kerns.ExpSquaredKernel(lguess, ndim=len(lguess), axes=range(X.shape[1])) #
            _gp = George.GP(kernel, mean=np.average(self._y, axis=0)) 
            if yerr is None: 
                _gp.compute(self._x)
            else: 
                _gp.compute(self._x, yerr=yerr[i])
            # optimize hyperparameters of the kernel
            p0 = _gp.get_parameter_vector()
            results = op.minimize(self._nll, p0, args=(_gp), jac=self._grad_nll, method="L-BFGS-B")
            _gp.set_parameter_vector(results.x)
            self.GPs.append(_gp)
            self.kernels.append(_gp.kernel) 
        return self.GPs 

    def kernels(self): 
        ''' select kernel for gaussian process. Unclear ATM how this should
        be structured.
        '''
        return None 

    def read_HODLHD(self, HODrange='sinha2017prior_narrow', method='nohl', samples=17): 
        ''' read in Halo Occupation Distribution (HOD) parameters determined by 
        Latin Hypercube Design (LHD). 
        '''
        self.HODrange = HODrange
        self.LHDmethod = method 
        self.LHDsamples = samples 
    
        # ***this is hardcoded for convenience*** 
        # needs to be consistent with lhd.HOD_LHD 
        if self.HODrange == 'sinha2017prior_narrow': # narrow alpha range 
            range_descrip = "Sinha et al (2017) prior with narrow alpha range"
            self.HODrange_min = [11., 0.001, 6., 12., 0.5]
            self.HODrange_max = [12.2, 1., 14., 14., 1.5]
        elif self.HODrange == 'sinha2017prior':  
            range_descrip = "Sinha et al (2017) prior"
            self.HODrange_min = [11., 0.001, 6., 12., 0.001]
            self.HODrange_max = [12.2, 1., 14., 14., 2.]
        self.HOD_params = ['logMmin', 'sigma_logM', 'logM0', 'logM1', 'alpha'] 
        self.HOD_labels = ['log $M_\mathrm{min}$', '$\sigma_{\mathrm{log}\,M}$', 'log $M_0$', 'log $M_1$', r'$\alpha$']

        # save HODLHD info 
        lhcube = lhd.HOD_LHD(HODrange=HODrange, samples=samples, method=method)
        self.lhcube = lhcube 
        return lhcube 
    
    def read_NeutObvs(self, obvs, mneut, nreal, nzbin, seed_hod, Nmesh=360, rsd=True, 
            krange=[0.01, 0.5], silent=False): 
        ''' Read in observables( theta_HODLHD ) 
        ** need more descriptions ** 
        '''
        if np.any([v is None for v in [self.HODrange, self.LHDmethod, self.LHDsamples]]): 
            raise ValueError("need to run self.read_HODLHD beforehand") 
        self.mneut = mneut 
        self.nreal = nreal 
        self.nzbin = nzbin
        if obvs == 'p0k':  # powerspectrum monopole 
            self.obvs = self._HODLHD_Pks(0, mneut, nreal, nzbin, seed_hod, 
                    HODrange=self.HODrange, method=self.LHDmethod, samples=self.LHDsamples, 
                    Nmesh=Nmesh, rsd=rsd, krange=krange, silent=silent)
        else: 
            raise NotImplementedError
        return self.obvs
     
    def _HODLHD_Pks(self, ell, mneut, nreal, nzbin, seed_hod, 
            HODrange='sinha2017prior_narrow', method='nohl', samples=17, Nmesh=360, rsd=True, 
            krange=None, silent=False): 
        ''' Read in powerspectrum measurements for HOD LHD 
        '''
        if ell not in [0,2,4]: 
            raise ValueError("only monopole, quadrupole, and hexadecapole") 
        plks = [] 
        for i_p in range(samples): 
            try: 
                plk_i = lhd.HODLHD_NeutObvs('plk', mneut, nreal, nzbin, seed_hod, i_p,
                        HODrange=HODrange, method=method, samples=samples, Nmesh=Nmesh, 
                        rsd=rsd, make=False, silent=silent)
            except ValueError: 
                continue 
            if krange is not None: 
                klim = np.where((plk_i['k'] > krange[0]) & (plk_i['k'] < krange[1]))
                plks.append(plk_i['p'+str(ell)+'k'][klim])
            else:
                plks.append(plk_i['p'+str(ell)+'k'])
        if krange is not None: 
            self.k = plk_i['k'][klim]
        else:
            self.k = plk_i['k']
        return np.array(plks) 
    
    def _nll(self, p, _gp) :
        # negative log likliehood 
        _gp.set_parameter_vector(p)
        ll = _gp.lnlikelihood(self._y, quiet=True)
        # The scipy optimizer doesn't play well with infinities.
        return -ll if np.isfinite(ll) else 1e25

    def _grad_nll(self, p, _gp):
        # grad log likelihood
        _gp.set_parameter_vector(p)
        return -_gp.grad_lnlikelihood(self._y, quiet=True)
