'''

emulator 

'''
import numpy as np 

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

    def read_HODLHD(self, HODrange='sinha2017prior_narrow', method='nohl', samples=17): 
        ''' read in Halo Occupation Distribution (HOD) parameters determined by 
        Latin Hypercube Design (LHD). 
        '''
        self.HODrange = HODrange
        self.LHDmethod = method 
        self.LHDsamples = samples 
        # save HODLHD info 
        lhcube = lhd.HOD_LHD(HODrange=HODrange, samples=samples, method=method)
        self.lhcube = lhcube 
        return lhcube 
    
    def read_NeutObvs(self, obvs, mneut, nreal, nzbin, seed_hod, Nmesh=360, rsd=True): 
        ''' Read in observables( theta_HODLHD ) 
        '''
        self.mneut = mneut 
        self.nreal = nreal 
        self.nzbin = nzbin
        if obvs == 'plk': 
            self.obvs = HODLHD_Pks(mneut, nreal, nzbin, seed_hod, 
                    HODrange=self.HODrange, method=self.LHDmethod, samples=self.LHDsamples, 
                    Nmesh=Nmesh, rsd=rsd)
        else: 
            raise NotImplementedError
        return self.obvs
     
    def trainGP(self):   
        '''
        '''
        lguess = (np.max(self.lhcube, 0) - np.min(self.lhcube, 0))/float(self.LHDsamples)
        print lguess
        kernel = 1.0 * Kerns.ExpKernel(metric=lguess, ndim=self.lhcube.shape[1])
        gp = George.GP(kernel) 
        gp.optimize(self.lhcube, pks)
        self.GP = gp 
        return gp


def HODLHD_Pks(mneut, nreal, nzbin, seed_hod, 
        HODrange='sinha2017prior_narrow', method='nohl', samples=17, Nmesh=360, rsd=True): 
    '''
    '''
    plks = [] 
    for i_p in range(samples): 
        plk_i = lhd.HODLHD_NeutObvs('plk', mneut, nreal, nzbin, seed_hod, i_p,
                HODrange=HODrange, method=method, samples=samples, Nmesh=Nmesh, rsd=rsd)
        plks.append(plk_i['p0k'][8])
    return np.array(plks) 


if __name__=="__main__":
    HODemulator(HODrange='sinha2017prior_narrow', samples=17, method='mdu')
    
