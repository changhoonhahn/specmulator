'''


'''
import bigfile as BF 
from nbodykit.lab import HODCatalog
#from nbodykit.lab import HaloCatalog


def Galaxies(halos, p_hod, seed=None): 
    ''' populate given halo catalog (halos) with galaxies
    based on HOD model with p_hod parameters 

    Parameters
    ----------
    halos : halotools.sim_manager.UserSuppliedHaloCatalog
        output of HaloCatalog.to_halotools. for our intents and purposes 
        this is the halo catalog
    p_hod : dict
        dictionary specifying the HOD parameters 
    '''
    # check halos 
    if 'alpha' not in p_hod.keys(): 
        raise ValueError
    if 'logMmin' not in p_hod.keys(): 
        raise ValueError
    if 'logM1' not in p_hod.keys(): 
        raise ValueError
    if 'logM0' not in p_hod.keys(): 
        raise ValueError
    if 'sigma_logM' not in p_hod.keys(): 
        raise ValueError

    # run HOD
    halocat = halos.to_halotools(halos.attrs['BoxSize'])
    hod = HODCatalog(halocat, seed=seed, **params)
    
    hod.save('%s/HOD'%(folder), ('Position', 'Velocity'))
    return None 
   

def Halos(): 
    '''
    '''
