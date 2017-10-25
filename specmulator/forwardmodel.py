'''


'''
import bigfile as BF 
import nbodykit.lab as NBlab
from nbodykit.lab import HODCatalog
from pmesh.pm import ParticleMesh
from scipy.interpolate import InterpolatedUnivariateSpline as interpolate


def Observables(gals, observable='plk'): 
    ''' Given galaxy catalog, measure specified observables 
    (e.g. powerspectrum multipoles). 
    '''
    if observable == 'plk': # powerspectrum multipole 
        # things are hardcoded below **for now**
        mesh = gals.to_mesh(window='tsc', Nmesh=256, compensated=True, position='RSDPosition')
        r = NBlab.FFTPower(mesh, mode='2d', dk=0.005, kmin=0.01, Nmu=5, los=[0,0,1], poles=[0,2,4])

        # think about the output some more! 
        poles = r.poles
        plk = {'k': poles['k']} 
        for ell in [0, 2, 4]:
            P = poles['power_%d' %ell].real
            if ell == 0: P = P - poles.attrs['shotnoise'] # add shotnoise to monopole 
            plk['p'+str(ell)+'k'] = P 
        return plk  


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

    # add RSD
    line_of_sight = [0,0,1]
    hod['RSDPosition'] = hod['Position'] + hod['VelocityOffset'] * line_of_sight
    
    hod.save('%s/HOD'%(folder), ('Position', 'Velocity', 'RSDPosition'))
    return None 


def RSDHalos(halos, LOS=[0,0,1]): 
    ''' Given halo catalog, return redshift space distorted halo catalog

    parameters
    ----------
    halos : 

    LOS : list, 3 elements
        list that specifies the line of sight 
    '''
    halo['RSDPosition'] = halo['Position'] + halo['VelocityOffset'] * LOS
    return halo 


def Halos(bs, nc, seed, nstep, seed_hod, Omega_m, p_alpha, p_logMin, p_logM1, p_logM0, p_sigma_logM):
    '''
    '''
    # setup initial conditions
    Omegacdm = Omega_m-0.049,
    cosmo   = NBlab.cosmology.Planck15.clone(Omega_cdm=Omegacdm, 
                                       h = 0.6711, Omega_b = 0.049)
    power   = NBlab.cosmology.LinearPower(cosmo, 0)
    klin    = np.logspace(-4,2,1000) 
    plin    = power(klin)
    pkfunc  = interpolate(klin, plin)
    
    # run the simulation
    pm      = ParticleMesh(BoxSize=bs, Nmesh=[nc, nc, nc])
    Q       = pm.generate_uniform_particle_grid()
    stages  = numpy.linspace(0.1, 1.0, nstep, endpoint=True)
    solver  = Solver(pm, cosmo, B=2)
    wn      = solver.whitenoise(seed)
    dlin    = solver.linear(wn, pkfunc)
    state   = solver.lpt(dlin, Q, stages[0])
    state   = solver.nbody(state, leapfrog(stages))
    
    # create the catalogue
    cat = ArrayCatalog({'Position' : state.X, 
    	  		'Velocity' : state.V, 
			'Displacement' : state.S,
			'Density' : state.RHO}, 
			BoxSize=pm.BoxSize, 
			Nmesh=pm.Nmesh,
			M0 = Omega_m * 27.75e10 * bs**3 / (nc/2.0)**3)
    cat['KDDensity'] = KDDensity(cat).density

    # run FOF to construct Halo catalog 
    fof = FOF(cat, linking_length=0.2, nmin=12)
    fofcat = fof.to_halos(particle_mass=cat.attrs['M0'], cosmo = cosmo, redshift = 0.0)      
    return fofcat  


