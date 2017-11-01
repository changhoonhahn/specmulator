'''
'''
import os
import numpy as np 
import nbodykit.lab as NBlab
from halotools.sim_manager import UserSuppliedHaloCatalog

import util as UT
import readsnap as RS
import readfof 

def NeutHalos(mneut, nreal, nzbin): 
    ''' Read halo catalogs generated by Paco
    
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
    '''
    if mneut == 0.1: 
        dir = ''.join([UT.dat_dir(), '0.10eV/', str(nreal)])
    else: 
        dir = ''.join([UT.dat_dir(), str(mneut), 'eV/', str(nreal)])
    # read in Gadget header
    header = RS.read_gadget_header(''.join([dir, '/snapdir_', str(nzbin).zfill(3), '/snap_', str(nzbin).zfill(3)]))
    
    # get cosmology from header 
    Omega_b = 0.049 # fixed baryon 
    cosmo = NBlab.cosmology.Planck15.clone(Omega_cdm=header['Omega_m']-Omega_b, h=header['h'], Omega_b=Omega_b)

    Fof = readfof.FoF_catalog(dir, nzbin, long_ids=False, swap=False, SFR=False)
    group_data = {}  
    group_data['Length']    = Fof.GroupLen
    group_data['Position']  = Fof.GroupPos/1e3
    group_data['Velocity']  = Fof.GroupVel
    group_data['Mass']      = Fof.GroupMass*1e10
    # calculate velocity offset
    rsd_factor = (1.+header['z']) / (100.*cosmo.efunc(header['z']))
    group_data['VelocityOffset']    = group_data['Velocity'] * rsd_factor
    #group_data['RSDPosition']       = group_data['Position'] + group_data['VelocityOffset'] * LOS
    
    # save to ArryCatalog for consistency
    cat = NBlab.ArrayCatalog(group_data, BoxSize=np.array([1000., 1000., 1000.])) 
    cat = NBlab.HaloCatalog(cat, cosmo=cosmo, redshift=header['z'], mdef='vir') 
    return cat


def NeutParticles(mneut, nreal, nzbin, clobber=False): 
    ''' Read particle catalog generated by Paco and return NBlab.ArrayCatalog

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

    clobber : bool, optional 
        if True, reconstructs the BigFile data 
    '''
    dir = ''.join([UT.dat_dir(), str(mneut), 'eV/', str(nreal), '/snapdir_', str(nzbin).zfill(3), '/'])
    dir_list = [dir+'/Matter'+sub for sub in ['', '/Position', '/Velocity', '/ID']]
    if (not np.all([os.path.isdir(dd) for dd in dir_list])) or clobber: 
        f = ''.join([dir, 'snap_', str(nzbin).zfill(3)]) 
        # read in Gadget header
        header = RS.read_gadget_header(f)

        # read in CDM particles (parttype = 1) and create catalogue
        particle_data = {} 
        particle_data['Position'] = RS.read_block(f, 'POS ', parttype=1)/1000. # Mpc/h
        particle_data['Velocity'] = RS.read_block(f, 'VEL ', parttype=1)
        particle_data['ID'] = RS.read_block(f, 'ID  ', parttype=1)
        cat = NBlab.ArrayCatalog(particle_data, BoxSize=np.array([header['boxsize'], header['boxsize'], header['boxsize']])) 
        #cat['KDDensity'] = KDDensity(cat).density
        cat.save(dir+'/Matter', ["Position", "Velocity", "ID"])
    else: 
        cat = NBlab.BigFileCatalog(dir+'/Matter', header='Header') 
    return cat 


def _NeutHalos(mneut, nreal, nzbin, clobber=False): 
    ''' Construct Friend-of-Friend halo catalog from Paco's neutrino 
    particle catalogs. 
    '''
    dir = ''.join([UT.dat_dir(), str(mneut), 'eV/', str(nreal), '/snapdir_', str(nzbin).zfill(3), '/'])
    f = ''.join([dir, 'snap_', str(nzbin).zfill(3)]) 
    # read in Gadget header
    header = RS.read_gadget_header(f)
    
    # get cosmology from header 
    Omega_b = 0.049 # fixed baryon 
    cosmo = NBlab.cosmology.Planck15.clone(Omega_cdm=header['Omega_m']-Omega_b, h=header['h'], Omega_b=Omega_b)
    
    if (not os.path.isdir(dir+'/FOF')) or clobber: 
        # DM particle mass (parttype=1)
        m_part = header['masses'][1]

        cat = NeutParticles(mneut, nreal, nzbin) # read in neutrino particles with 
        cat.attrs['Nmesh'] = [512, 512, 512]
        # calculate friend-of-friend with only CDM particles 
        fof = NBlab.FOF(cat, linking_length=0.2, nmin=20)
        # now make them into halos  
        fofcat = fof.to_halos(particle_mass=m_part, cosmo=cosmo, redshift=header['z'])      
        fofcat.save(dir+'/FOF', ['Position', 'Velocity', 'VelocityOffset', 'Mass', 'Radius'])
    else: 
        fofcat = NBlab.BigFileCatalog(dir+'/FOF', header='Header') 
    return fofcat 



"""
    dir = ''.join([UT.dat_dir(), str(mneut), 'eV/', str(nreal), '/groups_', str(nzbin).zfill(3), '/'])
    prefix = ''.join([dir, 'group_tab_', str(nzbin).zfill(3), '.']) 

    dt1 = np.dtype((np.float32,3))
    dt2 = np.dtype((np.float32,6))

    # initialize the data columns using the first file
    f = open(prefix+str(0), 'rb')
    NG         = np.fromfile(f, dtype=np.int32,  count=1)[0]
    TNG        = np.fromfile(f, dtype=np.int32,  count=1)[0]
    Nids       = np.fromfile(f, dtype=np.int32,  count=1)[0]
    TotNids    = np.fromfile(f, dtype=np.uint64, count=1)[0]
    Nfiles     = np.fromfile(f, dtype=np.uint32, count=1)[0]
    
    GroupLen    = np.empty(TNG, dtype=np.int32)
    GroupOffset = np.empty(TNG, dtype=np.int32)
    GroupMass   = np.empty(TNG, dtype=np.float32)
    GroupPos    = np.empty(TNG, dtype=dt1)
    GroupVel    = np.empty(TNG, dtype=dt1)
    GroupTLen   = np.empty(TNG, dtype=dt2)
    GroupTMass  = np.empty(TNG, dtype=dt2)
    
    # read in the rest of the snapshot group files
    istart = 0
    for i in range(0, Nfiles): 
        if i > 0: 
            f = open(prefix+str(i), 'rb')
            NG = np.fromfile(f, dtype=np.int32,  count=1)[0]
        if NG == 0: 
            continue 
        indx = slice(istart,istart+NG)
        GroupLen[indx]    = np.fromfile(f,dtype=np.int32,count=NG)
        GroupOffset[indx] = np.fromfile(f,dtype=np.int32,count=NG)
        GroupMass[indx]   = np.fromfile(f,dtype=np.float32,count=NG)
        GroupPos[indx]    = np.fromfile(f,dtype=dt1,count=NG)
        GroupVel[indx]    = np.fromfile(f,dtype=dt1,count=NG)
        GroupTLen[indx]   = np.fromfile(f,dtype=dt2,count=NG)
        GroupTMass[indx]  = np.fromfile(f,dtype=dt2,count=NG)
        istart += NG
        #curpos = f.tell()
        #f.seek(0,os.SEEK_END)
        #if curpos != f.tell(): print "Warning: finished reading before EOF for tab file",i
        f.close()

    # data columns
    sel = (GroupLen > 20)
    group_data = {}  
    group_data['Length'] = GroupLen[sel]
    group_data['Position'] = GroupPos[sel,:]/1e3
    group_data['Velocity'] = GroupVel[sel,:]
    group_data['Mass'] = GroupMass[sel]
    group_data['ID'] = np.arange(TNG)[sel]
        
    cat = NBlab.ArrayCatalog(group_data, BoxSize=np.array([1000., 1000., 1000.])) 
    return cat
"""

