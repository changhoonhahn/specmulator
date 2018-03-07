import numpy as np 
import env
import util as UT 
import data as Dat
import emulator as Emu
import compress as Comp
import forwardmodel as FM 

from ChangTools.plotting import prettycolors
import matplotlib as mpl 
import matplotlib.pyplot as plt
mpl.rcParams['text.usetex'] = True
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['axes.linewidth'] = 1.5
mpl.rcParams['axes.xmargin'] = 1
mpl.rcParams['xtick.labelsize'] = 'x-large'
mpl.rcParams['xtick.major.size'] = 5
mpl.rcParams['xtick.major.width'] = 1.5
mpl.rcParams['ytick.labelsize'] = 'x-large'
mpl.rcParams['ytick.major.size'] = 5
mpl.rcParams['ytick.major.width'] = 1.5
mpl.rcParams['legend.frameon'] = False



def HODemulator_PkPCA(mneut=0.0, nreal=1, nzbin=4, seed_hod=1, Nmesh=360, rsd=True, 
        HODrange='sinha2017prior_narrow', method='mdu', samples=40, krange=[0.01,0.5]):
    ''' Test simple gaussian process emulator by randomly 
    choosing an HOD parameter then comparing it to the observables
    of HOD parameters near it 
    '''
    emu = Emu.HODemulator() 
    lhcube = emu.read_HODLHD(HODrange=HODrange, method=method, samples=samples)
    plks = emu.read_NeutObvs('p0k', mneut, nreal, nzbin, seed_hod, Nmesh=Nmesh, rsd=rsd)

    # compress data with PCA 
    n_pca = 21  # number of PCA components 
    pca_fid = Comp.PCA_fid(n_pca, mneut=0.0, nreal=nreal, nzbin=nzbin, obvs='plk', poles=[0], 
            Nmesh=Nmesh, rsd=rsd, HODrange=HODrange, krange=krange)
    plks_pca = pca_fid.transform(plks)
    
    X_fid = Comp.X_fid(nreal, nzbin, obvs='plk', poles=[0], mneut=0.0, Nmesh=Nmesh, rsd=rsd, HODrange=HODrange, krange=krange)
    X_fid_pca = pca_fid.transform(X_fid)
    yerr = np.diag(np.cov(X_fid_pca.T))
    emu.trainGP(lhcube, plks_pca, yerr=yerr) # train GP 
    n_comp = len(emu.GPs) # number of GPs 
    print("%i GPs" % n_comp)

    # randomly select one theta in HOD parameter space 
    p_hod = np.zeros(len(emu.HOD_params))
    for ik, k in enumerate(emu.HOD_params): 
        p_hod[ik] = 0.5*(emu.HODrange_max[ik] + emu.HODrange_min[ik]) #+ \ np.random.uniform(-0.01, 0.33) * (emu.HODrange_max[ik] - emu.HODrange_min[ik])
    
    p0ks_gp_pca = np.zeros(plks_pca.shape[1])
    vars_gp_pca = np.zeros(n_comp)
    for i in range(n_comp): 
        mu, var = emu.GPs[i].predict(plks_pca[:,i], np.array((p_hod,)))
        p0ks_gp_pca[i] = mu
        vars_gp_pca[i] = np.diag(var)
    
    # select parameters in lhcube "closest" to 
    # the chosen parameter above
    rho_hod = np.sqrt(np.sum(((lhcube - p_hod)/(np.array(emu.HODrange_max) - np.array(emu.HODrange_min)))**2, axis=0))
    i_rhosort = np.argsort(rho_hod)
    i_NN = i_rhosort[:5]  # 5 nearest HOD LHD points 
    
    fig = plt.figure(figsize=(10,5)) 
    sub = fig.add_subplot(121)
    sub.plot(range(n_comp), p0ks_gp_pca, c='k', ls='--')  
    sub.fill_between(range(n_comp), p0ks_gp_pca - np.sqrt(vars_gp_pca), p0ks_gp_pca + np.sqrt(vars_gp_pca),linewidth=0, color='k', alpha=0.25) 
    for i in i_NN: 
        sub.plot(range(n_comp), plks_pca[i,:]) 
    # x-axis
    sub.set_xlim([0, n_comp]) 
    # y-axis
    sub.set_ylabel('PCA of $P(k)$', fontsize=25)

    sub = fig.add_subplot(122)
    sub.scatter(lhcube[:,0], lhcube[:,1], c='k')
    sub.scatter([lhcube[i,0] for i in i_NN], [lhcube[i,1] for i in i_NN], c='b')
    sub.scatter([p_hod[0]], [p_hod[1]], c='r') 
    # x-axis
    sub.set_xlim([emu.HODrange_min[0], emu.HODrange_max[0]]) 
    sub.set_xlabel(emu.HOD_labels[0], fontsize=20)
    # y-axis
    sub.set_ylim([emu.HODrange_min[1], emu.HODrange_max[1]]) 
    sub.set_ylabel(emu.HOD_labels[1], fontsize=20)
    f = ''.join([UT.fig_dir(), 'tests/test.HODemulator.trainGP_PCA_NN.png']) 
    fig.savefig(f)#, bbox_inches='tight') 
    return None 


def HODemulator_PkNN_compress(mneut=0.0, nreal=1, nzbin=4, seed_hod=1, Nmesh=360, rsd=True, 
        HODrange='sinha2017prior_narrow', method='mdu', samples=40, krange=[0.01,0.5]):
    ''' Test simple gaussian process emulator by randomly 
    choosing an HOD parameter then comparing it to the observables
    of HOD parameters near it 
    '''
    emu = Emu.HODemulator() 
    lhcube = emu.read_HODLHD(HODrange=HODrange, method=method, samples=samples)
    plks = emu.read_NeutObvs('p0k', mneut, nreal, nzbin, seed_hod, Nmesh=Nmesh, rsd=rsd)

    # compress data with PCA 
    n_pca = plks.shape[1]  # number of PCA components 
    pca_fid = Comp.PCA_fid(n_pca, mneut=0.0, nreal=nreal, nzbin=nzbin, obvs='plk', poles=[0], 
            Nmesh=Nmesh, rsd=rsd, HODrange=HODrange, krange=krange)
    plks_pca = pca_fid.transform(plks)
    
    X_fid = Comp.X_fid(nreal, nzbin, obvs='plk', poles=[0], mneut=0.0, Nmesh=Nmesh, rsd=rsd, HODrange=HODrange, krange=krange)
    X_fid_pca = pca_fid.transform(X_fid)
    yerr = np.diag(np.cov(X_fid_pca.T))
    emu.trainGP(lhcube, plks_pca)#, yerr=yerr) # train GP 
    n_comp = len(emu.GPs) # number of GPs 
    print("%i GPs" % n_comp)

    # randomly select one theta in HOD parameter space 
    p_hod = np.zeros(len(emu.HOD_params))
    for ik, k in enumerate(emu.HOD_params): 
        p_hod[ik] = 0.5*(emu.HODrange_max[ik] + emu.HODrange_min[ik]) #+ \ np.random.uniform(-0.01, 0.33) * (emu.HODrange_max[ik] - emu.HODrange_min[ik])
    
    p0ks_gp_pca = np.zeros(plks_pca.shape[1])
    vars_gp_pca = np.zeros(n_comp)
    for i in range(n_comp): 
        mu, var = emu.GPs[i].predict(plks_pca[:,i], np.array((p_hod,)))
        p0ks_gp_pca[i] = mu
        vars_gp_pca[i] = np.diag(var)
    
    p0ks_gp = pca_fid.inverse_transform(p0ks_gp_pca)

    # select parameters in lhcube "closest" to 
    # the chosen parameter above
    rho_hod = np.sqrt(np.sum(((lhcube - p_hod)/(np.array(emu.HODrange_max) - np.array(emu.HODrange_min)))**2, axis=0))
    i_rhosort = np.argsort(rho_hod)
    i_NN = i_rhosort[:5]  # 5 nearest HOD LHD points 
    
    fig = plt.figure(figsize=(10,5)) 
    sub = fig.add_subplot(121)
    sub.plot(emu.k, p0ks_gp, c='k', ls='--')  
    sub.fill_between(emu.k, 
            pca_fid.inverse_transform(p0ks_gp_pca - np.sqrt(vars_gp_pca)), 
            pca_fid.inverse_transform(p0ks_gp_pca + np.sqrt(vars_gp_pca)), 
            linewidth=0, color='k', alpha=0.25) 
    for i in i_NN: 
        sub.plot(emu.k, pca_fid.inverse_transform(plks_pca[i,:]))
    # x-axis
    sub.set_xscale('log') 
    sub.set_xlim([0.01, 0.5]) 
    sub.set_xlabel('k', fontsize=25)
    # y-axis
    sub.set_yscale('log') 
    sub.set_ylabel('$P(k)$', fontsize=25)

    sub = fig.add_subplot(122)
    sub.scatter(lhcube[:,0], lhcube[:,1], c='k')
    sub.scatter([lhcube[i,0] for i in i_NN], [lhcube[i,1] for i in i_NN], c='b')
    sub.scatter([p_hod[0]], [p_hod[1]], c='r') 
    # x-axis
    sub.set_xlim([emu.HODrange_min[0], emu.HODrange_max[0]]) 
    sub.set_xlabel(emu.HOD_labels[0], fontsize=20)
    # y-axis
    sub.set_ylim([emu.HODrange_min[1], emu.HODrange_max[1]]) 
    sub.set_ylabel(emu.HOD_labels[1], fontsize=20)
    f = ''.join([UT.fig_dir(), 'tests/test.HODemulator.trainGP_PkNN.compressed.png']) 
    fig.savefig(f)#, bbox_inches='tight') 
    return None 


def HODemulator_PkNN(mneut=0.0, nreal=1, nzbin=4, seed_hod=1, Nmesh=360, rsd=True, 
        HODrange='sinha2017prior_narrow', method='mdu', samples=40, krange=[0.01,0.5]):
    ''' Test simple gaussian process emulator by randomly 
    choosing an HOD parameter then comparing it to the observables
    of HOD parameters near it 
    '''
    emu = Emu.HODemulator() 
    lhcube = emu.read_HODLHD(HODrange=HODrange, method=method, samples=samples)
    plks = emu.read_NeutObvs('p0k', mneut, nreal, nzbin, seed_hod, 
            Nmesh=Nmesh, rsd=rsd, krange=krange)
    X_fid = Comp.X_fid(nreal, nzbin, obvs='plk', poles=[0], mneut=mneut, 
            Nmesh=Nmesh, rsd=rsd, HODrange=HODrange, krange=krange)
    yerr = np.diag(np.cov(X_fid.T))
    emu.trainGP(lhcube, plks, yerr=yerr) # train GP 
    n_comp = len(emu.GPs) # number of GPs 
    print("%i GPs" % n_comp)

    # randomly select one theta in HOD parameter space 
    p_hod = np.zeros(len(emu.HOD_params))
    for ik, k in enumerate(emu.HOD_params): 
        p_hod[ik] = 0.5*(emu.HODrange_max[ik] + emu.HODrange_min[ik]) #+ \ np.random.uniform(-0.01, 0.33) * (emu.HODrange_max[ik] - emu.HODrange_min[ik])
    
    p0ks_gp = np.zeros(n_comp)
    vars_gp = np.zeros(n_comp)
    for i in range(n_comp): 
        mu, var = emu.GPs[i].predict(plks[:,i], np.array((p_hod,)))
        p0ks_gp[i] = mu
        vars_gp[i] = np.diag(var)

    # select parameters in lhcube "closest" to 
    # the chosen parameter above
    rho_hod = np.sqrt(np.sum(((lhcube - p_hod)/(np.array(emu.HODrange_max) - np.array(emu.HODrange_min)))**2, axis=0))
    i_rhosort = np.argsort(rho_hod)
    i_NN = i_rhosort[:5]  # 5 nearest HOD LHD points 
    
    fig = plt.figure(figsize=(10,5)) 
    sub1 = fig.add_subplot(121)
    sub1.plot(emu.k, p0ks_gp, c='k', ls='--')  
    #sub.fill_between(emu.k, 
    #        p0ks_gp - np.sqrt(vars_gp), 
    #        p0ks_gp + np.sqrt(vars_gp), 
    #        linewidth=0, color='k', alpha=0.25) 
    # x-axis
    sub1.set_xscale('log') 
    sub1.set_xlim([0.01, 0.5]) 
    sub1.set_xlabel('k', fontsize=25)
    # y-axis
    sub1.set_yscale('log') 
    sub1.set_ylabel('$P(k)$', fontsize=25)

    sub2 = fig.add_subplot(122)
    sub2.scatter(lhcube[:,0], lhcube[:,1], c='k')
    for ii, i in enumerate(i_NN): 
        sub1.plot(emu.k, plks[i,:], c='C'+str(ii))
        sub2.scatter([lhcube[i,0]], [lhcube[i,1]], c='C'+str(ii))
    sub2.scatter([p_hod[0]], [p_hod[1]], c='r') 
    # x-axis
    sub2.set_xlim([emu.HODrange_min[0], emu.HODrange_max[0]]) 
    sub2.set_xlabel(emu.HOD_labels[0], fontsize=20)
    # y-axis
    sub2.set_ylim([emu.HODrange_min[1], emu.HODrange_max[1]]) 
    sub2.set_ylabel(emu.HOD_labels[1], fontsize=20)
    f = ''.join([UT.fig_dir(), 'tests/test.HODemulator.trainGP_PkNN.png']) 
    fig.savefig(f)#, bbox_inches='tight') 
    return None 


def HODemulator_trainGP_Pk(mneut=0.0, nreal=1, nzbin=4, seed_hod=1, Nmesh=360, rsd=True, 
        HODrange='sinha2017prior_narrow', method='mdu', samples=40, krange=[0.01, 0.5]):
    ''' Test simple gaussian process emulator by randomly 
    sampling HOD parameter space and compare it to GP predictions 
    '''
    emu = Emu.HODemulator() 
    lhcube = emu.read_HODLHD(HODrange=HODrange, method=method, samples=samples)
    plks = emu.read_NeutObvs('p0k', mneut, nreal, nzbin, seed_hod, 
            Nmesh=Nmesh, rsd=rsd, krange=krange, silent=True)
    mu_plks = np.sum(plks, axis=0)/float(plks.shape[0])
    dplks = plks - mu_plks

    X_fid = Comp.X_fid(nreal, nzbin, obvs='plk', poles=[0], mneut=0.0, 
            Nmesh=Nmesh, rsd=rsd, HODrange=HODrange, krange=krange, silent=True)
    yerr = np.diag(np.cov(X_fid.T))
    emu.trainGP(lhcube, dplks, yerr=yerr) # train GP 
    n_comp = len(emu.GPs) # number of GPs 
    print("%i GPs" % n_comp)

    # randomly sample HOD parameter space 
    p_hods = np.zeros((5, len(emu.HOD_params)))
    for i in range(p_hods.shape[0]): 
        for ik, k in enumerate(emu.HOD_params): 
            p_hods[i,ik] = 0.5*(emu.HODrange_max[ik] + emu.HODrange_min[ik]) + \
                    np.random.uniform(-0.1, 0.1) * (emu.HODrange_max[ik] - emu.HODrange_min[ik])

    p0ks_gp = np.zeros((p_hods.shape[0], n_comp))
    var_gp = np.zeros((p_hods.shape[0], n_comp))
    for i in range(n_comp): 
        mu, var = emu.GPs[i].predict(plks[:,i], p_hods)
        p0ks_gp[:,i] = mu + mu_plks[i]
        var_gp[:,i] = np.diag(var)
    
    p0ks_mock = []
    halos = Dat.NeutHalos(mneut, nreal, nzbin) 
    for i in range(p_hods.shape[0]): 
        # forward model at p_hods[i] 
        p_hod_i = {} 
        for t, k in zip(p_hods[i,:], emu.HOD_params): p_hod_i[k] = t 
        gals = FM.Galaxies(halos, p_hod_i, seed=seed_hod)  # populate the halo catalogs using HOD 
        gals['RSDPosition'] = FM.RSD(gals, LOS=[0,0,1]) # RSD position (hardcoded in the z direction) 
        plk = FM.Observables(gals, observable='plk', rsd=rsd, Nmesh=Nmesh, krange=krange)
        p0ks_mock.append(plk['p0k']) 
    k_mock = plk['k'] 

    # compare residual between the forwardmodeled P(k) to the Gaussian process
    # with the variance predicted by GP
    fig = plt.figure(figsize=(8,4)) 
    sub = fig.add_subplot(121)
    pretty_colors = prettycolors() 
    #sub.fill_between(emu.k, -np.sqrt(var_gp[i,:])/p0ks_mock[i], np.sqrt(var_gp[i,:])/p0ks_mock[i], color=pretty_colors[i], edgecolor="none", alpha=0.5)  
    # x-axis
    sub.set_xscale('log') 
    sub.set_xlim(krange) 
    sub.set_xlabel('k', fontsize=25)
    # y-axis
    sub.set_ylim([-1., 1.])
    #sub.set_yscale('log') 
    sub.set_ylabel('$\Delta P(k)/P(k)$', fontsize=25)

    sub2 = fig.add_subplot(122)
    sub2.scatter(lhcube[:,0], lhcube[:,1], c='k')
    for i in range(p_hods.shape[0]): 
        sub.plot(emu.k, (p0ks_gp[i,:]-p0ks_mock[i])/p0ks_mock[i], ls='--', c='C'+str(i)) 
        sub2.scatter([p_hods[i,0]], [p_hods[i,1]], c='C'+str(i))
    # x-axis
    sub2.set_xlim([emu.HODrange_min[0], emu.HODrange_max[0]]) 
    sub2.set_xlabel(emu.HOD_labels[0], fontsize=20)
    # y-axis
    sub2.set_ylim([emu.HODrange_min[1], emu.HODrange_max[1]]) 
    sub2.set_ylabel(emu.HOD_labels[1], fontsize=20)

    f = ''.join([UT.fig_dir(), 'tests/test.HODemulator.trainGP_Pk.png']) 
    fig.savefig(f, bbox_inches='tight') 
    return None 


def HODemulator_readNeutObvs(obvs='p0k', mneut=0.0, nreal=1, nzbin=4, seed_hod=1, Nmesh=360, rsd=True, 
        HODrange='sinha2017prior_narrow', method='mdu', samples=17): 
    ''' ***TESTED Nov 13, 2017***
    Test HODemulator.read_HODLHD
    '''
    emu = Emu.HODemulator() 
    lhcube = emu.read_HODLHD(HODrange=HODrange, method=method, samples=samples)
    plks = emu.read_NeutObvs(obvs, mneut, nreal, nzbin, seed_hod, Nmesh=Nmesh, rsd=rsd)
    fig = plt.figure() 
    sub = fig.add_subplot(111)
    for i in range(plks.shape[0]): 
        sub.plot(emu.k, plks[i,:]) 
    # x-axis
    sub.set_xscale('log') 
    sub.set_xlim([0.01, 0.5]) 
    sub.set_xlabel('k', fontsize=25)
    # y-axis
    sub.set_yscale('log') 
    sub.set_ylabel('$P(k)$', fontsize=25)
    f = ''.join([UT.fig_dir(), 'tests/test.HODemulator.', obvs, '.png']) 
    fig.savefig(f, bbox_inches='tight') 
    return None 


def HODemulator_readHODLHD(HODrange='sinha2017prior_narrow', method='mdu', samples=17): 
    ''' ***TESTED Nov 13, 2017*** 
    Test HODemulator.read_HODLHD
    '''
    emu = Emu.HODemulator() 
    lhcube = emu.read_HODLHD(HODrange=HODrange, method=method, samples=samples)
    fig = plt.figure(figsize=(9,9))
    for i in range(lhcube.shape[1]):
        for j in range(lhcube.shape[1]): 
            if i < j:
                sub = fig.add_subplot(lhcube.shape[1],lhcube.shape[1],lhcube.shape[1]*j+i+1)
                sub.scatter(lhcube[:,i], lhcube[:,j])
                sub.set_xlim([emu.HODrange_min[i], emu.HODrange_max[i]])
                if j == lhcube.shape[1]-1: 
                    sub.set_xlabel(emu.HOD_labels[i], fontsize=20)
                else: 
                    sub.set_xticks([])
                sub.set_ylim([emu.HODrange_min[j], emu.HODrange_max[j]])
                if i == 0: 
                    sub.set_ylabel(emu.HOD_labels[j], fontsize=20)
                else: 
                    sub.set_yticks([])            
            elif i == j:
                sub = fig.add_subplot(lhcube.shape[1],lhcube.shape[1],lhcube.shape[1]*j+i+1)
                sub.hist(lhcube[:,i], range=[emu.HODrange_min[i], emu.HODrange_max[i]], normed=True)
                sub.set_xlim([emu.HODrange_min[i], emu.HODrange_max[i]])
                if i != 0: 
                    sub.set_yticks([])            
                if i != lhcube.shape[1]-1: 
                    sub.set_xticks([])            
    f = ''.join([UT.fig_dir(), 'tests/test.HODemulator.LHD.', method, '.', str(samples), '.samples.png']) 
    fig.savefig(f, bbox_inches='tight') 
    return None 


if __name__=="__main__": 
    #HODemulator_PkNN_compress()
    #HODemulator_trainGP_Pk()
    #HODemulator_PkPCA()
    #HODemulator_PkNN()
    HODemulator_trainGP_Pk()
