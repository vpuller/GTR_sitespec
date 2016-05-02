# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 11:56:36 2016

@author: vpuller
"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import os

def mean_var(dist):
    dist_mean = dist.mean(axis=0)
    dist_var = (dist**2).mean(axis=0) - dist_mean**2 
    return dist_mean, dist_var
    
if __name__=="__main__":
    '''Generating sequences for a given tree using a GTR matrix'''
    
    plt.ioff()
    dir_name = '/ebio/ag-neher/share/users/vpuller/GTR_staggered/L100/' 
    outdir_name = dir_name + 'plots/'
    if not os.path.exists(outdir_name):
        os.makedirs(outdir_name)    
        
    with open(dir_name + 'seq0.txt','r') as seq0_file:
        seq0 = seq0_file.read()
#    with open(dir_name + 'alphabet.txt','r') as alpha_file:
#        alphabet = alpha_file.read()
    alphabet = 'ACGT'; q = len(alphabet); L = len(seq0)
    
    mmu = np.loadtxt(dir_name + 'mmu.txt')
    branch_lengths = np.loadtxt(dir_name + 'branch_lengths.txt')
    T = branch_lengths.sum()
    dist_to_root = np.loadtxt(dir_name + 'dist_to_root.txt')
    dist_to_root_mean = dist_to_root.mean()

    p0_a = np.loadtxt(dir_name + 'p0_a.txt')
    Wij0 = np.loadtxt(dir_name + 'Wij0.txt')

    data_names = ['aln','GTR','GTRanc','GTRtree', 'GTRtree_unique']
    Nsample = np.loadtxt(dir_name + 'Nsample.txt', dtype = 'int')
    p_a = {key: np.zeros((Nsample,mmu.shape[0], q, L)) for key in data_names}
    dist = {key: np.zeros((Nsample,mmu.shape[0])) for key in data_names}
    for jsample in xrange(Nsample):
        for jmu, mu in enumerate(mmu):
            for name in data_names:
                p_a[name][jsample,jmu,:,:] = np.loadtxt(dir_name + 'sample{}_mu{}_p_{}.txt'.format(jsample, jmu, name))
                dist[name][jsample,jmu] = np.sum((p_a[name][jsample,jmu,:,:] - p0_a)**2)  

    means = np.zeros((len(dist),len(mmu)))
    variances = np.zeros((len(dist),len(mmu)))
    for jname, name in enumerate(data_names):
        means[jname,:], variances[jname,:] = mean_var(dist[name])
    
    xx = mmu*dist_to_root_mean*Wij0[np.where(np.triu(Wij0) >0)].mean()
    leg = ('Tips alignment', 'GTR', 'GTR ancestors', 'GTR tree and ancestors')    
    plt.figure(30,figsize = (20,10)); plt.clf()
    plt.subplot(1,2,1)
    for jdata, mean_loc in enumerate(means):
        plt.errorbar(xx, np.log10(means[jdata,:]),\
        yerr = np.sqrt(variances[jdata,:])/means[jdata,:]/np.log(10))
    plt.xlabel(r'$\mu\bar{t}_{root}$', fontsize = 18)
    plt.ylabel(r'$\log_{10}\chi^2$', fontsize = 18)
    plt.legend(leg, fontsize = 18, loc = 0)

    plt.subplot(1,2,2)
    for jdata, mean_loc in enumerate(means):
        plt.errorbar(np.log10(xx), np.log10(means[jdata,:]),\
        yerr = np.sqrt(variances[jdata,:])/means[jdata,:]/np.log(10))
    plt.xlabel(r'$\log_{10} (\mu\bar{t}_{root})$', fontsize = 18)
    plt.ylabel(r'$\log_{10}\chi^2$', fontsize = 18)
    plt.legend(leg, loc = 0, fontsize = 18)
    plt.savefig(outdir_name + 'dist.pdf')
    plt.close(30)
    
    # studying reconstructed mutation rates
    mu_a = {}
    Wij = {}
    T_i = {}
    Tcons_i = {}
    for name in data_names[1:]:
        mu_arr = np.zeros((Nsample,mmu.shape[0],L))
        Wij_arr = np.zeros((Nsample,mmu.shape[0],q,q))
        T_i_arr = np.zeros((Nsample,mmu.shape[0], q, L))
        Tcons_i_arr = np.zeros((Nsample,mmu.shape[0], q, L))
        for jsample in xrange(Nsample):
            for jmu, mu in enumerate(mmu):
                mu_arr[jsample,jmu,:] = np.loadtxt(dir_name +\
                'sample{}_mu{}_mu_{}.txt'.format(jsample, jmu, name))
                Wij_arr[jsample,jmu,:,:] = np.loadtxt(dir_name +\
                'sample{}_mu{}_Wij_{}.txt'.format(jsample, jmu, name))
                T_i_arr[jsample,jmu,:,:] = np.mean(np.loadtxt(dir_name +\
                'sample{}_mu{}_T_i_{}.txt'.format(jsample,jmu, name)))
                Tcons_i_arr[jsample,jmu,:,:] = np.mean(np.loadtxt(dir_name +\
                'sample{}_mu{}_Tcons_i_{}.txt'.format(jsample,jmu, name)))
        mu_a[name] = mu_arr
        Wij[name] = Wij_arr
        T_i[name] = T_i_arr
        Tcons_i[name] = Tcons_i_arr

    dist_to_root_GTRtree = np.zeros((Nsample,mmu.shape[0]))
    for jsample in xrange(Nsample):
        for jmu, mu in enumerate(mmu):
            dist_to_root_GTRtree[jsample,jmu] = np.mean(np.loadtxt(dir_name +\
            'sample{}_mu{}_dist_to_root.txt'.format(jsample,jmu)))
    
    dist_to_root_GTRtree_unique = np.zeros((Nsample,mmu.shape[0]))
    for jsample in xrange(Nsample):
        for jmu, mu in enumerate(mmu):
            dist_to_root_GTRtree_unique[jsample,jmu] = np.mean(np.loadtxt(dir_name +\
            'sample{}_mu{}_dist_to_root_unique.txt'.format(jsample,jmu)))    
    
    factor = dist_to_root_mean*np.sum(Wij[name].mean(axis=(0,1)))/(q*(q-1))
    cols = ['g','r','c','m']
    plt.figure(30,figsize = (20,10)); plt.clf()
    for jname, name in enumerate([data_names[j] for j in [1,2,4]]):
        mu_mean, mu_var = mean_var((mu_a[name][0,:,:]*\
        T_i[name][0,:,:,:].sum(axis=1)/Tcons_i[name][0,:,:,:].sum(axis=1)).T)
        factor_tree = factor*T_i[name][0,:,:,0].sum(axis=1)/T_i['GTR'][0,:,:,0].sum(axis=1)
        if name == 'GTRtree':
            j0 = 10
            plt.errorbar(xx[j0:], mu_mean[j0:]*factor_tree, yerr = np.sqrt(mu_var[j0:])*factor_tree, color = cols[jname])
        elif name == 'GTRtree_unique':
            j0 = 10
            plt.errorbar(xx[j0:], mu_mean[j0:]*factor_tree, yerr = np.sqrt(mu_var[j0:])*factor_tree, color = cols[jname])
        else:
            plt.errorbar(xx, mu_mean*factor, yerr = np.sqrt(mu_var)*factor, color = cols[jname])
    plt.plot(xx,xx,':k')
    plt.xlabel(r'$\mu_0\bar{t}_{root}$', fontsize = 18)
    plt.ylabel(r'$\mu\bar{t}_{root}$', fontsize = 18)
    plt.legend(leg[1:], fontsize = 18, loc = 0)
    plt.savefig(outdir_name + 'mu_mu0.pdf')
    plt.close(30)
    
#    plt.figure(30,figsize = (20,10)); plt.clf()
#    for jname, name in enumerate([data_names[j] for j in [1,2,4]]):
#        mu_dist = ((mu_a['GTR'] - mmu[np.newaxis,:,np.newaxis])**2).sum(axis = 2)
#        mu_mean, mu_var = mean_var((mu_a[name][0,:,:]*\
#        T_i[name][0,:,:,:].sum(axis=1)/Tcons_i[name][0,:,:,:].sum(axis=1)).T)*\
#        T_i[name][0,:,:,0].sum(axis=1)/T_i['GTR'][0,:,:,0].sum(axis=1)
#       
#        if name == 'GTRtree':
#            j0 = 10
#            plt.errorbar(xx[j0:], mu_dist[j0:]*factor, yerr = np.sqrt(mu_var[j0:])*factor, color = cols[jname])
#        elif name == 'GTRtree_unique':
#            j0 = 10
#            plt.errorbar(xx[j0:], mu_dist[j0:]*factor, yerr = np.sqrt(mu_var[j0:])*factor, color = cols[jname])
#        else:
#            plt.errorbar(xx, mu_dist*factor, yerr = np.sqrt(mu_var)*factor, color = cols[jname])
#    plt.plot(xx,xx,':k')
#    plt.xlabel(r'$\mu_0\bar{t}_{root}$', fontsize = 18)
#    plt.ylabel(r'$\mu\bar{t}_{root}$', fontsize = 18)
#    plt.legend(leg[1:], fontsize = 18, loc = 0)
#    plt.savefig(outdir_name + 'mu_dist.pdf')
#    plt.close(30)
    
#    plt.figure(40, figsize = (24, 30)); plt.clf()
#    for jmu, mu in enumerate(mmu):
#        plt.subplot(5,4,jmu+1)
#        plt.hist(mu_a['GTR'][0,jmu,:]/mu)
#        plt.xlabel(r'$\mu$')
#        plt.title('mu = {}'.format(xx[jmu]))
#    plt.savefig(outdir_name + 'hist.pdf')
#    plt.close(40)
    
    