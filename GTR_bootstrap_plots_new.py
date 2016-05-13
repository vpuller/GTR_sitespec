# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 11:56:36 2016

@author: vpuller
"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import os

# Constants
h = 10**(-8)

def mean_var(dist):
    dist_mean = dist.mean(axis=0)
    dist_var = (dist**2).mean(axis=0) - dist_mean**2 
    return dist_mean, dist_var
        
class model_param_data(object):
    '''Loading model parametes and data'''
    def __init__(self,dir_name, model_name, data_names):
        '''loading model, simulation parameters and simulation data'''
        model_file_head = dir_name + model_name + '/'
        with open(model_file_head + 'seq0.txt','r') as seq0_file:
            seq0 = seq0_file.read()
        with open(model_file_head + 'alphabet.txt','r') as alpha_file:
            alphabet = alpha_file.read()
        self.q = len(alphabet); self.L = len(seq0)
        
        self.mmu = np.loadtxt(model_file_head + 'mmu.txt')
        self.branch_lengths = np.loadtxt(model_file_head + 'branch_lengths.txt')
        self.T = self.branch_lengths.sum()
        self.dist_to_root = np.loadtxt(model_file_head + 'dist_to_root.txt')
        self.dist_to_root_mean = self.dist_to_root.mean()
    
        self.p0_a = np.loadtxt(model_file_head + 'p0_a.txt')
        self.Wij0 = np.loadtxt(model_file_head + 'Wij0.txt')
        self.Nsample = np.loadtxt(model_file_head + 'Nsample.txt', dtype = 'int')
        
        self.mu0_a = np.zeros((self.mmu.shape[0], self.L))
        for jmu, mu in enumerate(self.mmu):
            self.mu0_a[jmu,:] = np.loadtxt(model_file_head + 'mu0_a{}.txt'.format(jmu))
            
        # loading data
        self.p_a = {}            
        self.mu_a = {}
        self.Wij = {}
        self.T_i = {}
        self.Tcons_i = {}
        for jname, name in enumerate(data_names):
            p_arr = np.zeros((self.Nsample,self.mmu.shape[0], self.q, self.L))
            mu_arr = np.zeros((self.Nsample,self.mmu.shape[0],self.L))
            Wij_arr = np.zeros((self.Nsample,self.mmu.shape[0],self.q,self.q))
            T_i_arr = np.zeros((self.Nsample,self.mmu.shape[0], self.q, self.L))
            Tcons_i_arr = np.zeros((self.Nsample,self.mmu.shape[0], self.q, self.L))
            for jsample in xrange(self.Nsample):
                for jmu, mu in enumerate(self.mmu):
                    p_arr[jsample,jmu,:,:] = np.loadtxt(dir_name + '{}/sample{}_mu{}_p_a.txt'.format(name, jsample, jmu))
                    if name != 'aln':
                        mu_arr[jsample,jmu,:] = np.loadtxt(dir_name +\
                        '{}/sample{}_mu{}_mu_a.txt'.format(name, jsample, jmu))
                        Wij_arr[jsample,jmu,:,:] = np.loadtxt(dir_name +\
                        '{}/sample{}_mu{}_Wij.txt'.format(name, jsample, jmu))
                        T_i_arr[jsample,jmu,:,:] = np.mean(np.loadtxt(dir_name +\
                        '{}/sample{}_mu{}_T_i.txt'.format(name, jsample,jmu)))
                        Tcons_i_arr[jsample,jmu,:,:] = np.mean(np.loadtxt(dir_name +\
                        '{}/sample{}_mu{}_Tcons_i.txt'.format(name, jsample,jmu)))
            self.p_a[name] = p_arr
            self.mu_a[name] = mu_arr
            self.Wij[name] = Wij_arr
            self.T_i[name] = T_i_arr
            self.Tcons_i[name] = Tcons_i_arr
            
                    
if __name__=="__main__":
    '''Generating sequences for a given tree using a GTR matrix'''
    
    plt.ioff()
    dir_name = '/ebio/ag-neher/share/users/vpuller/GTR_staggered/L100_corr1/' 
    outdir_name = dir_name + 'plots/'
    if not os.path.exists(outdir_name):
        os.makedirs(outdir_name)    

#    # loading data
#    data_names = ['aln','GTR','GTRanc','GTRtree', 'GTRtree_site','GTRtree_unique'] 
#    data_names = ['aln','GTR','GTRanc','GTRtree', 'GTRtree_unique'] 
    data_names = ['aln','GTR','GTRanc','GTRanc_rescaled', 'GTRtree', 'GTRtree_unique'] 
    MPD = model_param_data(dir_name, 'model', data_names)
    means = np.zeros((len(data_names),len(MPD.mmu)))
    variances = np.zeros((len(data_names),len(MPD.mmu)))
    for jname, name in enumerate(data_names):
        dist = np.sum((MPD.p_a[name] - MPD.p0_a[np.newaxis,np.newaxis,:,:])**2, axis = (2,3))/MPD.L 
        means[jname,:], variances[jname,:] = mean_var(dist)
    
    xx = MPD.mmu*MPD.dist_to_root_mean*MPD.Wij0[np.where(np.triu(MPD.Wij0) >0)].mean()
#    leg = ('Tips alignment', 'GTR', 'GTR ancestors', 'GTR tree and ancestors')    
    leg = tuple(data_names[:-1])
    plt.figure(30,figsize = (20,10)); plt.clf()
    plt.subplot(1,2,1)
    for jdata, mean_loc in enumerate(means[:-1,:]):
        plt.errorbar(xx, np.log10(means[jdata,:]),\
        yerr = np.sqrt(variances[jdata,:])/means[jdata,:]/np.log(10))
    plt.xlabel(r'$\mu\bar{t}_{root}$', fontsize = 18)
    plt.ylabel(r'$\log_{10}\chi^2_{\pi}$', fontsize = 18)
    plt.legend(leg, fontsize = 18, loc = 0)

    plt.subplot(1,2,2)
    for jdata, mean_loc in enumerate(means[:-1,:]):
        plt.errorbar(np.log10(xx), np.log10(means[jdata,:]),\
        yerr = np.sqrt(variances[jdata,:])/means[jdata,:]/np.log(10))
    plt.xlabel(r'$\log_{10} (\mu\bar{t}_{root})$', fontsize = 18)
    plt.ylabel(r'$\log_{10}\chi^2_{\pi}$', fontsize = 18)
    plt.legend(leg, loc = 0, fontsize = 18)
    plt.savefig(outdir_name + 'dist.pdf')
    plt.close(30)
    
    
    # studying reconstructed mutation rates
    dist_to_root_GTRtree = np.zeros((MPD.Nsample,MPD.mmu.shape[0]))
    for jsample in xrange(MPD.Nsample):
        for jmu, mu in enumerate(MPD.mmu):
            dist_to_root_GTRtree[jsample,jmu] = np.mean(np.loadtxt(dir_name +\
            'GTRtree/sample{}_mu{}_dist_to_root.txt'.format(jsample,jmu)))
    
    dist_to_root_GTRtree_unique = np.zeros((MPD.Nsample,MPD.mmu.shape[0]))
    for jsample in xrange(MPD.Nsample):
        for jmu, mu in enumerate(MPD.mmu):
            dist_to_root_GTRtree_unique[jsample,jmu] = np.mean(np.loadtxt(dir_name +\
            'GTRtree_unique/sample{}_mu{}_dist_to_root.txt'.format(jsample,jmu)))    
    
    factor = MPD.dist_to_root_mean*np.sum(MPD.Wij[name].mean(axis=(0,1)))/(MPD.q*(MPD.q-1))
    cols = ['g','r','c','m']
    plt.figure(30,figsize = (10,10)); plt.clf()
    for jname, name in enumerate([data_names[j] for j in [1,2,4]]):
        mu_mean, mu_var = mean_var((MPD.mu_a[name][0,:,:]*\
        MPD.T_i[name][0,:,:,:].sum(axis=1)/MPD.Tcons_i[name][0,:,:,:].sum(axis=1)).T)
        factor_tree = factor*MPD.T_i[name][0,:,:,0].sum(axis=1)/MPD.T_i['GTR'][0,:,:,0].sum(axis=1)
        if name in ['GTRtree', 'GTRtree_unique']:
            j0 = 10
            plt.errorbar(xx[j0:], (mu_mean*factor_tree)[j0:], yerr = (np.sqrt(mu_var)*factor_tree)[j0:], color = cols[jname])
        else:
            plt.errorbar(xx, mu_mean*factor_tree, yerr = np.sqrt(mu_var)*factor_tree, color = cols[jname])
    plt.plot(xx,xx,':k')
    plt.xlabel(r'$\mu_0\bar{t}_{root}$', fontsize = 18)
    plt.ylabel(r'$\mu\bar{t}_{root}$', fontsize = 18)
    plt.legend(leg[1:], fontsize = 18, loc = 0)
    plt.savefig(outdir_name + 'mu_mu0.pdf')
    plt.close(30)
    
    
    j0 = 10
    plt.figure(30,figsize = (10,10)); plt.clf()
#    plt.subplot(1,3,1)
#    for jname, name in enumerate([data_names[j] for j in [1,2,3]]):
##        factor_tree = T_i[name][0,:,:,0].sum(axis=1)/T_i['GTR'][0,:,:,0].sum(axis=1)
#        mu_adjusted = MPD.mu_a[name]*MPD.T_i[name].sum(axis=2)/MPD.Tcons_i[name].sum(axis=2)*\
#        MPD.T_i[name].sum(axis=2)/MPD.T_i['GTR'].sum(axis=2)
#        mu_dist = ((mu_adjusted - MPD.mmu[np.newaxis,:,np.newaxis])**2).sum(axis = 2)/MPD.L
#        dist_mean, dist_var = mean_var(mu_dist)        
##        plt.errorbar(xx, dist_mean*factor**2, yerr = np.sqrt(dist_var)*factor**2, color = cols[jname])
#        plt.errorbar(xx[j0:], np.log10(dist_mean*factor**2)[j0:],\
#        yerr = (np.sqrt(dist_var)/dist_mean/np.log(10))[j0:], color = cols[jname])
#    plt.xlabel(r'$\mu_0\bar{t}_{root}$', fontsize = 18)
#    plt.ylabel(r'$\log_{10}\chi^2_{\mu}$', fontsize = 18)
#    plt.legend(leg[1:], fontsize = 18, loc = 0)
#    
#    plt.subplot(1,3,2)
#    for jname, name in enumerate([data_names[j] for j in [1,2,3]]):
#        mu_adjusted = MPD.mu_a[name]*MPD.T_i[name].sum(axis=2)/MPD.Tcons_i[name].sum(axis=2)*\
#        MPD.T_i[name].sum(axis=2)/MPD.T_i['GTR'].sum(axis=2)
#        mu_dist = np.sqrt(((mu_adjusted - MPD.mmu[np.newaxis,:,np.newaxis])**2).sum(axis = 2))/MPD.L
#        dist_mean, dist_var = mean_var(mu_dist)        
##        plt.errorbar(xx[j0:], np.log10(dist_mean*factor)[j0:],\
##        yerr = (np.sqrt(dist_var)/dist_mean/np.log(10))[j0:], color = cols[jname])
#        plt.errorbar(xx[j0:], dist_mean[j0:]/MPD.mmu[j0:],\
#        yerr = np.sqrt(dist_var)[j0:]/MPD.mmu[j0:], color = cols[jname])
    for jname, name in enumerate([data_names[j] for j in [1,2,3]]):
        mu_adjusted = MPD.mu_a[name]*MPD.T_i[name].sum(axis=2)/MPD.Tcons_i[name].sum(axis=2)*\
        MPD.T_i[name].sum(axis=2)/MPD.T_i['GTR'].sum(axis=2)
        mu_dist = np.sqrt(((mu_adjusted - MPD.mu0_a[np.newaxis,:,:])**2).sum(axis = 2))/MPD.L
        dist_mean, dist_var = mean_var(mu_dist)        
        plt.errorbar(xx[j0:], dist_mean[j0:]/MPD.mmu[j0:],\
        yerr = np.sqrt(dist_var)[j0:]/MPD.mmu[j0:], color = cols[jname])
    plt.xlabel(r'$\mu_0\bar{t}_{root}$', fontsize = 18)
    plt.ylabel(r'$\log_{10}\chi^2_{\mu}/(2\mu)$', fontsize = 18)
    plt.legend(leg[1:], fontsize = 18, loc = 0)
    
#    plt.subplot(1,3,3)
#    for jname, name in enumerate([data_names[j] for j in [1,2,3]]):
#        mu_adjusted = MPD.mu_a[name]*MPD.T_i[name].sum(axis=2)/MPD.Tcons_i[name].sum(axis=2)*\
#        MPD.T_i[name].sum(axis=2)/MPD.T_i['GTR'].sum(axis=2)
#        mu_dist = np.sqrt(((mu_adjusted - MPD.mmu[np.newaxis,:,np.newaxis])**2).sum(axis = 2))/MPD.L
#        dist_mean, dist_var = mean_var(mu_dist)        
##        plt.errorbar(xx[j0:], np.log10(dist_mean*factor)[j0:],\
##        yerr = (np.sqrt(dist_var)/dist_mean/np.log(10))[j0:], color = cols[jname])
#        plt.errorbar(xx[j0:], np.log10(dist_mean[j0:]/MPD.mmu[j0:]),\
#        yerr = (np.sqrt(dist_var)[j0:]/dist_mean[j0:])/np.log(10), color = cols[jname])
#    plt.xlabel(r'$\mu_0\bar{t}_{root}$', fontsize = 18)
#    plt.ylabel(r'$\log_{10}\chi^2_{\mu}/(2\mu)$', fontsize = 18)
#    plt.legend(leg[1:], fontsize = 18, loc = 0)
    plt.savefig(outdir_name + 'mu_dist.pdf')
    plt.close(30)


#    mu_adjusted = MPD.mu_a['GTR']*MPD.T_i['GTR'].sum(axis=2)/MPD.Tcons_i['GTR'].sum(axis=2)
#    plt.figure(40, figsize = (24,30)); plt.clf()
#    for jmu, mu in enumerate(MPD.mmu):
#        plt.subplot(5,4,j+1)
##        plt.hist(MPD.mu_a['GTR'][0,jmu,:]/mu)
#        plt.hist(mu_adjusted[0,jmu,:]/mu)
#        plt.xlabel(r'$\mu$')
#    plt.savefig(outdir_name + 'mu_hist.pdf')
#    plt.close(40)
    

#    # comparing simulations for different sequence lengths    
#    LL = [100, 400, 800]
#    styles = ['-o','--s',':v']
#    colors = ['b','g','r','c','m']
#    plt.figure(40); plt.clf()
#    for jL, L in enumerate(LL):
#        print jL, L
#        dir_name_loc = '/ebio/ag-neher/share/users/vpuller/GTR_staggered/L{}/'.format(L) 
#        MPD_loc = model_param_data(dir_name_loc, 'model', data_names)
#        means = np.zeros((len(data_names),len(MPD_loc.mmu)))
#        variances = np.zeros((len(data_names),len(MPD_loc.mmu)))
#        for jname, name in enumerate(data_names):
#            dist = np.sum((MPD_loc.p_a[name] - MPD_loc.p0_a[np.newaxis,np.newaxis,:,:])**2, axis = (2,3))/MPD_loc.L 
#            means[jname,:], variances[jname,:] = mean_var(dist)
#        
#        xx = MPD_loc.mmu*MPD_loc.dist_to_root_mean*MPD_loc.Wij0[np.where(np.triu(MPD_loc.Wij0) >0)].mean()
#        for jdata, mean_loc in enumerate(means[:4,:]):
#            plt.errorbar(xx, np.log10(means[jdata,:]),\
#            yerr = np.sqrt(variances[jdata,:])/means[jdata,:]/np.log(10), fmt = styles[jL] + colors[jdata])
#    plt.xlabel(r'$\mu\bar{t}_{root}$', fontsize = 18)
#    plt.ylabel(r'$\log_{10}\chi^2_{\pi}$', fontsize = 18)
#    plt.legend(leg, fontsize = 18, loc = 0)
#    plt.savefig(outdir_name + 'distL.pdf')
#    plt.close(40)
    

    plt.figure(50, figsize = (20, 5*MPD.mmu.shape[0])); plt.clf()
    for jmu, mu in enumerate(MPD.mmu):
        plt.subplot(MPD.mmu.shape[0],1,jmu+1)
        plt.plot()
        for name in data_names:
            S_a = - np.sum((h+ MPD.p_a[name][0,jmu,:,:])*np.log(h+ MPD.p_a[name][0,jmu,:,:]), axis = 0)
            plt.plot(S_a/np.log(4))
        plt.xlabel('site'); plt.ylabel('entropy')
        plt.legend(data_names)
        S_a = -np.sum((h + MPD.p0_a)*np.log(h + MPD.p0_a), axis=0)
        plt.plot(S_a/np.log(4), '-k')
    plt.savefig(outdir_name + 'entropy.pdf')
    plt.close(50)
    