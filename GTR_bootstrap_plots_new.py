# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 11:56:36 2016

@author: vpuller
"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import os
import pickle

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
        self.p_a = {}; self.mu_a = {}; self.Wij = {}
        self.n_ij = {}; self.T_i = {}; self.Tcons_i = {}
        for jname, name in enumerate(data_names):
            print name
            p_arr = np.zeros((self.Nsample,self.mmu.shape[0], self.q, self.L))
            mu_arr = np.zeros((self.Nsample,self.mmu.shape[0],self.L))
            Wij_arr = np.zeros((self.Nsample,self.mmu.shape[0],self.q,self.q))
            n_ij_arr = np.zeros((self.Nsample, self.mmu.shape[0], self.q, self.q, self.L))
            T_i_arr = np.zeros((self.Nsample,self.mmu.shape[0], self.q, self.L))
            Tcons_i_arr = np.zeros((self.Nsample,self.mmu.shape[0], self.q, self.L))
            for jsample in xrange(self.Nsample):
                for jmu, mu in enumerate(self.mmu):
#                    if name != 'aln':
                    infile_head = dir_name + '{}/sample{}_mu{}_'.format(name, jsample, jmu)
                    if os.path.exists(infile_head + 'counts.pickle'):
                        with open(infile_head + 'counts.pickle','r') as infile:
                            (n_ij_arr[jsample, jmu,:,:,:], T_i_arr[jsample, jmu,:,:],\
                            Tcons_i_arr[jsample, jmu,:,:], root_states) = pickle.load(infile)
                        with open(infile_head + 'GTRinf.pickle','r') as infile:
                            (Wij_arr[jsample,jmu,:,:], p_arr[jsample,jmu,:,:], mu_arr[jsample,jmu,:]) = pickle.load(infile)
                    else:
                        p_arr[jsample,jmu,:,:] = np.loadtxt(dir_name + '{}/sample{}_mu{}_p_a.txt'.format(name, jsample, jmu))
            self.p_a[name] = p_arr
            self.mu_a[name] = mu_arr
            self.Wij[name] = Wij_arr
            self.n_ij[name] = n_ij_arr
            self.T_i[name] = T_i_arr
            self.Tcons_i[name] = Tcons_i_arr
    
    def chi2_pi(self, data_names):
#        MPD = model_param_data(dir_name, 'model', data_names)
        means = np.zeros((len(data_names),len(self.mmu)))
        variances = np.zeros((len(data_names),len(self.mmu)))
        for jname, name in enumerate(data_names):
            dist = np.sum((self.p_a[name] - self.p0_a[np.newaxis,np.newaxis,:,:])**2, axis = (2,3))/self.L 
            means[jname,:], variances[jname,:] = mean_var(dist)
        return means, variances
        
    def chi2_pi_plot(self, data_names, to_file, scale_by = 'branch_length'):
        if scale_by == 'dist_to_root':
            xx = self.mmu*self.dist_to_root_mean*self.Wij0.sum()/(self.q*(self.q-1))
            xlab = r'$\mu\bar{t}_{root}$'
        elif scale_by == 'tree_length':
            xx = self.mmu*self.dist_to_root_mean*self.Wij0.sum()/(self.q*(self.q-1))
            xlab = r'$\mu T_{tree}$'
        else:
            xx = self.mmu*self.branch_lengths.mean()*self.Wij0.sum()/(self.q*(self.q-1))
            xlab = r'$\mu\bar{t}_{branch}$' 
#        leg = tuple(data_names)
        means, variances = self.chi2_pi(data_names)
        plt.figure(30,figsize = (20,10)); plt.clf()
        plt.subplot(1,2,1)
        for jdata, mean_loc in enumerate(means):
            plt.errorbar(xx, np.log10(means[jdata,:]),\
            yerr = np.sqrt(variances[jdata,:])/means[jdata,:]/np.log(10))
        plt.xlabel(xlab, fontsize = 18)
        plt.ylabel(r'$\log_{10}\chi^2_{\pi}$', fontsize = 18)
        plt.legend(leg, fontsize = 18, loc = 0)
    
        plt.subplot(1,2,2)
        for jdata, mean_loc in enumerate(means):
            plt.errorbar(np.log10(xx), np.log10(means[jdata,:]),\
            yerr = np.sqrt(variances[jdata,:])/means[jdata,:]/np.log(10))
        plt.xlabel(xlab, fontsize = 18)
        plt.ylabel(r'$\log_{10}\chi^2_{\pi}$', fontsize = 18)
        plt.legend(leg, loc = 0, fontsize = 18)
        plt.savefig(to_file)
        plt.close(30)  
        return None
                    
if __name__=="__main__":
    '''Generating sequences for a given tree using a GTR matrix'''
    
    plt.ioff()
    dir_name = '/ebio/ag-neher/share/users/vpuller/GTR_staggered/L100_RRE/'
    outdir_name = dir_name + 'plots/'
    if not os.path.exists(outdir_name):
        os.makedirs(outdir_name)    

#    # loading data
#    data_names = ['aln','GTR','GTRanc','GTRanc_rescaled', 'GTRtree', 'GTRtree_unique', 'GTRtree_site'] 
    data_names = ['aln','GTR','GTRanc','GTRanc_rescaled', 'GTRtree', 'GTRtree_unique'] 
    MPD = model_param_data(dir_name, 'model', data_names)
    MPD.chi2_pi_plot(data_names, outdir_name + 'dist.pdf', scale_by = 'dist_to_root')
#    means = np.zeros((len(data_names),len(MPD.mmu)))
#    variances = np.zeros((len(data_names),len(MPD.mmu)))
#    for jname, name in enumerate(data_names):
#        dist = np.sum((MPD.p_a[name] - MPD.p0_a[np.newaxis,np.newaxis,:,:])**2, axis = (2,3))/MPD.L 
#        means[jname,:], variances[jname,:] = mean_var(dist)
#    
#    
##    xx = MPD.mmu*MPD.dist_to_root_mean*MPD.Wij0.sum()/(MPD.q*(MPD.q-1))
#    xx = MPD.mmu*MPD.branch_lengths.mean()*MPD.Wij0.sum()/(MPD.q*(MPD.q-1))
##    xx = np.zeros(MPD.mmu.shape)
##    T = MPD.T_i[name][0,:,:,:].sum(axis=1).mean()
##    for jmu, mu in enumerate(MPD.mmu):    
##        Q0mean = np.mean(MPD.mu0_a[jmu,:]*MPD.p0_a[:,np.newaxis,:]*MPD.Wij0[:,:,np.newaxis], axis=2)
##        xx[jmu] = T*Q0mean.sum()/Q0mean.shape[0]
##    leg = ('Tips alignment', 'GTR', 'GTR ancestors', 'GTR tree and ancestors')    
#    leg = tuple(data_names)
#    plt.figure(30,figsize = (20,10)); plt.clf()
#    plt.subplot(1,2,1)
#    for jdata, mean_loc in enumerate(means):
#        plt.errorbar(xx, np.log10(means[jdata,:]),\
#        yerr = np.sqrt(variances[jdata,:])/means[jdata,:]/np.log(10))
#    plt.xlabel(r'$\mu\bar{t}_{root}$', fontsize = 18)
##    plt.xlabel(r'$\mu T_{tree}$', fontsize = 18)
#    plt.ylabel(r'$\log_{10}\chi^2_{\pi}$', fontsize = 18)
#    plt.legend(leg, fontsize = 18, loc = 0)
#
#    plt.subplot(1,2,2)
#    for jdata, mean_loc in enumerate(means):
#        plt.errorbar(np.log10(xx), np.log10(means[jdata,:]),\
#        yerr = np.sqrt(variances[jdata,:])/means[jdata,:]/np.log(10))
#    plt.xlabel(r'$\log_{10} (\mu\bar{t}_{root})$', fontsize = 18)
##    plt.xlabel(r'$\log_{10} (\mu T_{tree})$', fontsize = 18)
#    plt.ylabel(r'$\log_{10}\chi^2_{\pi}$', fontsize = 18)
#    plt.legend(leg, loc = 0, fontsize = 18)
#    plt.savefig(outdir_name + 'dist.pdf')
#    plt.close(30)
    
    
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
    xx = MPD.mmu*MPD.dist_to_root_mean*MPD.Wij0[np.where(np.triu(MPD.Wij0) >0)].mean()
    cols = ['g','r','c','m']
    plt.figure(30,figsize = (10,10)); plt.clf()
    for jname, name in enumerate([data_names[j] for j in [1,2,3,4]]):
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
    plt.legend([data_names[j] for j in [1,2,3,4]], fontsize = 18, loc = 0)
    plt.savefig(outdir_name + 'mu_mu0.pdf')
    plt.close(30)
    
    
    j0 = 10
    plt.figure(30,figsize = (10,10)); plt.clf()
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
    

    plt.figure(50, figsize = (20, 6*MPD.mmu.shape[0])); plt.clf()
    for jmu, mu in enumerate(MPD.mmu):
        plt.subplot(MPD.mmu.shape[0],1,jmu+1)
        for name in data_names:
            S_a = - np.sum((h+ MPD.p_a[name][0,jmu,:,:])*np.log(h+ MPD.p_a[name][0,jmu,:,:]), axis = 0)
            plt.plot(S_a/np.log(4))
        plt.xlabel('site'); plt.ylabel('entropy')
        plt.legend(data_names)
        plt.title('jmu = {}, mu = {}'.format(jmu, mu))
        S_a = -np.sum((h + MPD.p0_a)*np.log(h + MPD.p0_a), axis=0)
        plt.plot(S_a/np.log(4), '-k')
    plt.savefig(outdir_name + 'entropy.pdf')
    plt.close(50)
    
    names_loc = ['GTR', 'GTRanc', 'GTRtree']
    plt.figure(50, figsize = (6*len(names_loc), 6*MPD.mmu.shape[0])); plt.clf()
    for jmu, mu in enumerate(MPD.mmu):
        for jname, name in enumerate(names_loc):
            plt.subplot(MPD.mmu.shape[0],len(names_loc),len(names_loc)*jmu + jname + 1)
            plt.hist(MPD.n_ij[name][0,jmu,:,:,:].sum(axis=(0,1)))
            plt.xlabel('n_ij.sum()'); plt.ylabel('jmu = {}'.format(jmu))
            plt.title(name)
    plt.savefig(outdir_name + 'n_ij_sum.pdf')
    plt.close(50)
    
    Ttree = {}
    plt.figure(60); plt.clf()
    for name in data_names[1:]: 
        Ttree[name] = MPD.T_i[name][0,:,:,:].sum(axis=1).mean(axis=1)
        if name in ['GTR', 'GTRanc', 'GTRtree_site']:
            factor = np.zeros(MPD.mmu.shape)
            for jmu, mu in enumerate(MPD.mmu):
                Q0mean = np.mean(MPD.mu0_a[jmu,:]*MPD.p0_a[:,np.newaxis,:]*MPD.Wij0[:,:,np.newaxis], axis=2)
                factor[jmu] = Q0mean.sum()/Q0mean.shape[0]
            plt.plot(MPD.mmu, Ttree[name]*factor)
        elif name in ['GTRanc_rescaled']:
            plt.plot(MPD.mmu, Ttree[name]*.75)
        else:
            plt.plot(MPD.mmu, Ttree[name])
    plt.legend(data_names[1:], loc = 0)
#    plt.legend(['GTRanc', 'GTRanc_rescaled', 'GTRtree'], loc = 0)
    plt.xlabel(r'$\mu$'); plt.ylabel('T') 
    plt.savefig(outdir_name + 'tree_lengths.pdf')
    plt.close(60)
    