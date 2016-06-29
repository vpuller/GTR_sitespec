# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 11:56:36 2016

@author: vpuller
"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import os
import pickle as pickle
from Bio import Phylo

# Constants
h = 10**(-8)

def mean_var(dist):
    dist_mean = dist.mean(axis=0)
    dist_var = (dist**2).mean(axis=0) - dist_mean**2 
    return dist_mean, dist_var
        
class model_param_data(object):
    '''Loading model parametes and data'''
    def __init__(self, dir_name, model_name, data_names, load_counts = False, load_Like = False):
        '''loading model, simulation parameters and simulation data'''
        model_file_head = dir_name + model_name + '/'
        with open(model_file_head + 'seq0.txt','r') as seq0_file:
            seq0 = seq0_file.read()
        with open(model_file_head + 'alphabet.txt','r') as alpha_file:
            self.alphabet = alpha_file.read()
        self.q = len(self.alphabet); self.L = len(seq0)
        
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
        self.W1_ij = {}; self.p1_i = {}
        self.logL = {}
#        self.logL1 = np.zeros((self.Nsample,self.mmu.shape[0]))
#        self.logL1exact = np.zeros((self.Nsample,self.mmu.shape[0]))
        for jname, name in enumerate(data_names):
#            print name
#            n_ij_arr = np.zeros((self.Nsample, self.mmu.shape[0], self.q, self.q, self.L))
#            T_i_arr = np.zeros((self.Nsample,self.mmu.shape[0], self.q, self.L))
#            Tcons_i_arr = np.zeros((self.Nsample,self.mmu.shape[0], self.q, self.L))            
            if name == 'aln':
                p_arr = np.zeros((self.Nsample,self.mmu.shape[0], self.q, self.L))
                for jsample in xrange(self.Nsample):
                    for jmu, mu in enumerate(self.mmu):
                        p_arr[jsample,jmu,:,:] = np.loadtxt(dir_name + '{}/sample{}_mu{}_p_a.txt'.format(name, jsample, jmu))
                self.p_a[name] = p_arr
            else:
                (Wij_arr, p_arr, mu_arr) = self.load_GTRmodels(dir_name, name)
                self.p_a[name] = p_arr
                self.mu_a[name] = mu_arr
                self.Wij[name] = Wij_arr
                
                if load_counts:
                    (n_ij_arr, T_i_arr, Tcons_i_arr) = self.load_counts(dir_name, name)
                    self.n_ij[name] = n_ij_arr
                    self.T_i[name] = T_i_arr
                    self.Tcons_i[name] = Tcons_i_arr
                
                if load_Like:
                    self.logL[name] = self.load_Likelihoods(dir_name, name)

#                    for jmu, mu in enumerate(self.mmu):
#                        infile_head = dir_name + '{}/sample{}_mu{}_'.format(name, jsample, jmu)
#                        if os.path.exists(infile_head + 'counts.pickle'):
#                            with open(infile_head + 'counts.pickle','r') as infile:
#                                (n_ij_arr[jsample, jmu,:,:,:], T_i_arr[jsample, jmu,:,:],\
#                                Tcons_i_arr[jsample, jmu,:,:], root_states) = pickle.load(infile)
#                            with open(infile_head + 'GTRinf.pickle','r') as infile:
#                                (Wij_arr[jsample,jmu,:,:], p_arr[jsample,jmu,:,:], mu_arr[jsample,jmu,:]) = pickle.load(infile)
                        
#            self.p_a[name] = p_arr
#            self.mu_a[name] = mu_arr
#            self.Wij[name] = Wij_arr
#            self.n_ij[name] = n_ij_arr
#            self.T_i[name] = T_i_arr
#            self.Tcons_i[name] = Tcons_i_arr
    
    def load_GTRmodels(self,dir_name, name):
        print 'loading reconstructed model: ' + name, self.q, self.L
        p_arr = np.zeros((self.Nsample,self.mmu.shape[0], self.q, self.L))
        mu_arr = np.zeros((self.Nsample,self.mmu.shape[0],self.L))
        Wij_arr = np.zeros((self.Nsample,self.mmu.shape[0],self.q,self.q))
        for jsample in xrange(self.Nsample):
            for jmu, mu in enumerate(self.mmu):
                infile_head = dir_name + '{}/sample{}_mu{}_'.format(name, jsample, jmu)
                if os.path.exists(infile_head + 'GTRinf.pickle'):
                    with open(infile_head + 'GTRinf.pickle','r') as infile:
                        (Wij_arr[jsample,jmu,:,:], p_arr[jsample,jmu,:,:], mu_arr[jsample,jmu,:]) = pickle.load(infile)
        return (Wij_arr, p_arr, mu_arr)
    
    def load_counts(self,dir_name, name):
        print 'loading counts: ' + name, self.q, self.L
        n_ij_arr = np.zeros((self.Nsample, self.mmu.shape[0], self.q, self.q, self.L))
        T_i_arr = np.zeros((self.Nsample,self.mmu.shape[0], self.q, self.L))
        Tcons_i_arr = np.zeros((self.Nsample,self.mmu.shape[0], self.q, self.L)) 
        for jsample in xrange(self.Nsample):
            for jmu, mu in enumerate(self.mmu):
                infile_head = dir_name + '{}/sample{}_mu{}_'.format(name, jsample, jmu)
                if os.path.exists(infile_head + 'GTRinf.pickle'):
                    with open(infile_head + 'counts.pickle','r') as infile:
                        (n_ij_arr[jsample, jmu,:,:,:], T_i_arr[jsample, jmu,:,:],\
                        Tcons_i_arr[jsample, jmu,:,:], root_states) = pickle.load(infile)
#                    with open(infile_head + 'GTRinf.pickle','r') as infile:
#                        (Wij_arr[jsample,jmu,:,:], p_arr[jsample,jmu,:,:], mu_arr[jsample,jmu,:]) = pickle.load(infile)
        return (n_ij_arr, T_i_arr, Tcons_i_arr)
    
    def load_Likelihoods(self,dir_name, name):
        print 'loading likelihoods model: ' + name, self.q, self.L
        logL_arr = np.zeros((self.Nsample,self.mmu.shape[0], 3))
        for jsample in xrange(self.Nsample):
            for jmu, mu in enumerate(self.mmu):
                infile_head = dir_name + '{}/sample{}_mu{}_'.format(name, jsample, jmu)
                if os.path.exists(infile_head + 'Likelihoods.txt'):
                    logL_arr[jsample,jmu,:] = np.loadtxt(infile_head + 'Likelihoods.txt')
#                    with open(infile_head + 'GTRinf.pickle','r') as infile:
#                        (Wij_arr[jsample,jmu,:,:], p_arr[jsample,jmu,:,:], mu_arr[jsample,jmu,:]) = pickle.load(infile)
        return logL_arr
        
    def chi2_pi(self, data_names):
        '''
        Calculating divergence between the model and the inferred
        nucleotide frequencies
        '''

        means = np.zeros((len(data_names),len(self.mmu)))
        variances = np.zeros((len(data_names),len(self.mmu)))
        for jname, name in enumerate(data_names):
            dist = np.sum((self.p_a[name] - self.p0_a[np.newaxis,np.newaxis,:,:])**2, axis = (2,3))/self.L 
            means[jname,:], variances[jname,:] = mean_var(dist)
        return means, variances
    
    def chi2_pi1(self, data_names):
        means = np.zeros((len(data_names),len(self.mmu)))
        variances = np.zeros((len(data_names),len(self.mmu)))
        for jname, name in enumerate(data_names):
            p_a = np.tile(self.p1_i[name][:,:,:,np.newaxis],(1,1,1,self.L))
            dist = np.sum((p_a - self.p0_a[np.newaxis,np.newaxis,:,:])**2, axis = (2,3))/self.L 
            means[jname,:], variances[jname,:] = mean_var(dist)
#            dist = np.sum((self.p1_i[name][:,:,:,np.newaxis] -\
#            self.p0_a[np.newaxis,np.newaxis,:,:])**2, axis = (2,3))/self.L 
#            means[jname,:], variances[jname,:] = mean_var(dist)
        return means, variances
    
    def KLmean_pi(self, data_names):
        '''
        Calculating KL divergence between the model and the inferred
        nucleotide frequencies
        '''
        means = np.zeros((len(data_names),len(self.mmu)))
        variances = np.zeros((len(data_names),len(self.mmu)))
        for jname, name in enumerate(data_names):
            KL = - np.sum(self.p0_a[np.newaxis,np.newaxis,:,:]*(np.log2(self.p_a[name] + h) -\
            np.log2(self.p0_a[np.newaxis,np.newaxis,:,:]+h)), axis = (2,3))/self.L
            means[jname,:], variances[jname,:] = mean_var(KL)
        return means, variances
        
    def scale_mu(self, scale_by):
        if scale_by == 'dist_to_root':
            factor = self.dist_to_root_mean*self.Wij0.sum()/(self.q*(self.q-1))
            xlab = r'$\mu\bar{t}_{root}$'
        elif scale_by == 'tree_length':
            factor = self.T*self.Wij0.sum()/(self.q*(self.q-1))
            xlab = r'$\mu T_{tree}$'
        else:
            factor = self.branch_lengths.mean()*self.Wij0.sum()/(self.q*(self.q-1))
            xlab = r'$\mu\bar{t}_{branch}$' 
        xx = self.mmu*factor
        return xx, xlab, factor
        
    def chi2_pi_plot(self, data_names, to_file, scale_by = 'branch_length', oneforall = False, data_legend = None):
        xx, xlab, factor = self.scale_mu(scale_by)

        if oneforall:
            means, variances = self.chi2_pi1(data_names)
        else:
            means, variances = self.chi2_pi(data_names)
        
        if data_legend is None:
            data_legend = list(data_names)
        plt.figure(30,figsize = (20,10)); plt.clf()
        plt.subplot(1,2,1)
        for jdata, mean_loc in enumerate(means):
            plt.errorbar(xx, np.log10(means[jdata,:]),\
            yerr = np.sqrt(variances[jdata,:])/means[jdata,:]/np.log(10))
        plt.xlabel(xlab, fontsize = 22)
        plt.ylabel(r'$\log_{10}\chi^2_{\pi}$', fontsize = 22)
        plt.legend(data_legend, fontsize = 22, loc = 0)
    
        plt.subplot(1,2,2)
        for jdata, mean_loc in enumerate(means):
            plt.errorbar(np.log10(xx), np.log10(means[jdata,:]),\
            yerr = np.sqrt(variances[jdata,:])/means[jdata,:]/np.log(10))
        plt.xlabel(xlab, fontsize = 22)
        plt.ylabel(r'$\log_{10}\chi^2_{\pi}$', fontsize = 22)
        plt.legend(data_legend, loc = 0, fontsize = 22)
        plt.savefig(to_file)
        plt.close(30)  
        return None
    
    def KLmean_plot(self, data_names, to_file, scale_by = 'branch_length', oneforall = False):
        xx, xlab, factor = self.scale_mu(scale_by)
        means, variances = self.KLmean_pi(data_names)
            
        plt.figure(30,figsize = (10,10)); plt.clf()
        for jdata, mean_loc in enumerate(means):
#            plt.errorbar(xx, np.log2(means[jdata,:]),\
#            yerr = np.sqrt(variances[jdata,:])/means[jdata,:]/np.log(2))
            plt.errorbar(xx, means[jdata,:], yerr = np.sqrt(variances[jdata,:]))
        plt.xlabel(xlab, fontsize = 18)
        plt.ylabel(r'$KL_{mean}$', fontsize = 18)
        plt.legend(data_names, fontsize = 18, loc = 0)
        plt.savefig(to_file)
        plt.close(30)  
        return None
        
    def chi2_mu_plot(self, data_names, to_file, j0 = 0, scale_by = 'branch_length', data_legend = None):
        xx, xlab, factor = self.scale_mu(scale_by)

        if data_legend is None:
            data_legend = list(data_names)
            
        plt.figure(30,figsize = (10,10)); plt.clf()
        for jname, name in enumerate(data_names):
            mu_adjusted = self.mu_a[name]*self.T_i[name].sum(axis=2)/self.Tcons_i[name].sum(axis=2)*\
            self.T_i[name].sum(axis=2)/self.T_i['GTR'].sum(axis=2)
            mu_dist = np.sqrt(((mu_adjusted - self.mu0_a[np.newaxis,:,:])**2).sum(axis = 2))/self.L
            dist_mean, dist_var = mean_var(mu_dist)        
            plt.errorbar(xx[j0:], dist_mean[j0:]/self.mmu[j0:],\
            yerr = np.sqrt(dist_var)[j0:]/self.mmu[j0:])
        plt.xlabel(xlab, fontsize = 18)
        plt.ylabel(r'$\chi^2_{\mu}/(2\mu)$', fontsize = 18)
        plt.legend(data_legend, fontsize = 18, loc = 0)
        plt.savefig(to_file)
        plt.close(30)
        return None

    def muave_mu_plot(self, data_names, to_file, j0 = 0, scale_by = 'branch_length', data_legend = None):
        xx, xlab, f = self.scale_mu(scale_by)
        if scale_by == 'dist_to_root':
            f0 = self.dist_to_root_mean
        elif scale_by == 'tree_length':
            f0 = self.T
        else:
            f0 = self.branch_lengths.mean()
        
        if data_legend is None:
            data_legend = list(data_names)
            
        plt.figure(30,figsize = (10,10)); plt.clf()
        for jname, name in enumerate(data_names):
#            mu_mean, mu_var = mean_var((self.mu_a[name][0,:,:]*\
#            self.T_i[name][0,:,:,:].sum(axis = 1)/self.Tcons_i[name][0,:,:,:].sum(axis = 1)).T)
#            print xx, mu_mean, mu_var
#            plt.errorbar(xx[j0:], f*mu_mean[j0:], yerr = f*np.sqrt(mu_var)[j0:])
            mu_mean, mu_var = mean_var((self.mu_a[name][0,:,:]*\
            self.T_i[name][0,:,:,:].sum(axis=1)/self.Tcons_i[name][0,:,:,:].sum(axis=1)).T)
##            factor_tree = np.ones(self.mmu.shape[0])*f
##            factor_tree = f0*np.sum(self.Wij[name][0,:,:,:], axis = (1,2))/(self.q*(self.q-1))
##            mu_mean, mu_var = mean_var((self.mu_a[name][0,:,:]*\
##            self.T_i[name][0,:,:,:].sum(axis=1)/self.Tcons_i[name][0,:,:,:].sum(axis=1)).T)
#            factor = f0*np.sum(self.Wij[name][0,:,:,:], axis = (1,2))/(self.q*(self.q-1))
#            factor_tree = factor*self.T_i[name][0,:,:,:].sum(axis = (1,2))/self.T_i['GTR'][0,:,:,:].sum(axis = (1,2))
#            plt.errorbar(xx[j0:], (mu_mean*factor_tree)[j0:], yerr = (np.sqrt(mu_var)*factor_tree)[j0:])
            fact = f0*np.sum(self.Wij[name][0,:,:,:], axis = (1,2))/(self.q*(self.q-1))
            factor = fact*self.T_i[name][0,:,:,:].sum(axis = (1,2))/self.T_i['GTR'][0,:,:,:].sum(axis = (1,2))         
            plt.errorbar(xx[j0:], (mu_mean*factor)[j0:], yerr = (np.sqrt(mu_var)*factor)[j0:])
        plt.plot(xx,xx,':k')
        plt.xlabel(xlab, fontsize = 18)
        plt.ylabel(xlab, fontsize = 18)
        plt.legend(data_legend, fontsize = 18, loc = 0)
        plt.savefig(to_file)
        plt.close(30)
        return None
    
    def site_entropy_plot(self, data_names, to_file):
        leg = list(data_names)
        leg.append('model')
        plt.figure(50, figsize = (20, 6*self.mmu.shape[0])); plt.clf()
        for jmu, mu in enumerate(self.mmu):
            plt.subplot(self.mmu.shape[0],1,jmu+1)
            for name in data_names:
                S_a = - np.sum((h+ self.p_a[name][0,jmu,:,:])*np.log(h+ self.p_a[name][0,jmu,:,:]), axis = 0)
                plt.plot(S_a/np.log(4))
            S_a = -np.sum((h + self.p0_a)*np.log(h + self.p0_a), axis=0)
            plt.plot(S_a/np.log(4), '-k')
            plt.xlabel('site'); plt.ylabel('entropy')
            plt.legend(leg)
            plt.title('jmu = {}, mu = {:.2f}'.format(jmu, mu))
        plt.savefig(to_file)
        plt.close(50)
        return None
    
    def inference_vs_model_freq_plot(self, name, to_file):
        plt.figure(55, figsize = (6*self.q, 6*self.mmu.shape[0])); plt.clf()
        for jmu, mu in enumerate(self.mmu):
            for jnuc, nuc in enumerate(self.alphabet):
                plt.subplot(self.mmu.shape[0], self.q, jmu*self.q + jnuc + 1)
                p_mean, p_var = mean_var(self.p_a[name][:,jmu,jnuc,:])
                plt.errorbar(self.p0_a[jnuc,:], p_mean, yerr = np.sqrt(p_var), fmt = '')
                plt.xlabel('model'); plt.ylabel('inferred')
                plt.title(r'$\pi_{}, jmu = {}$'.format(nuc, jmu))
        plt.savefig(to_file)
        plt.close(55)
        return None
        
    def nij_hist(self, data_names, to_file):
        plt.figure(50, figsize = (6*len(data_names), 6*self.mmu.shape[0])); plt.clf()
        for jmu, mu in enumerate(self.mmu):
            for jname, name in enumerate(data_names):
                plt.subplot(self.mmu.shape[0],len(data_names),len(data_names)*jmu + jname + 1)
                plt.hist(self.n_ij[name][0,jmu,:,:,:].sum(axis=(0,1)))
                plt.xlabel('n_ij.sum()')
                plt.ylabel('jmu = {}, mu = {:.2f}'.format(jmu, mu))
                plt.title(name)
        plt.savefig(to_file)
        plt.close(50)
        return None
        
#    def Akaike_plot(self, data_names, to_file, scale_by = 'branch_length'):
#        xx, xlab, f = self.scale_mu(scale_by)
#        k = self.q*(self.q-1)//2 - 1 + self.q*self.L
#        k1 = self.q*(self.q+1)//2
#        
#        leg = list(data_names)
#        leg.extend(['GTR1','GTR1exact'])
#        plt.figure(30,figsize = (20,10)); plt.clf()
#        plt.subplot(1,2,1)
#        for name in data_names:
#            plt.plot(xx, self.logL[name].sum(axis=0))
#        plt.plot(xx, self.logL1.sum(axis=0))
#        plt.plot(xx, self.logL1exact.sum(axis=0))
#        plt.xlabel(xlab, fontsize = 18)
#        plt.ylabel(r'$\log{L}$', fontsize = 18)
#        plt.legend(leg, fontsize = 18, loc = 0)
#    
#        plt.subplot(1,2,2)
#        for name in data_names:
#            plt.plot(xx, k - self.logL[name].sum(axis=0))
#        plt.plot(xx, k1 - self.logL1.sum(axis=0))
#        plt.plot(xx, k1 - self.logL1exact.sum(axis=0))
#        plt.xlabel(xlab, fontsize = 18)
#        plt.ylabel(r'$k - \log{L}$', fontsize = 18)
#        plt.legend(leg, loc = 0, fontsize = 18)
#        plt.savefig(to_file)
#        plt.close(30)  
#        return None
        
    def site_freq_hist(self, to_file):
        plt.figure(20, figsize = (6*self.q, 6)); plt.clf()
        for jnuc, nuc in enumerate(self.alphabet):
            plt.subplot(1, self.q, jnuc + 1)
            plt.hist(self.p0_a[jnuc,:])
        plt.savefig(to_file)
        plt.close(20) 
        return None
    
    def site_freq_dist(self, data_names, to_file, jmu = 1):
        plt.figure(20, figsize = (6*(len(data_names) + 1), 6)); plt.clf()
        plt.subplot(1, len(data_names) + 1, 1)
        for jnuc, nuc in enumerate(self.alphabet):
            nn_bin, x_bin = np.histogram(self.p0_a[jnuc,:],bins = 10, range = (0,1))
            plt.plot(x_bin[:-1] + .5*x_bin[1], nn_bin/self.L)
            plt.xlabel(r'$\pi_{i\alpha}$', fontsize = 18)
            plt.legend(self.alphabet, fontsize = 18)
            plt.title('Model', fontsize = 18)
        for jname, name in enumerate(data_names):
            plt.subplot(1, len(data_names) + 1, jname + 2)
            for jnuc, nuc in enumerate(self.alphabet):
                nn_bin, x_bin = np.histogram(self.p_a[name][0,jmu,jnuc,:],bins = 10, range = (0,1))
                plt.plot(x_bin[:-1] + .5*x_bin[1], nn_bin/self.L)
            plt.xlabel(r'$\pi_{i\alpha}$', fontsize = 18)
            plt.legend(self.alphabet, fontsize = 18)
            plt.title(name, fontsize = 18)
        plt.savefig(to_file)
        plt.close(20) 
        return None
    
    def site_freq_scatter(self, data_names, to_file, jmu = 1):
        plt.figure(20, figsize = (7*len(data_names), 7)); plt.clf()
        for jname, name in enumerate(data_names):
            plt.subplot(1, len(data_names), jname + 1)
            for jnuc, nuc in enumerate(self.alphabet):                
                p_mean, p_var = mean_var(self.p_a[name][:,jmu,jnuc,:])
                plt.errorbar(self.p0_a[jnuc,:], p_mean, yerr = np.sqrt(p_var), fmt = None)
            plt.xlabel(r'$\pi_{i\alpha}^{(0)}$', fontsize = 18)
            plt.ylabel(r'$\pi_{i\alpha}$', fontsize = 18)
            plt.legend(self.alphabet, fontsize = 18, loc = 0)
            plt.title(name, fontsize = 18)
        plt.savefig(to_file)
        plt.close(20)
        return None
        
    def Akaike_plot(self, data_names, to_file, scale_by = 'branch_length', data_legend = None):
        xx, xlab, factor = self.scale_mu(scale_by)

        if data_legend is None:
            data_legend = list(data_names)
            
        kk = np.array([self.q*self.L + self.q*(self.q-1)//2 - 1,\
        self.q-1 + self.q*(self.q-1)//2, (self.q-1 + self.q*(self.q-1)//2)*self.L])
        plt.figure(30,figsize = (8*len(data_names), 16)); plt.clf()
        for jname, name in enumerate(data_names):
            plt.subplot(2,len(data_names),jname + 1)
            for jlogL in xrange(3):
                logL_mean, logL_var = mean_var(self.logL[name][:,:,jlogL])        
                plt.errorbar(xx, logL_mean, yerr = np.sqrt(logL_var))
            plt.xlabel(xlab, fontsize = 24)
            plt.ylabel(r'$\log L$', fontsize = 24)
            plt.legend(('single W_ij', 'single GTR', 'site specific GTR'), fontsize = 18, loc = 0)
            plt.title(data_legend[jname], fontsize = 24)
            
            plt.subplot(2,len(data_names), len(data_names) + jname + 1)
            for jlogL in xrange(3):
                logL_mean, logL_var = mean_var(self.logL[name][:,:,jlogL])        
                plt.errorbar(xx, kk[jlogL] - logL_mean, yerr = np.sqrt(logL_var))
            plt.xlabel(xlab, fontsize = 24)
            plt.ylabel(r'$AIC = k - \log L$', fontsize = 24)
            plt.legend(('single W_ij', 'single GTR', 'site specific GTR'), fontsize = 18, loc = 0)
        plt.savefig(to_file)
        plt.close(30)
        return None
        
        
if __name__=="__main__":
    '''Loading data and making plots'''
    
    plt.ioff()
    plt.close('all')
    dir_name = '/ebio/ag-neher/share/users/vpuller/GTR_staggered/L100_poltree/'
    outdir_name = dir_name + 'plots/'
    if not os.path.exists(outdir_name):
        os.makedirs(outdir_name)    

    # loading data
    data_names = ['aln','GTR','GTRanc_vadim', 'GTRanc_pavel', 'GTRtree_vadim', 'GTRtree_pavel', 'GTRsite_vadim']
#    data_names = ['aln','GTR','GTRanc_vadim', 'GTRtree_vadim','GTRsite_vadim']
    data_legend = list(data_names)
#    data_legend = ['Alignment', 'Inferred', 'Inferred, ancestral', 'Inferred, tree, ancestral']
#    data_names = ['aln', 'GTR']
    MPD = model_param_data(dir_name, 'model', data_names, load_counts = False, load_Like = True)

    
    # divergence in pi
    MPD.chi2_pi_plot(data_names[:-1], outdir_name + 'dist.pdf', scale_by = 'branch_length',data_legend = data_legend)
    MPD.chi2_pi_plot(['GTRanc_vadim', 'GTRsite_vadim'], outdir_name + 'dist_sitespec.pdf',\
    scale_by = 'branch_length',data_legend = ['JC anc. rec.', 'site-specific anc. rec.'])
    
#    MPD.KLmean_plot(data_names, outdir_name + 'KL.pdf', scale_by = 'branch_length')
#    MPD.chi2_pi_plot(data_names, outdir_name + 'dist1.pdf', scale_by = 'branch_length', oneforall = False)
    
#    # site entropies
#    MPD.site_entropy_plot(data_names, outdir_name + 'entropy.pdf')
    
    # plot site frequency distributions
#    for jmu, mu in enumerate(MPD.mmu): 
    jmu = 14
    MPD.site_freq_dist(['aln', 'GTR'], outdir_name + 'site_freq_dist{}.pdf'.format(jmu), jmu = jmu)
    
    # frequency scatter-plots
#    for jmu, mu in enumerate(MPD.mmu): 
    MPD.site_freq_scatter(['aln', 'GTR'], outdir_name + 'site_freq_scatter{}.pdf'.format(jmu), jmu = jmu)
        
#    for name in data_names:
#        MPD.inference_vs_model_freq_plot(name, outdir_name + 'freq_{}.pdf'.format(name))
        
#    # n_ij histograms
#    MPD.nij_hist(data_names, outdir_name + 'n_ij_sum.pdf')
    
#    # divergence in mu
#    MPD.chi2_mu_plot(data_names[1:], outdir_name + 'mu_dist.pdf', j0 = 10,\
#    scale_by = 'branch_length', data_legend = data_legend[1:])
#      
#    
#    # average over sequence mutation rate
#    MPD.muave_mu_plot(data_names[1:], outdir_name + 'mu_mu0.pdf', j0 = 5,\
#    scale_by = 'branch_length', data_legend = data_legend[1:])

#
#    # Akaike information criterion
    MPD.Akaike_plot(data_names[1:-1], outdir_name + 'Akaike.pdf',\
    scale_by = 'branch_length',data_legend = data_legend[1:])

    
    
    # comparing simulations for different sequence lengths    
    scalingplot = False
    if scalingplot:
        LL = [100, 400]
    #    data_names = ['aln','GTR']
        styles = ['-o','--s',':v']
        colors = ['b','g','r','c','m', 'y','k']
        plt.figure(40); plt.clf()
        for jL, L in enumerate(LL):
            print jL, L
            dir_name_loc = '/ebio/ag-neher/share/users/vpuller/GTR_staggered/L{}_simplex1/'.format(L) 
            MPD_loc = model_param_data(dir_name_loc, 'model', data_names[:-1])
            xx, xlab, factor = MPD_loc.scale_mu('branch_length')
            means, variances = MPD_loc.chi2_pi(data_names)
            for jdata, mean_loc in enumerate(means):
                plt.errorbar(xx, np.log10(means[jdata,:]),\
                yerr = np.sqrt(variances[jdata,:])/means[jdata,:]/np.log(10), fmt = styles[jL] + colors[jdata])        
        plt.xlabel(xlab, fontsize = 18)
        plt.ylabel(r'$\log_{10}\chi^2_{\pi}$', fontsize = 18)
        plt.legend(data_legend, fontsize = 18, loc = 0)
        plt.savefig(outdir_name + 'distL.pdf')
        plt.close(40)

    
#    xx, xlab, f = MPD.scale_mu('branch_length')
#    Ttree = {}
#    plt.figure(60); plt.clf()
#    for name in data_names[1:]: 
#        Ttree[name] = MPD.T_i[name][0,:,:,:].sum(axis=1).mean(axis=1)
#        if name in ['GTR', 'GTRanc', 'GTRtree_site']:
#            factor = np.zeros(MPD.mmu.shape)
#            for jmu, mu in enumerate(MPD.mmu):
#                Q0mean = np.mean(MPD.mu0_a[jmu,:]*MPD.p0_a[:,np.newaxis,:]*MPD.Wij0[:,:,np.newaxis], axis=2)
#                factor[jmu] = Q0mean.sum()/Q0mean.shape[0]
#            plt.plot(xx, Ttree[name]*factor)
#        elif name in ['GTRanc_rescaled']:
#            plt.plot(xx, Ttree[name]*.75)
#        else:
#            plt.plot(xx, Ttree[name])
#    plt.legend(data_names[1:], loc = 0)
##    plt.legend(['GTRanc', 'GTRanc_rescaled', 'GTRtree'], loc = 0)
#    plt.xlabel(xlab); plt.ylabel('T') 
#    plt.savefig(outdir_name + 'tree_lengths.pdf')
#    plt.close(60)
    
#    for jsample in xrange(MPD.Nsample):
#        for jmu, mu in enumerate(MPD.mmu):
#            tree_file_name = dir_name + 'sample{}_mu{}_rectree_tmp.nwk'.format(jsample, jmu) 
#            bio_tree = Phylo.read(tree_file_name, 'newick')
#            print bio_tree.root.branch_length            
#            bio_tree.ladderize()
#            
#            fig10=plt.figure(10,figsize=(20,20))    
#            plt.clf()
#            ax10=fig10.add_subplot(1,1,1)
#            Phylo.draw(bio_tree, do_show = False, axes = ax10)
#            plt.savefig(tree_file_name[:-4] + '.pdf')
#            plt.close(10)
#    