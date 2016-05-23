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
from Bio import Phylo

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
        self.W1_ij = {}; self.p1_i = {}
        self.logL = {}
        self.logL1 = np.zeros((self.Nsample,self.mmu.shape[0]))
        self.logL1exact = np.zeros((self.Nsample,self.mmu.shape[0]))
        for jname, name in enumerate(data_names):
            print name
            p_arr = np.zeros((self.Nsample,self.mmu.shape[0], self.q, self.L))
            mu_arr = np.zeros((self.Nsample,self.mmu.shape[0],self.L))
            Wij_arr = np.zeros((self.Nsample,self.mmu.shape[0],self.q,self.q))
            n_ij_arr = np.zeros((self.Nsample, self.mmu.shape[0], self.q, self.q, self.L))
            T_i_arr = np.zeros((self.Nsample,self.mmu.shape[0], self.q, self.L))
            Tcons_i_arr = np.zeros((self.Nsample,self.mmu.shape[0], self.q, self.L))
            W1_arr = np.zeros((self.Nsample,self.mmu.shape[0],self.q,self.q))
            p1_arr = np.zeros((self.Nsample,self.mmu.shape[0], self.q))
            logL = np.zeros((self.Nsample,self.mmu.shape[0]))
            
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
                        with open(infile_head + 'GTRone.pickle','r') as infile:
                            (W1_arr[jsample,jmu,:,:], p1_arr[jsample,jmu,:]) = pickle.load(infile)
                        logL[jsample, jmu] = np.loadtxt(infile_head + 'logL.txt')
                    else:
                        p_arr[jsample,jmu,:,:] = np.loadtxt(dir_name + '{}/sample{}_mu{}_p_a.txt'.format(name, jsample, jmu))
                    if name == 'GTR':
                        self.logL1[jsample, jmu] = np.loadtxt(infile_head + 'logL1.txt')
                        self.logL1exact[jsample, jmu] = np.loadtxt(infile_head + 'logLexact.txt')
                        
            self.p_a[name] = p_arr
            self.mu_a[name] = mu_arr
            self.Wij[name] = Wij_arr
            self.n_ij[name] = n_ij_arr
            self.T_i[name] = T_i_arr
            self.Tcons_i[name] = Tcons_i_arr
            self.W1_ij[name] = W1_arr
            self.p1_i[name] = p1_arr
            self.logL[name] = logL
    
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
            p_a = np.tile(MPD.p1_i[name][:,:,:,np.newaxis],(1,1,1,self.L))
            dist = np.sum((p_a - self.p0_a[np.newaxis,np.newaxis,:,:])**2, axis = (2,3))/self.L 
            means[jname,:], variances[jname,:] = mean_var(dist)
#            dist = np.sum((self.p1_i[name][:,:,:,np.newaxis] -\
#            self.p0_a[np.newaxis,np.newaxis,:,:])**2, axis = (2,3))/self.L 
#            means[jname,:], variances[jname,:] = mean_var(dist)
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
        
    def chi2_pi_plot(self, data_names, to_file, scale_by = 'branch_length', oneforall = False):
        xx, xlab, factor = self.scale_mu(scale_by)

        if oneforall:
            means, variances = self.chi2_pi1(data_names)
        else:
            means, variances = self.chi2_pi(data_names)
        plt.figure(30,figsize = (20,10)); plt.clf()
        plt.subplot(1,2,1)
        for jdata, mean_loc in enumerate(means):
            plt.errorbar(xx, np.log10(means[jdata,:]),\
            yerr = np.sqrt(variances[jdata,:])/means[jdata,:]/np.log(10))
        plt.xlabel(xlab, fontsize = 18)
        plt.ylabel(r'$\log_{10}\chi^2_{\pi}$', fontsize = 18)
        plt.legend(data_names, fontsize = 18, loc = 0)
    
        plt.subplot(1,2,2)
        for jdata, mean_loc in enumerate(means):
            plt.errorbar(np.log10(xx), np.log10(means[jdata,:]),\
            yerr = np.sqrt(variances[jdata,:])/means[jdata,:]/np.log(10))
        plt.xlabel(xlab, fontsize = 18)
        plt.ylabel(r'$\log_{10}\chi^2_{\pi}$', fontsize = 18)
        plt.legend(data_names, loc = 0, fontsize = 18)
        plt.savefig(to_file)
        plt.close(30)  
        return None
    
    def chi2_mu_plot(self, data_names, to_file, j0 = 0, scale_by = 'branch_length'):
        xx, xlab, factor = self.scale_mu(scale_by)

        plt.figure(30,figsize = (10,10)); plt.clf()
        for jname, name in enumerate([data_names[j] for j in [1,2,3]]):
            mu_adjusted = self.mu_a[name]*self.T_i[name].sum(axis=2)/self.Tcons_i[name].sum(axis=2)*\
            self.T_i[name].sum(axis=2)/self.T_i['GTR'].sum(axis=2)
            mu_dist = np.sqrt(((mu_adjusted - self.mu0_a[np.newaxis,:,:])**2).sum(axis = 2))/self.L
            dist_mean, dist_var = mean_var(mu_dist)        
            plt.errorbar(xx[j0:], dist_mean[j0:]/self.mmu[j0:],\
            yerr = np.sqrt(dist_var)[j0:]/self.mmu[j0:])
        plt.xlabel(xlab, fontsize = 18)
        plt.ylabel(r'$\chi^2_{\mu}/(2\mu)$', fontsize = 18)
        plt.legend(data_names[1:], fontsize = 18, loc = 0)
        plt.savefig(to_file)
        plt.close(30)
        return None

    def muave_mu_plot(self, data_names, to_file, j0 = 0, scale_by = 'branch_length'):
        xx, xlab, f = self.scale_mu(scale_by)
        if scale_by == 'dist_to_root':
            f0 = self.dist_to_root_mean
        elif scale_by == 'tree_length':
            f0 = self.T
        else:
            f0 = self.branch_lengths.mean()
            
        plt.figure(30,figsize = (10,10)); plt.clf()
        for jname, name in enumerate(data_names):
            mu_mean, mu_var = mean_var((self.mu_a[name][0,:,:]*\
            self.T_i[name][0,:,:,:].sum(axis=1)/self.Tcons_i[name][0,:,:,:].sum(axis=1)).T)
#            factor_tree = np.ones(self.mmu.shape[0])*f
#            factor_tree = f0*np.sum(self.Wij[name][0,:,:,:], axis = (1,2))/(self.q*(self.q-1))
#            mu_mean, mu_var = mean_var((self.mu_a[name][0,:,:]*\
#            self.T_i[name][0,:,:,:].sum(axis=1)/self.Tcons_i[name][0,:,:,:].sum(axis=1)).T)
            factor = f0*np.sum(self.Wij[name][0,:,:,:], axis = (1,2))/(self.q*(self.q-1))
            factor_tree = factor*self.T_i[name][0,:,:,:].sum(axis = (1,2))/self.T_i['GTR'][0,:,:,:].sum(axis = (1,2))
            plt.errorbar(xx[j0:], (mu_mean*factor_tree)[j0:], yerr = (np.sqrt(mu_var)*factor_tree)[j0:])
        plt.plot(xx,xx,':k')
        plt.xlabel(xlab, fontsize = 18)
        plt.ylabel(xlab, fontsize = 18)
        plt.legend(data_names, fontsize = 18, loc = 0)
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
        
    def Akaike_plot(self, data_names, to_file, scale_by = 'branch_length'):
        xx, xlab, f = self.scale_mu(scale_by)
        k = self.q*(self.q-1)//2 - 1 + self.q*self.L
        k1 = self.q*(self.q+1)//2
        
        leg = list(data_names)
        leg.extend(['GTR1','GTR1exact'])
        plt.figure(30,figsize = (20,10)); plt.clf()
        plt.subplot(1,2,1)
        for name in data_names:
            plt.plot(xx, self.logL[name].sum(axis=0))
        plt.plot(xx, self.logL1.sum(axis=0))
        plt.plot(xx, self.logL1exact.sum(axis=0))
        plt.xlabel(xlab, fontsize = 18)
        plt.ylabel(r'$\log{L}$', fontsize = 18)
        plt.legend(leg, fontsize = 18, loc = 0)
    
        plt.subplot(1,2,2)
        for name in data_names:
            plt.plot(xx, k - self.logL[name].sum(axis=0))
        plt.plot(xx, k1 - self.logL1.sum(axis=0))
        plt.plot(xx, k1 - self.logL1exact.sum(axis=0))
        plt.xlabel(xlab, fontsize = 18)
        plt.ylabel(r'$k - \log{L}$', fontsize = 18)
        plt.legend(leg, loc = 0, fontsize = 18)
        plt.savefig(to_file)
        plt.close(30)  
        return None
        
if __name__=="__main__":
    '''Generating sequences for a given tree using a GTR matrix'''
    
    plt.ioff()
    dir_name = '/ebio/ag-neher/share/users/vpuller/GTR_staggered/L100_test1/'
    outdir_name = dir_name + 'plots/'
    if not os.path.exists(outdir_name):
        os.makedirs(outdir_name)    

#    # loading data
#    data_names = ['aln','GTR','GTRanc','GTRanc_rescaled', 'GTRtree', 'GTRtree_unique', 'GTRtree_site'] 
    data_names = ['aln','aln_all','GTR','GTRanc','GTRanc_rescaled', 'GTRtree', 'GTRtree_unique'] 
#    data_names = ['aln','GTR','GTRanc','GTRtree'] 
    MPD = model_param_data(dir_name, 'model', data_names)
    
    # divergence in pi
    MPD.chi2_pi_plot(data_names, outdir_name + 'dist.pdf', scale_by = 'branch_length')
    MPD.chi2_pi_plot(data_names, outdir_name + 'dist1.pdf', scale_by = 'branch_length', oneforall = True)
    
    # site entropies
    MPD.site_entropy_plot(['aln','GTR','GTRanc', 'GTRtree'], outdir_name + 'entropy.pdf')
    
    # n_ij histograms
    MPD.nij_hist(['GTR', 'GTRanc', 'GTRanc_rescaled', 'GTRtree'], outdir_name + 'n_ij_sum.pdf')
    
    # divergence in mu
    MPD.chi2_mu_plot(data_names, outdir_name + 'mu_dist.pdf', j0 = 10, scale_by = 'branch_length')
      
    
    # average over sequence mutation rate
#    MPD.muave_mu_plot([data_names[j] for j in [1,2,3,4]], outdir_name + 'mu_mu0.pdf', j0 = 10, scale_by = 'branch_length')
    MPD.muave_mu_plot(['GTR'], outdir_name + 'mu_mu0.pdf', j0 = 10, scale_by = 'branch_length')


    # Akaike information criterion
    MPD.Akaike_plot(data_names[2:], outdir_name + 'Akaike.pdf', scale_by = 'branch_length' )
    
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

    
    xx, xlab, f = MPD.scale_mu('branch_length')
    Ttree = {}
    plt.figure(60); plt.clf()
    for name in data_names[1:]: 
        Ttree[name] = MPD.T_i[name][0,:,:,:].sum(axis=1).mean(axis=1)
        if name in ['GTR', 'GTRanc', 'GTRtree_site']:
            factor = np.zeros(MPD.mmu.shape)
            for jmu, mu in enumerate(MPD.mmu):
                Q0mean = np.mean(MPD.mu0_a[jmu,:]*MPD.p0_a[:,np.newaxis,:]*MPD.Wij0[:,:,np.newaxis], axis=2)
                factor[jmu] = Q0mean.sum()/Q0mean.shape[0]
            plt.plot(xx, Ttree[name]*factor)
        elif name in ['GTRanc_rescaled']:
            plt.plot(xx, Ttree[name]*.75)
        else:
            plt.plot(xx, Ttree[name])
    plt.legend(data_names[1:], loc = 0)
#    plt.legend(['GTRanc', 'GTRanc_rescaled', 'GTRtree'], loc = 0)
    plt.xlabel(xlab); plt.ylabel('T') 
    plt.savefig(outdir_name + 'tree_lengths.pdf')
    plt.close(60)
    
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