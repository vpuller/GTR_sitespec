# -*- coding: utf-8 -*-
"""
Created on Thu Apr 21 10:59:11 2016

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
    
    mmu = np.loadtxt(dir_name + 'mmu.txt')
    branch_lengths = np.loadtxt(dir_name + 'branch_lengths.txt')
    dist_to_root = np.loadtxt(dir_name + 'dist_to_root.txt')
    dist_to_root_mean = dist_to_root.mean()
    p0_a = np.loadtxt(dir_name + 'p0_a.txt')
    Wij0 = np.loadtxt(dir_name + 'Wij0.txt')
    
    data_names = ['tips','GTR','GTRanc','GTRtree']
    parallel = True
    if parallel:
        Nsample = np.loadtxt(dir_name + 'Nsample.txt', dtype = 'int')
        dist_GTR = np.zeros((Nsample,mmu.shape[0]))
        dist_freqs = np.zeros((Nsample,mmu.shape[0]))
        dist_GTR_anc = np.zeros((Nsample,mmu.shape[0]))
        dist_GTR_tree = np.zeros((Nsample,mmu.shape[0]))    
        for jsample in xrange(Nsample):
            for jmu, mu in enumerate(mmu):
                dist_freqs[jsample, jmu], dist_GTR[jsample,jmu],\
                dist_GTR_anc[jsample,jmu], dist_GTR_tree[jsample,jmu] =\
                np.loadtxt(dir_name + 'sample{}_mu{}.txt'.format(jsample, jmu))
#                os.remove(dir_name + 'sample{}_mu{}.txt'.format(jsample, jmu))   

        data_all = [dist_freqs, dist_GTR, dist_GTR_anc, dist_GTR_tree]
        for jdata, data in enumerate(data_all):
            np.savetxt(dir_name + data_names[jdata] + '_data.txt',data)
    else:
        data_all = []
        for jdata, data_name in enumerate(data_names):
            data_all.append(np.loadtxt(dir_name + data_names[jdata] + '_data.txt'))

    means = np.zeros((len(data_all),len(mmu)))
    variances = np.zeros((len(data_all),len(mmu)))
    for jdata, data in enumerate(data_all):
        means[jdata,:], variances[jdata,:] = mean_var(data)
    
    leg = ('Tips alignment', 'GTR', 'GTR ancestors', 'GTR tree and ancestors')    
    xx = mmu*dist_to_root_mean*Wij0[np.where(np.triu(Wij0) >0)].mean()
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
    plt.savefig(dir_name + 'dist.pdf')
    plt.close(30)
    
    
    LL = [100, 500, 1000]
    styles = ['-o','--s',':v']
    colors = ['b','g','r','c']
    plt.figure(40); plt.clf()
    for jL, L in enumerate(LL):
        dir_name_loc = '/ebio/ag-neher/share/users/vpuller/GTR_staggered/L{}/'.format(L) 
        mmu = np.loadtxt(dir_name_loc + 'mmu.txt')
#        branch_lengths = np.loadtxt(dir_name_loc + 'branch_lengths.txt')
        dist_to_root = np.loadtxt(dir_name_loc + 'dist_to_root.txt')
        dist_to_root_mean = dist_to_root.mean()
    
        data_names = ['tips','GTR','GTRanc','GTRtree']
        data_all = []
        for jdata, data_name in enumerate(data_names):
            data_all.append(np.loadtxt(dir_name_loc + data_names[jdata] + '_data.txt'))
    
        means = np.zeros((len(data_all),len(mmu)))
        variances = np.zeros((len(data_all),len(mmu)))
        for jdata, data in enumerate(data_all):
            means[jdata,:], variances[jdata,:] = mean_var(data/L)
        
        for jdata, mean_loc in enumerate(means):
            plt.errorbar(xx, np.log10(means[jdata,:]),\
            yerr = np.sqrt(variances[jdata,:])/means[jdata,:]/np.log(10), fmt = styles[jL] + colors[jdata])
            plt.xlabel(r'$\mu\bar{t}_{root}$', fontsize = 18)
            plt.ylabel(r'$\log_{10}(\chi^2/L)$', fontsize = 18)
            plt.legend(leg, fontsize = 18, loc = 0)
    plt.savefig(dir_name + 'distL.pdf')
    plt.close(40)