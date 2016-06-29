# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 12:32:36 2016

@author: vpuller
"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import pickle as pickle

#Constants
h= 10**(-8)
alphabet = 'ACGT'

def site_entropy(freqs):
    return -np.sum((freqs + h)*np.log(freqs + h), axis = 0)
    
if __name__=="__main__":
    '''GTR reconstruction for HIV data: plots'''
    
    plt.ioff()
    plt.close('all')
    dir_name = '/ebio/ag-neher/share/users/vpuller/GTR_staggered/HIV_test/'
    gene = 'pol'
    
    freqs = np.loadtxt(dir_name + '{}_aln_p_a.txt'.format(gene))
    S_a = site_entropy(freqs)
    
    data_names = ['aln','GTRanc_pavel', 'GTRanc_vadim']
#    with open(dir_name + '{}_GTRanc_pavel/'.format(gene) + 'counts.pickle','r') as infile:
#        (n2_ij, T2_i, Tcons2_i, root_states2) = pickle.load(infile)
    with open(dir_name + '{}_GTRanc_pavel/'.format(gene) + 'GTRinf.pickle','r') as infile:
        (W2_ij, p2_a, mu2_a) = pickle.load(infile)
    branch_lengths2 = np.loadtxt(dir_name + '{}_GTRanc_pavel/'.format(gene) + 'branch_lengths.txt')
    
#    with open(dir_name + '{}_GTRanc_vadim/'.format(gene) + 'GTRinf.pickle','r') as infile:
#        (W6_ij, p6_a, mu6_a) = pickle.load(infile)
#    branch_lengths6 = np.loadtxt(dir_name + '{}_GTRanc_vadim/'.format(gene) + 'branch_lengths.txt')
    
#    np.savetxt(output_file_head + 'dist_to_root.txt', dist_to_root)  
#    np.savetxt(output_file_head + 'logL.txt',np.array([logL]))
#    if oneforall:
#        with open(output_file_head + 'GTRone.pickle','w') as outfile:
#            pickle.dump((W1_ij, p1_i), outfile)
#    if oneforall and sitespec:
#        np.savetxt(output_file_head + 'Likelihoods.txt', np.array([logL,logL_one,logL_sitespec]))
    
    S2_a = site_entropy(p2_a)
    plt.figure(10, figsize = (18,6)); plt.clf()
    plt.subplot(1,3,1)
    plt.hist(S_a/np.log(len(alphabet)))
    plt.xlabel(r'$S_{\alpha}$')
    plt.title('Alighnment')

    plt.subplot(1,3,2)
    plt.hist(S2_a/np.log(len(alphabet)))
    plt.xlabel(r'$S_{\alpha}$')
    plt.title('GTR')

    plt.subplot(1,3,3)
    plt.hist(mu2_a*branch_lengths2.mean())
    plt.xlabel(r'$\mu_{\alpha}\bar{t}_{branch}$' )
    plt.title('GTR')
    plt.savefig(dir_name + '{}_S_mu.pdf'.format(gene))
    plt.close(10)
    
    plt.figure(20, figsize = (6*len(alphabet), 6)); plt.clf()
    for jnuc, nuc in enumerate(alphabet):
        plt.subplot(1, len(alphabet), jnuc + 1)
        plt.scatter(freqs[jnuc,:], p2_a[jnuc,:])
        plt.ylabel(r'$\pi_{\alpha}$, GTR')
        plt.xlabel(r'$\pi_{\alpha}$, aln')
        plt.title(nuc)
    plt.savefig(dir_name + '{}_freqs.pdf'.format(gene))
    plt.close(20)
    
#    S2_a = site_entropy(p2_a)
#    S6_a = site_entropy(p6_a)
#    entropies = [S_a, S2_a, S6_a]
#    plt.figure(10, figsize = (6*len(entropies),6)); plt.clf()
#    for j, S in enumerate(entropies):
#        plt.subplot(1, 3, j + 1)
#        plt.hist(S/np.log(len(alphabet)))
#        plt.xlabel('S_a')
#        plt.title(data_names[j])
#    plt.savefig(dir_name + '{}_S_a.pdf'.format(gene))
#    plt.close(10)
#    
#    plt.figure(10, figsize = (12,6)); plt.clf()
#    plt.subplot(1, 2, 1)
#    plt.hist(mu2_a*branch_lengths2.mean())
#    plt.xlabel(r'$\mu\bar{t}_{branch}$' )
#    plt.title('pavel')
#    plt.subplot(1, 2, 2)
#    plt.hist(mu6_a*branch_lengths6.mean())
#    plt.xlabel(r'$\mu\bar{t}_{branch}$' )
#    plt.title('vadim')
#    plt.savefig(dir_name + '{}_mu_a.pdf'.format(gene))
#    plt.close(10)
    