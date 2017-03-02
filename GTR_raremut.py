# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 16:29:49 2017

@author: vpuller
"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from Bio import Phylo
import os, sys
import create_model as createGTR

import dress_tree as DressTree
import GTR_inference as GTRinf
#sys.path.append('/ebio/ag-neher/share/users/vpuller/CLUSTER/GTR_bootstrap')
#import dress_tree 

np.random.seed(seed=134)

#Constants
h= 10**(-8)
Hw = 8
fs = 24
cols = ['b', 'g','r','c', 'm']

def nuc_freq(aln_array, alphabet = 'ACGT'):
    '''calculating frequency of different nucleotides in array'''
    return np.array([(aln_array == nuc).mean(axis=0) for nuc in alphabet])


def mu_bayes(n_ija, T_ija, W_ij, p_a):
#    T_ja = T_ija[range(q), range(q),:]
    T_ja = np.sum(T_ija + np.transpose(T_ija, axes = (1,0,2)), axis = 1)/2.
    bbeta = np.sum(p0_a*W_ij.dot(T_ja), axis =0)
    kk = np.sum(n_ija, axis = (0,1))
    return kk/bbeta
    
if __name__=="__main__":
    '''GTR reconstruction: dealing with rare mutations'''
    plt.ioff()
    plt.close('all')
    
    # create GTR model --> distribution of mu_a, pi_ia
    # evolve sequences
    # reconstruct GTR model --> distribution of mu_a, pi_ia
    
    alphabet = 'ACGT'
    q = len(alphabet)
    L = 900 # sequence length
    Nsample = 3
    outdir_name = '/ebio/ag-neher/share/users/vpuller/GTR_staggered/Rare_mutations/'
    if not os.path.exists(outdir_name):
        os.makedirs(outdir_name)    
    
    tree_file_name = '/ebio/ag-neher/share/users/vpuller/GTR_staggered/trees/test_tree500.nwk'
    bio_tree = Phylo.read(tree_file_name, 'newick')
    for clade in bio_tree.find_clades():
        clade.branch_length = 1.
    Phylo.write(bio_tree, tree_file_name[:-4] + '_unit.nwk', 'newick')

    tree_file_name = '/ebio/ag-neher/share/users/vpuller/GTR_staggered/trees/test_tree500_unit.nwk'
    bio_tree = Phylo.read(tree_file_name, 'newick')
    bio_tree.root.branch_length = h
    branch_lengths = np.array([clade.branch_length for clade in bio_tree.find_clades()])
#    dist_to_root = np.array([bio_tree.distance(clade) for clade in bio_tree.get_terminals()])

    # Creating and saving GTR model
#    model_spec = {'mu': ('gamma', 1., 5),
#                  'pi': 'simplex',
#                  'Wij': ('JC',),
#                  'seq0': 'equilibrated'} 
    model_spec = {'mu': ('uniform', 1., 5),
                  'pi': 'equal',
                  'Wij': ('JC',),
                  'seq0': 'equilibrated'} 
    if not os.path.exists(outdir_name + 'model/'):
        os.makedirs(outdir_name + 'model/')   
    Wij0, p0_a, mu0_a, seq0 = createGTR.create_GTR_model(L, model_spec,\
    model_file_head = outdir_name + 'model/')
    
    Hp0 = - np.sum(p0_a*np.log(p0_a+h)/np.log(q), axis=0)
    pseq = np.array([seq0.count(nuc) for nuc in alphabet], dtype = 'float')
    pseq = pseq/np.sum(pseq)
    Hseq = - np.sum(pseq*np.log(pseq+h)/np.log(q))
    
    
    # simulating tree
    dress = DressTree.dress_tree(tree_file_name, seq0, mu_a = mu0_a, Wij = Wij0,\
    p0_a = p0_a, alphabet = alphabet, use_eigenvalues = True)
    
    # leaves aln
    arr = np.array([list(clade.seq) for clade in dress.tree.get_terminals()])
    paln = nuc_freq(arr, alphabet = alphabet)
    Hpaln = - np.sum(paln*np.log(paln+h)/np.log(q), axis=0)


    # counts and GTR reconstruction
    n_ij, T_ij, root_states = GTRinf.nuc_sub_count(dress.tree, alphabet = alphabet)   
    W_ij, p_a, mu_a = GTRinf.GTR_simult(n_ij,T_ij,root_states) 
    Hp = - np.sum(p_a*np.log(p_a+h)/np.log(q), axis=0)
    
#    mmu_th0 = mu_bayes(n_ij, T_ij, Wij0, p0_a)
    mmu_th1 = mu_bayes(n_ij, T_ij, W_ij, p_a)
    kk = np.sum(n_ij, axis = (0,1))
    mmu_unit = kk*(q**2)/branch_lengths.sum()
#    T_ja = T_ij[range(q), range(q),:]
    T_ja = np.sum(T_ij + np.transpose(T_ij, axes = (1,0,2)), axis = 1)/2.
    bbeta = np.sum(p0_a*W_ij.dot(T_ja), axis =0)
    mmu_th0 = kk/bbeta
    
    fig, ax = plt.subplots(1,4, figsize = (4*Hw, Hw))
    for jmu, mu in enumerate([mu0_a, mu_a, mmu_th0, mmu_th1, mmu_unit]):
        nn, edges = np.histogram(mu)
        ax[0].plot((edges[:-1] + edges[1:])/2., nn, label = str(jmu))
        ax[0].axvline(mu.mean(), lw=2, color=cols[jmu])
    ax[0].set_xlabel('mu_a', fontsize = .8*fs)
    ax[0].legend(fontsize = fs, loc = 0)

    for H in [Hp0, Hpaln, Hp]:
        nn, edges = np.histogram(H)
        ax[1].plot((edges[:-1] + edges[1:])/2., nn)
    ax[1].set_xlabel('Hp', fontsize = .8*fs)
    ax[1].legend(['model', 'aln', 'GTR'], fontsize = fs, loc = 0)

    nsum = n_ij.sum(axis=(0,1), dtype = 'int')
    nn = np.array([np.count_nonzero(nsum == j) for j in xrange(np.max(nsum)+1)])
    ax[2].plot(range(np.max(nsum)+1), nn)
    ax[2].set_xlabel('# mutations', fontsize = .8*fs)
    
    nn, edges = np.histogram(branch_lengths, 20)
    ax[3].plot((edges[:-1] + edges[1:])/2., nn)
    ax[3].set_xlabel('branch length', fontsize = .8*fs)
    plt.savefig(outdir_name + 'hist1.pdf')
    plt.close()
    
##    ncr = 1
##    jj = np.where(nsum > ncr)[0]
##    a = np.array([np.count_nonzero(n_ij[:,:,j]) for j in xrange(L)])
##    [np.count_nonzero(a == k) for k in xrange(13)]
#    m_i = root_states + n_ij.sum(axis=1)
##    [np.count_nonzero((m_i > 0).sum(axis=0) == j) for j in xrange(5)]
#    jj = np.where((m_i>0).sum(axis=0) > 1)[0]
#    W_ij, p_a, mu_a = GTRinf.GTR_simult(n_ij[:,:,jj],T_ij[:,:,jj],root_states[:,jj]) 
#    Hp = - np.sum(p_a*np.log(p_a+h)/np.log(q), axis=0)
#    
#    fig, ax = plt.subplots(1,3, figsize = (3*Hw, Hw))
#    for mu in [mu0_a[jj], mu_a]:
#        nn, edges = np.histogram(mu)
#        ax[0].plot((edges[:-1] + edges[1:])/2., nn)
#    ax[0].set_xlabel('mu_a', fontsize = .8*fs)
#    ax[0].legend(['model', 'GTR'], fontsize = fs, loc = 0)
#
#    for H in [Hp0[jj], Hpaln[jj], Hp]:
#        nn, edges = np.histogram(H)
#        ax[1].plot((edges[:-1] + edges[1:])/2., nn)
#    ax[1].set_xlabel('Hp', fontsize = .8*fs)
#    ax[1].legend(['model', 'aln', 'GTR'], fontsize = fs, loc = 0)
#
#    nsum = n_ij.sum(axis=(0,1), dtype = 'int')
#    nn = np.array([np.count_nonzero(nsum == j) for j in xrange(np.max(nsum)+1)])
#    ax[2].plot(range(np.max(nsum)+1), nn)
#    ax[2].set_xlabel('# mutations', fontsize = .8*fs)
#    plt.savefig(outdir_name + 'hist_jj.pdf')
#    plt.close()
