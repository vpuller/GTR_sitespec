# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 13:33:49 2016

@author: vpuller
"""
from __future__ import division
import numpy as np
#from scipy import linalg as LA
#from scipy import misc
from Bio import Phylo #, AlignIO
import matplotlib.pyplot as plt
#from scipy import optimize
import sys, os
import GTR_twoseq
import generate_tree_GTR as gentree
#import copy as copy
import time

#sys.path.append('/ebio/ag-neher/share/users/vpuller/myTOOLS/') 
#import Vadim_toolbox_file as vp
#import GTR_class

#Constants
h= 10**(-8)

def GTR_wcat(n_ij, T_i, root_states, idx0, idx_cat):
    '''
    Solve GTR model with specified categories
    
    Input arguments:
    n_ij, T_i, root_i - transition counts, times in states and root states for every site
    idx0 - indices of sites not to be categized
    idx_cat - indices for sites to be categorized
    '''
    ncat_ij = n_ij[:,:,idx0]
    Tcat_i = T_i[:,idx0]
    root_states_cat = root_states[:,idx0]
    for idx in idx_cat:
        ncat_ij = np.append(ncat_ij,n_ij[:,:,idx].sum(axis=2, keepdims = True),axis=2)
        Tcat_i = np.append(Tcat_i,T_i[:,idx].sum(axis=1, keepdims = True),axis=1)
        root_states_cat = np.append(root_states_cat,root_states[:,idx].sum(axis=1, keepdims = True),axis=1)
    
    Wcat_ij, pcat, mucat = gentree.GTR_simult(ncat_ij,Tcat_i,root_states_cat) 
    pcat_a = np.zeros(p_a.shape)
    pcat_a[:,idx0] = pcat[:,:len(idx0)]
    mucat_a = np.zeros(mu_a.shape)
    mucat_a[idx0] = mucat[:len(idx0)]
    for jcat, idx in enumerate(idx_cat):
        pcat_a[:,idx] = np.array([pcat[:,len(idx0)+ jcat]]).T
        mucat_a[idx] = mucat[len(idx0)+ jcat]
    
    return Wcat_ij, pcat_a, mucat_a

def Tcategories(T_ia, Nbins):
    '''Split residence times into categories
    
    T_ia - residence times for nucs for every site
    Nbins - number of time bins
    '''
    cats = []
    Ttot = np.max(T_ia.sum(axis=0))
    q = T_ia.shape[0]
    dig = np.array([np.digitize(T_ia[jnuc,:],np.linspace(0.,Ttot+h,num = Nbins + 1)) for jnuc in xrange(q)])
    for jjbin in set([tuple(row) for row in dig.T]):
        idx = np.where(np.array([dig[jnuc,:] == jjbin[jnuc] for jnuc in xrange(q)]).sum(axis=0) == q)[0]
        cats.append(idx_q[idx])
    return cats
    
if __name__=="__main__":
    '''Clustering sites to improve GTR performance'''
    
    plt.ioff()
    outdir_name = '/ebio/ag-neher/share/users/vpuller/GTR_staggered/tmp_new/'
    if not os.path.exists(outdir_name):
        os.makedirs(outdir_name)    
    
    tree_file_name = '/ebio/ag-neher/share/users/vpuller/Maximum_likelihood/ML_test_sequences/RRE/HIV1B500_RRE_nuc_tree.nwk'    
    bio_tree = Phylo.read(tree_file_name, 'newick')
    bio_tree.root.branch_length = h
    Ttot = bio_tree.total_branch_length()
    
    # generate root sequence
    alphabet = 'ACGT'
    q = len(alphabet)
    L = 2*10**2 # sequence length
    seq0 = GTR_twoseq.random_seq_p0(L,alphabet = alphabet)
    P0 = gentree.seq_to_P(seq0, alphabet = alphabet)
    
    #dressing tree with sequences
    mu0 = 1.
    mu0_a = mu0*np.ones(L)
    W0_ij = q*(np.ones((q,q)) - np.eye(q))
#    p0_a = np.ones((q,L))/q
    p0_a = np.random.exponential(size = (q,L)); p0_a = p0_a/np.sum(p0_a,axis=0)
    Snuc0 = - ((p0_a + h)*np.log(p0_a + h)).sum(axis = 0)
    
    t0 = time.time()
    dress = gentree.dress_tree(bio_tree, seq0, mu_a = mu0_a, Wij = W0_ij, p0_a = p0_a, alphabet = alphabet)
    t1 = time.time()
    print 't = ', t1-t0

    
    # array of all sequences on the tree
#    arr = np.array([list(clade.seq) for clade in dress.tree.find_clades()])
    arr = np.array([list(clade.seq) for clade in dress.tree.find_clades()])
    freqs = gentree.nuc_freq(arr, alphabet = alphabet)
    Snuc = - ((freqs + h)*np.log(freqs + h)).sum(axis = 0)    
    
    # GTR fitting
    n_ij, T_i, root_states = gentree.nuc_sub_count(dress.tree, alphabet = alphabet)   
    W_ij, p_a, mu_a = gentree.GTR_simult(n_ij,T_i,root_states) 
    S_a = -((h+p_a)*np.log(h+p_a)).sum(axis=0)
    
    dp_a = p_a/np.sqrt(n_ij.sum(axis=1))
    dp_a[np.where(np.sqrt(n_ij.sum(axis=1)) == 0)] = 1.
#    dS_a = -(dp_a*np.log(h+p_a)).sum(axis=0)
    
    dist_freqs = np.sum((freqs - p0_a)**2)            
    dist_GTR = np.sum((p_a - p0_a)**2)
    print 'distance to alignment = {}\ndistance to GTR = {}'.format(dist_freqs,dist_GTR)     

    ncr = 9
    
#    idx0 = np.where((n_ij.sum(axis=1) < ncr).sum(axis=0) < q)[0]
#    idx_q = np.where((n_ij.sum(axis=1) < ncr).sum(axis=0) == q)[0]
#    ncat_ij = np.zeros((q,q,len(idx0) + 1))
#    ncat_ij = n_ij[:,:,idx0]
#    ncat_ij = np.append(ncat_ij,n_ij[:,:,idx_q].sum(axis=2, keepdims = True),axis=2)
#    Tcat_i = T_i[:,idx0]
#    Tcat_i = np.append(Tcat_i,T_i[:,idx_q].sum(axis=1, keepdims = True),axis=1)
#    root_states_cat = root_states[:,idx0]
#    root_states_cat = np.append(root_states_cat,root_states[:,idx_q].sum(axis=1, keepdims = True),axis=1)
    
    idx0 = np.where((n_ij.sum(axis=1) < ncr).sum(axis=0) < q)[0]
    idx_q = np.where((n_ij.sum(axis=1) < ncr).sum(axis=0) == q)[0]
    Nbins = 10
#    idx_cat = []
#    dig = np.array([np.digitize(T_i[jnuc,idx_q],np.linspace(0.,Ttot+h,num = Nbins + 1)) for jnuc in xrange(q)])
#    for jjbin in set([tuple(row) for row in dig.T]):
#        idx = np.where(np.array([dig[jnuc,:] == jjbin[jnuc] for jnuc in xrange(q)]).sum(axis=0) == q)[0]
#        idx_cat.append(idx_q[idx])
    
    idx_cat = Tcategories(T_i[:,idx_q], Nbins)
#    if q == 2:
#        idx_cat = []
#        nnbins, Tbins = np.histogram(T_i[0,idx_q], bins = Nbins)
#        dig = np.digitize(T_i[0,idx_q],np.linspace(0.,Ttot+h,num = Nbins + 1))
#        idx_cat.extend([idx_q[np.where(dig == jbin + 1)[0]] for jbin in xrange(Nbins) if np.count_nonzero(dig == jbin + 1) > 0])
#    elif q == 4:
#        idx_cat = []
#        nnbins, Tbins = np.histogram(T_i[0,idx_q], bins = Nbins)

    
#    ncat_ij = n_ij[:,:,idx_cat[0]]
#    Tcat_i = T_i[:,idx_cat[0]]
#    root_states_cat = root_states[:,idx_cat[0]]
#    for idx in idx_cat[1:]:
#        ncat_ij = np.append(ncat_ij,n_ij[:,:,idx].sum(axis=2, keepdims = True),axis=2)
#        Tcat_i = np.append(Tcat_i,T_i[:,idx].sum(axis=1, keepdims = True),axis=1)
#        root_states_cat = np.append(root_states_cat,root_states[:,idx].sum(axis=1, keepdims = True),axis=1)
#    
#    Wcat_ij, pcat, mucat = gentree.GTR_simult(ncat_ij,Tcat_i,root_states_cat) 
#    pcat_a = np.zeros(p_a.shape)
#    pcat_a[:,idx_cat[0]] = pcat[:,:len(idx_cat[0])]
#    mucat_a = np.zeros(mu_a.shape)
#    mucat_a[idx_cat[0]] = mucat[:len(idx_cat[0])]
#    for jcat, idx in enumerate(idx_cat[1:]):
#        pcat_a[:,idx] = np.array([pcat[:,len(idx_cat[0])+ jcat]]).T
#        mucat_a[idx] = mucat[len(idx_cat[0])+ jcat]
    
    Wcat_ij, pcat_a, mucat_a = GTR_wcat(n_ij,T_i,root_states, idx0, idx_cat) 
    
    Scat_a = -((h+pcat_a)*np.log(h+pcat_a)).sum(axis=0)
    dist_GTRcat = np.sum((pcat_a - p0_a)**2)   
    print 'distance to GTR cat. = {}'.format(dist_GTRcat) 
#    dist_GTRcat1 = np.sum((pcat_a1 - p0_a)**2)   
#    print 'distance to GTR cat1. = {}'.format(dist_GTRcat1) 
    
    Like0 = gentree.dress_tree(dress.tree, seq0, mu_a = mu0_a, Wij = W0_ij, p0_a = p0_a, alphabet = alphabet, logL = True)
    Like_GTR = gentree.dress_tree(dress.tree, seq0, mu_a = mu_a, Wij = W_ij, p0_a = p_a, alphabet = alphabet, logL = True)
    Like_GTRcat = gentree.dress_tree(dress.tree, seq0, mu_a = mucat_a, Wij = Wcat_ij, p0_a = pcat_a, alphabet = alphabet, logL = True)
    print 'logL0 = {}\nlogL_GTR = {}\nlogL_GTRcat = {}'.format(Like0.logL, Like_GTR.logL, Like_GTRcat.logL) 
    
    print 'AIC_GTR = {}\nAIC_GTRcat = {}'.format(q*L - Like_GTR.logL, q*(len(idx0) + len(idx_cat)) - Like_GTRcat.logL)
    #plots
    plt.figure(10,figsize = (20,6*(q+1))); plt.clf()
    for jnuc, nuc in enumerate(alphabet):
        plt.subplot(q+1,1,jnuc+1)
        plt.plot(p0_a[jnuc,:])
        plt.plot(freqs[jnuc,:])
        plt.errorbar(np.arange(L),p_a[jnuc,:],yerr = dp_a[jnuc,:])
        plt.xlabel('site'); plt.ylabel('p_' + nuc)
        plt.legend(('model','tree seq.','GTR'))
    plt.subplot(q+1,1,q+1)
    plt.plot(Snuc0/np.log(q))
    plt.plot(Snuc/np.log(q))
    plt.plot(S_a/np.log(q))
    plt.plot(Scat_a/np.log(q))
    for jcat, idx in enumerate(idx_cat[1:]):
        plt.scatter(idx,.95*np.ones(len(idx)))
#    plt.scatter(idx_q,.95*np.ones(len(idx_q)))
    plt.xlim(0,L); plt.ylim(0.,1.)
    plt.xlabel('site'); plt.ylabel('entropy')
    plt.legend(('model','tree seq.','GTR','GTR, cat.'),loc = 0)
    plt.savefig(outdir_name + 'freqs.pdf')
    plt.close(10)
    
    
    plt.figure(20, figsize = (6*q,6)); plt.clf()
    for jnuc, nuc in enumerate(alphabet):
        plt.subplot(1,q,jnuc+1)
        plt.hist(T_i[jnuc,:]/Ttot, bins = 20)
        plt.xlabel(r'$T_'+nuc +'$'); plt.title(nuc)
    plt.savefig(outdir_name + 'Times_hist.pdf')
    plt.close(20)
    
    plt.figure(20, figsize = (6*q,6)); plt.clf()
    for jnuc, nuc in enumerate(alphabet):
        plt.subplot(1,q,jnuc+1)
        idx = np.where(n_ij[jnuc,:,:].sum(axis = 0) < ncr)[0]
        for jnuc0, nuc0 in enumerate(alphabet):   
            plt.hist(T_i[jnuc0,idx]/Ttot, bins = 20,alpha = 0.5)
        plt.xlabel(r'$T$'); plt.title('-->' + nuc)
        plt.legend([r'$T_' + nuc0 + '$' for nuc0 in alphabet],loc = 0)
    plt.savefig(outdir_name + 'Times_hist1.pdf')
    plt.close(20)
    
#    #categorizing sites by the number of nuc. states with insufficient counts 
#    idx_cat = []
#    for jcat in xrange(q + 1):
#        idx_cat.append(np.where((n_ij.sum(axis=1) < ncr).sum(axis=0) == jcat)[0])    
    
#    # categorizing by the number and identity of nuc. states with insufficient counts
#    idx_cat = []
#    titles_cat = []
#    for jcat in xrange(q+1):
#        if jcat == 0 or jcat == q:
#            idx_cat.append(np.where((n_ij.sum(axis=1) < ncr).sum(axis=0) == jcat)[0])
#            titles_cat.append('n = {}'.format(jcat))
#        elif jcat == 1:
#            for jnuc, nuc in enumerate(alphabet):
#                idx_cat.append(np.where(((n_ij.sum(axis=1) < ncr).sum(axis=0) == jcat)*\
#                (n_ij[jnuc,:,:].sum(axis=0) < ncr))[0])
#                titles_cat.append('n = {}, nuc = {}'.format(jcat,nuc))
#        else:
#            for jnuc in xrange(misc.comb(q,jcat)):
#                idx_cat.append(np.where((n_ij.sum(axis=1) < ncr).sum(axis=0) == jcat)[0])
#                titles_cat.append('n = {}'.format(jcat))
#                
#    plt.figure(20, figsize = (6*len(idx_cat),6)); plt.clf()
#    for jcat, idx in enumerate(idx_cat):
#        plt.subplot(1,len(idx_cat),jcat+1)
##        idx = np.where((n_ij.sum(axis=1) < ncr).sum(axis=0) == jcat)[0]
#        if len(idx > 0):
#            for jnuc0, nuc0 in enumerate(alphabet):   
#                plt.hist(T_i[jnuc0,idx]/Ttot, bins = 20,alpha = 0.5)
#        plt.xlabel(r'$T$'); plt.title(titles_cat[jcat])
#        plt.legend([r'$T_' + nuc0 + '$' for nuc0 in alphabet],loc = 0)
#    plt.savefig(outdir_name + 'Times_bycategory_hist.pdf')
#    plt.close(20)
    
#    plt.figure(20, figsize = (6*(q+1),6)); plt.clf()
#    for jnuc, nuc in enumerate(alphabet):
#        plt.subplot(1,q+1,jnuc+1)
#        plt.hist(freqs[jnuc,:])
#        plt.xlabel('frequency'); plt.title(nuc)
#    plt.subplot(1,q+1,q+1)
#    plt.hist(Snuc/np.log(q))
#    plt.xlabel('entropy'); plt.title('Snuc')
#    plt.savefig(outdir_name + 'freq_hist.pdf')
#    plt.close(20)
#    
#    plt.figure(20, figsize = (6*(q+1),6)); plt.clf()
#    for jnuc, nuc in enumerate(alphabet):
#        plt.subplot(1,q+1,jnuc+1)
#        plt.hist(p_a[jnuc,:])
#        plt.xlabel('fraction of time'); plt.title(nuc)
#    plt.subplot(1,q+1,q+1)
#    plt.hist(S_a/np.log(q))
#    plt.xlabel('entropy'); plt.title('Snuc')
#    plt.savefig(outdir_name + 'freqGTR_hist.pdf')
#    plt.close(20)
    