# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 18:09:38 2016

@author: vpuller
"""

from __future__ import division
import numpy as np
from scipy import linalg as LA
from Bio import AlignIO, Phylo
import matplotlib.pyplot as plt
from scipy import optimize
import sys, os
import GTR_twoseq

sys.path.append('/ebio/ag-neher/share/users/vpuller/myTOOLS/') 
import Vadim_toolbox_file as vp
import GTR_class

#Constants
h= 10**(-8)
#
#def assemble_Qij(mu_a, Wij, p0_a):
#    q = len(alphabet); L = len(seq)
#    #defining mutation matrix
#    if mu_a is None:
#        mu_a = np.ones(L)
#    
#    if Wij is None:
#        Wij = q*(np.ones((q,q)) - np.eye(q))
#        
#    if p0_a is None:
#       p0_a = np.ones((q,L))

def seq_to_P(seq,alphabet = 'ACGT'):
    q = len(alphabet); L = len(seq)
    seq_num=np.zeros(len(seq),dtype = 'int')
    seq_arr = np.array(list(seq))
    for jbase, base in enumerate(alphabet):
        seq_num[np.where(seq_arr == base)] = jbase   
#    seq_num = vp.string_to_numbers(seq)
    P = np.zeros((q,L),dtype = 'int')
    P[seq_num,range(L)] = 1.
    return P

def P_to_seq(P,alphabet = 'ACGT'):
    seq_num = np.zeros(P.shape[1],dtype = 'int')
    seq_num[np.where(P > 0)[1]] = np.where(P >0)[0]
    return ''.join(vp.numbers_to_string(seq_num,alphabet = alphabet))
    
#def evolve_seq_Qija(seq, t, Qija, alphabet = 'ACGT'):
#    '''
#    Evolve sequence using site - specific GTR model
#    
#    Input parameters:
#    seq - sequence
#    t - time
#    Qij - site-specific GTR matrices
#    '''
#    q = len(alphabet); L = len(seq)
#    
#    # converting sequence to vector form
#    P0 = seq_to_P(seq,alphabet = alphabet)
#    
#    # evolving sequence
#    Pt = np.array([LA.expm(Qija[:,:,ja]*t).dot(P0[:,ja]) for ja in xrange(L)]).T
##    Pt = Qija.dot(P0)[range(L),:,range(L)].T
#    
#    # measuring sequence
#    MC = np.random.uniform(size = L)
#    P1 = np.zeros((q,L))
#    for jq in xrange(q):
#        P1[jq,np.where((Pt[:jq,:].sum(axis=0) < MC)*(MC < Pt[:jq+1,:].sum(axis=0)))[0]] = 1.
#    
#    # converting result to sequence form
#    seq1_num = np.zeros(L,dtype = 'int')
#    seq1_num[np.where(P1 > 0)[1]] = np.where(P1 >0)[0]
#    return ''.join(vp.numbers_to_nucs(seq1_num))


def evolve_seq_Qija_P(P0, t, Qija):
    '''
    Evolve sequence using site - specific GTR model
    '''
    (q,L) = P0.shape
    
    # evolving sequence
    Pt = np.array([LA.expm(Qija[:,:,ja]*t).dot(P0[:,ja]) for ja in xrange(L)]).T
#    Pt = Qija.dot(P0)[range(L),:,range(L)].T
    
    # measuring sequence
    MC = np.random.uniform(size = L)
    P1 = np.zeros((q,L))
    for jq in xrange(q):
        P1[jq,np.where((Pt[:jq,:].sum(axis=0) < MC)*(MC < Pt[:jq+1,:].sum(axis=0)))[0]] = 1.
    
    return P1
    
def dress_tree_seq(tree, seq_root, mu_a = None, Wij = None, p0_a = None, alphabet = 'ACGT'):
    '''
    Dress tree withs equences evolved according to the specified GTR model
    '''
    q = len(alphabet); L = len(seq_root)
    #defining mutation matrix
    if mu_a is None:
        mu_a = np.ones(L)
    
    if Wij is None:
        Wij = q*(np.ones((q,q)) - np.eye(q))
        
    if p0_a is None:
       p0_a = np.ones((q,L))
       
#    mu_a = np.ones(L)
#    Wij = q*(np.ones((q,q)) - np.eye(q))
#    p0_a = np.ones((q,L))/q
    Qija = np.array([mu_a[ja]*np.diag(p0_a[:,ja]).dot(Wij) - 
    np.diag(np.sum(mu_a[ja]*np.diag(p0_a[:,ja]).dot(Wij), axis = 0)) for ja in xrange(L)])
    Qija = np.transpose(Qija,(1,2,0))
    
#    n = 0
    tree.root.seq = seq0
#    print seq0
    tree.root.P = seq_to_P(seq0,alphabet = alphabet)
    for parent in tree.get_nonterminals():
        for child in parent:
#            print n
#            n +=1
#            child.seq = GTR_twoseq.evolve_seq_GTR(parent.seq,child.branch_length, Qij = Qij)
#            child.P = seq_to_P(child.seq)
#            child.seq = evolve_seq_Qija(clade.seq,child.branch_length, Qija)
            child.P = evolve_seq_Qija_P(parent.P,child.branch_length, Qija)
            child.seq = P_to_seq(child.P,alphabet = alphabet)
#            if parent == tree.root:
#                print child.seq
            
    return tree
    
def nuc_sub_count(tree,alphabet = 'ACGT', fpar = 0.5):        
        '''
        Function that counts the numbers of nucleotide substitutions on tree
        
        Input arguments:
        tree - tree dressed withs equences
        alphabet - nucleotide alphabet
        fpar - fraction of the branch length assigned to the parent nucletide state
        '''
        print alphabet
        q = len(alphabet); L = len(tree.root.seq)
        substitutions = np.zeros((q,q,L),dtype = 'int')
        T_i=np.zeros((q,L))
        for parent in tree.get_nonterminals():
            seq_parent_num=vp.string_to_numbers(parent.seq,alphabet = alphabet)
            for child in parent:
                if child.branch_length >0:
                    seq_child_num=vp.string_to_numbers(child.seq,alphabet = alphabet)
                    T_i[seq_parent_num,range(L)] += fpar*child.branch_length
                    T_i[seq_child_num,range(L)] += (1-fpar)*child.branch_length
                    substitutions[seq_child_num,seq_parent_num,range(L)] +=1  
                    
        root_states = GTR_twoseq.seq_to_P(tree.root.seq,alphabet = alphabet)        
        n_ij = np.copy(substitutions)     
        n_ij[range(q),range(q),:] = 0        
        return n_ij, T_i, root_states

def nuc_freq(aln_array, alphabet = 'ACGT'):
    '''calculating frequency of different nucleotides in array'''
#    freq_single=np.zeros((len(alphabet),aln_array.shape[1]))
#    for nuc_ind, nuc in enumerate(alphabet):
#        freq_single[nuc_ind,:] = (aln_array==nuc).mean(axis=0)
    return np.array([(aln_array == nuc).mean(axis=0) for nuc in alphabet])

def GTR_simult(n_ija, T_ia,root_states, dp = 10**(-4), Nit = 10**4):
    '''
    Calculating components of GTR matrix for multiple sites
    
    Input arguments:
    n_ija - transition counts matrix
    T_ia - times spent in each nucleotide state
    root_states - the nucleotide states of the tree root (for pseudocount correction)
    idx_cum - indices of the nucleotides to use
    
    Returns:
    attempt matrix, equilibrium populations by site, mutation rates
    '''
#    q = n_ija.shape[0]; L = n_ija.shape[2]
    (q,L) = T_ia.shape
    n_ija[range(q),range(q),:] = 0
    Lambda = np.sum(root_states,axis=0)
    
    p_ia_old=np.zeros((q,L))    
    mu_a = np.ones(L)
    p_ia = np.ones((q,L))/q
    
    count = 0
    while (LA.norm(p_ia_old-p_ia) > dp or np.abs(1-np.max(p_ia.sum(axis=0))) > dp) and count < Nit:
        count += 1
        p_ia_old = np.copy(p_ia)
        W_ij = (n_ija+np.transpose(n_ija,axes=(1,0,2))).sum(axis = 2)/((mu_a*p_ia).dot(T_ia.T)+T_ia.dot((mu_a*p_ia).T)+h)
        W_ij = W_ij/np.sum(W_ij)
        p_ia = (np.sum(n_ija,axis=1)+root_states)/(mu_a*np.dot(W_ij,T_ia)+Lambda)
        mu_a = n_ija.sum(axis=1).sum(axis=0)/(h+np.diag(p_ia.T.dot(W_ij).dot(T_ia)))
    
    if count >= Nit:
        print 'WARNING: maximum number of iterations has been reached in GTR_simult',\
        np.min(p_ia.sum(axis=0)), np.max(p_ia.sum(axis=0))
        if LA.norm(p_ia_old-p_ia) > dp:
            print '    the iterative scheme has not converged'
        elif np.abs(1-np.max(p_ia.sum(axis=0))) > dp:
            print '    the iterative scheme has converged, but proper normalization was not reached'
            
    return W_ij, p_ia, mu_a
        
if __name__=="__main__":
    '''Generating sequences for a given tree using a GTR matrix'''
    
    plt.ioff()
    outdir_name = '/ebio/ag-neher/share/users/vpuller/GTR_staggered/tmp/'
    if not os.path.exists(outdir_name):
        os.makedirs(outdir_name)    
    
    tree_file_name = '/ebio/ag-neher/share/users/vpuller/Maximum_likelihood/ML_test_sequences/RRE/HIV1B500_RRE_nuc_tree.nwk'    
    bio_tree = Phylo.read(tree_file_name, 'newick')
    bio_tree.root.branch_length = h
    
    fig10=plt.figure(10,figsize = (120,240))
    plt.clf()
    ax10=fig10.add_subplot(1,1,1)
    Phylo.draw(bio_tree,do_show = False,axes=ax10)
    plt.draw()
    plt.savefig(outdir_name+'tree.pdf')
    plt.close(10)
    
    
    # generate root sequence
    alphabet = 'AG' #'ACGT'
    q = len(alphabet)
    L = 10**2 # sequence length
    seq0 = GTR_twoseq.random_seq_p0(L,alphabet = alphabet)
    P0 = seq_to_P(seq0, alphabet = alphabet)
    
    #dressing tree with sequences
    mu = 20.
#    Qij = mu*(np.ones((q,q)) - q*np.eye(q))
#    Wij0 = q*(np.ones((q,q)) - np.eye(q))
#    Qij = np.diag(pi0).dot(Wij0); Qij = Qij - np.diag(np.sum(Qij,axis=0))
    
    mu_a = mu*np.ones(L)
    Wij = q*(np.ones((q,q)) - np.eye(q))
#    p0_a = np.ones((q,L))/q
    p0_a = np.random.exponential(size = (q,L)); p0_a = p0_a/np.sum(p0_a,axis=0)
    Snuc0 = - ((p0_a + h)*np.log(p0_a + h)).sum(axis = 0)
    
#    Qija = np.array([mu_a[ja]*np.diag(p0_a[:,ja]).dot(Wij) - 
#    np.diag(np.sum(mu_a[ja]*np.diag(p0_a[:,ja]).dot(Wij), axis = 0)) for ja in xrange(L)])
#    Qija = np.transpose(Qija,(1,2,0))
#    
#    bio_tree.root.seq = seq0
#    for clade in bio_tree.get_nonterminals():
#        for child in clade:
##            child.seq = GTR_twoseq.evolve_seq_GTR(clade.seq,child.branch_length, Qij = Qij)
#            child.seq = evolve_seq_Qija(clade.seq,child.branch_length, Qija)
    
    bio_tree = dress_tree_seq(bio_tree, seq0, mu_a = mu_a, Wij = Wij, p0_a = p0_a, alphabet = alphabet)
    n_ij, T_i, root_states = nuc_sub_count(bio_tree, alphabet = alphabet)   
    
#    [clade.seq for clade in bio_tree.find_clades()]  
    arr = np.array([list(clade.seq) for clade in bio_tree.find_clades()])
    
    freqs = nuc_freq(arr)
    Snuc = - ((freqs + h)*np.log(freqs + h)).sum(axis = 0)

    fracT = T_i/np.sum(T_i, axis = 0)
    ST =  - ((fracT + h)*np.log(fracT + h)).sum(axis = 0)
    
    
#    GTR = GTR_class.GTR_class(bio_tree, warnings = 'off',Nit_max = 10**4)
#    GTR.substitution_count()            
#    W_ij, p_a, mu_a = GTR.GTR_simultaneous(n_ij,T_i,root_states,range(L)) 
#    S_a1 = -((h+p_a)*np.log(h+p_a)).sum(axis=0)
    
    W_ij, p_a, mu_a = GTR_simult(n_ij,T_i,root_states) 
    S_a = -((h+p_a)*np.log(h+p_a)).sum(axis=0)
    
    p_appx = np.zeros((q,L))
    for jnuc in xrange(q):
        idx = range(q)
        idx.remove(jnuc)
        p_appx[jnuc,:] = n_ij[idx,:,:][:,idx,:].sum(axis = (0,1))/T_i[idx,:].sum(axis = 0)
    p_appx = p_appx/np.sum(p_appx, axis =0)
    S_appx = -((h+p_appx)*np.log(h+p_appx)).sum(axis=0)
    
    plt.figure(10,figsize = (20,12)); plt.clf()
    plt.subplot(2,1,1)
    plt.plot(freqs.T)
    plt.xlabel('site'); plt.ylabel('frequency'); plt.legend(list(alphabet))
    plt.subplot(2,1,2)
    plt.plot(Snuc0/np.log(q))
    plt.plot(Snuc/np.log(q))
    plt.plot(ST/np.log(q))
    plt.plot(S_a/np.log(q))
    plt.plot(S_appx/np.log(q))
    plt.xlabel('site'); plt.ylabel('entropy')
    plt.legend(('model','nucs', 'times', 'GTR','estimate'),loc = 0)
    plt.savefig(outdir_name + 'freqs.pdf')
    plt.close(10)
    
    plt.figure(20, figsize = (6*(q+1),6)); plt.clf()
    for jnuc, nuc in enumerate(alphabet):
        plt.subplot(1,q+1,jnuc+1)
        plt.hist(freqs[jnuc,:])
        plt.xlabel('frequency'); plt.title(nuc)
    plt.subplot(1,q+1,q+1)
    plt.hist(Snuc/np.log(q))
    plt.xlabel('entropy'); plt.title('Snuc')
    plt.savefig(outdir_name + 'freq_hist.pdf')
    plt.close(20)
    
    plt.figure(20, figsize = (6*(q+1),6)); plt.clf()
    for jnuc, nuc in enumerate(alphabet):
        plt.subplot(1,q+1,jnuc+1)
        plt.hist(fracT[jnuc,:])
        plt.xlabel('fraction of time'); plt.title(nuc)
    plt.subplot(1,q+1,q+1)
    plt.hist(ST/np.log(q))
    plt.xlabel('entropy'); plt.title('Snuc')
    plt.savefig(outdir_name + 'fracT_hist.pdf')
    plt.close(20)
    
    plt.figure(20, figsize = (6*(q+1),6)); plt.clf()
    for jnuc, nuc in enumerate(alphabet):
        plt.subplot(1,q+1,jnuc+1)
        plt.hist(p_a[jnuc,:])
        plt.xlabel('fraction of time'); plt.title(nuc)
    plt.subplot(1,q+1,q+1)
    plt.hist(S_a/np.log(q))
    plt.xlabel('entropy'); plt.title('Snuc')
    plt.savefig(outdir_name + 'freqGTR_hist.pdf')
    plt.close(20)
    
    