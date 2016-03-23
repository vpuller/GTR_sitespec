# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 16:32:39 2016

@author: vpuller
"""

from __future__ import division
import numpy as np
from scipy import linalg as LA
import matplotlib.pyplot as plt
from scipy import optimize
import sys, os

sys.path.append('/ebio/ag-neher/share/users/vpuller/myTOOLS/') 
import Vadim_toolbox_file as vp

#Constants
h= 10**(-8)

def random_seq(L,alphabet = 'ACGT'):
    '''Generate a random sequence of length L'''
    q = len(alphabet)
    seq0_num = np.random.randint(0,q,size = L)
    return vp.numbers_to_string(seq0_num,alphabet = alphabet)

def random_seq_p0(L,p0 = None, alphabet = 'ACGT'):
    '''Generate a random sequence of length L 
    with probabilities of nucleotides given by vector p0'''
    q = len(alphabet)
    if p0 is None:
        p0 = np.ones(q)/q
    else:
        p0 = p0/p0.sum()
    MC = np.random.uniform(size = L)
    seq_arr = np.zeros(L,dtype = 'S1')
    for jnuc, nuc in enumerate(alphabet):
        seq_arr[np.where((p0[:jnuc].sum() < MC)*(MC < p0[:jnuc+1].sum()))[0]] = nuc
    return ''.join(seq_arr)
    
def dist_simple(seq1,seq2,q = 4):
    '''
    Estimating distance between two sequences relying on a simple GTR model
    '''
    seq1_arr = np.array(list(seq1))
    seq2_arr = np.array(list(seq2))
    n = np.count_nonzero(seq1_arr == seq2_arr)
    m = np.count_nonzero(seq1_arr != seq2_arr)
    x = (n*(q-1) - m)/((n + m)*(q-1))
#    print n, m, x, q, (n*(q-1) - m),((n + m)*(q-1))
    return - np.log(x)/q

def evolve_seq_GTR(seq, t, Qij = None, mu = 1., alphabet = 'ACGT'):
    '''
    Evolve sequence according to GTR model
    
    Input parameters:
    seq - sequence
    t - time
    Qij - GTR matrix
    '''
    L = len(seq)
    #defining mutation matrix
    if Qij is None:
        q = len(alphabet)
        qt = (1. - np.exp(-mu*q*t))/q
        pt = (1. + (q-1.)*np.exp(-mu*q*t))/q
        pij = np.ones((q,q))*qt + np.eye(q)*(pt - qt)
    else:
        q = Qij.shape[0]
        pij = LA.expm(Qij*t)
    
    # converting sequence to vector form
    P0 = seq_to_P(seq,alphabet = alphabet)
    
    # evolving sequence
    Pt = pij.dot(P0)
    
    # measuring sequence
    MC = np.random.uniform(size = L)
    P1 = np.zeros((q,L))
    for jq in xrange(q):
        P1[jq,np.where((Pt[:jq,:].sum(axis=0) < MC)*(MC < Pt[:jq+1,:].sum(axis=0)))[0]] = 1.
    
    # converting result to sequence form
    seq1_num = np.zeros(L,dtype = 'int')
    seq1_num[np.where(P1 > 0)[1]] = np.where(P1 >0)[0]
    return ''.join(vp.numbers_to_nucs(seq1_num))

def inferGTR_simple(seq1,seq2,t,alphabet = 'ACGT'):
    '''
    Inferring GTR matrix from two sequences
    '''
    def logLike_GTR_simple(mu):
        Qij = mu*(np.ones((q,q)) - q*np.eye(q))
        pij = LA.expm(Qij*t)
        logL = - np.log((P2*pij.dot(P1)).sum(axis=0)).sum()
        return logL
        
    q = len(alphabet)
    P1 = seq_to_P(seq1,alphabet = alphabet)
    P2 = seq_to_P(seq2,alphabet = alphabet)

    res = optimize.minimize_scalar(logLike_GTR_simple)
    return res.x
        
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

#Iterative scheme implementation
def mutations_count(seq1,seq2,alphabet = 'ACGT'):
    seq1_arr = np.array(list(seq1))
    seq2_arr = np.array(list(seq2))
    nij = np.zeros((len(alphabet),len(alphabet)),dtype = 'int')
    for jnuc1, nuc1 in enumerate(alphabet):
        for jnuc2, nuc2 in enumerate(alphabet):
            nij[jnuc2,jnuc1] = np.count_nonzero((seq1_arr == nuc1)*(seq2_arr == nuc2))
    return nij

def inferGTR_iter(nij, t = 1., Nit = 10**4, tol = 10**(-8)):
    mij = nij - np.diag(np.diag(nij))
    nii = np.diag(nij)
    q = nij.shape[0]
    p0 = np.ones(q)/q
    
    for s in xrange(Nit):
        Wij = (mij + mij.T)/(p0[:,np.newaxis]*nii + nii[:,np.newaxis]*p0 + h)
        p = mij.sum(axis=1)/Wij.dot(nii)
#        print '\ns = {}, p = {},\nWij = \n{}'.format(s,p, Wij)
#        print Wij.dot(nii)
        if np.isnan(p).any():
            print 's = {}, p = {},\nnij = \n{},\nWij = \n{},\nmij = \n{}'.format(s,p, nij, Wij, mij)
#            print (mij + mij.T)
#            print (p0[:,np.newaxis]*nii + nii[:,np.newaxis]*p0)
        if LA.norm(p - p0)/LA.norm(p) < tol:
            break
        else:
            p0 = p/p.sum()
    return Wij/t, p

if __name__=="__main__":
    '''testing GTR inference for distance between two sequences'''
    
    plt.ioff()
    outdir_name = '/ebio/ag-neher/share/users/vpuller/GTR_staggered/tmp/'
    if not os.path.exists(outdir_name):
        os.makedirs(outdir_name)
    
    q = 4 #alphabet size
    L = 10**3
    pi0 = np.array([0.7,.1,.1,.1])
#    pi0 = np.array([.3,.3,.3,.1])
    seq0 = random_seq_p0(L,pi0)
    seq0_arr = np.array(list(seq0))
    P0 = seq_to_P(seq0)
    

    Nseq = 10**2
    tt = np.arange(.02,.1,.01)
    Qij = np.ones((q,q)) - q*np.eye(q)
    Wij0 = q*(np.ones((q,q)) - np.eye(q))
    
#    pi0 = np.array([1. - 3*h,h,h,h])
#    pi0 = np.ones(q)/q
    Qij = np.diag(pi0).dot(Wij0); Qij = Qij - np.diag(np.sum(Qij,axis=0))
    
#    beta0 = 1.; alpha0 = 1.
#    beta = 3*beta0/(alpha0 + 2*beta0); alpha = 3*alpha0/(alpha0 + 2*beta0)
#    Qij = beta*np.ones((q,q)) + (alpha - beta)*(np.diag(np.ones(2),2) + np.diag(np.ones(2),k = -2))
#    Qij = Qij - np.diag(np.sum(Qij,axis=0))
    
    Qtype = 'JC' #'random_pi' #'JC'
    tt_MC = np.zeros((tt.shape[0],Nseq))
    tt_iter = np.zeros((tt.shape[0],Nseq))
    n_overfit = np.zeros(tt.shape[0], dtype = 'int')
    mut = np.zeros((tt.shape[0],Nseq))
    for jt, t in enumerate(tt):
        for jseq in xrange(Nseq):
            if Qtype == 'random_pi':
                pi0 = np.random.uniform(size = q); pi0 = pi0/pi0.sum()
#                print pi0
#                pi0 = np.ones(q)/q
#                pi0 = np.array([.7,.1,.1,.1])
                Qij = np.diag(pi0).dot(Wij0); Qij = Qij - np.diag(np.sum(Qij,axis=0))
            
            seq1 = evolve_seq_GTR(seq0, t, Qij = Qij)
            tt_MC[jt,jseq] = dist_simple(seq0,seq1)
            
            mut[jt,jseq] = np.count_nonzero(np.array(list(seq1)) != seq0_arr)
            nij = mutations_count(seq0,seq1)
            if np.count_nonzero(nij + nij.T == 0) > 0:
                n_overfit[jt] += 1
                continue
            else:    
                Wij, p0 = inferGTR_iter(nij,t = 1.)
#                print p0
#                tt_iter[jt,jseq] = np.mean(np.diag(p0).dot(Wij)[np.triu_indices(q,1)])
                tt_iter[jt,jseq] = np.mean(Wij[np.triu_indices(q,1)])/4

    tt_mean = np.array([tt_MC[jt,np.where(np.isfinite(tt_MC[jt,:]))[0]].mean() for jt, t in enumerate(tt)]) 
    tt_var = np.array([(tt_MC[jt,np.where(np.isfinite(tt_MC[jt,:]))[0]]**2).mean() for jt, t in enumerate(tt)]) - tt_mean**2
    
    tt_mean1 = np.array([tt_iter[jt,np.where(tt_iter[jt,:] > 0)[0]].mean() for jt, t in enumerate(tt)]) 
    tt_var1 = np.array([(tt_iter[jt,np.where(tt_iter[jt,:] > 0)[0]]**2).mean() for jt, t in enumerate(tt)]) - tt_mean1**2

    plt.figure(20); plt.clf()
    plt.errorbar(tt,tt_mean, yerr = np.sqrt(tt_var), fmt = ':d')
    plt.errorbar(tt,tt_mean1, yerr = np.sqrt(tt_var1), fmt = ':p')
    plt.plot(tt,tt,'--')
    plt.xlabel('branch length')
    plt.ylabel('estimated branch length')
    plt.xlim([tt[0] - .01, tt[-1] + .01])
    plt.legend(('Simple','iter','True'),loc = 0)
    plt.savefig(outdir_name + 'time_divergence_pi.pdf'); plt.close(20)
    
    
    mut_mean = mut.mean(axis=1)
    mut_var = (mut**2).mean(axis=1) - mut_mean**2
    plt.figure(20); plt.clf()
    plt.errorbar(tt,mut_mean, yerr = np.sqrt(mut_var), fmt = ':d')
    plt.xlabel('branch length')
    plt.ylabel('# mutations')
    plt.xlim([tt[0] - .01, tt[-1] + .01])
#    plt.legend(('Simple','iter','True'),loc = 0)
    plt.savefig(outdir_name + 'mutations.pdf'); plt.close(20)
    