# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 16:32:39 2016

@author: vpuller
"""

from __future__ import division
import numpy as np
from scipy import linalg as LA
#import scipy.stats as stats
import matplotlib.pyplot as plt
#import pandas as pd
from scipy import optimize
import sys, os
#import ANOVA_functions as ANOVA

sys.path.append('/ebio/ag-neher/share/users/vpuller/myTOOLS/') 
import Vadim_toolbox_file as vp

#Constants
h= 10**(-8)

def random_seq(L,alphabet = 'ACGT'):
    '''Generate a randoms equence of length L'''
    q = len(alphabet)
    seq0_num = np.random.randint(0,q,size = L)
    return vp.numbers_to_string(seq0_num,alphabet = alphabet)
    
def evolve_seq(seq, t, mu = 1.,alphabet = 'ACGT'):
    '''
    Evolve sequence according to the symmetric GTR model
    '''
    #defining mutation matrix
    q = len(alphabet)
    qt = (1. - np.exp(-mu*q*t))/q
    pt = (1. + (q-1.)*np.exp(-mu*q*t))/q
    pij = np.ones((q,q))*qt + np.eye(q)*(pt - qt)
    
    # converting sequence to vector form
#    seq_arr = np.array(list(seq))
    seq_num = vp.string_to_numbers(seq,alphabet = alphabet)
    P0 = np.zeros((q,L))
    P0[seq_num,range(L)] = 1.
    
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
#    seq1_arr = vp.numbers_to_nucs(seq1_num)
#    seq1 = ''.join(seq1_arr)
    return ''.join(vp.numbers_to_nucs(seq1_num))


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
    Evolve sequence according to the symmetric GTR model
    
    Input parameters:
    seq - sequence
    t - time
    Qij - GTR matrix
    '''
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
#    seq_arr = np.array(list(seq))
    seq_num = vp.string_to_numbers(seq,alphabet = alphabet)
    P0 = np.zeros((q,L))
    P0[seq_num,range(L)] = 1.
    
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
#    seq1_arr = vp.numbers_to_nucs(seq1_num)
#    seq1 = ''.join(seq1_arr)
    return ''.join(vp.numbers_to_nucs(seq1_num))

#def muWpi_to_Qij(muWpi,q):
#    mu = muWpi[0]
#    ww = muWpi[1:-(q-1)]
#    pi0 = muWpi[-(q-1):]; pi0 = np.append(pi0,1. - pi0.sum())
#    Wij = np.zeros((q,q))
#    l = 0
#    for jq in xrange(1,q):
#        Wij += np.diag(ww[l:l+q-jq],jq)
#        l += q-jq
#    Wij = Wij + Wij.T
#    Qij = mu*np.diag(pi0).dot(Wij)
#    Qij = Qij - np.diag(Qij.sum(axis=0))
#    return mu, Wij, pi0, Qij

def muWpi_to_Qij(muWpi,q):
    mu = muWpi[0]
    ww = muWpi[1:-q]
    pi0 = muWpi[-q:]
    Wij = np.zeros((q,q))
    l = 0
    for jq in xrange(1,q):
        Wij += np.diag(ww[l:l+q-jq],jq)
        l += q-jq
    Wij = Wij + Wij.T
    Qij = mu*np.diag(pi0).dot(Wij)
    Qij = Qij - np.diag(Qij.sum(axis=0))
    return mu, Wij, pi0, Qij
    
def Wij_to_ww(Wij):
    q = Wij.shape[0]
    l = 0
    ww = np.zeros(int(q*(q-1)/2))
    for jq in xrange(1,q):
        ww[l:l + q -jq] = np.diag(Wij,jq)
        l += q - jq
    return ww
    
#def inferGTR(seq1,seq2,t,alphabet = 'ACGT'):
#    '''
#    Inferring GTR matrix from two sequences
#    '''
#    def logLike_GTR(muWpi):
#        #Reconstruct GTR matrix from muWpi 
##        mu, Wij, pi0, Qij = muWpi_to_Qij(muWpi**2,q)
#        mu = muWpi**2
#        Qij = mu*(np.ones((q,q)) - q*np.eye(q))
#        Pij = LA.expm(Qij*t)
#        logL = - np.log((P2*Pij.dot(P1)).sum(axis=0)).sum()
#        if np.abs(pi0.sum() - 1.) > 10**(-2):
#            logL = (pi0.sum() - 1.)**2
#        elif np.abs(Wij.sum() -1.) > 10**(-2):
#            logL = (Wij.sum() -1.)**2
#        else:
#            logL = - np.log((P2*Pij.dot(P1)).sum(axis=0)).sum()
#            
##        print mu, Wij.sum(), pi0.sum(), logL,\
##        (P2*P1).sum(), P1.sum()
##        np.log((P2*np.eye(q).dot(P1)).sum(axis=0)).sum()
#        if np.isnan(logL):
#            print '\nmu = ',mu,'\nWij =\n', Wij,Wij.sum(),\
#            '\npi0 = ', pi0, pi0.sum(), '\nQij = ', Qij,\
#            '\nL = ', logL,\
#            Pij.dot(P1)
#            raw_input('continue?')
#        return logL
#        
#    q = len(alphabet)
#    P1 = seq_to_P(seq1,alphabet = alphabet)
#    P2 = seq_to_P(seq2,alphabet = alphabet)
#
#    mu = 1.
#    Wij = np.triu(np.ones((q,q)),1); Wij = Wij + Wij.T; Wij = Wij/Wij.sum()
#    ww = Wij_to_ww(Wij)
#    pi0 = np.ones(q)/q
##    muWpi0 = [mu]; muWpi0.extend(ww); muWpi0.extend(pi0)
##    muWpi0 = np.sqrt(np.array(muWpi0))
##    step = .1*muWpi0
##    tol = 10**(-8)
##    muWpi = vp.amoeba_vp(logLike_GTR,muWpi0,args=(),a = step,tol_x = tol,tol_f = tol, Nit = 10**4)
##    return muWpi_to_Qij(muWpi**2,q)[:-1]
#    res = optimize.minimize_scalar(logLike_GTR)
#    return res.x**2, res.x,res.x

def inferGTR_simple(seq1,seq2,t,alphabet = 'ACGT'):
    '''
    Inferring GTR matrix from two sequences
    '''
    def logLike_GTR_simple(mu):
#        qt = (1. - np.exp(-mu*q*t))/q
#        pt = (1. + (q-1.)*np.exp(-mu*q*t))/q
#        pij = np.ones((q,q))*qt + np.eye(q)*(pt - qt)
        
        Qij = mu*(np.ones((q,q)) - q*np.eye(q))
        pij = LA.expm(Qij*t)
        logL = - np.log((P2*pij.dot(P1)).sum(axis=0)).sum()
        
#        seq1_arr = np.array(list(seq1))
#        seq2_arr = np.array(list(seq2))
#        n = np.count_nonzero(seq1_arr == seq2_arr)
#        m = np.count_nonzero(seq1_arr != seq2_arr)
#        logL0 = - n*np.log(pt) - m*np.log(qt)
#        print mu, logL, logL0
#        if np.isnan(logL):
#            print P1.T,P2.T,pij,mu,logL
#            raw_input('continue?')
        return logL
        
    q = len(alphabet)
    P1 = seq_to_P(seq1,alphabet = alphabet)
    P2 = seq_to_P(seq2,alphabet = alphabet)

    res = optimize.minimize_scalar(logLike_GTR_simple)
    return res.x

#def logLike_GTR_tmp(seq1,seq2,t,mu = 1.,alphabet = 'ACGT'):
#    P1 = seq_to_P(seq1,alphabet = alphabet)
#    P2 = seq_to_P(seq2,alphabet = alphabet)
#    qt = (1. - np.exp(-mu*q*t))/q
#    pt = (1. + (q-1.)*np.exp(-mu*q*t))/q
#    pij = np.ones((q,q))*qt + np.eye(q)*(pt - qt)
#    logL = np.log((P2*pij.dot(P1)).sum(axis=0)).sum()
#    return logL
        
def seq_to_P(seq,alphabet = 'ACGT'):
    seq_num=np.zeros(len(seq),dtype = 'int')
    seq_arr = np.array(list(seq))
    for jbase, base in enumerate(alphabet):
        seq_num[np.where(seq_arr == base)] = jbase   
#    seq_num = vp.string_to_numbers(seq)
    P = np.zeros((q,L))
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
        Wij = (mij + mij.T)/(p0[:,np.newaxis]*nii + nii[:,np.newaxis]*p0)
        p = mij.sum(axis=1)/Wij.dot(nii)
#        print 's = {}, p = {}, Wij = \n{}'.format(s,p, Wij)
        if np.isnan(p).any():
            print 's = {}, p = {},\nnij = \n{}'.format(s,p, nij)
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
    seq0 = random_seq(L)
    seq0_num = vp.string_to_numbers(seq0)
    seq0_arr = np.array(list(seq0))
    P0 = np.zeros((q,L))
    P0[seq0_num,range(L)] = 1.
    
#    t = .01
#    seq1 = evolve_seq(seq0,t)
##    seq0 = 'GGAGCAATCATCCACTGAGT'
##    seq1 = 'GGAGCAATCATCCACTGCGT'
#    nij = mutations_count(seq0,seq1)
#    print 'n_ij = \n{}\n'.format(nij)
#    
#    Wij, p0 = inferGTR_iter(nij,t = 1.)
#    print '\nmu = ', np.sum(Wij), '\nWij = \n', Wij, '\np0 = ',p0  
#    print 't = {}'.format(np.mean(np.diag(p0).dot(Wij)[np.triu_indices(q,1)]))
#    print inferGTR_simple(seq0,seq1,1.), dist_simple(seq0,seq1)
    
#    t = .2 #branch length
    Nseq = 10**2
    tt = np.arange(.01,.1,.01)
    
    Qij = np.ones((q,q)) - q*np.eye(q)
#    beta0 = 1.; alpha0 = 1.
#    beta = 3*beta0/(alpha0 + 2*beta0); alpha = 3*alpha0/(alpha0 + 2*beta0)
#    Qij = beta*np.ones((q,q)) + (alpha - beta)*(np.diag(np.ones(2),2) + np.diag(np.ones(2),k = -2))
#    Qij = Qij - np.diag(np.sum(Qij,axis=0))
    
    tt_MC = np.zeros((tt.shape[0],Nseq))
    tt_iter = np.zeros((tt.shape[0],Nseq))
    n_overfit = np.zeros(tt.shape[0], dtype = 'int')
    for jt, t in enumerate(tt):
        for jseq in xrange(Nseq):            
            seq1 = evolve_seq_GTR(seq0, t, Qij = Qij)
            tt_MC[jt,jseq] = dist_simple(seq0,seq1)
            
            nij = mutations_count(seq0,seq1)
            if np.count_nonzero(nij == 0) > 0:
                n_overfit[jt] += 1
                continue
            else:    
                Wij, p0 = inferGTR_iter(nij,t = 1.)
                tt_iter[jt,jseq] = np.mean(np.diag(p0).dot(Wij)[np.triu_indices(q,1)])

#    tt_MC = np.array([[dist_simple(seq0,evolve_seq_GTR(seq0, t, Qij = Qij)) for j in xrange(Nseq)] for t in tt])
    tt_mean = np.array([tt_MC[jt,np.where(np.isfinite(tt_MC[jt,:]))[0]].mean() for jt, t in enumerate(tt)]) 
    tt_var = np.array([(tt_MC[jt,np.where(np.isfinite(tt_MC[jt,:]))[0]]**2).mean() for jt, t in enumerate(tt)]) - tt_mean**2
    
    tt_mean1 = np.array([tt_iter[jt,np.where(tt_iter[jt,:] > 0)[0]].mean() for jt, t in enumerate(tt)]) 
    tt_var1 = np.array([(tt_iter[jt,np.where(tt_iter[jt,:] > 0)[0]]**2).mean() for jt, t in enumerate(tt)]) - tt_mean1**2

    
    plt.figure(20); plt.clf()
#    plt.errorbar(tt,tt_mean0, yerr = np.sqrt(tt_var0), fmt = ':o')
    plt.errorbar(tt,tt_mean, yerr = np.sqrt(tt_var), fmt = ':d')
    plt.errorbar(tt,tt_mean1, yerr = np.sqrt(tt_var1), fmt = ':p')
    plt.plot(tt,tt,'--')
    plt.xlabel('branch length')
    plt.ylabel('estimated branch length')
    plt.xlim([tt[0] - .01, tt[-1] + .01])
    plt.legend(('Simple','iter','True'),loc = 0)
    plt.savefig(outdir_name + 'time_divergence.pdf'); plt.close(20)
    