# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 16:32:39 2016

@author: vpuller
"""

from __future__ import division
import numpy as np
#import scipy.stats as stats
import matplotlib.pyplot as plt
import pandas as pd
#from scipy import optimize
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
    
    # convertings equence to vector form
#    seq_arr = np.array(list(seq))
    seq_num = vp.string_to_numbers(seq,alphabet = alphabet)
    P0 = np.zeros((q,L))
    P0[seq_num,range(L)] = 1.
    
    # evolving sequence
    Pt = pij.dot(P0)
    
    # measurings equence
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

if __name__=="__main__":
    '''testing GTR inference for distance between two sequences'''
    
    plt.ioff()
    outdir_name = '/ebio/ag-neher/share/users/vpuller/GTR_staggered/tmp/'
    if not os.path.exists(outdir_name):
        os.makedirs(outdir_name)
    
    q = 4 #alphabet size
    L = 10**3
#    seq0_num = np.random.randint(0,q,size = L)
#    seq0 = ''.join(vp.numbers_to_nucs(seq0_num))
    seq0 = random_seq(L)
    seq0_num = vp.string_to_numbers(seq0)
    seq0_arr = np.array(list(seq0))
    P0 = np.zeros((q,L))
    P0[seq0_num,range(L)] = 1.
    
    
    t = .2 #branch length
    Nseq = 10**3
    tt = np.array([dist_simple(seq0,evolve_seq(seq0, t)) for j in xrange(Nseq)])
    
    plt.figure(10); plt.clf()
    plt.hist(tt[np.where(np.isfinite(tt))[0]],bins = 20)
    plt.xlabel('branch length')
    plt.savefig(outdir_name + 'hist.pdf'); plt.close(10)
        
    