# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 16:53:18 2017

@author: vpuller
"""
from __future__ import division
import numpy as np

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

def create_GTR_model(L, model_spec, alphabet = 'ACGT', model_file_head = None):
    '''
    Create site-specific GTR model for a sequence
    
    Input arguments:
    L - sequence length
    alphabet - sequence alphabet
    mu_model - model for mutation rates: uniform, exp or gamma
    gamma_k - shape of gamma distribution, if gamma model is used
    pi_model - model for equilibrated frequencies: equal or simplex
    W_model - model for attempt matrix, if None - use Jukes-Cantor
    seq0_model - model for the root sequence: equal or equilibrated
    model_file_head - location for saving the model parameters
    '''
    q = len(alphabet)

    # defining mutation rates
#    mu0 = 2.
    if model_spec['mu'][0] == 'uniform':
        mu0 = model_spec['mu'][1]
        mu0_a = mu0*np.ones(L)
    elif model_spec['mu'][0] == 'exp':
        mu0 = model_spec['mu'][1]
#        mu0_a = mu0*np.random.exponential(size = L)
        mu0_a = np.random.exponential(scale = mu0, size = L)
    elif model_spec['mu'][0] == 'gamma':
        mu0, gamma_k = model_spec['mu'][1:3]
        mu0_a = mu0*np.random.gamma(gamma_k, mu0/gamma_k, size = L)
    else:
        print 'WARNING: unknown model for mu'
        
    # defining attempt matrix
    if model_spec['Wij'][0] == 'JC':
        Wij0 = q*(np.ones((q,q)) - np.eye(q))
#        beta = 1; alpha = beta
#        Wij0 = beta*(np.ones((q,q)) - np.eye(q)) +\
#        (alpha - beta)*(np.diag(np.ones(2),2) + np.diag(np.ones(2),k = -2))
    else:
        print 'WARNING: unknown model for Wij'
    Wij0 = Wij0/Wij0.sum()
    
    # defining equilibrium populations
    if model_spec['pi'] == 'equal':
        p0_a = np.ones((q,L))/q
    elif model_spec['pi'] == 'simplex':
        p0_a = np.random.exponential(size = (q,L)); p0_a = p0_a/np.sum(p0_a,axis=0)
    else:
        print 'WARNING: unknown model for p_a'

    # generate root sequence
    if model_spec['seq0'] == 'equal':
        seq0 = random_seq_p0(L,alphabet = alphabet)
    elif model_spec['seq0'] == 'equilibrated':
        seq0 = random_seq_p0(L,p0 = p0_a, alphabet = alphabet)
    else:
        print 'WARNING: unknown model for seq0'
        
    if model_file_head is not None:
        # Saving GTR model
        np.savetxt(model_file_head + 'mu0_a.txt', mu0_a)
        np.savetxt(model_file_head + 'p0_a.txt', p0_a)
        np.savetxt(model_file_head + 'Wij0.txt', Wij0)
        with open(model_file_head + 'seq0.txt','w') as seq0_file:
            seq0_file.write(seq0)
        
    return Wij0, p0_a, mu0_a, seq0
