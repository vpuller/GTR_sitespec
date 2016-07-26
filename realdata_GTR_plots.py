# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 12:32:36 2016

@author: vpuller
"""
from __future__ import division
import numpy as np
from scipy import linalg as LA
import matplotlib.pyplot as plt
import pickle as pickle
from scipy.stats import spearmanr    
    
#from Bio import AlignIO
import os

#Constants
h= 10**(-8)
alphabet = 'ACGT'

def nuc_freq(aln_array, alphabet = 'ACGT'):
    '''calculating frequency of different nucleotides in array'''
    return np.array([(aln_array == nuc).mean(axis=0) for nuc in alphabet])
    
#def site_entropy(freqs, correct_bias = True):
#    if correct_bias:
#        (q, L) = freqs.shape
#        return -np.sum((freqs + h)*np.log(freqs + h), axis = 0) + (q-1)/(2*L)
#    else:
#        return -np.sum((freqs + h)*np.log(freqs + h), axis = 0)

def site_entropy(freqs):
    return -np.sum((freqs + h)*np.log(freqs + h), axis = 0)
        
def convert(string):
    if string[0] in ['>', '<']:
        return float(string[1:])
    elif string == 'nan':
        return 0.001
    else:
        return float(string)

def bin_scatter(X, Y, q):
    nn, edges = np.histogram(X, q)
    jjq = [np.where((edges[jq] < X)*(X <= edges[jq+1] ))[0] for jq in xrange(q)]
    return np.array([X[jj].mean() for jj in jjq]), np.array([Y[jj].mean() for jj in jjq])

def percentile_scatter(X, Y, q):
#    nn, edges = np.histogram(X, q)
    edges = np.array([np.percentile(X, 100*jq/q) for jq in xrange(q+1)])
    edges[0] = 0.
    jjq = [np.where((edges[jq] < X)*(X <= edges[jq+1] ))[0] for jq in xrange(q)]
    return np.array([X[jj].mean() for jj in jjq]), np.array([Y[jj].mean() for jj in jjq])    

def eigvals_GTR(GTRmodel):
        (Wij, p_a, mu_a) = GTRmodel
        w_a = np.zeros(p_a.shape)
        for ja in xrange(p_a.shape[1]):
            Qij = mu_a[ja]*np.diag(p_a[:,ja]).dot(Wij)
            Qij = Qij - np.diag(Qij.sum(axis=0))
            (w,vr) = LA.eig(Qij)
            w_a[:,ja] = np.sort(w.real)
        return w_a
        
if __name__=="__main__":
    '''GTR reconstruction for HIV data: plots'''
    
    plt.ioff()
    plt.close('all')
    subtype = 'ALL' #'ALL' # 'B'
    gene = 'vpu' #'vpr' #'vpu' #'nef' #'env' #'gag' # 'vif' #'pol'
    dir_name = '/ebio/ag-neher/share/users/vpuller/GTR_staggered/HIV1{}_test/'.format(subtype)
    outdir_name = dir_name + gene + '/'
    if not os.path.exists(outdir_name):
        os.makedirs(outdir_name)  

    

    # alignment statistics    
    freqs = np.loadtxt(dir_name + '{}_aln_p_a.txt'.format(gene))
    idx_HXB2 = np.loadtxt(dir_name + '{}_idx_HXB2.txt'.format(gene), dtype = 'int')
    S_a = site_entropy(freqs)
    
    with open(dir_name + '{}_GTRanc_pavel/'.format(gene) + 'GTRinf.pickle','r') as infile:
        (W2_ij, p2_a, mu2_a) = pickle.load(infile)
    branch_lengths2 = np.loadtxt(dir_name + '{}_GTRanc_pavel/'.format(gene) + 'branch_lengths.txt')
    S2_a = site_entropy(p2_a)

    # loading fitness coefficients data
#    datafile_name = '/ebio/ag-neher/share/users/vpuller/HIV_fitness_landscape/data/nuc_{}_selection_coeffcients_B.tsv'.format(gene)
    datafile_name = '/ebio/ag-neher/share/users/rneher/HIV_fitness_landscape/data/fitness_pooled/nuc_{}_fitness_costs_B.tsv'.format(gene)
    cons_arr, HXB2_arr = np.loadtxt(datafile_name, skiprows = 2, usecols = (1,2), dtype = 'str', unpack = True)
    cons_seq = ''.join(cons_arr)
    HXB2_seq = ''.join(HXB2_arr)
    
#    lower, median, upper = np.loadtxt(datafile_name, skiprows = 2, usecols = (3,4,5), dtype = 'str', unpack = True)
    lower, median, upper = np.loadtxt(datafile_name, skiprows = 2, usecols = (3,4,5), \
    converters = {3: convert, 4: convert, 5: convert},unpack = True)
    
#    idx_syn = range(0, len(idx_HXB2),3)
#    idx_HXB2 = idx_HXB2[idx_syn]
#    median = median[idx_syn]

#    idx_finite = np.where((.001 < median)*(median < .1))[0]    
#    idx_HXB2 = idx_HXB2[idx_finite]
#    median = median[idx_finite]

#    idx1 = np.where(np.max(p2_a[:,idx_HXB2], axis = 0) > 0.6)[0]
#    idx_HXB2 = idx_HXB2[idx1]
#    median = median[idx1]


    plt.figure(10, figsize = (24,6)); plt.clf()
    plt.subplot(1,4,1)
    plt.hist(S_a/np.log(len(alphabet)))
    plt.xlabel(r'$S_{\alpha}$')
    plt.title('aln')

    plt.subplot(1,4,2)
    plt.hist(S2_a/np.log(len(alphabet)))
    plt.xlabel(r'$S_{\alpha}$')
    plt.title(r'$GTR, \pi_{i\alpha}$')

    plt.subplot(1,4,3)
    plt.hist(mu2_a*branch_lengths2.mean())
    plt.xlabel(r'$\mu_{\alpha}\bar{t}_{branch}$' )
    plt.title(r'$GTR, \mu$')
    
    plt.subplot(1,4,4)
    plt.hist(median)
    plt.xlabel(r'$s$' )
    plt.title(r'$fitness$')
    plt.savefig(outdir_name + '{}_S_mu.pdf'.format(gene))
    plt.close(10)
    
    # correlating fitness with entropy     
    xx = [S_a[idx_HXB2]/np.log(len(alphabet)), S2_a[idx_HXB2]/np.log(len(alphabet)), mu2_a[idx_HXB2]*branch_lengths2.mean()]
    titls = ['aln', r'GTR, $\pi_{i\alpha}$', r'GTR, $\mu$']
    labls = [r'$S_{\alpha}$', r'$S_{\alpha}$', r'$\mu_{\alpha}\bar{t}_{branch}$' ]
    plt.figure(40, figsize = (18,24)); plt.clf()
    for jx, x in enumerate(xx):
        plt.subplot(4,3,jx + 1)
        plt.scatter(x, median)
        plt.xlabel(labls[jx], fontsize = 18)
        plt.ylabel('s', fontsize = 18)
        plt.title(titls[jx], fontsize = 18)
        print '\n', spearmanr(x, median)
        print np.corrcoef(x, median)
    
        plt.subplot(4,3,3 + jx + 1)
        plt.scatter(x, np.log10(median))
        plt.xlabel(labls[jx], fontsize = 18)
        plt.ylabel(r'$\log_{10}s$', fontsize = 18)
        print np.corrcoef(x, np.log10(median))
        
        plt.subplot(4,3,6 + jx + 1)
        plt.scatter(np.log10(x+h), median)
        plt.xlabel(r'$\log_{10}$' + labls[jx], fontsize = 18)
        plt.ylabel('s', fontsize = 18)
        print np.corrcoef(np.log10(x+h), median)
#        plt.xlim(-6., 0.)
        
        plt.subplot(4,3,9 + jx + 1)
        plt.scatter(np.log10(x+h), np.log10(median), c=np.log(xx[(jx-1)%3]))
        plt.xlabel(r'$\log_{10}$' + labls[jx], fontsize = 18)
        plt.ylabel(r'$\log_{10}s$', fontsize = 18)
        print np.corrcoef(np.log10(x+h), np.log10(median))
#        plt.xlim(-3., 0.)
    plt.savefig(outdir_name + '{}_corr_new.pdf'.format(gene))
    plt.close(40) 

    
    # inferring fitness coefficient from reconstructed populations data
#    seff = np.diff(np.sort(np.log(p2_a[:, idx_HXB2] + h), axis = 0), axis =0)[-1,:]
#    pmax_a = np.max(p2_a[:, idx_HXB2], axis = 0)
    pmax_a = np.max(freqs[:, idx_HXB2], axis = 0)
    seff = np.log(pmax_a) - np.log(1. - pmax_a + h)
#    F_a = np.log(p2_a[:, idx_HXB2]+h)
#    seff = np.sort(F_a - np.max(F_a, axis = 0), axis = 0)[]
    plt.figure(50, figsize = (24,6)); plt.clf()
    plt.subplot(1,4,1)
    plt.scatter(seff, median)
    plt.xlabel('s_{eff}'); plt.ylabel('s')
    
    plt.subplot(1,4,2)
    plt.scatter(seff, np.log10(median))
    
    plt.subplot(1,4,3)
    plt.scatter(np.log10(seff), median)
    
    plt.subplot(1,4,4)
    plt.scatter(np.log10(seff), np.log10(median))
    plt.savefig(outdir_name + '{}_seff.pdf'.format(gene))
    plt.close(50)
    print '\n', spearmanr(seff, median)
    print np.corrcoef(seff, median)[0,1]
    print np.corrcoef(seff, np.log10(median))[0,1]
    print np.corrcoef(np.log10(seff), median)[0,1]
    print np.corrcoef(np.log10(seff), np.log10(median))[0,1]
    
    
    plt.figure(60, figsize = (18,6)); plt.clf()
    plt.subplot(1,3,1)
    plt.scatter(seff, mu2_a[idx_HXB2])
    plt.xlabel(r'$s_{eff}$')
    plt.ylabel(r'$\mu_{\alpha}$')
    
    plt.subplot(1,3,2)
    plt.scatter(seff, np.log10(mu2_a[idx_HXB2]))
    plt.xlabel(r'$s_{eff}$')
    plt.ylabel(r'$\log_{10}\mu_{\alpha}$')
    
    plt.subplot(1,3,3)
    plt.scatter(np.log10(seff), np.log10(mu2_a[idx_HXB2]))
    plt.xlabel(r'$\log_{10}s_{eff}$')
    plt.ylabel(r'$\log_{10}\mu_{\alpha}$')
    plt.savefig(outdir_name + '{}_mu_v_seff.pdf'.format(gene))
    plt.close(60)
    print '\n', spearmanr(seff, mu2_a[idx_HXB2])
    
#>>> plt.figure(); plt.scatter(xx[-1]+.0001, xx[-2]+.0001, c=np.log(median))
#<matplotlib.figure.Figure object at 0x7f499c19cfd0>
#<matplotlib.collections.PathCollection object at 0x7f49a2f80c10>
#>>> plt.ylim(0.00001,2)
#(1e-05, 2)
#>>> plt.xlim(0.00001,2)
#(1e-05, 2)
#>>> 
    