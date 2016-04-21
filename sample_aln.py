# -*- coding: utf-8 -*-
"""
Created on Tue Apr  5 12:02:49 2016

@author: vpuller
"""
from __future__ import division
import numpy as np
from scipy import misc
from Bio import Phylo, AlignIO
import matplotlib.pyplot as plt
import sys, os
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
import time
#import subprocess as sp

import GTR_twoseq
import generate_tree_GTR as gentree

#sys.path.append('/ebio/ag-neher/share/users/vpuller/myTOOLS/') 
#import Vadim_toolbox_file as vp
#import GTR_class

sys.path.append('/ebio/ag-neher/share/users/vpuller')
#import time_tree
#from time_tree import io
from time_tree.treetime.gtr import GTR
from time_tree.treetime import io

#Constants
h= 10**(-8)

def aln_from_tree(tree, all_nodes = True):   
    aln = MultipleSeqAlignment([], annotations={})
    if all_nodes:
        for jclade,clade in enumerate(tree.find_clades()):
            seq_record=SeqRecord(Seq(clade.seq, generic_dna), id=str(clade.name))
            aln.append(seq_record)
    else:
        for jclade,clade in enumerate(tree.get_terminals()):
            seq_record=SeqRecord(Seq(clade.seq, generic_dna), id=str(clade.name))
            aln.append(seq_record)
    return aln
        
if __name__=="__main__":
    '''
    given a tree topology:
    - generate root sequence
    - generate GTR model
    - evolve root sequence along the tree
    - sample from leaves sequences
    - generate a tree from the sample
    - infer ancestral sequences
    '''
    plt.ioff()
    outdir_name = '/ebio/ag-neher/share/users/vpuller/GTR_staggered/report_sampling/'
    if not os.path.exists(outdir_name):
        os.makedirs(outdir_name)    
    
    file_name_head = outdir_name + ''
    
    tree_file_name = '/ebio/ag-neher/share/users/vpuller/Maximum_likelihood/ML_test_sequences/RRE/HIV1B2000_RRE_nuc_tree.nwk'    
    bio_tree = Phylo.read(tree_file_name, 'newick')
    bio_tree.root.branch_length = h
  
    alphabet = 'ACGT'
    q = len(alphabet)
    L = 10**2 # sequence length
    
    
    #specifying GTR model
    # defining mutation rates
    mu0 = .2
    mu_a0 = mu0*np.ones(L)
#    mu_a0 = mu0*np.random.exponential(size = L)
    
    # defining attempt frequencies
#    Wij = q*(np.ones((q,q)) - np.eye(q))
    beta = q; alpha = 3.*beta
    Wij = beta*(np.ones((q,q)) - np.eye(q)) +\
    (alpha - beta)*(np.diag(np.ones(2),2) + np.diag(np.ones(2),k = -2))
    
    # defining equilibrium populations
#    p0_a = np.ones((q,L))/q
    p0_a = np.random.exponential(size = (q,L)); p0_a = p0_a/np.sum(p0_a,axis=0)
    Snuc0 = - ((p0_a + h)*np.log(p0_a + h)).sum(axis = 0)

    # generate root sequence
    seq0 = GTR_twoseq.random_seq_p0(L,alphabet = alphabet)
#    seq0 = 'A'*L
#    seq0 = vp.numbers_to_string(np.argmax(p0_a, axis = 0))
    P0 = gentree.seq_to_P(seq0, alphabet = alphabet)
    
    #dressing tree with sequences
    t0 = time.time()
    dress = gentree.dress_tree(bio_tree, seq0, mu_a = mu_a0, Wij = Wij, p0_a = p0_a, alphabet = alphabet)
    t1 = time.time()
    print 't = ', t1-t0
    
    np.savetxt(outdir_name + 'Wij.txt', Wij)
    np.savetxt(outdir_name + 'mu_a.txt',mu_a0)
    np.savetxt(outdir_name + 'p0_a.txt',p0_a)
    aln = aln_from_tree(dress.tree, all_nodes = False)
    with open(outdir_name + 'aln_terminals.fasta','w') as aln_file:
        AlignIO.write(aln,aln_file,'fasta')    
    
    n_ij, T_i, root_states = gentree.nuc_sub_count(dress.tree, alphabet = alphabet)  
    W_ij, p_a, mu_a = gentree.GTR_simult(n_ij,T_i,root_states) 
    dist_GTR0 = np.sum((p_a - p0_a)**2)
    
    arr = np.array(aln)
    freqs = gentree.nuc_freq(arr, alphabet = alphabet)
    dist_freqs0 = np.sum((freqs - p0_a)**2)
    
    #sampling from the alignment
#    os.system(" ".join(call))
    
    Nsamples = 10
    NNseq = range(100,len(aln),100) # number of sequences in a sample
    dist_GTR = np.zeros((Nsamples,len(NNseq)))
    dist_freqs = np.zeros((Nsamples,len(NNseq)))
    Ttree = np.zeros((Nsamples,len(NNseq)))
    
    for jN, Nseq in enumerate(NNseq):
        jjsamples = np.random.randint(0,high = len(aln),size = (Nsamples,Nseq))
        for jsample, jj in enumerate(jjsamples):
            print jN, jsample
            aln_sample = MultipleSeqAlignment([], annotations={})
            for j in np.unique(jj):
                aln_sample.append(aln[j])
            aln_file_name = outdir_name + 'aln_sample_tmp.fasta'
            with open(aln_file_name,'w') as aln_file:
                AlignIO.write(aln_sample,aln_file,'fasta')  
            
            tree_file_name = outdir_name + 'tree_sample_tmp.fasta'
            call = ['/ebio/ag-neher/share/programs/bin/fasttree', '-nt','-quiet',\
            aln_file_name, ' > ', tree_file_name]
            os.system(' '.join(call))
            
            gtr = GTR.standard()
            t = io.treetime_from_newick(gtr, tree_file_name)
            io.set_seqs_to_leaves(t, AlignIO.read(aln_file_name, 'fasta'))
    #        t.optimize_seq_and_branch_len()
            t.reconstruct_anc('ml')
            Ttree[jsample,jN] = t.tree.total_branch_length()
            
            for clade in t.tree.find_clades():
                clade.seq = ''.join(clade.sequence)
                
            n_ij, T_i, root_states = gentree.nuc_sub_count(t.tree, alphabet = alphabet)   
            W_ij, p_a, mu_a = gentree.GTR_simult(n_ij,T_i,root_states) 
            dist_GTR[jsample,jN] = np.sum((p_a - p0_a)**2)
            
    #        arr = np.array([list(clade.seq) for clade in dress.tree.get_terminals()])
            arr = np.array(aln_sample)
            freqs = gentree.nuc_freq(arr, alphabet = alphabet)
    #        Snuc = - ((freqs + h)*np.log(freqs + h)).sum(axis = 0)
            dist_freqs[jsample,jN] = np.sum((freqs - p0_a)**2)
            
            os.remove(aln_file_name); os.remove(tree_file_name)
    
    dist_freqs_mean = dist_freqs.mean(axis=0)
    dist_freqs_var = (dist_freqs**2).mean(axis=0) - dist_freqs_mean**2 
    dist_GTR_mean = dist_GTR.mean(axis=0)
    dist_GTR_var = (dist_GTR**2).mean(axis=0) - dist_GTR_mean**2    
    Ttree_mean = Ttree.mean(axis=0)
    Ttree_var = (Ttree**2).mean(axis=0) - Ttree_mean**2
    
    plt.figure(30,figsize = (30,10)); plt.clf()
    plt.subplot(1,3,1)
    plt.errorbar(NNseq, dist_freqs_mean, yerr = np.sqrt(dist_freqs_var))
    plt.errorbar(NNseq, dist_GTR_mean, yerr = np.sqrt(dist_GTR_var))
    plt.plot(NNseq, dist_GTR0*np.ones(len(NNseq)),'--')
    plt.plot(NNseq, dist_freqs0*np.ones(len(NNseq)),'--')
    plt.xlabel(r'$N_{seq}$', fontsize = 18)
    plt.ylabel(r'$\chi^2$', fontsize = 18)
    plt.legend(('tree seq.', 'GTR'), fontsize = 18)

    plt.subplot(1,3,2)
    plt.errorbar(NNseq, np.log10(dist_freqs_mean), yerr = np.sqrt(dist_freqs_var)/dist_freqs_mean/np.log(10))
    plt.errorbar(NNseq, np.log10(dist_GTR_mean), yerr = np.sqrt(dist_GTR_var)/dist_GTR_mean/np.log(10))
    plt.plot(NNseq, np.log10(dist_GTR0)*np.ones(len(NNseq)),'--')
    plt.plot(NNseq, np.log10(dist_freqs0)*np.ones(len(NNseq)),'--')
    plt.xlabel(r'$N_{seq}$', fontsize = 18);
    plt.ylabel(r'$\log_{10}\chi^2$', fontsize = 18)
    plt.legend(('tree seq.', 'GTR'), fontsize = 18)

    plt.subplot(1,3,3)
    plt.errorbar(np.log10(NNseq), np.log10(dist_freqs_mean), yerr = np.sqrt(dist_freqs_var)/dist_freqs_mean/np.log(10))
    plt.errorbar(np.log10(NNseq), np.log10(dist_GTR_mean), yerr = np.sqrt(dist_GTR_var)/dist_GTR_mean/np.log(10))
    plt.plot(np.log10(NNseq), np.log10(dist_GTR0)*np.ones(len(NNseq)),'--')
    plt.plot(np.log10(NNseq), np.log10(dist_freqs0)*np.ones(len(NNseq)),'--')
    plt.xlabel(r'$\log_{10} N_{seq}$', fontsize = 18)
    plt.ylabel(r'$\log_{10}\chi^2$', fontsize = 18)
    plt.legend(('tree seq.', 'GTR'), loc = 0, fontsize = 18)
    plt.savefig(file_name_head + 'dist_vs_N.pdf')
    plt.close(30)

    
    plt.figure(30,figsize = (30,10)); plt.clf()
    plt.subplot(1,3,1)
    plt.errorbar(Ttree_mean*mu0, dist_freqs_mean, yerr = np.sqrt(dist_freqs_var))
    plt.errorbar(Ttree_mean*mu0, dist_GTR_mean, yerr = np.sqrt(dist_GTR_var))
    plt.plot(Ttree_mean*mu0, dist_GTR0*np.ones(len(Ttree_mean)),'--')
    plt.plot(Ttree_mean*mu0, dist_freqs0*np.ones(len(Ttree_mean)),'--')
    plt.xlabel(r'$\mu T_{tree}$', fontsize = 18)
    plt.ylabel(r'$\chi^2$', fontsize = 18)
    plt.legend(('tree seq.', 'GTR'), fontsize = 18)

    plt.subplot(1,3,2)
    plt.errorbar(Ttree_mean*mu0, np.log10(dist_freqs_mean), yerr = np.sqrt(dist_freqs_var)/dist_freqs_mean/np.log(10))
    plt.errorbar(Ttree_mean*mu0, np.log10(dist_GTR_mean), yerr = np.sqrt(dist_GTR_var)/dist_GTR_mean/np.log(10))
    plt.plot(Ttree_mean*mu0, np.log10(dist_GTR0)*np.ones(len(Ttree_mean)),'--')
    plt.plot(Ttree_mean*mu0, np.log10(dist_freqs0)*np.ones(len(Ttree_mean)),'--')
    plt.xlabel(r'$\mu T_{tree}$', fontsize = 18)
    plt.ylabel(r'$\log_{10}\chi^2$', fontsize = 18)
    plt.legend(('tree seq.', 'GTR'), fontsize = 18)

    plt.subplot(1,3,3)
    plt.errorbar(np.log10(Ttree_mean*mu0), np.log10(dist_freqs_mean), yerr = np.sqrt(dist_freqs_var)/dist_freqs_mean/np.log(10))
    plt.errorbar(np.log10(Ttree_mean*mu0), np.log10(dist_GTR_mean), yerr = np.sqrt(dist_GTR_var)/dist_GTR_mean/np.log(10))
    plt.plot(np.log10(Ttree_mean*mu0), np.log10(dist_GTR0)*np.ones(len(Ttree_mean)),'--')
    plt.plot(np.log10(Ttree_mean*mu0), np.log10(dist_freqs0)*np.ones(len(Ttree_mean)),'--')
    plt.xlabel(r'$\log_{10} (\mu T_{tree})$', fontsize = 18)
    plt.ylabel(r'$\log_{10}\chi^2$', fontsize = 18)
    plt.legend(('tree seq.', 'GTR'), loc = 0, fontsize = 18)
    plt.savefig(file_name_head + 'dist_vs_T.pdf')
    plt.close(30)    
        