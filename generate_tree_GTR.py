# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 18:09:38 2016

@author: vpuller
"""

from __future__ import division
import numpy as np
from scipy import linalg as LA
from Bio import Phylo, AlignIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
import matplotlib.pyplot as plt
#from scipy import optimize
import sys, os
import GTR_twoseq
import copy as copy
import time

sys.path.append('/ebio/ag-neher/share/users/vpuller/myTOOLS/') 
import Vadim_toolbox_file as vp

sys.path.append('/ebio/ag-neher/share/users/vpuller')
from time_tree.treetime.gtr import GTR
from time_tree.treetime import io

#Constants
h= 10**(-8)

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

def eigLR(M):
    (w,vr) = LA.eig(M)
    vl = np.zeros(M.shape, dtype = 'complex')
    for j in xrange(M.shape[0]):
        b = np.zeros(M.shape[0]); b[j] = 1.
        vl[:,j] = LA.solve(vr.T,b)
    
    return w.real, vl.real, vr.real
    
class dress_tree(object):
    def __init__(self,tree,seq_root,mu_a = None, Wij = None, p0_a = None,\
    alphabet = 'ACGT', use_eigenvalues = True, logL = False):
        '''
        TODO
        '''
        self.q = len(alphabet); self.L = len(seq_root)
        self.tree = copy.deepcopy(tree)
        self.alphabet = alphabet
        #defining mutation matrix
        if mu_a is None:
            mu_a = np.ones(L)
        
        if Wij is None:
            Wij = self.q*(np.ones((self.q,self.q)) - np.eye(self.q))
            
        if p0_a is None:
           p0_a = np.ones((self.q,L))

        self.tree.root.seq = seq_root
        self.tree.root.P = self.seq_to_P(seq_root)
          
        if use_eigenvalues:
            self.Qija = np.zeros((self.L,self.q,self.q))
            self.w_a = np.zeros((self.L,self.q))
            self.vl_a = np.zeros((self.L,self.q,self.q))
            self.vr_a = np.zeros((self.L,self.q,self.q))
            for ja in xrange(self.L):
                Qij = mu_a[ja]*np.diag(p0_a[:,ja]).dot(Wij)
                Qij = Qij - np.diag(Qij.sum(axis=0))
#                (w,vl,vr) = LA.eig(Qij, left = True)
                (w,vl,vr) = eigLR(Qij)
                self.w_a[ja,:] = w
                self.vl_a[ja,:,:] = vl
                self.vr_a[ja,:,:] = vr/np.diag(vl.T.dot(vr))
                vr = vr/np.diag(vl.T.dot(vr))
                self.Qija[ja,:,:] = Qij
#                self.Qija[ja,:,:] = vr.dot(np.diag(w)).dot(vl.T)
#                print '\n',ja, np.sum((Qij - vr.dot(np.diag(w)).dot(vl.T))**2), (vl.T.dot(vr)).sum()/q
            if logL:
                self.logL = 0.
                for parent in self.tree.get_nonterminals():
                    P_parent = self.seq_to_P(parent.seq)
                    for child in parent:
                        P_child = self.seq_to_P(child.seq)
                        self.logL += np.log(np.array([P_child[:,ja].dot(self.vr_a[ja,:,:])\
                        .dot(np.diag(np.exp(self.w_a[ja,:]*child.branch_length)))\
                        .dot(self.vl_a[ja,:,:].T).dot(P_parent[:,ja]) 
                        for ja in xrange(self.L)])).sum()
#                        print parent.name, child.name, logL
            else:
                for parent in self.tree.get_nonterminals():
                    for child in parent:
                        child.P = self.evolve_seq_eig(parent.P,child.branch_length)
                        child.seq = self.P_to_seq(child.P)
        else:
            Qija = np.array([mu_a[ja]*np.diag(p0_a[:,ja]).dot(Wij) - 
            np.diag(np.sum(mu_a[ja]*np.diag(p0_a[:,ja]).dot(Wij), axis = 0)) for ja in xrange(L)])
            self.Qija = np.transpose(Qija,(1,2,0))
        
            for parent in self.tree.get_nonterminals():
                for child in parent:
                    child.P = self.evolve_seq_Qija_P(parent.P,child.branch_length)
                    child.seq = self.P_to_seq(child.P)
                
    def seq_to_P(self,seq):
#        q = len(alphabet); L = len(seq)
        seq_num=np.zeros(len(seq),dtype = 'int')
        seq_arr = np.array(list(seq))
        for jbase, base in enumerate(self.alphabet):
            seq_num[np.where(seq_arr == base)] = jbase   
        P = np.zeros((self.q,self.L),dtype = 'int')
        P[seq_num,range(self.L)] = 1.
        return P
    
    def P_to_seq(self,P):
        seq_num = np.zeros(self.L,dtype = 'int')
        seq_num[np.where(P > 0)[1]] = np.where(P >0)[0]
        return ''.join(self.numbers_to_string(seq_num))
    
    def numbers_to_string(self,seq_num):
        '''converting sequence of nucleotides represented by numbers 
        into a char array of nucleotides'''
        seq=np.zeros(seq_num.shape[0],dtype='|S1')
        for jbase, base in enumerate(self.alphabet):
            seq[np.where(seq_num == jbase)] = base
        return ''.join(seq)
    
#    def string_to_numbers(self,seq):
#        '''
#        Function for converting string to numbers according to the symbols's place in the alphabet
#        '''
#        seq_array = np.array(list(seq))
#        seq_num=np.zeros(len(seq_array),dtype = 'int')
#        for jbase, base in enumerate(self.alphabet):
#            seq_num[np.where(seq_array == base)] = jbase    
#        return seq_num
    
    def evolve_seq_Qija_P(self,P0, t):
        '''
        Evolve sequence using site - specific GTR model
        '''
        # evolving sequence
        Pt = np.array([LA.expm(self.Qija[:,:,ja]*t).dot(P0[:,ja]) for ja in xrange(self.L)]).T
        
        # measuring sequence
        MC = np.random.uniform(size = self.L)
        P1 = np.zeros((self.q,self.L))
        for jq in xrange(self.q):
            P1[jq,np.where((Pt[:jq,:].sum(axis=0) < MC)*(MC < Pt[:jq+1,:].sum(axis=0)))[0]] = 1
        return P1    

    def evolve_seq_eig(self,P0, t):
        '''
        Evolve sequence using site - specific GTR model
        '''
        # evolving sequence
        Pt = np.array([self.vr_a[ja,:,:].dot(np.diag(np.exp(self.w_a[ja,:]*t))).dot(self.vl_a[ja,:,:].T).dot(P0[:,ja]) 
        for ja in xrange(self.L)]).T

        # measuring sequence
        MC = np.random.uniform(size = self.L)
        P1 = np.zeros((self.q,self.L))
        for jq in xrange(self.q):
            P1[jq,np.where((Pt[:jq,:].sum(axis=0) < MC)*(MC < Pt[:jq+1,:].sum(axis=0)))[0]] = 1
        return P1              
    
def nuc_sub_count(tree, alphabet = 'ACGT', fpar = 0.5):        
        '''
        Function that counts the numbers of nucleotide substitutions on tree
        
        Input arguments:
        tree - tree dressed withs equences
        alphabet - nucleotide alphabet
        fpar - fraction of the branch length assigned to the parent nucletide state
        '''
#        print alphabet
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
                    
        root_states = seq_to_P(tree.root.seq,alphabet = alphabet)        
        n_ij = np.copy(substitutions)     
        n_ij[range(q),range(q),:] = 0   
        return n_ij, T_i, root_states

def nuc_freq(aln_array, alphabet = 'ACGT'):
    '''calculating frequency of different nucleotides in array'''
    return np.array([(aln_array == nuc).mean(axis=0) for nuc in alphabet])

def GTR_simult(sub_ija, T_ia,root_states, dp = 10**(-4), Nit = 10**4):
    '''
    Calculating components of GTR matrix for multiple sites
    
    Input arguments:
    sub_ija - transition counts matrix
    T_ia - times spent in each nucleotide state
    root_states - the nucleotide states of the tree root (for pseudocount correction)
    idx_cum - indices of the nucleotides to use
    
    Returns:
    attempt matrix, equilibrium populations by site, mutation rates
    '''
#    q = n_ija.shape[0]; L = n_ija.shape[2]
    (q,L) = T_ia.shape
    n_ija = np.copy(sub_ija)
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

def mean_var(dist):
    dist_mean = dist.mean(axis=0)
    dist_var = (dist**2).mean(axis=0) - dist_mean**2 
    return dist_mean, dist_var

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
    '''Generating sequences for a given tree using a GTR matrix'''
    
    plt.ioff()
    outdir_name = '/ebio/ag-neher/share/users/vpuller/GTR_staggered/report2/'
    if not os.path.exists(outdir_name):
        os.makedirs(outdir_name)    
    
    file_name_head = outdir_name + ''
    
    tree_file_name = '/ebio/ag-neher/share/users/vpuller/Maximum_likelihood/ML_test_sequences/RRE/HIV1B500_RRE_nuc_tree.nwk'    
    bio_tree = Phylo.read(tree_file_name, 'newick')
    bio_tree.root.branch_length = h
  
    alphabet = 'ACGT'
    q = len(alphabet)
    L = 10**2 # sequence length
    
    
    #specifying GTR model
    # defining mutation rates
    mu0 = 2.
    mu_a0 = mu0*np.ones(L)
#    mu_a0 = mu0*np.random.exponential(size = L)
    
    # defining attempt frequencies
#    Wij = q*(np.ones((q,q)) - np.eye(q))
    beta = q; alpha = beta
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
    P0 = seq_to_P(seq0, alphabet = alphabet)
    
    #dressing tree with sequences
    t0 = time.time()
    dress = dress_tree(bio_tree, seq0, mu_a = mu_a0, Wij = Wij, p0_a = p0_a, alphabet = alphabet)
    t1 = time.time()
    print 't = ', t1-t0
    
    # array of all sequences on the tree
#    arr = np.array([list(clade.seq) for clade in dress.tree.find_clades()])
    arr = np.array([list(clade.seq) for clade in dress.tree.get_terminals()])
    freqs = nuc_freq(arr, alphabet = alphabet)
    Snuc = - ((freqs + h)*np.log(freqs + h)).sum(axis = 0)
    
    # GTR fitting
    n_ij, T_i, root_states = nuc_sub_count(dress.tree, alphabet = alphabet)   
    W_ij, p_a, mu_a = GTR_simult(n_ij,T_i,root_states) 
    S_a = -((h+p_a)*np.log(h+p_a)).sum(axis=0)
    
    dp_a = p_a/np.sqrt(n_ij.sum(axis=1))
    dp_a[np.where(np.sqrt(n_ij.sum(axis=1)) == 0)] = 1.
    dS_a = -(dp_a*np.log(h+p_a)).sum(axis=0)
    
    
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
#    plt.plot(ST/np.log(q))
    plt.plot(S_a/np.log(q))
#    plt.errorbar(np.arange(L),S_a/np.log(q), yerr = dS_a/np.log(q))
#    plt.plot(S_appx/np.log(q))
    plt.xlabel('site'); plt.ylabel('entropy')
    plt.legend(('model','tree seq.','GTR'),loc = 0)
    plt.savefig(file_name_head + 'freqs.pdf')
    plt.close(10)
    
    plt.figure(20, figsize = (6*(q+1),6)); plt.clf()
    for jnuc, nuc in enumerate(alphabet):
        plt.subplot(1,q+1,jnuc+1)
        plt.hist(freqs[jnuc,:])
        plt.xlabel('frequency'); plt.title(nuc)
    plt.subplot(1,q+1,q+1)
    plt.hist(Snuc/np.log(q))
    plt.xlabel('entropy'); plt.title('Snuc')
    plt.savefig(file_name_head + 'freq_hist.pdf')
    plt.close(20)
    
    plt.figure(20, figsize = (6*(q+1),6)); plt.clf()
    for jnuc, nuc in enumerate(alphabet):
        plt.subplot(1,q+1,jnuc+1)
        plt.hist(p_a[jnuc,:])
        plt.xlabel('fraction of time'); plt.title(nuc)
    plt.subplot(1,q+1,q+1)
    plt.hist(S_a/np.log(q))
    plt.xlabel('entropy'); plt.title('Snuc')
    plt.savefig(file_name_head + 'freqGTR_hist.pdf')
    plt.close(20)

    plt.figure(20, figsize = (18,6)); plt.clf()
    plt.subplot(1,3,1)
    for jnuc in xrange(q):
        nn_bin, x_bin = np.histogram(p0_a[jnuc,:],bins = 10, range = (0,1))
        plt.plot(x_bin[:-1] + .5*x_bin[1], nn_bin/L)
        plt.legend(alphabet); plt.title('Model')
    plt.subplot(1,3,2)
    for jnuc in xrange(q):
        nn_bin, x_bin = np.histogram(freqs[jnuc,:],bins = 10, range = (0,1))
        plt.plot(x_bin[:-1] + .5*x_bin[1], nn_bin/L)
        plt.legend(alphabet); plt.title('Alignment')
    plt.subplot(1,3,3)
    for jnuc in xrange(q):
        nn_bin, x_bin = np.histogram(p_a[jnuc,:],bins = 10, range = (0,1))
        plt.plot(x_bin[:-1] + .5*x_bin[1], nn_bin/L)
        plt.legend(alphabet); plt.title('GTR inference')
    plt.savefig(file_name_head + 'hist.pdf')
    
    #testing dependence of the quality of GTR reconstruction on mutation rate
    #dressing tree with sequences
    branch_lengths = np.array([clade.branch_length for clade in dress.tree.find_clades()])
    branch_mean = branch_lengths.mean()
    branch_var = (branch_lengths**2).mean() - branch_mean**2
    dist_to_root = np.array([dress.tree.distance(clade) for clade in dress.tree.get_terminals()])
    dist_to_root_mean = dist_to_root.mean()
    
    mu_max = 1/branch_lengths.max()
#    mmu = np.linspace(0.1,10,num = 20)
    mmu = 0.1*10**np.linspace(0,2,num=10)

##    Wij = q*(np.ones((q,q)) - np.eye(q))
#    beta = q; alpha = 3.*beta
#    Wij = beta*(np.ones((q,q)) - np.eye(q)) +\
#    (alpha - beta)*(np.diag(np.ones(2),2) + np.diag(np.ones(2),k = -2))
#    
##    p0_a = np.ones((q,L))/q
#    p0_a = np.random.exponential(size = (q,L)); p0_a = p0_a/np.sum(p0_a,axis=0)
##    Snuc0 = - ((p0_a + h)*np.log(p0_a + h)).sum(axis = 0)
    
    Nsample = 5
    dist_GTR = np.zeros((Nsample,mmu.shape[0]))
    dist_freqs = np.zeros((Nsample,mmu.shape[0]))
    dist_GTR_anc = np.zeros((Nsample,mmu.shape[0]))
    dist_GTR_tree = np.zeros((Nsample,mmu.shape[0]))
#    Hamdist_GTR = np.zeros((Nsample,mmu.shape[0]))
#    Hamdist_freqs = np.zeros((Nsample,mmu.shape[0]))
    for jsample in xrange(Nsample):
        print 'sample #{}'.format(jsample)
        for jmu, mu in enumerate(mmu):
            mu_a = mu*np.ones(L)
#            mu_a = mu*mu_a0/mu0
            
            dress = dress_tree(bio_tree, seq0, mu_a = mu_a, Wij = Wij, p0_a = p0_a, alphabet = alphabet, use_eigenvalues = True)
            
            # array of all sequences on the tree
            arr = np.array([list(clade.seq) for clade in dress.tree.find_clades()])
            freqs = nuc_freq(arr, alphabet = alphabet)
    #        Snuc = - ((freqs + h)*np.log(freqs + h)).sum(axis = 0)
            dist_freqs[jsample,jmu] = np.sqrt(np.sum((freqs - p0_a)**2))
#            Hamdist_freqs[jsample,jmu] = L - np.sum(freqs*p0_a)
            
            
            # GTR fitting
            n_ij, T_i, root_states = nuc_sub_count(dress.tree, alphabet = alphabet)   
            W_ij, p_a, mu_a = GTR_simult(n_ij,T_i,root_states) 
    #        S_a = -((h+p_a)*np.log(h+p_a)).sum(axis=0)
            dist_GTR[jsample,jmu] = np.sqrt(np.sum((p_a - p0_a)**2))
#            Hamdist_GTR[jsample,jmu] = L - np.sum(p_a*p0_a)
            
            
            # GTR fitting with ancestral reconstruction
            aln = aln_from_tree(dress.tree, all_nodes = False)
            aln_file_name = outdir_name + 'aln_tmp.fasta'
            with open(aln_file_name,'w') as aln_file:
                AlignIO.write(aln,aln_file,'fasta')    
            
            gtr = GTR.standard()
            t = io.treetime_from_newick(gtr, tree_file_name)
            io.set_seqs_to_leaves(t, AlignIO.read(aln_file_name, 'fasta'))
    #        t.optimize_seq_and_branch_len()
            t.reconstruct_anc('ml')
#            Ttree[jsample,jN] = t.tree.total_branch_length()
            
            for clade in t.tree.find_clades():
                clade.seq = ''.join(clade.sequence)
                
            n_ij, T_i, root_states = nuc_sub_count(t.tree, alphabet = alphabet)   
            W_ij, p_a, mu_a = GTR_simult(n_ij,T_i,root_states) 
            
            dist_GTR_anc[jsample,jmu] = np.sum((p_a - p0_a)**2)
            
            
            # GTR fitting with tree reconstruction
            rec_tree_file_name = outdir_name + 'tree_tmp.fasta'
            call = ['/ebio/ag-neher/share/programs/bin/fasttree', '-nt','-quiet',\
            aln_file_name, ' > ', rec_tree_file_name]
            os.system(' '.join(call))
            
            gtr = GTR.standard()
            t = io.treetime_from_newick(gtr, rec_tree_file_name)
            io.set_seqs_to_leaves(t, AlignIO.read(aln_file_name, 'fasta'))
    #        t.optimize_seq_and_branch_len()
            t.reconstruct_anc('ml')
#            Ttree[jsample,jN] = t.tree.total_branch_length()
            
            for clade in t.tree.find_clades():
                clade.seq = ''.join(clade.sequence)
                
            n_ij, T_i, root_states = nuc_sub_count(t.tree, alphabet = alphabet)   
            W_ij, p_a, mu_a = GTR_simult(n_ij,T_i,root_states) 
            dist_GTR_tree[jsample,jmu] = np.sum((p_a - p0_a)**2)
            
            os.remove(aln_file_name); #os.remove(rec_tree_file_name)   
    
    data_all = [dist_freqs, dist_GTR, dist_GTR_anc, dist_GTR_tree]
    means = np.zeros((len(data_all),len(mmu)))
    variances = np.zeros((len(data_all),len(mmu)))
    for jdata, data in enumerate(data_all):
        means[jdata,:], variances[jdata,:] = mean_var(data)
    
#    dist_freqs_mean, dist_freqs_var = mean_var(dist_freqs)
#    dist_GTR_mean, dist_GTR_var = mean_var(dist_GTR)
#    dist_GTR_anc_mean, dist_GTR_anc_var = mean_var(dist_GTR_anc)
#    dist_GTR_tree_mean, dist_GTR_tree_var = mean_var(dist_GTR_tree)
    
    leg = ('Tips alignment', 'GTR', 'GTR ancesors', 'GTR tree and ancestors')
#    means = [dist_freqs_mean, dist_GTR_mean, dist_GTR_anc_mean, dist_GTR_tree_mean]
#    variances = [dist_freqs_var, dist_GTR_var, dist_GTR_anc_var, dist_GTR_tree_var]
    
    plt.figure(30,figsize = (20,10)); plt.clf()
    plt.subplot(1,2,1)
    for jdata, mean_loc in enumerate(means):
        plt.errorbar(mmu*dist_to_root_mean, np.log10(means[jdata,:]),\
        yerr = np.sqrt(variances[jdata,:])/means[jdata,:]/np.log(10))
#    plt.errorbar(mmu*dist_to_root_mean, np.log10(dist_freqs_mean), yerr = np.sqrt(dist_freqs_var)/dist_freqs_mean/np.log(10))
#    plt.errorbar(mmu*dist_to_root_mean, np.log10(dist_GTR_mean), yerr = np.sqrt(dist_GTR_var)/dist_GTR_mean/np.log(10))
#    plt.errorbar(mmu*dist_to_root_mean, np.log10(dist_GTR_anc_mean), yerr = np.sqrt(dist_GTR_anc_var)/dist_GTR_anc_mean/np.log(10))
    plt.xlabel(r'$\mu\bar{t}_{root}$', fontsize = 18)
    plt.ylabel(r'$\log_{10}\chi^2$', fontsize = 18)
    plt.legend(leg, fontsize = 18)

    plt.subplot(1,2,2)
    for jdata, mean_loc in enumerate(means):
        plt.errorbar(np.log10(mmu*dist_to_root_mean), np.log10(means[jdata,:]),\
        yerr = np.sqrt(variances[jdata,:])/means[jdata,:]/np.log(10))
#    plt.errorbar(np.log10(mmu*dist_to_root_mean), np.log10(dist_freqs_mean), yerr = np.sqrt(dist_freqs_var)/dist_freqs_mean/np.log(10))
#    plt.errorbar(np.log10(mmu*dist_to_root_mean), np.log10(dist_GTR_mean), yerr = np.sqrt(dist_GTR_var)/dist_GTR_mean/np.log(10))
#    plt.errorbar(np.log10(mmu*dist_to_root_mean), np.log10(dist_GTR_anc_mean), yerr = np.sqrt(dist_GTR_anc_var)/dist_GTR_anc_mean/np.log(10))
    plt.xlabel(r'$\log_{10} (\mu\bar{t}_{root})$', fontsize = 18)
    plt.ylabel(r'$\log_{10}\chi^2$', fontsize = 18)
    plt.legend(leg, loc = 0, fontsize = 18)
    plt.savefig(file_name_head + 'dist.pdf')
    plt.close(30)
    
    