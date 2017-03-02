#!/usr/bin/env python
"""
Created on Wed May 11 11:25:13 2016

@author: vpuller
"""
from __future__ import division
import numpy as np
from scipy import linalg as LA
from Bio import Phylo #, AlignIO

#Constants
h= 10**(-8)

class dress_tree(object):
    '''Class that evolves sequences along the tree, givern tree topology'''
    def __init__(self,tree_file_name,seq_root, mu_a = None, Wij = None, p0_a = None,\
    alphabet = 'ACGT', use_eigenvalues = True, logL = False):
        '''
        TODO
        '''
        self.q = len(alphabet); self.L = len(seq_root)
#        self.tree = copy.deepcopy(tree)
        self.tree = Phylo.read(tree_file_name, 'newick')
        self.tree.root.branch_length = h  
        self.alphabet = alphabet
        #defining mutation matrix
        if mu_a is None:
            mu_a = np.ones(self.L)
        
        if Wij is None:
            Wij = self.q*(np.ones((self.q,self.q)) - np.eye(self.q))
            
        if p0_a is None:
           p0_a = np.ones((self.q,self.L))

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
                (w,vl,vr) = self.eigLR(Qij)
                self.w_a[ja,:] = w
                self.vl_a[ja,:,:] = vl
                self.vr_a[ja,:,:] = vr/np.diag(vl.T.dot(vr))
                vr = vr/np.diag(vl.T.dot(vr))
                self.Qija[ja,:,:] = Qij
#                if ja==8:
#                    import pdb; pdb.set_trace()
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
#                        jcheck = 32 #[ 7, 32, 74, 84]
#                        if child.seq[jcheck] != parent.seq[jcheck]:
#                            print parent.seq[jcheck] + '-->' + child.seq[jcheck] 
        else:
            Qija = np.array([mu_a[ja]*np.diag(p0_a[:,ja]).dot(Wij) - 
            np.diag(np.sum(mu_a[ja]*np.diag(p0_a[:,ja]).dot(Wij), axis = 0)) for ja in xrange(self.L)])
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

    def eigLR(self,M):
        (w,vr) = LA.eig(M)
        ij = np.where(np.iscomplex(vr).sum(axis=0))[0]
        if len(ij) ==2:
            vr[:,ij] = np.array([vr[:,ij[0]].real , vr[:,ij[0]].imag]).T
    
        vl = np.zeros(M.shape, dtype = 'complex')
        for j in xrange(M.shape[0]):
            b = np.zeros(M.shape[0]); b[j] = 1.
            vl[:,j] = LA.solve(vr.T,b)
        
        return w.real, vl.real, vr.real