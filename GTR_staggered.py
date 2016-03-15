# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 12:06:56 2015

@author: vpuller
"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import sys
from Bio import AlignIO, Phylo
import os
#import itertools

sys.path.append('/ebio/ag-neher/share/users/vpuller/myTOOLS/') 
#sys.path.append('./') 
import ancestor_reconstruction_class1 as ARClass
import structures_on_tree1 as StrucTree 
#import Vadim_toolbox as Vtools
import Vadim_toolbox_file as vp
#import Struct_inv_file1 as Sinv

def exchange(seq):
    '''Obtain all sequences obtained by purine-purine or 
    pyrimidine-pyrimidine exchange in seq'''
    seqAG = seq.replace('A','a').replace('G','A').replace('a','G')
    seqCT = seq.replace('C','c').replace('T','C').replace('c','T')
    seqAG_CT = seqCT.replace('A','a').replace('G','A').replace('a','G')
    
    return [seq,seqAG,seqCT,seqAG_CT]
            
def exchange_purpyr(seq):
    '''exchange purines with pyrimidines in seq'''
    seq_new = seq.replace('A','a').replace('G','g').replace('C','A').replace('T','G').replace('a','C').replace('g','T')
    return seq_new

def GTR_compact(idx_groups,n_ij,T_ia,root_states):
    '''Running simultaneous GTR calculation by considering groups of 
    nucleotides as a single nucleotide
    
    idx_groups - list of indices correspodning to each group    
    '''
    Ngroups = len(idx_groups)
    n_ij_groups = np.zeros((5,5,Ngroups))
    T_i_groups = np.zeros((5,Ngroups))
    root_states_groups = np.zeros((5,Ngroups))
    for jgroup,idx_group in enumerate(idx_groups):
        n_ij_groups[:,:,jgroup] = n_ij[:,:,idx_group].sum(axis =2)
        T_i_groups[:,jgroup] = T_ia[:,idx_group].sum(axis=1) 
        root_states_groups[:,jgroup] = root_states[:,idx_group].mean(axis=1) 
    
    W_ij_groups, p_i_groups, mu_groups = GTR.GTR_simultaneous(n_ij_groups,T_i_groups,root_states_groups,range(Ngroups))
    p_i_rec = np.zeros((5,n_ij.shape[2]))
    mu_rec = np.zeros(n_ij.shape[2])
    for jgroup, idx_group in enumerate(idx_groups):
        p_i_rec[:,idx_group] = np.tile(p_i_groups[:,jgroup],(len(idx_group),1)).T
        mu_rec[idx_group] = mu_groups[jgroup]
        
    return W_ij_groups, p_i_rec, mu_rec
    
if __name__=="__main__":
    '''Accelerating GTR calculation by dividing nucleotide frequencies into 
    classes and using nucleotide re-ordering'''
    h= 10**(-8)
#    vp=Vtools.Vadim_toolbox()
    np.set_printoptions(precision=4,suppress=True)
#    plt.ioff()
    
    New_gene = False
    if New_gene:                
        '''loading alignment'''    
        subtype= 'HIV1BC_' #'HIV1ALL_' #
        gene_folder = 'vif' #'gp41' #'env' #'vpr' #'gag' #'nef' #'vif' #'RRE' #'gag' #'pol1'
        gene = 'vif_nuc' #'gp41_nuc' #'envcut1_nuc'# 'vpr_nuc' #'gag_nuc' #'nef_nuc' #'vif_nuc' #'RRE_nuc' #'gag_nuc' # 'pol1' #
        print subtype+gene
        
        with open('/ebio/ag-neher/share/users/vpuller/Maximum_likelihood/ML_test_sequences/'+gene_folder + '/'+subtype+gene + '.fasta','r') as sequence_file:
            aln=AlignIO.read(sequence_file, 'fasta')
        
        aln_array = np.array(aln)
        nuc_freq=vp.nucleotide_frequency(aln_array)     
        nuc_freq_max=np.max(nuc_freq,axis=0)
        cons_seq_array=vp.numbers_to_nucs(np.argmax(nuc_freq,axis=0))
        cons_seq_string=''.join(cons_seq_array) 
        
        '''importing shape values'''
        SHAPE_data_file = '/ebio/ag-neher/share/users/vpuller/Maximum_likelihood/ML_test_sequences/structures_from_SHAPE/SHAPEs_2014.txt'
        SHAPE_seq_array = np.loadtxt(SHAPE_data_file, dtype = '|S1', skiprows=1,usecols = (1,))
        SHAPE_seq_string = ''.join(SHAPE_seq_array)
        SHAPE_data = np.loadtxt(SHAPE_data_file, skiprows=1,usecols = (2,))
    
        idx = np.where(nuc_freq[4,:] < .8)[0]
        if subtype == 'HIV1BC_':
            if gene == 'gag_nuc':
                idx_to_delete = [351,352,353,370,371,373]
                idx = np.delete(idx,idx_to_delete)
            elif gene == 'nef_nuc':
                idx_to_delete = [66,67,68,195,196,197]
                idx = np.delete(idx,idx_to_delete)
            elif gene == 'pol_nuc':
                idx_to_delete = [81,82,83,126,127,128]
                idx = np.delete(idx,idx_to_delete)
            elif gene == 'gp41_nuc':
                idx_to_delete =range(825,846)
                idx = np.delete(idx,idx_to_delete)
        
        elif subtype == 'HIV1ALL_':
            if gene == 'gag_nuc':
                idx_to_delete = []
                idx_to_delete.extend(range(1494,1497))
                idx_to_delete.extend(range(1479,1491))
                idx_to_delete.extend(range(1377,1383))
                idx_to_delete.extend(range(1146,1149))
                idx_to_delete.extend(range(381,384))
                idx_to_delete.extend(range(372,378))
                idx_to_delete.extend(range(363,366))
                idx_to_delete.extend(range(327,345))
                idx = np.delete(idx,idx_to_delete)
            elif gene == 'pol_nuc':
                idx_to_delete = []
                idx_to_delete.extend(range(78,84))
                idx = np.delete(idx,idx_to_delete)
            elif gene == 'nef_nuc':
                idx_to_delete = []
                idx_to_delete.extend(range(189,192))
                idx = np.delete(idx,idx_to_delete)
            
        idx_weeks = np.copy(idx)
        
        jnuc0, score = vp.seq_match(SHAPE_seq_string,''.join(cons_seq_array[idx]))
        jnucN = jnuc0 + len(idx)
        
        SHAPE = np.zeros(len(cons_seq_string))
        SHAPE[idx] = SHAPE_data[jnuc0:jnucN]
        SHAPE_seq_gapped_array = np.repeat('-',len(cons_seq_string))
        SHAPE_seq_gapped_array[idx] = SHAPE_seq_array[jnuc0:jnucN]
        SHAPE_codons, SHAPE_aas = vp.nucseq_to_codons(SHAPE_seq_gapped_array,0,len(SHAPE_seq_gapped_array),0)
        
    
        '''loading tree'''
        new_tree_file_name='/ebio/ag-neher/share/users/vpuller/Maximum_likelihood/ML_test_sequences/'+gene_folder + '/'+subtype+gene+'_new_tree0mid.nwk'
        bio_tree = Phylo.read(new_tree_file_name, 'newick')
        AR=ARClass.ancestor_reconstruction_class(bio_tree,aln,use_branch_lengths=True)
     
        '''loading reconstructed sequences'''
        tree_seq_file_name='/ebio/ag-neher/share/users/vpuller/Maximum_likelihood/ML_test_sequences/'+gene_folder + '/'+subtype+gene+'_tree0mid_seq.fasta'
            
        with open(tree_seq_file_name,'r') as tree_seq_file:
            aln_all = AlignIO.read(tree_seq_file,'fasta')    
    
        aln_all_array=np.array(aln_all)
        
        
        '''initializing the structures_on_tree class'''    
        rf=0
        nuc0=0
        nucN=aln_all_array.shape[1]
        struct_on_tree = StrucTree.structures_on_tree(AR.tree,aln_all,[],structure=[nuc0,nucN],reference_frame = rf, warnings = 'off')     
        top_node = struct_on_tree.tree.root
            
            
        import GTR_class
        GTR = GTR_class.GTR_class(top_node, warnings = 'off',Nit_max = 10**4)
        AA_array, Codons_array, Nstop = GTR.nucs_to_aas_codons(aln_array)
        AA_freq=vp.codon_frequency(AA_array,GTR.aa_bases)
        AA_freq_max=np.max(AA_freq,axis=0)    
        AA_cons_array = GTR.aa_bases[np.argmax(AA_freq,axis=0)]
        Codon_freq=vp.codon_frequency(Codons_array,GTR.codons)
        
        p_aln = vp.nucleotide_frequency(aln_array)
        p_aln_all = vp.nucleotide_frequency(aln_all_array)
#        sequence_statistics(AA_cons_array,GTR.aa_bases)
        GTR.substitution_count()
            
        W_ij, p_a, mu_a = GTR.GTR_nucs(np.where(cons_seq_array != '-')[0])
        ssigma_a = p_a**2/(GTR.n_ij.sum(axis = 1)+GTR.root_states+h)  
        S_a = -((h+p_a)*np.log(h+p_a)).sum(axis=0)
    
    dir_name = '/ebio/ag-neher/share/users/vpuller/GTR_staggered/'
#    plt.figure(10,figsize=(20,10))
#    plt.clf()        
#    plt.plot(S_a)
#    plt.plot(np.sort(S_a))
#    plt.ylabel('Entropy')   
#    plt.legend(('unsorted','sorted'))

    p_sort = np.sort(p_a,axis=0)
    
    p_aln_sort = np.sort(p_aln,axis=0)
    p_aln_all_sort = np.sort(p_aln_all,axis=0)
    plt.figure(40,figsize=(20,10))
    plt.clf()
    for j in xrange(5):    
        plt.plot(p_sort[-j-1,np.argsort(p_sort[-1,:])])
    plt.ylabel('p_sort') 
    plt.savefig(dir_name + 'nuc_freq_sorted.pdf')
    plt.close(40)
    
    plt.figure(45,figsize=(20,10))
    plt.clf()
    plt.plot(p_sort[-1,np.argsort(p_sort[-1,:])])
    plt.plot(p_aln_sort[-1,np.argsort(p_sort[-1,:])])
    plt.plot(p_aln_all_sort[-1,np.argsort(p_sort[-1,:])])
    plt.title('Dominant nuc. frequency')
    plt.legend(('GTR','alignment','tree'))
    plt.ylabel('p_sort') 
    plt.savefig(dir_name + 'Dominant_nuc_freq.pdf')
    plt.close(45)

    p_pur = p_a[0,:] + p_a[2,:]
    p_pyr = p_a[1,:] + p_a[3,:]    
    jjpur = np.argsort(p_pur)
    plt.figure(50,figsize=(20,10))
    plt.clf()        
    plt.plot(p_pur[jjpur])
    plt.plot(p_pyr[jjpur])
    plt.plot(p_a[4,jjpur])
    plt.ylabel('p_sort') 
    plt.legend(('purines: A,G', 'pyrimidines: C,T','gaps'))
    plt.savefig(dir_name + 'Pur_Pyr.pdf')
    plt.close(50)
    
    
    '''Splitting into classes by nucleotide order'''
    nuc_order = np.array([''.join(GTR.bases[p[::-1]]) for p in np.argsort(p_a,axis=0).T])
    Perms = vp.permute_str('ACGT-')     
    counts = np.array([np.count_nonzero(nuc_order == perm) for perm in Perms],dtype = 'int')
    for j in np.argsort(counts)[::-1]:
        if counts[j]>0:
            print '\n',Perms[j],counts[j],'\n', p_sort[::-1,np.where(nuc_order == Perms[j])[0]].T

    Perms1 = []
    Perms2 = []
    for perm in Perms:
        if perm not in Perms1:
            Perms1.extend(exchange(perm))
            Perms2.append(exchange(perm))
            
    counts2 = np.array([np.array([np.count_nonzero(nuc_order == perm) for perm in perm2]).sum() for perm2 in Perms2],dtype = 'int')
    for j in np.argsort(counts2)[::-1]:
        print Perms2[j],counts2[j]
        


    '''Splitting into "natural" classes'''
    idx0 = np.where(np.max(p_aln,axis=0) == 1)[0] # Class in which no mutations occur
    idx_else = np.where(np.max(p_aln,axis=0) < 1)[0]
    n_ij_groups = np.zeros((5,5,2))
    n_ij_groups[:,:,0] = GTR.n_ij[:,:,idx0].sum(axis =2)
    n_ij_groups[:,:,1] = GTR.n_ij[:,:,idx_else].sum(axis =2)
    T_i_groups = np.zeros((5,2))
    T_i_groups[:,0] = GTR.T_ia[:,idx0].sum(axis=1) 
    T_i_groups[:,1] = GTR.T_ia[:,idx_else].sum(axis=1)
    root_states_groups = np.zeros((5,2))
    root_states_groups[:,0] = GTR.root_states[:,idx0].mean(axis=1) 
    root_states_groups[:,1] = GTR.root_states[:,idx_else].mean(axis=1)

    GTR.Nit_max = 10**4    
    W_ij_groups, p_i_groups, mu_groups = GTR.GTR_simultaneous(n_ij_groups,T_i_groups,root_states_groups,range(2))
    print 'W_ij =\n',W_ij_groups,'\np_i=\n',p_i_groups,'\nmu=\n',mu_groups,'\n'
    
    pcut = .999    
#    nuc_order_GTR = np.array([''.join(GTR.bases[p[::-1]]) for p in np.argsort(p_a,axis=0).T])
    nuc_order_aln = np.array([''.join(GTR.bases[p[::-1]]) for p in np.argsort(p_aln,axis=0).T])
    
    
    '''Splitting sequence into classes WITHOUT nucleotide permutations'''
#    Perms_aln = [perm for perm in Perms if np.count_nonzero(nuc_order_aln == perm) > 0]
#    idx_groups = [np.where(nuc_order_aln == perm)[0] for perm in Perms if np.count_nonzero(nuc_order_aln == perm) > 0]
#    Perms_aln = [perm for perm in Perms if np.count_nonzero(nuc_order == perm) > 0]
    idx_groups = [np.where(nuc_order == perm)[0] for perm in Perms if np.count_nonzero(nuc_order == perm) > 0]    
    W_ij_groups, p_i_rec, mu_rec = GTR_compact(idx_groups,GTR.n_ij,GTR.T_ia,GTR.root_states)
    
    p_i_rec_sort = np.sort(p_i_rec,axis=0)
    
    
    '''Splitting sequence into classes WITH nucleotide permutations'''
    group_swap = np.array([np.where(np.array(Perms2) == order) for order in nuc_order]).squeeze()
    idx_swap = [np.where(group_swap[:,0] == j)[0] for j in xrange(len(Perms2)) if np.count_nonzero(group_swap[:,0] == j) >0]
    
    n_ij_swap = np.zeros(GTR.n_ij.shape, dtype = 'int')
    T_ia_swap = np.zeros(GTR.T_ia.shape)
    root_states_swap = np.zeros(GTR.root_states.shape)
    
    idx_0 = np.where(group_swap[:,1] == 0)[0]
    if len(idx_0) >0:
        n_ij_swap[:,:,idx_0] = GTR.n_ij[:,:,idx_0]
        T_ia_swap[:,idx_0] = GTR.T_ia[:,idx_0]
        root_states_swap[:,idx_0] = GTR.root_states[:,idx_0]
        
    idx_AG = np.where(group_swap[:,1] == 1)[0]
    if len(idx_AG) >0:
        n_ij_swap[:,:,idx_AG] = GTR.n_ij[[2,1,0,3,4],:,:][:,[2,1,0,3,4],:][:,:,idx_AG]
        T_ia_swap[:,idx_AG] = GTR.T_ia[[2,1,0,3,4],:][:,idx_AG]
        root_states_swap[:,idx_AG] = GTR.root_states[[2,1,0,3,4],:][:,idx_AG]
        
    idx_CT = np.where(group_swap[:,1] == 2)[0]
    if len(idx_CT) >0:
        n_ij_swap[:,:,idx_CT] = GTR.n_ij[[0,3,2,1,4],:,:][:,[0,3,2,1,4],:][:,:,idx_CT]
        T_ia_swap[:,idx_CT] = GTR.T_ia[[0,3,2,1,4],:][:,idx_CT]
        root_states_swap[:,idx_CT] = GTR.root_states[[0,3,2,1,4],:][:,idx_CT]
        
    idx_AG_CT = np.where(group_swap[:,1] == 3)[0]
    if len(idx_AG_CT) >0:
        n_ij_swap[:,:,idx_AG_CT] = GTR.n_ij[[2,3,0,1,4],:,:][:,[2,3,0,1,4],:][:,:,idx_AG_CT]
        T_ia_swap[:,idx_AG_CT] = GTR.T_ia[[2,3,0,1,4],:][:,idx_AG_CT]
        root_states_swap[:,idx_AG_CT] = GTR.root_states[[2,3,0,1,4],:][:,idx_AG_CT]
    
    GTR.Nit_max = 10**5      
    W_ij_swap, p_i_swap, mu_swap = GTR_compact(idx_swap,n_ij_swap,T_ia_swap,root_states_swap)
    
    if len(idx_AG) >0:
        p_i_swap[:,idx_AG] = p_i_swap[[2,1,0,3,4],:][:,idx_AG]
        
    if len(idx_CT) >0:
        p_i_swap[:,idx_CT] = p_i_swap[[0,3,2,1,4],:][:,idx_CT]
        
    if len(idx_AG_CT) >0:
        p_i_swap[:,idx_AG_CT] = p_i_swap[[2,3,0,1,4],:][:,idx_AG_CT]

    p_i_swap_sort = np.sort(p_i_swap,axis=0)   
    
    idx_sort = np.argsort(p_sort[-1,:])
    plt.figure(100,figsize=(20,10))
    plt.clf()        
    plt.plot(p_sort[-1,idx_sort])
    plt.plot(p_i_rec_sort[-1,idx_sort])
    plt.plot(p_i_swap_sort[-1,idx_sort])
    plt.ylabel('p_sort')   
    plt.legend(('GTR','simple classes','swap classes'),loc=4)
    plt.savefig(dir_name + 'simple_ACGT_classes.pdf')    
    
    idx_sort = np.argsort(p_i_rec_sort[-1,:])
    plt.figure(110,figsize=(20,10))
    plt.clf()        
    plt.plot(p_sort[-1,idx_sort])
    plt.plot(p_i_rec_sort[-1,idx_sort])
    plt.plot(p_i_swap_sort[-1,idx_sort])
    plt.ylabel('p_sort')   
    plt.legend(('GTR','simple classes','swap classes'),loc=4)
#    plt.savefig(dir_name + 'simple_ACGT_classes.pdf')  

    idx_sort = np.argsort(p_i_swap_sort[-1,:])
    plt.figure(120,figsize=(20,10))
    plt.clf()        
    plt.plot(p_sort[-1,idx_sort])
    plt.plot(p_i_rec_sort[-1,idx_sort])
    plt.plot(p_i_swap_sort[-1,idx_sort])
    plt.ylabel('p_sort')   
    plt.legend(('GTR','simple classes','swap classes'),loc=4)
#    plt.savefig(dir_name + 'simple_ACGT_classes.pdf')  