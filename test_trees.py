# -*- coding: utf-8 -*-
"""
Created on Tue May 17 18:16:32 2016

@author: vpuller
"""
from __future__ import division
#import numpy as np
from Bio import Phylo
import sys, os
import matplotlib.pyplot as plt

sys.path.append('/ebio/ag-neher/share/users/vpuller/betatree/')
from betatree import betatree as bt

if __name__=="__main__":
    '''Generate test trees'''
    plt.ioff()
    outdir_name = '/ebio/ag-neher/share/users/vpuller/GTR_staggered/trees/'
    if not os.path.exists(outdir_name):
        os.makedirs(outdir_name)    
    
    nleaves = 1000
    myT = bt.betatree(nleaves,alpha = 2)
    myT.yule = True
    myT.coalesce()
    myT.BioTree.ladderize()


    for jclade, clade in enumerate(myT.BioTree.get_terminals()):
        clade.name = 'A' + clade.name
            
    Phylo.write(myT.BioTree, outdir_name + 'test_tree{}.nwk'.format(nleaves), 'newick')
    
    fig10=plt.figure(10,figsize=(20,20))    
    plt.clf()
    ax10=fig10.add_subplot(1,1,1)
    Phylo.draw(myT.BioTree,do_show = False,axes = ax10)
    plt.savefig(outdir_name + 'test_tree.pdf')
    plt.close(10)
    