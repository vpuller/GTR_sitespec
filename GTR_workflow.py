# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 14:52:16 2017

@author: vpuller
"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import pydot
    
if __name__=="__main__":
    '''Making graph of workflow in GTR reconstruction project'''
    plt.ioff()
    plt.close('all')
    
    dir_name = '/ebio/ag-neher/share/users/vpuller/GTR_staggered/'
    G = pydot.Dot(graph_type = 'digraph')
    Gsimul = pydot.Subgraph(graph_type = 'digraph', rank = 'min')
    Gsimul.add_node(pydot.Node(name = 'GTR', label = 'GTR model\n Wij, mu_a, p_ia', shape = 'box'))
    Gsimul.add_node(pydot.Node(name = 'tr0', label = 'Tree', shape = 'box'))
    G.add_subgraph(Gsimul)
    G.add_node(pydot.Node(name = 'evol', label = 'evolving sequences'))
    
    Gstrip = pydot.Subgraph(graph_type = 'digraph', rank = 'same')    
    Gstrip.add_node(pydot.Node(name = 'tr_dress', label = 'Dressed tree', shape = 'box'))
    Gstrip.add_node(pydot.Node(name = 'tr_strip', label = 'Stripped tree', shape = 'box'))
    Gstrip.add_node(pydot.Node(name = 'aln', label = 'Alignment', shape = 'box'))
    G.add_subgraph(Gstrip)
    
    G.add_node(pydot.Node(name = 'tr_rec', label = 'tree reconstruction'))
    G.add_node(pydot.Node(name = 'anc_rec', label = 'ancestral reconstruction'))
    G.add_node(pydot.Node(name = 'gtr_rec', label = 'GTR reconstruction'))
    
    Gresult = pydot.Subgraph(graph_type = 'digraph', rank = 'max')
    G.add_node(pydot.Node(name = 'GTR_rec', label = 'Reconstructed GTR model', shape = 'box'))
    G.add_subgraph(Gresult)
    
    G.add_edge(pydot.Edge('GTR','evol'))
    G.add_edge(pydot.Edge('tr0','evol'))
    G.add_edge(pydot.Edge('evol','tr_dress'))
    G.add_edge(pydot.Edge('tr_dress','gtr_rec'))
    G.add_edge(pydot.Edge('tr_dress','tr_strip'))
    G.add_edge(pydot.Edge('tr_strip','anc_rec'))
    G.add_edge(pydot.Edge('anc_rec','gtr_rec'))
    G.add_edge(pydot.Edge('tr_strip','aln'))
    G.add_edge(pydot.Edge('aln','tr_rec'))
    G.add_edge(pydot.Edge('tr_rec','anc_rec'))
    G.add_edge(pydot.Edge('gtr_rec','GTR_rec'))
    G.write_pdf(dir_name + 'GTR_workflow.pdf')
    