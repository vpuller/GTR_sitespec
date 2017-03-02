# -*- coding: utf-8 -*-
"""
Created on Wed May 18 11:17:35 2016

@author: vpuller
"""
from __future__ import division
import numpy as np
from scipy import linalg as LA
from scipy.optimize import minimize
import sys

sys.path.insert(0, '/ebio/ag-neher/share/users/vpuller/myTOOLS/') 
import Vadim_toolbox_file as vp

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
    
def nuc_sub_count(tree, alphabet = 'ACGT', fpar = 0.5, use_pro = False):        
        '''
        Function that counts the numbers of nucleotide substitutions on tree
        
        Input arguments:
        tree - tree dressed withs equences
        alphabet - nucleotide alphabet
        fpar - fraction of the branch length assigned to the parent nucletide state
        '''
#        print alphabet
        q = len(alphabet); L = len(tree.root.seq)
        substitutions = np.zeros((q,q,L))
        T_ij = np.zeros((q,q,L))
        for parent in tree.get_nonterminals():
            seq_parent_num=vp.string_to_numbers(parent.seq,alphabet = alphabet)
            for child in parent:
                if child.branch_length >0:
                    if not use_pro:
                        seq_child_num=vp.string_to_numbers(child.seq,alphabet = alphabet)
                        T_ij[seq_child_num,seq_parent_num,range(L)] += child.branch_length 
                        substitutions[seq_child_num,seq_parent_num,range(L)] +=1  
                    else:
                        T_ij += child.corr*child.branch_length                        
                        substitutions += child.corr
        
        if use_pro:
            root_states = np.copy(tree.root.P)   
        else:
            root_states = seq_to_P(tree.root.seq,alphabet = alphabet)        
        n_ij = np.copy(substitutions)     
        n_ij[range(q),range(q),:] = 0
        return n_ij, T_ij, root_states

        
def GTR_simult(sub_ija, T_ija,root_states, dp = 10**(-4), Nit = 10**4):
    '''
    Calculating components of GTR matrix for multiple sites
    
    Input arguments:
    sub_ija - transition counts matrix
    T_ia - times spent in each nucleotide state
    root_states - the nucleotide states of the tree root (for pseudocount correction)
    idx_cum - indices of the nucleotides to use
    
    Returns:
    attempt matrix (which sums to 1), equilibrium populations by site, mutation rates
    '''
#    q = n_ija.shape[0]; L = n_ija.shape[2]
    (q,L) = root_states.shape
    n_ija = np.copy(sub_ija)
    n_ija[range(q),range(q),:] = 0
    Lambda = np.sum(root_states,axis=0)
    T_ia = np.sum(T_ija + np.transpose(T_ija, axes = (1,0,2)), axis = 1)/2.
#    T_a = T_ia.sum(axis=0)
    dT_ija = np.copy(T_ija); dT_ija[range(q),range(q),:] = 0. 
    
    p_ia_old=np.zeros((q,L))    
    mu_a = np.ones(L)
    p_ia = np.ones((q,L))/q
    
    n_ij = n_ija.sum(axis=2)
    n_ia = n_ija.sum(axis=1)
    n_a = n_ija.sum(axis=1).sum(axis=0)
    
    W_ij = np.ones(q) - np.eye(q)
    count = 0
    while (LA.norm(p_ia_old-p_ia) > dp or np.abs(1-np.max(p_ia.sum(axis=0))) > dp) and count < Nit:
#        import ipdb; ipdb.set_trace()
        count += 1
        p_ia_old = np.copy(p_ia)
#        W_ij = (n_ija+np.transpose(n_ija,axes=(1,0,2))).sum(axis = 2)/((mu_a*p_ia).dot(T_ia.T)+T_ia.dot((mu_a*p_ia).T)+h)
#        W_ij = W_ij/np.sum(W_ij)
#        p_ia = (np.sum(n_ija,axis=1)+root_states)/(mu_a*np.dot(W_ij,T_ia)+Lambda)
#        mu_a = n_ija.sum(axis=1).sum(axis=0)/(h+np.diag(p_ia.T.dot(W_ij).dot(T_ia)))
        
#        S_ij = (mu_a*p_ia).dot(T_ia.T)
        TW = np.transpose(dT_ija/(W_ij[:,:,np.newaxis]+h),(2,0,1))
        TWW = np.transpose(dT_ija/(W_ij[:,:,np.newaxis]**2+h),(2,0,1))
        TT_ija = T_ia[np.newaxis,:,:] -\
        0.5*(np.transpose(TW.dot(W_ij),(2,1,0)) + np.transpose(W_ij.dot(TW),(0,2,1)) -\
        np.transpose(W_ij.dot(TWW),(0,2,1))*W_ij[:,:,np.newaxis])
        
        TT_ija = T_ia[np.newaxis,:,:]
        S_ij = np.sum(mu_a *p_ia[:,np.newaxis,:]*TT_ija, axis = 2)
        W_ij = (n_ij+n_ij.T)/(S_ij + S_ij.T+h)
        W_ij = W_ij/np.sum(W_ij)
        
        Q_ia = np.dot(W_ij,T_ia) - .5*W_ij.dot(TW).dot(W_ij)[range(q),:,range(q)]
        p_ia = (n_ia+root_states)/(mu_a*Q_ia + Lambda)

#        R_a = np.diag(p_ia.T.dot(W_ij).dot(T_ia))
#        R_a = np.sum(p_ia*W_ij.dot(T_ia),axis = 0)
        R_a = np.sum(p_ia*Q_ia,axis = 0)   
        mu_a = n_a/(h + R_a)
#        mu_a = (n_a + h/(np.max(W_ij)np.sqrt())/(h + R_a)
        
    if count >= Nit:
        print 'WARNING: maximum number of iterations has been reached in GTR_simult',\
        np.min(p_ia.sum(axis=0)), np.max(p_ia.sum(axis=0))
        if LA.norm(p_ia_old-p_ia) > dp:
            print '    the iterative scheme has not converged'
        elif np.abs(1-np.max(p_ia.sum(axis=0))) > dp:
            print '    the iterative scheme has converged, but proper normalization was not reached'
            
    return W_ij, p_ia, mu_a

def GTR_oneforall(sub_ija, T_ija,root_states, dp = 10**(-4), Nit = 10**4):
    '''
    Calculating GTR matrix assuming that it is the same for all sites
    Input arguments:
    sub_ija - transition counts matrix
    T_ia - times spent in each nucleotide state
    root_states - the nucleotide states of the tree root (for pseudocount correction)
    
    Returns:
    attempt matrix, equilibrium populations by site
    '''
    (q,L) = root_states.shape
    n_ij = sub_ija.sum(axis = 2)
    n_ij[range(q),range(q)] = 0
    
    m_i = root_states.sum(axis = 1)
    Lambda = m_i.sum()
#    T_ia = np.sum(T_ija + np.transpose(T_ija, axes = (1,0,2)), axis = 1)/2.
#    T_i = T_ia.sum(axis = 1)
    T_ij = T_ija.sum(axis = 2)
    T_i = (T_ij + T_ij.T).sum(axis=0)/2.
    dT_ij = T_ij - np.diag(np.diag(T_ij))

    p_old=np.zeros(q)    
    p_i = np.ones(q)/q
    
    W_ij = np.ones(q) - np.eye(q)
    count = 0
    while (LA.norm(p_old-p_i) > dp or np.abs(1-np.max(p_i.sum(axis=0))) > dp) and count < Nit:
#        print p_i
#        raw_input('continue?')
        count += 1
        p_old = np.copy(p_i)
#        W_ij = (n_ij + n_ij.T)/(np.array([p_i]).T*T_i + np.array([T_i]).T*p_i+h)
#        p_i = (n_ij.sum(axis=1) + m_i)/(W_ij.dot(T_i)+Lambda)

        S_ij = np.array([p_i]).T*T_i -.5*(dT_ij/(W_ij+h)).dot(W_ij).dot(np.diag(p_i)) -\
        0.5*np.diag(p_i).dot(W_ij).dot(dT_ij/(W_ij+h)) + \
        .5*W_ij.dot(np.diag(p_i)).dot(W_ij)*dT_ij/(W_ij+h)**2
        W_ij = (n_ij + n_ij.T)/(S_ij + S_ij.T+h)
        Q_i = W_ij.dot(T_i) - .5*np.diag(W_ij.dot(dT_ij/(W_ij+h)).dot(W_ij))
        p_i = (n_ij.sum(axis=1) + m_i)/(Q_i+Lambda)
    
    if count >= Nit:
        print 'WARNING: maximum number of iterations has been reached in GTR_oneforall,\np_i = ',\
        p_i, p_i.sum()
        if LA.norm(p_old-p_i) > dp:
            print '    the iterative scheme has not converged, ', LA.norm(p_old-p_i)
        elif np.abs(1-np.max(p_i.sum(axis=0))) > dp:
            print '    the iterative scheme has converged, but proper normalization was not reached', np.sum(p_i)
            
    return W_ij, p_i

def GTR_sitespec(n_ija, T_ija,root_states, dp = 10**(-4), Nit = 10**4):
    '''
    Calculating separate GTR matrix for every site
    Input arguments:
    sub_ija - transition counts matrix
    T_ia - times spent in each nucleotide state
    root_states - the nucleotide states of the tree root (for pseudocount correction)
    
    Returns:
    attempt matrix, equilibrium populations by site
    '''
    def GTR_singlesite(n_ij, T_ij, m_i):
        n_ij[range(q),range(q)] = 0
        Lambda = m_i.sum()
        T_i = (T_ij + T_ij.T).sum(axis=0)/2.
#        T_ii = np.diag(T_ij)
#        dT_ij = T_ij - np.diag(np.diag(T_ij))
    
        p_old=np.zeros(q)    
        p_i = np.ones(q)/q
        p_i = m_i
        
        W_ij = (np.ones(q) - np.eye(q))/(q*(q-1))
        count = 0
        while (LA.norm(p_old-p_i) > dp or np.abs(1-np.max(p_i.sum(axis=0))) > dp) and count < Nit:
            count += 1
            p_old = np.copy(p_i)
            W_ij = (n_ij + n_ij.T)/(np.array([p_i]).T*T_i + np.array([T_i]).T*p_i+h)
#            mu = W_ij.sum()/(q*(q-1))
#            W_ij = mu*(np.ones(q) - np.eye(q))
            p_i = (n_ij.sum(axis=1) + m_i)/(W_ij.dot(T_i)+Lambda)
            
#            S_ij = np.array([p_i]).T*T_i -.5*(dT_ij/(W_ij+h)).dot(W_ij).dot(np.diag(p_i)) -\
#            0.5*np.diag(p_i).dot(W_ij).dot(dT_ij/(W_ij+h)) + \
#            0.5*W_ij.dot(np.diag(p_i)).dot(W_ij)*(dT_ij/(W_ij+h)**2)
##            S_ij = np.array([p_i]).T*T_i
##            W_ij = (n_ij + n_ij.T + h*W_JC)/(S_ij + S_ij.T+h)
#            W_ij = (n_ij + n_ij.T)/(S_ij + S_ij.T+h)
#            
#            Q_i = W_ij.dot(T_i) - .5*np.diag(W_ij.dot(dT_ij/(W_ij+h)).dot(W_ij))
##            Q_i = W_ij.dot(T_i)
#            p_i = (n_ij.sum(axis=1) + m_i)/(Q_i+Lambda)
            
#            mu = n_ij.sum()/(.5*dT_ij.sum() + np.sum(p_i*(T_ii.sum() - T_ii)))
#            p_i = (n_ij.sum(axis=1) + m_i)/(mu*(T_ii.sum() - T_ii + 0.5*dT_ij.sum())+ Lambda)
#            W_ij = mu*(np.ones(q) - np.eye(q))

        if count >= Nit:
            print 'WARNING: maximum number of iterations has been reached in GTR_oneforall,\np_i = ',\
            p_i, p_i.sum()
            if LA.norm(p_old-p_i) > dp:
                print '    the iterative scheme has not converged, ', LA.norm(p_old-p_i)
#                print '    total counts: \n', n_ij
#                print '    W_ij = \n', W_ij
#                import ipdb; ipdb.set_trace()
            elif np.abs(1-np.max(p_i.sum(axis=0))) > dp:
                print '    the iterative scheme has converged, but proper normalization was not reached', np.sum(p_i)
                
        return W_ij, p_i
    
    (q,L) = root_states.shape
#    n_ij = sub_ija.sum(axis = 2)
#    m_i = root_states.sum(axis = 1)
#    T_i = T_ia.sum(axis = 1)
#    (W_ij, p_i) = GTR_singlesite(n_ij, T_i, m_i)
#    T_ia = np.sum(T_ija + np.transpose(T_ija, axes = (1,0,2)), axis = 1)/2.
    W_ija = np.zeros((q,q,L))
    p_ia = np.zeros((q,L))
#    W_JC = (np.ones(q) - np.eye(q))/(q*(q-1))
    for ja in xrange(L):
#        print ja
        (W_ija[:,:,ja], p_ia[:,ja]) = GTR_singlesite(n_ija[:,:,ja], T_ija[:,:,ja], root_states[:,ja])

    return W_ija, p_ia
    
    
def GTR_exact1(tree, alphabet = 'ACGT', amoeba = True):
    '''
    Calculating GTR matrix (same for all sites) using the exact likelihood
    
    Input arguments
    tree - tree dressed with sequences
    alphabet - nucleotide alphabet
    '''
    import tree_likelihood as TL
    tl = TL.tree_likelihood(tree, alphabet = alphabet)

    def logLike(x):
        Wij = np.zeros((q,q)); Wij[idx] = x[:-q]; Wij = Wij + Wij.T
        p_i = x[-q:]/x[-q:].sum();  
        p_a = np.tile(p_i[:,np.newaxis],(1,L))
        logL = tl.logL((Wij, p_a, mu_a))
#        print '\n{}, {}, {}'.format(Wij, p_i, logL)
        return -logL
        
    q = len(alphabet); L = len(tree.root.seq)
    Wij0 = np.ones((q,q)) - np.eye(q)
    p0_i = np.ones(q)/q
    mu_a = np.ones(L)
    idx = np.triu_indices(q, 1)
    x0 = np.zeros(q*(q+1)//2)
    x0[:-q] = Wij0[idx]
    x0[-q:] = p0_i
    
    if amoeba:
        xx = vp.amoeba_vp(logLike,x0,args=(),a = x0, tol_x = h, tol_f = h, Nit = 10**4)
    else:
        res = minimize(logLike, x0)
        xx = res.x
#        jit = 0; Nit = 100
#        dfdx = np.zeros(x0.shape)
#        f_old = logLike(x0)
#        f_new = np.copy(f_old)
#        xx = np.copy(x0)
#        nu = .01
#        for j in xrange(Nit):
#            for jx, x in enumerate(xx):
#                dx = np.zeros(xx.shape); dx[jx] = 0.01*x
#                dfdx[jx] = (logLike(xx + dx) - f_old)/dx[jx]
#            xx = xx - nu*dfdx
#            f_new = logLike(xx)
#            print j, xx, f_new, nu
#            if f_new < f_old:
#                nu = nu/2.
#            if abs(f_new/f_old - 1) < h:
#                break
#
#        if jit == Nit-1:
#            print 'WARNING from GTR_exact1: Maximum number of iterations has been reached, ', abs(f_old/f_new - 1)
            
    Wij = np.zeros((q,q)); Wij[idx] = xx[:-q]; Wij = Wij + Wij.T
    p_i = xx[-q:]
    return Wij, p_i

def GTR_exact1_MC(tree, alphabet = 'ACGT', tol = 10**(-4), Nit = 10**4, nu = 0.01, VERBOSE = False, T = 10.):
    '''
    Calculating GTR matrix (same for all sites) using the exact likelihood 
    and Monte-Carlo minimization 
    
    Input arguments
    tree - tree dressed with sequences
    alphabet - nucleotide alphabet
    '''
    import tree_likelihood as TL
    tl = TL.tree_likelihood(tree, alphabet = alphabet)

#    def logLike(Wij, p_i):
##        Wij = np.zeros((q,q)); Wij[idx] = x[:-q]; Wij = Wij + Wij.T
##        p_i = x[-q:]/x[-q:].sum();  
#        p_a = np.tile(p_i[:,np.newaxis],(1,L))
#        logL = tl.logL((Wij, p_a, mu_a))
##        print '\n{}, {}, {}'.format(Wij, p_i, logL)
#        return -logL
        
    q = len(alphabet); L = len(tree.root.seq)
    Wij = np.ones((q,q)) - np.eye(q)
    p_i = np.ones(q)/q
    mu_a = np.ones(L)
    idx = np.triu_indices(q, 1)
    
#    f_old = logLike(Wij, p_i)
    p_a = np.tile(np.array([p_i]),(L,1)).T
    f_old = tl.logL((Wij, p_a, mu_a))
    for j in xrange(Nit):
        dW = nu*np.sum(Wij)*np.random.random(size = Wij.shape)
        dW = dW + dW.T
        np.fill_diagonal(dW, 0.)
        dp0 = np.random.exponential(size = (q,2)); dp0 = dp0/np.sum(dp0,axis=0)
        dp = nu*(dp0[:,0] - dp0[:,1])
        p_a = np.tile(np.array([p_i + dp]),(L,1)).T
        f_new = tl.logL((Wij + dW, p_a, mu_a))
        
#        r = np.exp((f_new - f_old)/T)
#        r = r/(1. + r)
#        eps = np.random.uniform()
#        if VERBOSE:
#            print j, f_old, abs(f_old/f_new - 1), r, eps
#        if abs(f_old/f_new - 1) < tol:
#            break
#        elif r > eps:
#            f_old = np.copy(f_new)
#            Wij = Wij + dW
#            p_i = p_i + dp
#            print 'accept step'
        
        i = np.random.randint(0, high = q*(q+1)//2)
        dW = np.zeros((q,q)); dp = np.zeros(q)
        if i < q*(q-1)//2:
            dW[idx[0][i], idx[1][i]] += Wij[idx[0][i], idx[1][i]]*np.random.uniform()*nu
            dW = dW + dW.T
        else:
            dp0 = np.random.exponential(size = (q,2)); dp0 = dp0/np.sum(dp0,axis=0)
            dp = nu*(dp0[:,0] - dp0[:,1])
        
#        f_new = logLike(Wij + dW, p_i + dp)
        p_a = np.tile(np.array([p_i + dp]),(L,1)).T
        f_new = tl.logL((Wij + dW, p_a, mu_a))
        if VERBOSE:
            print j, f_old, abs(f_new/f_old - 1), p_i
        if abs(f_new/f_old - 1) < tol:
            break
        elif f_new > f_old:
            f_old = np.copy(f_new)
            Wij = Wij + dW
            p_i = p_i + dp
    
    if j == Nit + 1:
        print 'maximum number of iterations has been reached, ', abs(f_old/f_new - 1)
    return Wij, p_i
