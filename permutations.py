# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 17:17:17 2015

@author: vpuller
"""

def permute_list(L):
#    print L
    '''Permutations of the elements of list L of different elements'''
    if len(L) == 2:
        perms = [L,L[::-1]] 
    else:
        perms = []
        for element in L:
            Lred = list(L)
            Lred.remove(element)
            for Lred_perm in permute_list(Lred):
                Lnew = [element]
                Lnew.extend(Lred_perm)
                perms.append(Lnew)
    
    return perms
            
def permute_str(L):
#    print L
    '''Permutations of the elements of string L of different elements'''
    if len(L) == 2:
        perms = [L,L[::-1]] 
    else:
        perms = []
        for element in L:
            Lred = L.replace(element,'')
            for Lred_perm in permute_str(Lred):
                Lnew = element + Lred_perm
                perms.append(Lnew)
    
    return perms
    
if __name__=="__main__":
#    L = ['a','b']
#    permutations(L)
#    L = ['a','b','c']
#    print permute_list(L)
    L = 'ACGT-'
    Perms = permute_str(''.join(L))
    print Perms
    