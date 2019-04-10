#!/usr/bin/python
# codes to perform nma/alad + sod/water interaction energy calculation
# process the Q-CHEM BSSE results
# Oct 2014, Jing Huang

import sys, string, copy
# import numpy as np
import math
from math import sqrt, pi, cos

def ptinker(s1,s2,r,a,d):
    filename=s1+"."+s2+"."+str(r)+"_"+str(a)+"_"+str(d)
    try:
        tfile = open(filename+".out",'r')
    except IOError:                     #catch exception
        print ('File '+filename+' do not exist!\n')

    for line in tfile.readlines():
        ll=string.split(line)
        if len(ll)==4:
            if(ll[0]=="DE," and ll[1]=="kJ/mol"): 
                e=float(ll[3])/4.184 # from kJ/mol to kcal/mol

    print r, a, d, e


if __name__ == "__main__":
    resi="nma"
    resi2="sod"
    
    r1=2.4  
    ptinker(resi,resi2,r1,180,180)
    for a1 in range(165,75,-15):
        for d1 in range(180,-180,-15):
            ptinker(resi,resi2,r1,a1,d1)

    
    
    
