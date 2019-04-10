#! /usr/bin/python
# compute Hval value (in kcal/mol)
# Usage; hval.py [data file for gas phase potential energy] [data file for condense phase potential energy] [number of molecules in the box] [simulation temperature]
# heat of vaporization dHvap = <gas phase ENERgy> - <condense phase ENERgy> / # molecules in the box + R*T where T is the simulation temperature
# Jing Huang, Oct. 2014

import string,sys,math

try:
        xfile = open(sys.argv[1],'r')
except IOError:                     #catch exception
        print ('Files do not exist!\n')

try:
        yfile = open(sys.argv[2],'r')
except IOError:                     #catch exception
        print ('Files do not exist!\n')

number=int(sys.argv[3])
temp=float(sys.argv[4])


x1=[]
for line in xfile.readlines():
    x1.append(float(string.split(line)[0]))

y1=[]
for line in yfile.readlines():
    y1.append(float(string.split(line)[0]))

gasene=sum(x1)/len(x1)
boxene=sum(y1)/len(y1)

R = 8.31441/4184.0
print gasene-boxene/float(number)+R*temp


