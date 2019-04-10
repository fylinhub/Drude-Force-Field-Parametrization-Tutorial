#! /usr/bin/python
# compute molecular volume
# Vm = <VOLUme> / number of molecules in the box
# Jing Huang, Oct. 2014

import string,sys,math

try:
        xfile = open(sys.argv[1],'r')
except IOError:                     #catch exception
        print ('Files do not exist!\n')


number=int(sys.argv[2])


x1=[]
for line in xfile.readlines():
    x1.append(float(string.split(line)[0])**3)

vbox=sum(x1)/len(x1)

print vbox/float(number)


