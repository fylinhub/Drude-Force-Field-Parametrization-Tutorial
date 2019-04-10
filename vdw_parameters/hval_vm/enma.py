#! /usr/bin/python
# compute dielectric constant for NMA
# e = 
# Jing Huang, Oct. 2014

import string,sys,math

einf=1.72

try:
        xfile = open(sys.argv[1],'r')
except IOError:                     #catch exception
        print ('Files do not exist!\n')

try:
        yfile = open(sys.argv[2],'r')
except IOError:                     #catch exception
        print ('Files do not exist!\n')

temp=float(sys.argv[3])


mx=[]
my=[]
mz=[]
for line in xfile.readlines():
    mx.append(float(string.split(line)[1]))
    my.append(float(string.split(line)[2]))
    mz.append(float(string.split(line)[3]))

y1=[]
for line in yfile.readlines():
    y1.append((float(string.split(line)[0]))**3)

volume=sum(y1)/len(y1)

kb=1.3806503 #*10^(-23), Boltzman constant   
ke=8.987551787368 #*10^9, Coulomb's constant
C=1.60217646 #*10^(-19), electron charge to coulombs
debye2ea=0.20819434
convt=10000*C**2*ke/kb
convt=debye2ea*debye2ea*convt*4.0/3.0*3.1415926

# convt=30339.2504016 # slightly different from 30339.05 that epsilon_combine.pl uses

a=0.0
aa=0.0
for x in mx:
    aa += x*x
    a += x
aax=aa/len(mx)
ax=a/len(mx)

a=0.0
aa=0.0
for x in my:
    aa += x*x
    a += x
aay=aa/len(my)
ay=a/len(my)

a=0.0
aa=0.0
for x in mz:
    aa += x*x
    a += x
aaz=aa/len(mz)
az=a/len(mz)

#print aax, ax, aay, ay, aaz, az
dm=aax-ax*ax+aay-ay*ay+aaz-az*az

e=einf+dm*convt/volume/temp
print e




