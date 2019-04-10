#!/usr/bin/python
# codes to perform nma/alad + sod/water interaction energy calculation
# part1: prepare for calculations using Q-chem and tinker
# Oct 2014, Jing Huang

import sys, string, copy
# import numpy as np
import math
from math import sqrt, pi, cos

def parse_crd(crdname):
    try:
        crd = open(crdname+".crd",'r')
    except IOError:                     #catch exception
        print ('Files do not exist!\n')
        
        
    buffer = crd.readline()
    while buffer[0:1] == "*":
        buffer = crd.readline()
    atomnum = int(buffer)
    atoms = []
    for i in range(atomnum):
        line=string.split(crd.readline())
        atoms.append([line[3][0],float(line[4]),float(line[5]),float(line[6])])
    return atoms

def convert(s1,s2,r,a,d):
    filename=s1+"."+s2+"."+str(r)+"_"+str(a)+"_"+str(d)
    coor= parse_crd(filename)
    if s1=="nma" and s2=="sod":
        write_nma_sod(coor,filename)
    if s1=="nma" and s2=="tip3":
        write_nma_tip3(coor,filename)
    if s1=="nma" and s2=="pot":
        write_nma_pot(coor,filename)

def write_nma_tip3(atoms,crdname):
    qchem = open("opt."+crdname+".inp",'w')
    subq = open("run."+crdname+".sh",'w')
    xyz = open(crdname+".xyz",'w')

    # write qchem input files:
    qchem.write("$molecule\n0 1\n")
    for i in range(15):
        s=' %s%12.5f%12.5f%12.5f\n' % (atoms[i][0],atoms[i][1], atoms[i][2], atoms[i][3])
        qchem.write(s)
    qchem.write("$end\n")
    qchem.write("\n")
    qchem.write("$rem\n")
    qchem.write("JOBTYPE                OPT\n")
    qchem.write("EXCHANGE               HF\n")
    qchem.write("CORRELATION            RIMP2\n")
    qchem.write("BASIS                  CC-PVQZ\n")
    qchem.write("AUX_BASIS              RIMP2-CC-PVQZ\n")
    qchem.write("TIME_MP2               TRUE\n")
    qchem.write("PURECART               1111\n")
    qchem.write("N_FROZEN_CORE          FC\n")
    qchem.write("MEM_STATIC             1000\n")
    qchem.write("SCF_CONVERGENCE        8\n")
    qchem.write("THRESH                 12\n")
    qchem.write("symmetry = false\n")
    qchem.write("max_scf_cycles = 50\n")
    qchem.write("xc_grid = 000075000302\n")
    qchem.write("$end\n\n")
    qchem.write("$opt\n")
    qchem.write("FIXED\n")
    for i in range(1,14):
        s=str(i)+" XYZ\n"
        qchem.write(s)
    qchem.write("ENDFIXED\n")
    qchem.write("$end\n")
    

# write qchem submission files on kickout:
    subq.write("#$ -S /bin/bash\n")
    subq.write("#$ -cwd\n")
    subq.write("#$ -V\n")
    subq.write("#$ -e /home/huangj/jobreport\n")
    subq.write("#$ -o /home/huangj/jobreport\n")
    subq.write("#$ -j y\n")
    s="#$ -N "+crdname+"\n"
    subq.write(s)
    subq.write("#$ -pe mpirun1 1\n")
    subq.write("#$ -l h_data=1600M,h_rt=300:00:00\n")
    subq.write("#$ -R y\n")
    subq.write("ulimit -c 4\n\n")
    subq.write("QC=/opt/mackerell/apps/qchem/4.0.1\n")
    subq.write(". $QC/qcenv.sh\n")
    subq.write("QCSCRATCH=`pwd`\n")
    subq.write("QCLOCALSCR=/tmp/$USER/tmp$$\n")
    subq.write("if [ ! -d $QCLOCALSCR ] ; then mkdir -p $QCLOCALSCR; fi\n")
    subq.write("export QC QCSCRATCH QCLOCALSCR\n\n")
    s="$QC/bin/qchem opt."+crdname+".inp opt."+crdname+".out\n"
    subq.write(s)
    
    #write tinker xyz files
    xyz.write("    15\n")
    s="     1  CL %12.6f%12.6f%12.6f   130     2     3     4     5\n"  % (atoms[0][1], atoms[0][2], atoms[0][3])
    xyz.write(s)
    s="     2  HL1%12.6f%12.6f%12.6f   131     1\n"                    % (atoms[1][1], atoms[1][2], atoms[1][3])
    xyz.write(s)
    s="     3  HL2%12.6f%12.6f%12.6f   131     1\n"                    % (atoms[2][1], atoms[2][2], atoms[2][3])
    xyz.write(s)
    s="     4  HL3%12.6f%12.6f%12.6f   131     1\n"                    % (atoms[3][1], atoms[3][2], atoms[3][3])
    xyz.write(s)
    s="     5  C  %12.6f%12.6f%12.6f   124     1     6     7\n"        % (atoms[4][1], atoms[4][2], atoms[4][3])
    xyz.write(s)
    s="     6  O  %12.6f%12.6f%12.6f   125     5\n"                    % (atoms[5][1], atoms[5][2], atoms[5][3])
    xyz.write(s)
    s="     7  N  %12.6f%12.6f%12.6f   126     5     8     9\n"        % (atoms[6][1], atoms[6][2], atoms[6][3])
    xyz.write(s)
    s="     8  H  %12.6f%12.6f%12.6f   127     7\n"                    % (atoms[7][1], atoms[7][2], atoms[7][3])
    xyz.write(s)
    s="     9  CR %12.6f%12.6f%12.6f   128     7    10    11    12\n"  % (atoms[8][1], atoms[8][2], atoms[8][3])
    xyz.write(s)
    s="    10  HR1%12.6f%12.6f%12.6f   129     9\n"                    % (atoms[9][1], atoms[9][2], atoms[9][3])
    xyz.write(s)
    s="    11  HR2%12.6f%12.6f%12.6f   129     9\n"                    % (atoms[10][1], atoms[10][2], atoms[10][3])
    xyz.write(s)
    s="    12  HR3%12.6f%12.6f%12.6f   129     9\n"                    % (atoms[11][1], atoms[11][2], atoms[11][3])
    xyz.write(s)
    s="    13  OH2%12.6f%12.6f%12.6f    36    14    15\n"                          % (atoms[12][1], atoms[12][2], atoms[12][3])
    xyz.write(s)
    s="    14  H1 %12.6f%12.6f%12.6f    37    13\n"                          % (atoms[13][1], atoms[13][2], atoms[13][3])
    xyz.write(s)
    s="    15  H2 %12.6f%12.6f%12.6f    37    13\n"                          % (atoms[14][1], atoms[14][2], atoms[14][3])
    xyz.write(s)    

def write_nma_pot(atoms,crdname):
# note that we are using Gaussian to compute this
    qchem = open(crdname+".com",'w')
    subq = open("run."+crdname+".sh",'w')
    xyz = open(crdname+".xyz",'w')

    # write gaussian input files:
    qchem.write("%mem=300MW\n%NProc=2\n#MP2/6-311++G** Counterpoise=2\n\nNMA K+\n\n1,1 0,1 1,1\n")
    for i in range(12):
        s='%s 0%12.5f%12.5f%12.5f 1\n' % (atoms[i][0],atoms[i][1], atoms[i][2], atoms[i][3])
        qchem.write(s)
    s='K 0%12.5f%12.5f%12.5f 2\n\n' % (atoms[12][1], atoms[12][2], atoms[12][3])
    qchem.write(s)

    # write gaussian submission files on kickout:
    subq.write("#$ -S /bin/bash\n")
    subq.write("#$ -cwd\n")
    subq.write("#$ -V\n")
    subq.write("#$ -e /home/huangj/jobreport\n")
    subq.write("#$ -o /home/huangj/jobreport\n")
    subq.write("#$ -j y\n")
    s="#$ -N "+crdname+"\n"
    subq.write(s)
    subq.write("#$ -pe smp 2\n")
    subq.write("#$ -l h_data=800M,h_rt=300:00:00\n")
    subq.write("#$ -R y\n\n")

    subq.write("ulimit -c 4\n")
    subq.write("G03ROOT=/opt/mackerell/apps/gaussian/g03\n")
    subq.write("GAUSS_EXEDIR=$G03ROOT\n")
    subq.write("GAUSS_ARCHDIR=$HOME\n")
    subq.write("GAUSS_SCRDIR=/tmp/$USER/gauss$$\n")
    subq.write("GMAIN=$GAUSS_EXEDIR\n")
    subq.write("LD_LIBRARY_PATH=$GAUSS_EXEDIR:/opt/pgi/2012/lib\n")
    subq.write("PATH=$GAUSS_EXEDIR:$PATH\n")
    subq.write("export G03ROOT GAUSS_EXEDIR GAUSS_ARCHDIR GAUSS_SCRDIR PATH GMAIN LD_LIBRARY_PATH\n")
    subq.write("cwd=`pwd`\n")
    subq.write("if [ ! -d $GAUSS_SCRDIR ] ; then mkdir -p $GAUSS_SCRDIR; fi\n")
    subq.write("cd $GAUSS_SCRDIR\n")
    s="g03 < $cwd/"+crdname+".com > $cwd/"+crdname+".log\n"
    subq.write(s)
    subq.write("rm -f *.rwf\n")


    #write tinker xyz files
    xyz.write("    13\n")
    s="     1  CL %12.6f%12.6f%12.6f   130     2     3     4     5\n"  % (atoms[0][1], atoms[0][2], atoms[0][3])
    xyz.write(s)
    s="     2  HL1%12.6f%12.6f%12.6f   131     1\n"                    % (atoms[1][1], atoms[1][2], atoms[1][3])
    xyz.write(s)
    s="     3  HL2%12.6f%12.6f%12.6f   131     1\n"                    % (atoms[2][1], atoms[2][2], atoms[2][3])
    xyz.write(s)
    s="     4  HL3%12.6f%12.6f%12.6f   131     1\n"                    % (atoms[3][1], atoms[3][2], atoms[3][3])
    xyz.write(s)
    s="     5  C  %12.6f%12.6f%12.6f   124     1     6     7\n"        % (atoms[4][1], atoms[4][2], atoms[4][3])
    xyz.write(s)
    s="     6  O  %12.6f%12.6f%12.6f   125     5\n"                    % (atoms[5][1], atoms[5][2], atoms[5][3])
    xyz.write(s)
    s="     7  N  %12.6f%12.6f%12.6f   126     5     8     9\n"        % (atoms[6][1], atoms[6][2], atoms[6][3])
    xyz.write(s)
    s="     8  H  %12.6f%12.6f%12.6f   127     7\n"                    % (atoms[7][1], atoms[7][2], atoms[7][3])
    xyz.write(s)
    s="     9  CR %12.6f%12.6f%12.6f   128     7    10    11    12\n"  % (atoms[8][1], atoms[8][2], atoms[8][3])
    xyz.write(s)
    s="    10  HR1%12.6f%12.6f%12.6f   129     9\n"                    % (atoms[9][1], atoms[9][2], atoms[9][3])
    xyz.write(s)
    s="    11  HR2%12.6f%12.6f%12.6f   129     9\n"                    % (atoms[10][1], atoms[10][2], atoms[10][3])
    xyz.write(s)
    s="    12  HR3%12.6f%12.6f%12.6f   129     9\n"                    % (atoms[11][1], atoms[11][2], atoms[11][3])
    xyz.write(s)
    s="    13  POT%12.6f%12.6f%12.6f     8\n"                          % (atoms[12][1], atoms[12][2], atoms[12][3])
    xyz.write(s)
    

def write_nma_sod(atoms,crdname):
    qchem = open(crdname+".inp",'w')
    subq = open("run."+crdname+".sh",'w')
    xyz = open(crdname+".xyz",'w')

    # write qchem input files:
    qchem.write("$molecule\n1 1\n--\n0 1\n")
    for i in range(12):
        s=' %s%12.5f%12.5f%12.5f\n' % (atoms[i][0],atoms[i][1], atoms[i][2], atoms[i][3])
        qchem.write(s)
    qchem.write("--\n1 1\n")
    s='Na%12.5f%12.5f%12.5f\n' % (atoms[12][1], atoms[12][2], atoms[12][3])
    qchem.write(s)
    qchem.write("$end\n")
    qchem.write("\n")
    qchem.write("$rem\n")
    qchem.write("JOBTYPE                BSSE\n")
    qchem.write("EXCHANGE               HF\n")
    qchem.write("CORRELATION            RIMP2\n")
    qchem.write("BASIS                  CC-PVQZ\n")
    qchem.write("AUX_BASIS              RIMP2-CC-PVQZ\n")
    qchem.write("TIME_MP2               TRUE\n")
    qchem.write("PURECART               1111\n")
    qchem.write("N_FROZEN_CORE          FC\n")
    qchem.write("MEM_STATIC             1000\n")
    qchem.write("SCF_CONVERGENCE        8\n")
    qchem.write("THRESH                 12\n")
    qchem.write("symmetry = false\n")
    qchem.write("max_scf_cycles = 50\n")
    qchem.write("xc_grid = 000075000302\n")
    qchem.write("$end\n")

    # write qchem submission files on kickout:
    subq.write("#$ -S /bin/bash\n")
    subq.write("#$ -cwd\n")
    subq.write("#$ -V\n")
    subq.write("#$ -e /home/huangj/jobreport\n")
    subq.write("#$ -o /home/huangj/jobreport\n")
    subq.write("#$ -j y\n")
    s="#$ -N "+crdname+"\n"
    subq.write(s)
    subq.write("#$ -pe mpirun1 1\n")
    subq.write("#$ -l h_data=1600M,h_rt=300:00:00\n")
    subq.write("#$ -R y\n")
    subq.write("ulimit -c 4\n\n")
    subq.write("QC=/opt/mackerell/apps/qchem/4.0.1\n")
    subq.write(". $QC/qcenv.sh\n")
    subq.write("QCSCRATCH=`pwd`\n")
    subq.write("QCLOCALSCR=/tmp/$USER/tmp$$\n")
    subq.write("if [ ! -d $QCLOCALSCR ] ; then mkdir -p $QCLOCALSCR; fi\n")
    subq.write("export QC QCSCRATCH QCLOCALSCR\n\n")
    s="$QC/bin/qchem "+crdname+".inp "+crdname+".out\n"
    subq.write(s)
    
    #write tinker xyz files
    xyz.write("    13\n")
    s="     1  CL %12.6f%12.6f%12.6f   130     2     3     4     5\n"  % (atoms[0][1], atoms[0][2], atoms[0][3])
    xyz.write(s)
    s="     2  HL1%12.6f%12.6f%12.6f   131     1\n"                    % (atoms[1][1], atoms[1][2], atoms[1][3])
    xyz.write(s)
    s="     3  HL2%12.6f%12.6f%12.6f   131     1\n"                    % (atoms[2][1], atoms[2][2], atoms[2][3])
    xyz.write(s)
    s="     4  HL3%12.6f%12.6f%12.6f   131     1\n"                    % (atoms[3][1], atoms[3][2], atoms[3][3])
    xyz.write(s)
    s="     5  C  %12.6f%12.6f%12.6f   124     1     6     7\n"        % (atoms[4][1], atoms[4][2], atoms[4][3])
    xyz.write(s)
    s="     6  O  %12.6f%12.6f%12.6f   125     5\n"                    % (atoms[5][1], atoms[5][2], atoms[5][3])
    xyz.write(s)
    s="     7  N  %12.6f%12.6f%12.6f   126     5     8     9\n"        % (atoms[6][1], atoms[6][2], atoms[6][3])
    xyz.write(s)
    s="     8  H  %12.6f%12.6f%12.6f   127     7\n"                    % (atoms[7][1], atoms[7][2], atoms[7][3])
    xyz.write(s)
    s="     9  CR %12.6f%12.6f%12.6f   128     7    10    11    12\n"  % (atoms[8][1], atoms[8][2], atoms[8][3])
    xyz.write(s)
    s="    10  HR1%12.6f%12.6f%12.6f   129     9\n"                    % (atoms[9][1], atoms[9][2], atoms[9][3])
    xyz.write(s)
    s="    11  HR2%12.6f%12.6f%12.6f   129     9\n"                    % (atoms[10][1], atoms[10][2], atoms[10][3])
    xyz.write(s)
    s="    12  HR3%12.6f%12.6f%12.6f   129     9\n"                    % (atoms[11][1], atoms[11][2], atoms[11][3])
    xyz.write(s)
    s="    13  SOD%12.6f%12.6f%12.6f     7\n"                          % (atoms[12][1], atoms[12][2], atoms[12][3])
    xyz.write(s)
    



if __name__ == "__main__":
    resi="nma"
    resi2="pot"
    
    r1=2.6  
    convert(resi,resi2,r1,180,180)
    for a1 in range(90,180,15):
        #for d1 in range(180,-180,-15):
        for d1 in range(180,-15,-15): # symmetry
            convert(resi,resi2,r1,a1,d1)

#    for a1 in range(90,180,15):
#        for d1 in [90,-90]:
#            convert(resi,resi2,r1,a1,d1)



    
    
    
