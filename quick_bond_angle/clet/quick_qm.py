#!/usr/bin/python
import sys, re, os

#####################################################################
# Author: Xibing He                                                 #
# cretated :  12/20/2010                                            #
# latest modified: 02/14/2012                                       #
# E-mail: xibing@outerbanks.umaryland.edu                           #
# ###################################################################

# ###################################################################
# Purpose: extract the values of bonds, angles, dihedrals from the  #
#          gaussian output file <RESI_opt_freq_mp2.log>, extract    #
#          the atoms names from the file <RESI_water_hf.xyz>, and   #
#          put them together into file <quick_qm_mm.str>            #
# ###################################################################

#print sys.argv
if len(sys.argv) < 2:
   print "Usage: ./quick_qm_mm.py <residue>"
   print "Example: ./quick_qm_mm.py meso"
   print "Purpose: extract the values of bonds, angles, dihedrals from the gaussian output file,"
   print "         and generate the CHARMM input for quick measurements"
   print "Caution: The working directory should have the following two files exist:"
   print "              '<residue>_water_mp2.xyz' and '<residue>_opt_freq_mp2.log'"
   sys.exit()
else:
   residue = sys.argv[1]
   print "Residue name: %s"%residue

file1 = residue+"_water_mp2.xyz"
if os.path.exists(file1):
   print "Open file %s"%file1
   inp1 = open (file1, "r")
else:
   print "Error: file %s does not exist."%file1
   print "Is the residue name correct?? Lower case or capital case??"
   sys.exit(-1)

lines1 = inp1.readlines()
names = []
atom = 0
for line in lines1:
   words = line.split()
   atom = atom + 1
   if (atom == len(names) + 1):
      names.append(words[0])
   else:
      print "ERROR, the sequence number (%d) of new atom (%s) is NOT equal to len(names) + 1."%(atom, words[3])
      sys.exit(-1)
print names
inp1.close()

file2 = residue+"_opt_freq_mp2.log"
if os.path.exists(file2):
   print "Open file %s"%file2
   inp2 = open (file2, "r")
else:
   print "Error: file %s does not exist."%file1
   sys.exit(-1)

txt = inp2.read()
begin = "                           !   Optimized Parameters   !\n                           ! (Angstroms and Degrees)  !\n --------------------------                            --------------------------\n ! Name  Definition              Value          Derivative Info.                !"
end = """ --------------------------------------------------------------------------------
 GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad"""
#end = "Normal termination of Gaussian"
s1 = txt.find(begin, 0)
print "Find begin string at %d"%s1
s2 = txt.find(end,0)
print "Find end string at %d"%s2

file3 = "quick_qm_"+residue+".str"

out = open (file3, "w")
out.write("* quick bond and angle\n")
out.write("* \n")
out.write("!Bonds\n")

pattern = r" ! (R.*)-DE/DX =    0.0                 !\n"
result2 = re.findall(pattern, txt[s1:s2])
for line in result2:
   words = line.split()
   numbers = re.match(r"R\((\d+),(\d+)\)", words[1])
   atom1 = int(numbers.group(1))
   atom2 = int(numbers.group(2))
   value = float(words[2])
   out.write("!%s-%s %6.4f\n"%(names[atom1-1],names[atom2-1],value))
out.write("\n")

out.write("!Angles\n")
pattern = r" ! (A.*)-DE/DX =    0.0                 !\n"
result2 = re.findall(pattern, txt[s1:s2])
for line in result2:
   words = line.split()
   if words[1][0:2] == "A(":
      numbers = re.match(r"A\((\d+),(\d+),(\d+)\)", words[1])
      atom1 = int(numbers.group(1))
      atom2 = int(numbers.group(2))
      atom3 = int(numbers.group(3))
      value = float(words[2])
   elif words[1][0:2] == "L(":
      numbers = re.match(r"L\((\d+),(\d+),(\d+),(\d+),-(\d+)\)", words[1])
      atom1 = int(numbers.group(1))
      atom2 = int(numbers.group(2))
      atom3 = int(numbers.group(3))
      atom4 = int(numbers.group(4))
      atom5 = int(numbers.group(5))
      if atom5 == 1:
         value = 360.0 - float(words[2])
      else:
         continue
   else:
      print "***** The following line is a new situation. Contact author Xibing He."
      print line
      raise SystemExit
   out.write("!%s-%s-%s %8.4f\n"%(names[atom1-1],names[atom2-1],names[atom3-1],value))
out.write("\n")

#print "Info of bonds, angles, and dihedrals are extracted from <%s> and put into <%s>"%(file2, file3)
#print "Next: You should copy the content of <%s> to <%s_charmm.inp>"%(file3, residue) 

inp2.close()
out.close()
