#!/usr/bin/python
import sys, re, os

#####################################################################
# Author: Xibing He                                                 #
# Created: 09/30/2010                                               #
# Latest modified: 02/09/2012                                       #
# E-mail: xibing@outerbanks.umaryland.edu                           #
# ###################################################################

# ###################################################################
# Purpose: read the topology of a molecule, generate the CHARMM     # 
#          input file for the calculation of vibration modes        #
# Principle: This is based on the Table III & Figure 1 in the paper #
#          JACS 1979, 101 (10), 2550-2560, by P. Pulay et al.       #
# ###################################################################

# ###################################################################
# Users please cite the following paper:                            #
# Xibing He, Pedro E. M. Lopes and Alexander D. MacKerell, Jr.      #
# Polarizable Empirical Force Field for Acyclic Poly-Alcohols Based #
# on the Classical Drude Oscillator.                                #
# Submitted to Biopolymers, 2013                                    #
# ###################################################################


###
def Print_Caution():
   print "****Caution: If the molecule has ring structure, for each ring you need to add a line like \"!RING planar 5 C1 C2 C3 C4 C5\" in topology file. The word after 'RING' (eg. planar) doesn't matter to this script, but it matters to other scripts which you may use."

max_length_user_name = 12    ## max length of a username is 12 characters

###
def Print_Usage():
   print "****Usage: ./molvib_mm.py <residue_name> <topology_keyword> [auxiliary_ arguments]"
   print "   <residue_name>: eg. mesn, etoh (usually lower cases)"
   print "   <topology_keyword>: if this keyword is longer than %d letters, it will be used as the name (and path) of the topology file."%max_length_user_name
   print "                       if this keyword is \"full\", the topology file would be \"/raid/kenno/param/general/toppar/top_cgenff_all_2a.inp\""
   print "                       otherwise, this keyword is considered as a user name, and the topology file would be \"/raid/<user>/param/general/toppar/toppar_<user>_2a.str\""
   print "   <auxiliary_ arguments>: optional arguments for special cases, so that you can manually define the symmetry type. For instance, an atom is connected to several atoms, among which none is hydrogen. See examples below."
   print "   Example: ./molvib_mm.py psna xibing"
   print "   Example: ./molvib_mm.py mesn the_path/your_topology_file_name"
   print "   Example: ./molvib_mm.py pbsm full \"sp3 C:S1 H1:O11 H2:O12 X:CZ Y:N2\""
   Print_Caution()
   sys.exit()
###

# Check the arguments
if len(sys.argv) < 3:
   Print_Usage()
else:
   residue = sys.argv[1]
   print "Residue name: %s"%residue
   if len(sys.argv[2]) > max_length_user_name:
      print "****Caution: the 2nd argument is longer than %d characters, therefore the script consider this as the name of the topology file."%max_length_user_name
      topology_file = sys.argv[2]
   elif sys.argv[2] == "full":
      #user_name = sys.argv[2]
      #print "User name: %s"%user_name
      topology_file="/raid/kenno/param/general/toppar/top_cgenff_all_2a.inp"
   else:
      user_name = sys.argv[2]
      print "User name: %s"%user_name
      topology_file="/raid/%s/param/general/toppar/toppar_%s_2a.str"%(user_name, user_name)
   print "Topology file: %s"%topology_file


### Initialize lists and indices
ic_index = 0;       ic_list = [];      ic_bfly = []
umat_index = 0;     umat_list = []
ped_index = 0;      ped_list = []
dic_ring_flag = {}; list_rings = []
ring_IC = [];       ring_UMAT = [];    ring_PED = []

# Open the topology file
if os.path.exists(topology_file):
   print "Open file %s"%topology_file
   inp = open (topology_file, "r")
else:
   print "****Error: file %s does not exist."%topology_file
   sys.exit(-1)

txt = inp.read()
inp.close()
begin = "\nRESI %s "%residue.upper() ## First string of the residue definition
end = "\nRESI"  #  ## First string of the next residue definition
s1 = txt.find(begin, 0)              ## position of the beginning string
if (s1 > 0 ):
  print "Find beginning string of the residue definition at %d"%s1
else:
  print "****Error: cannot find '\\nRESI %s' in %s"%(residue.upper(), topology_file)
  sys.exit(-1)
s2 = txt.find(end,s1+4)                ## position of the ending string
if (s2 > s1 ):
  print "Find ending string of the residue definition at %d"%s2
else:
  print "Cannot find string '\\nRESI' after '\\nRESI %s', which means '\\nRESI %s' is the last residue definition in topology file."%(residue.upper(), residue.upper())
  end = "\nEND"  # if it is the last residue in topology file
  s2 = txt.find(end,s1+4)              ## position of the ending string
  if (s2 > s1 ):
     print "Find ending string of the residue definition at %d"%s2
  else:
     print "****Error: in topology file, %s is the last residue defined, but the script cannot find the appropriate ending line '\\nEND'."%residue.upper()
     sys.exit(-1)

# Search the definitions of atoms
print "Read in the lines of \"\\nATOMS \" of this residue from file %s:"%topology_file
pattern = "\n(ATOM .*)"
result = re.findall(pattern, txt[s1:s2])
# Get all of the atom names
list_atoms = []  ## all of the atoms in this residue
dic_atoms = {}   ## { atom:i} for all atoms
for line in result:
   words = line.split()
   list_atoms.append(words[1])
n_atoms = len(list_atoms) # total number of atoms in this residue
# Check if there are two atoms with the same name
for i in range(n_atoms-1):
   for j in range(i+1,n_atoms):
      if list_atoms[i] == list_atoms[j]:
	 print "****Error: two atoms have the same name: %s. Check file %s"%(list_atoms[i], topology_file)
# Create a dictionary for pairs: {atom_name : atom_index}
for i in range(n_atoms):
   atom = list_atoms[i]
   dic_atoms[atom] = i+1  ### i starts from 0, atom sequence starts from 1
print "Totally %d atoms in molecule %s:"%(n_atoms, residue.upper())
# Print the atom_names and corresponding atom_index
from operator import itemgetter                    # debug
print sorted(dic_atoms.items(), key=itemgetter(1)) # debug

# Search the definition of the bonds
pattern = "\n([BD]O[NU][DB].*)" # findall ^BOND or ^DOUB(LE)
result = re.findall(pattern, txt[s1:s2])
connections = {}
for atom in list_atoms:
   connections[atom] = []
for line in result:
   words = line.split()
   for i in range(1, len(words)-1, 2):
      connections[words[i]].append(words[i+1])
      connections[words[i+1]].append(words[i])

# Check if there are lines about rings
pattern = "(!RING .*)"
result = re.findall(pattern, txt[s1:s2])
if len(result) == 0:
   print "*****CAUTION: No line starting with '!RING ...' is found. The script think this residue has no ring structure."
   print "*****If this residue does have ring structure, please add a line like '!RING p 6 C1 C2 C3 C4 C5 C6' for each ring in the topology file."
else:
   print "This residue has %d ring(s):"%len(result)
   ### Get the ring atoms
   for line in result:   ## eg. "!RING p 5 C6 C7 N8 C9 C14"
      list_rings.append([])
      print line[:] ## debug
      words = line.split()
      try:
         natoms_ring = int(words[2])
      except ValueError:
         print "*****Error: the string '%s' in above line is not an integer, which is supposed to show the number of atoms in one ring."%words[2]
         print "*****Correct example: '!RING * 5 C1 C2 C3 C4 N5'"
         print "*****Check your topology file."
         sys.exit(-1)
      list_rings[-1].append(natoms_ring)  # type of the ring, eg. 5, or 6
      if len(words[3:]) != natoms_ring:
         print "*****Error: the total number of atoms (%d) in above line is not equal to the given value (%d)."%(len(words[3:]), natoms_ring)
         print "*****Check your topology file."
         sys.exit(-1)
      for atom in words[3:]:
         dic_ring_flag[atom] = 1
         list_rings[-1].append(atom)
   #print "list_rings", list_rings  ## debug
   ### Generate the <residue>.t2 file
   t2 = open("%s.t2"%residue, "w")
   for ring in list_rings:
      t2.write("RING %d\n"%ring[0])
      for atom in ring[1:]:
         #pass
         t2.write("%s %d"%(atom, dic_atoms[atom]))
         for substituent in connections[atom]:
            if not substituent in dic_ring_flag:
               t2.write(" %s %d"%(substituent, dic_atoms[substituent]))
         t2.write("\n")
   t2.close()
   print "File %s.t2 has been generated."%residue
   ### Generate the <residue>.rso file
   command = ""
   if os.path.exists("./ringsys"):
      command = "./ringsys < %s.t2 > %s.rso"%(residue, residue)
   elif os.path.exists("/raid/kenno/devl/ringsys"):
      command = "/raid/kenno/devl/ringsys < %s.t2 > %s.rso"%(residue, residue)
   else:
      if len(sys.argv) > 3:
         for item in sys.argv[3:]:
            #print item ## debug
            words = item.split()
            if len(words) < 2: continue
            if words[0] == "ringsys":
               if os.path.exists(words[1]):
                  command = "%s < %s.t2 > %s.rso"%(words[1], residue, residue)
      else:
         print "*****Error: Cannot find the script 'ringsys', please copy this script to current directory,"
         print "            or use an argument to show the path."
         print "     Example: ./molvib_mm.py meso xibing \"ringsys /raid/kenno/devl/ringsys\""
         sys.exit()
   if len(command) > 2:
      print "Try to generate file %s.rso"%residue
      print "     %s"%command
      os.system(command)
   ### Read the <residue>.rso file
   rso_file = "%s.rso"%residue
   if os.path.exists(rso_file):
      print "Open file %s"%rso_file
      inp1 = open(rso_file, "r")
   else:
      print "****Error: file %s does not exist."%rso_file
      print "Because there is (are) ring structure(s), you need to generate a %s.t2 file like:"%residue
      print "        /raid/kenno/param/general/models/b57a/b57a.t2"
      print "        /raid/kenno/param/general/models/b57b/b57b.t2"
      print "        /raid/kenno/param/general/models/3hin/3hin.t2"
      print "And then run the following command:"
      print "        /raid/kenno/devl/ringsys < %s.t2 > %s.rso"%(residue, residue)
      sys.exit(-1)
   info = inp1.read()
   inp1.close()
   s_IC = info.find("\nIC", 0)
   s_UMAT = info.find("\nUMAT", 0)
   s_PED = info.find("\nPED", 0)
   s_END = info.find("\nEND", 0)
   if (s_IC < 0 or s_UMAT < 0 or s_PED < 0 or s_END < 0):
      print "****Error: cannot find one of the following in %s file: '^IC', '^UMAT', '^PED', '^END'"%rso_file 
      sys.exit(-1)
   ring_IC = info[s_IC+1:s_UMAT].split("\n")[1:]
   ring_UMAT = info[s_UMAT+1:s_PED].split("\n")[1:-1]
   ring_PED = info[s_PED+1:s_END].split("\n")[1:-1]
   ic_index = len(ring_IC)
   umat_index = int(ring_UMAT[-1].split()[-3])
   ped_index = int(ring_PED[-1].split()[-2])
   print "IC, UMAT & PED info of ring(s) has been read from file %s."%rso_file
   #print "After reading %s, ic_index = %d, umat_index = %d."%(rso_file, ic_index, umat_index) ## debug

### If two rings share one edge (two atoms), add the butterfly torsions to IC, UMAT and PED
#print "the length of list_rings is: %d: "%len(list_rings)
if len(list_rings) > 1:  # more than one ring
   print "Check if two rings share one edge(two atoms)"
   for ring1 in list_rings[:-1]:
      #print ring1 ## debug
      for ring2 in list_rings[1:]:
         #print ring2 ## debug
         if ring1 == ring2: continue
         edge = [x for x in ring1[1:] if x in ring2[1:]]  ## atoms belong to both ring1 and ring2
         print "edge is:"; print edge ## debug
         if len(edge) > 2:
            print "****Error: two rings share more than 2 atoms. New situation. Contact author Xibing He (xhe@rx.umaryland.edu)."
            sys.exit(-1)
         elif len(edge) == 2:
            umat_index = umat_index + 1
            i_butterfly_count = 0
            for atom1 in connections[edge[0]]:
               if atom1 != edge[1]:
                  for atom4 in connections[edge[1]]:
                     if atom4 != edge[0]:
                        same_ring = 0
                        for ring in list_rings:
                           if (atom1 in ring) and (atom4 in ring): same_ring = 1
                        if same_ring == 0:
                           ic_index = ic_index + 1
                           ic_bfly.append(" 4%3d%3d%3d%3d  ! %d bfly %s-%s-%s-%s\n"%(dic_atoms[atom1], dic_atoms[edge[0]], dic_atoms[edge[1]], dic_atoms[atom4], ic_index, atom1, edge[0], edge[1], atom4))
                           ##umat_list.append("%3d %2d  1."%(umat_index, ic_index))
                           i_butterfly_count += 1
                           if i_butterfly_count % 2 == 1: umat_list.append("%3d %2d  1."%(umat_index, ic_index))
                           else: umat_list.append("%3d %2d -1."%(umat_index, ic_index))
	                   umat_list[-1] = umat_list[-1] + "        "
         #else:  ## len(edge) < 2
   #if (umat_index%4 != 0): umat_list[-1] = umat_list[-1] + "\n"
   umat_list[-1] = umat_list[-1] + "\n"
   ped_index = ped_index + 1
   ped_list.append("%3d  %-10s"%(umat_index, "BFLY"))
#for item in ic_bfly: ## debug
#   print item   ## debug

# Generate the items about bonds in the IC table
# Generate the items about bonds in the UMAT table
list_names1 = []  ## atoms of rings
list_names2 = []  ## atoms directly connected to rings
list_names3 = []  ## atoms in the chain parts 
dic_names1 = {}
dic_names2 = {}
dic_names3 = {}
for i in range(n_atoms):
   if (atom in dic_ring_flag): list_names1.append(atom)
print "Read in the lines of \"BOND \" and \"DOUBLE \" of this residue:"
pattern = "\n([BD]O[NU][DB].*)" # findall ^BOND or ^DOUB(LE)
result = re.findall(pattern, txt[s1:s2])
for line in result:
   print line ## debug
   words = line.split()
   for i in range(1, len(words)-1, 2):
      k1 = 0; k2 = 0  ## the atom sequence numbers of atom1 and atom2 in a bond
      for j in range(len(list_atoms)):
         if (list_atoms[j] == words[i]): k1 = j+1 ## j starts from 0
         if (list_atoms[j] == words[i+1]): k2 = j+1 ## j starts from 0
      #print "words[i] = %s, n1 = %d, words[i+1] = %s, n2 = %d"%(words[i], k1, words[i+1], k2) # debug
      if (k1 > 0 and k2 > 0): ## find the two atoms of a bond in the previous atom list 
         if ((words[i] in dic_ring_flag) and (words[i+1] in dic_ring_flag)):
            dic_names1[words[i]] = dic_atoms[words[i]]
            dic_names1[words[i+1]] = dic_atoms[words[i+1]]
         elif not ((words[i] in dic_ring_flag) or (words[i+1] in dic_ring_flag)):
            dic_names3[words[i]] = dic_atoms[words[i]]
            dic_names3[words[i+1]] = dic_atoms[words[i+1]]
            ic_index = ic_index + 1
            if words[i][0] == 'H' and words[i+1][0] != 'H':
              ic_list.append(" 1%3d%3d  0  0  ! %d %s-%s\n"%(k2, k1, ic_index, words[i+1], words[i]))
            elif words[i][0] == 'H' and words[i+1][0] == 'H':
              print "*****Error: Find a bond between two atoms: %s and %s, both of them have names starting with 'H'."
              print "            If an atom is not hydrogen, don't name it like 'H***'"
              sys.exit()
            else:
              ic_list.append(" 1%3d%3d  0  0  ! %d %s-%s\n"%(k1, k2, ic_index, words[i], words[i+1]))

            umat_index = umat_index + 1
            umat_list.append("%3d %2d  1."%(umat_index, ic_index))
            if (umat_index%4 == 0):
               umat_list[-1] = umat_list[-1] + "\n"
            else:
	       umat_list[-1] = umat_list[-1] + "        "
         else:
            if (words[i] in dic_ring_flag): dic_names1[words[i]] = dic_atoms[words[i]]
            else: dic_names2[words[i]] = dic_atoms[words[i]]
            if (words[i+1] in dic_ring_flag): dic_names1[words[i+1]] = dic_atoms[words[i+1]]
            else: dic_names2[words[i+1]] = dic_atoms[words[i+1]]
      else:  ## not find
         print "****Error: atom %s or %s in BONDS list is not in the ATOM list. Check %s"%(words[i], words[i+1],topology_file)
         sys.exit(-1)
if (umat_index%4 != 0 and len(umat_list) > 0): umat_list[-1] = umat_list[-1] + "\n"

#print sorted(dic_atoms.items(), key=itemgetter(1)) # debug
#print sorted(dic_ring_flag.items(), key=itemgetter(1)) # debug
#print sorted(dic_names1.items(), key=itemgetter(1)) # debug
#print sorted(dic_names2.items(), key=itemgetter(1)) # debug
print "atoms belonging to chain parts:"
print sorted(dic_names3.items(), key=itemgetter(1)) # debug
for item in sorted(dic_names3.items(), key=itemgetter(1)):
   list_names3.append(item[0])
#print list_names3 ## debug
#for item in ic_list: # debug
#   print item        # debug
#for item in umat_list: # debug
#   print item        # debug
#for item in ped_list: # debug
#   print item        # debug
#sys.exit() # debug

# Generate the items about bonds in the PED table
### Define a function Get_PED_bonds
def Get_PED_bonds(ic_list):
   global ped_index, ped_list
   list_temp = []
   for item in ic_list:
      s = item.find("!")
      words = item[s+1:].split()
      index = words[0]
      atom1 = words[1].split("-")[0]
      atom2 = words[1].split("-")[1]
      list_temp.append([index, atom1, atom2])
   #print list_temp ## debug
   list_temp_bak = list_temp[:]
   #print list_temp_bak ## debug [['1', 'S', 'O1'], ['2', 'S', 'O2'], ['3', 'S', 'O3'], ... ]
   n = len(list_temp)
   list_similar = []
   while (n > 0):
      #print n ## debug
      first = list_temp[0] ## ['1', 'S', 'O1']
      to_del = [0]
      similar = [first[0]] ## ['1']
      for i in range(1,len(list_temp)):
	 item = list_temp[i] ## ['2', 'S', 'O2']
	 if first[1] == item[1] and first[2][0] == item[2][0]: ## 's' == 's' and 'O' == 'O'
	    similar.append(item[0]) ## ['1', '2']
	    to_del.append(i) ## [2]
      k = 0
      for j in to_del: # [2,3]
	 list_temp.pop(j-k)
	 k = k + 1
      list_similar.append(similar) ## [ ['1', '2', '3'], [...] ]
      n = len(list_temp)
   #print list_similar ## debug [['1', '2', '3'], ['4'], ['5', '6'], ['7'], ['8', '9'], ['10'], ['11', '12', '13']]
                      ## debug b57b [['63', 'C5', 'H51'], ['64', 'C5', 'H52']]
   for item in list_similar: ## ['1', '2', '3'], or b57b [['63', '64']]
      index = int(item[0])   ## 1 or 63
      for temp in list_temp_bak: ## ['1', 'S', 'O1'], or ['63', 'C5', 'H51']
         if index == int(temp[0]):
            atom1 = temp[1]
            atom2 = temp[2]
      #print "index: %d, atom1: %s, atom2: %s"%(index, atom1, atom2) ## debug #index: 1, atom1: S, atom2: O1
      if atom2[0] == "H": word = "s%s-%s"%(atom1, atom2[0])
      else: word = "s%s-%s"%(atom1, atom2)
      ##if atom2[0] == "H": word = "str%s-%s"%(atom1, atom2[0])
      ##else: word = "str%s-%s"%(atom1, atom2)
      ped_index = ped_index + 1
      for ttt in ic_list:
         #print ttt ## debug
         ppp = "  0  0  ! \d+ %s-%s"%(atom1, atom2)
         rrr = re.findall(ppp, ttt)
         if len(rrr) > 0:
            #print rrr
            iii = int(rrr[0].split()[-2])
      for ttt in umat_list:
         sss = ttt.split()
         for jjj in range(1, len(sss), 3):
            if int(sss[jjj]) == iii:
               kkk = int(sss[jjj-1])
      ped_list.append("%3d  %-10s"%(kkk, word))
      #print "ped_index = %d, word = %s"%(ped_index, word)  ## debug
   return list_similar
### End of definition of function Get_PED_bonds
list_similar = Get_PED_bonds(ic_list) ## debug

#for item in ic_list: # debug
#   print item        # debug
#for item in umat_list: # debug
#   print item        # debug
#for item in ped_list: # debug
#   print item        # debug
#sys.exit()  ## debug

###
def Print_Symmetry(center, others, symmetry_type):
   print "The symmetry type of atom %4s is %13s, which is connected to %d atom(s)"%(center, symmetry_type, len(others)) ## debug
###

### Define a function, which judge the symmetry type of a center atom
def Get_symmetry_type (center, others):
### center is the name of the center atom
### others is a list of names of atoms connecting to the center atom
   others.sort()
   list_repeat = []
   n_H = 0  # number of Hydrogen atoms in 'others'
   n_O = 0  # number of Oxygen atoms in 'others'
   n_C = 0  # number of Carbon atoms in 'others'
   for i in range(len(others)):
      list_repeat.append(0)
      if others[i][0] == "H": n_H = n_H + 1
      if others[i][0] == "O": n_O = n_O + 1
      if others[i][0] == "C": n_C = n_C + 1
      for j in range(len(others)):
	 if others[i][0] == others[j][0]:  ## The first letter of atom names are the same
            list_repeat[i] = list_repeat[i] + 1
   #print "Analyze the atom %4s, which is connected to %d atom(s)"%(center, len(others)) ## debug
   if len(others) > 4:
      symmetry_type = "unknown"
      Print_Symmetry(center, others, symmetry_type)
      return symmetry_type, []
   elif len(others) == 4:
      if (n_H == 3):
	 symmetry_type = "methyl_sp3"
	 for i in range(len(others)):
	    if others[i][0] == "H": index_H1 = i; atom_H1 = others[i]; break
	 for i in range(index_H1+1,len(others)):
	    if others[i][0] == "H": index_H2 = i; atom_H2 = others[i]; break
	 for i in range(index_H2+1,len(others)):
	    if others[i][0] == "H": index_H3 = i; atom_H3 = others[i]; break
         for i in range(len(others)):
	    if (i != index_H1 and i != index_H2 and i != index_H3): index_X = i; atom_X = others[i]; break
         Print_Symmetry(center, others, symmetry_type)
	 return symmetry_type, [atom_H1, atom_H2, atom_H3, atom_X];
      elif (n_H == 2):
	 symmetry_type = "methylene_sp3"
	 for i in range(len(others)):
	    if others[i][0] == "H": index_H1 = i; atom_H1 = others[i]; break
	 for i in range(index_H1+1,len(others)):
	    if others[i][0] == "H": index_H2 = i; atom_H2 = others[i]; break
         for i in range(len(others)):
	    if (i != index_H1 and i != index_H2): index_X = i; atom_X = others[i]; break
         for i in range(len(others)):
	    if (i != index_H1 and i != index_H2 and i != index_X): index_Y = i; atom_Y = others[i]; break
         Print_Symmetry(center, others, symmetry_type)
	 return symmetry_type, [atom_H1, atom_H2, atom_X, atom_Y]
      elif (n_H == 1):
	 symmetry_type = "methine_sp3"
	 for i in range(len(others)):
	    if others[i][0] == "H": index_H1 = i; atom_H1 = others[i]; break
         for i in range(len(others)):
	    if (i != index_H1): index_X = i; atom_X = others[i]; break
         for i in range(len(others)):
	    if (i != index_H1 and i != index_X): index_Y = i; atom_Y = others[i]; break
         for i in range(len(others)):
	    if (i != index_H1 and i != index_X and i != index_Y): index_Z = i; atom_Z = others[i]; break
         Print_Symmetry(center, others, symmetry_type)
	 return symmetry_type, [atom_H1, atom_X, atom_Y, atom_Z]
      else:
	 #print "Caution: no Hydrogen attached to center atom %s"%center
         symmetry_type = "sp3"
         if len(sys.argv) > 3:
	    match_flag = 0
	    for item in sys.argv[3:]:
		#print item ## debug
		words = item.split()
		if len(words) < 6: continue
		if words[0] == "sp3":
	           if words[1].split(":")[1] == center:
                      match_flage = 1
		      if words[5].split(":")[0] == "Z":
			 symmetry_type = "methine_sp3"
			 atom_H1 = words[2].split(":")[1]
			 atom_X = words[3].split(":")[1]
			 atom_Y = words[4].split(":")[1]
			 atom_Z = words[5].split(":")[1]
                         Print_Symmetry(center, others, symmetry_type)
			 return symmetry_type, [atom_H1, atom_X, atom_Y, atom_Z]
		      elif words[5].split(":")[0] == "Y":
			 symmetry_type = "methylene_sp3"
			 atom_H1 = words[2].split(":")[1]
			 atom_H2 = words[3].split(":")[1]
			 atom_X = words[4].split(":")[1]
			 atom_Y = words[5].split(":")[1]
                         Print_Symmetry(center, others, symmetry_type)
			 return symmetry_type, [atom_H1, atom_H2, atom_X, atom_Y]
		      else:
			 symmetry_type = "methyl_sp3"
			 atom_H1 = words[2].split(":")[1]
			 atom_H2 = words[3].split(":")[1]
			 atom_H3 = words[4].split(":")[1]
			 atom_X = words[5].split(":")[1]
                         Print_Symmetry(center, others, symmetry_type)
			 return symmetry_type, [atom_H1, atom_H2, atom_H3, atom_X]
	    if (match_flag == 0):
	       print "****Error: Need a correct argument for atom %s with symmetry type %s."%(center, symmetry_type)
	       Print_Usage()
         else:
	    print "****Error: atom %s is connected to 4 atoms, and none of them is hydrogen."%center
	    print "       You need to give one more argument like one of the following:"
            print "             \"sp3 C:<atom0> H1:<atom1> H2:<atom2> H3:<atom3> X:<atom4>\""
            print "             \"sp3 C:<atom0> H1:<atom1> H2:<atom2> X:<atom3> Y:<atom4>\""
            print "             \"sp3 C:<atom0> H1:<atom1> X:<atom2> Y:<atom3> Z:<atom4>\""
	    sys.exit(-1)
         Print_Symmetry(center, others, symmetry_type)
         return symmetry_type, []
   elif len(others) == 3:
      if (n_H == 2):
	 if center[0] == "C": symmetry_type = "methylene_sp2"
	 elif center[0] == "N": symmetry_type = "amino"
         else: symmetry_type = "methylene_sp2"  #######################!!!!!!!!!!!!!!!!!!!!!!
	 for i in range(len(others)):
	    if others[i][0] == "H": index_H1 = i; atom_H1 = others[i]; break
	 for i in range(index_H1+1,len(others)):
	    if others[i][0] == "H": index_H2 = i; atom_H2 = others[i]; break
         for i in range(len(others)):
	    if (i != index_H1 and i != index_H2): index_X = i; atom_X = others[i]; break
         Print_Symmetry(center, others, symmetry_type)
	 return symmetry_type, [atom_H1, atom_H2, atom_X]
      elif (n_H == 1):
	 if center[0] == "C": symmetry_type = "methine_sp2"
	 elif center[0] == "N": symmetry_type = "imino"
         else: symmetry_type = "methine_sp2"
	 for i in range(len(others)):
	    if others[i][0] == "H": index_H1 = i; atom_H1 = others[i]; break
         for i in range(len(others)):
	    if (i != index_H1): index_X = i; atom_X = others[i]; break
         for i in range(len(others)):
	    if (i != index_H1 and i != index_X): index_Y = i; atom_Y = others[i]; break
         Print_Symmetry(center, others, symmetry_type)
	 return symmetry_type, [atom_H1, atom_X, atom_Y]
      else:
	 #print "Caution: no Hydrogen attached to center atom %s"%center
         symmetry_type = "sp2"
         if len(sys.argv) > 3:
	    match_flag = 0
	    for item in sys.argv[3:]:
		#print item ## debug
		words = item.split()
		if len(words) < 5: continue
		if words[0] == "sp2":
	           if words[1].split(":")[1] == center:
                      match_flage = 1
		      if words[4].split(":")[0] == "Y":
			 symmetry_type = "methine_sp2"
			 atom_H1 = words[2].split(":")[1]
			 atom_X = words[3].split(":")[1]
			 atom_Y = words[4].split(":")[1]
                         Print_Symmetry(center, others, symmetry_type)
			 return symmetry_type, [atom_H1, atom_X, atom_Y]
		      else:
			 symmetry_type = "methylene_sp2"
			 atom_H1 = words[2].split(":")[1]
			 atom_H2 = words[3].split(":")[1]
			 atom_X = words[4].split(":")[1]
                         Print_Symmetry(center, others, symmetry_type)
			 return symmetry_type, [atom_H1, atom_H2, atom_X]
	    if (match_flag == 0):
	       print "****Error: Need a correct argument for atom %s with symmetry type %s."%(center, symmetry_type)
	       Print_Usage()
         else:
	    print "****Error: atom %s is connected to 3 atoms, and none of them is hydrogen."%center
	    print "       You need to give one more argument like \"sp2 C:<atom0> H1:<atom1> H2:<atom2> X:<atom3>\""
	    print "       or like \"sp2 C:<atom0> H1:<atom1> X:<atom2> Y:<atom3>\""
	    sys.exit(-1)
         Print_Symmetry(center, others, symmetry_type)
         return symmetry_type, []
   elif len(others) == 2:
      symmetry_type = "sp"
      Print_Symmetry(center, others, symmetry_type)
      return symmetry_type, others
   elif len(others) == 1:
      symmetry_type = "end" ##### The atom is only connected to one atom, therefore it is an end atom.
      Print_Symmetry(center, others, symmetry_type)
      return symmetry_type, others
   else:
      symmetry_type = "unknown"
      Print_Symmetry(center, others, symmetry_type)
      return symmetry_type, []
### End of definition of function Get_symmetry_type 

### Generate the items about angles in the IC table
### Generate the items about angles in the UMAT table
for atom in list_names3:
   symmetry_type, other_atoms = Get_symmetry_type(atom, connections[atom])
   ###
   if len(other_atoms) < 1:
      print "****Error: the symmetry type of atom %s is unknown or hard to determine. Contact author Xibing He (xhe@rx.umaryland.edu.)"%atom
   ###
   if symmetry_type == "unknown":
      print "****Error: the symmetry type of atom %s is unknown. Contact author Xibing He (xhe@rx.umaryland.edu.)"%atom
      sys.exit(-1)
   ###
   elif symmetry_type == "end":
      continue  #######################
   elif symmetry_type == "sp2":
      print "****Error: the symmetry type of atom %s is sp2. Contact author Xibing He (xhe@rx.umaryland.edu)."%atom
      sys.exit(-1)
   ###
   elif symmetry_type == "methyl_sp3":
      atom_H1 = other_atoms[0]; atom_H2 = other_atoms[1]; atom_H3 = other_atoms[2]; atom_X = other_atoms[3]
      ic_index = ic_index + 1
      ic_list.append(" 2%3d%3d%3d  0  ! %d %s-%s-%s %s\n"%(dic_atoms[atom_H2], dic_atoms[atom], dic_atoms[atom_H3], ic_index, atom_H2, atom, atom_H3, 'a1'))
      ic_index = ic_index + 1
      ic_list.append(" 2%3d%3d%3d  0  ! %d %s-%s-%s %s\n"%(dic_atoms[atom_H3], dic_atoms[atom], dic_atoms[atom_H1], ic_index, atom_H3, atom, atom_H1, 'a2'))
      ic_index = ic_index + 1
      ic_list.append(" 2%3d%3d%3d  0  ! %d %s-%s-%s %s\n"%(dic_atoms[atom_H1], dic_atoms[atom], dic_atoms[atom_H2], ic_index, atom_H1, atom, atom_H2, 'a3'))
      ic_index = ic_index + 1
      ic_list.append(" 2%3d%3d%3d  0  ! %d %s-%s-%s %s\n"%(dic_atoms[atom_X], dic_atoms[atom], dic_atoms[atom_H1], ic_index, atom_X, atom, atom_H1, 'b1'))
      ic_index = ic_index + 1
      ic_list.append(" 2%3d%3d%3d  0  ! %d %s-%s-%s %s\n"%(dic_atoms[atom_X], dic_atoms[atom], dic_atoms[atom_H2], ic_index, atom_X, atom, atom_H2, 'b2'))
      ic_index = ic_index + 1
      ic_list.append(" 2%3d%3d%3d  0  ! %d %s-%s-%s %s\n"%(dic_atoms[atom_X], dic_atoms[atom], dic_atoms[atom_H3], ic_index, atom_X, atom, atom_H3, 'b3'))
      umat_index = umat_index + 1
      umat_list.append("%3d %2d  1.        %3d %2d  1.        %3d %2d  1.        %3d %2d -1.\n%3d %2d -1.        %3d %2d -1.\n"%(umat_index, ic_index-5, umat_index, ic_index-4, umat_index, ic_index-3, umat_index, ic_index-2, umat_index, ic_index-1, umat_index, ic_index-0))
      ped_index = ped_index + 1
      word = "ds%s%s"%(atom, atom_H1[0])
      ##word = "sym%s-%s"%(atom, atom_H1[0])
      ped_list.append("%3d  %-10s"%(umat_index, word))
      umat_index = umat_index + 1
      umat_list.append("%3d %2d  2.        %3d %2d -1.        %3d %2d -1.\n"%(umat_index, ic_index-5, umat_index, ic_index-4, umat_index, ic_index-3))
      ped_index = ped_index + 1
      word = "da%s%s"%(atom, atom_H1[0])
      ##word = "asy%s-%s"%(atom, atom_H1[0])
      ped_list.append("%3d  %-10s"%(umat_index, word))
      umat_index = umat_index + 1
      umat_list.append("%3d %2d  1.        %3d %2d -1.\n"%(umat_index, ic_index-4, umat_index, ic_index-3))
      ped_index = ped_index + 1
      word = "da%s%s'"%(atom, atom_H1[0])
      ##word = "asp%s-%s"%(atom, atom_H1[0])
      ped_list.append("%3d  %-10s"%(umat_index, word))
      umat_index = umat_index + 1
      umat_list.append("%3d %2d  2.        %3d %2d -1.        %3d %2d -1.\n"%(umat_index, ic_index-2, umat_index, ic_index-1, umat_index, ic_index-0))
      ped_index = ped_index + 1
      word = "r%s%s"%(atom, atom_H1[0])
      ##word = "roc%s-%s"%(atom, atom_H1[0])
      ped_list.append("%3d  %-10s"%(umat_index, word))
      umat_index = umat_index + 1
      umat_list.append("%3d %2d  1.        %3d %2d -1.\n"%(umat_index, ic_index-1, umat_index, ic_index-0))
      ped_index = ped_index + 1
      word = "r%s%s'"%(atom, atom_H1[0])
      ##word = "rop%s-%s"%(atom, atom_H1[0])
      ped_list.append("%3d  %-10s"%(umat_index, word))
   ###
   elif symmetry_type == "methylene_sp3":
      atom_H1 = other_atoms[0]; atom_H2 = other_atoms[1]; atom_X = other_atoms[2]; atom_Y = other_atoms[3]
      ic_index = ic_index + 1
      ic_list.append(" 2%3d%3d%3d  0  ! %d %s-%s-%s %s\n"%(dic_atoms[atom_H1], dic_atoms[atom], dic_atoms[atom_H2], ic_index, atom_H1, atom, atom_H2, 'a'))
      ic_index = ic_index + 1
      ic_list.append(" 2%3d%3d%3d  0  ! %d %s-%s-%s %s\n"%(dic_atoms[atom_X], dic_atoms[atom], dic_atoms[atom_Y], ic_index, atom_X, atom, atom_Y, 'g'))
      ic_index = ic_index + 1
      ic_list.append(" 2%3d%3d%3d  0  ! %d %s-%s-%s %s\n"%(dic_atoms[atom_H1], dic_atoms[atom], dic_atoms[atom_Y], ic_index, atom_H1, atom, atom_Y, 'b1'))
      ic_index = ic_index + 1
      ic_list.append(" 2%3d%3d%3d  0  ! %d %s-%s-%s %s\n"%(dic_atoms[atom_H2], dic_atoms[atom], dic_atoms[atom_Y], ic_index, atom_H2, atom, atom_Y, 'b2'))
      ic_index = ic_index + 1
      ic_list.append(" 2%3d%3d%3d  0  ! %d %s-%s-%s %s\n"%(dic_atoms[atom_H1], dic_atoms[atom], dic_atoms[atom_X], ic_index, atom_H1, atom, atom_X, 'b3'))
      ic_index = ic_index + 1
      ic_list.append(" 2%3d%3d%3d  0  ! %d %s-%s-%s %s\n"%(dic_atoms[atom_H2], dic_atoms[atom], dic_atoms[atom_X], ic_index, atom_H2, atom, atom_X, 'b4'))
      umat_index = umat_index + 1
      umat_list.append("%3d %2d  5.        %3d %2d  1.\n"%(umat_index, ic_index-5, umat_index, ic_index-4))
      ped_index = ped_index + 1
      word = "c%s-%s"%(atom, atom_H1[0])
      ##word = "sci%s-%s"%(atom, atom_H1[0])
      ped_list.append("%3d  %-10s"%(umat_index, word))
      umat_index = umat_index + 1
      umat_list.append("%3d %2d  1.        %3d %2d  5.\n"%(umat_index, ic_index-5, umat_index, ic_index-4))
      ped_index = ped_index + 1
      word = "d%s%s%s"%(atom, atom_X[0], atom_Y[0])
      ##word = "sci%s-%s-%s"%(atom, atom_X[0], atom_Y[0])
      ped_list.append("%3d  %-10s"%(umat_index, word))
      umat_index = umat_index + 1
      umat_list.append("%3d %2d  1.        %3d %2d -1.        %3d %2d  1.        %3d %2d -1.\n"%(umat_index, ic_index-3, umat_index, ic_index-2, umat_index, ic_index-1, umat_index, ic_index-0))
      ped_index = ped_index + 1
      word = "r%s%s"%(atom, atom_H1[0])
      ##word = "roc%s-%s"%(atom, atom_H1[0])
      ped_list.append("%3d  %-10s"%(umat_index, word))
      umat_index = umat_index + 1
      umat_list.append("%3d %2d  1.        %3d %2d  1.        %3d %2d -1.        %3d %2d -1.\n"%(umat_index, ic_index-3, umat_index, ic_index-2, umat_index, ic_index-1, umat_index, ic_index-0))
      ped_index = ped_index + 1
      word = "w%s%s"%(atom, atom_H1[0])
      ##word = "wag%s-%s"%(atom, atom_H1[0])
      ped_list.append("%3d  %-10s"%(umat_index, word))
      umat_index = umat_index + 1
      umat_list.append("%3d %2d  1.        %3d %2d -1.        %3d %2d -1.        %3d %2d  1.\n"%(umat_index, ic_index-3, umat_index, ic_index-2, umat_index, ic_index-1, umat_index, ic_index-0))
      ped_index = ped_index + 1
      word = "i%s%s"%(atom, atom_H1[0])
      ##word = "twi%s-%s"%(atom, atom_H1[0])
      ped_list.append("%3d  %-10s"%(umat_index, word))
   ###
   elif symmetry_type == "methine_sp3":
      atom_H1 = other_atoms[0]; atom_X = other_atoms[1]; atom_Y = other_atoms[2]; atom_Z = other_atoms[3]
      ic_index = ic_index + 1
      ic_list.append(" 2%3d%3d%3d  0  ! %d %s-%s-%s %s\n"%(dic_atoms[atom_H1], dic_atoms[atom], dic_atoms[atom_X], ic_index, atom_H1, atom, atom_X, 'b1'))
      ic_index = ic_index + 1
      ic_list.append(" 2%3d%3d%3d  0  ! %d %s-%s-%s %s\n"%(dic_atoms[atom_H1], dic_atoms[atom], dic_atoms[atom_Y], ic_index, atom_H1, atom, atom_Y, 'b2'))
      ic_index = ic_index + 1
      ic_list.append(" 2%3d%3d%3d  0  ! %d %s-%s-%s %s\n"%(dic_atoms[atom_H1], dic_atoms[atom], dic_atoms[atom_Z], ic_index, atom_H1, atom, atom_Z, 'b3'))
      ic_index = ic_index + 1
      ic_list.append(" 2%3d%3d%3d  0  ! %d %s-%s-%s %s\n"%(dic_atoms[atom_X], dic_atoms[atom], dic_atoms[atom_Y], ic_index, atom_X, atom, atom_Y, 'axcy'))
      ic_index = ic_index + 1
      ic_list.append(" 2%3d%3d%3d  0  ! %d %s-%s-%s %s\n"%(dic_atoms[atom_Y], dic_atoms[atom], dic_atoms[atom_Z], ic_index, atom_Y, atom, atom_Z, 'aycz'))
      ic_index = ic_index + 1
      ic_list.append(" 2%3d%3d%3d  0  ! %d %s-%s-%s %s\n"%(dic_atoms[atom_Z], dic_atoms[atom], dic_atoms[atom_X], ic_index, atom_Z, atom, atom_X, 'azcx'))
      umat_index = umat_index + 1
      umat_list.append("%3d %2d  2.        %3d %2d -1.        %3d %2d -1.\n"%(umat_index, ic_index-5, umat_index, ic_index-4, umat_index, ic_index-3))
      ped_index = ped_index + 1
      word = "r%s%s"%(atom, atom_H1[0])
      ##word = "roc%s-%s"%(atom, atom_H1[0])
      ped_list.append("%3d  %-10s"%(umat_index, word))
      umat_index = umat_index + 1
      umat_list.append("%3d %2d  1.        %3d %2d -1.\n"%(umat_index, ic_index-4, umat_index, ic_index-3))
      ped_index = ped_index + 1
      word = "r%s%s'"%(atom, atom_H1[0])
      ##word = "rop%s-%s"%(atom, atom_H1[0])
      ped_list.append("%3d  %-10s"%(umat_index, word))
      umat_index = umat_index + 1
      umat_list.append("%3d %2d  4.        %3d %2d  1.        %3d %2d  1.\n"%(umat_index, ic_index-2, umat_index, ic_index-1, umat_index, ic_index-0))
      ped_index = ped_index + 1
      word = "d%s%s%s"%(atom_X, atom, atom_Y)
      ##word = "%s%s%s"%(atom_X, atom, atom_Y)
      ped_list.append("%3d  %-10s"%(umat_index, word))
      umat_index = umat_index + 1
      umat_list.append("%3d %2d  4.        %3d %2d  1.        %3d %2d  1.\n"%(umat_index, ic_index-1, umat_index, ic_index-0, umat_index, ic_index-2))
      ped_index = ped_index + 1
      word = "d%s%s%s"%(atom_Y, atom, atom_Z)
      ##word = "%s%s%s"%(atom_Y, atom, atom_Z)
      ped_list.append("%3d  %-10s"%(umat_index, word))
      umat_index = umat_index + 1
      umat_list.append("%3d %2d  4.        %3d %2d  1.        %3d %2d  1.\n"%(umat_index, ic_index-0, umat_index, ic_index-2, umat_index, ic_index-1))
      ped_index = ped_index + 1
      word = "d%s%s%s"%(atom_Z, atom, atom_X)
      ##word = "%s%s%s"%(atom_Z, atom, atom_X)
      ped_list.append("%3d  %-10s"%(umat_index, word))
   ###
   elif symmetry_type == "methylene_sp2" or symmetry_type == "amino":
      atom_H1 = other_atoms[0]; atom_H2 = other_atoms[1]; atom_X = other_atoms[2]
      ic_index = ic_index + 1
      ic_list.append(" 2%3d%3d%3d  0  ! %d %s-%s-%s %s\n"%(dic_atoms[atom_H1], dic_atoms[atom], dic_atoms[atom_H2], ic_index, atom_H1, atom, atom_H2, 'a'))
      ic_index = ic_index + 1
      ic_list.append(" 2%3d%3d%3d  0  ! %d %s-%s-%s %s\n"%(dic_atoms[atom_H1], dic_atoms[atom], dic_atoms[atom_X], ic_index, atom_H1, atom, atom_X, 'b1'))
      ic_index = ic_index + 1
      ic_list.append(" 2%3d%3d%3d  0  ! %d %s-%s-%s %s\n"%(dic_atoms[atom_H2], dic_atoms[atom], dic_atoms[atom_X], ic_index, atom_H2, atom, atom_X, 'b2'))
      ic_index = ic_index + 1
      ic_list.append(" 3%3d%3d%3d%3d  ! %d %s improper\n"%(dic_atoms[atom_H1], dic_atoms[atom_H2], dic_atoms[atom_X], dic_atoms[atom], ic_index, atom))
      umat_index = umat_index + 1
      umat_list.append("%3d %2d  2.        %3d %2d -1.        %3d %2d -1.\n"%(umat_index, ic_index-3, umat_index, ic_index-2, umat_index, ic_index-1))
      ped_index = ped_index + 1
      if symmetry_type == "methylene_sp2": word = "d%s%s"%(atom, atom_H1)
      ##if symmetry_type == "methylene_sp2": word = "sym%s-%s"%(atom, atom_H1[0])
      else: word = "d%s%s"%(atom, atom_H1)
      ##else: word = "sci%s-%s"%(atom, atom_H1[0])
      ped_list.append("%3d  %-10s"%(umat_index, word))
      umat_index = umat_index + 1
      umat_list.append("%3d %2d  1.        %3d %2d -1.\n"%(umat_index, ic_index-2, umat_index, ic_index-1))
      ped_index = ped_index + 1
      word = "d%s%s"%(atom, atom_X)  ### ???
      ##word = "roc%s-%s"%(atom, atom_H1[0])
      ped_list.append("%3d  %-10s"%(umat_index, word))
      umat_index = umat_index + 1
      umat_list.append("%3d %2d  1.\n"%(umat_index, ic_index-0))
      ped_index = ped_index + 1
      word = "w%s%s"%(atom, atom_X)  ### ???
      ##word = "wag%s"%(atom, )
      ped_list.append("%3d  %-10s"%(umat_index, word))
   ###
   elif symmetry_type == "methine_sp2" or symmetry_type == "imino":
      atom_H1 = other_atoms[0]; atom_X = other_atoms[1]; atom_Y = other_atoms[2]
      ic_index = ic_index + 1
      ic_list.append(" 2%3d%3d%3d  0  ! %d %s-%s-%s %s\n"%(dic_atoms[atom_X], dic_atoms[atom], dic_atoms[atom_Y], ic_index, atom_X, atom, atom_Y, 'a'))
      ic_index = ic_index + 1
      ic_list.append(" 2%3d%3d%3d  0  ! %d %s-%s-%s %s\n"%(dic_atoms[atom_X], dic_atoms[atom], dic_atoms[atom_H1], ic_index, atom_X, atom, atom_H1, 'b1'))
      ic_index = ic_index + 1
      ic_list.append(" 2%3d%3d%3d  0  ! %d %s-%s-%s %s\n"%(dic_atoms[atom_Y], dic_atoms[atom], dic_atoms[atom_H1], ic_index, atom_Y, atom, atom_H1, 'b2'))
      ic_index = ic_index + 1
      ic_list.append(" 3%3d%3d%3d%3d  ! %d %s improper\n"%(dic_atoms[atom_H1], dic_atoms[atom_X], dic_atoms[atom_Y], dic_atoms[atom], ic_index, atom))
      umat_index = umat_index + 1
      umat_list.append("%3d %2d  2.        %3d %2d -1.        %3d %2d -1.\n"%(umat_index, ic_index-3, umat_index, ic_index-2, umat_index, ic_index-1))
      ped_index = ped_index + 1
      word = "d%s%s%s"%(atom_X, atom, atom_Y)
      ##word = "%s%s%s"%(atom_X, atom, atom_Y)
      ped_list.append("%3d  %-10s"%(umat_index, word))
      umat_index = umat_index + 1
      umat_list.append("%3d %2d  1.        %3d %2d -1.\n"%(umat_index, ic_index-2, umat_index, ic_index-1))
      ped_index = ped_index + 1
      word = "d%s%s"%(atom, atom_H1[0])
      ##word = "roc%s-%s"%(atom, atom_H1[0])
      ped_list.append("%3d  %-10s"%(umat_index, word))
      umat_index = umat_index + 1
      umat_list.append("%3d %2d  1.\n"%(umat_index, ic_index-0))
      ped_index = ped_index + 1
      word = "w%s%s"%(atom, atom_H1)  ### ???
      ##word = "wag%s"%(atom, )
      ped_list.append("%3d  %-10s"%(umat_index, word))
   elif symmetry_type == "sp":
      atom_X = other_atoms[0]; atom_Y = other_atoms[1]
      ic_index = ic_index + 1
      ic_list.append(" 2%3d%3d%3d  0  ! %d %s-%s-%s\n"%(dic_atoms[atom_X], dic_atoms[atom], dic_atoms[atom_Y], ic_index, atom_X, atom, atom_Y))
      umat_index = umat_index + 1
      umat_list.append("%3d %2d  1.\n"%(umat_index, ic_index-0))
      ped_index = ped_index + 1
      word = "i%s"%(atom, ) ### ???
      ped_list.append("%3d  %-10s"%(umat_index, word))
###

#for item in ic_list: # debug
#   print item        # debug
#for item in umat_list: # debug
#   print item        # debug
#for item in ped_list: # debug
#   print item        # debug
#sys.exit()  ## debug

### Generate the items about dihedrals in the IC table
list_dihedrals = [['atom1','atom2','atom3','atom4']]  ## make the list un-empty, so that "for dihedral in list_dihedrals" can loop at least once
for atom in list_names3:
   if len(connections[atom]) < 2: continue
   atom2 = atom
   for atom3 in connections[atom2]:
      umat_index = umat_index + 1
      umat_list.append("")
      count = 0
      for atom1 in connections[atom2]:
         if atom1 != atom3:
            for atom4 in connections[atom3]:
               if atom4 != atom2:
                  flag_already = 0
                  for dihedral in list_dihedrals:
                     if ((atom1 == dihedral[0] and atom2 == dihedral[1] and atom3 == dihedral[2] and atom4 == dihedral[3]) or (atom4 == dihedral[0] and atom3 == dihedral[1] and atom2 == dihedral[2] and atom1 == dihedral[3] )):
                        flag_already = 1
                  if flag_already == 0:
                        list_dihedrals.append([atom1, atom2, atom3, atom4])
	 	        ic_index = ic_index + 1
		        ic_list.append(" 4%3d%3d%3d %2d  ! %d %s-%s-%s-%s\n"%(dic_atoms[atom1], dic_atoms[atom2], dic_atoms[atom3], dic_atoms[atom4], ic_index, atom1, atom2, atom3, atom4))
		        umat_list[-1] = umat_list[-1] + "%3d %2d  1."%(umat_index, ic_index)
		        count = count + 1
		        if (count % 4 == 0 and count > 0): umat_list[-1] = umat_list[-1] + "\n"
		        elif (count % 4 != 0): umat_list[-1] = umat_list[-1] + "        "
      if count == 0:
         umat_index = umat_index -1; umat_list = umat_list[:-1]
      else:
         if (count % 4 != 0):
            umat_list[-1] = umat_list[-1][:-8] + "\n"
         ped_index = ped_index + 1
         word = "tor%s-%s"%(atom2, atom3)
         ped_list.append("%3d  %-10s"%(umat_index, word))
###

### debug
#for item in ic_bfly: ## debug
#   print item   ## debug
#for item in ic_list: # debug
#   print item        # debug
#for item in umat_list: # debug
#   print item        # debug
#for item in ped_list: # debug
#   print item        # debug
#sys.exit()

## Open the output file
file3 = "molvib_mm.str"
out = open (file3, "w")
# Write number of freedoms
out.write("MOLVIB NDI1 %d NDI2 %d NDI3 %d SECO\n"%(3*n_atoms-6, 3*n_atoms-6, ic_index))
out.write("GFX\nPRNT    0\n")
out.write("DIM %5d%5d%5d\n"%(3*n_atoms-6, 3*n_atoms-6, ic_index))
# Write the IC table
out.write("IC\n")
for item in ring_IC:
   out.write(item+"\n")
for item in ic_bfly:
   out.write(item)
for item in ic_list:
   out.write(item)
print "IC table has been written."
# Write the UMAT table
out.write("UMAT     0    1    0           ! row normalization\n")
for item in ring_UMAT:
   out.write(item+"\n")
for item in umat_list:
   out.write(item)
out.write("-1\n")
print "UMAT table has been written."
# Write the PED table
count = 0
for item in list_similar:
   if len(item) > 1:
      count = count + 1
out.write("PED   %4d   15\n"%count)
for item in list_similar:
   if len(item) > 1:
      out.write("%3d"%len(item))
      for i in item:
         for umat in umat_list:
            words = umat.split()
            for j in range(0,len(words),3):
               if words[j+1] == i:
	          out.write("%3s"%words[j])
      out.write("\n")
for item in ring_PED:
   out.write(item+"\n")
for i in range(len(ped_list)):
   out.write(ped_list[i])
   if (i % 4 == 0): out.write("\n")
if (i % 4 != 0): out.write("\n")
out.write("-1\nEND\n")
print "PED table has been written."

print "***** MOLVIB_inp_generation is done"
print "***** CAUTION: check molvib_mm.str (especially the PED table) before using it!"
out.close()
sys.exit()
