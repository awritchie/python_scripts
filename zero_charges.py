#! /usr/bin/env python
# Zero the charges in a gromacs .top topology file of atoms given by an index file

import sys
import os
import numpy as np  

atoms_to_zero = "CB HB1 HB2 SG CD NE"
a2z = atoms_to_zero.split()


def usage() :
	print "USAGE: %s <top> <ndx>"%(sys.argv[0])
	sys.exit()

try :
    top = sys.argv[1]
    topfile = open(top,'r')
    toplines = topfile.readlines()
    topfile.close()
except : 
    print "Error reading top file"
    usage()
try :
    ndx = sys.argv[2] 
    ndxi = np.loadtxt(ndx, dtype="int") 
except : 
    print "Error reading index file"
    usage()

# zero charge on SCN
ztop = open(top.replace("TCHG","TCHG-SCN"),'w')
inAtoms = False
for line in toplines :
    if "[ atoms ]" in line : 
        inAtoms = True
        ztop.write(line)
        continue
    elif "[" in line :
        inAtoms = False
        ztop.write(line)
        continue
    else :
        ztop.write(line)
	continue
    if inAtoms and not line.startswith(";") : 
        ll = line.split()
        index = int(ll[0])
        if index in ndxi :
            charge = ll[6]
            ztop.write(line.replace(charge,"%11.6f"%(0.0)))
	    continue
        else : 
	    ztop.write(line)
	    continue
    else :
        ztop.write(line)
	continue
ztop.close()
print "Done writing to %s"%(os.path.abspath(top.replace("TCHG","TCHG-SCN")))
## zero charge on CH2SCN
#ztop = open(top.replace("TCHG","TCHG-CH2SCN"),'w')
#for line in toplines :
#    if "CNC" in line :
#        same = True
#        for each in ["CB","HB1","HB2","SG","CD","NE"] :
#            if " %s "%each in line :
#                charge = line[45:56]
#                ztop.write(line.replace(charge,"%11.6f"%0.0))
#                same = False
#                break
#        if same :
#            ztop.write(line)
#    else :
#        ztop.write(line)
#ztop.close()
## zero charge on Ral, residues 1 to 100
#ztop = open(top.replace("TCHG","TCHG-Ral"),'w')
#inatoms = False
#for line in toplines :
#    same = True
#    if line.startswith("[ ") :
#        inatoms = False
#    if inatoms and line.split() :
#        if not line.startswith(";") :
#            resid = int(line[20:24])
#            if resid > 0 and resid <= 100 :
#                charge = line[45:56]
#                ztop.write(line.replace(charge,"%11.6f"%0.0))
#                same = False
#    if line.startswith("[ atoms ]") :
#        inatoms = True
#    if same :
#        ztop.write(line)
#ztop.close()
## zero charge on solvent and ions
#ztop = open(top.replace("TCHG","TCHG-Sol"),'w')
#inatoms = False
#for line in toplines :
#    same = True
#    if line.startswith("[ ") :
#        inatoms = False
#    if inatoms and line.split() :
#        if not line.startswith(";") :
#            if "Na+" in line or "Cl-" in line :
#                charge = line[36:47]
#               ztop.write(line.replace(charge,"%11.6f"%0.0))
#                same = False
#    if line.startswith("[ atoms ]") :
#        inatoms = True
#    if same :
#        ztop.write(line.replace("tip3p","zero_tip3p"))
#ztop.close()
