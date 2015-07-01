#! /usr/bin/env python

from optparse import OptionParser
import sys

parser=OptionParser()
parser.add_option("-f", \
                  "--structure", \
                  dest="groname", \
                  metavar="GRO",\
                  help="A GROMACS .gro protein structure"\
                  )
parser.add_option("-r",\
                  "--rotationangle",\
                  dest="rotstep",\
                  metavar="ROTSTEP",\
                  help="Rotation step size. Default: 60"\
                  )
parser.add_option("-s",\
                  "--sidechain",\
                  dest="residue",\
                  help="Side chain identifier (residue number in the .gro file)"\
                  )
parser.add_option("-c",\
                  "--fix",\
                  dest="fix_zero",\
                  default=True,\
                  action="store_false",\
                  help="Do not fix it so structure <groname>-0.gro has a dihedral of 0"\
                  )
parser.add_option("-x",\
                  "--chi",\
                  dest="chi",\
                  metavar="X",\
                  help="Dihedral angle to look at (chi: 1, 2, ... etc.).  For residues \
having multiple heavy delta or gama atoms, the last one listed in the .gro is used",\
                  default=1\
                  )
parser.add_option("-t",\
                  "--threshold",\
                  dest="threshold",\
                  metavar="THRESHOLD",\
                  help="Minimum allowable distance (coordinate file units) between atoms.  Default=0.15 nm)",\
                  default=0.15\
                  )

(options, args) = parser.parse_args()
groname = options.groname
rotstep = int(options.rotstep)
resid = options.residue # leaving as a string since it will be read as a string
fix_zero = options.fix_zero
chi = int(options.chi) # converting to int since it's later alled as an integers
threshold=float(options.threshold)

if fix_zero != True :
        fix_zero = False


"""
Taking a GRO file:
1. center atoms around CB (only S=C=N needs to be moved)
2. rotate around that bond in steps of 30 degrees = pi/6
3. move atoms back by centering vector

0         1         2         3         4
01234567890123456789012345678901234567890123456789
   19CNC     SG  270   7.235   8.444   3.168

residue name = [5:8] (3 characters)
atom name = [13:15] (2 characters)
x coords = [20:28] (8)
y coords = [28:36]
z coords = [36:44]

"""
 

import sys
from numpy import pi, array, dot, linalg, cross
from math import cos, sin, atan2

def mag( v ) :
	return dist( v, array((0,0,0)) )

def diff( a1, a2 ):
	return a1-a2

def dist( a1, a2 ) :
	a = diff( a1, a2 )
	return ( dot(a,a) )**.5
def inCoord(line) :
        inline=False
        try :
            float(line.split()[-1])
            inline=True
        except :
            return False
        try :
            float(line.split()[-2])
            inline=True
        except :
            return False
        try :
            float(line.split()[-3])
            inline=True
        except :
            return False
        try :
            float(line.split()[-4])
            inline=True
        except :
            return False
        if inline :
            return True
def getR(atom1, atom2) :
        a1= (atom2[0]-atom1[0])**2
        a2= (atom2[1]-atom1[1])**2
        a3= (atom2[2]-atom1[2])**2
        d=(a1+a2+a3)**.5
        return d

try:
	grofile = open( groname )
except :
	print "usage: %s <groname> <resid>" % sys.argv[0]
	sys.exit()

grolines = grofile.readlines()
grofile.close()
    
names2coords = {}
move_these=False
atomsToMove = []
atomsToKeep = []
print "Looking at the following residue: "
for line in grolines :
    resnum = line[2:5].split()[0]
    if resnum == resid :
        print line.strip()
        resname = line[5:8]
        x = line[20:28]
        y = line[28:36]
        z = line [36:44]
        x = float(x)
        y = float(y)
        z = float(z)
        atomname = line[12:15].split()[0]
        names2coords[atomname]=array(( x,y,z ))
        if atomname == "CB" :
                move_these = True
        elif atomname == "C" :
                move_these = False
        if move_these :
                atomsToMove.append(atomname)
        if "G" in atomname[-2:] and "H" not in atomname :
                gama=atomname
        elif "D" in atomname[-2:] and "H" not in atomname :
                delta=atomname
    elif inCoord(line) and line[5:8] != "SOL" :
        try :
                x = line[20:28]
                y = line[28:36]
                z = line [36:44]
                x = float(x)
                y = float(y)
                z = float(z)
                atomsToKeep.append(array((x,y,z)))
        except :
                print line.strip(), "not included as atom coordinates"


print "\nConverting to array:\n",names2coords
print "\nSide chain atoms are:\n",atomsToMove

atomlists = []
atomlists += atomsToMove

# centering
if chi == 1 :
        centeringVector = -names2coords[ "CA" ]
elif chi == 2 :
        centeringVector = -names2coords[ "CB" ]

newcoords = {}
for atom in atomsToMove :
    oldcoords = names2coords[atom]
    newcoords[atom] = oldcoords+centeringVector
#    print "\n",newcoords[atom]

if chi == 1 :
        # unit vector along the CA-CB axis
        rotationAxis = newcoords["CB"]
elif chi == 2:
        # unit vector along the CB-?G axis
        rotationAxis = newcoords[gama]

( ux, uy, uz ) = rotationAxis/mag(rotationAxis)
ux2 = ux**2 
uxuy = ux*uy
uxuz=ux*uz
uy2 = uy**2
uyuz=uy*uz
uz2 = uz**2

# calculate the chi n angle
if fix_zero:
        if chi == 1 :
                b1 = names2coords["CA"]-names2coords["N"]
                b2 = names2coords["CB"]-names2coords["CA"]
                b3 = names2coords[gama]-names2coords["CB"]
        elif chi == 2 :
                b1 = names2coords["CB"]-names2coords["CA"]
                b2 = names2coords[gama]-names2coords["CB"]
                b3 = names2coords[delta]-names2coords[gama]
        b1 = b1/linalg.norm(b1)
        b2 = b2/linalg.norm(b2)
        b3 = b3/linalg.norm(b3)
        n1 = cross(b1,b2)
        n2 = cross(b2,b3)
        m1 = cross(n1,b2)
        x = dot(n1,n2)
        y = dot(m1,n2)
        dihedral=-atan2(y,x)*180/pi

        print "\n\nDihedral angle is currently: \n<",dihedral,">\nAccomodating...\n"
else :
        dihedral = 0
        
log = groname.replace(".gro",".log")
logfile=open(log,"w")
for rotdegrees in range(0,360,rotstep) :
    distances=[]
    distances.append(0)
    degree_shift=rotdegrees
    tries=0
    change=1
    while min(distances) <= threshold : 
        #print "%d * pi/6" % i
        #print rotdegrees
        distances = []
        angle = (degree_shift-dihedral)*(pi/180) # convert to radians
        #print degree_shift, dihedral, angle*180/pi
        s = sin( angle )
        c = cos( angle )
        rotmat = [ [ ux2 + (1-ux2)*c, uxuy*(1-c) - uz*s, uxuz*(1-c)+uy*s ],
                [ uxuy*(1-c) + uz*s, uy2+(1-uy2)*c, uyuz*(1-c) - ux*s ],
                [ uxuz*(1-c) - uy*s, uyuz*(1-c)+ux*s, uz2+(1-uz2)*c ] ]
        rotmat = array(rotmat)
        newatomcoords = (range(len(atomsToMove)))
        for t in range(1, len(atomsToMove)) :
            atomlists[t]=newcoords[atomsToMove[t]]
#            print "printing old coordinates", atomsToMove[t] 
#            print atomlists[t]
            newatomcoords[t] = dot ( rotmat, atomlists[t] ) - centeringVector
            for coordinate in range(len(atomsToKeep)) :
                this_dist=getR(newatomcoords[t], atomsToKeep[coordinate])
                distances.append(this_dist)
                if this_dist <= threshold :
                        break
        tries += 1
        if tries%40 == 0 and tries < 1000 : 
            print "The rotamer, %i, has been shifted %i times.  It is now in a region being examined by a different rotamer angle.  Consider decreasing the threshold criteria from %.5f." %(rotdegrees, tries, threshold)
        elif tries == 1000 : 
            print "The rotamer, %i, has been rotated 1000 times and is still not above the threshold minimumum distance of %.5f from every atom.  Consider changing some parameters..." %(rotdegrees, threshold)
            sys.exit()
##        if tries >= 40 :
##                print "Error: ",
##                print "The rotamer, %i, has been shifted %i times.  It is now in a region being examined by a different rotamer angle.  Consider changing the threshold criteria." %(rotdegrees,tries)
##                logfile.write("The rotamer, %i, has been shifted %i times.  It is now in a region being examined by a different rotamer angle.  Consider changing the threshold criteria.\n" %(rotdegrees,tries))
##                logfile.close()
##                exit()
#        print min(distances), rotdegrees, degree_shift, change
        this_shift=degree_shift
        if bool(tries & 1) : # odd numbers
                degree_shift = rotdegrees + rotstep/20. * change
        else : # even numbers
                degree_shift = rotdegrees - rotstep/20. * change
                change += 1

#        print "printing new coordinates", atomsToMove[t]
#        print newatomcoords[t]

    print "\nNew coordinates for rotamer %s are:\n"%rotdegrees,newatomcoords
    print "Centered on %.2f (%.2f) after %i steps with a minimum distance of %.3f."%(this_shift,rotdegrees,tries-1,min(distances))
    logfile.write("Centered on %.2f (%.2f) after %i steps with a minimum distance of %.3f.\n"%(this_shift,rotdegrees,tries-1,min(distances)))    
    filename = groname.replace( ".gro", "-%d.gro" % rotdegrees )
    print "Writing to", filename
    file = open( filename, "w" )
    startRotation = False
    for line in grolines :
        resnum = line[2:5].split()[0]
        if resnum == resid :
            atomname = line[12:15].split()[0]
            for t in range(0, len(atomsToMove)):
                if atomname == atomsToMove[t] and startRotation :
                    print "writing", atomsToMove[t], "coordinates to", filename
                    funkyzeit = newatomcoords[t]
                    newline = line [0:20]
                    newline += "%8.3f" % funkyzeit[0]
                    newline += "%8.3f" % funkyzeit[1]
                    newline += "%8.3f" % funkyzeit[2]
                    newline += line[44:]
                    line = newline
                elif chi == 1 and atomname == "CB" or \
                   chi == 2 and atomname == gama :
                        startRotation = True
                        
        file.write( line.rstrip() +"\n" )
    file.close()

logfile.close()
    
#exit()
