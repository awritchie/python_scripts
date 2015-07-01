#! /usr/bin/env python

import numpy
from math import pi, atan2, cos, sin, isnan
from structure_apbs import *
from fileUtilities import *

maxwidth = 79


    
class sidechains_2_change() : 
    def __init__(self, name=None,resid1=None, mutation1=None, resid2=None, mutation2=None) :
        self.name = name
        self.resid1 = resid1
        self.mutation1 = mutation1
        self.resid2 = resid2
        self.mutation2 = mutation2

def pdb_coords( fileline ) :
    return numpy.array((float(fileline[30:38]),float(fileline[38:46]),float(fileline[46:54])))  
    
def pdb_resid( fileline ) : 
    return int(fileline[23:26].split()[0])

def pdb_atomid( fileline ) : 
    return int(fileline[6:11].split()[0])

def pdb_resname( fileline ) : 
    return fileline[16:21].split()[0]
    
def pdb_atomname( fileline ) : 
    return fileline[12:16].split()[0]

def gro_coords( fileline ) :
    x=float(fileline[20:28])
    y=float(fileline[28:36])
    z=float(fileline[36:44])
    return numpy.array((x,y,z))
    
def gro_resid( fileline ) :
    return int(fileline[0:5].split()[0])

def gro_atomid( fileline ) : 
    return int(fileline[15:20].split()[0])
    
def gro_resname( fileline ) : 
    name = fileline[5:9].split()[0]
    if len(name) > 3 : 
        name = name[1:]
    return name

def gro_atomname( fileline ) : 
    return fileline[11:15].split()[0]
 
def pdb_name( filename ) : 
    filename = os.path.abspath(filename)
    if filename[:-4] != '.pdb' and '.' in filename : 
        return filename.replace(filename.split('.')[-1],'pdb')
    else :
        return filename+'.pdb'
    return filename
    
def gro_name( filename ) : 
    filename = os.path.abspath(filename)
    if filename[:-4] != '.gro' and '.' in filename : 
        return filename.replace(filename.split('.')[-1],'gro')
    else :
        return filename+'.gro'
    return filename
        
def gmx_input_names( structureFile ) : 
    structureFile = os.path.abspath(structureFile)
    structureFile = structureFile.replace(".gro","").replace(".pdb","")
    return structureFile, structureFile+".top",structureFile+".itp",structureFile+'.tpr'

def run_tleap( inname, outname ) : 
    cmd = ""
    leapin = inname.replace('.gro','.in').replace('.pdb','.in')
    outname = pdb_name(outname)
    backup_outname(outname)
    lines = "x = loadpdb %s\nsavepdb x %s\nquit\n"%(pdb_name(inname),outname)
    writeLines(leapin,lines)
    deffnm, top, itp, tpr = gmx_input_names( outname )
    if '.gro' in inname : 
        printw("Must use a .pdb file for tleap!")
        sys.exit()
        return None
    cmd += 'tleap -f leaprc.ff03CNC -f %s;\n'%leapin
    cmd += 'echo 1 | pdb2gmx -f %s -o %s -water None -ignh -merge all -p %s -i %s;\n'%(outname,gro_name(outname),top,itp)
    cmd += 'rm %s %s\n'%(top,itp)
    os.system(cmd)
    if isFile(gro_name(outname)) : 
        return gro_name(outname)
    return None

def minimize( structure, mdp, ndx = None, rerun = False, ignh = True, dihres = False ) :
    ignh_string = "-ignh"
    if not ignh :
        ignh_string = ""
    deffnm, top, itp, tpr = gmx_input_names(structure)
    deffnm += '.min'
    outgro = '%s.gro'%deffnm
    if isFile(structure) and isFile(mdp) :
        tmpname = "%s/pdb2gmx.%s"%(os.path.dirname(structure),os.path.basename(structure))
        if not isFile(top) or not isFile(tmpname) or rerun :
            cmd = 'echo 1 | pdb2gmx -f %s -o %s -water None -merge all -p %s -i %s %s;\n'%(structure,tmpname,top,itp,ignh_string)
            if ndx != None :
                cmd += 'genrestr -n %s -o %s -f %s;\n'%(ndx,itp,tmpname)
            os.system(cmd)
        if isFile(top) and isFile(tmpname) :
            if dihres :
                newlines = ""
                hasDih = False
                lines = readFile(top)
                for line in lines :
                    if "; Include Position restraint file" in line :
                        newlines += dihres
                        hasDih = True
                    newlines += line
                if not hasDih :
                    printbox( "Dihedral restraint not added...")
                    sys.exit()
                writeLines(top,newlines)
            if not isFile(tpr) or rerun :
                cmd = 'grompp -f %s -o %s -c %s -p %s;\n '%(mdp,tpr,tmpname,top)
                os.system(cmd)
        if isFile(tpr) :
            cmd = 'mdrun -s %s -deffnm %s;\n'%(tpr,deffnm)
            if not isFile(outgro) or rerun :
                os.system(cmd)
            if isFile(outgro,kill=True) :
                return (outgro,tpr,top)
            else :
                s="Structure file "+isFile(Structure)+" and mdp file "+isFile(mdp)
                printw(s)
    return (outgro,tpr,top)

def mdEnergy( gro, trr, top, mdp, rerun = False) :
    deffnm = nameExt(gro,'.rerun')
    if not isFile("%s.tpr"%deffnm) or rerun :
        cmd = "grompp -f %s -po %s -c %s -t %s -p %s -o %s;\n"%(mdp,deffnm,gro,trr,top,deffnm)
        os.system(cmd)
    if isFile("%s.tpr"%deffnm) :
        if not isFile("%s.edr"%deffnm) :
            cmd = "mdrun -s %s -rerun %s -deffnm %s -maxh 0.05;\n"%(deffnm,gro,deffnm)
            os.system(cmd)
        if isFile("%s.edr"%deffnm) :
            cmd = "echo Potential | g_energy -o %s -f %s -xvg none;\n"%(deffnm,deffnm)
            os.system(cmd)
    line = readFile("%s.xvg"%deffnm)[0]
    energy = float(line.split()[1])
    return (deffnm,energy)
    

def outdated_align_sequence( a, b ) : 
    """
    This is not a rigorous sequence alignment!  It assumes the one structure is 
    a modified version of the other and the only gaps are at the ends!
    """
    m = a.nresidues
    n = b.nresidues
    dif = abs(m-n)
    if m < n : 
        s = a
        l = b
        LongShort = False
    else : 
        l = a
        s = b
        LongShort = True
    shorter = s.sequence
    longer = l.sequence
    bestscore = -1*(len(shorter)+len(longer))
    bestLongoffset = 0
    bestShortoffset = 0
    shortResid = []
    longResid = []
    for k in range(len(shorter)) :
        for i in range(dif) :  
            score = 0
            for j in range(min(n,m)-k) : 
                if longer[i+j] == shorter[j+k] :
                    score += 1
            if score > bestscore : 
                bestLongoffset = i
                bestShortoffset = k
                bestscore = score
    
    for i in range(bestLongoffset) :
        shorter.insert(i,'-')
    for i in range(bestShortoffset) : 
        longer.insert(i,'-')
    for i in range(abs(len(shorter)-len(longer))) : 
        shorter.insert(i+len(shorter),'-')
    for i in range(max(m,n)) : 
        print "%3i %5s %5s"%(i+1,shorter[i],longer[i]),
        if shorter[i][:2] != longer[i][:2] : 
            print " <\n",
        else : 
            print "\n",
    print "Best Score = %i, Best Long Offset = %i, Best Short Offset = %i"%(bestscore,bestLongoffset,bestShortoffset)
    
    if len(longer) != len(shorter) : 
        printw("The lengths of the aligned sequences are off... wtf?")
        return 0
        
    lstop = 0
    sstop = 0
    llist = []
    slist = []
    for i in range(len(longer)) :
        if longer[i] != '-' and shorter[i] != '-' : 
            for j in range(lstop,l.natoms) : 
                if longer[i] == l.resname[j] and l.atom[j] == 'CA' : 
#                    print i,longer[i],j,l.resid[j],l.resname[j],l.atom[j],l.index[j]
                    llist.append(l.resid[j])
                    lstop = j+1
                    break
            for k in range(sstop,s.natoms) : 
                if shorter[i] == s.resname[k] and s.atom[k] == 'CA' : 
#                    print i, shorter[i],k,s.resid[k],s.resname[k],s.atom[k]
                    slist.append(s.resid[k])
                    sstop = k+1
                    break
##    print len(llist),len(slist)
#    for i in range(len(llist)) : 
#        print llist[i],slist[i]

        
    if LongShort : 
        return longer, shorter, llist, slist
    else : 
        return shorter, longer, slist, llist

def diff_ndx( new, original, ndxName ) : 
    newseq, oseq, newresid, oresid  = align_sequence(new,original)
    index = []
    for i in range(new.natoms) :
        resid = new.resid[i]
        for j in range(len(newresid)) :
            if resid == newresid[j] :
                # This tells us that the residue was part of the aligned sequence
                # and we still need to know if it's a changed or unchanged residue
                if oseq[resid-1][:2] != newseq[resid-1][:2] or (oseq[resid-2] == "-" and oseq[resid] == "-") or (newseq[resid-2] == "-" and newseq[resid] == "-") :
                    pass
                    #print new.atom[i],resid,new.resname[i],newseq[resid-1],oseq[resid-1]
                elif 'H' != new.atom[i][0] :
                    index.append(new.index[i])
    s = "[ unperturbed ]\n"
    for i in index : 
        s += "%i\n"%(i) # gros start couting at one, writeIndex starts counting at zero
    backup_outname(ndxName)
    writeLines(ndxName,s)
    return ndxName

def kabsch_alignment( pose1, pose2 , pose1sel = [], pose2sel = [], ndx=False ):
    """
    Performs the Kabsch alignment algorithm upon the N CA C O of two selections.
   
    This is a hacky modification for usage with class structure,
 
    @param pose1: First protein
    @param pose2: Second protein
    @param pose1sel: List of residues to align in the first protein
    @param pose2sel: List of residues to align in the second protein
    """	

    # default to all the residues
    if not pose1sel:
        pose1sel = range(1,pose1.nresidues+1)
    if not pose2sel:
        pose2sel = range(1,pose2.nresidues+1)

    # obtain a list of the atom positions to align
#    atoms_to_fit = ['CA','N','C','O','OC1']
    atoms_to_fit = ['CA','N','C']
    stsel1 = pose1.extract_coordinates(pose1sel,atoms_to_fit,ndx=ndx)
    stsel2 = pose2.extract_coordinates(pose2sel,atoms_to_fit,ndx=ndx)

    # remember all of each pose's atom positions, to apply transformation
    #    later for updating
    molsel1 = pose1.coord
    molsel2 = pose2.coord
    
    # check for consistency, same number of selected residues
    assert len(stsel1) == len(stsel2)
    L = len(stsel1)
    assert L > 0
    
    # must alway center the two proteins to avoid affine transformations
    #     Center the two proteins to their selections
    COM1 = numpy.sum(stsel1,axis=0) / float(L)
    COM2 = numpy.sum(stsel2,axis=0) / float(L)
    stsel1 -= COM1
    stsel2 -= COM2

    # Initial residual, see Kabsch.
    E0 = numpy.sum( numpy.sum(stsel1 * stsel1,axis=0),axis=0) + numpy.sum( numpy.sum(stsel2 * stsel2,axis=0),axis=0)

    # This beautiful step provides the answer.  V and Wt are the orthonormal
    # bases that when multiplied by each other give us the rotation matrix, U.
    # S, (Sigma, from SVD) provides us with the error!  Isn't SVD great!
    V, S, Wt = numpy.linalg.svd( numpy.dot( numpy.transpose(stsel2), stsel1))
    # we already have our solution, in the results from SVD.
    # we just need to check for reflections and then produce
    # the rotation.  V and Wt are orthonormal, so their det's
    # are +/-1.
    reflect = float(str(float(numpy.linalg.det(V) * numpy.linalg.det(Wt))))
    if reflect == -1.0:
        S[-1] = -S[-1]
        V[:,-1] = -V[:,-1]

    RMSD = E0 - (2.0 * sum(S))
    RMSD = numpy.sqrt(abs(RMSD / L))

    #U is simply V*Wt
    U = numpy.dot(V, Wt)
    # rotate and translate the molecule
    stsel2 = numpy.dot((molsel2 - COM2), U)
    stsel2 = stsel2.tolist()
    # center the molecule
#    stsel1 = molsel1 - COM1
#    stsel1 = stsel1.tolist()

    center_on_COM1 = COM1
    if pose1.ext == "pdb" and pose2.ext == "gro" :
        center_on_COM1 /= 10.
    elif pose1.ext == "gro" and pose2.ext == "pdb" :
        center_on_COM1 *= 10.

    # apply the changes to both poses
#    for i in range(pose1.natoms) :
#        pose1.coord[i]=stsel1[i]           # centered on COM
#        pose1.coord[i]=stsel1[i] + COM1    # return to original position
    for i in range(pose2.natoms) :
        #pose2.coord[i]=stsel2[i]           # centered on COM
        pose2.coord[i]=stsel2[i] + center_on_COM1
    
    pose1base = os.path.basename(pose1.name.replace('.%s'%pose1.ext,''))
    pose2base = os.path.basename(pose2.name.replace('.%s'%pose2.ext,''))
    pose1name = pose1.name.replace('.%s'%pose1.ext,'.fit2.%s.%s'%(pose2base,pose1.ext))
    pose2name = pose2.name.replace('.%s'%pose2.ext,'.fit2.%s.%s'%(pose1base,pose2.ext))
#    backup_outname(pose1name)
    backup_outname(pose2name)
#    if pose1.ext == "gro" :
#        pose1name = pose1.write_gro(pose1name)
#    elif pose1.ext == "pdb" :
#        pose1name = pose1.write_pdb(pose1name)
    if pose2.ext == "gro" :
        pose2name = pose2.write_gro(pose2name)
    elif pose2.ext == "pdb" :
        pose2name = pose2.write_pdb(pose2name)
    print 'RMSD=%f' % RMSD
    return pose1name, pose2name

    
class structure() : 
    def __init__(self, name = None, ignore_solvent = False):
        printbox("\n\nReading %s\n\n"%name)
        self.name = name
        self.ext = name.split('.')[-1]
        self.line = []
        self.atom = []
        self.coord = []
        self.index = []
        self.resid = []
        self.resname = []
        self.perturbed = []
        self.writeIndex = []
        self.natoms = 0
        self.sequence = []
        self.sequence_0 = []
        self.nterm = []
        self.cterm = []
        lines = readFile(self.name)
        at_nterm = True
        for line in lines :
            if ("SOL" in line or "HOH" in line) and ignore_solvent :
                break
            if ("Na+" in line or "Cl-" in line) and ignore_solvent :
                break
            if "TER" not in line and "END" not in line and "ENDMDL" not in line :
                if self.ext == 'gro' :
                    try :
                        #print line,
                        #print "0123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 "
                        #float(line.split()[5])
                        float(line[36:41])
                        #print line.strip()
                        self.line.append(line)
                        self.atom.append(gro_atomname(line))
                        self.coord.append(gro_coords(line))
                        self.index.append(gro_atomid(line))
                        self.resid.append(gro_resid(line))
                        self.resname.append(gro_resname(line))
                        self.writeIndex.append(self.natoms)
                        #print gro_atomname(line),
                        #print gro_coords(line),
                        #print gro_atomid(line),
                        #print gro_resid(line),
                        #print gro_resname(line),
                        #print self.natoms
                        if "OC2" == gro_atomname(line) and not at_nterm : 
                            self.cterm.append(gro_resid(line))
                            at_nterm = True
                        if at_nterm and gro_atomname(line) == "CA" : 
                            self.nterm.append(gro_resid(line))
                            at_nterm = False
                        self.natoms += 1
                    except (ValueError, IndexError) :
                        pass
                elif self.ext == 'pdb' :
                    if ('ATOM' in line or 'HETATM' in line) : 
                        try :
                            l=[line[31:38],line[38:46],line[46:54],line[54:60],line[60:66]]
                            for value in l :
                                """print "value=",value
                                print "float-",float(value)
                                print "\n"""
                                float(value)
                            self.line.append(line)
                            self.atom.append(pdb_atomname(line))
                            self.coord.append(pdb_coords(line))
                            self.index.append(pdb_atomid(line))
                            self.resid.append(pdb_resid(line))
                            self.resname.append(pdb_resname(line))
                            self.writeIndex.append(self.natoms)
                            if ("OC2" == pdb_atomname(line) or "OXT" == pdb_atomname(line)) and not at_nterm : 
                                self.cterm.append(pdb_resid(line))
                                at_nterm = True
                            if at_nterm and pdb_atomname(line) == "CA" : 
                                self.nterm.append(pdb_resid(line))
                                at_nterm = False
                            self.natoms += 1
                        except (ValueError, IndexError) : 
                            pass
                    if line.startswith('END') : 
                        # Guarantee we only get the first structure if multiple present
                        break
                else : 
                    printw('%s does not have a .gro or .pdb file type extension'%self.name)
                    break
                    sys.exit()
                if self.natoms > 1 and self.resid[self.natoms-2] != self.resid[self.natoms-1] : 
                    self.sequence.append(self.resname[self.natoms-1])
                elif self.natoms == 1 : 
                    self.sequence.append(self.resname[0])
        if self.ext == 'gro' : 
            self.groHeader = lines[0]
            self.groFooter = lines[-1]
        else : 
            self.groHeader = "TITLE Based on %s\n"%self.name
            self.groFooter = "%10.5f%10.5f%10.5f\n"%(0.,0.,0.)

        self.line = numpy.array(self.line)
        self.atom = numpy.array(self.atom)     
        self.index = numpy.array(self.index)
        self.coord = numpy.array(self.coord)
        self.resid = numpy.array(self.resid)
        self.resname = numpy.array(self.resname)
        
        # Copy the original structure 
        self.line_0 = list(self.line)
        self.atom_0 = list(self.atom)
        self.index_0 = list(self.index)
        self.coord_0 = list(self.coord)
        self.resid_0 = list(self.resid)
        self.resname_0 = list(self.resname)
        self.sequence_0 = list(self.sequence)
        self.extreme_coords()
        self.nresidues = len(self.sequence)
        
    def newline(self,n,index=None,resid=None,resname=None,atom=None,coord=(float('nan'),float('nan'),float('nan')),newext=None) :    
        if index == None : index = self.index[n]        
        if resid == None : resid = self.resid[n]
        if resname == None : resname = self.resname[n]
        if atom == None : atom = self.atom[n]
        if isnan(coord[0]) : coord = self.coord[n]
        if newext == None : newext = self.ext
        if newext == 'gro' : 
            if self.ext == 'pdb' : 
                coord = coord/10.
            if resid in self.cterm : 
                if resname[0] != "C" and len(resname) < 4 : 
                    resname = "C%s"%resname
            elif resid in self.nterm : 
                if resname[0] != "N" and len(resname) < 4 : 
                    resname = "N%s"%resname
            elif len(resname) > 3 : 
                resname = resname[1:]
            newline = "%5i%-5s%5s%5i%8.3f%8.3f%8.3f\n"%(resid,resname,atom,index,coord[0],coord[1],coord[2])
            return newline
        elif newext == 'pdb' : 
            if self.ext == 'gro' : 
                coord = coord*10.
            if len(resname) > 3 : 
                resname = resname[1:]
            newline = "%-6s%5i %4s %-4s  %3i    %8.3f%8.3f%8.3f%6.2f%6.2f\n"%('ATOM',index,atom,resname,resid,coord[0],coord[1],coord[2],1.,0.)
            return newline
        else : 
            printw('%s does not have a .gro or .pdb file type extension'%self.name)
        return 0
    
    def dihedral(self,index) : 
        if len(index) != 4 : 
            printw("There are %i elements in the dihedral index.  There should only be 4!"%len(index))
            return False
        index = numpy.array(index) - 1 
        b1 = self.coord[index[1]] - self.coord[index[0]]
        b2 = self.coord[index[2]] - self.coord[index[1]]
        b3 = self.coord[index[3]] - self.coord[index[2]]
        b1 /= numpy.linalg.norm(b1)
        b2 /= numpy.linalg.norm(b2)
        b3 /= numpy.linalg.norm(b3)
        n1 = numpy.cross(b1,b2)
        n2 = numpy.cross(b2,b3)
        m1 = numpy.cross(n1,b2)
        x = numpy.dot(n1,n2)
        y = numpy.dot(m1,n2)
        dihedral = -atan2(y,x)*180/pi
        return dihedral
        
    def rotate(self,new_dihedral,dih_ndx,index) :
        print dih_ndx,
        dih_ndx = sorted(dih_ndx) # I think the relevant dihedrals wll ALWAYS be indexed smallest to largest
        print dih_ndx
        index = sorted(index)
        cur_dih = self.dihedral(dih_ndx)
        if not cur_dih : 
            return False
        for i in dih_ndx[1:-1] : # If the middle 2 atoms are in the index, remove them
            try : 
                index.remove(i)
            except ValueError :
                pass
        if dih_ndx[0] in index and dih_ndx[-1] in index : 
            printw("Both ends of the dihedral definition atom are in the index to-be rotated.  Only 1 should be present...")
            return False
        angle = ( new_dihedral - cur_dih )*(pi/180)
        printbox("From %.3f to %.3f by %.3f"%(cur_dih,new_dihedral,angle*180/pi))
        s = sin( angle )
        c = cos( angle )
        if 0 in index : 
            index.remove(0)
        index = numpy.array(index) - 1
        centeringVector = -1 * self.coord[dih_ndx[1]-1]
        rotationAxis = self.coord[dih_ndx[2]-1] + centeringVector
        printcenter("Current Dihedral angle about \n<%s(%i%s-%i) - %s(%i%s-%i) - %s(%i%s-%i) - %s(%i%s-%i)>\nis %.3f degrees"\
        %(self.atom[dih_ndx[0]-1],self.resid[dih_ndx[0]-1],self.resname[dih_ndx[0]-1],dih_ndx[0],\
        self.atom[dih_ndx[1]-1],self.resid[dih_ndx[1]-1],self.resname[dih_ndx[1]-1],dih_ndx[1],\
        self.atom[dih_ndx[2]-1],self.resid[dih_ndx[2]-1],self.resname[dih_ndx[2]-1],dih_ndx[2],\
        self.atom[dih_ndx[3]-1],self.resid[dih_ndx[3]-1],self.resname[dih_ndx[3]-1],dih_ndx[3],cur_dih))
        printw("\nCentering about %s(%i%s)"%(self.atom[dih_ndx[1]-1],self.resid[dih_ndx[1]-1],self.resname[dih_ndx[1]-1]))
        printw("Rotation about %s(%i%s) - %s(%i%s)\n"%(self.atom[dih_ndx[1]-1],self.resid[dih_ndx[1]-1],self.resname[dih_ndx[1]-1],\
        self.atom[dih_ndx[2]-1],self.resid[dih_ndx[2]-1],self.resname[dih_ndx[2]-1]))
        coords = [] 
        for i in range(len(index)) : 
            self.perturbed.append(index[i]+1)
            coords.append( numpy.array(self.coord[index[i]]) + numpy.array(centeringVector))
        ( ux, uy, uz ) = rotationAxis/numpy.linalg.norm(rotationAxis)
        ux2 = ux**2
        uxuy = ux*uy
        uxuz = ux*uz
        uy2 = uy**2
        uyuz = uy*uz
        uz2 = uz**2
        rotmat = [ [ ux2 + (1-ux2)*c, uxuy*(1-c) - uz*s, uxuz*(1-c) + uy*s ],
                   [ uxuy*(1-c) + uz*s, uy2 + (1-uy2)*c, uyuz*(1-c) - ux*s ],
                   [ uxuz*(1-c) - uy*s, uyuz*(1-c) + ux*s, uz2 + (1-uz2)*c ] ]
        rotmat = numpy.array(rotmat)
        for i in self.writeIndex : 
            for j in range(len(index)) :
                if self.index[i] == index[j]+1 : 
                    coords[j] = numpy.dot(rotmat,coords[j]) - centeringVector
                    self.line[i] = self.newline(index[j],coord=coords[j])
                    self.coord[i] = coords[j]
                    break
        return True
        
    def replace_residue(self,resnum,replacement) : 
        if replacement != 'GLY' : 
            backbone=numpy.array(( "N"\
                ,"CA"\
                ,"C"\
                ,"O"\
                ,"OXT"\
                ,"CB"\
                ))
        else : 
            backbone=numpy.array(( "N"\
                ,"CA"\
                ,"C"\
                ,"O"\
                ,"OXT"\
                ))
        mutated = False
        for i in range(self.natoms) : 
            if self.resid[i] == resnum : 
                if self.atom[i] in backbone : 
                    self.resname[i] = replacement
                    self.line[i] = self.newline(i,resname=replacement)
                    printw("Mutated <%s(%i%s) to %s(%i%s)>"%(self.atom_0[i],self.resid_0[i],self.resname_0[i],self.atom[i],self.resid[i],self.resname[i]))
                    mutated = True
                    #print self.line_0[i],self.line[i]
                else :
                    self.writeIndex.remove(i)
        return mutated
    
    def met_2_cnc(self,resnum) : 
        """Specialized version of replace_residue intended to turn the MET into a CNC"""
        self.sequence[resnum - 1 ] = 'CNC'
        sameNames = numpy.array(("N"\
                          ,"H"\
                          ,"CA"\
                          ,"HA"\
                          ,"CB"\
                          ,"HB2"\
                          ,"HB3"\
                          ,"C"\
                          ,"O"\
                          ))
        remove = []
        sg = False
        cd = False
        if self.ext == 'gro' : 
            cf = 0.1
        for i in range(self.natoms) : 
            if self.resid[i] == resnum : 
                self.perturbed.append(i+1)
                self.resname[i] = 'CNC'
                self.line[i] = self.newline(i)
                if self.atom[i] == "CG" : 
                    self.atom[i] = "SG"
                    self.line[i] = self.newline(i)
                    sg = True
                    SG = self.coord[i]
                elif self.atom[i] == "SD" :
                    self.atom[i] = "CD"
                    self.line[i] = self.newline(i)
                    if sg : 
                        cd = True
                        CD = SG + (self.coord[i] - SG)*1.679/numpy.linalg.norm(self.coord[i] - SG) * cf # 1.679 Angstrom is SG-CD bond lngth in Amber03CNC
                        self.coord[i] = CD
                elif self.atom[i] == "CE" : 
                    self.atom[i] = "NE"
                    self.line[i] = self.newline(i)
                    if sg and cd : 
                        NE = CD + (CD - SG)*1.138/numpy.linalg.norm(CD - SG) * cf # 1.138 Angstrom is CD-NE triple bond length in Amber03CNC
                        self.coord[i] = NE
                elif self.atom[i] not in sameNames : 
                    remove.append(self.writeIndex[i])
        for i in remove : 
            if i in self.writeIndex : 
                self.writeIndex.remove(i)
        
    def write_gro(self,filename, maxresnum = []) :
        if not maxresnum :
            maxresnum = 2*self.natoms # in case atoms have been added for some reason
        filename = gro_name(filename)
        backup_outname(filename)
        filelines = self.groHeader
        filelines += " NATOMS\n"
        natoms = 0
        for i in self.writeIndex :
            if self.resid[i] <= maxresnum :
                filelines += self.newline(i,newext='gro')
                natoms += 1
        filelines = filelines.replace("NATOMS","%i"%natoms)
        filelines += self.groFooter
        file = open(filename,'w')
        file.write(filelines)
        file.close()
        printw("Done writing structure to %s"%filename)
        return filename
    
    def write_pdb(self,filename, maxresnum = []):
        """--ALWAYS-- use pdb2gmx -ignh or tleap on these pdb fies!!"""
        if not maxresnum :
            maxresnum = 2*self.natoms # in case atoms have been added for some reason
        filename = pdb_name(filename)
        backup_outname(filename) 
        filelines = "TITLE\t%s\n"%filename
        for i in self.writeIndex :
            if self.resid[i] <= maxresnum :
                if self.atom[i][0] != "H" :
                    filelines += self.newline(i,newext='pdb')
                if self.atom[i] == 'OC2' :
                    filelines += 'TER\n'
        filelines += "TER\nENDMDL"
        file = open(filename,'w')
        file.write(filelines)
        file.close()
        printw("Done writing structure to %s"%filename)
        return filename
    
    def extreme_coords(self) : 
        xs = []
        ys = []
        zs = []
        for i in self.writeIndex : 
            coord = self.coord[i]
            xs.append(coord[0])
            ys.append(coord[1])
            zs.append(coord[2])
        self.min=numpy.array((min(xs),min(ys),min(zs)))
        self.max=numpy.array((max(xs),max(ys),max(zs)))

    def extract_coordinates(self,residues,atomnames=[],ndx=False) :
        coords = []
        if ndx :
            for i in range(self.natoms) :
                if self.index[i] in residues :
                    coords.append(self.coord[i])
        else :
            if not atomnames :
                atomnames = ['CA']
            for i in range(self.natoms) :
                if self.resid[i] in residues :
                    if self.atom[i] in atomnames :
                        coords.append(self.coord[i])
        return numpy.array(coords)

    def extract_index(self,residues,atomnames=[]) :
        ndxs = []
        if not atomnames :
            atomnames = ['CA']
        for i in range(self.natoms) :
            if self.resid[i] in residues :
                if self.atom[i] in atomnames :
                    ndxs.append(self.index[i])
        return numpy.array(ndxs)

def merge_structures(pose1, pose2, outname, residue1sel = [], residue2sel = [], new_cterm = False, position = 0 ) :
    """
    Use new_nterm = True if you are adding residues to the n-terminus.
    For adding new chain, use new_cterm = True, as the new chain has a new c-terminus.
	
	Use position to insert the structure in front of resid position+1.  
    IE:
	pose1 = GAMGS
            12345

	pose2 = VTRKL IIVST EDEKY
            12345 67890 12345
	Desired outcome : VTRKL GAMGS IIVST EDEKY
                      12345 67890 12345 67890
	usage: merge_structures(pose1,pose2,outname,range(5),range(20),new_nterm=False,position=5)

    BUG:
    > new_cterm is too vague.  It will change ALL OC1 and OC2 in pose2 prior to the insertion 
      location of pose1, ALL OC1 and OC2 in pose2, and NONE after the insertion.  Probably 
      good enough?
    """
    ext='gro'
    if not residue1sel:
        residue1sel = range(1,pose1.nresidues+1)   
    else : 
        residue1sel = numpy.sort(residue1sel)
    if not residue2sel:
        residue2sel = range(1,pose2.nresidues+1)
    else : 
        reisdue2sel = numpy.sort(residue2sel)
    line = ""
    n = 1
    resid = 0 
    prev_resid = 0  
    for i in range(pose2.natoms) : 
        if pose2.resid[i] in residue2sel and pose2.resid[i] < position : 
            if pose2.resid[i] != prev_resid :
                prev_resid = pose2.resid[i]
                resid += 1
            if new_cterm : 
                if pose2.atom[i] == 'OC2' : 
                    line += pose2.newline(i,index=n,resid=resid,atom='O',newext=ext)
                    n += 1
                elif pose2.atom[i] != 'OC1' : 
                    line += pose2.newline(i,index=n,resid=resid,newext=ext)
                    n += 1
            else : 
                line += pose2.newline(i,index=n,resid=resid,newext=ext)
                n += 1
    prev_resid = 0 
    for i in range(pose1.natoms) : 
        if pose1.resid[i] in residue1sel : 
            if pose1.resid[i] != prev_resid : 
                prev_resid = pose1.resid[i]
                resid += 1
            if not new_cterm: 
                if pose1.atom[i] == 'OC2' : 
                    line += pose1.newline(i,index=n,resid=resid,atom='O',newext=ext)
                    n += 1
                elif pose1.atom[i] != 'OC1' : 
                    line += pose1.newline(i,index=n,resid=resid,newext=ext)
                    n += 1
            else : 
                line += pose1.newline(i,index=n,resid=resid,newext=ext)
                n += 1
    prev_resid = 0
    for i in range(pose2.natoms) : 
        if pose2.resid[i] in residue2sel and pose2.resid[i] >= position : 
            if pose2.resid[i] != prev_resid : 
                prev_resid = pose2.resid[i]
                resid += 1
            line += pose2.newline(i,index=n,resid=resid,newext=ext)
            n += 1

    if ext == 'gro' :
        outname = gro_name(outname)
        line = 'TITLE %s\n %i\n%s'%(outname,n-1,line)
        line +="%10f%10f%10f\n"%(5.,5.,5.)
    elif ext == 'pdf' : 
        outname = pdb_name(outname)
    writeLines(outname,line)
    return outname
    
    
def align_sequence( struct1, struct2, gscore=-3, match_HIS2HIx = True ) : 
    """
    A Local alignment of 2 structures.  A gap score (gscore) = -3 appears
    to work fairly well for my limited test structures. If the chain lengths 
    are very similar, I imagine you could use a lower gscore, but I don't 
    see why it's necessary.
    """
    if match_HIS2HIx : 
        his = 2
    else : 
        his = 3

    aseq = struct1.sequence
    bseq = struct2.sequence
    for i in range(1) : 
        aseq.insert(i,'-')
        bseq.insert(i,'-')
    aseq = numpy.array(aseq)
    bseq = numpy.array(bseq)
    ai = len(aseq)
    bi = len(bseq) 

    mat = numpy.zeros(shape=(ai,bi))
    for a in range(ai) : 
        for b in range(bi) : 
            if aseq[a] == '-' : 
                score = b*gscore
            elif bseq[b] == '-' : 
                score = a*gscore
            elif aseq[a][:his] == bseq[b][:his] : 
                score = 2
            else: 
                score = -1
            mat[a][b] = score
    """
    This is the dynamic programming matrix.
    """
    """
    print "%4s,"%"",
    for i in range(ai) :
        print "%4s,"%aseq[i],
    print "\n",
    for j in range(bi) :
        print "%4s,"%(bseq[j]),
        for i in range(ai) :
            print "%4i,"%mat[i][j],
        print "\n",
    print "\n"
    """
    fillmat = numpy.zeros(shape=(ai,bi)) 
    for a in range(1,ai) : 
        for b in range(1,bi) : 
            score = ( mat[a][b]+fillmat[a-1][b-1], fillmat[a][b-1]+gscore, fillmat[a-1][b]+gscore, 0 )
            fillmat[a][b] = max(score)
    del(mat)
    """
    This is the scoring matrix that you follow the diagnals of
    <http://www.avatar.se/molbioinfo2001/dynprog/dynamic.html>
    """
    """
    print "%3s,"%"",
    for i in range(ai) :
        print "%3s,"%aseq[i],
    print "\n",
    for j in range(bi) :
        print "%3s,"%(bseq[j]),
        for i in range(ai) :
            print "%3i,"%(fillmat[i][j]),
#            if i > 20 : break
        print "\n",
#        if j > 20 : break
    print "\n"
    """
    # I need to know the location of the highest score; thats where we begin
    # the alignment
    imax = 0
    jmax = 0
    maxmat = 1e-10
    for i in range(ai) : 
        for j in range(bi) : 
            if fillmat[i][j] > maxmat : 
                imax = i
                jmax = j
                maxmat = fillmat[i][j]
#    print imax,jmax

    i=imax
    j=jmax
    aligned_a = []
    aligned_b = []
    while (i > 1) and (j > 1) : 
#        while j > 1 :   
        up = fillmat[i][j-1]
        left = fillmat[i-1][j]
        diag = fillmat[i-1][j-1]
        best = max(up,left,diag)
        if best == diag : 
            i = i-1
            j = j-1
            aligned_a.append(aseq[i])
            aligned_b.append(bseq[j])
        elif best == up : 
            i = i
            j = j-1
            aligned_a.append('-')
            aligned_b.append(bseq[j])
        elif best == left : 
            i = i-1
            j = j
            aligned_a.append(aseq[i])
            aligned_b.append('-')
    del(fillmat)
    # Keeping these for generating a list of aligned residue numbers
    aFirst = i
    bFirst = j
    aligned_struct1 = list(aligned_a)
    aligned_struct2 = list(aligned_b)
    # Add back in the residues from the beginning of the sequence
    for k in range(i-1,0,-1) : 
        aligned_struct1.append(aseq[k])
        aligned_struct2.append('-')
    for k in range(j-1,0,-1) : 
        aligned_struct1.append('-')
        aligned_struct2.append(bseq[k])
    # Reverse the order so it starts on the lowest resid 
    aligned_a = aligned_a[::-1]
    aligned_b = aligned_b[::-1]
    aligned_struct1 = aligned_struct1[::-1]
    aligned_struct2 = aligned_struct2[::-1]
    # Add back in the unalined residues
    for k in range(imax,ai) : 
        aligned_struct1.append(aseq[k])
    for k in range(jmax,bi) : 
        aligned_struct2.append(bseq[k])
    # Now add gaps to the end of the shorter thing
    for k in range(len(aligned_struct1),len(aligned_struct2)) : 
        aligned_struct1.append('-')
    for k in range(len(aligned_struct2),len(aligned_struct1)) :
        aligned_struct2.append('-')
    total = len(aligned_struct1)
    if len(aligned_struct2) != total: 
        printw("The sequence lengths are different, wtf? %s = %i, %s = %i"%(struct1.name,total,struct2.name,len(aligned_struct2)))
        return 0
    for k in range(total) : 
        if aligned_struct1[k][:his] != aligned_struct2[k][:his] : 
            printw( "%4i%6s%6s <"%(k+1,aligned_struct1[k],aligned_struct2[k]))
        else :
            printw( "%4i%6s%6s"%(k+1,aligned_struct1[k],aligned_struct2[k]))
            
    start = aFirst
    struct1_resids = []
    for k in range(len(aligned_a)) :
        for l in range(start,struct1.natoms) : 
            if struct1.atom[l] == 'CA' and struct1.resname[l] == aligned_a[k] : 
                struct1_resids.append(struct1.resid[l])
                start = l+1
                break
    start = bFirst
    struct2_resids = []
    for k in range(len(aligned_b)) :
        for l in range(start,struct2.natoms) : 
            if struct2.atom[l] == 'CA' and struct2.resname[l] == aligned_b[k] : 
                struct2_resids.append(struct2.resid[l])
                start = l+1
                break
    return aligned_struct1, aligned_struct2, struct1_resids, struct2_resids

def setup_nitrile_dih_rotation( structure, chi = 1) :
    dih_ndx = []
    rotation_ndx = []
    if chi == 1 :
        x_dih = [ 'N', 'CA', 'CB', 'SG' ]
        x_rot = [ 'HB1', 'HB2', 'SG', 'CD', 'NE' ]
    elif chi == 2 :
        x_dih = [ 'CA', 'CB', 'SG', 'CD' ]
        x_rot = [ 'CD', 'NE']
    else :
        printw("ERROR, chi %i was selected for CNC.  We can only do chi 1 and chi 2!"%chi)
        sys.exit()
    for i in range(structure.natoms) :
        if structure.resname[i] == 'CNC' :
            if structure.atom[i] in x_dih :
                dih_ndx.append(structure.index[i])
            if structure.atom[i] in x_rot :
                rotation_ndx.append(structure.index[i])
    if len(dih_ndx) == 4 :
        return dih_ndx,rotation_ndx
    else :
        printw("Too many atoms in the dihedral! Exiting...")
        sys.exit()

