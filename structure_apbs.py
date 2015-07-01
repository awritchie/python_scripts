#! /usr/local/bin/python

from fileUtilities import *
import os
import numpy
import glob

bctemplate="\
read \n\
    mol pqr PQRNAME \n\
end \n\
\n\
# stage 4 \n\
# based on suggestions from Gernot Kierseritzky and stage4 described \n\
# in stages-logic.txt; changed to lpbe for boundary and convergence \n\
# issues when we go to finer grids, and spl2 because of fewer off- \n\
# the-grid warnings. \n\
elec \n\
    mg-manual \n\
\n\
    dime 97 97 97 \n\
    glen 250.0 250.0 250.0 # 2.5 Angstroms \n\
\n\
    gcent mol 1 \n\
    mol 1 \n\
    lpbe \n\
    bcfl sdh # faster than mdh and probably good enough if the coarse grid is large enough \n\
\n\
    ion charge 1 conc 0.150 radius 2.0 \n\
    ion charge -1 conc 0.150 radius 2.0 \n\
    pdie PDIE \n\
    sdie SDIE \n\
    chgm spl2 \n\
    srfm mol \n\
    srad 1.4 \n\
    sdens 10.0 \n\
    temp 300.0 \n\
    calcenergy total \n\
    calcforce no \n\
end \n"

template="\
# stage2 \n\
# 193 193 193 \n\
# 10 10 10 \n\
elec \n\
    mg-manual \n\
\n\
    glen 10 10 10 \n\
    dime 193 193 193 \n\
    gcent CDCOORDS \n\
    mol 1 \n\
    lpbe \n\
    bcfl focus \n\
\n\
    ion charge 1 conc 0.150 radius 2.0 \n\
    ion charge -1 conc 0.150 radius 2.0 \n\
    pdie PDIE \n\
    sdie SDIE \n\
    chgm spl2 \n\
    srfm mol \n\
    srad 1.4 \n\
    sdens 10.0 \n\
    temp 300.0 \n\
    calcenergy total \n\
    calcforce no \n\
    write atompot flat OUTPUT \n\
end \n"

def makeXYZ( trr, tpr, top, ndx, out, parm) :
    cmd = "/Users/ritchie/Utilities/new_gromacs/gmx2xyz-4.6.5/gmx2xyz-4.6.5/gmx2xyz -s %s -f %s -p %s -n %s -a %s -o %s"%(tpr,trr,top,ndx,parm,out)
    os.system(cmd)
    globsearch = out.replace('.xyz','*.xyz')
    xyznames = glob.glob(globsearch)
    return (xyznames)

def makePQR( trr, tpr, top, ndx, outpqr, outfield, pdie, dat) :
    lines = readFile(top)
    atoms_to_zero = []
    for line in lines :
        if "CNC" in line :
            l = line.split()
            if l[4] == "CD" :
                CD = int(l[5])
            elif l[4] == "NE" :
                NE = int(l[5])
            elif l[4] in ['CB','HB1','HB2','HB3','SG'] :
                atoms_to_zero.append(int(l[5]))
    zerostring = '\"'
    for atom in atoms_to_zero :
        zerostring += "%i "%atom
    zerostring += '\"'
    cmd = "/Users/ritchie/Utilities/new_gromacs/gmx2pqr-4.6.5/gmx2pqr-4.6.5/gmx2pqr -s %s -f %s -d %s -o %s -of %s -n %s -a1 %i -a2 %i -a3 %s -pdie %s -xvg none"%(tpr,trr,dat,outpqr,outfield,ndx,CD,NE,zerostring,pdie)
    os.system(cmd)
    globsearch = outpqr.replace('.pqr','*.pqr')
    pqrnames = glob.glob(globsearch)
    return (pqrnames)

def makeAPBSin( pqr, pdie, sdie=78 ) :
    # get the midpoint and bondlength
    lines = readFile(pqr)
    haveCD = False
    haveNE = False
    nDUM = 0
    for line in lines :
        if "DUM" in line :
            nDUM += 1
        if "CNC" in line :
            l=line.split()
            if l[2] == "CD" :
                cd = numpy.array([ float(l[5]), float(l[6]), float(l[7]) ])
                haveCD = True
            elif l[2] == "NE" :
                ne = numpy.array([ float(l[5]), float(l[6]), float(l[7]) ])
                haveNE = True
        if haveCD and haveNE :
            break
    midpoint = cd/2. + ne/2.
    bondlength = numpy.linalg.norm(ne-cd)
    do0 = True
    do1 = True
    # Configure file names for SDIE != PDIE
    inputname0 = nameExt(pqr,'.%d.in'%sdie)
    epot = os.path.basename(nameExt(pqr,'.%d_pot'%sdie))
    atompot0 = os.path.basename(nameExt(pqr,'.%d'%sdie))
    if os.path.exists('%s.dx'%atompot0) : do0 = False
    # Boundary condition template for SDIE != PDIE
    inputstring = bctemplate.replace("PQRNAME",os.path.basename(pqr))
    # template for SDIE != PDIE
    inputstring += template.replace("CDCOORD","%8.3f %8.3f %8.3f"%tuple(midpoint))
    inputstring = inputstring.replace("SDIE","%d"%sdie)
    inputstring = inputstring.replace("PDIE","%d"%pdie)
    inputstring = inputstring.replace("OUTPUT","%s"%atompot0)
    inputstring += "\nquit\n"
    writeLines(inputname0,inputstring)
    # Configure file names for SDIE == PDIE
    inputname1 = nameExt(pqr,'.%d.in'%pdie)
    epot = os.path.basename(nameExt(pqr,'.%d_pot'%pdie))
    atompot1 = os.path.basename(nameExt(pqr,'.%d'%pdie))
    if os.path.exists('%s.dx'%atompot1) : do1 = False
    # Boundary condition template for SDIE != PDIE
    inputstring = bctemplate.replace("PQRNAME",os.path.basename(pqr))
    # template for SDIE != PDIE
    inputstring += template.replace("CDCOORD","%8.3f %8.3f %8.3f"%tuple(midpoint))
    inputstring = inputstring.replace("SDIE","%d"%pdie)
    inputstring = inputstring.replace("PDIE","%d"%pdie)
    inputstring = inputstring.replace("OUTPUT","%s"%atompot1)
    inputstring += "\nquit\n"
    writeLines(inputname1,inputstring)
    return (inputname0,inputname1),(do0,do1),bondlength,nDUM,("%s.dx"%atompot0,"%s.dx"%atompot1)

def extractField(outputs,bondlength,nDUM) :
    outname = outputs[0].replace("78_atompot.txt","field-err.xvg")
    for output in outputs :
        lines = readFile(output)[4:nDUM+4]
        pot = []
        for line in lines :
            pot.append(float(line))
        m1 = int(numpy.ceil((nDUM-1)/2-1))
        m2 = int(numpy.ceil((nDUM-1)/2))
        m3 = int(numpy.floor((nDUM-1)/2))
        m4 = int(numpy.floor((nDUM-1)/2+1))
        if "78_atompot" in output :
            field78 = (pot[m1]-pot[m2])/(bondlength/(nDUM-1))/2 + \
                (pot[m3]-pot[m4])/(bondlength/(nDUM-1))/2
            outname = output.replace("78_atompot.txt","f78-f1-rf.xvg")
        elif "1_atompot" in output :
            field1 = (pot[m1]-pot[m2])/(bondlength/(nDUM-1))/2 + \
                      (pot[m3]-pot[m4])/(bondlength/(nDUM-1))/2
    fields = [ field78, field1, field78-field1 ]
    string = "%20.6e %20.6e %20.6e\n"%(field78,field1,field78-field1)
    writeLines(outname,string)
    return fields

def nameExt( filename, ext ) :
    if ext.startswith('.') :
        ext = ext[1:]
    filename = os.path.abspath(filename)
    if filename[:-len(ext.split()[-1])] != '.%s'%ext and '.' in filename :
        return filename.replace(filename.split('.')[-1],ext)
    else :
        return filename+'.%s'%ext
    return filename

def makeAPBScoords( trr, tpr, top, ndx, out, pdie=1, parm="/Users/ritchie/Desktop/tmp_amoeba/amoeba.prm", dat="/Users/ritchie/Utilities/apbs/AMBER.DAT", rerun = False ) :
    # Since the box is a grid, being consistent about how the protein is oriented
    # is probably a good idea.
    if ".gro" in trr or ".pdb" in trr :
        cmd = "echo 0 | editconf -f %s -o %s -center 0 0 0 -princ"%(trr,trr)
    xyzname = nameExt(out,'.xyz')
    pqrname = nameExt(out,'.pqr')
    xvgname = nameExt(out,'.coulomb.xvg')
    printbox('%s\n%s\n%s'%(xyzname,pqrname,xvgname))
    xyznames = makeXYZ( trr, tpr, top, ndx, xyzname, parm )
    pqrnames = makePQR( trr, tpr, top, ndx, pqrname, xvgname, pdie, dat )
    for pqr in pqrnames :
        input,missing,bondlength,nDUM,output = makeAPBSin(pqr,pdie)
#        for i in range(len(input)) :
#            if missing[i] or rerun :
#                cmd="apbs %s"%input[i]
#                os.system(cmd)
#        extractField(output,bondlength,nDUM)

