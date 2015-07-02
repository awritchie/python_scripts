#! /usr/bin/env python

import sys
import numpy as np

try :
    coords = sys.argv[1]
    force = sys.argv[2]
    out = sys.argv[3]
except :
    print "Usage: %s <CN_coords.xvg> <TCHG_force.xvg> <outname.xvg>"%sys.argv[0]
    sys.exit()

def project( a, b ) :
    return np.dot(a,b)/np.sqrt(np.dot(b,b))

cnlines = np.loadtxt(coords,"float64")
tclines = np.loadtxt(force,"float64")

maxiter = len(tclines)
fields = np.zeros(maxiter,"float64")
for i in range(maxiter) :
    cd = cnlines[i][1:4]
    ne = cnlines[i][4:7]
    tc = tclines[i][1:]
    
    cn = ne - cd
    # Units are kJ/(mol nm e-)
    # 1 kJ/(mol nm e-) = 0.04012 kbT/(A e-)
    fields[i] = project(tc,cn)*0.04012

np.savetxt(out, fields)