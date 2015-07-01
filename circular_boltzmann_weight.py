#!/usr/bin/env python

#nbins = 72

from numpy import array, sqrt, zeros, loadtxt, pi, sin, cos, sum
from math import atan2, isnan
import sys

def rad( x ) :
    return x*pi/180

def deg( x ) :
    return x*180/pi

# get_averages returns the weighted cartesian vector sum and weighted sum of
# the sins and cos of the angles.  Everything is converted to radians.  The
# single probability of the bin being looked at is entered.  By summing all
# of the xs, all of the ys, all of the sumcos**2, and all of the sumsin**2,
# the the weighted mean angle and weighted variance is produced. 
def get_averages( angles, probability ) :
    vectors = [ array(( cos( rad(angle) ), sin( rad(angle) ) )) for
                angle in angles ]
    xs, ys, sumsin, sumcos = 0, 0, 0, 0
    sample_size = len( vectors )
    for i in range( sample_size ) :
        if not isnan(vectors[i][0])  and not isnan(vectors[i][1]) :
            x, y = float( vectors[i][0] ), float( vectors[i][1] )
            xs += x * probability/sample_size
            ys += y * probability/sample_size
            sumsin += y * probability/sample_size
            sumcos += x * probability/sample_size
    R = sqrt( sumcos**2 + sumsin**2 )
    average = deg( atan2( ys, xs ) )
    variance = deg( 1-R )
    return xs, ys, sumcos, sumsin    

try :
    probfile = sys.argv[1]
    datafile = sys.argv[2]
    probabilities = loadtxt( probfile, "float64" )
    nbins=len(probabilities)
    if sum(probabilities) > 1 or sum(probabilities) < 0.999 :
        print "\nWARNING\n>>>>Probabilities sum to %.30f<<<<\n" %sum(probabilities)
except:
    print "usage: ./script <prob file> <data file>"
    print "The probability file is the one-dimensional list of probabilities."
    print "The data file lists files:"
    print " <bin file> <data file>"
    print "which will be used for additional weighting of the result."
    sys.exit()

# the print i, counts[i],"\n"  of the degree of freedom in bin i which was sampled ncounts times is
# prob(i) * ( d[1] + d[2] + ... ) / ncounts

file = open( datafile )
datafilelines = file.readlines()
file.close()


binfile0, dofdata0 = datafilelines[0].split()
try :
    doffile0 = open( dofdata0 )
    dof0=doffile0.readlines()[0]
    values = len(dof0.split())
except IOError :
    print "file '%s' or '%s' not read, skipping" % ( binfile, dofdata )
    values = 1
grammar = "s" if values > 1 else ""
print "Looking at %i column%s of data" %(values,grammar)

for value in range(values) :
    counts = [ [] for i in range(nbins) ]
    nread = 0
    for line in datafilelines :
        binfile, dofdata = line.split()
        print "reading", dofdata, "using bins", binfile, ";",

        try : 
            binsfile = open( binfile )
            bins = binsfile.readlines()
            doffile = open( dofdata )
            dof=doffile.readlines()

            nread += 1
        except IOError :
            print "file '%s' or '%s' not read, skipping" % ( binfile, dofdata )
            continue

        lenbins = len(bins)
        lendof = len(dof)
        if lenbins != lendof :
            print "WARNING! %d points in bin file, %d points in dof file" % ( lenbins, lendof )
        else :
            print "%d points" % lenbins
        
        for i in range(len(bins)) :
            try :
                thisbin = int( bins[i].strip() )
                thisdof = float( dof[i].split()[value] )
                if isnan(thisdof) :
                    print "NaN Error in file %s, line %i!" %(dofdata,i+1)
            except ValueError :
                print "Incomprehensible '", dof[i].split()[value], "' at point", i
                continue
            except IndexError:
                print "IndexError in dof file %s at line %d" % ( dofdata, i )
                print thisdof
                continue
    
            try:
                counts[thisbin].append( thisdof )
            except IndexError :
                print "IndexError"
                counts[0].append( thisdof ) # close

            
    vartotal = 0
    xs, ys, sumcos, sumsin = 0, 0, 0, 0
    average = []
    for i in range( nbins ) :
        meansq = 0
#       print "bin", i, "counts", len(counts[i])
        if len(counts[i]) > 0 :
            angles = array( counts[i] )
            x, y, sc, ss = get_averages( angles, probabilities[i] )
            xs += x
            ys += y
            sumcos += sc
            sumsin += ss
    mean = deg( atan2( ys, xs ) )
    R = sqrt( sumcos**2 + sumsin**2 )
    variance = deg( 1-R )

    print value,":",probfile, nread, mean, variance 
