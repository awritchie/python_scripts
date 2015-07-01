#! /usr/bin/env python

import os
import sys
import h5py as h5
import numpy as np
import itertools

def closest_value( array, cutoff, minimum = 1 ) :
    avg = np.average(array)
    stdev = np.std(array)
    length = len(array)
    for i in range(min(length,minimum),length + 1) :
        for each in list(itertools.combinations(array,i)) :
            if ( abs(np.average(each) - avg) <= cutoff ) :
                index = [ array.index(j) for j in each ]
                return index, np.average(each)

def cluster( h5name, ref_group, cutoff = 0.1, min_frames_per_bin = 3 ) :
    plural = "s"
    if min_frames_per_bin == 1 : plural = ""
    print >> sys.stderr,"*"*75
    print >> sys.stderr,"Reading %s"%h5name
    print >> sys.stderr,"Finding subarrays with averages within +/- %.6f of bin averages and at\n least %d frame%s"%(cutoff,min_frames_per_bin,plural)
    try :
        wham = h5.File(h5name,"r")
    except :
        print >> sys.stderr,"Error opening %s\n"%h5name
        sys.exit()
    try :
        probdat = wham['/Ensemble/Probability']
        chi1, chi2, prob = probdat[:,0], probdat[:,1], probdat[:,2]
    except :
        print >> sys.stderr,"Error reading /Ensemble/Probability\n"
        sys.exit()
    try :
        full_avg = wham['/Ensemble/%s'%ref_group][-1][0]
    except :
        print >> sys.stderr,"Error reading /Ensemble/%s\n"%(ref_group)
        sys.exit()
    try :
        bincounts = wham['/Ensemble/BinCounts'][:]
        nbins = len(bincounts)
    except :
        print >> sys.stderr,"Error reading /Ensemble/BinCounts\n"
        sys.exit()
    try :
        b, e = wham['/Ensemble'].attrs.get('First frame read, last frame read')
    except :
        print >> sys.stderr,"Error finding the span of frames included in average\n"
        sys.exit()
    try :
        nexp = len(wham['Trajectories'].values())
    except :
        print >> sys.stderr,"Error finding the number of trajectories\n"
        sys.exit()
    ref_frames = []
    bin_frames = []
    for i in range(nexp) :
        item = "/Trajectories/Traj-%d/%s"%(i,ref_group)
        bin_item = "/Trajectories/Traj-%d/Bin"%(i)
        try :
            ref_frames.append(wham[item][:])
        except :
            print >> sys.stderr,"Error opening %s in %s\n"%(item,h5name)
            sys.exit()
        try :
            bin_frames.append(wham[bin_item][:])
        except :
            print >> sys.stderr,"Error opening %s in %s\n"(bin_item,h5name)

    # Prepare the bin_array as a 2D list
    bin_array = []
    for i in range(nbins) :
        bin_array.append([])

    # Map each frame to the bin number; i=experiment, j=frame number within experiment i
    frames_to_keep = []
    for i in range(nexp) :
        frames_to_keep.append([])
        for j in range( len(ref_frames[i]) ) :
            if ( j >= b and j<= e ):
                bin_ij = bin_frames[i][j]
                bin_array[bin_ij].append([i,j])

    # Find the average reference value for each bin
    total = 0
    total_frames = 0
    total_keep = 0
    kept_size = []
    for i in range(nbins) :
        if len(bin_array[i]) > 0 :
            dat_array = []
            for j in bin_array[i] :
                m,n = j
                dat_array.append(ref_frames[m][n])
            index, sub_array_avg = closest_value(dat_array, cutoff, minimum = min_frames_per_bin )
            kept_size.append(len(index))
            total_keep += len(index)
            total_frames += len(dat_array)
            total += sub_array_avg * prob[i]
            for j in index :
                m,n = bin_array[i][j]
                frames_to_keep[m].append(n)

    # Print some statistics information to stderr so that the output of this
    # can be piped using bash
    print >> sys.stderr,"Keeping %d out of %d frames (%.3f%%)"%(total_keep, total_frames, 100. * total_keep / total_frames )
    print >> sys.stderr,"Average number of frames per bin: %.3f +/- %.3f"%(np.average(kept_size),np.std(kept_size))
    print >> sys.stderr,"Maximum number of frames in a bin: %d"%(max(kept_size))
    print >> sys.stderr,"%-25s: %12.6f"%("Clustered Average",total)
    print >> sys.stderr,"%-25s: %12.6f"%("Full Average",full_avg)
    print >> sys.stderr, "%-25s: %12.6f"%("dAverage",abs(total-full_avg))
    print >> sys.stderr,"*"*75

    for i in range(nexp) :
        nframes = len(frames_to_keep[i])
        outname = Andrew_Save_frame_subarray(h5name, wham, ref_group, i, 2,cutoff)
        frame_index = open(outname,"w")
        if nframes > 0 :
            sorted_keep = sorted(frames_to_keep[i])
            for j in range(nframes) :
                frame_index.write("%d\n"%sorted_keep[j])
        frame_index.close()

def Andrew_Save_frame_subarray( h5name, opened_h5file, ref_group, trajn, ndof, cutoff ) :
    """
        This is incredibly specific to my personal simulation directory heirarchy.
    """
    outpath = os.path.dirname(os.path.abspath(h5name))
    for dof in range(ndof) :
        outpath += "/%d" %(int(opened_h5file['/Trajectories/Traj-%d/BiasingCoordinate/DoF-%d'%(trajn,dof)].attrs.get('phi0 (degrees), dphi (degrees), force constant (kJ/mol/degree^2), temperature (K)')[0]))
    outpath += "/clustered_%s.%.0e.fndx"%(ref_group,cutoff)
    return outpath

try :
    h5name = sys.argv[1]
    ref_group = sys.argv[2]
except :
    print "Usage: <H5 File Name> <Reference Group>"
    sys.exit()

cluster( h5name, ref_group, cutoff = 1, min_frames_per_bin = 1 )