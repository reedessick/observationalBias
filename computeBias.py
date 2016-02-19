#!/usr/bin/python
usage="computeBias.py [--options] psdfile.pkl psdfile.pkl psdfile.pkl ..."
description = "compute the observational bias associated with the psd files"
author = "reed.essick@ligo.org"

#-------------------------------------------------

import os

import pickle

from collections import defaultdict

import healpy as hp
import numpy as np

import utils as bayesburst_utils ### may get conflicts between snr and bayesburst modules...
import detector_cache

from lal.lal import GreenwichMeanSiderealTime as GMST

from optparse import OptionParser

#-------------------------------------------------

parser = OptionParser(usage=usage, description=description)

parser.add_option("-v", "--verbose", default=False, action="store_true")

parser.add_option("-n", "--nside", default=128, type="int", help="must be a power of 2")

parser.add_option("-f", "--min-freq", default=32.0, type="float", help="the minimum frequency at which we evaluate the PSD")
parser.add_option("-F", "--max-freq", default=512.0, type="float", help="the maximum frequency at which we evaluate the PSD")
parser.add_option("-d", "--dfreq", default=10.0, type="float", help="the frequency spacing for averaging the PSDs")

parser.add_option("-o", "--output-dir", default=".", type="string")
parser.add_option("-t", "--tag", default="", type="string")

parser.add_option("", "--saveIndividualMaps", default=False, action="store_true")
parser.add_option("", "--alsoEarthFixed", default=False, action="store_true")

opts, args = parser.parse_args()

if opts.tag:
    opts.tag = "_"+opts.tag

if not os.path.exists(opts.output_dir):
    os.makedirs(opts.output_dir)

#-------------------------------------------------

### sort  psdfiles by timestamps
if opts.verbose:
    print "sorting %d psdfiles"%(len(args))
psdfiles = defaultdict( list )
ifos = set()
for psdfilename in args:
    filename = os.path.basename(psdfilename).strip(".pkl").split("-")
    s, d = (int(_) for _ in filename[-2:])
    psdfiles[(s, s+d)].append( psdfilename )
    ifos.add( filename[0] )
segs = psdfiles.keys()
segs.sort(key=lambda l:l[0]) ### assume non-overlapping segments...
ifos = sorted(ifos)

#-------------------------------------------------

### instantiate eigval maps
npix = hp.nside2npix(opts.nside)
tracemap = np.zeros((npix,), dtype=float) ### array storing the trace of the individual eigenvalues
if opts.alsoEarthFixed:
    eftracemap = np.zeros((npix,), dtype=float)

### the cooridnates corresponding to our pixels
thetas, ras = hp.pix2ang( opts.nside, np.arange(npix) )

### set up frequencies at which we sample the eigenvalue maps
smpl_freqs = np.arange( opts.min_freq, opts.max_freq+opts.dfreq, opts.dfreq)
nsmpl = len(smpl_freqs)
freqWeights = smpl_freqs**(-7./3)

### instantiate detector and network objects
detectors = dict( (ifo, detector_cache.detectors[ifo]) for ifo in ifos )
network = bayesburst_utils.Network( detectors=detectors.values(), freqs=smpl_freqs )

#-------------------------------------------------

### iterate through segments
if opts.verbose:
    print "processing : "
for ind, seg in enumerate(segs):
    if opts.verbose:
        print "    %d - %d"%seg
    gps = 0.5*np.sum(seg)
    gmst = GMST( gps )

    ### update PSD objects
    for filename in psdfiles[seg]:
        ifo = os.path.basename(filename).split("-")[0]
        if opts.verbose:
            print "      %s -> %s"%(ifo, filename)
        file_obj = open(filename, "r")
        psd = pickle.load(file_obj)
        file_obj.close()
        detectors[ifo].set_psd( psd.get_psd(), freqs=psd.get_freqs() )

    ### compute phis from gps
    phis = (ras - gmst)%(2*np.pi)

    ### compute network sensitivity and weight relative parts of the sky accordingly
    if opts.verbose:
        print "      computing tracemap in celestial coordinates"
    A = network.A( thetas, phis, 0.0, no_psd=False ) ### compute sensitivity matrix
    trace = A[:,:,0,0] + A[:,:,1,1] ### take the trace of A ~ sensitivity
    trace = np.sum( trace * freqWeights, axis=1 ) ### weight relative freuqencies and sum over them
    trace = trace**(1.5) ### weight appropriately for p(d) ~ D^2 dD
    trace /= np.sum(trace) ### normalize so this is a probability distribution        

    if opts.saveIndividualMaps:
        outfits = "%s/celestial-freqWeightTrace%s-%d-%d.fits"%(opts.output_dir, opts.tag, seg[0], seg[1]-seg[0])
        if opts.verbose:
            print "        writing : %s"%(outfits)
        hp.write_map(outfits, trace)

    ### compute distance weight
    invDsqr = 0.0
    for detector in detectors.values():
        invDsqr += np.sum( freqWeights / detector.psd.interpolate( smpl_freqs ) )
    distWeight = invDsqr**(-1.5) ### weight appropriately for p(d) ~ D^2 dD

    ### add to the stack
    tracemap += distWeight * trace

    ### sanity check by plotting in Earth-Fixed coords
    if opts.alsoEarthFixed:
        if opts.verbose:
            print "      computing tracemap in Earth-Fixed coords"
        phis = ras ### plot in earth-fixed so this is independent of time
        A = network.A( thetas, phis, 0.0, no_psd=False ) ### compute sensitivity matrix
        trace = A[:,:,0,0] + A[:,:,1,1] ### take the trace of A ~ sensitivity
        trace = np.sum( trace * freqWeights, axis=1 ) ### weight relative freuqencies and sum over them
        trace = trace**(1.5) ### weight appropriately for p(d) ~ D^2 dD
        trace /= np.sum(trace) ### normalize so this is a probability distribution        

        if opts.saveIndividualMaps:
            outfits = "%s/earthFixed-freqWeightTrace%s-%d-%d.fits"%(opts.output_dir, opts.tag, seg[0], seg[1]-seg[0])
            if opts.verbose:
                print "        writing : %s"%(outfits)
            hp.write_map(outfits, trace)

        eftracemap += distWeight * trace

### normalize map
tracemap /= np.sum(tracemap)

### write output
outfits = "%s/celestial-freqDistWeightTrace%s-%d-%d.fits"%(opts.output_dir, opts.tag, segs[0][0], segs[-1][1]-segs[0][0])
if opts.verbose:
    print "writing : %s"%(outfits)
hp.write_map(outfits, tracemap)

if opts.alsoEarthFixed:
    eftracemap /= np.sum(eftracemap)

    outfits = "%s/earthFixed-freqDistWeightTrace%s-%d-%d.fits"%(opts.output_dir, opts.tag, segs[0][0], segs[-1][1]-segs[0][0])
    if opts.verbose:
        print "writing : %s"%(outfits)
    hp.write_map(outfits, tracemap)
