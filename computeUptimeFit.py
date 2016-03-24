#!/usr/bin/python
usage="computeUptimeFit.py [--options] file.pkl file.pkl file.pkl"
description="reads in the bins and counts from pickled files and fits them, plotting the overlays of raw histograms, normalized histograms, and fitted curves"
author = "reed.essick@ligo.org"

import pickle

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

import numpy as np
import scipy

from optparse import OptionParser

#-------------------------------------------------

parser=OptionParser(usage=usage, description=description)

parser.add_option("-v", "--verbose", default=False, action="store_true")

parser.add_option("-o", "--output-dir", default=".", type="string")
parser.add_option("-t", "--tag", default="", type="string")

opts, args = parser.parse_args()

if opts.tag:
    opts.tag = "_"+opts.tag

#-------------------------------------------------

rfig = plt.figure()
rax = rfig.gca()

nfig = plt.figure()
nax = nfig.gca()

X = np.linspace(0, 86400., 1001)
Y = np.zeros_like(X)
S = np.zeros_like(X)

### overlay raw counts and normalized counts
for pkl in args:
    if opts.verbose:
        print "reading : %s"%pkl
    file_obj = open(pkl, "r")
    bins = pickle.load(file_obj)
    counts = pickle.load(file_obj)
    file_obj.close()

    x = 0.5*(bins[1:]+bins[:-1])/3600.
    rax.plot( x, counts )

    counts /= np.sum(counts)
    nax.plot( x, counts )

    Y += np.mean(counts)
    S += np.var(counts)

Y /= len(args)
S /= len(args)

nax.plot( X/3600., Y, 'k--', linewidth=2 )
nax.plot( X/3600., Y+(2*S)**0.5, 'k--', linewidth=2 )
nax.plot( X/3600., Y-(2*S)**0.5, 'k--', linewidth=2 )

dphi = np.pi/8
phase = -dphi
for i in xrange(6):
    nax.plot( X/3600., Y+(2*S)**0.5*np.sin(2*np.pi*X/86400. - phase), 'k-', linewidth=2 )
    phase += dphi

rax.set_xlabel('hours')
nax.set_xlabel('hours')

rax.set_ylabel('seconds')
nax.set_ylabel('fraction of livetime')

rax.set_xlim(xmin=0, xmax=24)
nax.set_xlim(xmin=0, xmax=24)

rax.grid(True, which="both")
nax.grid(True, which="both")


figname = "%s/rawCounts%s.png"%(opts.output_dir, opts.tag)
if opts.verbose:
    print "saving : %s"%figname
rfig.savefig(figname)
plt.close(rfig)

figname = "%s/normCounts%s.png"%(opts.output_dir, opts.tag)
if opts.verbose:
    print "saving : %s"%figname
nfig.savefig(figname)
plt.close(nfig)
