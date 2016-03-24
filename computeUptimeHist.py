#!/usr/bin/python

usage="computeUptimeHist.py [--options] start end"
description="queries data and computes histograms based on uptime"
author = "reed.essick@ligo.org"

#-------------------------------------------------

import os
import sys

import numpy as np

#import snr_utils as psd_utils

from glue.ligolw import ligolw
from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw import utils as ligolw_utils

from laldetchar.idq import event

import subprocess as sp

import pickle

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

from optparse import OptionParser 

#-------------------------------------------------

parser = OptionParser(usage=usage, description=description)

parser.add_option("-v", "--verbose", default=False, action="store_true")

parser.add_option("-o", "--observatory", default=[], action="append", type="string")

#parser.add_option("-c", "--channel", default="GDS-CALIB_STRAIN", type="string")
#parser.add_option("-F", "--frame-type", default="HOFT_C00", type="string")

parser.add_option("-f", "--flag", default="DMT-ANALYSIS_READY:1", type="string")

parser.add_option("-u", "--segdb-url", default="https://segments.ligo.org", type="string")

parser.add_option("-O", "--output-dir", default=".", type="string")

#parser.add_option("-d", "--psd-dur", default=64, type="int")
#parser.add_option("-n", "--num-segs", default=16, type="int")

opts, args = parser.parse_args()

if len(args) != 2:
    raise ValueError("must supply two arguments : %s"%(usage))
start, end = [float(_) for _ in args]

if not os.path.exists(opts.output_dir):
    os.makedirs(opts.output_dir)

opts.observatory.sort()

#-------------------------------------------------

### query segments
segs = [[start, end]]
for ifo in opts.observatory:

    segfilename = "%s/%s1_%s-%d-%d.xml"%(opts.output_dir, ifo, opts.flag.replace(":","_"), start, end-start)
    cmd = "ligolw_segment_query_dqsegdb -t %s -q -a %s1:%s -s %d -e %d -o %s"%(opts.segdb_url, ifo, opts.flag, start, end, segfilename)
    if opts.verbose:
        print "querying %s segments for %s\n    %s"%(opts.flag, ifo, cmd)

    output = sp.Popen( cmd.split(), stdout=sp.PIPE, stderr=sp.PIPE ).communicate()

    ### iterate over segments, computing PSDs
    xmldoc = ligolw_utils.load_filename(segfilename, contenthandler=lsctables.use_in(ligolw.LIGOLWContentHandler))
    segs = event.andsegments( [segs, [[row.start_time, row.end_time] for row in table.get_table(xmldoc, lsctables.SegmentTable.tableName)]] ) ### take intersection of segments

### write intersection to ascii file
segfilename = "%s/intersection-%d-%d.seg"%(opts.output_dir, start, end-start)
if opts.verbose:
    print "found %d sec of joint livetime"%(event.livetime(segs))
    print "writing : %s"%(segfilename)
file_obj = open(segfilename, 'w')
for s, e in segs:
    print >> file_obj, s, e
file_obj.close()

#-------------------------------------------------

### build histograms
### will have to do this by hand (unfortunately)

### set up bins
day = 86400.
Nbins = 60*12 ### one bin per minute
bins = np.linspace(0, day, Nbins+1)
counts = np.zeros((Nbins,), dtype=float)
binDur = bins[1]-bins[0]

### determine reference start time
tref = 1125964817 ###Fri Sep 11 00:00:00 GMT 2015

if opts.verbose:
    print "computing binning"

for s, e in segs:
    if opts.verbose:
        print "processing : %d - %d"%(s, e)
    relative_start = (s - tref)%(day) ### where the segment starts in the day
    dur = e - s

    for i in xrange(Nbins): ### find the bin we start in
        if bins[i+1] > relative_start:
            break
    else:
        raise ValueError("could not find relative start's position in binning")

    seg = bins[i+1]-relative_start ### add the initial bit to the bin we start in
    if seg > dur:
        seg = dur
    counts[i] += seg
    dur -= seg
    i = (i+1)%(Nbins)
    
    while dur > 0: ### iterate until we run out of this segment, adding time to each bin as needed.
        if dur < binDur:
            counts[i] += dur
            break
        else:
            counts[i] += binDur
            dur -= binDur
        i = (i+1)%(Nbins)

filename = "%s/UptimeHist-%s-%d-%d.pkl"%(opts.output_dir, opts.flag.replace(":","_"), start, end-start)
if opts.verbose:
    print filename
file_obj = open(filename, "w")
pickle.dump(bins, file_obj)
pickle.dump(counts, file_obj)
file_obj.close()

counts /= np.sum(counts) ### compute fractional occupation

counts /= binDur/day ### normalize by the expected amount

### plot histogram
if opts.verbose:
    print "plotting the histogram"
fig = plt.figure()
ax = fig.gca()

x = []
y = []
Y = []
for i in xrange(Nbins):
    x += list(bins[i:i+2])
    y += [0]*2
    Y += [counts[i]]*2
ax.fill_between( x, y, Y )

ax.plot([0, 86400], [1, 1], 'k:', alpha=0.5)

ax.set_title('%d sec joint livetime out of %d sec'%(event.livetime(segs), end-start))

ax.set_xlabel('hours relative to Sept 11 00:00:00 GMT 2015')
ax.set_ylabel('fraction of joint livetime / uniform distribution')

ax.xaxis.set_ticks([i*3600 for i in xrange(25)])
ax.xaxis.set_ticklabels(["%d"%i for i in xrange(25)])

ax.set_xlim(xmin=0, xmax=86400)

figname = "%s/UptimeHist-%s-%d-%d.png"%(opts.output_dir, opts.flag.replace(":","_"), start, end-start)
if opts.verbose:
    print "saving : %s"%(figname)
fig.savefig( figname )
plt.close( fig )
