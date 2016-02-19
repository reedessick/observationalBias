#!/usr/bin/python

usage="computePSDs.py [--options] start end"
description="queries data and computes PSDs"
author = "reed.essick@ligo.org"

#-------------------------------------------------

import os
import sys

import snr_utils as psd_utils

from glue.ligolw import ligolw
from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw import utils as ligolw_utils

from laldetchar.idq import event

import subprocess as sp

import pickle

from optparse import OptionParser 

#-------------------------------------------------

parser = OptionParser(usage=usage, description=description)

parser.add_option("-v", "--verbose", default=False, action="store_true")

parser.add_option("-o", "--observatory", default=[], action="append", type="string")

parser.add_option("-c", "--channel", default="GDS-CALIB_STRAIN", type="string")
parser.add_option("-F", "--frame-type", default="HOFT_C00", type="string")

parser.add_option("-f", "--flag", default="DMT-ANALYSIS_READY:1", type="string")

parser.add_option("-u", "--segdb-url", default="https://segments.ligo.org", type="string")

parser.add_option("-O", "--output-dir", default=".", type="string")

parser.add_option("-d", "--psd-dur", default=64, type="int")
parser.add_option("-n", "--num-segs", default=16, type="int")

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

### iterate through segments and estimte PSD for each IFO
if opts.verbose:
    print "  processing:"
for s, e in segs:
    while s < e:
        _e = min(e, s+opts.psd_dur)
        if opts.verbose:
            print "    %d - %d"%(s, _e)

        dirtmplate = "%s/"%(opts.output_dir) + "%s" + "1-%s_%s-%d"%(opts.frame_type.replace("-","_"), opts.channel.replace("-","_"), s/100000)
        tmplate = "%s/%s" + "1-%s_%s-%d-%d.pkl"%(opts.frame_type.replace("-","_"), opts.channel.replace("-","_"), s, _e-s)

        for ifo in opts.observatory:
            chan = "%s1:%s"%(ifo, opts.channel)

            cmd = "gw_data_find -o %s --type %s1_%s -u file -s %d -e %d"%(ifo, ifo, opts.frame_type, s, _e)
            if opts.verbose:
                print "      %s"%(cmd)
            frames = [_.strip().replace("file://localhost","") for _ in sp.Popen( cmd.split(), stdout=sp.PIPE, stderr=sp.PIPE ).communicate()[0].split()]
            if len(frames):
                if opts.verbose:
                    for frame in frames:
                        print "        %s"%(frame)
                    print "      computing PSD"
                psd = psd_utils.frames2PSD( frames, chan, start=s, stop=_e, num_segs=opts.num_segs, overlap=0.0 )

                outdir = dirtmplate%(ifo)
                if not os.path.exists(outdir):
                    os.makedirs( outdir )

                psdfilename = tmplate%(outdir, ifo)
                if opts.verbose:
                    print "      writing : %s"%(psdfilename)
                file_obj = open(psdfilename, "w")
                pickle.dump(psd, file_obj)
                file_obj.close()
            else:
                print >> sys.stderr, "no frames found for %d - %d"%(s, _e)

        s += opts.psd_dur
