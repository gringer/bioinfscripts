#!/usr/bin/env python

'''
reads event data from ONT fast5 file, writes event data matrix to output.

Copyright 2016, David Eccles (gringer) <bioinformatics@gringene.org>

Permission to use, copy, modify, and/or distribute this software for
any purpose with or without fee is hereby granted. The software is
provided "as is" and the author disclaims all warranties with regard
to this software including all implied warranties of merchantability
and fitness. The parties responsible for running the code are solely
liable for the consequences of code excecution.
'''

import sys
import h5py
import numpy

def generate_event_matrix(fileName):
    '''write out event matrix from fast5, return False if not present'''
    try:
        h5File = h5py.File(fileName, 'r')
        h5File.close()
    except:
        return False
    with h5py.File(fileName, 'r') as h5File:
      runMeta = h5File['UniqueGlobalKey/tracking_id'].attrs
      runID = '%s_%s' % (runMeta["device_id"],runMeta["run_id"][0:16])
      eventBase = "/Analyses/EventDetection_000/Reads/"
      readNames = h5File[eventBase]
      for readName in readNames:
        eventLocation = "/Analyses/EventDetection_000/Reads/%s/Events" % readName
        outData = h5File[eventLocation]
        headers = outData.dtype
        sys.stdout.write(",".join(headers.names)+"\n")
        # There *has* to be an easier way to do this while preserving
        # precision. Reading element by element seems very inefficient
        for line in outData:
          res=map(str,line)
          sys.stdout.write(readName + "," + ",".join(res) + "\n")

if len(sys.argv) < 2:
    sys.stderr.write('Usage: sys.argv[0] <fast5 file name>')
    sys.exit(1)

generate_event_matrix(sys.argv[1]);
