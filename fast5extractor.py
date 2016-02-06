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

import os
import sys
import h5py
import numpy

def generate_event_matrix(fileName, header=True):
    '''write out event matrix from fast5, return False if not present'''
    try:
        h5File = h5py.File(fileName, 'r')
        h5File.close()
    except:
        return False
    with h5py.File(fileName, 'r') as h5File:
      runMeta = h5File['UniqueGlobalKey/tracking_id'].attrs
      channelMeta = h5File['UniqueGlobalKey/channel_id'].attrs
      runID = '%s_%s' % (runMeta["device_id"],runMeta["run_id"][0:16])
      eventBase = "/Analyses/EventDetection_000/Reads/"
      readNames = h5File[eventBase]
      for readName in readNames:
        readMetaLocation = "/Analyses/EventDetection_000/Reads/%s" % readName
        eventLocation = "/Analyses/EventDetection_000/Reads/%s/Events" % readName
        outMeta = h5File[readMetaLocation].attrs
        channel = str(channelMeta["channel_number"])
        mux = str(outMeta["start_mux"])
        headers = h5File[eventLocation].dtype
        outData = h5File[eventLocation][()] # load entire array into memory
        if(header):
            sys.stdout.write("runID,channel,mux,read,"+",".join(headers.names)+"\n")
        # There *has* to be an easier way to do this while preserving
        # precision. Reading element by element seems very inefficient
        for line in outData:
          res=map(str,line)
          # data seems to be normalised, but just in case it isn't, here's the formula for
          # future reference: pA = (raw + offset)*range/digitisation
          # (using channelMeta[("offset", "range", "digitisation")])
          sys.stdout.write(",".join((runID,channel,mux,readName)) + "," + ",".join(res) + "\n")

if len(sys.argv) < 2:
    sys.stderr.write('Usage: sys.argv[0] <fast5 file name>')
    sys.exit(1)

fileArg = sys.argv[1]
seenHeader = False

if(os.path.isdir(fileArg)):
    sys.stderr.write("Processing directory '%s':\n" % fileArg)
    for dirPath, dirNames, fileNames in os.walk(fileArg):
        fc = len(fileNames)
        for fileName in fileNames:
            if(fileName.endswith(".fast5")): # only process fast5 files
                sys.stderr.write("  Processing file '%s'..." % fileName)
                generate_event_matrix(os.path.join(dirPath, fileName), not seenHeader)
                fc -= 1
                seenHeader = True
                if(fc == 1):
                    sys.stderr.write(" done (%d more file to process)\n" % fc)
                else:
                    sys.stderr.write(" done (%d more files to process)\n" % fc)
elif(os.path.isfile(fileArg)):
    generate_event_matrix(fileArg)
