#!/usr/bin/env python

'''
reads elements from ONT fast5 file and writes them to standard output.

Copyright 2016, David Eccles (gringer) <bioinformatics@gringene.org>

Permission to use, copy, modify, and distribute this software for any
purpose with or without fee is hereby granted. The software is
provided "as is" and the author disclaims all warranties with regard
to this software including all implied warranties of merchantability
and fitness. In other words, the parties responsible for running the
code are solely liable for the consequences of code execution.
'''

import os
import sys
import h5py
import numpy
from collections import deque, Counter, OrderedDict
from itertools import islice
from bisect import insort, bisect_left
from struct import pack
from array import array

def generate_consensus_matrix(fileName, header=True):
    '''write out 2D consensus matrix from fast5, return False if not present'''
    try:
        h5File = h5py.File(fileName, 'r')
        h5File.close()
    except:
        return False
    with h5py.File(fileName, 'r') as h5File:
        runMeta = h5File['UniqueGlobalKey/tracking_id'].attrs
        channelMeta = h5File['UniqueGlobalKey/channel_id'].attrs
        runID = '%s_%s' % (runMeta["device_id"],runMeta["run_id"][0:16])
        eventBaseTemp = "/Analyses/Basecall_1D_000/BaseCalled_template/Events/"
        eventBaseComp = "/Analyses/Basecall_1D_000/BaseCalled_complement/Events/"
        evtTemp = h5File[eventBaseTemp][()]
        evtComp = h5File[eventBaseComp][()]
        channelRate = channelMeta["sampling_rate"]
        evtTempStart = map(int,evtTemp["start"] * channelRate)
        evtTempLen = map(int,evtTemp["length"] * channelRate)
        evtCompStart = map(int,evtComp["start"] * channelRate)
        evtCompLen = map(int,evtComp["length"] * channelRate)
        tempRawStart = evtTempStart[0]
        compRawStart = evtCompStart[0]
        alignmentBase = "/Analyses/Basecall_2D_000/BaseCalled_2D/Alignment/"
        if(not alignmentBase in h5File):
            return False
        readName = ""
        mux = ""
        rawReadBase = "/Raw/Reads/"
        for tReadName in h5File[rawReadBase]:
            readName = tReadName
            readMeta = h5File['%s%s' % (rawReadBase, readName)].attrs
            mux = str(readMeta["start_mux"])
            channel = str(channelMeta["channel_number"])
        alnHeaders = h5File[alignmentBase].dtype
        outAlnData = h5File[alignmentBase][()] # load entire array into memory
        if(header):
            sys.stdout.write("runID,channel,mux,read,"+
                             "tempStart,tempEnd,compStart,compEnd,bpPos," +
                             ",".join(alnHeaders.names) + "\n")
        lastkmer = ""
        bpPos = -3
        lastbpPos = 1
        tempStart = -1
        tempEnd = -1
        compStart = -1
        compEnd = -1
        for line in outAlnData:
            nextkmer = line["kmer"]
            moved = False
            if(lastkmer != nextkmer):
                if(lastkmer[1:] == nextkmer[:-1]):
                    bpPos += 1
                elif(lastkmer[2:] == nextkmer[:-2]):
                    bpPos += 2
                elif(lastkmer[3:] == nextkmer[:-3]):
                    bpPos += 3
                elif(lastkmer[4:] == nextkmer[:-4]):
                    bpPos += 4
                else:
                    bpPos += 5
                moved = True
                lastkmer = nextkmer
            if(line["template"] != -1):
                tempStart = (evtTempStart[line["template"]] if moved
                             else tempStart)
                tempEnd = (evtTempStart[line["template"]] +
                           evtTempLen[line["template"]])
            if(line["complement"] != -1):
                compStart = evtCompStart[line["complement"]]
                compEnd = ((evtCompStart[line["complement"]] +
                            evtCompLen[line["complement"]]) if moved
                           else compEnd)
            res=map(str,line)
            if(moved):
                sys.stdout.write(",".join((runID,channel,mux,readName,
                                           str(tempStart-tempRawStart),str(tempEnd-tempRawStart),
                                           str(compStart-compRawStart),str(compEnd-compRawStart),
                                           str(bpPos))) +
                                 "," + ",".join(res) + "\n")

def generate_eventdir_matrix(fileName, header=True, direction=None):
    '''write out event matrix from fast5, return False if not present'''
    try:
        h5File = h5py.File(fileName, 'r')
        h5File.close()
    except:
        return False
    with h5py.File(fileName, 'r') as h5File:
      runMeta = h5File['UniqueGlobalKey/tracking_id'].attrs
      channelMeta = h5File['UniqueGlobalKey/channel_id'].attrs
      channel = str(channelMeta["channel_number"])
      runID = '%s_%s' % (runMeta["device_id"],runMeta["run_id"][0:16])
      dir = "complement" if (direction=="r") else "template"
      eventBase = "/Analyses/EventDetection_000/Reads/"
      readNames = h5File[eventBase]
      mux = 0
      # get mux for the read
      for readName in readNames:
        readMetaLocation = "/Analyses/EventDetection_000/Reads/%s" % readName
        outMeta = h5File[readMetaLocation].attrs
        mux = str(outMeta["start_mux"])
      eventLocation = "/Analyses/Basecall_1D_000/BaseCalled_%s/Events" % (dir)
      if(not eventLocation in h5File):
          return False
      readNames = h5File[eventLocation]
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
        # - might also be useful to know start_time from outMeta["start_time"]
        #   which should be subtracted from event/start
        sys.stdout.write(",".join((runID,channel,mux,readName)) + "," + ",".join(res) + "\n")

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
          # - might also be useful to know start_time from outMeta["start_time"]
          #   which should be subtracted from event/start
          sys.stdout.write(",".join((runID,channel,mux,readName)) + "," + ",".join(res) + "\n")

def generate_fastq(fileName, callID="000"):
    '''write out fastq sequence(s) from fast5, return False if not present'''
    try:
        h5File = h5py.File(fileName, 'r')
        h5File.close()
    except:
        return False
    with h5py.File(fileName, 'r') as h5File:
        rowData = get_telemetry(h5File, callID)
        seqBase1D = "/Analyses/Basecall_1D_%s" % callID
        seqBase2D = "/Analyses/Basecall_2D_%s" % callID
        v1_2File = False
        if(not (seqBase1D in h5File) and (seqBase2D in h5File)):
            seqBase1D = seqBase2D
            v1_2File = True
        if( (rowData["templateCalledBases"] > 0) and
            (rowData["templateRawLength"] / rowData["templateCalledBases"] <= 25)):
            baseTemp = "%s/BaseCalled_template/Fastq" % seqBase1D
            sys.stdout.write("@1Dtemp_" +
                             "_".join(map(lambda x: str(rowData[x]),
                                ("runID","channel","mux","read"))) +
                                " ")
            sys.stdout.write(str(h5File[baseTemp][()][1:]))
        if( (rowData["complementCalledBases"] > 0) and
            (rowData["complementRawLength"] / rowData["complementCalledBases"] <= 25)):
            baseComp = "%s/BaseCalled_complement/Fastq" % seqBase1D
            sys.stdout.write("@1Dcomp_"+
                             "_".join(map(lambda x: str(rowData[x]),
                                ("runID","channel","mux","read"))) +
                                " ")
            sys.stdout.write(str(h5File[baseComp][()][1:]))
        if(seqBase2D in h5File):
            base2D = "%s/BaseCalled_2D/Fastq" % seqBase2D
            if((base2D in h5File)):
                sys.stdout.write("@2Dcons_"+
                             "_".join(map(lambda x: str(rowData[x]),
                                ("runID","channel","mux","read")))+
                                " ")
                sys.stdout.write(str(h5File[base2D][()][1:]))

## Running median
## See [http://code.activestate.com/recipes/578480-running-median-mean-and-mode/]
def runningMedian(seq, M):
    if(M % 2 == 0):
        sys.stderr.write("Error: median window size must be odd")
        sys.exit(1)
    seq = iter(seq)
    s = []
    m = M // 2
    s = [item for item in islice(seq,M)]
    d = deque(s)
    median = lambda : s[m] # if bool(M&1) else (s[m-1]+s[m])/2
    s.sort()
    medians = [median()] * (m+1) # set initial m samples to median of first M
    for item in seq:
        old = d.popleft()          # pop oldest from left
        d.append(item)             # push newest in from right
        del s[bisect_left(s, old)] # locate insertion point and then remove old
        insort(s, item)            # insert newest such that new sort is not required
        medians.append(median())
    medians.extend([median()] * (m))
    return medians

def get_telemetry(h5File, callID):
    runMeta = h5File['UniqueGlobalKey/tracking_id'].attrs
    channelMeta = h5File['UniqueGlobalKey/channel_id'].attrs
    useRaw = False
    rowData = OrderedDict(
        [('runID','%s_%s' % (runMeta["device_id"],runMeta["run_id"][0:16])),
         ('channel',channelMeta["channel_number"]),
         ('mux',''),('read',''),
         ('offset',channelMeta["offset"]),
         ('range',channelMeta["range"]),
         ('digitisation',channelMeta["digitisation"]),
         ('sampleRate',channelMeta["sampling_rate"]),
         ('rawStart',''),('rawLength',-1),
         ('templateRawStart',''),('templateRawLength',-1),
         ('templateCalledEvents',''),('templateCalledBases',-1),
         ('complementRawStart',''),('complementRawLength',-1),
         ('complementCalledEvents',''),('complementCalledBases',-1)
        ])
    callBase = "/Analyses/Basecall_1D_%s/Summary" % (callID)
    eventBase = "/Analyses/EventDetection_%s/Reads" % (callID)
    if(not eventBase in h5File):
        useRaw = True
        eventBase = "/Raw/Reads"
    readNames = h5File[eventBase]
    # get mux for the read
    for readName in readNames:
        readMetaLocation = "%s/%s" % (eventBase,readName)
        outMeta = h5File[readMetaLocation].attrs
        rowData["mux"] = outMeta["start_mux"]
        rowData["read"] = readName
        rowData["rawStart"] = outMeta["start_time"]
        rowData["rawLength"] = outMeta["duration"]
    for dir in ('template','complement'):
        callBase = "/Analyses/Basecall_1D_%s" % (callID)
        metaLoc = ("%s/Summary/basecall_1d_%s" % (callBase,dir) if useRaw
                   else "%s/BaseCalled_%s/Events" % (callBase,dir))
        if(metaLoc in h5File):
            dirMeta = h5File[metaLoc].attrs
            rowData["%sRawStart" % dir] = int(dirMeta["start_time"] *
                                              rowData["sampleRate"])
            rowData["%sRawLength" % dir] = int(dirMeta["duration"] *
                                               rowData["sampleRate"])
        metaLoc = ("%s/Summary/basecall_1d_%s" % (callBase,dir))
        if(metaLoc in h5File):
            dirMeta = h5File[metaLoc].attrs
            rowData["%sCalledEvents" % dir] = dirMeta["called_events"]
            rowData["%sCalledBases" % dir] = dirMeta["sequence_length"]
    return(rowData)

def generate_telemetry(fileName, callID="000", header=True):
    '''Create telemetry matrix from read files; any per-read summary
       statistics that would be useful to know'''
    try:
        h5File = h5py.File(fileName, 'r')
        h5File.close()
    except:
        return False
    with h5py.File(fileName, 'r') as h5File:
        rowData = get_telemetry(h5File, callID)
        if(header):
            sys.stdout.write(",".join(rowData.keys()) + "\n")
            # here's the raw to pA formula for future reference:
            # pA = (raw + offset)*range/digitisation
            # (using channelMeta[("offset", "range", "digitisation")])
        sys.stdout.write(",".join(map(str,rowData.values())) + "\n")

def generate_raw(fileName, callID="000", medianWindow=21):
    '''write out raw sequence from fast5, with optional running median
       smoothing, return False if not present'''
    try:
        h5File = h5py.File(fileName, 'r')
        h5File.close()
    except:
        return False
    with h5py.File(fileName, 'r') as h5File:
      runMeta = h5File['UniqueGlobalKey/tracking_id'].attrs
      channelMeta = h5File['UniqueGlobalKey/channel_id'].attrs
      runID = '%s_%s' % (runMeta["device_id"],runMeta["run_id"][0:16])
      eventBase = "/Raw/Reads"
      if(not eventBase in h5File):
          return False
      readNames = h5File[eventBase]
      readNameStr = ""
      for readName in readNames:
        readRawLocation = "%s/%s/Signal" % (eventBase, readName)
        outData = h5File[readRawLocation][()] # load entire raw data into memory
        if(medianWindow==1):
            sys.stdout.write(outData)
        else:
            array("H",runningMedian(outData, M=medianWindow)).tofile(sys.stdout)

def generate_dir_raw(fileName, callID="000", medianWindow=1, direction=None):
    '''write out directional raw sequence from fast5, return False if not present'''
    try:
        h5File = h5py.File(fileName, 'r')
        h5File.close()
    except:
        return False
    with h5py.File(fileName, 'r') as h5File:
      runMeta = h5File['UniqueGlobalKey/tracking_id'].attrs
      channelMeta = h5File['UniqueGlobalKey/channel_id'].attrs
      runID = '%s_%s' % (runMeta["device_id"],runMeta["run_id"][0:16])
      eventBase = "/Raw/Reads"
      if(not eventBase in h5File):
          return False
      seqBase1D = "/Analyses/Basecall_1D_%s" % callID
      dir = "complement" if (direction=="r") else "template"
      eventMetaBase = "%s/BaseCalled_%s/Events" % (seqBase1D, dir)
      if(not eventMetaBase in h5File):
          return False
      eventMeta = h5File[eventMetaBase].attrs
      absRawStart = eventMeta["start_time"] * channelMeta["sampling_rate"]
      absRawEnd = (eventMeta["start_time"]
                + eventMeta["duration"]) * channelMeta["sampling_rate"]
      readNames = h5File[eventBase]
      readNameStr = ""
      for readName in readNames:
        readRawMeta = h5File["%s/%s" % (eventBase, readName)].attrs
        relRawStart = int(absRawStart - readRawMeta["start_time"])
        relRawEnd = int(absRawEnd - readRawMeta["start_time"])
        sys.stderr.write("Writing (%d..%d) from %s\n" %
                         (relRawStart, relRawEnd, readName))
        readRawLocation = "%s/%s/Signal" % (eventBase, readName)
        signal = h5File[readRawLocation][relRawStart:relRawEnd] # subset for direction
        ## Remove extreme values from signal
        meanSig = sum(signal) / len(signal)
        madSig = sum(map(lambda x: abs(x - meanSig), signal)) / len(signal)
        minSig = meanSig - madSig * 6
        maxSig = meanSig + madSig * 6
        rangeFilt = numpy.vectorize(lambda x: meanSig if
                                    ((x < minSig) or (x > maxSig)) else x);
        signal = rangeFilt(signal)
        sys.stdout.write(signal) # write to file

def usageQuit():
    sys.stderr.write('Error: No file or directory provided in arguments\n\n')
    sys.stderr.write('Usage: %s <dataType> <fast5 file name>\n' % sys.argv[0])
    sys.stderr.write(' where <dataType> is one of the following:\n')
    sys.stderr.write('  fastq     - extract base-called fastq data\n')
    sys.stderr.write('  event     - extract uncalled model event matrix\n')
    sys.stderr.write('  consensus - extract consensus alignment matrix\n')
    sys.stderr.write('  eventfwd  - extract model event matrix (template)\n')
    sys.stderr.write('  eventrev  - extract model event matrix (complement)\n')
    sys.stderr.write('  telemetry - extract read statistics matrix\n')
    sys.stderr.write('  raw       - extract raw data without smoothing\n')
    sys.stderr.write('  rawfwd    - extract raw data from template\n')
    sys.stderr.write('  rawrev    - extract raw data from complement\n')
    sys.stderr.write('  rawsmooth - raw data, running-median smoothing\n')
    sys.exit(1)

if len(sys.argv) < 3:
    usageQuit()

dataType = sys.argv[1]
if(not dataType in ("fastq", "fasta", "event", "consensus", "eventfwd",
                    "eventrev", "telemetry", "raw", "rawfwd", "rawrev", "rawsmooth")):
    sys.stderr.write('Error: Incorrect dataType\n\n')
    usageQuit()

fileArg = sys.argv[2]
seenHeader = False

if(os.path.isdir(fileArg)):
    sys.stderr.write("Processing directory '%s':\n" % fileArg)
    for dirPath, dirNames, fileNames in os.walk(fileArg):
        fc = len(fileNames)
        for fileName in fileNames:
            if(fileName.endswith(".fast5")): # only process fast5 files
                if((fc == 2) or ((fc-1) % 100 == 0)):
                    sys.stderr.write("  Processing file '%s'..." % fileName)
                if(dataType == "event"):
                    generate_event_matrix(os.path.join(dirPath, fileName), header=not seenHeader)
                elif(dataType == "consensus"):
                    generate_consensus_matrix(os.path.join(dirPath, fileName), header=not seenHeader)
                elif(dataType == "telemetry"):
                    generate_telemetry(os.path.join(dirPath, fileName), header=not seenHeader)
                elif(dataType == "fastq"):
                    generate_fastq(os.path.join(dirPath, fileName))
                elif(dataType == "raw"):
                    sys.stderr.write(" Error: raw output only works for single files!\n")
                    usageQuit()
                fc -= 1
                seenHeader = True
                if(fc == 1):
                    sys.stderr.write(" done (%d more file to process)\n" % fc)
                elif(fc % 100 == 0):
                    sys.stderr.write(" done (%d more files to process)\n" % fc)
elif(os.path.isfile(fileArg)):
    if(dataType == "event"):
        generate_event_matrix(fileArg)
    elif(dataType == "consensus"):
        generate_consensus_matrix(fileArg)
    elif(dataType == "eventfwd"):
        generate_eventdir_matrix(fileArg, direction="f")
    elif(dataType == "eventrev"):
        generate_eventdir_matrix(fileArg, direction="r")
    elif(dataType == "telemetry"):
        generate_telemetry(fileArg)
    elif(dataType == "fastq"):
        generate_fastq(fileArg)
    elif(dataType == "rawsmooth"):
        generate_raw(fileArg, medianWindow=21)
    elif(dataType == "raw"):
        generate_raw(fileArg)
    elif(dataType == "rawfwd"):
        generate_dir_raw(fileArg, direction="f")
    elif(dataType == "rawrev"):
        generate_dir_raw(fileArg, direction="r")
else:
    usageQuit()
