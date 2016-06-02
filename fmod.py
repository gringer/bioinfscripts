#!/usr/bin/python

import sys
import array
import struct
import wave
from math import sin, asin, pi, log, exp, sqrt

## Carry out frequency modulation; vary the frequency of a wave based
## on a signal
##
## Here is the rough process:
# 1) find out the phase of the current frequency at the previous position
# 2) calculate the new amplitude by adding one time slice to that phase

def fmod(outFile, signal, minFreq, maxFreq, oldRate, newRate,
         speed=1.0, volume=0.8, logScale=True):
    oldRate = oldRate * speed
    logRange = log(maxFreq) - log(minFreq)
    linRange = maxFreq - minFreq
    meanSig = sum(signal) / len(signal)
    madSig = sum(map(lambda x: abs(x - meanSig), signal)) / len(signal)
    minSig = meanSig - madSig * 4;
    maxSig = meanSig + madSig * 4;
    if(min(signal) > minSig):
        minSig = min(signal)
    if(max(signal) < maxSig):
        maxSig = max(signal)
    ## limit signal to within MAD range, scale to 0..100
    signal = map(lambda x: float(0) if (x < minSig) else
                 (float(99) if x > maxSig else
                  float(99) * (x - minSig) / (maxSig - minSig)), signal)
    newFreqs = (map(
        lambda x: exp((log(x + 1) / log(100)) * logRange + log(minFreq)),
        signal)) if logScale else (
            map(lambda x: (x/100) * linRange + minFreq, signal))
    # number of sound samples per signal sample
    newPerOld  = (float(newRate) / (float(oldRate)))
    amp = [0] * int(len(signal) * newPerOld)
    # start in first quadrant [0..pi/2)
    quadrant = 0
    oldPhase = 0
    oldAmplitude = 0
    for s in xrange(len(signal)-1):
        ## work out phase of previous signal
        if(quadrant == 0):
            oldPhase = asin(oldAmplitude)
        if(quadrant == 1):
            oldPhase = pi - asin(oldAmplitude)
        if(quadrant == 2):
            oldPhase = pi - asin(oldAmplitude)
        if(quadrant == 3):
            oldPhase = 2 * pi + asin(oldAmplitude)
        newFreq = newFreqs[s+1]
        ## determine phase step (for input)
        sigPhaseStep = (newFreq * 2 * pi / oldRate)
        ## determine phase step (for output)
        outPhaseStep = (newFreq * 2 * pi / newRate)
        ## determine new phase
        newPhase = (oldPhase + sigPhaseStep) % (2 * pi)
        ## calculate new quadrant
        quadrant = int(newPhase / (pi/2))
        ## determine new amplitude
        oldAmplitude = sin(newPhase)
        ## write signal to file
        sStart = int(s*newPerOld);
        sEnd = int((s+1)*newPerOld);
        packedSamples = map(
            lambda x: struct.pack('h',int(sin(oldPhase + x*outPhaseStep) *
                                          volume * 32767)),
            xrange(sEnd-sStart))
        fmodOut.writeframes(''.join(packedSamples))

rate=int(sys.argv[2])

with open(sys.argv[1], "rb") as f:
    outRate = rate*15
    data = array.array('H', f.read())
    #print(",".join(map(str,newData[:50])))
    fmodOut = wave.open('out.wav', 'w')
    fmodOut.setparams((1, 2, outRate, 0, 'NONE', 'not compressed'))
    fmod(outFile=fmodOut, signal=data, minFreq=50, maxFreq=rate/4,
         speed=0.0625,
         oldRate=rate, newRate=outRate, volume=0.03)
    fmodOut.close()

    #with open ("out.raw", "wb") as soundFile:
    #    soundFile.write(newData)
