#!/usr/bin/python

import sys
import array
from math import sin, asin, pi, log, exp

## Carry out frequency modulation; vary the frequency of a wave based
## on a signal
##
## Here is the rough process:
# 1) find out the phase of the current frequency at the previous position
# 2) calculate the new amplitude by adding one time slice to that phase

def fmod(signal, minFreq, maxFreq, rate):
    logRange = log(maxFreq) - log(minFreq)
    newFreqs = map(lambda x: exp((log(x+1) / log(65536)) *
                                logRange + log(minFreq)), signal)
    amp = [0] * len(signal)
    # start in first quadrant [0..pi/2)
    quadrant = 0
    for s in xrange(len(signal)-1):
        ## work out phase of previous signal
        oldPhase = asin(amp[s]) + quadrant * (pi/2)
        newFreq = newFreqs[s+1]
        ## determine phase step
        phaseStep = (newFreq * 2 * pi / rate)
        ## determine new phase
        newPhase = (oldPhase + phaseStep) % (2 * pi)
        ## calculate new quadrant
        quadrant = int(newPhase / (pi/2))
        ## determine new amplitude
        amp[s+1] = sin(newPhase)
    return map(lambda x: int(x * 32767), amp);
#    return(newAmp)

rate=int(sys.argv[2])

with open(sys.argv[1], "rb") as f:
    data = array.array('H', f.read())
    newData = array.array('h', fmod(data, 220, 1500, rate))
    print(",".join(map(str,newData[:50])))
    with open ("out.raw", "wb") as soundFile:
        soundFile.write(newData)
