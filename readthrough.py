#!/usr/bin/python
import HTSeq
import string
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import re
import sys

def writeSeqTrans(name, exons, seq, interval):
    if(interval.strand == "-"):
        #print(">%s (rc)\n%s" % (name, seq))
        seq = string.replace(seq,"|","")
        translation = str(Seq(seq, generic_dna).reverse_complement().translate())
        endText = " from RC strand"
    else:
        #print(">%s\n%s" % (name, seq))
        seq = string.replace(seq,"|","")
        translation = str(Seq(seq, generic_dna).translate())
        endText = ""
    translation = string.replace(translation, "*", "X", 1)
    translation = re.sub("\\*.*$","*",translation)
    if(translation.endswith("X")):
        endText += ", no sequence beyond stop codon"
        translation = string.replace(translation, "X", "*", 1)
        throughDist = 0
    elif(not translation.endswith("*")):
        endText += ", no additional stop codons found"
        throughDist = 0
    else:
        throughDist = translation.find("*") - translation.find("X")
        endText += ", stop codon distance: %d" % throughDist
    name = "%05d_%s" % (throughDist, name)
    if(numExons == 1):
        print(">%s [translation of %d exon%s]\n%s" % (name, exons, endText, translation))
    else:
        print(">%s [translation of %d exons%s]\n%s" % (name, exons, endText, translation))

def writeSeqMrna(name, exons, seq, interval):
    if(interval.strand == "-"):
        #print(">%s (rc)\n%s" % (name, seq))
        seq = string.replace(seq,"|","")
        translation = str(Seq(seq, generic_dna).reverse_complement().translate())
        endText = " from RC strand"
    else:
        #print(">%s\n%s" % (name, seq))
        seq = string.replace(seq,"|","")
        translation = str(Seq(seq, generic_dna).translate())
        endText = ""
    translation = string.replace(translation, "*", "X", 1)
    translation = re.sub("\\*.*$","*",translation)
    if(translation.endswith("X")):
        endText += ", no sequence beyond stop codon"
        translation = string.replace(translation, "X", "*", 1)
        throughDist = 0
    elif(not translation.endswith("*")):
        endText += ", no additional stop codons found"
        throughDist = 0
    else:
        throughDist = translation.find("*") - translation.find("X")
        endText += ", stop codon distance: %d" % throughDist
    name = "%05d_%s" % (throughDist, name)
    if(numExons == 1):
        print(">%s [translation of %d exon%s]\n%s" % (name, exons, endText, translation))
    else:
        print(">%s [translation of %d exons%s]\n%s" % (name, exons, endText, translation))


# load Scer reference genome
scerFile = open('saccharomyces_cerevisiae.gff', 'r')
gffLines = ()
fastaLines = ()
hitFasta = False
for(line in scerFile):
    if(line.startswith('##FASTA')):
        hitFasta = True
    if(not hitFasta):
        gfflines.append(line)
    else if(not line.startswith('#')):
        fastaLines.append(line)


gtfFile = HTSeq.GFF_Reader(gffLines)
# load all sequences into memory
sequences = dict()
for s in HTSeq.FastaReader(fastaLines):
    s.seq = string.replace(s.seq,"\r","")
    sequences[s.name] = s

lastSequence = ""
lastName = ""
lastInterval = None
extend = 50 # number of base pairs to extend
readDistance = 5000
numExons = 0
numFeatures = 0
sys.stderr.write("\n")
for feature in gtfFile:
    sys.stderr.write("\rFeatures read: %d" % numFeatures)
    numFeatures += 1
    if(feature.type != "CDS"):
        continue
    if(feature.name != lastName):
        if((lastInterval != None) and (lastInterval.strand != "-")):
            # add readDistance bases to end of sequence
            lastSequence += sequences[lastInterval.chrom].seq[lastInterval.end:lastInterval.end+readDistance]
        # write out sequence / translation
        if(lastName != ""):
            writeSeqTrans(lastName, numExons, lastSequence, lastInterval)
        lastSequence = ""
        numExons = 0
        if(feature.iv.strand == "-"):
            # add readDistance bases to start of sequence
            startPos = lastInterval.start - readDistance
            if(startPos < 0):
                startPos = 0
            lastSequence += sequences[feature.iv.chrom].seq[startPos:lastInterval.end]
        lastName = feature.name
        lastInterval = feature.iv
    # add current interval to sequence
    lastSequence += sequences[feature.iv.chrom].seq[feature.iv.start:feature.iv.end]
    lastSequence += "|"
    numExons += 1

# write out sequence / translation
writeSeqTrans(lastName, numExons, lastSequence, lastInterval)

sys.stderr.write("\n")
