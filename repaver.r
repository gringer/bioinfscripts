#!/usr/bin/Rscript

## repaver.r - Repetitive pattern visualiser for extremely long reads

library(reticulate);

dnaSeqFile <- if(length(commandArgs(TRUE) > 0)){
                  commandArgs(TRUE)[1];
              } else "data/circ-Nb-ec3-mtDNA.fasta";

## Create python function to quickly generate kmer location dictionary
py_run_string("
from Bio import SeqIO
from collections import defaultdict
from itertools import groupby
def getKmerLocs(seqFile, kSize=17):
   fileKmers = dict()
   for record in SeqIO.parse(seqFile, \"fasta\"):
      kmers = defaultdict(list)
      seq = str(record.seq)
      for k, v in zip([seq[d:d+kSize] for d in xrange(len(seq)-kSize+1)], xrange(len(seq)-kSize+1)):
         kmers[k].append(v)
      fileKmers[record.id] = kmers
   return(fileKmers)
");

## Generate kmer location dictionary
system.time(res <- py$getKmerLocs(dnaSeqFile, as.integer(17)));
