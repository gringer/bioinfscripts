#!/usr/bin/Rscript
## repaver.r - REpetitive PAttern Visualiser for Extremely-long Reads

library(reticulate);

dnaSeqFile <- if(length(commandArgs(TRUE) > 0)){
                  commandArgs(TRUE)[1];
              } else "data/circ-Nb-ec3-mtDNA.fasta";

## Create python function to quickly generate kmer location dictionary
py_run_string("
from string import maketrans
from Bio import SeqIO
from collections import defaultdict
compTransTable = maketrans('ACGTUYRSWMKDVHBXNacgtuyrswmkdvhbxn',
                           'TGCAARYSWKMHBDVXNtgcaaryswkmhbdvxn');
def comp(seq):
  return(seq.translate(compTransTable))
def rev(seq):
  return(seq[::-1])
def rc(seq):
  return(seq.translate(compTransTable)[::-1])
def getKmerLocs(seqFile, kSize=17):
   fileKmers = dict()
   for record in SeqIO.parse(seqFile, \"fasta\"):
      kmers = defaultdict(list)
      seq = str(record.seq)
      for k, v in zip([seq[d:d+kSize] for d in
            xrange(len(seq)-kSize+1)], xrange(len(seq)-kSize+1)):
         kmers[k].append(v)
      kmers = {k:v for k, v in kmers.iteritems() if (
        (not 'N' in k) and (
           (len(v) > 1) or             ## repeated
           (k[::-1] in kmers) or       ## reverse
           (comp(k) in kmers) or       ## complement
           (comp(k[::-1]) in kmers)))} ## reverse complement
      kmers['length'] = len(seq) ## add in length, because it's cheap
      fileKmers[record.id] = kmers
   return(fileKmers)
");

kmerLength <- 10;

## Generate kmer location dictionary
system.time(res <- py$getKmerLocs(dnaSeqFile, as.integer(kmerLength)));

for(dnaSeqMapName in names(res)){
    dnaSeqMap <- res[[dnaSeqMapName]];
    sLen <- dnaSeqMap$length;
    par(mgp=c(2,0.5,0));
    plot(NA, xlim=c(0,sLen), ylim=c(sLen,0),
         xlab=ifelse(sLen >= 10^6, "Base Location (Mb)", "Base Location (kb)"),
         ylab=ifelse(sLen >= 10^6, "Base Location (Mb)", "Base Location (kb)"),
         axes=FALSE,
         main=sprintf("%s (k=%d)", dnaSeqMapName, kmerLength));
    if(sLen >= 10^6){
        axis(1, at=axTicks(1), labels=pretty(axTicks(1))/10^6);
        axis(2, at=rev(axTicks(2)), labels=pretty(axTicks(2))/10^6);
    } else {
        axis(1, at=axTicks(1), labels=pretty(axTicks(1))/1000);
        axis(2, at=rev(axTicks(2)), labels=pretty(axTicks(2))/1000);
    }
    revNames <- sapply(names(dnaSeqMap), py$rev);
    revCNames <- sapply(names(dnaSeqMap), py$rc);
    compNames <- sapply(names(dnaSeqMap), py$comp);
    repeatedKmers <- sapply(dnaSeqMap, function(x){length(x) > 1});
    rcKmers <- revCNames %in% names(dnaSeqMap);
    rKmers <- revNames %in% names(dnaSeqMap);
    cKmers <- compNames %in% names(dnaSeqMap);
    ## f,c,rc,r : red, orange, blue, green
    plotPoints <- NULL;
    for(kposs in dnaSeqMap[repeatedKmers]){
        plotPoints <- rbind(plotPoints, 
                            data.frame(x=rep(kposs, length(kposs)),
                                       y=rep(kposs, each=length(kposs))));
    }
    points(plotPoints, pch=15, col="#8b000040", cex=0.5);
    plotPoints <- NULL;
    for(kmer in names(dnaSeqMap)[cKmers]){
        kposs <- dnaSeqMap[[kmer]];
        oposs <- dnaSeqMap[[py$comp(kmer)]];
        plotPoints <- rbind(plotPoints, 
                            data.frame(x=rep(kposs, length(oposs)),
                                       y=rep(oposs, each=length(kposs))));
    }
    points(plotPoints, pch=15, col="#FF7F0040", cex=0.5);
    plotPoints <- NULL;
    for(kmer in names(dnaSeqMap)[rcKmers]){
        kposs <- dnaSeqMap[[kmer]];
        oposs <- dnaSeqMap[[py$rc(kmer)]];
        plotPoints <- rbind(plotPoints, 
                            data.frame(x=rep(kposs, length(oposs)),
                                       y=rep(oposs, each=length(kposs))));
    }
    points(plotPoints, pch=15, col="#0000FF40", cex=0.5);
    plotPoints <- NULL;
    for(kmer in names(dnaSeqMap)[rKmers]){
        kposs <- dnaSeqMap[[kmer]];
        oposs <- dnaSeqMap[[py$rev(kmer)]];
        plotPoints <- rbind(plotPoints, 
                            data.frame(x=rep(kposs, length(oposs)),
                                       y=rep(oposs, each=length(kposs))));
    }
    points(plotPoints, pch=15, col="#00A00040", cex=0.5);
    legend("bottomleft",
           legend=c("Forward","Complement","RevComp","Reverse"),
           fill=c("#8b000040","#FF7F0040","#0000FF40","#00A00040"),
           bg="#FFFFFFE0", inset=0.05);
}
