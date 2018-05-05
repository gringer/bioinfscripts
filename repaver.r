#!/usr/bin/Rscript
## repaver.r - REpetitive PAttern Visualiser for Extremely-long Reads

library(reticulate);

outputStyle <- "circular";
kmerLength <- 10;

dnaSeqFile <- if(length(commandArgs(TRUE) > 0)){
                  commandArgs(TRUE)[1];
              } else "data/circ-Nb-ec3-mtDNA.fasta";

## helper R functions
fib.divs <- round(10^((0:4)/5) * 2) * 0.5; ## splits log decades into 5

valToSci <- function(val, unit = ""){
    sci.prefixes <- c("", "k", "M", "G", "T", "P", "E", "Z", "Y");
    units <- rep(paste(sci.prefixes,unit,sep=""), each=3);
    logRegion <- floor(log10(val))+1;
    conv.units <- units[logRegion];
    conv.div <- 10^rep(0:(length(sci.prefixes)-1) * 3, each = 3)[logRegion];
    conv.val <- val / conv.div;
    conv.val[val == 0] <- 0;
    conv.units[val == 0] <- unit;
    return(sprintf("%s %s",conv.val,conv.units));
}

## Create python function to quickly generate a kmer location
## dictionary, then filter out unique kmers
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

## Generate filtered kmer location dictionary
system.time(res <- py$getKmerLocs(dnaSeqFile, as.integer(kmerLength)));

for(dnaSeqMapName in names(res)){
    dnaSeqMap <- res[[dnaSeqMapName]];
    if(outputStyle == "dotplot"){
        png(width=1920, height=1080, pointsize=18);
    } else {
        png(width=1920, height=1920, pointsize=24);
    }
    sLen <- dnaSeqMap$length;
    if(outputStyle == "dotplot"){
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
    } else if(outputStyle == "logplot"){
        par(mgp=c(2.5,1,0), mar=c(4,6,3,0.5),
            cex.axis=1.5, cex.lab=1.5, cex.main=2);
        plot(NA, xlim=c(0,sLen), ylim=c(1,sLen), log="y",
             xlab=ifelse(sLen >= 10^6, "Base Location (Mb)", "Base Location (kb)"),
             ylab="",
             axes=FALSE,
             main=sprintf("%s (k=%d)", dnaSeqMapName, kmerLength));
        if(sLen >= 10^6){
            axis(1, at=axTicks(1), labels=pretty(axTicks(1))/10^6, lwd=3);
        } else {
            axis(1, at=axTicks(1), labels=pretty(axTicks(1))/1000, lwd=3);
        }
        drMax <- ceiling(log10(sLen));
        axis(2, at= 10^(0:drMax), las=2, lwd=3, cex.axis=1.5,
             labels=valToSci(10^(0:drMax)));
        axis(2, at= rep(1:9, each=drMax+1) * 10^(0:drMax), labels=FALSE);
        abline(h=10^(0:drMax), col="#80808050", lwd = 3);
        mtext("Feature distance (bp)", 2, line=4.5, cex=1.5);
    } else if(outputStyle == "circular"){
        par(mgp=c(2.5,1,0), mar=c(3,3,3,3),
            cex.axis=1.5, cex.lab=1.5, cex.main=2);
        plot(NA, xlim=c(-1,1), ylim=c(-1,1), axes=FALSE, ann=FALSE,
             main=sprintf("%s (k=%d)", dnaSeqMapName, kmerLength));
    }
    revNames <- sapply(names(dnaSeqMap), py$rev);
    revCNames <- sapply(names(dnaSeqMap), py$rc);
    compNames <- sapply(names(dnaSeqMap), py$comp);
    repeatedKmers <- sapply(dnaSeqMap, function(x){length(x) > 1});
    rcKmers <- revCNames %in% names(dnaSeqMap);
    rKmers <- revNames %in% names(dnaSeqMap);
    cKmers <- compNames %in% names(dnaSeqMap);
    ## f,c,rc,r : red, orange, blue, green
    plotPointsF <- NULL;
    plotPointsC <- NULL;
    plotPointsRC <- NULL;
    plotPointsR <- NULL;
    for(kposs in dnaSeqMap[repeatedKmers]){
        plotPointsF <- rbind(plotPointsF, 
                            data.frame(x=rep(kposs, length(kposs)),
                                       y=rep(kposs, each=length(kposs))));
    }
    for(kmer in names(dnaSeqMap)[cKmers]){
        kposs <- dnaSeqMap[[kmer]];
        oposs <- dnaSeqMap[[py$comp(kmer)]];
        plotPointsC <- rbind(plotPointsC,
                            data.frame(x=rep(kposs, length(oposs)),
                                       y=rep(oposs, each=length(kposs))));
    }
    for(kmer in names(dnaSeqMap)[rcKmers]){
        kposs <- dnaSeqMap[[kmer]];
        oposs <- dnaSeqMap[[py$rc(kmer)]];
        plotPointsRC <- rbind(plotPointsRC,
                            data.frame(x=rep(kposs, length(oposs)),
                                       y=rep(oposs, each=length(kposs))));
    }
    for(kmer in names(dnaSeqMap)[rKmers]){
        kposs <- dnaSeqMap[[kmer]];
        oposs <- dnaSeqMap[[py$rev(kmer)]];
        plotPointsR <- rbind(plotPointsR,
                            data.frame(x=rep(kposs, length(oposs)),
                                       y=rep(oposs, each=length(kposs))));
    }
    if(outputStyle == "dotplot"){
        points(plotPointsF, pch=15, col="#8b000040", cex=0.5);
        points(plotPointsC, pch=15, col="#FF7F0040", cex=0.5);
        points(plotPointsRC, pch=15, col="#0000FF40", cex=0.5);
        points(plotPointsR, pch=15, col="#00A00040", cex=0.5);
        legend("bottomleft",
               legend=c("Forward","Complement","RevComp","Reverse"),
               fill=c("#8b000040","#FF7F0040","#0000FF40","#00A00040"),
               bg="#FFFFFFE0", inset=0.05);
    } else if(outputStyle == "logplot"){
        plotPointsF$dist <-  plotPointsF$y -  plotPointsF$x;
        plotPointsC$dist <-  plotPointsC$y -  plotPointsC$x;
        plotPointsRC$dist <- plotPointsRC$y - plotPointsRC$x;
        plotPointsR$dist <-  plotPointsR$y -  plotPointsR$x;
        plotPointsF <-  subset(plotPointsF, dist > 0);
        plotPointsC <-  subset(plotPointsC, dist > 0);
        plotPointsRC <- subset(plotPointsRC, dist > 0);
        plotPointsR <-  subset(plotPointsR, dist > 0);
        points(plotPointsF$x,  plotPointsF$dist,
               pch=15, col="#8b000040", cex=0.5); # red
        points(plotPointsC$x,  plotPointsC$dist,
               pch=15, col="#FDC08640", cex=0.5); # salmon
        points(plotPointsRC$x, plotPointsRC$dist,
               pch=15, col="#0000FF40", cex=0.5); # blue
        points(plotPointsR$x,  plotPointsR$dist,
               pch=15, col="#00A00040", cex=0.5); # green
        points(plotPointsF$y,  plotPointsF$dist,
               pch=15, col="#9000A040", cex=0.5); # magenta
        points(plotPointsC$y,  plotPointsC$dist,
               pch=15, col="#FF7F0040", cex=0.5); # orange
        points(plotPointsRC$y, plotPointsRC$dist,
               pch=15, col="#00A09040", cex=0.5); # cyan
        points(plotPointsR$y,  plotPointsR$dist,
               pch=15, col="#A0900040", cex=0.5); # yellow
        legend(x = "bottom",
               fill=c("#9000a0","#8b0000",
                      "#fdc086","#ff7f00",
                      "#00a090","#0000ff",
                      "#a09000","#00a000"),
               legend=c("Repeat (L)",  "Repeat (R)",
                        "Comp (L)",    "Comp (R)",
                        "RevComp (L)", "RevComp (R)",
                        "Reverse (L)", "Reverse (R)"),
               bg="#FFFFFFE0", horiz=FALSE, inset=0.01, ncol=4);
    } else if(outputStyle == "circular"){
        plotPointsF$dist <- ifelse(plotPointsF$x > plotPointsF$y,
                                   sLen + (plotPointsF$y - plotPointsF$x),
                                   pmin(plotPointsF$y - plotPointsF$x));
        plotPointsC$dist <- ifelse(plotPointsC$x > plotPointsC$y,
                                   sLen + (plotPointsC$y - plotPointsC$x),
                                   pmin(plotPointsC$y - plotPointsC$x));
        plotPointsRC$dist <- ifelse(plotPointsRC$x > plotPointsRC$y,
                                   sLen + (plotPointsRC$y - plotPointsRC$x),
                                   pmin(plotPointsRC$y - plotPointsRC$x));
        plotPointsR$dist <- ifelse(plotPointsR$x > plotPointsR$y,
                                   sLen + (plotPointsR$y - plotPointsR$x),
                                   pmin(plotPointsR$y - plotPointsR$x));
        plotPointsF <-  subset(plotPointsF, dist <= sLen/2);
        plotPointsC <-  subset(plotPointsC, dist <= sLen/2);
        plotPointsRC <- subset(plotPointsRC, dist <= sLen/2);
        plotPointsR <-  subset(plotPointsR, dist <= sLen/2);
        plotPointsF$r <- 1 * (1 - sqrt(plotPointsF$dist) / sqrt(sLen));
        plotPointsC$r <- 1 * (1 - sqrt(plotPointsC$dist) / sqrt(sLen));
        plotPointsRC$r <- 1 * (1 - sqrt(plotPointsRC$dist) / sqrt(sLen));
        plotPointsR$r <- 1 * (1 - sqrt(plotPointsR$dist) / sqrt(sLen));
        points(plotPointsF$r*cos(plotPointsF$x/sLen*2*pi),
               plotPointsF$r*sin(plotPointsF$x/sLen*2*pi),
               pch=15, col="#8b000040", cex=0.5); # red
        points(plotPointsF$r*cos(plotPointsF$y/sLen*2*pi),
               plotPointsF$r*sin(plotPointsF$y/sLen*2*pi),
               pch=15, col="#9000A040", cex=0.5); # magenta
        points(plotPointsC$r*cos(plotPointsC$x/sLen*2*pi),
               plotPointsC$r*sin(plotPointsC$x/sLen*2*pi),
               pch=15, col="#FDC08640", cex=0.5); # salmon
        points(plotPointsC$r*cos(plotPointsC$y/sLen*2*pi),
               plotPointsC$r*sin(plotPointsC$y/sLen*2*pi),
               pch=15, col="#FF7F0040", cex=0.5); # orange
        points(plotPointsRC$r*cos(plotPointsRC$x/sLen*2*pi),
               plotPointsRC$r*sin(plotPointsRC$x/sLen*2*pi),
               pch=15, col="#0000FF40", cex=0.5); # blue
        points(plotPointsRC$r*cos(plotPointsRC$y/sLen*2*pi),
               plotPointsRC$r*sin(plotPointsRC$y/sLen*2*pi),
               pch=15, col="#00A09040", cex=0.5); # cyan
        points(plotPointsR$r*cos(plotPointsR$x/sLen*2*pi),
               plotPointsR$r*sin(plotPointsR$x/sLen*2*pi),
               pch=15, col="#00A00040", cex=0.5); # green
        points(plotPointsR$r*cos(plotPointsR$y/sLen*2*pi),
               plotPointsR$r*sin(plotPointsR$y/sLen*2*pi),
               pch=15, col="#A0900040", cex=0.5); # yellow
        legend(x = "bottom",
               fill=c("#9000a0","#8b0000",
                      "#fdc086","#ff7f00",
                      "#00a090","#0000ff",
                      "#a09000","#00a000"),
               legend=c("Repeat (L)",  "Repeat (R)",
                        "Comp (L)",    "Comp (R)",
                        "RevComp (L)", "RevComp (R)",
                        "Reverse (L)", "Reverse (R)"),
               bg="#FFFFFFE0", horiz=FALSE, inset=0.01, ncol=4);
    }
    invisible(dev.off());
}
