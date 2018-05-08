#!/usr/bin/Rscript
## repaver.r - REpetitive PAttern Visualiser for Extremely-long Reads

outputStyle <- "semicircular"; ## dotplot, circular, semicircular, profile
outputType <- "png";
kmerLength <- 17;
filePrefix <- "repaver";
fileName <- "";

dnaSeqFile <- "data/circ-Nb-ec3-mtDNA.fasta";

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

usage <- function(){
  cat("usage: ./repaver.r",
      "<fasta/fastq file> [options]\n");
  cat("\nOther Options:\n");
  cat("-help            : Only display this help message\n");
  cat("-k <int>         : Set kmer length\n");
  cat("-style <string>  : Output file style",
      "(dotplot|profile|circular|semicircular)\n");
  cat("-type <string>   : Output file type (png|svg)\n");
  cat("-prefix <string> : Output file prefix (default: 'repaver')\n");
  cat("\n");
}

argLoc <- 1;

if(length(commandArgs(TRUE)) < 1){
      usage();
      quit(save = "no", status=0);
}

while(!is.na(commandArgs(TRUE)[argLoc])){
  if(file.exists(commandArgs(TRUE)[argLoc])){ # file existence check
      dnaSeqFile <- commandArgs(TRUE)[argLoc];
  } else {
    if(commandArgs(TRUE)[argLoc] == "-help"){
      usage();
      quit(save = "no", status=0);
    }
    else if(commandArgs(TRUE)[argLoc] == "-k"){
      kmerLength <- as.numeric(commandArgs(TRUE)[argLoc+1]);
      argLoc <- argLoc + 1;
    }
    else if(commandArgs(TRUE)[argLoc] == "-type"){
      outputType <- commandArgs(TRUE)[argLoc+1];
      argLoc <- argLoc + 1;
    }
    else if(commandArgs(TRUE)[argLoc] == "-style"){
      outputStyle <- commandArgs(TRUE)[argLoc+1];
      argLoc <- argLoc + 1;
    }
    else if(commandArgs(TRUE)[argLoc] == "-prefix"){
      filePrefix <- commandArgs(TRUE)[argLoc+1];
      argLoc <- argLoc + 1;
    }
    else {
      cat("Error: Argument '",commandArgs()[argLoc],
          "' is not understood by this program\n\n", sep="");
      usage();
      quit(save = "no", status=0);
    }
  }
  argLoc <- argLoc + 1;
}

if(!(outputStyle %in% c("dotplot","profile","circular","semicircular"))){
    cat("Error: Plot style should be one of ",
        "'(dotplot|profile|circular|semicircular)'\n");
    usage();
    quit(save = "no", status=0);
}

if(!(outputType %in% c("svg","png"))){
    cat("Error: File output type should be one of ",
        "'(svg|png)'\n");
    usage();
    quit(save = "no", status=0);
}

library(reticulate);

## Create python function to quickly generate a kmer location
## dictionary and determine the distance from kmers that have been seen before
## This is chunked into no more than ~5000 pieces so that the final image
## generation is slightly quicker.
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
   fileChunks = dict()
   for record in SeqIO.parse(seqFile, \"fasta\"):
      kmers = defaultdict(set)
      chunks = dict()
      chunks['F'] = defaultdict(set)
      chunks['R'] = defaultdict(set)
      chunks['RC'] = defaultdict(set)
      chunks['C'] = defaultdict(set)
      seq = str(record.seq)
      seqLen = len(seq) ## add in length, because it's cheap
      baseBlockSize = int(seqLen / 5000) ## limit to 5000 position slots
      if(baseBlockSize < 1):
         baseBlockSize = 1
      for k, v in zip([seq[d:d+kSize] for d in
            xrange(len(seq)-kSize+1)], xrange(len(seq)-kSize+1)):
         krev = k[::-1]
         kcomp = comp(k)
         krc = kcomp[::-1]
         checkstr = {'F': k, 'R': krev, 'C': kcomp, 'RC':krc}
         chunkID = 'b' + str(int(v / baseBlockSize) * baseBlockSize)
         for type,ko in checkstr.iteritems():
            if(ko in kmers):
               for p in kmers[ko]:
                  chunks[type][chunkID].add(v-p)
         kmers[k].add(v)
      for k, v in chunks.iteritems():
          chunks[k] = {kv:list(vv) for kv,vv in chunks[k].iteritems()}
      fileChunks[record.id] = dict({'len':seqLen, 'blockSize':baseBlockSize,
                                   'chunks':chunks})
   return(fileChunks)
");

## Generate filtered kmer location dictionary
cat("Generating chunk difference dictionary... ");
my.time <- Sys.time();
res <- py$getKmerLocs(dnaSeqFile, as.integer(kmerLength));
cat(sprintf("done in %0.2f %s\n",
            Sys.time() - my.time, attr(Sys.time() - my.time, "units")));


fileName <- sprintf("%s_k%d.%s",filePrefix, kmerLength, outputType);
if(length(names(res)) > 1){
    fileName <- sprintf("%s_k%d_%%03d.%s",filePrefix, kmerLength, outputType);
}
fileCounter <- 1;

#print(str(res[[1]]));

if(outputStyle %in% c("profile", "semicircular")){
    if(outputType == "svg"){
        svg(filename=fileName, width=20, height=11.25, pointsize=22);
    } else {
        png(filename=fileName, width=1920, height=1080,
            pointsize=22, antialias="gray");
    }
} else {
    if(outputType == "svg"){
        svg(filename=fileName, width=11.25, height=11.25, pointsize=18);
    } else {
        png(filename=fileName, width=1080, height=1080,
            pointsize=18, antialias="gray");
    }
}
for(dnaSeqMapName in names(res)){
    dnaSeqMap <- res[[dnaSeqMapName]];
    sLen <- dnaSeqMap$len;
    sBS <- dnaSeqMap$b;
    cat(sprintf("Processing %s [length: %d; %d bases per block]\n",
                dnaSeqMapName, sLen, sBS));
    if(outputStyle == "dotplot"){
        par(mgp=c(2,0.5,0));
        plot(NA, xlim=c(0,sLen), ylim=c(sLen,0),
             xlab=ifelse(sLen >= 10^6, "Base Location (Mb)",
                         "Base Location (kb)"),
             ylab=ifelse(sLen >= 10^6, "Base Location (Mb)",
                         "Base Location (kb)"),
             axes=FALSE,
             main=sprintf("%s (k=%d)", dnaSeqMapName, kmerLength));
        if(sLen >= 10^6){
            axis(1, at=axTicks(1), labels=pretty(axTicks(1))/10^6);
            axis(2, at=rev(axTicks(2)), labels=pretty(axTicks(2))/10^6);
        } else {
            axis(1, at=axTicks(1), labels=pretty(axTicks(1))/1000);
            axis(2, at=rev(axTicks(2)), labels=pretty(axTicks(2))/1000);
        }
    } else if(outputStyle == "profile"){
        par(mgp=c(2.5,1,0), mar=c(4,6,3,0.5),
            cex.axis=1.5, cex.lab=1.5, cex.main=2);
        plot(NA, xlim=c(0,sLen), ylim=c(1,sLen), log="y",
             xlab=ifelse(sLen >= 10^6,
                         "Base Location (Mb)", "Base Location (kb)"),
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
        par(mgp=c(2.5,1,0), mar=c(2.5,2,1.5,2),
            cex.axis=1.5, cex.lab=1.5, cex.main=2);
        plot(NA, xlim=c(-1.1,1.1), ylim=c(-1.2,1),
             axes=FALSE, xlab="", ylab="",
             main=sprintf("%s (k=%d)", dnaSeqMapName, kmerLength));
    } else if(outputStyle == "semicircular"){
        par(mgp=c(2.5,1,0), mar=c(2.5,2,2.5,2),
            cex.axis=1.5, cex.lab=1.5, cex.main=2);
        xMul <- (1.2 / par()$pin[2]) * par()$pin[1];
        plot(NA, xlim=c(-xMul/2, xMul/2), ylim=c(-0.2, 1),
             axes=FALSE, xlab="", ylab="",
             main=sprintf("%s (k=%d)", dnaSeqMapName, kmerLength));
    }
    ## f,c,rc,r : red, orange, blue, green
    plotPoints <- NULL;
    dc <- dnaSeqMap$chunks;
    for(type in as.character(names(dnaSeqMap$chunks)[sapply(dc,length)>0])){
        if(length(dc[[type]]) == 0){
            next;
        }
        my.time <- Sys.time();
        cat(sprintf("Processing %s... ",
            c(F="repeats", C="complements",
              RC="reverse complements", R="reverses")[type]));
        collectFunc <- function(kposs){
            vals=dc[[type]][[kposs]];
            data.frame(y=rep(as.numeric(substring(kposs,2)),
                             length(vals)),
                       dist=vals,
                       type=rep(type, length(vals)),
                       stringsAsFactors=FALSE);
        }
        plotPoints <-
            rbind(plotPoints,
                  Reduce(rbind,
                         sapply(names(dc[[type]]), simplify=FALSE,
                                collectFunc)));
        cat(sprintf(" done in %0.2f %s\n", Sys.time() - my.time,
                    attr(Sys.time() - my.time, "units")));
    }
    plotPoints$x <- plotPoints$y - plotPoints$dist;
    my.time <- Sys.time();
    cat("Drawing plot... ");
    if(outputStyle == "dotplot"){
        points(plotPoints, pch=15,
               col=c(F="#8b000040",RC="#FF7F0040",
                     C="#0000FF40",R="#00A00040")[plotPoints$type], cex=0.5);
        legend("bottomleft",
               legend=c("Forward","Complement","RevComp","Reverse"),
               fill=c("#8b000040","#0000FF40","#FF7F0040","#00A00040"),
               bg="#FFFFFFE0", inset=0.05);
    } else if(outputStyle == "profile"){
        ## left symbols
        points(x=plotPoints$x, y=plotPoints$dist, pch=15,
               col=c(F="#9000A040",RC="#FF7F0040",
                     C="#00A09040",R="#A0900040")[plotPoints$type], cex=0.5);
        ## right symbols
        points(x=plotPoints$y, y=plotPoints$dist, pch=15,
               col=c(F="#8b000040",RC="#FDC08640",
                     C="#0000FF40",R="#00A00040")[plotPoints$type], cex=0.5);
        legend(x = "bottom",
               fill=c("#9000a0","#8b0000",
                      "#00a090","#0000ff",
                      "#fdc086","#ff7f00",
                      "#a09000","#00a000"),
               legend=c("Repeat (L)",  "Repeat (R)",
                        "Comp (L)",    "Comp (R)",
                        "RevComp (L)", "RevComp (R)",
                        "Reverse (L)", "Reverse (R)"),
               bg="#FFFFFFE0", horiz=FALSE, inset=0.01, ncol=4);
    } else if(outputStyle == "circular"){
        distFlips <- (plotPoints$dist > sLen/2);
        plotPoints[distFlips,c("x","y","dist")] <-
            plotPoints[distFlips,c("y","x","dist")];
        plotPoints$dist[distFlips] <- (sLen - plotPoints$dist[distFlips]);
        ## Convert distance to radius. This is a piecewise function
        ## with the following properties:
        ## * Starts off as a log function
        ## * Remainder is a linear function
        ## * The transition point is the point where the slope is equal
        ## * The transition point is 1/3 along the radius
        ## * The plot ends at (sLen/2, 1)
        ## * The base of the log is sLen/12
        ## Note: slope of log[b](x) = 1/(x*log(b))
        pwFun <- function(d){
            a <- sLen/50; ## changes to linear approximately 1/3 of the way in
            aProp <- (sLen / 2) / a;
            ifelse(d < a,
                   log(d) / log(a),
                   d/(a * log(a)) + (1 - 1/log(a))) /
                ((aProp-1) / log(a) + 1) * 0.75 + 0.25;
        }
        cat("converting distances... ");
        plotPoints$r <- pwFun(plotPoints$dist);
        cat("drawing points... ");
        drMax <- ceiling(log10(sLen));
        scalePtsMajor <- rep(1, each=drMax+1) * 10^(0:drMax);
        scalePtsMajor <- c(scalePtsMajor[scalePtsMajor < sLen/2], sLen/2);
        scalePts <- rep(1:9, each=drMax+1) * 10^(0:drMax);
        scalePts <- scalePts[scalePts <= sLen/2];
        for(p in scalePtsMajor){ # rings for log scale
            points(x=pwFun(p)*cos(seq(0,2*pi, length.out=360)),
                   y=pwFun(p)*sin(seq(0,2*pi, length.out=360)),
                   type="l", lwd=3, col="#808080A0");
        }
        distPts <- (1:99)*10^(drMax-2);
        distPts <- c(head(distPts[distPts < sLen], -1), sLen);
        if(length(distPts) > 20){
            distPts <- (1:9)*10^(drMax-1);
            distPts <- signif(c(distPts[distPts < sLen], sLen),3);
        }
        segments(x0=0.18*cos(distPts / sLen * 2*pi), # tick marks
                 x1=0.2*cos(distPts / sLen * 2*pi),
                 y0=0.18*sin(distPts / sLen * 2*pi),
                 y1=0.2*sin(distPts / sLen * 2*pi),
                 lwd=2, col="#000000A0");
        for(dpi in seq_along(distPts)){ # tick labels for base location
            text(x=0.14*cos(distPts[dpi] / sLen * 2*pi),
                 y=0.14*sin(distPts[dpi] / sLen * 2*pi),
                 labels=valToSci(signif(distPts[dpi],3)), cex=0.5,
                 srt=if((distPts[dpi]/sLen * 360 >= 90) &&
                        (distPts[dpi]/sLen * 360 < 270)){
                         (distPts[dpi]/sLen * 360 + 180);
                     } else {
                         (distPts[dpi]/sLen * 360);
                     },
                 col="black");
        }
        points(x=0.2*cos(seq(10^(drMax-2)/sLen * 2*pi/2,
                             2*pi, length.out=360)),
               y=0.2*sin(seq(10^(drMax-2)/sLen * 2*pi/2,
                             2*pi, length.out=360)),
               type="l", lwd=3, col="#000000A0"); # tick circle
        ## left symbols
        points(plotPoints$r*cos(plotPoints$x/sLen*2*pi),
               plotPoints$r*sin(plotPoints$x/sLen*2*pi),
               pch=ifelse(outputType=="png",20,"•"),
               col=c(F="#9000A040",RC="#FF7F0040",
                     C="#00A09040",R="#A0900040")[plotPoints$type],
               cex=ifelse(outputType=="png",0.5,1));
        ## right symbols
        points(plotPoints$r*cos(plotPoints$y/sLen*2*pi),
               plotPoints$r*sin(plotPoints$y/sLen*2*pi),
               pch=ifelse(outputType=="png",20,"•"),
               col=c(F="#8b000040",RC="#FDC08640",
                     C="#0000FF40",R="#00A00040")[plotPoints$type],
               cex=ifelse(outputType=="png",0.5,1));
        rect(xleft=pwFun(head(scalePtsMajor,1))-0.025,
             xright=pwFun(tail(scalePtsMajor,1))+0.05,
             ytop=0.13, ybottom=-0.13, col="#FFFFFFA0", border=NA);
        arrows(x0=pwFun(head(scalePts,-1)), x1=pwFun(tail(scalePts,-1)),
               y0=0, angle=90, code=3, length=0.1, lwd=2, col="#80808080");
        arrows(x0=pwFun(head(scalePtsMajor,-1)),
               x1=pwFun(tail(scalePtsMajor,-1)),
               y0=0, angle=90, code=3, length=0.15, lwd=3, col="#00000080");
        text(x=pwFun(scalePtsMajor), y=0, col="black",
             labels=valToSci(signif(scalePtsMajor,2)), pos=1, offset=1,
             cex=0.5);
        text(x=mean(range(pwFun(scalePtsMajor))), y=0, col="black",
             labels="Feature Distance (bases)", pos=3, offset=1, cex=0.75);
        text(x=0, y=0, labels="Sequence\nLocation\n(bases)", col="black",
             cex=0.75);
        legend(x = "bottom",
               fill=c("#9000a0","#8b0000",
                      "#00a090","#0000ff",
                      "#fdc086","#ff7f00",
                      "#a09000","#00a000"),
               legend=c("Repeat (L)",  "Repeat (R)",
                        "Comp (L)",    "Comp (R)",
                        "RevComp (L)", "RevComp (R)",
                        "Reverse (L)", "Reverse (R)"),
               bg="#FFFFFFE0", horiz=FALSE, inset=0.01, ncol=4);
    } else if(outputStyle == "semicircular"){
        ## Convert distance to radius. This is a piecewise function
        ## with the following properties:
        ## * Starts off as a log function
        ## * Remainder is a linear function
        ## * The transition point is the point where the slope is equal
        ## * The transition point is 1/3 along the radius
        ## * The plot ends at (sLen/2, 1)
        ## * The base of the log is sLen/12
        ## Note: slope of log[b](x) = 1/(x*log(b))
        pwFun <- function(d){
            a <- sLen/25; ## changes to linear approximately 1/3 of the way in
            aProp <- sLen / a;
            ifelse(d < a,
                   log(d) / log(a),
                   d/(a * log(a)) + (1 - 1/log(a))) /
                ((aProp-1) / log(a) + 1) * 0.75 + 0.25;
        }
        cat("converting distances... ");
        plotPoints$r <- pwFun(plotPoints$dist);
        cat("drawing points... ");
        drMax <- ceiling(log10(sLen));
        scalePtsMajor <- rep(1, each=drMax+1) * 10^(0:drMax);
        scalePtsMajor <- c(scalePtsMajor[scalePtsMajor < sLen], sLen);
        scalePts <- rep(1:9, each=drMax+1) * 10^(0:drMax);
        scalePts <- scalePts[scalePts <= sLen];
        for(p in scalePtsMajor){ # rings for log scale
            points(x=pwFun(p)*cos(seq(0,pi, length.out=180)),
                   y=pwFun(p)*sin(seq(0,pi, length.out=180)),
                   type="l", lwd=3, col="#808080A0");
        }
        distPts <- (0:99)*10^(drMax-2);
        distPts <- c(head(distPts[distPts < sLen], -1), sLen);
        if(length(distPts) > 25){
            distPts <- (1:9)*10^(drMax-1);
            distPts <- signif(c(distPts[distPts < sLen], sLen),3);
        }
        segments(x0=-0.18*cos(distPts / sLen * pi), # tick marks for circle
                 x1=-0.2*cos(distPts / sLen * pi),
                 y0=0.18*sin(distPts / sLen * pi),
                 y1=0.2*sin(distPts / sLen * pi),
                 lwd=2, col="#000000A0");
        for(dpi in seq_along(distPts)){ # tick labels for base location
            text(x=-0.14*cos(distPts[dpi] / sLen * pi),
                 y=0.14*sin(distPts[dpi] / sLen * pi),
                 labels=valToSci(signif(distPts[dpi],3)), cex=0.5,
                 srt=-if((distPts[dpi]/sLen * 180 >= 90) &&
                        (distPts[dpi]/sLen * 180 < 270)){
                         (distPts[dpi]/sLen * 180 + 180);
                     } else {
                         (distPts[dpi]/sLen * 180);
                     },
                 col="black");
        }
        points(x=0.2*cos(seq(0, pi, length.out=180)),
               y=0.2*sin(seq(0, pi, length.out=180)),
               type="l", lwd=3, col="#000000A0"); # tick circle
        ## left symbols
        points(-plotPoints$r*cos(plotPoints$x/sLen*pi),
               plotPoints$r*sin(plotPoints$x/sLen*pi),
               pch=ifelse(outputType=="png",20,"•"),
               col=c(F="#9000A040",RC="#FF7F0040",
                     C="#00A09040",R="#A0900040")[plotPoints$type],
               cex=ifelse(outputType=="png",0.5,1));
        ## right symbols
        points(-plotPoints$r*cos(plotPoints$y/sLen*pi),
               plotPoints$r*sin(plotPoints$y/sLen*pi),
               pch=ifelse(outputType=="png",20,"•"),
               col=c(F="#8b000040",RC="#FDC08640",
                     C="#0000FF40",R="#00A00040")[plotPoints$type],
               cex=ifelse(outputType=="png",0.5,1));
        ## tick marks (left)
        arrows(x0=-pwFun(head(scalePts,-1)),
               x1=-pwFun(tail(scalePts,-1)),
               y0=0, angle=90, code=3, length=0.1, lwd=2, col="#80808080");
        arrows(x0=-pwFun(head(scalePtsMajor,-1)),
               x1=-pwFun(tail(scalePtsMajor,-1)),
               y0=0, angle=90, code=3, length=0.15, lwd=3, col="#00000080");
        ## tick marks (right)
        arrows(x0=pwFun(head(scalePts,-1)),
               x1=pwFun(tail(scalePts,-1)),
               y0=0, angle=90, code=3, length=0.1, lwd=2, col="#80808080");
        arrows(x0=pwFun(head(scalePtsMajor,-1)),
               x1=pwFun(tail(scalePtsMajor,-1)),
               y0=0, angle=90, code=3, length=0.15, lwd=3, col="#00000080");
        text(x=pwFun(scalePtsMajor), y=0, col="black",
             labels=valToSci(signif(scalePtsMajor,2)), pos=1, offset=1,
             cex=0.5); # tick labels for distance axis
        text(x=-pwFun(scalePtsMajor), y=0, col="black",
             labels=valToSci(signif(scalePtsMajor,2)), pos=1, offset=1,
             cex=0.5); # tick labels for distance axis
        text(x=-mean(range(pwFun(scalePtsMajor))), y=-0.05, col="black",
             labels="Feature Distance (bases)", pos=1, offset=1, cex=0.75);
        text(x=mean(range(pwFun(scalePtsMajor))), y=-0.05, col="black",
             labels="Feature Distance (bases)", pos=1, offset=1, cex=0.75);
        text(x=0, y=0, labels="Sequence\nLocation\n(bases)", col="black",
             cex=0.75);
        legend(x = "bottom",
               fill=c("#9000a0","#8b0000",
                      "#00a090","#0000ff",
                      "#fdc086","#ff7f00",
                      "#a09000","#00a000"),
               legend=c("Repeat (L)",  "Repeat (R)",
                        "Comp (L)",    "Comp (R)",
                        "RevComp (L)", "RevComp (R)",
                        "Reverse (L)", "Reverse (R)"),
               bg="#FFFFFFE0", horiz=FALSE, inset=0.01, ncol=4);
    }
    cat(sprintf(" done in %0.2f %s\n", Sys.time() - my.time,
                attr(Sys.time() - my.time, "units")));
    cat("Written to '",sprintf(fileName,fileCounter),"'\n",sep="");
    fileCounter <- fileCounter+1;
}
invisible(dev.off());
