#!/usr/bin/Rscript
## Signal viewer for nanopore sensor traces
## expects a file containing only raw signal data, as produced by
## 'porejuicer.py raw <file.fast5>'

## Default parameters
sigFileName <- "signal.bin";
channel <- -1;
read <- -1;
imageName <- "signal_out.pdf";
title <- "";
doPlot <- TRUE;

usage <- function(){
  cat("usage: ./signal_viewer.r",
      "<signal.bin> [options]\n");
  cat("\nOther Options:\n");
  cat("-out <file> : Write image to <file> [default: signal_out.pdf]\n");
  cat("-noplot     : Don't plot any data (just carry out calculations)\n");
  cat("-title      : Title for additional annotations\n");
  cat("-help       : Only display this help message\n");
  cat("\n");
}

argLoc <- grep("--args",commandArgs()) + 1;
if((length(argLoc) == 0) || (is.na(argLoc))){
      usage();
      quit(save = "no", status=0);
}

while(!is.na(commandArgs()[argLoc])){
  if(file.exists(commandArgs()[argLoc])){ # file existence check
    sigFileName <- commandArgs()[argLoc];
  } else {
    if(commandArgs()[argLoc] == "-help"){
      usage();
      quit(save = "no", status=0);
    }
    else if(commandArgs()[argLoc] == "-noplot"){
        doPlot <- FALSE;
    }
    else if(commandArgs()[argLoc] == "-out"){
      imageName <- commandArgs()[argLoc+1];
      argLoc <- argLoc + 1;
    }
    else if(commandArgs()[argLoc] == "-title"){
      title <- commandArgs()[argLoc+1];
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

fileLen <- file.size(sigFileName);
data.sig <- readBin(sigFileName, what=integer(), size=2, signed=FALSE,
                    n=fileLen/2);

orig.sig.len <- length(data.sig);

dMed <- median(data.sig);
dMad <- mad(data.sig);
dMin <- max(min(data.sig),dMed-4*dMad,0);
dMax <- min(max(data.sig),dMed+4*dMad,65535);

if(doPlot){
    png("untrimmed.png", width=1280, height=720, pointsize=16);
    if(orig.sig.len > 100000){
        #smoothScatter(1:length(data.sig), data.sig, main="Untrimmed raw signal (smoothed scatter plot)", xlab="Samples",
        #              ylab="Unadjusted raw signal", nbin=c(512,256), nrpoints=1000, bandwidth = c(orig.sig.len/512,20));
        plot(1:length(data.sig), data.sig, pch=".", main="Untrimmed raw signal", xlab="Samples", ylab="Unadjusted raw signal",
             ylim=c(min(data.sig),min(max(data.sig[data.sig < 2000]),2000)));
    } else {
        plot(1:length(data.sig), data.sig, type="l", main="Untrimmed raw signal", xlab="Samples", ylab="Unadjusted raw signal",
             ylim=c(min(data.sig),min(max(data.sig[data.sig < 2000]),2000)));
    }
    if(title != ""){
        mtext(title, line=0.5);
    }
}


library(caTools); ## for runmean

##rangeRLE <- rle((runmed(data.sig,11) > dMin) & (runmed(data.sig,11) < dMax));
rangeRLE <- rle(runmean(abs(runmed(data.sig,31)-runmed(data.sig,101)),1001) > 10);
startPoint <- 1;
if(length(rangeRLE$lengths) > 1){
    ## subset on longest region that fits the expected range
    indexRLE <- order(-rangeRLE$values, -rangeRLE$lengths)[1];
    startPoint <- ifelse(indexRLE == 1, 1,
                         cumsum(rangeRLE$lengths)[indexRLE-1] + 50);
    dataLen <- rangeRLE$lengths[indexRLE];
    ## adjust start point to be the first point that dips below the median value
    data.sig <- head(tail(data.sig, -startPoint), dataLen);
    ## re-adjust statistics
    if(length(data.sig) > 1){
        dMed <- median(data.sig);
        dMad <- mad(data.sig);
        dMin <- max(min(data.sig),dMed-4*dMad,0);
        dMax <- min(max(data.sig),dMed+4*dMad,65535);
    }
}

## adjust start point to be the first point that dips below the median value
startPoint <- startPoint + which.min(data.sig > dMed);
data.sig <- tail(data.sig, -which.min(data.sig > dMed));
## re-adjust statistics
if(length(data.sig) > 1){
    dMed <- median(data.sig);
    dMad <- mad(data.sig);
    dMin <- max(min(data.sig),dMed-4*dMad,0);
    dMax <- min(max(data.sig),dMed+4*dMad,65535);
}
if(doPlot){
    abline(v=c(0,length(data.sig))+startPoint, col="red", lwd=2);
    abline(h=dMed, col="red", lwd=2);
    text(x=startPoint, y=dMax, pos=4, labels=startPoint);
    text(x=startPoint + length(data.sig), y=dMax, pos=2, labels=startPoint + length(data.sig));
}

if(length(data.sig) / orig.sig.len < 0.25){
    cat("Error: signal data reduced to less than 25% of original size after noise/plateau trimming\n");
    cat(sprintf("[Remaining proportion: %0.2f%%]\n", 100 * length(data.sig) / orig.sig.len));
    quit(save="no", status=1);
    if(doPlot){
        dummy <- dev.off();
    }
}

## filter out huge signal spikes
data.sig[data.sig > dMax] <- dMax;

hpFlats <- rle(data.sig > (dMed + 1*dMad)); ## find marginal crossover points
hpFlats$values[hpFlats$lengths < 64] <- FALSE; ## exclude small gaps
hpFlats <- rle(inverse.rle(hpFlats)); ## regenerate RLE-encoding
lowPoints <- which(!hpFlats$values);
lowPoints <- lowPoints[(lowPoints > 1) & (lowPoints < length(hpFlats$values))];
hpflat.searchTable <- rbind(cumsum(hpFlats$lengths)[lowPoints-1], hpFlats$lengths[lowPoints-1],
                            hpFlats$lengths[lowPoints], hpFlats$lengths[lowPoints+1], cumsum(hpFlats$lengths)[lowPoints]);
## hairpin "low bit" length should be similar to adjacent bits
hp.ratio <- hpflat.searchTable[3,,drop=FALSE] / colSums(hpflat.searchTable[c(2,4),,drop=FALSE]);
hp.likely <- which((hp.ratio < 3) & (hp.ratio > 1/3));
hp.breakPoints <- NULL;
if(length(hp.likely) > 0){
    hp.breakPoints <- round((hpflat.searchTable[1,hp.likely] + hpflat.searchTable[5,hp.likely])/2);
}

if(doPlot){
    if(length(hp.breakPoints) >= 1){
        abline(v=startPoint + hp.breakPoints, col="#00FF0080", lwd = 3);
        text(x=startPoint + hp.breakPoints, y=dMax, pos=3, labels=startPoint + hp.breakPoints);
    }
    dummy <- dev.off();
}

## data.sig <- (data.sig + 3) * (1479.8 / 8192);

rml <- round(length(data.sig)/50) * 2 + 1; ## running median length
if(doPlot){
    png("drift.png", width=1280, height=720, pointsize=24);
    par(mar=c(4,4,0.5,0.5));
    if(orig.sig.len > 100000){
        plot((1:length(data.sig))/4000, data.sig, pch=19, cex=0.25,
             xlab="time (s)", ylab="Unadjusted raw signal", col="#80808020");
    } else {
        plot((1:length(data.sig))/4000, data.sig, type="l",
             xlab="time (s)", ylab="Unadjusted raw signal", col="grey");
    }
    if(orig.sig.len > 100000){
        psamp <- seq(1,length(data.sig), length.out=20000);
        points((1:length(data.sig))[psamp]/4000,
               runmed(data.sig, rml, endrule="constant")[psamp],
               type="l", lwd=3, col="black");
    } else {
        points((1:length(data.sig))/4000,
               runmed(data.sig, rml, endrule="constant"),
               type="l", lwd=3, col="black");
    }
}
glm.res <- glm(y ~ x,
               data=data.frame(x=(1:length(data.sig))/4000,
                               y=runmed(data.sig,rml, endrule="constant")));
if(doPlot){
    abline(glm.res, col="#FF000040", lty="dashed", lwd=3);
    text(length(data.sig)/8000, min(data.sig), pos=3,
         sprintf("Running median drift (k=%d): %0.1f units per second",
                 rml, glm.res$coefficients[2]), col="darkred");
} else {
    cat(sprintf("Running median drift (k=%d): %0.1f units per second\n",
                rml, glm.res$coefficients[2]));
}
if(doPlot){
    glm2.res <- glm(y ~ x,
                    data=data.frame(x=(1:length(data.sig))/4000,
                                    y=data.sig));
    abline(glm2.res, col="#0000FF40", lty="dashed", lwd=3);
    text(length(data.sig)/8000, min(data.sig)+dMad, pos=3,
         sprintf("Unadjusted drift: %0.1f units per second",
                 glm2.res$coefficients[2]), col="darkblue");
    dummy <- dev.off();
    if(length(data.sig) < 100000){
        pdf("trimmed.pdf", width=12, height=8, pointsize=16);
        par(mar=c(4.5,6,1,1), mfrow=c(3,1));
        xsplit <- length(data.sig)/(3*4000);
        for(start in 0:2){
            plot((1:length(data.sig)+startPoint)/4000, data.sig, pch=19, type="l",
                 xlim=c(xsplit*start,xsplit*(start+1)) + startPoint/4000,
                 cex=0.25, xlab="", las=1, ylab="");
            if(start == 1){
                mtext("Unadjusted Raw Signal", side=2, line=4, xpd=NA);
            }
            if(start == 2){
                mtext("Time (s)", side=1, line=3, xpd=NA);
            }
        }
        dummy <- dev.off();
    }
}


sampleRate <- 4000;
dRange <- dMax-dMin;

if(doPlot && (length(data.sig) > 100000) && grepl("PDF$", imageName, ignore.case=TRUE)){
    cat("Data length too long for PDF plot, changing to PNG\n");
    imageName <- sub(".pdf$", "_%02d.png", imageName, ignore.case=TRUE);
}

hp.breakPoints <- c(1, hp.breakPoints, length(data.sig));

if(doPlot){
    sigAspect <- (dRange / (length(hp.breakPoints) - 1)) / length(data.sig);
    sw <- 8; ## signal plot width
    sh <- 11; ## signal plot height
    if(grepl("\\.pdf$", imageName)){
        sigLines <- min(15,round(sh / (sigAspect * sw * 2)));
        pdf(imageName, paper="a4", width=sw, height=sh);
    } else if (grepl("\\.png$", imageName)){
        sw <- 1600; ## signal plot width
        sh <- 1600; ## signal plot height
        sigLines <- min(20,round(sh / (sigAspect * sw * 2)));
        png(imageName, width=sw, height=sh, pointsize=24);
    }
    for(sigStart in 2:length(hp.breakPoints)-1){
        ## subselect raw sequence between breakpoints (including hairpin)
        sub.data.sig <- data.sig[hp.breakPoints[sigStart]:hp.breakPoints[sigStart+1]];
        ## recalculate parameters for current segment
        dMed <- median(sub.data.sig);
        dMad <- mad(sub.data.sig);
        dMin <- max(min(sub.data.sig),dMed-4*dMad,0);
        dMax <- min(max(sub.data.sig),dMed+4*dMad,65535);
        dRange <- max(dMed-dMin,dMax-dMed)*2;
        sub.data.sig <- sub.data.sig - dMed + dRange/2;
        sigAspect <- dRange / length(sub.data.sig);
        sigLines <- min(15,round(sh / (sigAspect * sw * 2)));
        par(mar=c(0.5,0.5,0.5,0.5));
        width <- ceiling(length(sub.data.sig) / sigLines);
        plot(NA, xlim=c(0,width), ylim=c(0,dRange*sigLines),
             axes=FALSE, ann=FALSE);
        for(x in 1:sigLines){
            startPoint <- (x-1) * width;
            yPoints <- if(startPoint == 0){
                           head(sub.data.sig, width);
                       } else {
                           head(tail(sub.data.sig,-startPoint),width);
                       }
            if(length(yPoints) > 20000){
                points(x=1:length(yPoints),
                       y=yPoints + (sigLines - x) * dRange, pch=".");
            } else {
                points(x=1:length(yPoints),
                       y=yPoints + (sigLines - x) * dRange, type="l");
            }
            segments(x0=1, x1=length(yPoints), y0=(sigLines - x + 0.5) * dRange,
                     col="#6495EDA0", lwd=3);
            segments(x0=1, x1=length(yPoints), y0=(sigLines - x + 0.5) * dRange,
                     col="#FFF8DCA0", lwd=1);
            for(t in seq(1,width,length.out=5)){
                tVal=round((startPoint+hp.breakPoints[sigStart]+t-2) / 4000,2);
                cVal=floor((tVal*100) %% 100);
                text(t,(sigLines - x) * dRange, tVal, col=rainbow(100)[cVal+1],
                     adj=ifelse(t==1,0,ifelse(t==width,1,0.5)), cex=0.71);
            }
        }
    }
    dummy <- dev.off();
    cat(sprintf("Done... written to '%s'\n", imageName));
}
