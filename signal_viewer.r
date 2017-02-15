#!/usr/bin/Rscript
## expects a file containing only raw signal data, as produced by
## 'porejuicer.py raw <file.fast5>'

## Default parameters
sigFileName <- "signal.bin";
channel <- -1;
read <- -1;
imageName <- "signal_out.pdf";
doPlot <- TRUE;

usage <- function(){
  cat("usage: ./signal_viewer.r",
      "<signal.bin> [options]\n");
  cat("\nOther Options:\n");
  cat("-out <file> : Write image to <file> [default: signal_out.pdf]\n");
  cat("-noplot     : Don't plot any data (just carry out calculations)\n");
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

if(doPlot){
    png("untrimmed.png", width=1280, height=720, pointsize=16);
    plot(data.sig, type="l", main="Untrimmed raw signal", xlab="Samples",
         ylab="Unadjusted raw signal");
    dummy <- dev.off();
}

dMed <- median(data.sig);
dMad <- mad(data.sig);
dMin <- max(min(data.sig),dMed-4*dMad,0);
dMax <- min(max(data.sig),dMed+4*dMad,65535);

rangeRLE <- rle((runmed(data.sig,11) > dMin) & (runmed(data.sig,11) < dMax));
print(rangeRLE);
if(length(rangeRLE$lengths) > 1){
    ## subset on longest TRUE region
    indexRLE <- order(-rangeRLE$values, -rangeRLE$lengths)[1];
    print(indexRLE);
    startPoint <- ifelse(indexRLE == 1, 1,
                         cumsum(rangeRLE$lengths)[startPoint-1] + 50);
    data.sig <- tail(data.sig, -startPoint);
    data.sig <- head(data.sig, rangeRLE$lengths[indexRLE]);
    if(length(data.sig) > 1){
        dMed <- median(data.sig);
        dMad <- mad(data.sig);
        dMin <- max(min(data.sig),dMed-4*dMad,0);
        dMax <- min(max(data.sig),dMed+4*dMad,65535);
    }
} else { ## trim off 5 samples to deal with short-length initial peaks
    data.sig <- tail(data.sig, -5);
}

## filter out huge signal spikes
data.sig[data.sig > dMax] <- dMax;

## data.sig <- (data.sig + 3) * (1479.8 / 8192);

if(length(data.sig) == 0){
    cat("Warning: no signal data found after noise trimming\n");
    quit(save="no", status=1);
}

rml <- round(length(data.sig)/50) * 2 + 1; ## running median length
if(doPlot){
    png("drift.png", width=1280, height=720, pointsize=24);
    par(mar=c(4,4,0.5,0.5));
    plot((1:length(data.sig))/4000, data.sig, type="l",
         xlab="time (s)", ylab="Unadjusted raw signal", col="grey");
    points((1:length(data.sig))/4000, runmed(data.sig, rml, endrule="constant"),
           type="l", lwd=3, col="black");
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
glm2.res <- glm(y ~ x,
    data=data.frame(x=(1:length(data.sig))/4000,
                    y=data.sig));
if(doPlot){
    abline(glm2.res, col="#0000FF40", lty="dashed", lwd=3);
    text(length(data.sig)/8000, min(data.sig)+dMad/2, pos=3,
         sprintf("Unadjusted drift: %0.1f units per second",
                 glm2.res$coefficients[2]), col="darkblue");
    dummy <- dev.off();
} else {
    cat(sprintf("Unadjusted drift: %0.1f units per second\n",
                glm2.res$coefficients[2]));
}


sampleRate <- 4000;
dRange <- dMax-dMin;
data.sig <- data.sig - dMin;

sigAspect <- dRange / length(data.sig);

sw <- 8; ## signal plot width
sh <- 11; ## signal plot height

if(doPlot){
    if(grepl("\\.pdf$", imageName)){
        sigLines <- min(20,round(sh / (sigAspect * sw * 2)));
        pdf(imageName, paper="a4", width=sw, height=sh);
    } else if (grepl("\\.png$", imageName)){
        sw <- 1600; ## signal plot width
        sh <- 1600; ## signal plot height
        sigLines <- min(20,round(sh / (sigAspect * sw * 2)));
        png(imageName, width=sw, height=sh, pointsize=24);
    }
    par(mar=c(0.5,0.5,0.5,0.5));
    width <- ceiling(length(data.sig) / sigLines);
    plot(NA, xlim=c(0,width), ylim=c(0,dRange*sigLines),
         axes=FALSE, ann=FALSE);
    for(x in 1:sigLines){
        startPoint <- (x-1) * width;
        yPoints <- if(startPoint == 0){
                       head(data.sig, width);
                   } else {
                       head(tail(data.sig,-startPoint),width);
                   }
        points(x=1:length(yPoints),
               y=yPoints + (sigLines - x) * dRange, type="l");
        for(t in seq(1,width,length.out=5)){
            tVal=round((startPoint+t-1) / 4000,2);
            cVal=floor((tVal*10) %% 10);
            text(t,(sigLines - x) * dRange, tVal, col=rainbow(10)[cVal+1],
                 adj=ifelse(t==1,0,ifelse(t==width,1,0.5)), cex=0.71);
        }
    }
    dummy <- dev.off();
    cat(sprintf("Done... written to '%s'\n", imageName));
}

