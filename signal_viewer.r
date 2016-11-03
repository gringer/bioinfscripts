#!/usr/bin/Rscript
## expects a file containing only raw signal data, as produced by
## fast5extractor.py
myFileName <- "out.bin";
channel <- 389;
read <- 490;
createImage <- FALSE;
imageName <- "";

if(length(commandArgs(TRUE)) > 0){
    myFileName <- commandArgs(TRUE)[1];
    channel <- as.numeric(commandArgs(TRUE)[2]);
    read <- as.numeric(commandArgs(TRUE)[3]);
}

if(length(commandArgs(TRUE)) > 3){
    imageName <- commandArgs(TRUE)[4];
    createImage <- TRUE;
}

fileLen <- file.size(myFileName);
data.sig <- readBin(myFileName, what=integer(), size=2, signed=FALSE,
                    n=fileLen/2);

if(length(data.sig) < 4000){
    quit(save="no");
}

dMed <- median(data.sig);
dMad <- mad(data.sig);
dMin <- max(min(data.sig),dMed-4*dMad,0);
dMax <- min(max(data.sig),dMed+4*dMad,65535);
rangeRLE <- rle((runmed(data.sig,11) > dMin) & (runmed(data.sig,11) < dMax));
if(length(rangeRLE$lengths) > 1){
    startPoint <- head(tail(cumsum(rangeRLE$lengths),2),1) + 50;
    data.sig <- tail(data.sig, -startPoint);
    if(length(data.sig) > 1){
        dMed <- median(data.sig);
        dMad <- mad(data.sig);
        dMin <- max(min(data.sig),dMed-4*dMad,0);
        dMax <- min(max(data.sig),dMed+4*dMad,65535);
    }
}

if(length(data.sig) < 8000){
    quit(save="no");
}

sampleRate <- 4000;
dRange <- dMax-dMin;
#X11(width=14.25, height=7.5);

if(createImage){
    png(imageName);
} else {
    X11(width=7, height=5);
}
par(mar=c(1,1,2,1));
startPTM <- proc.time()["elapsed"];
sp <- 0;
while(sp < (length(data.sig) - 8000)){
    sp <- (proc.time()["elapsed"] - startPTM) * sampleRate;
    plot(NA, ylim=c(0,(dRange)*3), xlim=c(0,4000),
         ann=FALSE, axes=FALSE);
    mtext(sprintf("Raw Signal from Channel %d, Read %d", channel, read),
          side=3, cex=2);
    points(data.sig[(sp+00001):(sp+04000)]-dMin+dRange*2, type="l");
    points(data.sig[(sp+04001):(sp+08000)]-dMin+dRange*1, type="l");
    points(data.sig[(sp+08001):(sp+12000)]-dMin+dRange*0, type="l");
    spTLabs <- floor(sp/(4000));
    spTLabs <- spTLabs+(seq(0.25, 6, by=0.25));
    spTLabs <- spTLabs[(spTLabs * sampleRate) < length(data.sig)];
    spTLposSamp <- (spTLabs*4000 - sp);
    spTLabs <- spTLabs[(spTLposSamp > 0) & (spTLposSamp < 12000)]
    spTLposSamp <- spTLposSamp[(spTLposSamp > 0) & (spTLposSamp < 12000)]
    spTLposX <- (spTLposSamp) %% 4000;
    spTLposY <- (2-floor((spTLabs*4000 - sp) / 4000)) * dRange;
    text(x=spTLposX, y=spTLposY, labels=spTLabs,
         col=rainbow(10, alpha=0.5)[(spTLabs - 1) %% 10 + 1], cex=2);
    Sys.sleep(1/10);
}
Sys.sleep(2);
dummy <- dev.off();

