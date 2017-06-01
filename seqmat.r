#!/usr/bin/Rscript

fileName <- commandArgs(TRUE)[1];
rptSize <- as.numeric(commandArgs(TRUE)[2]);

inLines <- readLines(fileName);

inName <- substring(inLines[1],2);
inSeq <- c(A=1, C=2, G=3, T=4)[unlist(strsplit(paste(inLines[-1],collapse=""),""))];

type <- "png";
if(any(commandArgs(TRUE) == "pdf")){
    type <- "pdf";
}

if(type == "png"){
    png("sequence_matrix.png", width=1200, height=1200, pointsize=24);
} else if(type == "pdf"){
    pdf("sequence_matrix.pdf", width=12, height=12, pointsize=16);
}
lis <- length(inSeq);
numLines <- floor(lis/rptSize + 1);
inSeq <- c(inSeq,rep(5,rptSize));
subSeq <- inSeq[1:(numLines*rptSize)];
par(mar=c(0.5,5,1,0.5), mgp=c(3.5,1,0));
image(x=1:rptSize, y=1:numLines-1, matrix(subSeq,nrow=rptSize),
      main=sprintf("%s (%0.3f kb, %d bases / line)", sub(" .*$","",inName),
                   lis/1000, rptSize), ylab="Base location",
      cex.main=0.8, xaxt="n", yaxt="n", useRaster=TRUE,
      col=c("green","blue","yellow","red","lightgrey"));
axis(2, at=round(seq(0, numLines, length.out=20)),
     labels=round(seq(0, numLines, length.out=20)) * rptSize+1,
     las=2, cex.axis=1);
dummy <- dev.off();

numLoops <- ceiling(length(inSeq) / rptSize);
##startCount <- rptSize / 1.2;
startCount <- rptSize;
startRadius <- 0.2;
endRadius <- 1.0;
##loopIncrement <- ((rptSize * 1.2) - (rptSize / 1.2)) / numLoops;
loopIncrement <- 0;


if(type == "png"){
    png("sequence_circle.png", width=1200, height=1200, pointsize=24);
} else if(type == "pdf"){
    pdf("sequence_circle.pdf", width=12, height=12, pointsize=16);
}
par(mar=c(3,3,3,3));
plot(NA, xlim=c(-1,1), ylim=c(-1,1), axes=FALSE, ann=FALSE);
mtext(sprintf("%s (%0.3f kb, %d bases / ring)", sub(" .*$","",inName),
              lis/1000, rptSize));
## Pre-population plot variables
radiusFactor <- (endRadius - startRadius) / numLoops;
radius <- unlist(sapply(1:numLoops,function(x){
  loopCount <- startCount + loopIncrement * (x-1);
  ((x-1) + seq(0, 1-(1/loopCount), length.out=loopCount-1)) *
    radiusFactor + startRadius}));
angle <- unlist(sapply(1:numLoops,function(x){
  loopCount <- startCount + loopIncrement * (x-1);
  seq(0, 2*pi-(2*pi/loopCount), length.out=loopCount-1)}));
## draw the spiral
points(rev(radius) * cos(rev(angle)), rev(radius) * sin(rev(angle)),
       col=c("#00FF0060","#0000FF60","#FFFF0060","#FF000060")[rev(inSeq)],
       pch=16, cex=rev(sqrt(radius)) * 2);
dummy <- dev.off();
