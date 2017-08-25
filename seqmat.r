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
      col=c("#006400","#0000FF","#FFD700","#FF6347","lightgrey"));
axis(2, at=round(seq(0, numLines, length.out=20)),
     labels=round(seq(0, numLines, length.out=20)) * rptSize+1,
     las=2, cex.axis=1);
dummy <- dev.off();

numLoops <- (length(inSeq) / rptSize);
##startCount <- rptSize / 1.2;
startCount <- rptSize;
startRadius <- ifelse(numLoops < 10, 0.5, 0.2);
endRadius <- 1.0;
##loopIncrement <- ((rptSize * 1.2) - (rptSize / 1.2)) / numLoops;


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
## integrate(2*pi*r,r=startRadius..endRadius)
## => pi((endRadius)²-(startRadius)²)
dTot <- pi*(endRadius^2 - startRadius^2); ## total "distance" travelled
theta <- seq(0, numLoops * 2*pi, length.out=length(inSeq)); ## traversed angle
deg <- (theta / (2*pi)) * 360;
r <- seq(sqrt(startRadius), sqrt(endRadius),
         length.out=length(inSeq))^2; ## path radius
s <- pi * (r^2 - startRadius^2); ## traversed distance at each pos
ds <- c(s[2],diff(s)); ## distance difference at each pos

## par(mfrow=c(4,1));
## plot(s, main="s");
## plot(r, main="r");
## plot(deg, main="deg");
## plot(theta, main="theta");

## draw the spiral
points(rev(r) * cos(rev(theta)), rev(r) * sin(rev(theta)),
       col=c("#00640080","#0000FF80","#FFD70080","#FF634780")[rev(inSeq)],
       pch=16, cex=rev(sqrt(r)) * (7/log(numLoops)));
legend("center", legend=c("A","C","G","T"), inset=0.2,
       fill=c("#006400","#0000FF","#FFD700","#FF6347"),
       cex=ifelse(numLoops < 10, 1, 0.71));
invisible(dev.off());
