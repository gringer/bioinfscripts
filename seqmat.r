#!/usr/bin/Rscript

fileName <- commandArgs(TRUE)[1];
rptSize <- as.numeric(commandArgs(TRUE)[2]);

inLines <- readLines(fileName);
inLines[2] <- paste(inLines[-1],collapse="");

inName <- substring(inLines[1],2);
inSeq <- c(A=1, C=2, G=3, T=4)[unlist(strsplit(inLines[2],""))];

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
