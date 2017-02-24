#!/usr/bin/Rscript

rptSize <- as.numeric(commandArgs(TRUE)[1]);

inLines <- readLines("stdin");
inName <- substring(inLines[1],2);
inSeq <- c(A=1, C=2, G=3, T=4)[unlist(strsplit(inLines[2],""))];

png("sequence_matrix.png", width=1280, height=1280,
    pointsize=24);
subSeq <- inSeq[1:(floor(length(inSeq)/rptSize)*rptSize)];
par(mar=c(0.5,0.5,1,0.5));
image(matrix(subSeq,nrow=rptSize),
      main=sprintf("%s (%0.3f kb, %d bp repeat)",inName,
                   length(inSeq)/1000, rptSize),
      cex.main=0.8, xaxt="n", yaxt="n", useRaster=TRUE,
      col=c("green","blue","yellow","red"));
dummy <- dev.off();
