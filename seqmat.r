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
subSeq <- inSeq[1:(floor(length(inSeq)/rptSize)*rptSize)];
par(mar=c(0.5,0.5,1,0.5));
image(matrix(subSeq,nrow=rptSize),
      main=sprintf("%s (%0.3f kb, %d bases / line)",inName,
                   length(inSeq)/1000, rptSize),
      cex.main=0.8, xaxt="n", yaxt="n", useRaster=TRUE,
      col=c("green","blue","yellow","red"));
dummy <- dev.off();
