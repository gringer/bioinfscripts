#!/usr/bin/Rscript

if(length(commandArgs(TRUE)) < 1){
    cat("Error: no fasta file specified\n");
    cat("syntax: ./fouriervar.r <input.fa>\n");
    quit(save="no");
}

library(Biostrings, quietly=TRUE, warn.conflicts=FALSE, verbose=FALSE);
library(bitops, quietly=TRUE, warn.conflicts=FALSE, verbose=FALSE);
library(caTools, quietly=TRUE, warn.conflicts=FALSE, verbose=FALSE);

fuzz <- 2;

fileName <- commandArgs(TRUE)[1];
seqs <- readDNAStringSet(fileName)[[1]];
slen <- length(seqs);
fLimit <- 1000;

baseBin <- sapply(c("A","C","G","T"), function(x){c((as.vector(seqs) == x)+0,rep(0,fuzz))});
basePoss <-
    sapply(c("A","C","G","T"), function(x){
        res <- which(as.vector(seqs) == x);
    });

spectrum <- matrix(0, slen, fLimit);
for(base in c("A","C","G","T")){
    cat("Base:",base,"\n");
    possB <- basePoss[[base]];
    poss <- possB;
    for(cycle in 1:ncol(spectrum)){
        poss <- possB[possB > (cycle + fuzz)] - cycle - fuzz;
        spectrum[poss, cycle] <-
            baseBin[poss + 0, base] *  1 +
            baseBin[poss + 1, base] *  4 +
            baseBin[poss + 2, base] * 16 +
            baseBin[poss + 3, base] *  4 +
            baseBin[poss + 4, base] *  1;
    }
}

colnames(spectrum) <- 1:ncol(spectrum);

cat("Removing noise\n");
## Remove transient repeats
spectrum[spectrum < 5] <- 0;
binSpec <- sign(spectrum);
filtSpec <- apply(binSpec, 2, function(x){
    rx <- rle(x);
    rx$values[rx$lengths <= (fuzz*3)] <- 0;
    rx$values <- (rx$values == 1);
    inverse.rle(rx);
});

spectrum[!filtSpec] <- 0;

cat("Writing to file\n");

pdf(paste0("DNA_spectrogram_",sub("\\.[^\\.]$","",fileName),".pdf"),
    height=8, width=11);
image(x=1:length(seqs), y=1:ncol(spectrum), spectrum,
      col=(colorRampPalette(c("white",blues9)))(100),
      xlab="Sequence location (bp)",
      ylab="Repeat width (bp)", useRaster=TRUE);
abline(h=seq(0,fLimit, by=10), lty="dashed", lwd=1, col="#FF5000A0");
invisible(dev.off());
