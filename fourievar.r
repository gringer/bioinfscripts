#!/usr/bin/Rscript

if(length(commandArgs(TRUE)) < 1){
    cat("Error: no fasta file specified\n");
    cat("syntax: ./fouriervar.r <input.fa>\n");
    quit(save="no");
}

library(Biostrings, quietly=TRUE, warn.conflicts=FALSE, verbose=FALSE);
library(bitops, quietly=TRUE, warn.conflicts=FALSE, verbose=FALSE);
library(caTools, quietly=TRUE, warn.conflicts=FALSE, verbose=FALSE);

fileName <- commandArgs(TRUE)[1];
seqs <- readDNAStringSet(fileName)[[1]];
slen <- length(seqs);

baseBin <- sapply(c("A","C","G","T"), function(x){(as.vector(seqs) == x)+0});
basePoss <-
    sapply(c("A","C","G","T"), function(x){which(as.vector(seqs) == x)});

fLimit <- 600;

spectrum <- matrix(0, slen, fLimit);
colnames(spectrum) <- 1:ncol(spectrum);
for(cycle in 1:ncol(spectrum)){
    for(base in c("A","C","G","T")){
        for(pshift in -2:2){
            poss <- basePoss[[base]] - cycle + pshift;
            poss <- poss[(poss >= 1) & (poss <= slen)];
            spectrum[poss,cycle] <- spectrum[poss,cycle] +
                bitShiftL(baseBin[poss,base], 4-abs(pshift*2));
        }
    }
}

## Remove transient repeats
spectrum[spectrum < 8] <- 0;

##filtSpec <- apply(spectrum[10:11,], 1, function(x){
##    rx <- rle(x);
##    rs$values[rs$values < 10];
##    rs$values[rs$lengths == 1];
##    print(rle(x));
##});

smoothSpec <- apply(spectrum, 2, runmean, 20);

png(paste0("DNA_spectrogram_",sub("\\.[^.]$","",fileName),".png"),
    height=1800, width=1800, pointsize=36);
image(x=1:length(seqs), y=1:ncol(spectrum), spectrum,
      col=(colorRampPalette(c("white",blues9)))(100),
      xlab="Sequence location (bp)",
      ylab="Repeat width (bp)");
abline(h=seq(0,fLimit, by=10), lty="dashed", lwd=2, col="#FF5000A0");
invisible(dev.off());
