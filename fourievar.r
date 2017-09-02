#!/usr/bin/Rscript

library(Biostrings);
library(bitops);
library(caTools);

seqs <- readDNAStringSet("tig00022132_rpt0.fa")[[1]];
slen <- length(seqs);

baseBin <- sapply(c("A","C","G","T"), function(x){(as.vector(seqs) == x)+0});
basePoss <-
    sapply(c("A","C","G","T"), function(x){which(as.vector(seqs) == x)});

fLimit <- 300;

spectrum <- matrix(0, slen, fLimit);
colnames(spectrum) <- 1:ncol(spectrum);
for(cycle in 1:ncol(spectrum)){
    cat(cycle);
    for(base in c("A","C","G","T")){
        for(pshift in -2:2){
            poss <- basePoss[[base]] - cycle + pshift;
            poss <- poss[(poss >= 1) & (poss <= slen)];
            spectrum[poss,cycle] <- spectrum[poss,cycle] +
                bitShiftL(baseBin[poss,base], 4-abs(pshift*2));
        }
    }
}

smoothSpec <- apply(spectrum, 2, runmean, 20);

png("Spectral_decomposition.png", height=1800, width=1800, pointsize=36);
image(x=1:length(seqs), y=1:ncol(spectrum), spectrum,
      col=(colorRampPalette(c("white",blues9)))(100),
      xlab="Sequence location (bp)",
      ylab="Repeat width (bp)");
invisible(dev.off());
