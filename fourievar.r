#!/usr/bin/Rscript

if(length(commandArgs(TRUE)) < 1){
    cat("Error: no fasta file specified\n");
    cat("syntax: ./fouriervar.r <input.fa>\n");
    quit(save="no");
}

library(Biostrings, quietly=TRUE, warn.conflicts=FALSE, verbose=FALSE);

fuzz <- 2;
doPlot <- FALSE;

#fileName <- "/bioinf/presentations/2017-Sep-03/tig00022132_rpt0.fa";
fileName <- commandArgs(TRUE)[1];
seqs <- as.vector(readDNAStringSet(fileName)[[1]]);
slen <- length(seqs);
fLimit <- 1000;
fLimit <- min(slen-1, fLimit);

baseCplx <- complex(real=c(A=1, C=0, G=-1, T=0)[seqs],
                    imaginary=c(A=0, C=1, G=0, T=-1)[seqs]);

spectrum <- matrix(0, slen, fLimit);
for(cycle in 2:fLimit){
    spectrum[1:(slen-cycle), cycle] <-
        (head(baseCplx, -cycle) == tail(baseCplx, -cycle)) + 0;
}

colnames(spectrum) <- 1:ncol(spectrum);

## cat("Removing noise\n");
## ## Remove transient repeats
## spectrum[spectrum < 4] <- 0;
## binSpec <- sign(spectrum);
## filtSpec <- apply(binSpec, 2, function(x){
##     rx <- rle(x);
##     rx$values[rx$lengths <= (fuzz*3)] <- 0;
##     rx$values <- (rx$values == 1);
##     inverse.rle(rx);
## });

## spectrum[!filtSpec] <- 0;
## spectrum[,c(1:2)] <- NA;

scols <- colSums(sign(spectrum)) / (slen - 1:ncol(spectrum));
cat(fileName, which(scols == max(scols, na.rm=TRUE)), max(scols, na.rm=TRUE),
    median(scols, na.rm=TRUE), mad(scols, na.rm=TRUE), "\n");

scols[scols==0] <- NA;

if(doPlot){
    cat("Writing to file\n");
    pdf(paste0("DNA_spectrogram_",sub("\\.[^\\.]$","",fileName),".pdf"),
        height=8, width=11);
                                        #pdf(paste0("DNA_spectrogram_test.pdf"),
                                        #    height=8, width=11);
    image(x=1:length(seqs), y=1:ncol(spectrum), spectrum,
          col=(colorRampPalette(c("white",blues9)))(100),
          xlab="Sequence location (bp)",
          ylab="Repeat width (bp)", useRaster=TRUE);
    plot(x=1:ncol(spectrum), scols);
    abline(h=median(scols, na.rm=TRUE)+(mad(scols, na.rm=TRUE)*3*c(-1,1)));
    invisible(dev.off());
}
