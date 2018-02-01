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
seqs <- readDNAStringSet(fileName);

print(seqs);

for(si in 1:length(seqs)){
    mySeq <- as.vector(seqs[[si]]);
    slen <- length(mySeq);
    fLimit <- 1000;
    fLimit <- min(slen-1, fLimit);
    
    baseCplx <- complex(real=c(A=1, C=0, G=-1, T=0)[mySeq],
                        imaginary=c(A=0, C=1, G=0, T=-1)[mySeq]);
    
    spectrum <- sapply(2:fLimit, function(cycle){
        sum(head(baseCplx, -cycle) == tail(baseCplx, -cycle)) / (slen-cycle);
    });

    png("out.png");
    plot(spectrum);
    invisible(dev.off());

    names(spectrum) <- 2:fLimit;
    spectrum <- (spectrum - mean(spectrum)) / sd(spectrum);
    spectrum <- spectrum[order(-spectrum)];

    print(head(cbind(spectrum, qProb = -10 * log10(1-pnorm(spectrum))),20));
}
