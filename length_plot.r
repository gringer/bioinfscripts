#!/usr/bin/Rscript

fileNames <- commandArgs(TRUE);

if(length(fileNames) == 0){
    fileNames <- list.files(pattern="lengths_.*\\.txt(\\.gz)?");
}

valToSci <- function(val, unit = ""){
    sci.prefixes <- c("", "k", "M", "G", "T", "P", "E", "Z", "Y");
    units <- rep(paste(sci.prefixes,unit,sep=""), each=3);
    logRegion <- floor(log10(val))+1;
    conv.units <- units[logRegion];
    conv.div <- 10^rep(0:(length(sci.prefixes)-1) * 3, each = 3)[logRegion];
    conv.val <- val / conv.div;
    conv.val[val == 0] <- 0;
    conv.units[val == 0] <- unit;
    return(sprintf("%s %s",conv.val,conv.units));
}

sequence.hist <- function(lengths, invert = TRUE, ...){
    fib.divs <- round(10^((0:4)/5) * 2) * 0.5; ## splits log decades into 5
    histBreaks <- round(rep(10^(0:16),each=5) * fib.divs);
    lengthRange <- range(lengths);
    ## filter on actual data range
    histBreaks <- histBreaks[(which.min(histBreaks < lengthRange[1])-1):
                             which.max(histBreaks > lengthRange[2])];
    seqd.bases <- tapply(lengths,cut(lengths, breaks=histBreaks), sum);
    seqd.counts <- tapply(lengths,cut(lengths, breaks=histBreaks), length);
    xBreaks <- round(rep(10^(0:16),each=5) * fib.divs);
    axRange <- range(seqd.bases);
    xBreaks <- xBreaks[(which.min(xBreaks < axRange[1])):
                       (which.max(xBreaks > axRange[2])-1)];
    xBreaksMajor <- xBreaks[log10(xBreaks) - floor(log10(xBreaks)) < 0.001];
    xBreaksMinor <- rep(xBreaksMajor,each=9) * 1:9;
    xBreaksMinor <- xBreaksMinor[(which.min(xBreaksMinor < axRange[1])):
                                 (which.max(xBreaksMinor > axRange[2])-1)];
    barPos <- barplot(if(invert){rev(seqd.bases)} else {seqd.bases},
                      log = "x", las = 1, axes = FALSE, col = "steelblue",
                      horiz = TRUE, names.arg = rep("",length(seqd.bases)),
                      ylab = "",
                      xlab = "Aggregate sequence length (number of sequences)",
                      ...);
    barGap <- diff(barPos)[1];
    barOffset <- barPos[1] - barGap/2;
    axis(2, at = if(invert){rev(seq(barOffset,by=barGap,
                                    length.out = length(histBreaks)))}
         else {seq(barOffset,by=barGap,length.out = length(histBreaks))},
         labels = valToSci(histBreaks,"b"), las = 2);
    mtext("Fragment size", side=2, line=5);
    axis(1, at = xBreaksMajor, labels = valToSci(xBreaksMajor, "b"), lwd=3);
    axis(1, at = xBreaksMinor, labels = FALSE);
    text.poss <- ((log10(seqd.bases) < mean(par("usr")[1:2]))+1)*2;
    text.poss[is.na(text.poss)] <- 4;
    text.col <- c("white","black")[((log10(seqd.bases) <
                                     mean(par("usr")[1:2]))+1)];
    text(seqd.bases,if(invert){rev(barPos)} else {barPos},
         paste(valToSci(signif(seqd.bases,4)),
               " (", seqd.counts, ")", sep = ""),
         pos=text.poss, col=text.col, cex = 0.65);
}

plain.hist <- function(lengths, invert = TRUE, ...){
    fib.divs <- round(10^((0:4)/5) * 2) * 0.5; ## splits log decades into 5
    histBreaks <- round(rep(10^(0:16),each=5) * fib.divs);
    lengthRange <- range(lengths);
    ## filter on actual data range
    histBreaks <- histBreaks[(which.min(histBreaks < lengthRange[1])-1):
                             which.max(histBreaks > lengthRange[2])];
    seqd.bases <- tapply(lengths,cut(lengths, breaks=histBreaks), sum);
    seqd.counts <- tapply(lengths,cut(lengths, breaks=histBreaks), length);
    xBreaks <- round(rep(10^(0:16),each=5) * fib.divs);
    axRange <- range(seqd.counts);
    xBreaks <- xBreaks[(which.min(xBreaks < axRange[1])):
                       (which.max(xBreaks > axRange[2])-1)];
    xBreaksMajor <- xBreaks[log10(xBreaks) - floor(log10(xBreaks)) < 0.001];
    xBreaksMinor <- rep(xBreaksMajor,each=9) * 1:9;
    xBreaksMinor <- xBreaksMinor[(which.min(xBreaksMinor < min(xBreaksMajor))):
                                 (which.max(xBreaksMinor > max(xBreaksMajor))-1)];
    barPos <- barplot(log10(if(invert){rev(seqd.counts)} else {seqd.counts})+0.025,
                      las = 1, axes = FALSE, col = "steelblue",
                      horiz = TRUE, names.arg = rep("",length(seqd.counts)),
                      ylab = "", xlim=c(log10(max(seqd.counts, na.rm=TRUE))+0.025,
                                        log10(min(seqd.counts, na.rm=TRUE))-0.025),
                      xlab = "Number of sequences (Aggregate length)",
                      ...);
    barGap <- diff(barPos)[1];
    barOffset <- barPos[1] - barGap/2;
    axis(4, at = if(invert){rev(seq(barOffset,by=barGap,
                                    length.out = length(histBreaks)))}
         else {seq(barOffset,by=barGap,length.out = length(histBreaks))},
         labels = valToSci(histBreaks,"b"), las = 2, pos=0);
    mtext("Fragment size", side=4, line=5);
    axis(1, at = log10(xBreaksMajor)+0.025, labels = valToSci(xBreaksMajor), lwd=3);
    axis(1, at = log10(xBreaksMinor)+0.025, labels = FALSE);
    text.poss <- ((log10(seqd.counts) > mean(par("usr")[1:2]))+1)*2;
    text.poss[is.na(text.poss)] <- 2;
    text.col <- c("white","black")[((log10(seqd.counts) <
                                     mean(par("usr")[1:2]))+1)];
    text(log10(seqd.counts)+0.025,if(invert){rev(barPos)} else {barPos},
         paste(seqd.counts,
               " (", valToSci(signif(seqd.bases,4),"b"), ")", sep = ""),
         pos=text.poss, col=text.col, cex = 0.65);
}


pdf("MinION_Reads_SequenceHist.pdf", paper="a4r",
    width=11, height=8);
par(mar=c(5.5,6.5,2.5,1.5));
dens.mat <- sapply(fileNames, function(x){
    cat(x,"...");
    data <- scan(x, comment.char=" ", quiet=TRUE);
    cat(" done\n");
    res <- density(log10(data), from=2, to=6, bw=0.1);
    res.out <- res$y;
    names(res.out) <- round(res$x,3);
    subName <- sub("lengths_(.*)\\.txt(\\.gz)?","\\1", x);
    png(sprintf("MinION_Reads_SequenceHist_%s.png",
                subName), pointsize=24,
        width=1280, height=960);
    par(mar=c(5.5,6.5,2.5,1.5));
    sequence.hist(data,
                  main=sprintf("Read Length Distribution Plot (%s)",subName));
    dummy <- dev.off();
    png(sprintf("MinION_Reads_PlainHist_%s.png",
                subName), pointsize=24,
        width=1280, height=960);
    par(mar=c(5.5,1.5,2.5,6.5));
    plain.hist(data,
               main=sprintf("Read Count Distribution Plot (%s)",subName));
    dummy <- dev.off();
    sequence.hist(data,
                  main=sprintf("Read Length Distribution Plot (%s)",subName));
    plain.hist(data,
               main=sprintf("Read Count Distribution Plot (%s)",subName));
    res.out;
});
dummy <- dev.off();
colnames(dens.mat) <- sub("lengths_(.*)\\.txt(\\.gz)?","\\1",colnames(dens.mat));

bpdens.mat <- dens.mat * 10^as.numeric(rownames(dens.mat));

{
    png("MinION_Bases_DigElec_white.png", pointsize=24,
        width=240*ncol(bpdens.mat), height=960);
    par(mar=c(3,5,0.5,0.5), mgp=c(4,1,0));
    image(x=seq_len(ncol(bpdens.mat)), ann=TRUE, axes=FALSE,
          y=as.numeric(rownames(bpdens.mat)),
          z=t(bpdens.mat), col=colorRampPalette(hsv(h=27/360,s=0,v=seq(1,0,by=-0.001)), bias=1)(100),
          ylab = "Read Length");
    abline(v=seq_len(ncol(bpdens.mat)+1)-0.5);
    abline(h=log10(c(1,2,5)) + rep(0:5, each=3),
           lty="dashed", col="#00000020");
    axis(1,at=seq_len(length(fileNames)), labels=colnames(bpdens.mat),
         lwd=0);
    axis(2,at=log10(c(1,2,5)) + rep(0:5, each=3), las=2,
         labels=valToSci(as.numeric(paste0(c(1,2,5),
             rep(substring("00000",first=0,last=0:5),each=3)))));
    dummy <- dev.off();
}

{
    png("MinION_Reads_DigElec_white.png", pointsize=24,
        width=240*ncol(dens.mat), height=960);
    par(mar=c(3,5,0.5,0.5), mgp=c(4,1,0));
    image(x=seq_len(ncol(dens.mat)), ann=TRUE, axes=FALSE,
          y=as.numeric(rownames(dens.mat)),
          z=t(dens.mat), col=colorRampPalette(hsv(h=27/360,s=0,v=seq(1,0,by=-0.001)), bias=1)(100),
          ylab = "Read Length");
    abline(v=seq_len(ncol(dens.mat)+1)-0.5);
    abline(h=log10(c(1,2,5)) + rep(0:5, each=3),
           lty="dashed", col="#00000020");
    axis(1,at=seq_len(length(fileNames)),
         labels=colnames(dens.mat), lwd=0);
    axis(2,at=log10(c(1,2,5)) + rep(0:5, each=3), las=2,
         labels=valToSci(as.numeric(paste0(c(1,2,5),
             rep(substring("00000",first=0,last=0:5),each=3)))));
    dummy <- dev.off();
}

{
    png("MinION_Bases_DigElec_black.png", pointsize=24,
        width=240*ncol(bpdens.mat), height=960);
    par(mar=c(3,5,0.5,0.5), mgp=c(4,1,0));
    image(x=seq_len(ncol(bpdens.mat)), ann=TRUE, axes=FALSE,
          y=as.numeric(rownames(bpdens.mat)),
          z=t(bpdens.mat), col=colorRampPalette(hsv(h=27/360,s=1,v=seq(0,1,by=0.001)), bias=1.25)(100),
          ylab = "Read Length");
    abline(v=seq_len(ncol(bpdens.mat)+1)-0.5, lwd=3);
    abline(h=log10(c(1,2,5)) + rep(0:5, each=3),
           lty="dashed", col="#FFFFFF40");
    axis(1,at=seq_len(length(fileNames)),
         labels=colnames(bpdens.mat),
         lwd=0);
    axis(2,at=log10(c(1,2,5)) + rep(0:5, each=3), las=2,
         labels=valToSci(as.numeric(paste0(c(1,2,5),
             rep(substring("00000",first=0,last=0:5),each=3)))));
    dummy <- dev.off();
}

{
    png("MinION_Reads_DigElec_black.png", pointsize=24,
        width=240*ncol(dens.mat), height=960);
    par(mar=c(3,5,0.5,0.5), mgp=c(4,1,0));
    image(x=seq_len(ncol(dens.mat)), ann=TRUE, axes=FALSE,
          y=as.numeric(rownames(dens.mat)),
          z=t(dens.mat), col=colorRampPalette(hsv(h=27/360,s=1,v=seq(0,1,by=0.001)), bias=1.25)(100),
          ylab = "Read Length");
    abline(v=seq_len(ncol(dens.mat)+1)-0.5, lwd=3);
    abline(h=log10(c(1,2,5)) + rep(0:5, each=3),
           lty="dashed", col="#FFFFFF40");
    axis(1,at=seq_len(length(fileNames)),
         labels=colnames(dens.mat),
         lwd=0);
    axis(2,at=log10(c(1,2,5)) + rep(0:5, each=3), las=2,
         labels=valToSci(as.numeric(paste0(c(1,2,5),
             rep(substring("00000",first=0,last=0:5),each=3)))));
    dummy <- dev.off();
}

{
    png("MinION_Bases_Cumulative.png", pointsize=24,
        width=1600, height=960);
    par(mgp=c(4,1,0), mar=c(5.5,5.5,0.5,0.5));
    plot(NA, xlim=range(as.numeric(rownames(bpdens.mat))), ylim=c(0,1),
         type="l", xaxt="n", xlab = "Read Length",
         ylab = "");
    for(col in seq_len(ncol(bpdens.mat))){
        points(as.numeric(rownames(bpdens.mat)),1-cumsum(bpdens.mat[,col]) /
                   sum(bpdens.mat[,col]), type="l",
               col=hcl(h=col/ncol(bpdens.mat)*360, l=70, c=80), lwd=3);
    }
    legend("bottomleft", cex=0.71,
           legend=colnames(bpdens.mat), ncol=ceiling(ncol(bpdens.mat)/16),
           fill=hcl(h=(1:ncol(bpdens.mat))/ncol(bpdens.mat)*360, l=70, c=80),
           inset=0.05);
    mtext("Cumulative Base Proportion",2,3);
    axis(1,at=log10(c(1,2,5)) + rep(0:5, each=3), las=2,
         labels=valToSci(as.numeric(paste0(c(1,2,5),
             rep(substring("00000",first=0,last=0:5),each=3)))));
    dummy <- dev.off();
}

{
    png("MinION_Reads_Cumulative.png", pointsize=24,
        width=1600, height=960);
    par(mgp=c(4,1,0), mar=c(5.5,5.5,0.5,0.5));
    plot(NA, xlim=range(as.numeric(rownames(dens.mat))), ylim=c(0,1),
         type="l", xaxt="n", xlab = "Read Length",
         ylab = "");
    for(col in seq_len(ncol(dens.mat))){
        points(as.numeric(rownames(dens.mat)),1-cumsum(dens.mat[,col]) /
                   sum(dens.mat[,col]), type="l",
               col=hcl(h=col/ncol(dens.mat)*360, l=70, c=80), lwd=3);
    }
    legend("bottomleft", cex=0.71,
           legend=colnames(dens.mat), ncol=ceiling(ncol(bpdens.mat)/16),
           fill=hcl(h=(1:ncol(dens.mat))/ncol(dens.mat)*360, l=70, c=80),
           inset=0.05);
    mtext("Cumulative Read Proportion",2,3);
    axis(1,at=log10(c(1,2,5)) + rep(0:5, each=3), las=2,
         labels=valToSci(as.numeric(paste0(c(1,2,5),
             rep(substring("00000",first=0,last=0:5),each=3)))));
    dummy <- dev.off();
}

{
    png("MinION_Reads_Density.png", pointsize=24,
        width=1600, height=960);
    par(mgp=c(4,1,0));
    ymax <- max(dens.mat);
    plot(NA, xlim=range(as.numeric(rownames(dens.mat))), ylim=c(0,ymax),
         type="l", xaxt="n", xlab = "Read Length",
         ylab = "");
    for(col in seq_len(ncol(dens.mat))){
        points(as.numeric(rownames(dens.mat)),dens.mat[,col], type="l",
               col=hcl(h=col/ncol(dens.mat)*360, l=70, c=80), lwd=3);
    }
    legend("topright", legend=colnames(dens.mat),
           fill=hcl(h=(1:ncol(dens.mat))/ncol(dens.mat)*360, l=70, c=80),
           inset=0.05);
    mtext("Read Density",2,3);
    axis(1,at=log10(c(1,2,5)) + rep(0:5, each=3), las=2,
         labels=valToSci(as.numeric(paste0(c(1,2,5),
             rep(substring("00000",first=0,last=0:5),each=3)))));
    dummy <- dev.off();
}

{
    png("MinION_Bases_Density.png", pointsize=24,
        width=1600, height=960);
    ymax <- max(t(bpdens.mat) / colSums(bpdens.mat));
    par(mgp=c(4,1,0));
    plot(NA, xlim=range(as.numeric(rownames(bpdens.mat))), ylim=c(0,ymax),
         type="l", xaxt="n", xlab = "Read Length",
         ylab = "");
    for(col in seq_len(ncol(bpdens.mat))){
        points(as.numeric(rownames(bpdens.mat)),bpdens.mat[,col] /
                   sum(bpdens.mat[,col]), type="l", lwd=3,
               col=hcl(h=col/ncol(bpdens.mat)*360, l=70, c=80));
    }
    legend("topleft", cex=0.71, legend=colnames(bpdens.mat), ncol=ceiling(ncol(bpdens.mat)/16),
           fill=hcl(h=(1:ncol(bpdens.mat))/ncol(bpdens.mat)*360, l=70, c=80),
           inset=0.05);
    mtext("Base Density (arbitrary scale)",2,3);
    axis(1,at=log10(c(1,2,5)) + rep(0:5, each=3), las=2,
         labels=valToSci(as.numeric(paste0(c(1,2,5),
             rep(substring("00000",first=0,last=0:5),each=3)))));
    dummy <- dev.off();
}


{
    pdf("MinION_Reads_DigElec.pdf", width=12, height=8);
    par(mar=c(3,5,0.5,0.5), mgp=c(4,1,0));
    image(x=seq_len(ncol(bpdens.mat)), ann=TRUE, axes=FALSE,
          y=as.numeric(rownames(bpdens.mat)),
          z=t(bpdens.mat), col=colorRampPalette(hsv(h=27/360,s=0,v=seq(1,0,by=-0.001)), bias=1)(100),
          ylab = "Read Length");
    abline(v=seq_len(ncol(bpdens.mat)+1)-0.5);
    abline(h=log10(c(1,2,5)) + rep(0:5, each=3),
           lty="dashed", col="#00000020");
    axis(1,at=seq_len(length(fileNames)),labels=colnames(bpdens.mat),
         lwd=0);
    axis(2,at=log10(c(1,2,5)) + rep(0:5, each=3), las=2,
         labels=valToSci(as.numeric(paste0(c(1,2,5),
             rep(substring("00000",first=0,last=0:5),each=3)))));
    image(x=seq_len(ncol(bpdens.mat)), ann=TRUE, axes=FALSE,
          y=as.numeric(rownames(bpdens.mat)),
          z=t(bpdens.mat), col=colorRampPalette(hsv(h=27/360,s=1,v=seq(0,1,by=0.001)), bias=1.25)(100),
          ylab = "Read Length");
    abline(v=seq_len(ncol(bpdens.mat)+1)-0.5, col="#000000", lwd=5);
    abline(h=log10(c(1,2,5)) + rep(0:5, each=3),
           lty="dashed", col="#FFFFFF40");
    axis(1,at=seq_len(length(fileNames)),labels=colnames(bpdens.mat),
         lwd=0);
    axis(2,at=log10(c(1,2,5)) + rep(0:5, each=3), las=2,
         labels=valToSci(as.numeric(paste0(c(1,2,5),
             rep(substring("00000",first=0,last=0:5),each=3)))));
    ## Cumulative Reads
    par(mar=c(5.5,5,0.5,0.5), mgp=c(4,1,0));
    plot(NA, xlim=range(as.numeric(rownames(dens.mat))), ylim=c(0,1),
         type="l", xaxt="n", xlab = "Read Length",
         ylab = "");
    for(col in seq_len(ncol(dens.mat))){
        points(as.numeric(rownames(dens.mat)),1-cumsum(dens.mat[,col]) /
                   sum(dens.mat[,col]), type="l", lwd=3,
               col=hcl(h=col/ncol(dens.mat)*360, l=70, c=80));
    }
    legend("topright", legend=colnames(dens.mat),
           fill=hcl(h=(1:ncol(dens.mat))/ncol(dens.mat)*360, l=70, c=80),
           inset=0.05);
    mtext("Cumulative Read Proportion",2,3);
    axis(1,at=log10(c(1,2,5)) + rep(0:5, each=3), las=2,
         labels=valToSci(as.numeric(paste0(c(1,2,5),
             rep(substring("00000",first=0,last=0:5),each=3)))));
    ## Cumulative Bases
    plot(NA, xlim=range(as.numeric(rownames(bpdens.mat))), ylim=c(0,1),
         type="l", xaxt="n", xlab = "Read Length",
         ylab = "");
    for(col in seq_len(ncol(bpdens.mat))){
        points(as.numeric(rownames(bpdens.mat)),1-cumsum(bpdens.mat[,col]) /
                   sum(bpdens.mat[,col]), type="l", lwd=3,
               col=hcl(h=col/ncol(bpdens.mat)*360, l=70, c=80));
    }
    legend("topright", legend=colnames(bpdens.mat),
           fill=hcl(h=(1:ncol(bpdens.mat))/ncol(bpdens.mat)*360, l=70, c=80),
           inset=0.05);
    mtext("Cumulative Base Proportion",2,3);
    axis(1,at=log10(c(1,2,5)) + rep(0:5, each=3), las=2,
         labels=valToSci(as.numeric(paste0(c(1,2,5),
             rep(substring("00000",first=0,last=0:5),each=3)))));
    ## Base call density plot
    ymax <- max(t(bpdens.mat) / colSums(bpdens.mat));
    par(mgp=c(4,1,0), mar=c(5.5,5,1,1));
    plot(NA, xlim=range(as.numeric(rownames(bpdens.mat))), ylim=c(0,ymax),
         type="l", xaxt="n", xlab = "Read Length",
         ylab = "");
    for(col in seq_len(ncol(bpdens.mat))){
        points(as.numeric(rownames(bpdens.mat)),bpdens.mat[,col] /
                   sum(bpdens.mat[,col]), type="l", lwd=3,
               col=hcl(h=col/ncol(bpdens.mat)*360, l=70, c=80));
    }
    legend("topleft", cex=0.71, legend=colnames(bpdens.mat), ncol=ceiling(ncol(bpdens.mat)/16),
           fill=hcl(h=(1:ncol(bpdens.mat))/ncol(bpdens.mat)*360, l=70, c=80),
           inset=0.05);
    mtext("Base Density (arbitrary scale)",2,3);
    axis(1,at=log10(c(1,2,5)) + rep(0:5, each=3), las=2,
         labels=valToSci(as.numeric(paste0(c(1,2,5),
             rep(substring("00000",first=0,last=0:5),each=3)))));
    dummy <- dev.off();
}

