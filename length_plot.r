#!/usr/bin/Rscript

## Read length distribution plotter; expects input files named 'lengths_<id>.txt[.gz]'
## The first [space-separated] column will be processed as sequence lengths

scriptArgs <- commandArgs(TRUE);

fileNames <-
    c("/mnt/gg_nanopore/gringer/reads_LC_LN18_2017-Oct-31/called_LC_LN18_20h/workspace/lengths_1_C_BC01.txt",
      "/mnt/gg_nanopore/gringer/reads_LC_LN18_2017-Oct-31/called_LC_LN18_20h/workspace/lengths_2_Gy_BC09.txt",
      "/mnt/gg_nanopore/gringer/reads_LC_LN18_2017-Oct-31/called_LC_LN18_20h/workspace/lengths_3_TMZ_BC03.txt",
      "/mnt/gg_nanopore/gringer/reads_LC_LN18_2017-Oct-31/called_LC_LN18_20h/workspace/lengths_4_Dox_BC05.txt",
      "/mnt/gg_nanopore/gringer/reads_LC_LN18_2017-Oct-31/called_LC_LN18_20h/workspace/lengths_5_H2O2_BC06.txt",
      "/mnt/gg_nanopore/gringer/reads_LC_LN18_2017-Oct-31/called_LC_LN18_20h/workspace/lengths_6_UV_BC07.txt",
      "/mnt/gg_nanopore/gringer/reads_LC_LN18_2017-Oct-31/called_LC_LN18_20h/workspace/lengths_8_+ve_BC08.txt",
      "/mnt/gg_nanopore/gringer/reads_LC_LN18_2017-Oct-31/called_LC_LN18_20h/workspace/lengths_unclassified.txt");










plotVals <- TRUE;
plotHoriz <- TRUE;
plotCombined <- FALSE;

fileNames <- scriptArgs[file.exists(scriptArgs)];
scriptArgs <- scriptArgs[!file.exists(scriptArgs)];

if("-n" %in% scriptArgs){
    plotVals <- FALSE;
}

if("-h" %in% scriptArgs){
    plotHoriz <- FALSE;
}

if("-c" %in% scriptArgs){
    plotCombined <- TRUE;
}

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

sequence.hist <- function(lengths, horiz = TRUE, barValues=TRUE, invert = TRUE, main = "", ...){
    fib.divs <- round(10^((0:4)/5) * 2) * 0.5; ## splits log decades into 5
    histBreaks <- round(rep(10^(0:16),each=5) * fib.divs);
    lengthRange <- range(lengths);
    ## filter on actual data range
    histBreaks <- unique(histBreaks[(which.min(histBreaks < lengthRange[1])-1):
                                    which.max(histBreaks > lengthRange[2])]);
    seqd.bases <- seqd.na.bases <- tapply(lengths,cut(lengths, breaks=histBreaks), sum);
    seqd.counts <- seqd.na.counts <- tapply(lengths,cut(lengths, breaks=histBreaks), length);
    seqd.bases[is.na(seqd.bases)] <- 0;
    seqd.counts[is.na(seqd.counts)] <- 0;
    xBreaks <- round(rep(10^(0:16),each=5) * fib.divs);
    axRange <- range(seqd.bases);
    xBreaks <- xBreaks[(which.min(xBreaks < axRange[1])):
                       (which.max(xBreaks > axRange[2])-1)];
    xBreaksMajor <- xBreaks[log10(xBreaks) - floor(log10(xBreaks)) < 0.001];
    xBreaksMinor <- rep(xBreaksMajor,each=9) * 1:9;
    xBreaksMinor <- xBreaksMinor[(which.min(xBreaksMinor < axRange[1])):
                                 (which.max(xBreaksMinor > axRange[2])-1)];
    barPos <- barplot(if(invert){seqd.na.bases} else {rev(seqd.na.bases)}, main = "",
                      log = ifelse(horiz, "x", "y"), axes = FALSE, col = "steelblue",
                      horiz = horiz, names.arg = rep("",length(seqd.bases)), ann=FALSE,
                      ...);
    barGap <- diff(barPos)[1];
    barOffset <- barPos[1] - barGap/2;
    axis(ifelse(horiz,2,1), at = if(!invert){rev(seq(barOffset,by=barGap,
                                    length.out = length(histBreaks)))}
         else {seq(barOffset,by=barGap,length.out = length(histBreaks))},
         labels = valToSci(histBreaks,"b"), las = 2, cex.axis=ifelse(horiz,1,1.41));
    mtext(main, side=3, line=ifelse(horiz,0,1.5), cex=2);
    mtext("Fragment size", side=ifelse(horiz,2,1), line=4);
    mtext("Aggregate sequence length (number of sequences)", side=ifelse(horiz,1,2), line=ifelse(horiz,4,5));
    axis(ifelse(horiz,1,2), at = xBreaksMajor, labels = valToSci(xBreaksMajor, "b"), lwd=3, las=1, cex.axis=1.41);
    axis(ifelse(horiz,1,2), at = c(head(xBreaksMinor,1),tail(xBreaksMinor,1)), labels = FALSE, lwd=3);
    axis(ifelse(horiz,1,2), at = xBreaksMinor, labels = FALSE);
    if(barValues && horiz){
        text.poss <- ((log10(seqd.bases) < mean(par("usr")[1:2]))+1)*2;
        text.poss[is.na(text.poss)] <- 4;
        text.col <- c("white","black")[((log10(seqd.bases) <
                                         mean(par("usr")[1:2]))+1)];
        text(seqd.bases,if(!invert){rev(barPos)} else {barPos},
             paste(valToSci(signif(seqd.bases,4)),
                   " (", seqd.counts, ")", sep = ""),
             pos=text.poss, col=text.col, cex = 0.65);
    }
}

plain.hist <- function(lengths, horiz=TRUE, barValues=TRUE, invert = TRUE, main = "", ...){
    fib.divs <- round(10^((0:4)/5) * 2) * 0.5; ## splits log decades into 5
    histBreaks <- round(rep(10^(0:16),each=5) * fib.divs);
    lengthRange <- range(lengths);
    ## filter on actual data range
    histBreaks <- unique(histBreaks[(which.min(histBreaks < lengthRange[1])-1):
                                    which.max(histBreaks > lengthRange[2])]);
    seqd.bases <- seqd.na.bases <- tapply(lengths,cut(lengths, breaks=histBreaks), sum);
    seqd.counts <- seqd.na.counts <- tapply(lengths,cut(lengths, breaks=histBreaks), length);
    seqd.bases[is.na(seqd.bases)] <- 0;
    seqd.counts[is.na(seqd.counts)] <- 0;
    xBreaks <- round(rep(10^(0:16),each=5) * fib.divs);
    axRange <- range(seqd.counts);
    xBreaks <- xBreaks[(which.min(xBreaks < axRange[1])):
                       (which.max(xBreaks > axRange[2])-1)];
    xBreaksMajor <- xBreaks[log10(xBreaks) - floor(log10(xBreaks)) < 0.001];
    xBreaksMinor <- rep(xBreaksMajor,each=9) * 1:9;
    xBreaksMinor <- xBreaksMinor[(which.min(xBreaksMinor < min(xBreaksMajor))):
                                 (which.max(xBreaksMinor > max(xBreaksMajor))-1)];
    barPos <- barplot(log10(if(!invert){rev(seqd.na.counts)} else {seqd.na.counts})+0.025,
                      las = 1, axes = FALSE, ann=FALSE, col = "steelblue", main = "",
                      horiz = horiz, names.arg = rep("",length(seqd.counts)),
                      xlim=if(!horiz){NULL} else {c(log10(max(seqd.na.counts, na.rm=TRUE))+0.025,
                                                   log10(min(seqd.na.counts, na.rm=TRUE))-0.025)},
                      ylim=if(horiz){NULL} else {c(log10(max(seqd.na.counts, na.rm=TRUE))+0.025,
                                                   log10(min(seqd.na.counts, na.rm=TRUE))-0.025)},
                      ...);
    mtext(main, side=ifelse(horiz,3,1), line=ifelse(horiz,0,1.5), cex=2);
    mtext("Number of sequences (Aggregate length)", side=ifelse(horiz,1,2), line=ifelse(horiz,4,5));
    barGap <- diff(barPos)[1];
    barOffset <- barPos[1] - barGap/2;
    axis(ifelse(horiz,4,3), at = if(!invert){rev(seq(barOffset,by=barGap,
                                                    length.out = length(histBreaks)))}
                                 else {seq(barOffset,by=barGap,length.out = length(histBreaks))},
         labels = valToSci(histBreaks,"b"), las = 2, pos=0, cex.axis=ifelse(horiz,1,1.41));
    mtext("Fragment size", side=ifelse(horiz,4,1), line=5);
    axis(ifelse(horiz,1,2), at = log10(xBreaksMajor)+0.025, labels = valToSci(xBreaksMajor),
         lwd=3, cex.axis=1.41, las=1);
    axis(ifelse(horiz,1,2), at = c(head(xBreaksMinor,1),tail(xBreaksMinor,1)), labels = FALSE, lwd=3);
    axis(ifelse(horiz,1,2), at = log10(xBreaksMinor)+0.025, labels = FALSE);
    if(barValues && horiz){
        text.poss <- ((log10(seqd.counts) > mean(par("usr")[1:2]))+1)*2;
        text.poss[is.na(text.poss)] <- 2;
        text.col <- c("white","black")[((log10(seqd.counts) <
                                         mean(par("usr")[1:2]))+1)];
        text(log10(seqd.counts)+0.025,if(!invert){rev(barPos)} else {barPos},
             paste(seqd.counts,
                   " (", valToSci(signif(seqd.bases,4),"b"), ")", sep = ""),
             pos=text.poss, col=text.col, cex = 0.65);
    }
}

## generate fibonacci bins for logarithmic data storage
sampled.lengths <- NULL;
fib.divs <- round(10^(seq(0,1-1/5,length.out=5)) * 2) * 0.5; ## splits log decades into 5
histBreaks <- unique(c(0,round(rep(10^(0:16),each=5) * fib.divs)));
histCentres <- 10^((log10(head(histBreaks,-1)) + log10(tail(histBreaks,-1)))/2);
## fetch data and store in a table, binned as above
countData <- data.frame(row.names=tail(histBreaks, -1));
baseData <- data.frame(row.names=tail(histBreaks, -1));
dens.mat <- sapply(fileNames, function(x){
    cat(x,"...");
    file.data <- scan(x, comment.char=" ", quiet=TRUE);
    file.data <- file.data[file.data>0];
    cat(" done\n");
    cat("Number of sequences:",length(file.data),"\n");
    cat("Length quantiles:\n");
    print(quantile(file.data, probs=seq(0,1,by=0.1)));
    lengthRange <- range(file.data);
    seqd.bases <- tapply(file.data, cut(file.data, breaks=histBreaks), sum);
    seqd.counts <- tapply(file.data, cut(file.data, breaks=histBreaks), length);
    seqd.bases[is.na(seqd.bases)] <- 0;
    seqd.counts[is.na(seqd.counts)] <- 0;
    names(seqd.bases) <- valToSci(tail(histBreaks,-1));
    names(seqd.counts) <- valToSci(tail(histBreaks,-1));
    seqd.label <- sub(".txt(.gz)?$","",
                      sub("^lengths_(called_)?","",basename(x)));
    baseData[[seqd.label]] <<- seqd.bases;
    countData[[seqd.label]] <<- seqd.counts;
});

baseRange <- c(min(which(rowSums(countData) > 1)),
               max(which(rowSums(countData) > 1)));

library(RColorBrewer);

pdf("Sequence_curves.pdf", width=16, height=8);
#### Create density plots ####
## read counts
par(mar=c(4.5, 6.5, 1, 1), las=1, mgp=c(5,1,0), cex.lab=1.5, bg="white");
plot(NA, type="l", log="x", xaxt="n", yaxt="n",
     xlim=c(histCentres[baseRange[1]], histCentres[baseRange[2]]),
     ylim=c(0,max(countData)*1.1),
     ylab="Sequenced Reads", xlab = "");
drMax <- max(log10(histCentres[baseRange[2]]));
axis(1, at= (10^(0:(drMax+1))), lwd=3, cex.axis=1.5,
     labels=valToSci(10^(0:(drMax+1))));
axis(1, at= (rep(1:9, each=drMax+1) * 10^(0:drMax)), labels=FALSE);
axis(2, at= pretty(unlist(countData)),
     labels=valToSci(pretty(unlist(countData))), cex.axis = 1.5 );
mtext("Read Length", side=1, line=3, cex=1.5);
smoothed.count.mat <- sapply(1:ncol(countData), function(pcol){
    smoothed.data <- spline(x=log10(histCentres[baseRange[1]:baseRange[2]]),
                            y=countData[baseRange[1]:baseRange[2],pcol],
                            n=10*length(histBreaks));
    smoothed.data$x <- 10^smoothed.data$x;
    points(smoothed.data, type="l", lwd=3,
           col=brewer.pal(ncol(countData), "Set2")[pcol]);
    res <- smoothed.data$y;
    names(res) <- signif(smoothed.data$x,4);
    res;
});
legend("topright", fill=brewer.pal(ncol(baseData), "Set2")[1:ncol(baseData)],
       legend=colnames(baseData), inset=c(0.025,0.05), cex=1.5,
       bg="#FFFFFFD0");
#### Create density plots ####
## sequenced lengths
par(mar=c(4.5, 6.5, 1, 1), las=1, mgp=c(5,1,0), cex.lab=1.5);
plot(NA, type="l", log="x", xaxt="n", yaxt="n",
     xlim=c(histCentres[baseRange[1]], histCentres[baseRange[2]]),
     ylim=c(0,max(baseData)*1.1),
     ylab="Sequenced Bases", xlab = "");
drMax <- max(log10(histCentres[baseRange[2]]));
axis(1, at= (10^(0:(drMax+1))), lwd=3, cex.axis=1.5,
     labels=valToSci(10^(0:(drMax+1))));
axis(1, at= (rep(1:9, each=drMax+1) * 10^(0:drMax)), labels=FALSE);
axis(2, at= pretty(unlist(baseData)),
     labels=valToSci(pretty(unlist(baseData))), cex.axis = 1.5 );
mtext("Read Length", side=1, line=3, cex=1.5);
smoothed.base.mat <- sapply(1:ncol(baseData), function(pcol){
    smoothed.data <- spline(x=log10(histCentres[baseRange[1]:baseRange[2]]),
                            y=baseData[baseRange[1]:baseRange[2],pcol],
                            n=10*length(histBreaks));
    smoothed.data$x <- 10^smoothed.data$x;
    points(smoothed.data, type="l", lwd=3,
           col=brewer.pal(ncol(baseData), "Set2")[pcol]);
    res <- smoothed.data$y;
    names(res) <- signif(smoothed.data$x,4);
    res;
});
legend("topright", fill=brewer.pal(ncol(baseData), "Set2")[1:ncol(baseData)],
       legend=colnames(baseData), inset=c(0.025,0.05), cex=1.5,
       bg="#FFFFFFD0");
## sequenced base proportion
par(mar=c(4.5, 6.5, 1, 1), las=1, mgp=c(5,1,0), cex.lab=1.5);
plot(NA, type="l", log="x", xaxt="n", yaxt="n",
     xlim=c(histCentres[baseRange[1]], histCentres[baseRange[2]]),
     ylim=c(0,100),
     ylab="Cumulative Sequenced Bases (%)", xlab = "");
drMax <- max(log10(histCentres[baseRange[2]]));
axis(1, at= (10^(0:(drMax+1))), lwd=3, cex.axis=1.5,
     labels=valToSci(10^(0:(drMax+1))));
axis(1, at= (rep(1:9, each=drMax+1) * 10^(0:drMax)), labels=FALSE);
axis(2, at= (0:10 * 10), cex.axis = 1.5 );
mtext("Read Length", side=1, line=3, cex=1.5);
for(pcol in 1:ncol(baseData)){
    smoothed.data <-
        spline(x=log10(histCentres[baseRange[1]:baseRange[2]]),
               y=100 - cumsum(baseData[baseRange[1]:baseRange[2],pcol]) /
                   sum(baseData[baseRange[1]:baseRange[2],pcol]) * 100,
                            n=10*length(histBreaks));
    smoothed.data$x <- 10^smoothed.data$x;
    points(smoothed.data, type="l", lwd=3,
           col=brewer.pal(ncol(baseData), "Set2")[pcol]);
    points(c(rep(smoothed.data$x[order(abs((smoothed.data$y) - 90))[1]],2),
             rep(smoothed.data$x[order(abs((smoothed.data$y) - 50))[1]],2),
             rep(smoothed.data$x[order(abs((smoothed.data$y) - 10))[1]],2)),
           c(90,0,50,0,10,0), col=brewer.pal(ncol(baseData), "Set2")[pcol]);
    segments(x0=c(rep(smoothed.data$x[order(abs((smoothed.data$y) - 90))[1]]),
                  rep(smoothed.data$x[order(abs((smoothed.data$y) - 50))[1]]),
                  rep(smoothed.data$x[order(abs((smoothed.data$y) - 10))[1]])),
             y0=c(90,50,10), y1=c(0,0,0), lty=c("dotted","dashed"),
             col=brewer.pal(ncol(baseData), "Set2")[pcol]);
}
abline(h=c(10,50,90), lty="dashed", col="#80808080", lwd=2);
legend("topright", fill=brewer.pal(ncol(baseData), "Set2")[1:ncol(baseData)],
       legend=colnames(baseData), inset=c(0.025,0.05), cex=1.5,
       bg="#FFFFFFD0");
invisible(dev.off());


colnames(smoothed.base.mat) <- colnames(baseData);
smoothed.baseNorm.mat <- apply(smoothed.base.mat,2,function(x){x/sum(x)});

quit("no");

## Create density matrix and histogram plots
pdf("MinION_Reads_SequenceHist.pdf", paper="a4r",
    width=11, height=8);
par(mar=c(5.5,6.5,2.5,4.5));
dens.mat <- sapply(fileNames, function(x){
    cat(x,"...");
    data <- scan(x, comment.char=" ", quiet=TRUE);
    data <- data[data>0];
    cat(" done\n");
    cat("Number of sequences:",length(data),"\n");
    cat("Length quantiles:\n");
    print(quantile(data, probs=seq(0,1,by=0.1)));

    fib.divs <- round(10^(seq(0,1-1/5,length.out=5)) * 2) * 0.5; ## splits log decades into 5
    histBreaks <- unique(c(0,round(rep(10^(0:15),each=5) * fib.divs)));
    seqd.bases <- seqd.na.bases <- tapply(lengths,cut(lengths, breaks=histBreaks), sum);
    seqd.counts <- seqd.na.counts <- tapply(lengths,cut(lengths, breaks=histBreaks), length);
    seqd.bases[is.na(seqd.bases)] <- 0;
    names(seqd.bases) <- tail(histBreaks,-1);
    seqd.counts[is.na(seqd.counts)] <- 0;
    names(seqd.counts) <- tail(histBreaks,-1);

    res.out <- density(log10(data), bw=0.1);
    res.out <- res$y;
    names(res.out) <- round(res$x,3);
    subName <- sub("lengths_(.*)\\.txt(\\.gz)?","\\1", x);
    sampled.lengths <<- cbind(sampled.lengths, sample(data, 10000, replace=TRUE));
    str(sampled.lengths);
    if(plotCombined){
        png(sprintf("MinION_Reads_SequenceHist_%s.png",
                    subName), pointsize=24,
            width=2274, height=1046);
        if(plotHoriz){
            par(mfrow=c(1,2));
        } else {
            par(mfrow=c(2,1));
        }
    } else {
        png(sprintf("MinION_Reads_SequenceHist_%s.png",
                    subName), pointsize=24,
            width=1280, height=960);
    }
    par(mar=c(5.5,6.5,3.5,4.5));
    sequence.hist(data, horiz=plotHoriz, barValues=plotVals,
                  main=sprintf("Read Length Distribution Plot (%s)",subName));
    if(!plotCombined){
        dummy <- dev.off();
        png(sprintf("MinION_Reads_PlainHist_%s.png",
                    subName), pointsize=24,
            width=1280, height=960);
    }
    if(plotHoriz){
        par(mar=c(5.5,6.5,2.5,4.5));
    } else {
        par(mar=c(5,6.5,4.5,4.5));
    }
    plain.hist(data, horiz=plotHoriz, barValues=plotVals,
               main=sprintf("Read Count Distribution Plot (%s)",subName));
    dummy <- dev.off();
    par(mar=c(5.5,6.5,3.5,4.5));
    sequence.hist(data, horiz=plotHoriz, barValues=plotVals,
                  main=sprintf("Read Length Distribution Plot (%s)",subName));
    if(plotHoriz){
        par(mar=c(5.5,6.5,2.5,4.5));
    } else {
        par(mar=c(5,6.5,4.5,4.5));
    }
    plain.hist(data, horiz=plotHoriz, barValues=plotVals,
               main=sprintf("Read Count Distribution Plot (%s)",subName));
    res.out;
});
dummy <- dev.off();
colnames(dens.mat) <- sub("lengths_(.*)\\.txt(\\.gz)?","\\1",colnames(dens.mat));


quit(save="n");

bpdens.mat <- dens.mat * 10^as.numeric(rownames(dens.mat));

{
    png("DigElec_Bases_white.png", pointsize=24,
        width=80*ncol(smoothed.baseNorm.mat), height=960);
    par(mar=c(6,5,0.5,0.5), mgp=c(4,1,0));
    image(x=seq_len(ncol(smoothed.baseNorm.mat)), ann=TRUE, axes=FALSE,
          y=log10(as.numeric(rownames(smoothed.baseNorm.mat))),
          z=t(smoothed.baseNorm.mat), col=colorRampPalette(hsv(h=27/360,s=0,v=seq(1,0,by=-0.001)), bias=1)(100),
          ylab = "Read Length", xlab="");
    abline(v=seq_len(ncol(smoothed.baseNorm.mat)+1)-0.5, lwd=5, col="white");
    abline(h=log10(c(1,2,5)) + rep(0:5, each=3),
           lty="dashed", col="#00000020");
    axis(1,at=seq_len(length(fileNames)), labels=colnames(smoothed.baseNorm.mat),
         lwd=0, las=2, cex.axis=0.75);
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
                   sum(bpdens.mat[,col]), type="l", lwd=5,
               col=hcl(h=col/ncol(bpdens.mat)*360, l=70, c=80));
    }
    legend("topleft", cex=1, legend=colnames(bpdens.mat), ncol=ceiling(ncol(bpdens.mat)/16),
           fill=hcl(h=(1:ncol(bpdens.mat))/ncol(bpdens.mat)*360, l=70, c=80),
           inset=0.05);
    mtext("Base Density (arbitrary scale)",2,3);
    axis(1,at=log10(c(1,2,5)) + rep(0:5, each=3), las=2,
         labels=valToSci(as.numeric(paste0(c(1,2,5),
             rep(substring("00000",first=0,last=0:5),each=3)))));
    dummy <- dev.off();
}

{
    qqpfun <- function(A, B){
        points(qqplot(A, B, plot.it=FALSE));
    }
    if(ncol(dens.mat)>1){
        png("MinION_QQPlot.png", pointsize=24,
            width=1600, height=1600);
        pairs(sampled.lengths, labels=colnames(dens.mat), panel=qqpfun);
        invisible(dev.off());
    }
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
                   sum(bpdens.mat[,col]), type="l", lwd=4,
               col=hcl(h=col/ncol(bpdens.mat)*360, l=70, c=80));
    }
    legend("topleft", cex=1.5, legend=colnames(bpdens.mat), ncol=ceiling(ncol(bpdens.mat)/16),
           fill=hcl(h=(1:ncol(bpdens.mat))/ncol(bpdens.mat)*360, l=70, c=80),
           inset=0.05);
    mtext("Base Density (arbitrary scale)",2,3);
    axis(1,at=log10(c(1,2,5)) + rep(0:5, each=3), las=2,
         labels=valToSci(as.numeric(paste0(c(1,2,5),
             rep(substring("00000",first=0,last=0:5),each=3)))));
    dummy <- dev.off();
}

