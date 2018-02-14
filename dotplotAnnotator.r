#!/usr/bin/Rscript

library(png);

fname <- commandArgs(TRUE);
img <- readPNG(fname);

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

fib.divs <- round(10^((0:4)/5) * 2) * 0.5; ## splits log decades into 5

seqName <- sub("^.*?_(.*).png","\\1",fname);
seqLen <- sub("_.*$","",fname);
seqStart <- 0;
seqEnd <- 0;
if(grepl("-", seqLen)){
    seqStart <- as.numeric(sub("-.*$","",seqLen));
    seqEnd <- as.numeric(sub("^.*?-","",seqLen));
    seqLen <- seqEnd - seqStart + 1;
} else {
    seqLen <- as.numeric(seqLen);
    seqStart <- 1;
    seqEnd <- seqLen;
}

print(c(seqName, seqLen, seqStart, seqEnd));

svg(sprintf("featurePlot_%s.svg", seqName),
    width=12.8, height=7.2, pointsize=12);
par(mar=c(5,6,3,1.5), cex.axis=1.5, cex.lab=1.5, cex.main=2);
plot(NA, main=sprintf("Feature profile (%s)", seqName),
     xlab="Sequence Location (kbp)",
     ylab="", log="y", xlim=c(seqStart,seqEnd)/1000,
     ylim=c(1,seqLen), yaxt="n");
rasterImage(img, xleft=seqStart/1000, xright=seqEnd/1000,
            ybottom=1, ytop=seqLen);
drMax <- ceiling(log10(seqLen));
axis(2, at= 10^(0:drMax), las=2, lwd=3, cex.axis=1.5, labels=valToSci(10^(0:drMax)));
axis(2, at= rep(1:9, each=drMax+1) * 10^(0:drMax), labels=FALSE);
abline(h=10^(0:drMax), col="#80808050", lwd = 3);
mtext("Feature distance (bp)", 2, line=4.5, cex=1.5);
legend(x = "bottom",
    fill=c("#9000a0","#8b0000","#00a090","#0000FF","#a09000","#00A000"),
    legend=c("Repeat (L)",  "Repeat (R)",
             "RevComp (L)", "RevComp (R)",
             "Reverse (L)", "Reverse (R)"),
    bg="#FFFFFF",
    horiz=TRUE, inset=0.01);
invisible(dev.off());
