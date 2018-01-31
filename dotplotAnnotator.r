#!/usr/bin/Rscript

valToSci <- function(val, unit = ""){
    sci.prefixes <- c("", "k", "M", "G", "T", "P", "E", "Z", "Y");
    units <- rep(paste(sci.prefixes,unit,sep=""), each=3);
    logRegion <- floor(log10(val))+1;
    conv.units <- units[logRegion];
    conv.div <- 10^rep(0:(length(sci.prefixes)-1) * 3, each = 3)[logRegion];
    conv.val <- val / conv.div;
    conv.val[val == 0] <- 0;
    conv.units[val == 0] <- unit;
    return(ifelse(conv.units == "", conv.val,
                  sprintf("%s %s",conv.val,conv.units)));
}

library(png);
img.alt <- readPNG("alt_tig7744.png");
seqLen <- 119722;

png("tig7744_thinkable.png", width=1280, height=720, pointsize=16);
par(mar=c(4.5,5.5,0.5,0.5));
plot(NA, xlim=c(1,seqLen/1000), ylim=c(1,seqLen),
     log="y", ann=FALSE, axes=FALSE,
     frame.plot=TRUE);
rasterImage(img.alt, xleft=1, xright=seqLen/1000, ybottom=1, ytop=seqLen);
legend("bottom", horiz=TRUE, inset=0.025, cex=1.2,
       legend=c("Repeat(L)","Repeat(R)",
                "Reverse(L)","Reverse(R)","RevComp(L)","RevComp(R)"),
       fill=c("#9000a0","#8b0000","#00a090","#0000ff","#a09000","#00a000"));
mtext("tig7744 location (kbp)", side=1, line=3, cex=1.5);
mtext("Feature distance (bp)", side=2, line=3.5, cex=1.5);
axis(1, cex.axis=1.5, lwd=3);
drMax <- ceiling(log10(seqLen));
axis(2, at= 10^(0:drMax), lwd=3, cex.axis=1.5, labels=valToSci(10^(0:drMax)),
     las=2);
axis(2, at= rep(1:9, each=drMax+1) * 10^(0:drMax), labels=FALSE);
abline(h=10^(1:drMax), col="#80808040",
       lwd=3);
invisible(dev.off());
