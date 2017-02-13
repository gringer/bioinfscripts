#!/usr/bin/Rscript

fileNames <- commandArgs(TRUE);

if(length(fileNames) == 0){
    fileNames <- list.files(pattern="lengths_.*\\.txt(\\.gz)?");
}

dens.mat <- sapply(fileNames, function(x){
    cat(x,"...");
    data <- scan(x, comment.char=" ", quiet=TRUE);
    cat(" done\n");
    res <- density(log10(data), from=2, to=6, bw=0.1);
    res.out <- res$y;
    names(res.out) <- round(res$x,3);
    res.out;
});
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
    axis(1,at=seq_len(length(fileNames)),labels=colnames(bpdens.mat),
         lwd=0);
    axis(2,at=log10(c(1,2,5)) + rep(0:5, each=3), las=1,
         labels=paste0(c(1,2,5),
                       rep(substring("00000",first=0,last=0:5),each=3)));
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
    axis(1,at=seq_len(length(fileNames)),labels=colnames(dens.mat),
         lwd=0);
    axis(2,at=log10(c(1,2,5)) + rep(0:5, each=3), las=1,
         labels=paste0(c(1,2,5),
                       rep(substring("00000",first=0,last=0:5),each=3)));
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
    axis(1,at=seq_len(length(fileNames)),labels=colnames(bpdens.mat),
         lwd=0);
    axis(2,at=log10(c(1,2,5)) + rep(0:5, each=3), las=1,
         labels=paste0(c(1,2,5),
                       rep(substring("00000",first=0,last=0:5),each=3)));
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
    axis(1,at=seq_len(length(fileNames)),labels=colnames(dens.mat),
         lwd=0);
    axis(2,at=log10(c(1,2,5)) + rep(0:5, each=3), las=1,
         labels=paste0(c(1,2,5),
                       rep(substring("00000",first=0,last=0:5),each=3)));
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
         labels=paste0(c(1,2,5),
             rep(substring("00000",first=0,last=0:5),each=3)));
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
         labels=paste0(c(1,2,5),
             rep(substring("00000",first=0,last=0:5),each=3)));
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
         labels=paste0(c(1,2,5),
             rep(substring("00000",first=0,last=0:5),each=3)));
    dummy <- dev.off();
}

{
    png("MinION_Bases_Density.png", pointsize=24,
        width=1600, height=960);
    ymax <- max(bpdens.mat / colSums(bpdens.mat));
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
         labels=paste0(c(1,2,5),
             rep(substring("00000",first=0,last=0:5),each=3)));
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
    axis(2,at=log10(c(1,2,5)) + rep(0:5, each=3), las=1,
         labels=paste0(c(1,2,5),
                       rep(substring("00000",first=0,last=0:5),each=3)));
    image(x=seq_len(ncol(bpdens.mat)), ann=TRUE, axes=FALSE,
          y=as.numeric(rownames(bpdens.mat)),
          z=t(bpdens.mat), col=colorRampPalette(hsv(h=27/360,s=1,v=seq(0,1,by=0.001)), bias=1.25)(100),
          ylab = "Read Length");
    abline(v=seq_len(ncol(bpdens.mat)+1)-0.5, col="#000000", lwd=5);
    abline(h=log10(c(1,2,5)) + rep(0:5, each=3),
           lty="dashed", col="#FFFFFF40");
    axis(1,at=seq_len(length(fileNames)),labels=colnames(bpdens.mat),
         lwd=0);
    axis(2,at=log10(c(1,2,5)) + rep(0:5, each=3), las=1,
         labels=paste0(c(1,2,5),
                       rep(substring("00000",first=0,last=0:5),each=3)));
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
         labels=paste0(c(1,2,5),
             rep(substring("00000",first=0,last=0:5),each=3)));
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
         labels=paste0(c(1,2,5),
             rep(substring("00000",first=0,last=0:5),each=3)));
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
         labels=paste0(c(1,2,5),
             rep(substring("00000",first=0,last=0:5),each=3)));
    dummy <- dev.off();
}

