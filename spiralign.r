#!/usr/bin/Rscript

library(msa);

## Nucleotides
input.seqs <- readDNAStringSet("notUSCO_EOG091H04CB.tran.fa");
names(input.seqs) <- sub(" .*$","",names(input.seqs));
msa.df <- data.frame(t(as.matrix(msa(input.seqs))), stringsAsFactors=FALSE);
colnames(msa.df) <- 1:3;

efg.cols <- c("G" = "gold", "C" = "blue", "A" = "darkgreen",
              "T" = "red", "-" = "grey90");

png("msa.png", width=1024, height=1024, pointsize=24);
par(mar=c(0.5,0.5,0.5,0.5));
loops <- 5;
lstt <- 3;
lend <- loops+lstt;
msa.df$r <- seq(0,loops,length.out=nrow(msa.df))+lstt;
msa.df$theta <- seq(0,loops,length.out=nrow(msa.df)) * 2*pi;
msa.df$deg <- seq(0,loops,length.out=nrow(msa.df)) * 360+90;
msa.df$x <- msa.df$r * cos(msa.df$theta);
msa.df$y <- msa.df$r * sin(msa.df$theta);
plot(NA,xlim=c(-lend,lend), ylim=c(-lend,lend), ann=FALSE, axes=FALSE);
for(p in seq(1,nrow(msa.df))){
    pr <- msa.df$r[p];
    pt <- msa.df$theta[p];
    text(x=(pr-0.25)*cos(pt), y=(pr-0.25)*sin(pt), labels="▅",
         srt=msa.df$deg[p],
         cex=0.5, col=efg.cols[msa.df[p,1]]);
    text(x=pr*cos(pt), y=pr*sin(pt), labels="▅",
         srt=msa.df$deg[p],
         cex=0.5, col=efg.cols[msa.df[p,2]]);
    text(x=(pr+0.25)*cos(pt), y=(pr+0.25)*sin(pt), labels="▅",
         srt=msa.df$deg[p],
         cex=0.5, col=efg.cols[msa.df[p,3]]);
}
text(0,0.5, expression(italic(Nippostrongylus)));
text(0,0, expression(italic(brasiliensis)));
text(0,-0.5, expression(italic(fbp)));
invisible(dev.off());

## Amino Acids
input.seqs <- readAAStringSet("notUSCO_EOG091H04CB.prot.fa");
names(input.seqs) <- sub(" .*$","",names(input.seqs));
msa.df <- data.frame(t(as.matrix(msa(input.seqs))), stringsAsFactors=FALSE);
colnames(msa.df) <- 1:3;

rasmol.cols <- c("D" = "#E60A0A", "E" = "#E60A0A",
                 "C" = "#E6E600", "M" = "#E6E600",
                 "K" = "#145AFF", "R" = "#145AFF",
                 "S" = "#FA9600", "T" = "#FA9600",
                 "F" = "#3232AA", "Y" = "#3232AA",
                 "N" = "#00DCDC", "Q" = "#00DCDC",
                 "G" = "#EBEBEB",
                 "L" = "#0F820F", "V" = "#0F820F", "I" = "#0F820F",
                 "A" = "#C8C8C8",
                 "W" = "#B45AB4",
                 "H" = "#8282D2",
                 "P" = "#DC9682",
                 "-" = "grey20");

png("msa_aa.png", width=1024, height=1024, pointsize=24);
par(mar=c(0.5,0.5,0.5,0.5), bg="black");
loops <- 5;
lstt <- 3;
lend <- loops+lstt;
msa.df$r <- seq(0,loops,length.out=nrow(msa.df))+lstt;
msa.df$theta <- seq(0,loops,length.out=nrow(msa.df)) * 2*pi;
msa.df$deg <- seq(0,loops,length.out=nrow(msa.df)) * 360+90;
msa.df$x <- msa.df$r * cos(msa.df$theta);
msa.df$y <- msa.df$r * sin(msa.df$theta);
plot(NA,xlim=c(-lend,lend), ylim=c(-lend,lend), ann=FALSE, axes=FALSE);
pcex <- 0.75;
for(p in seq(1,nrow(msa.df))){
    pr <- msa.df$r[p];
    pt <- msa.df$theta[p];
    text(x=-(pr-0.25)*cos(pt), y=(pr-0.25)*sin(pt), labels="▅",
         srt=-msa.df$deg[p],
         cex=pcex, col=rasmol.cols[msa.df[p,1]]);
    text(x=-pr*cos(pt), y=pr*sin(pt), labels="▅",
         srt=-msa.df$deg[p],
         cex=pcex, col=rasmol.cols[msa.df[p,2]]);
    text(x=-(pr+0.25)*cos(pt), y=(pr+0.25)*sin(pt), labels="▅",
         srt=-msa.df$deg[p],
         cex=pcex, col=rasmol.cols[msa.df[p,3]]);
}
text(0,0.5, expression(italic(Nippostrongylus)), col="white");
text(0,0, expression(italic(brasiliensis)), col="white");
text(0,-0.5, "Fructose-1,6-bisphosphatase", col="white", cex=0.75);
invisible(dev.off());
