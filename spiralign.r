#!/usr/bin/Rscript

setwd("/bioinf/MIMR-2017-Jul-01-GBIS/GLG/spiralign/");

## Nucleotides
library(msa);
input.seqs <- readDNAStringSet("notUSCO_EOG091H04CB.tran.fa");
names(input.seqs) <- sub(" .*$","",names(input.seqs));
msa.df <- data.frame(t(as.matrix(msa(input.seqs))), stringsAsFactors=FALSE);
colnames(msa.df) <- 1:3;

efg.cols <- c("G" = "gold", "C" = "blue", "A" = "darkgreen",
              "T" = "red", "-" = "grey20");

png("msa.png", width=1024, height=1024, pointsize=24);
par(mar=c(0.5,0.5,0.5,0.5), bg="black");
loops <- 5;
lstt <- 3;
lend <- loops+lstt;
## integrate(2*pi*r,r=lstt..x)
## => pi(x²-(lstt)²)
dTot <- pi*((lstt + loops)^2 - (lstt)^2); ## total "distance" travelled
## s = pi(x²-(lstt)²)
## => s/pi = x² - (lstt)²
## => x = sqrt((lstt)² + s/pi) 
msa.df$s <- seq(0,dTot, length.out=nrow(msa.df)); ## distance at each pos
msa.df$r <- sqrt(lstt^2 + msa.df$s/pi); ## path radius at each pos
msa.df$theta <- msa.df$r * 2*pi; ## traversed angle at each pos
msa.df$deg <- (msa.df$theta / (2*pi)) * 360;
msa.df$x <- msa.df$r * cos(msa.df$theta);
msa.df$y <- msa.df$r * sin(msa.df$theta);
plot(NA,xlim=c(-lend,lend), ylim=c(-lend,lend), ann=FALSE, axes=FALSE);
pcex <- 0.45;
for(p in seq(1,nrow(msa.df))){
    pr <- msa.df$r[p];
    pt <- msa.df$theta[p];
    text(x=-(pr-0.25)*cos(pt), y=(pr-0.25)*sin(pt), labels="▅",
         srt=-msa.df$deg[p],
         cex=pcex*0.9, col=efg.cols[msa.df[p,1]]);
    text(x=-pr*cos(pt), y=pr*sin(pt), labels="▅",
         srt=-msa.df$deg[p],
         cex=pcex, col=efg.cols[msa.df[p,2]]);
    text(x=-(pr+0.25)*cos(pt), y=(pr+0.25)*sin(pt), labels="▅",
         srt=-msa.df$deg[p],
         cex=pcex*1.05, col=efg.cols[msa.df[p,3]]);
}
text(0,0.5, expression(italic(Nippostrongylus)), col="white");
text(0,0, expression(italic(brasiliensis)), col="white");
text(0,-0.5, "Fructose-1,6-bisphosphatase", col="white", cex=0.75);
invisible(dev.off());

## Amino Acids
library(msa);
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
loops <- 2.75;
lstt <- 3;
lend <- loops+lstt;
## integrate(2*pi*r,r=lstt..x)
## => pi(x²-(lstt)²)
dTot <- pi*((lstt + loops)^2 - (lstt)^2); ## total "distance" travelled
## s = pi(x²-(lstt)²)
## => s/pi = x² - (lstt)²
## => x = sqrt((lstt)² + s/pi) 
msa.df$s <- seq(0,dTot, length.out=nrow(msa.df)); ## distance at each pos
msa.df$r <- sqrt(lstt^2 + msa.df$s/pi); ## path radius at each pos
msa.df$theta <- msa.df$r * 2*pi; ## traversed angle at each pos
msa.df$deg <- (msa.df$theta / (2*pi)) * 360;
msa.df$x <- msa.df$r * cos(msa.df$theta);
msa.df$y <- msa.df$r * sin(msa.df$theta);
plot(NA,xlim=c(-lend,lend), ylim=c(-lend,lend), ann=FALSE, axes=FALSE);
pcex <- 0.9;
for(p in seq(1,nrow(msa.df))){
    pr <- msa.df$r[p];
    pt <- msa.df$theta[p];
    text(x=-(pr-0.25)*cos(pt), y=(pr-0.25)*sin(pt), labels="▅",
         srt=-msa.df$deg[p],
         cex=pcex*0.9, col=rasmol.cols[msa.df[p,1]]);
    text(x=-pr*cos(pt), y=pr*sin(pt), labels="▅",
         srt=-msa.df$deg[p],
         cex=pcex, col=rasmol.cols[msa.df[p,2]]);
    text(x=-(pr+0.25)*cos(pt), y=(pr+0.25)*sin(pt), labels="▅",
         srt=-msa.df$deg[p],
         cex=pcex*1.05, col=rasmol.cols[msa.df[p,3]]);
}
text(0,0.5, expression(italic(Nippostrongylus)), col="white");
text(0,0, expression(italic(brasiliensis)), col="white");
text(0,-0.5, "Fructose-1,6-bisphosphatase", col="white", cex=0.75);
invisible(dev.off());
