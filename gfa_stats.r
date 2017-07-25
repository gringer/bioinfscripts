#!/usr/bin/Rscript

library(igraph);
library(Biostrings);
library(RColorBrewer);

ghost.data <- list();

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

sequence.hist <- function(lengths, lengthRange=NULL){
    fib.divs <- round(10^((0:4)/5) * 2) * 0.5; ## splits log decades into 5
    histBreaks <- round(rep(10^(0:16),each=5) * fib.divs);
    if(is.null(lengthRange)){
        lengthRange <- range(lengths);
    }
    ## filter on actual data range
    histBreaks <- histBreaks[(which.min(histBreaks < lengthRange[1])-1):
                             which.max(histBreaks > lengthRange[2])];
    seqd.bases <- seqd.na.bases <-
        tapply(lengths,cut(lengths, breaks=histBreaks), sum);
    seqd.counts <- seqd.na.counts <-
        tapply(lengths,cut(lengths, breaks=histBreaks), length);
    seqd.bases[is.na(seqd.bases)] <- 0;
    seqd.counts[is.na(seqd.counts)] <- 0;
    names(seqd.bases) <- paste0(head(valToSci(histBreaks),-1));
    names(seqd.bases)[length(seqd.bases)] <-
        paste0(tail(names(seqd.bases),1),";",
               tail(valToSci(histBreaks),1));
    seqd.bases;
}

wtsi.lengths <- read.table("~/db/fasta/nippo/lengths_Nb_WTSI.txt",
                           col.names=c("length","contig"));
ghost.data[["WTSI"]] <- sequence.hist(wtsi.lengths$length,
                                      lengthRange=myXRange);

setwd("/mnt/gg_nanopore/gringer/ONT_Jan17/GFA_stats");

#gfaName <- "Nb_ONTCFED_65bpTrim_t1.contigs.gfa";
gfaName <- "NbL5_ONTA.contigs.gfa";
#gfaName <- "Nb_ONTA_65bpTrim_t3.contigs.gfa";
#gfaName <- "NbL5_ONTDECAF_t1.contigs.flipped.gfa";
data.lines <- readLines(gfaName);

data.lengths <- sapply(strsplit(grep("^S", data.lines, value=TRUE),"\t"),
                          function(x){
                              val <- as.numeric(substring(x[4],6));
                              names(val) <- x[2];
                              val;
                          });
data.linklist <- strsplit(grep("^L", data.lines, value=TRUE),"\t");
data.link.df <- data.frame(t(sapply(data.linklist,c)), stringsAsFactors=FALSE);
colnames(data.link.df) <- c("type","from","fromDir",
                            "to","toDir","cigar");
data.link.df <- data.link.df[,1:5];
invertLines <- (data.link.df$fromDir == "-") & (data.link.df$toDir == "-");
data.link.df[invertLines,c("from","to")] <- data.link.df[invertLines,c("to","from")];
data.link.df[invertLines,"fromDir"] <- "+";
data.link.df[invertLines,"toDir"] <- "+";
data.link.df <- unique(data.link.df);

dup.fromDir <- table(paste0(data.link.df$from,data.link.df$fromDir));
dup.fromDir <- dup.fromDir[dup.fromDir > 1];
dup.fromDir <- sub(".$","",names(dup.fromDir));
dup.toDir <- table(paste0(data.link.df$to,data.link.df$toDir));
dup.toDir <- dup.toDir[dup.toDir > 1];
dup.toDir <- sub(".$","",names(dup.toDir));
dup.contigs <- union(dup.fromDir, dup.toDir);
link.lines.df <- data.link.df[(data.link.df$from %in% dup.contigs) |
                              (data.link.df$to   %in% dup.contigs),];
link.lines.contigs <- union(link.lines.df$from,link.lines.df$to);
paste(sub("tig0+","",link.lines.contigs),collapse=",");
length(link.lines.contigs);

#data.bilink.df <- unique(rbind(data.link.df[,c("from","to")],
#                               data.link.df[,c("to","from")]));
data.graph <- graph.data.frame(data.link.df[,c("from","to")]);
data.clusters <- clusters(data.graph);
data.clLengths <- tapply(data.clusters$membership, data.clusters$membership,
                         function(x){data.lengths[names(x)]});
data.clSizes <- data.clusters$csize;
data.unlinked <-
    names(data.lengths)[!(names(data.lengths) %in%
                          c(data.link.df$from,data.link.df$to))];
tigs.clSizes <- data.clusters$membership[names(data.lengths)];
names(tigs.clSizes) <- names(data.lengths);
tigs.clSizes[is.na(tigs.clSizes)] <- 0;

## barplot for length of contigs
pdf(sprintf("%s.pdf",gfaName), width=12, height=4);
par(mar=c(4.5,4.5,2,0.5));
options(scipen=15);
clBreaks <- c(0,paste(tail(sort(unique(floor((tigs.clSizes-1)/30))),-1)*30+1,
                      tail(sort(unique(floor((tigs.clSizes-1)/30))),-1)*30+30,
                      sep="-"));
bcols <- colorRampPalette(brewer.pal(11,"Spectral"))(length(clBreaks));
bcols <- colorRampPalette(brewer.pal(11,"Spectral"))(10);
myXRange <- c(1672,2048309); #range(data.lengths);
myYRange <- c(0,60); #range(data.lengths);
bar.data <- t(sapply(sort(unique(floor((tigs.clSizes-1)/30))), function(x){
    sequence.hist(data.lengths[floor((tigs.clSizes-1)/30) == x],
                  lengthRange=myXRange);
}));
ghost.data[[sprintf("%s",gfaName)]] <- colSums(bar.data);
ghost.data[[sprintf("%s",gfaName)]];
b.res <- barplot(bar.data/1000000, col=bcols, border=NA,
                 main=sprintf("%s",gfaName),
                 legend.text=clBreaks, las=2, ylab = "Aggregate length (Mb)",
                 xlab="Contig length", xaxt="n", ylim=c(0,60),
                 args.legend=list(x="topright", inset=0.05,
                                  title="Linked Contigs", cex=0.6));
for(gi in 1:length(ghost.data)){
    points(spline(x=b.res, y=ghost.data[[gi]]/1000000, n=100),
           type="l", lty=gi);
}
legend("topleft", lty=1:length(ghost.data), legend=names(ghost.data),
       inset=0.05, cex=0.6);
b.int <- b.res[2]-b.res[1];
axis(1, at=seq(head(b.res,1)-b.int/2, tail(b.res,1)+b.int/2, by=b.int),
     labels = unlist(strsplit(colnames(bar.data),";")), las=2);
#abline(h=1:100 * 5);
invisible(dev.off());
