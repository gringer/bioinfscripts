#!/usr/bin/Rscript

inFileName <- commandArgs(TRUE)[1];

valToSci <- function(val, unit = ""){
    sci.prefixes <- c("", "k", "M", "G", "T", "P", "E", "Z", "Y");
    units <- rep(paste(sci.prefixes,unit,sep=""), each=3);
    logRegion <- floor(log10(val))+1;
    conv.units <- units[logRegion];
    conv.div <- 10^rep(0:(length(sci.prefixes)-1) * 3, each = 3)[logRegion];
    conv.val <- val / conv.div;
    conv.val[val == 0] <- 0;
    conv.units[val == 0] <- unit;
    return(ifelse(conv.units == "", conv.val, sprintf("%s %s",conv.val,conv.units)));
}

rankPlot <- function(data){
    options(scipen = 6);
    par(mar=c(6,6,1,1)+0.1, mgp = c(4.5,1,0), las = 1, mfrow=c(1,2),
        cex.lab=1.5, bg="#FFFFFFFF");
    dataRange <- c(1,max(data$maxrank));
    plot(maxrank ~ minrank, data=data,
         log="xy", xlim=dataRange, ylim=dataRange, ylab="Maximum Rank",
         cex=2, xlab="Minimum Rank (All Markers)", pch=21,
         axes=FALSE, bg="#4682B480", col="#4682B4");
    abline(h=10^(0:max(log10(dataRange))), v=10^(0:max(log10(dataRange))),
           col="#00000020", lty="dashed");
    abline(h=(nrow(data) * 0.05), col="red", lwd=2, lty="dotted");
    box(lwd=3);
    drMax <- max(log10(dataRange));
    axis(1, at= 10^(0:drMax), las=2, lwd=3, cex.axis=1.5, labels=valToSci(10^(0:drMax)));
    axis(2, at= 10^(0:drMax), lwd=3, cex.axis=1.5, labels=valToSci(10^(0:drMax)));
    axis(1, at= rep(1:9, each=drMax+1) * 10^(0:drMax), labels=FALSE);
    axis(2, at= rep(1:9, each=drMax+1) * 10^(0:drMax), labels=FALSE);
    plot(maxrank ~ minrank, data=subset(data,
                                maxrank < (nrow(data) * 0.05)),
         log="xy", xlim=dataRange, ylim=dataRange, ylab="Maximum Rank",
         cex=2, xlab=sprintf("Minimum Rank (%d Consistent Markers)",
                             sum(data$maxrank < (nrow(data) * 0.05))),
         pch=21, axes=FALSE, bg="#4682B480", col="#4682B4");
    abline(h=10^(0:max(log10(dataRange))), v=10^(0:max(log10(dataRange))),
           col="#00000020", lty="dashed");
    box(lwd=3);
    drMax <- max(log10(dataRange));
    axis(1, at= 10^(0:drMax), las=2, lwd=3, cex.axis=1.5, labels=valToSci(10^(0:drMax)));
    axis(2, at= 10^(0:drMax), lwd=3, cex.axis=1.5, labels=valToSci(10^(0:drMax)));
    axis(1, at= rep(1:9, each=drMax+1) * 10^(0:drMax), labels=FALSE);
    axis(2, at= rep(1:9, each=drMax+1) * 10^(0:drMax), labels=FALSE);
}


cs.stats <-
    c(1,247249719,224999719,
      2,242951149,237712649,
      3,199501827,194704827,
      4,191273063,187297063,
      5,180857866,177702766,
      6,170899992,167273992,
      7,158821424,154952424,
      8,146274826,142612826,
      9,140273252,120143252,
      10,135374737,131624737,
      11,134452384,131130853,
      12,132349534,130303534,
      13,114142980,95559980,
      14,106368585,88290585,
      15,100338915,81341915,
      16,88827254,78884754,
      17,78774742,77800220,
      18,76117153,74656155,
      19,63811651,55785651,
      20,62435964,59505253,
      21,46944323,34171998,
      22,49691432,34851332,
      23,154913754,151058754, # X
      24,57772954,25652954,   # Y
      25,16571,16571);        # MT
dim(cs.stats) <- c(3,25);
cs.stats <- as.data.frame(t(cs.stats));
colnames(cs.stats) <- c("Chromosome","Assembled","Sequenced");
rownames(cs.stats) <- cs.stats$Chromosome;
cs.stats$label <- cs.stats$Chromosome;
cs.stats$label[23:25] <- c("X","Y","MT");
cs.stats$startPoint <-
    c(0,cumsum(as.numeric(cs.stats$Assembled)))[-length(cs.stats$Assembled)];
cs.stats$colour <- rep(c("#4682B480","#5F9EA080"),13)[1:nrow(cs.stats)];

chrPlot <- function(parent, plotData = NULL){
    useLog <- ("maxrank" %in% colnames(parent));
    parent.ylim <- c(min(-log(parent$maxrank)),0);
    par(mar=c(2.5,6,1,1)+0.1, mgp = c(3.5,1,0), las = 1,
        cex.axis=1.5, cex.lab=1.5, bg="#FFFFFFFF", mfrow = c(2,1));
    assembled.half <-
        cs.stats$startPoint[which.min(cs.stats$startPoint <
                                      max(parent$adjPos)/2)];
    drawMarkers <-
        if(!is.null(plotData)){
            subset(parent, (adjPos < assembled.half) &
                           !(marker %in% plotData$marker));
        } else {
            subset(parent, (adjPos < assembled.half));
        }
    drawMarkers$col <-cs.stats[drawMarkers$chrom,"colour"];
    sampRows <- 1:nrow(drawMarkers);
    plot(-log(maxrank) ~ adjPos,
         data=drawMarkers[sampRows,],
         xlim=c(0,assembled.half),
         ylim=parent.ylim, xlab="",
         ylab=if(useLog){expression(-log(adjRank))},
         cex=1, pch=21, bg=drawMarkers$col[sampRows],
         col=sub("..$","",drawMarkers$col[sampRows]),
         axes=FALSE);
    if(!is.null(plotData)){
        points(-log(maxrank) ~ (adjPos),
               data=subset(plotData, adjPos <= max(drawMarkers$adjPos)),
               cex=1, pch=21, bg="#FA807280", col="#FA8072");
        box(lwd=1);
    }
    ## print tick marks
    axis(side = 1, at=cs.stats$startPoint[c(-25)], labels = FALSE, lwd = 2,
         cex.axis = 2);
    ## print labels
    axis(side = 1, at=(cs.stats$startPoint[c(-24,-25)] +
                       cs.stats$startPoint[c(-1,-25)]) / 2,
         labels = cs.stats$label[c(-24,-25)],
         las = 1, tick = FALSE, lwd = 2, cex.axis = 1);
    ## print y axis
    axis(2, las=2);
    drawMarkers <-
        if(!is.null(plotData)){
            subset(parent, (adjPos >= assembled.half) &
                           !(marker %in% plotData$marker));
        } else {
            subset(parent, (adjPos >= assembled.half));
        }
    drawMarkers$col <-cs.stats[drawMarkers$chrom,"colour"];
    sampRows <- 1:nrow(drawMarkers);
    plot(-log(maxrank) ~ adjPos,
         data=drawMarkers[sampRows,],
         xlim=c(assembled.half, assembled.half*2),
         ylim=parent.ylim, xlab="",
         ylab=expression(-log(adjRank)),
         cex=1, pch=21, bg=drawMarkers$col[sampRows],
         col=sub("..$","",drawMarkers$col[sampRows]),
         axes=FALSE);
    if(!is.null(plotData)){
        points(-log(maxrank) ~ (adjPos),
               data=subset(plotData, adjPos >= min(drawMarkers$adjPos)),
               cex=1, pch=21, bg="#FA807280", col="#FA8072");
        box(lwd=1);
    }
    ## print tick marks
    axis(side = 1, at=cs.stats$startPoint[c(-25)], labels = FALSE, lwd = 2,
         cex.axis = 2);
    ## print labels
    axis(side = 1, at=(cs.stats$startPoint[c(-24,-25)] +
                       cs.stats$startPoint[c(-1,-25)]) / 2,
         labels = cs.stats$label[c(-24,-25)],
         las = 1, tick = FALSE, lwd = 2, cex.axis = 1);
    ## print y axis
    axis(2, las=2);
}

inFileDir <- dirname(inFileName);
inFileBase <- sub(".csv(.gz)?$","",basename(inFileName));

cat("Loading data...");
data.df <- read.csv(inFileName);
cat(" done\n");

cat("Getting chromosome parameters...");
data.df$chrom <- as.numeric(sub("^.*?_([^_]*?)_[^_]*$","\\1",data.df$marker));
data.df$pos <- as.numeric(sub("^.*_","",data.df$marker));
data.df <- subset(data.df, chrom != 0);
data.df$startPoint <- cs.stats[data.df$chrom,]$startPoint;
data.df <- subset(data.df, !is.na(startPoint));
data.df <- subset(data.df, (chrom > 0) & (chrom < 24)); ## exclude unknown, Y, MT chromosome
data.df$adjPos <- data.df$startPoint + data.df$pos;
consistent.data.df <- subset(data.df, (maxrank < (nrow(data.df) * 0.05)));
cat(" done\n");

cat("Creating rank plot...");
rankFileName <- sprintf("%s/rankPlot_%s.png", inFileDir, inFileBase);
png(rankFileName,
    width=2400, height=1200,
    pointsize=24, antialias="gray");
rankPlot(data = data.df);
dummy <- graphics.off();
cat(sprintf(" done [%s]\n", rankFileName));

cat("Creating genome plot...");
genomeFileName <- sprintf("%s/genomePlot_%s.png", inFileDir, inFileBase);
png(genomeFileName,
    width=3500, height=1500,
    pointsize=36, antialias="gray");
chrPlot(data.df, consistent.data.df);
dummy <- dev.off();
cat(sprintf(" done [%s]\n", genomeFileName));
