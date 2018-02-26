#!/usr/bin/Rscript

argLoc <- 1;

usage <- function(){
  cat("usage: ./circularAnnotator.r <gff3 file> <genome size> [options]\n");
  cat("\nOther Options:\n");
  cat("-help                  : Only display this help message\n");
  cat("-size <width>x<height> : Size of the output image (in pixels)\n");
  cat("-type (png/pdf)        : Image type (default 'png')\n");
 cat("\n");
}

gffFileName <- "";
genomeSize <- -1;
type <- "png";
sizeX <- -1;
sizeY <- -1;

if(length(commandArgs(TRUE)) == 0){
      usage();
      quit(save = "no", status=0);
}

while(!is.na(commandArgs(TRUE)[argLoc])){
    arg <- commandArgs(TRUE)[argLoc];
    argLoc <- argLoc + 1;
    if(arg == "-help"){
      usage();
      quit(save = "no", status=0);
    } else if(arg == "-size"){
        arg <- unlist(strsplit(commandArgs(TRUE)[argLoc], "x"));
        argLoc <- argLoc + 1;
        sizeX <- as.numeric(arg[1]);
        sizeY <- as.numeric(arg[2]);
    } else if(arg == "-type"){
        arg <- commandArgs(TRUE)[argLoc];
        argLoc <- argLoc + 1;
        type <- arg;
    } else {
        if(file.exists(arg)){
            gffFileName <- arg;
        } else if ((genomeSize == -1) && grepl("^[0-9]+$", arg)){
            genomeSize <- as.numeric(arg);
        } else {
            cat("Error: Argument '", arg,
                "' is not understood by this program\n\n", sep="");
            usage();
            quit(save = "no", status=0);
        }
    }
}

if(sizeX == -1){
    if((type == "pdf") || (type == "svg")){
        sizeX = 8;
        sizeY = 8;
    } else {
        sizeX = 1600;
        sizeY = 1600;
    }
}

if(gffFileName == ""){
      usage();
      quit(save = "no", status=0);
}

if(genomeSize == -1){
    gffLines <- readLines(gffFileName, n=100);
    if(any(grepl("^..sequence-region", gffLines))){
        genomeSize <-
            as.numeric(sub("^.* ","",
                           gffLines[grep("^..sequence-region", gffLines)[1]]));
    } else {
        cat("No genome size found\n");
        usage();
        quit(save = "no", status=0);
    }
}

## Circular plot
library(plotrix);

draw.wedge <- function (x=0, y=0, rad1=0.5, rad2=1,
                       angle1 = deg1 * pi/180,
    angle2 = deg2 * pi/180, deg1 = 0, deg2 = 45, n = 0.05, col = NA,
    lwd = NA, ...)
{
    if (all(is.na(col)))
        col <- par("col")
    if (all(is.na(lwd)))
        lwd <- par("lwd")
    xylim <- par("usr")
    ymult <- getYmult()
    devunits <- dev.size("px")
    draw.wedge.0 <- function(x, y, rad1, rad2, angle1, angle2, n, col,
        lwd, ...) {
        delta.angle <- (angle2 - angle1)
        if (n != as.integer(n))
            n <- as.integer(1 + delta.angle/n)
        delta.angle <- delta.angle/n
        angles <- angle1 + seq(0, length = (n+1)) * delta.angle
        if (n > 1) {
            half.lwd.user <- (lwd/2) * (xylim[2] - xylim[1])/devunits[1]
            adj.angle = delta.angle * half.lwd.user/(2 * (mean(rad1,rad2) +
                half.lwd.user))
            angles[2:n] = angles[2:n] - adj.angle
        }
        p1x <- x + rad1 * cos(angles)
        p1y <- y + rad1 * sin(angles) * ymult
        p2x <- rev(x + rad2 * cos(angles))
        p2y <- rev(y + rad2 * sin(angles) * ymult)
        polygon(x=c(p1x,p2x), y=c(p1y, p2y), col = col, lwd = lwd, ...)
    }
    xy <- xy.coords(x, y)
    x <- xy$x
    y <- xy$y
    a1 <- pmin(angle1, angle2)
    a2 <- pmax(angle1, angle2)
    angle1 <- a1
    angle2 <- a2
    args <- data.frame(x, y, rad1, rad2, angle1, angle2, n, col,
        lwd, stringsAsFactors = FALSE)
    for (i in 1:nrow(args)) do.call("draw.wedge.0", c(args[i, ],
        ...))
    invisible(args)
}


chrom.length <- genomeSize;

data.features.df <-
    read.delim(gffFileName, header=FALSE,
               comment.char="#",
               col.names=c("chrom","source","type","start", "end",
                           "score", "strand", "phase", "attributes"),
               stringsAsFactors=FALSE);
data.features.df$name <- sub("^.*Name=([^;]*).*$","\\1",
                             data.features.df$attributes);
data.features.df$attributes <- NULL;
data.features.df <- subset(data.features.df,grepl("gene$",type));
data.features.df$name <- sub("^mt-","",data.features.df$name);
data.features.df$plotName <- data.features.df$name;
data.features.df$mid <- (data.features.df$start + data.features.df$end)/2;
data.features.df$plotName <- sub("tRNA-(...)_...","\\1",
                                 data.features.df$plotName);
data.features.df$plotName[data.features.df$plotName == "AT-rich region"] <- "";
data.features.df$plotName[data.features.df$plotName == "UTR"] <- "";
data.features.df$plotPos <- 0;
data.features.df$plotPos[grepl("(tRNA|^T[a-z]$)",data.features.df$name)] <- -1:1;

dateStr <- format(Sys.Date(), "%Y-%b-%d");

if(type == "pdf"){
    pdf(sprintf("circular_diagram_%s.pdf", dateStr), width=sizeX, height=sizeY);
} else if(type == "svg") {
    svg(sprintf("circular_diagram_%s.svg", dateStr), width=sizeX, height=sizeY,
        pointsize = 16 / (8 / min(sizeX, sizeY)));
} else {
    png(sprintf("circular_diagram_%s.png", dateStr), width=sizeX, height=sizeY,
        pointsize = 48 / (1600 / min(sizeX, sizeY)));
}
par(mar=c(0,0,0,0));
radial.plot(NA, radial.lim=c(0,16), show.grid=FALSE);
draw.wedge(x=0,y=0,rad1=13.5, rad2=15.5,
           angle1=0,
           angle2=pi*2,
           col="grey");
for(i in 1:nrow(data.features.df)){
    plotName <- data.features.df[i,"plotName"];
    draw.wedge(x=0,y=0,rad1=13.5, rad2=15.5,
             angle1=-data.features.df[i,"start"]*2*pi/chrom.length+pi/2,
               angle2=-data.features.df[i,"end"]*2*pi/chrom.length+pi/2,
               col=ifelse(nchar(plotName)==0,"grey",
                   ifelse(nchar(plotName)==1,"cornflowerblue","cornsilk")));
    if(nchar(plotName) > 1){
        arctext(plotName,
                middle=-data.features.df[i,"mid"]*2*pi/chrom.length+pi/2,
                radius=ifelse(nchar(plotName) == 1, 10, 14.5)
                + 1*data.features.df[i,"plotPos"],
                cex=0.9,
                col=ifelse(nchar(plotName)<3, "darkgreen","black"));
    }
}
draw.circle(x=0, y=0, radius=13);
for(p in 0:(genomeSize / 1000)){
        arctext(paste0(p),
                middle=-p*1000*2*pi/chrom.length+pi/2,
                radius=12, cex=0.7);
        draw.radial.line(center=c(0,0), start=12.75, end=13,
                         angle=-p*1000*2*pi/chrom.length+pi/2);
}
dummy <- dev.off();

