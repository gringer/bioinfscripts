#!/usr/bin/Rscript

argLoc <- 1;

usage <- function(){
  cat("usage: ./seqmat.r <fasta file> <repeat unit size> [options]\n");
  cat("\nOther Options:\n");
  cat("-help               : Only display this help message\n");
  cat("-type (png/pdf)     : Image type (default 'png')\n");
  cat("-size <x>x<y>       : Image size (default 1200x1200 for png, 12x12 for PDF)\n");
  cat("-col (ef|cb)     : Colour palette (electrophoresis [default], colour-blind)\n");
 cat("\n");
}

fileName <- "";
rptSize <- -1;
type <- "png";
sizeX <- -1;
sizeY <- -1;
colType <- "ef";

seqRange <- FALSE;

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
    } else if(arg == "-col"){
        arg <- commandArgs(TRUE)[argLoc];
        argLoc <- argLoc + 1;
        colType <- arg;
    } else {
        if(grepl(":",arg)){
            subArgs <- unlist(strsplit(arg,"[:\\-]"));
            arg <- subArgs[1];
            seqRange <- subArgs[2:3];
        }
        if(file.exists(arg)){
            fileName <- arg;
        } else if ((rptSize == -1) && grepl("^[0-9]+$", arg)){
            rptSize <- as.numeric(arg);
        } else {
            cat("Error: Argument '", arg,
                "' is not understood by this program\n\n", sep="");
            usage();
            quit(save = "no", status=0);
        }
    }
}

if(rptSize == -1){
    cat("Error: repeat unit size has not been specified\n\n");
    usage();
    quit(save = "no", status=0);
}


inLines <- readLines(fileName);

inName <- substring(inLines[1],2);
inSeq <- c(A=1, C=2, G=3, T=4)[unlist(strsplit(paste(inLines[-1],collapse=""),""))];

if(seqRange[1] != FALSE){
    inSeq <- inSeq[seqRange[1]:seqRange[2]];
}

if(type == "png"){
    if(sizeX == -1){
        sizeX <- 1200;
    }
    if(sizeY == -1){
        sizeY <- 1200;
    }
} else if(type == "pdf"){
    if(sizeX == -1){
        sizeX <- 12;
    }
    if(sizeY == -1){
        sizeY <- 12;
    }
}

colPal <-
    if(colType == "cb"){
        c("#006400","#0000FF","#FFD700","#FF6347","#D3D3D3");
    } else {
        c("#006400","#0000FF","#FFD700","#8B0000","#D3D3D3");
    }

if(type == "png"){
    png("sequence_matrix.png", width=sizeX, height=sizeY, pointsize=24);
} else if(type == "pdf"){
    pdf("sequence_matrix.pdf", width=sizeX, height=sizeY, pointsize=16);
}
lis <- length(inSeq);
numLines <- floor(lis/rptSize + 1);
inSeq <- c(inSeq,rep(5,rptSize));
subSeq <- inSeq[1:(numLines*rptSize)];
par(mar=c(0.5,5,1,0.5), mgp=c(3.5,1,0));
image(x=1:rptSize, y=1:numLines-1, matrix(subSeq,nrow=rptSize),
      main=sprintf("%s%s (%0.3f kb, %d bases / line)", sub(" .*$","",inName),
                   ifelse(seqRange[1] == FALSE,"",
                          paste0(":",seqRange[1],"-",seqRange[2])),
                   lis/1000, rptSize), ylab="Base location",
      cex.main=0.8, xaxt="n", yaxt="n", useRaster=TRUE,
      col=colPal);
axis(2, at=round(seq(0, numLines, length.out=20)),
     labels=round(seq(0, numLines, length.out=20)) * rptSize+1,
     las=2, cex.axis=1);
dummy <- dev.off();

numLoops <- (length(subSeq) / rptSize);
##startCount <- rptSize / 1.2;
startCount <- rptSize;
startRadius <- 0.3;
endRadius <- 1.0;
##loopIncrement <- ((rptSize * 1.2) - (rptSize / 1.2)) / numLoops;


if(type == "png"){
    png("sequence_circle.png", width=max(sizeX,sizeY), height=max(sizeX,sizeY), pointsize=24);
} else if(type == "pdf"){
    pdf("sequence_circle.pdf", width=max(sizeX,sizeY), height=max(sizeX,sizeY), pointsize=16);
}
par(mar=c(1.5,1.5,0.5,1.5));
plot(NA, xlim=c(-1,1), ylim=c(-1,1), axes=FALSE, ann=FALSE);
mtext(sprintf("%s%s (%0.3f kb, %d bases / ring)", sub(" .*$","",inName),
              ifelse(seqRange[1] == FALSE,"",
                     paste0(":",seqRange[1],"-",seqRange[2])),
              lis/1000, rptSize), side=1, cex=1.5);
## Pre-population plot variables
## integrate(2*pi*r,r=startRadius..endRadius)
## => pi((endRadius)²-(startRadius)²)
dTot <- pi*(endRadius^2 - startRadius^2); ## total "distance" travelled
theta <- seq(0, (numLoops) * 2*pi, length.out=length(subSeq)); ## traversed angle
deg <- (theta / (2*pi)) * 360;
r <- seq(sqrt(startRadius), sqrt(endRadius),
         length.out=length(subSeq))^2; ## path radius
s <- pi * (r^2 - startRadius^2); ## traversed distance at each pos
ds <- c(s[2],diff(s)); ## distance difference at each pos

## par(mfrow=c(4,1));
## plot(s, main="s");
## plot(r, main="r");
## plot(deg, main="deg");
## plot(theta, main="theta");

## draw the spiral
points(rev(r) * cos(rev(theta)), rev(r) * sin(rev(theta)),
       col=paste0(colPal[rev(subSeq)],"A0"),
       pch=20, cex=rev(sqrt(r)) * (7/log(numLoops)));
legend("center", legend=c("A","C","G","T"), inset=0.2,
       fill=colPal[1:4],
       cex=ifelse(numLoops < 10, 1, 0.71));
invisible(dev.off());
