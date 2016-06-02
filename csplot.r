#!/usr/bin/Rscript

## csplot.r -- Generates a genome-wide plot of data from an input file
## that includes chromosome, location, and value.

## The output is similar to plots used in the WTCCC results paper
## (i.e. dark blue/light blue for different chromosomes), but is split
## over three lines.

## Author: David Eccles (gringer), 2009 <programming@gringer.org>

usage <- function(){
  cat("usage: ./csplot.r <file> [options]\n");
  cat("\nOther Options:\n");
  cat("-help              : Only display this help message\n");
  cat("-threshold <value> : Only display data greater than <value>\n");
  cat("-pointsize <value> : Multiplier for size of points in graph\n");
  cat("-invert            : Invert values (lowest value at top of graph)\n");
  cat("-transparent       : Use slightly transparent points\n");
  cat("-normlimit         : Limit value display to a reasonable normal distribution\n");
  cat("-limit <value>     : Trim values greater than <limit>\n");
  cat("-keep <value>      : Keep a proportion of values below the cutoff value\n");
  cat("-window <value>    : Calculate running median across window of size <value>\n");
  cat("-label <string>    : Label for Y axis\n");
  cat("\n");
}

dataFile <- FALSE; # marker chromosome location <other fields> value
valThreshold <- 0;
useNormDistForMax <- FALSE;
invertAxis <- FALSE;
trimLimit <- FALSE;
transparent <- FALSE;
sizeMul <- 2;
filterKeepProp <- 0; # keep this proportion of markers below threshold cutoff
filterRandom <- (filterKeepProp > 0); # keep a random sampling of markers below threshold cutoff
valueText <- expression(chi^2~Value);
windowSize <- 1;

argLoc <- 1;
while(!is.na(commandArgs(TRUE)[argLoc])){
  if(file.exists(commandArgs(TRUE)[argLoc])){ # file existence check
    if(dataFile == FALSE){
      dataFile <- commandArgs(TRUE)[argLoc];
    } else{
        cat("Error: More than one input file specified\n",dataFile);
        usage();
        quit(save = "no", status=1);
    }
  } else {
    parsed <- FALSE;
    if(commandArgs(TRUE)[argLoc] == "-help"){
      usage();
      parsed <- TRUE;
      quit(save = "no", status=0);
    }
    if(commandArgs(TRUE)[argLoc] == "-threshold"){
      valThreshold <- as.numeric(commandArgs(TRUE)[argLoc+1]);
      argLoc <- argLoc + 1;
      cat("Setting value threshold to ",valThreshold,"\n",sep="");
      parsed <- TRUE;
    }
    if(commandArgs(TRUE)[argLoc] == "-pointsize"){
      sizeMul <- as.numeric(commandArgs(TRUE)[argLoc+1]);
      argLoc <- argLoc + 1;
      cat("Setting point size multiplier to ",sizeMul,"\n",sep="");
      parsed <- TRUE;
    }
    if(commandArgs(TRUE)[argLoc] == "-invert"){
      invertAxis <- TRUE;
      cat("Using normal distribution to set upper graph limit\n");
      parsed <- TRUE;
    }
    if(commandArgs(TRUE)[argLoc] == "-transparent"){
      transparent <- TRUE;
      cat("Points will have partial transparency\n");
      parsed <- TRUE;
    }
    if(commandArgs(TRUE)[argLoc] == "-normlimit"){
      useNormDistForMax <- TRUE;
      cat("Using normal distribution to set upper graph limit\n");
      parsed <- TRUE;
    }
    if(commandArgs(TRUE)[argLoc] == "-limit"){
      trimLimit <- as.numeric(commandArgs(TRUE)[argLoc+1]);
      argLoc <- argLoc + 1;
      cat("Values greater than ",trimLimit," will be marked at graph extent\n",sep="");
      parsed <- TRUE;
    }
    if(commandArgs(TRUE)[argLoc] == "-keep"){
      filterKeepProp <- as.numeric(commandArgs(TRUE)[argLoc+1]);
      argLoc <- argLoc + 1;
      filterRandom <- (filterKeepProp > 0);
      cat("Including ",filterKeepProp,"x markers below threshold value\n",sep="");
      parsed <- TRUE;
    }
    if(commandArgs(TRUE)[argLoc] == "-window"){
      windowSize <- as.numeric(commandArgs(TRUE)[argLoc+1]);
      argLoc <- argLoc + 1;
      cat("Setting running median window size to ",windowSize,"\n",sep="");
      parsed <- TRUE;
    }
    if(commandArgs(TRUE)[argLoc] == "-label"){
      valueText <- commandArgs(TRUE)[argLoc+1];
      argLoc <- argLoc + 1;
      cat("Setting value text to ",valueText,"\n",sep="");
      parsed <- TRUE;
    }
    if(!parsed){
      cat("Error: unknown argument '",commandArgs(TRUE)[argLoc],"'\n",sep="");
    }
  }
  argLoc <- argLoc + 1;
}

if(dataFile == FALSE){
  cat("Error: No valid file given\n\n");
  usage();
  quit(save = "no", status=1);
}

## Chromosome Statistics -- from http://genome.ucsc.edu/goldenPath/stats.html#hg18
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
cs.stats$startPoint <- c(0,cumsum(as.numeric(cs.stats$Assembled)))[-length(cs.stats$Assembled)];

# takes about 10s on melinus with T1D data
cat("Reading data file...", file = stderr());
retVal = grep(pattern = "^\\D+$", x = readLines(dataFile, n = 1), perl = TRUE);
marker.statistics <- NULL;
markers.statistics <- read.table(dataFile, header = (length(retVal) > 0));
colnames(markers.statistics)[1] <- "Marker";
if(length(grep("/",markers.statistics[1,2])) > 0){
    colnames(markers.statistics)[2] <- "Mutation";
    colnames(markers.statistics)[3] <- "Chromosome";
    colnames(markers.statistics)[4] <- "Location";
} else {
    colnames(markers.statistics)[2] <- "Chromosome";
    colnames(markers.statistics)[3] <- "Location";
}
#rownames(markers.statistics) <- markers.statistics$Marker;
colnames(markers.statistics)[length(colnames(markers.statistics))] <- "Value";
markers.statistics$Chromosome <- sub("chr","",markers.statistics$Chromosome,
                                     ignore.case = TRUE);
markers.statistics$Chromosome <- sub("X","23",markers.statistics$Chromosome,
                                     ignore.case = TRUE);
markers.statistics$Chromosome <- sub("Y","24",markers.statistics$Chromosome,
                                     ignore.case = TRUE);
markers.statistics$Chromosome <- sub("MT?","25",markers.statistics$Chromosome,
                                     ignore.case = TRUE);
markers.statistics$Chromosome <- as.numeric(markers.statistics$Chromosome);
if(transparent){
    markers.statistics$Colour <- rep(c("#0000FF80","#00A0FF80"),13)[markers.statistics$Chromosome];
} else {
    markers.statistics$Colour <- rep(c("#0000FF","#00A0FF"),13)[markers.statistics$Chromosome];
}
markers.statistics$startPoint <- cs.stats[markers.statistics$Chromosome,]$startPoint;
cat("done!\n", file = stderr());

cat("Calculating limits...", file = stderr());
maxVal <- Inf;
if(useNormDistForMax){ # Filters out abnormally large spikes, redrawn as triangles in graph
  maxVal <- 4*qnorm(p = 1 - 1/dim(markers.statistics)[1], mean = mean(markers.statistics$Value), sd = sd(markers.statistics$Value));
}

if(!(trimLimit == FALSE)){
  maxVal <- trimLimit;
}
## valRange <- c(min(markers.statistics$Value, na.rm = TRUE),
##               min(maxVal, max(markers.statistics$Value, na.rm = TRUE)));
valRange <- c(min(c(0,markers.statistics$Value), na.rm = TRUE),
              min(maxVal, max(markers.statistics$Value, na.rm = TRUE)));
markers.statistics$randVal <- runif(dim(markers.statistics)[1]);
markers.filtered <- subset(markers.statistics, (Value >= valThreshold) | (filterRandom & (randVal < filterKeepProp)));
if(transparent){
    markers.filtered$Colour[markers.filtered$Value > valRange[2]] <-
        rep(c("#80000080","#A0404080"),13)[markers.filtered$Chromosome[markers.filtered$Value > valRange[2]]];
} else {
    markers.filtered$Colour[markers.filtered$Value > valRange[2]] <-
        rep(c("#800000","#A04040"),13)[markers.filtered$Chromosome[markers.filtered$Value > valRange[2]]];
}
markers.filtered$Value[markers.filtered$Value > valRange[2]] <- valRange[2]+0.01;
cat("done!\n", file = stderr());

cat("Generating png...", file = stderr());
X11.options(antialias="gray");
png("chromosome_plot.png", width = 1280, height = 720, pointsize=12,
    antialias="gray");
#svg("chromosome_plot.svg", width = 11, height = 8);
#Cairo_png("chromosome_plot.png", width = 11, height = 8);
if(!(is.expression(valueText)) && (valueText == FALSE)){
  par(mfrow = c(3,1), mar = c(4,5,1,1));
} else {
  par(mfrow = c(3,1), mar = c(4,7,1,1));
}
assembled.third <- cs.stats$startPoint[23]/3;
drawMarkers <- subset(markers.filtered,
                      ((startPoint + Location) >= assembled.third * 0) &
                      ((startPoint + Location) <= assembled.third * 1));
pchStyle <- (drawMarkers$Value > valRange[2]) + 16; # circles within range, triangles outside
plot(runmed(drawMarkers$startPoint+drawMarkers$Location, windowSize),
     drawMarkers$Value, col = drawMarkers$Colour,
     pch = pchStyle, cex = 0.5 * sizeMul, xaxt = "n", frame.plot = FALSE,
     ylab = "", ylim = valRange, xlim = c(0, assembled.third * 1),
     cex.axis = 2, lwd = 2, cex.lab = 2, xaxs = "i",
     xlab = "", las = 1);
## print tick marks
axis(side = 1, at=cs.stats$startPoint[c(-24,-25)], labels = FALSE, lwd = 2);
## print labels
axis(side = 1, at=(cs.stats$startPoint[c(-23,-24,-25)] +
       cs.stats$startPoint[c(-1,-24,-25)]) / 2, srt = 2,
     labels = seq(1,22), las = 1, tick = FALSE, lwd = 2, cex.axis = 2);
drawMarkers <- subset(markers.filtered,
                      ((startPoint + Location) >= assembled.third * 1) &
                      ((startPoint + Location) <= assembled.third * 2));
pchStyle <- (drawMarkers$Value > valRange[2]) + 16; # circles within range, triangles outside
plot(runmed(drawMarkers$startPoint+drawMarkers$Location, windowSize),
     drawMarkers$Value, col = drawMarkers$Colour,
     pch = pchStyle, cex = 0.5 * sizeMul, xaxt = "n", frame.plot = FALSE,
     ylab = valueText, ylim = valRange, xlim = c(assembled.third * 1, assembled.third * 2),
     cex.axis = 2, lwd = 5, cex.lab = 2, xaxs = "i",
     xlab = "", las = 1, mgp = c(5,1,0));
## print tick marks
axis(side = 1, at=cs.stats$startPoint[c(-24,-25)], labels = FALSE, lwd = 2,
     cex.axis = 2);
## print labels
axis(side = 1, at=(cs.stats$startPoint[c(-23,-24,-25)] +
       cs.stats$startPoint[c(-1,-24,-25)]) / 2,
     labels = seq(1,22), las = 1, tick = FALSE, lwd = 2, cex.axis = 2);
drawMarkers <- subset(markers.filtered,
                      ((startPoint + Location) >= assembled.third * 2) &
                      ((startPoint + Location) <= assembled.third * 3));
pchStyle <- (drawMarkers$Value > valRange[2]) + 16; # circles within range, triangles outside
plot(runmed(drawMarkers$startPoint+drawMarkers$Location, windowSize),
     drawMarkers$Value, col = drawMarkers$Colour,
     pch = pchStyle, cex = 0.5 * sizeMul, xaxt = "n", frame.plot = FALSE,
     ylab = "", ylim = valRange, xlim = c(assembled.third * 2, assembled.third * 3),
     cex.axis = 2, cex.lab = 2, xaxs = "i",
     xlab = paste("Chromosome location",
       if((valThreshold > 0) && (!filterRandom)){paste(", Value >", valThreshold)}, sep = ""), las = 1);
## print tick marks
axis(side = 1, at=cs.stats$startPoint[c(-24,-25)], labels = FALSE, lwd = 2,
     cex.axis = 2);
## print labels
axis(side = 1, at=(cs.stats$startPoint[c(-23,-24,-25)] +
       cs.stats$startPoint[c(-1,-24,-25)]) / 2,
     labels = seq(1,22), las = 1, tick = FALSE, lwd = 2, cex.axis = 2);
dummy <- dev.off();
cat("done!\n", file = stderr());
