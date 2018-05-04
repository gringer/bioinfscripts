#!/usr/bin/Rscript
genome <- "mm10";

features.meta <- read.csv("chrM_copies.csv");
features.meta$l2fc <- log(features.meta$qPct + 1);
features.meta$loc <- (features.meta$tS + features.meta$tE) / 2;

features.meta <- read.delim("geneAggregated_differential_coverage_4T1_WTvsp0.tsv");
features.meta$l2fc <- apply(features.meta[,c("fwd","rev")],1,function(x){head(max(abs(x)) * sign(x[abs(x) == max(abs(x))]),1)});
features.meta$loc <- (features.meta$Start + features.meta$End) / 2;

cytoURL <- paste(sprintf("http://genome.ucsc.edu/cgi-bin/hgTables?db=%s",
                         genome),
                 "hgta_outputType=primaryTable",
                 "hgta_regionType=genome",
                 "hgta_track=cytoBandIdeo",
                 "hgta_table=cytoBandIdeo",
                 "submit=submit",
                 "hgta_doTopSubmit=1", sep="&");

## Load in cytobands, strip "chr" from start, ignore mitochondrial sequence
cyto.df <- read.delim(cytoURL, stringsAsFactors=FALSE);
cyto.df$chr <- sub("^chr","",cyto.df$X.chrom);
cyto.df <- subset(cyto.df, !grepl("_", chr) & chr != "M");

## Convert chromosome names to integer version
chrs <- unique(cyto.df$chr);
chr.order <- 0:(length(chrs)-1);
chrs.numeric <- sort(as.numeric(chrs[grepl("[0-9]+",chrs)]));
chrs.alpha <- sort(chrs[!grepl("[0-9]+",chrs)]);
names(chr.order) <- c(chrs.numeric, chrs.alpha);
cyto.df$chrNum <- chr.order[cyto.df$chr];

## Set up scaling factors for chromosomes (number of bp per line)
nl <- 2;
chrs.per.line <- ceiling(length(chr.order) / nl);
bp.per.line <- ceiling(max(cyto.df$chromEnd/10^6*1.2))*10^6;

## Set up colours and hatching for cytobands
## stalk -- short arm of acrocentric chromosomes
## gvar -- heterochromatin
## acen -- centromeric region
cyto.df$col <-
    unlist(list(gpos100="black", gpos75="grey25", gpos66="grey33",
                gpos50="grey50",
                gpos25="grey75", gpos33="grey66", gneg="white",
                acen="blue", gvar="red", stalk="green")[cyto.df$gieStain]);
cyto.df$hang <- 0;
cyto.df$hang[cyto.df$gieStain == "acen"] <- 315;
cyto.df$hang[cyto.df$gieStain == "gvar"] <- 45;
cyto.df$dens <- NULL;
cyto.df$dens[cyto.df$gieStain == "acen"] <- 30;
cyto.df$dens[cyto.df$gieStain == "gvar"] <- 20;

## Determine extents of chromosomes
chrEnds <- tapply(cyto.df$chromEnd,cyto.df$chrNum,max);
chrStarts <- tapply(cyto.df$chromStart,cyto.df$chrNum,min);

## Determine adjusted locations for cytobands
cyto.df$x <- cyto.df$chrNum %% chrs.per.line;
cyto.df$yb <- (chrEnds[as.character(cyto.df$chrNum)]-cyto.df$chromStart) +
    bp.per.line * ((nl-1) - floor(cyto.df$chrNum / chrs.per.line));
cyto.df$yt <- (chrEnds[as.character(cyto.df$chrNum)]-cyto.df$chromEnd) +
    bp.per.line * ((nl-1) - floor(cyto.df$chrNum / chrs.per.line));

## Generate feature matrix with start/end locations of features
features.df <- features.meta;
#features.df$chr <- sub("^chr","",features.df$target);
features.df$chr <- sub("^chr","",features.df$Chr);
#features.df$start <- features.df$tS;
#features.df$end <- features.df$tE;
features.df$start <- features.df$Start;
features.df$end <- features.df$End;
features.df$col <- hsv(c(1/3,1/2,0)[-sign(features.df$l2fc)+2],
                       alpha=abs(features.df$l2fc)/
                           max(abs(range(features.df$l2fc))));

features.df <- subset(features.df, abs(l2fc) > 2);

## Create ideogram plot
svg("karyotype_mm10_Differential_Coverage.svg",width=19,height=10);
par(mar=c(0.1,0.1,0.1,0.1), cex = 2);
## set up plot extents
plot(NA, xlim = c(0,chrs.per.line),
     ylim = c(-bp.per.line*0.1,bp.per.line * nl), axes=FALSE,
     xlab="", ylab="");
## draw circles on chromosome ends
symbols(x=chr.order %% chrs.per.line+0.5,
        y=bp.per.line * ((nl-1) - floor(chr.order / chrs.per.line)),
        circles=rep(0.25,length(chrs)), add=TRUE, inches=FALSE, bg="grey");
symbols(x=chr.order %% chrs.per.line+0.5,
        y=bp.per.line * ((nl-1) - floor(chr.order / chrs.per.line))+chrEnds,
        circles=rep(0.25,length(chrs)), add=TRUE, inches=FALSE, bg="grey");
## Add in chromosome names
text(x=chr.order %% chrs.per.line+0.5,
     y=bp.per.line * ((nl-1.15) - floor(chr.order / chrs.per.line)),
     labels = names(chr.order));
## Cut out half of circles to leave two arcs
rect(xl=chr.order %% chrs.per.line+0.25, xr=chr.order %% chrs.per.line+0.75,
     yb=bp.per.line * ((nl-1) - floor(chr.order / chrs.per.line)),
     yt=bp.per.line * ((nl-1) - floor(chr.order / chrs.per.line)) + chrEnds,
     col="white", border=NA);
## Create cytobands
rect(xl=cyto.df$x+0.25, xr=cyto.df$x+0.75,
     yb=cyto.df$yb, yt=cyto.df$yt,
     density=cyto.df$dens, angle=cyto.df$hang,
     col = cyto.df$col, border=NA);
## Add border on ends of chromosomes
segments(x0=c(chr.order %% chrs.per.line+0.25,chr.order %% chrs.per.line+0.75),
     y0=bp.per.line * ((nl-1) - floor(chr.order / chrs.per.line)),
     y1=bp.per.line * ((nl-1) - floor(chr.order / chrs.per.line)) + chrEnds,
     col="black");
## ## Generate legend for hatched/coloured regions
## legend("topright", legend=c("acen","gvar","stalk"),
##        fill = c("blue","red","green"), inset=0.025);
## Add in feature information
for(l in 1:dim(features.df)[1]){
    xadj <- (chr.order[features.df$chr[l]] %% chrs.per.line);
    yadj <- bp.per.line * ((nl-1) - floor(chr.order[features.df$chr[l]] /
                                              chrs.per.line));
    y1 <- features.df$start[l]+yadj;
    y2 <- features.df$end[l]+yadj;
    if(features.df$l2fc[l] < 0){
        polygon(col=features.df$col[l],
                x=c(0.85,0.77,0.77,0.85) + xadj,
                y=c(y1-bp.per.line*0.02,y1,y2,y2+bp.per.line*0.02));
    } else {
        polygon(col=features.df$col[l],
                x=c(0.15,0.23,0.23,0.15) + xadj,
                y=c(y1-bp.per.line*0.02,y1,y2,y2+bp.per.line*0.02));
    }
}
legend("topright", legend=c("WT > ρ0", "ρ0 > WT"),
       fill = c("green","red"), inset=0.025);
graphics.off();
