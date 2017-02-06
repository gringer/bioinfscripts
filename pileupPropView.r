#!/usr/bin/Rscript

pileupFileName <- commandArgs(TRUE)[1];
pos.start <- as.numeric(commandArgs(TRUE)[2]);
pos.end <- as.numeric(commandArgs(TRUE)[3]);
                        
data.prop.df <- subset(
    read.csv(pileupFileName),
    (Position >= pos.start) & (Position <= pos.end));

prop.plot <- function(data.df){
    par(mar=c(5,5,0.5,1));
    res <- barplot(t(as.matrix(data.df[,c("A","C","G","T","d","pR","i")]) *
                     data.df$Coverage),
                   ylim=c(0,max(data.df$Coverage*1.4)),
                   xlim=c(0,nrow(data.df)*1.1),
                   xaxt="n", xlab = "Mitochondrial Genome Location",
                   ylab="Read Coverage", border=NA, space=0,
                   col=c("darkgreen","blue","black","red",
                         "steelblue","grey90","grey60"));
    legend("right", horiz=FALSE, legend=c("A","C","G","T","Del","Ref","Ins"),
           fill=c("darkgreen","blue","black","red",
                  "steelblue","grey90","grey60"));
    tckPoss <- pretty(data.df$Position);
    axis(1, at=res[match(tckPoss, data.df$Position)],
         labels=tckPoss);
}

png(paste0(sub("\\..*$","",pileupFileName),"_",pos.start,"-",pos.end,".png"),
    width=1366, height=718, pointsize=24);
prop.plot(data.prop.df);
dummy <- dev.off();
