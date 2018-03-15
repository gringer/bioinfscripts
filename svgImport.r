library(gridSVG);
library(grid);
library(XML);

x = c(0, 0.5, 1, 0.5)
y = c(0.5, 1, 0.5, 0)
grid.newpage()
grid.polygon(x,y, name="goodshape")
pat <- pattern(linesGrob(gp=gpar(col="black",lwd=3)),
  width = unit(5, "mm"), height = unit(5, "mm"),
  dev.width = 1, dev.height = 1)
# Registering pattern
registerPatternFill("pat", pat)
# Applying pattern fill
grid.patternFill("goodshape", label = "pat")
grid.export("test-pattern.svg")

svg.data <- xmlParse("dna_linear_DS_v4_backbone.svg");

xmlGetAttr(svg.data,"viewBox")


svginfo <- getNodeSet(svg.data, "//*[name()='svg']")[[1]];
svgVB <- xmlGetAttr(svginfo, "viewBox");
svgW <- xmlGetAttr(svginfo, "width");
svgH <- xmlGetAttr(svginfo, "height");
if(!is.null(svgVB)){
    svgVB <- as.numeric(strsplit(svgVB,"\\s+")[[1]]);
} else if(!is.null(svgW)){
    svgVB <- c(0,0,svgW,svgH);
}

svg.paths <- xpathApply(svg.data, '//svg:path', fun=function(x){
    list(d=xmlGetAttr(x, "d"), style=xmlGetAttr(x, "style"));
});

grid.newpage();

for(myPath in svg.paths){
    pathChunks <-
        unlist(regmatches(myPath$d,
                          gregexpr(paste0("([MmLl]?\\s*[0-9\\.]+",
                                          "[,\\s][0-9\\.]+|[zZ])"),
                                             myPath$d)));
    ## not yet implemented:
    ## # HhVv - Horiz / vert lines     # CcSs - Curve (cubic bezier)
    ## # QqTt - Curve (quadratic)      # Aa   - Elliptical arc 
    blankLineChunks <- grep("^[^MmLlHhVvCcSsQqTtAaZz]", pathChunks);
    pathChunks[blankLineChunks] <- paste0("#",pathChunks[blankLineChunks]);
    pathChunks <- sub("^(.)\\s+([0-9\\.]+)","\\1\\2", pathChunks);
    path.df <- data.frame(chunk = pathChunks, stringsAsFactors=FALSE);
    path.df$command <- substr(path.df$chunk,1,1);
    path.df$remainder <- sub("^.","",path.df$chunk);
    path.df$posX <- as.numeric(sub("^([\\s0-9\\.]+)(.*)$","\\1",
                                   path.df$remainder));
    path.df$posY <- as.numeric(sub("^([\\s0-9\\.]+)(,|\\s+)","",
                                   path.df$remainder));
    command.rle <- rle(path.df$command);
    anonPoss <- which(command.rle$values == "#");
    ## [8.3.2] "If a moveto is followed by multiple pairs of coordinates, the
    ## subsequent pairs are treated as implicit lineto commands."
    command.rle$values[anonPoss] <-
        ifelse(command.rle$values[anonPoss-1] == "M", "L",
        ifelse(command.rle$values[anonPoss-1] == "m", "l", "#"));
    path.df$command <- inverse.rle(command.rle);
    penX <- 0;  penY <- 0;  startX <- 0;  startY <- 0;
    absPoss <- matrix(NA,nrow=nrow(path.df), ncol=2);
    for(li in seq_along(path.df$command)){
        cmd <- path.df$command[li];
        px <- path.df$posX[li]; py <- path.df$posY[li];
        if((cmd == "m") || (cmd == "M")){
            startX <- px; startY <- py;
        } else if(cmd == "l"){
            px <- px + penX; py <- py + penY;
        } else if(cmd == "z"){
            px <- startX; py <- startY;
        }
        absPoss[li,] <- c(px, py);
        penX <- px; penY <- py;
    }
    path.df[,c("absX","absY")] <- absPoss;
    path.df$absX <- (path.df$absX - svgVB[1]) / (svgVB[3] - svgVB[1]);
    path.df$absY <- 1 - (path.df$absY - svgVB[2]) / (svgVB[4] - svgVB[2]);
    grid.polygon(x=path.df$absX, y=path.df$absY);
}
