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
    ## see https://stackoverflow.com/questions/4246077/matching-numbers-with-regular-expressions-only-digits-and-commas/4247184
    numRE <- "(([+-]?)(?=\\d|\\.\\d)\\d*(\\.\\d*)?([Ee]([+-]?\\d+))?)";
    pathChunks <-
        unlist(regmatches(myPath$d,
                          gregexpr(paste0("([MmLl]?\\s*",numRE,
                                          "[,\\s]",numRE,"|[zZ])"), perl=TRUE,
                                   myPath$d)));
    ## not yet implemented:
    ## # HhVv - Horiz / vert lines     # CcSs - Curve (cubic bezier)
    ## # QqTt - Curve (quadratic)      # Aa   - Elliptical arc 
    blankLineChunks <- grep("^[^MmLlHhVvCcSsQqTtAaZz]", pathChunks);
    pathChunks[blankLineChunks] <- paste0("#",pathChunks[blankLineChunks]);
     ## remove filler whitespace (if any)
    pathChunks <- sub(paste0("^(.)\\s+",numRE),"\\1\\2", pathChunks, perl=TRUE);
    path.df <- data.frame(chunk = pathChunks, stringsAsFactors=FALSE);
    path.df$command <- substr(path.df$chunk,1,1);
    path.df$remainder <- sub("^.","",path.df$chunk);
    path.df$posX <- as.numeric(sub(paste0("^",numRE,"(.*)$"),"\\1",
                                   path.df$remainder, perl=TRUE));
    path.df$posY <- as.numeric(sub(paste0("^",numRE,"(,|\\s+)"),"",
                                   path.df$remainder, perl=TRUE));
    #if(grepl("e",myPath$d)){
    #    break;
    #}
    command.rle <- rle(path.df$command);
    anonPoss <- which(command.rle$values == "#");
    ## [8.3.2] "If a moveto is followed by multiple pairs of coordinates, the
    ## subsequent pairs are treated as implicit lineto commands."
    command.rle$values[anonPoss] <-
        ifelse(command.rle$values[anonPoss-1] == "M", "L",
        ifelse(command.rle$values[anonPoss-1] == "m", "l", "#"));
    path.df$command <- inverse.rle(command.rle);
    #if(path.df$command[1] == "M"){
    #    next;
    #}
    penX <- 0;  penY <- 0;  startX <- 0;  startY <- 0;
    absPoss <- matrix(NA,nrow=nrow(path.df), ncol=2);
    for(li in seq_along(path.df$command)){
        cmd <- path.df$command[li];
        px <- path.df$posX[li]; py <- path.df$posY[li];
        if((cmd == "m") || (cmd == "M")){
            startX <- px; startY <- py;
        } else if(cmd == "l"){
            px <- px + penX; py <- py + penY;
        } else if((cmd == "z") || (cmd == "Z")){
            px <- startX; py <- startY;
        }
        absPoss[li,] <- c(px, py);
        penX <- px; penY <- py;
    }
    path.df[,c("absX","absY")] <- absPoss;
    path.df$absX <- (path.df$absX - svgVB[1]) / (svgVB[3] - svgVB[1]);
    path.df$absY <- 1 - (path.df$absY - svgVB[2]) / (svgVB[4] - svgVB[2]);
    pathStyle <- myPath$style;
    pathFill <- sub("^fill:(.*)[;$]","\\1",
                    regmatches(pathStyle,gregexpr("fill:.*?[;$]",pathStyle)));
    pathStroke <- sub("^stroke:(.*)[;$]","\\1",
                      regmatches(pathStyle,gregexpr("stroke:.*?[;$]",pathStyle)));
    if((tail(path.df$command,1) == "z") || (tail(path.df$command,1) == "Z")){
        grid.polygon(x=head(path.df$absX,-1), y=head(path.df$absY,-1),
                     gp=gpar(col=pathStroke,fill=pathFill));
    } else {
        grid.lines(x=path.df$absX, y=path.df$absY,
                     gp=gpar(col=pathStroke,fill=pathFill));
    }
}
