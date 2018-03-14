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

svg.paths <- xpathApply(svg.data, '//svg:path', fun=function(x){
    list(d=xmlGetAttr(x, "d"), style=xmlGetAttr(x, "style"));
});

cStart <- substr(svg.paths[[1]]$d,1,1);

pathChunks <- unlist(regmatches(svg.paths[[1]]$d,
                         gregexpr(paste0("([MmLl]?",
                                         "\\s*[0-9\\.]+",
                                         "[,\\s][0-9\\.]+|[zZ])"),
                                  svg.paths[[1]]$d)));
## not yet implemented:
## # HhVv - Horizontal and vertical lines
## # CcSs - Curve (cubic bezier)
## # QqTt - Curve (quadratic)
## # Aa   - Elliptical arc 
blankLineChunks <- grep("^[^MmLlHhVvCcSsQqTtAaZz]", pathChunks);
pathChunks[blankLineChunks] <- paste0("#",pathChunks[blankLineChunks]);
pathChunks <- sub("^(.)\\s+([0-9\\.]+)","\\1\\2", pathChunks);
path.df <- data.frame(chunk = pathChunks, stringsAsFactors=FALSE);
path.df$command <- substr(path.df$chunk,1,1);
path.df$remainder <- sub("^.","",path.df$chunk);
path.df$posX <- sub("^([\\s0-9\\.]+)(.*)$","\\1",path.df$remainder);
path.df$posY <- sub("^([\\s0-9\\.]+)(,|\\s+)","",path.df$remainder);
command.rle <- rle(path.df$command);
anonPoss <- which(command.rle$values == "#");
## [8.3.2] "If a moveto is followed by multiple pairs of coordinates, the
## subsequent pairs are treated as implicit lineto commands."
command.rle$values[anonPoss] <-
    ifelse(command.rle$values[anonPoss-1] == "M", "L",
           ifelse(command.rle$values[anonPoss-1] == "m", "l", "#"));
path.df$command <- inverse.rle(command.rle);
