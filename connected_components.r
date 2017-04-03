#!/usr/bin/Rscript

## Find clusters of contigs from a GFA produced by Canu that have a total length of >1Mb
## Identify sub-paths that have link identities >90%

library(igraph);
library(Biostrings);

idThreshold <- 0.9; ## identity threshold for overlap (don't trust anything below this)

gfaName <- "NbL5_ONTDECAF_t1.contigs.cleaned.gfa";
faName <- paste0(getwd(),"/",sub("\\.gfa$", ".fasta", gfaName));
gfa.sequences <- readDNAStringSet(faName);
names(gfa.sequences) <- sub(" .*$","",names(gfa.sequences));

data.lines <- readLines(gfaName);

data.lengths <- sapply(strsplit(grep("^S", data.lines, value=TRUE),"\t"),
                          function(x){
                              val <- as.numeric(substring(x[4],6));
                              names(val) <- x[2];
                              val;
                          });
data.linklist <- strsplit(grep("^L", data.lines, value=TRUE),"\t");

## rbind.fill doesn't work here

data.link.df <- data.frame(t(sapply(data.linklist,c)), stringsAsFactors=FALSE);
colnames(data.link.df) <- c("type","from","fromDir",
                            "to","toDir","cigar","flags");
invertLines <- (data.link.df$fromDir == "-") & (data.link.df$toDir == "-");
data.link.df[invertLines,c("from","to")] <- data.link.df[invertLines,c("to","from")];
data.link.df[invertLines,"fromDir"] <- "+";
data.link.df[invertLines,"toDir"] <- "+";
data.link.df <- unique(data.link.df);

data.graph <- graph.data.frame(data.link.df[,c("from","to")]);

data.clusters <- clusters(data.graph);
data.clLengths <- tapply(data.clusters$membership, data.clusters$membership,
                         function(x){data.lengths[names(x)]});

## pick clusters over 1Mb
data.largeClusters <- Filter(function(x){sum(x) >= 1e6}, data.clLengths);
data.largeClContigs <- unique(unlist(sapply(data.largeClusters,names)));
data.largeClLinks <- subset(data.link.df, (from %in% data.largeClContigs) | (to %in% data.largeClContigs));

for(clName in names(data.largeClusters)){
    oldwd <- getwd();
    cat(clName,"\n",sep="");
    dirName <- tempfile(sprintf("GFAcluster_%s_",clName));
    dir.create(dirName);
    setwd(dirName);
    clusterFastaName <- sprintf("cluster_%s_seqs.fa", clName);
    clusterComps <- names(data.largeClusters[[clName]]);
    system2("samtools", args=c("faidx", faName, clusterComps),
            stdout=clusterFastaName);
    cat(" Generating LAST index... ");
    system2("lastdb", args=rep(clusterFastaName,2));
    cat("done\n");
    cat(" Overlap mapping... ");
    overlapName <- sprintf("overlap_cluster_%s_seqs.tsv", clName);
    system2("lastal", args=c(unlist(strsplit("-P10 -f BlastTab+ -E 1e-50 -l 50 -k 50"," ")),
                             rep(clusterFastaName,2)), stdout=overlapName);
    cat("done\n");
    mapping.df <- read.delim(overlapName, comment.char="#", stringsAsFactors=FALSE,
                             col.names=c("qid","tid","pctid","alnLen","mismatch","gapOpen",
                                         "qStart","qEnd","tStart","tEnd", "eval", "bscore",
                                         "qLen", "tLen"));
    ## remove matches of a sequence with itself
    mapping.df <- subset(mapping.df, qid != tid);
    ## filter to select out only overlap sequences
    mapping.df <- subset(mapping.df,
                         (qStart == 1 | qEnd == 1 | qStart == qLen | qEnd == qLen) &
                         (tStart == 1 | tEnd == 1 | tStart == tLen | tEnd == tLen));
    ## filter to remove low percent identity
    mapping.df <- subset(mapping.df, pctid > (idThreshold * 100));
    ## use direction-relative notation (as per link graph)
    mapping.df$qDir <- ifelse(mapping.df$qStart > mapping.df$qEnd,"-","+");
    mapping.df$tDir <- ifelse(mapping.df$tStart > mapping.df$tEnd,"-","+");
    mapping.df$qStart[mapping.df$qDir == "-"] <-
        mapping.df$qLen[mapping.df$qDir == "-"] - mapping.df$qStart[mapping.df$qDir == "-"] + 1;
    mapping.df$qEnd[mapping.df$qDir == "-"] <-
        mapping.df$qLen[mapping.df$qDir == "-"] - mapping.df$qEnd[mapping.df$qDir == "-"] + 1;
    mapping.df$tStart[mapping.df$tDir == "-"] <-
        mapping.df$tLen[mapping.df$tDir == "-"] - mapping.df$tStart[mapping.df$tDir == "-"] + 1;
    mapping.df$tEnd[mapping.df$tDir == "-"] <-
        mapping.df$tLen[mapping.df$tDir == "-"] - mapping.df$tEnd[mapping.df$tDir == "-"] + 1;
    mapping.df$qMerged <- paste0(mapping.df$qid,mapping.df$qDir);
    mapping.df$tMerged <- paste0(mapping.df$tid,mapping.df$tDir);
    ## remove useless / unlinkable mappings (start-start, or end-end)
    mapping.df <- subset(mapping.df, !((qStart == 1 & tStart == 1) | (qEnd == qLen & tEnd == tLen)));
    ## collect mappable contig names
    mapping.tigs <- unique(c(mapping.df$qid, mapping.df$tid));
    ## add in inverse links
    mapping.df.rev <- mapping.df;
    mapping.df.rev$qDir <- c("-"="+", "+"="-")[mapping.df$qDir];
    mapping.df.rev$tDir <- c("-"="+", "+"="-")[mapping.df$tDir];
    mapping.df.rev$qStart <- mapping.df$qLen - mapping.df$qEnd + 1;
    mapping.df.rev$qEnd <- mapping.df$qLen - mapping.df$qStart + 1;
    mapping.df.rev$tStart <- mapping.df$tLen - mapping.df$tEnd + 1;
    mapping.df.rev$tEnd <- mapping.df$tLen - mapping.df$tStart + 1;
    ## regenerate merged names
    mapping.df.rev$qMerged <- paste0(mapping.df.rev$qid,mapping.df.rev$qDir);
    mapping.df.rev$tMerged <- paste0(mapping.df.rev$tid,mapping.df.rev$tDir);
    mapping.df <- unique(rbind(mapping.df, mapping.df.rev));
    cluster.links.df <- subset(data.link.df, ((from %in% clusterComps) | (to %in% clusterComps)) &
                                             (from %in% mapping.tigs) & (to %in% mapping.tigs));
    cluster.links.df.rev <- cluster.links.df;
    cluster.links.df.rev$from <- cluster.links.df$to;
    cluster.links.df.rev$to <- cluster.links.df$from;
    cluster.links.df.rev$fromDir <- c("-"="+", "+"="-")[cluster.links.df$toDir];
    cluster.links.df.rev$toDir <- c("-"="+", "+"="-")[cluster.links.df$fromDir];
    cluster.links.df <- rbind(cluster.links.df, cluster.links.df.rev);
    cluster.links.df$fromMerged <- paste0(cluster.links.df$from,cluster.links.df$fromDir);
    cluster.links.df$toMerged <- paste0(cluster.links.df$to,cluster.links.df$toDir);
    cluster.graph <- graph.data.frame(cluster.links.df[,c("fromMerged","toMerged")]);
    cluster.paths <- list();
    seen.vertexes <- NULL;
    seen.vertexes[V(cluster.graph)] <- FALSE;
    for(Vin in V(cluster.graph)){
        subClust.paths <- all_simple_paths(cluster.graph, Vin);
        if(!all(seen.vertexes[unlist(subClust.paths)])){
            seen.vertexes[unlist(subClust.paths)] <- TRUE;
            cluster.paths <- c(cluster.paths, subClust.paths);
        }
    }
    ## remove sub-paths that are entirely represented by another path
    pathSeqs <- sapply(cluster.paths, paste, collapse=";");
    pathSubSeqs <- sapply(pathSeqs, grepl, x=pathSeqs, fixed=TRUE);
    pathSubSeqs[upper.tri(pathSubSeqs, diag=TRUE)] <- FALSE;
    if(length(pathSubSeqs) > 1){
        cluster.paths <- cluster.paths[which(colSums(pathSubSeqs) == 0)];
    }
    ## generate path
    added.seqs <- NULL;
    for(clPath in cluster.paths){
        pNames <- names(clPath);
        pNames.basic <- sub("[-\\+]$","",names(clPath));
        if(paste(pNames.basic, collapse="_") %in% added.seqs){
            next;
        }
        pNames.dir <- sub("^.*([-\\+])$","\\1",names(clPath));
        path.sequences <- gfa.sequences[pNames.basic];
        path.sequences[pNames.dir == "-"] <- reverseComplement(path.sequences[pNames.dir == "-"]);
        path.linkNames <- paste(head(pNames,-1),tail(pNames,-1),sep=" -> ");
        mapping.linkNames <- paste(mapping.df$qMerged,mapping.df$tMerged,sep=" -> ");
        found.rle <- rle(path.linkNames %in% mapping.linkNames);
        found.start <- cumsum(c(1,head(found.rle$lengths,-1)))[found.rle$values];
        found.end <- cumsum(found.rle$lengths)[found.rle$values];
        for(subPath in seq_along(found.start)){
            subpath.linkNames <- path.linkNames[found.start[subPath]:found.end[subPath]];
            subpath.mapping <- mapping.df[match(subpath.linkNames, mapping.linkNames),];
            cat("found mapping for:",subpath.linkNames,"\n");
            sp.name <- paste(pNames.basic[found.start[subPath]:(found.end[subPath]+1)],collapse="_");
            sm <- subpath.mapping[1,];
            current.sequences <- as.character(path.sequences[[sm$qid]]);
            for(seqPos in seq_along(subpath.linkNames)){
                sm <- subpath.mapping[seqPos,];
                current.sequences <- c(current.sequences,as.character(path.sequences[[sm$tid]][(sm$tEnd+1):sm$tLen]));
            }
            cat(sprintf(">%s\n", sp.name),
                sprintf("%s\n", paste(current.sequences, collapse="")), sep="",
                file=sprintf("%s.fasta",sp.name));
        }
        added.seqs <- c(added.seqs, paste(pNames.basic, collapse="_"),
                        paste(rev(pNames.basic), collapse="_"));
    }
    setwd(oldwd);
    if(any(cluster.links.df$fromDir != cluster.links.df$toDir)){
        break;
    }
}
