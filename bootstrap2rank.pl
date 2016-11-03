#!/usr/bin/sh

pv bootstrap100_anzgene200_vs_1kg200.csv.gz | zcat | sort -t',' -k 2,2n -k 3,3rn | perl -F',' -lane 'if($run != $F[1]){$run = $F[1]; $rank=1; $lastVal=0; $lastRank=0} if($F[1] ne "bs.run"){$nextVal = $F[2]; $F[2] = ($F[2] == $lastVal) ? $lastRank : $rank; $lastRank = $F[2]; $lastVal = $nextVal; $rank++} print join(",",@F)' | sort -t',' -k 1,1 -k 3,3n | perl -F',' -lane 'if($marker ne $F[0]){if($marker){print $lastLine; print $_} $marker = $F[0]} {$lastLine = $_} END{ print $lastLine}' | gzip > maxminRank_bootstrap100_anzgene200_vs_1kg200.csv.gz
