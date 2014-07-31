#!/bin/sh
#
# Read a GAS or LINKAGE format pedigree, return a digraph in the dot language 
# call dot to make pedigree drawing
# 
AWK=$(which awk)
DOTEXE=$(which dot)
NEATOEXE=$(which neato)

for fil in $*
do
  for ped in `$AWK '!/^[!#]/ {print $1}' $fil | sort -u`
  do
     echo "Pedigree $ped"
     $AWK -v ped=$ped '
     BEGIN { shape["f"]="box,regular=1"
             shape["1"]="box,regular=1"
             shape["m"]="circle"
             shape["2"]="circle"
             shape["u"]="diamond"
             shape["0"]="diamond"
             shade["y"]="grey"
             shade["2"]="grey"
             shade["n"]="white"
             shade["1"]="white"
             shade["x"]="white"
             shade["0"]="white"
     }
     !/^[!#]/ && $1==ped {
             sex[$2]=$5
             aff[$2]="x" ; if ($6 ~ /[012nyx]/) aff[$2]=$6
             if($3!="x" && $3!="0") {
               marriage[$3,$4]++
               child[$3,$4,marriage[$3,$4]]=$2
             }
     }
     END   { print "digraph Ped_" ped " {"
#             print "# page =\"11,8.5\" ;"
             print "node [shape=diamond] ;"
             print "ratio =\"auto\" ;"
             print "mincross = 2.0 ;"
             print "label=\"Pedigree " ped "\" ;"
             print "rotate=0 ;"
             for(s in sex) {
               print "\"" s "\" [shape=" shape[sex[s]] ","  \
                     " style=filled,fillcolor=" shade[aff[s]] "] ;"
             }
             for(m in marriage) {
               n=split(m,par,"\034")
               mating_t="\"t_" par[1] "x" par[2] "\""
               mating_b="\"b_" par[1] "x" par[2] "\""
               print mating_t "[shape=diamond,style=filled," \
                     "label=\"\",height=.1,width=.1] ;"
               print mating_b "[shape=diamond,style=filled," \
                     "label=\"\",height=.1,width=.1] ;"
               print "\"" par[1] "\" -> " mating_t " [dir=none, weight=1, penwidth=3.0] ;"
               print "\"" par[2] "\" -> " mating_t " [dir=none, weight=1, penwidth=3.0] ;"
               print mating_t " -> " mating_b " [dir=none, weight=1, penwidth=3.0] ;"
               for(k=1;k<=marriage[par[1],par[2]];k++) {
                 print  mating_b " -> \"" child[par[1],par[2],k] "\"" \
                        " [dir=none, weight=2] ;"
               }
             }
             print "}"
     }' $fil > $ped.dot
#     echo "running ${DOTEXE} -Tsvg ${ped}.dot -o ${ped}_dot.svg"
     ${DOTEXE} -Tpdf ${ped}.dot -o ${ped}_dot.pdf
#     ${NEATOEXE} -Tsvg ${ped}.dot -o ${ped}_neato.svg
  done
done
