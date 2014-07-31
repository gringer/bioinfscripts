#!/bin/sh

if [ $# -lt 1 ]
    then echo "Error: no arguments have been specified"
    echo "usage <path1/sample1.fasta.gz> [path2/sample2.fasta] [-nomerge]"
    exit
fi

ffn1=$1
d1=$(dirname $(readlink -e ${ffn1}))
f1=$(basename ${ffn1})
y1=$(basename ${f1} .gz)
y1=$(basename ${y1} .gzip)
y1=$(basename ${y1} .fasta)
y1=$(basename ${y1} .fa)
mergeName=${y1}_mm2
refCount=$(zgrep -c '^>' ${d1}/${f1})
d2=""
f2=""
if [ $# -gt 1 ]
    then ffn2=$2
    d2=$(dirname $(readlink -e ${ffn2}))
    f2=$(basename ${ffn2})
    y2=$(basename ${f2} .gz)
    y2=$(basename ${y2} .gzip)
    y2=$(basename ${y2} .fasta)
    y2=$(basename ${y2} .fa)
    mergeName=${y1}_${y2}_mm2
    if [ $3 != "-nomerge" ]
        then echo "-- merging ${y2} into ${y1} --"
        else echo "-- carrying out all-vs-all overlap on ${y1} + ${y2} --"
    fi
    else echo "-- carrying out all-vs-all overlap on ${y1} --"
fi
mkdir -p ${d1}/${mergeName}
cd ${d1}/${mergeName}
zcat -f ${d1}/${f1} > ${mergeName}.seq
if [ $# -gt 1 ]
    then zcat -f ${d2}/${f2} >> ${mergeName}.seq
fi
~/bin/amos/toAmos -s ${mergeName}.seq -o ${mergeName}.afg
if [ $3 = "-nomerge" ]
    then ~/bin/amos/minimus2 ${mergeName}
    else ~/bin/amos/minimus2 ${mergeName} -D REFCOUNT=${refCount}
fi
echo "-- DONE --"
