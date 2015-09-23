#!/bin/sh

readType=$1; # IonTorrent, Nanopore, Illumina
bamFile=$2;
regFile=$3; # /data/all/david/diagnostic/design/IAD31669_Submitted.bed
refFile=$(samtools view -H ${bamFile} | grep '\.fa' | perl -pe 's#^.*-f #/mnt/ihbi_ngs/iontorrent-kgq4#;s/(\.fa(sta)?).*$/$1/');

echo "BAM File: '${bamFile}'" > /dev/stderr;
echo "Reference File: '${refFile}'" > /dev/stderr;
echo "Region File: '${regFile}'" > /dev/stderr;

pileupOpts="";

if [ ${readType} = "Nanopore" ] ; then
    # initial assumption: similar to IonTorrent, with lower base quality and lower INDEL coefficient
    pileupOpts="-d 10000 -L 10000 -Q 5 -h 20 -o 10 -e 17 -m 10 -f ${refFile} ${bamFile}";
fi

if [ ${readType} = "Illumina" ] ; then
    # Illumina should work with default options
    pileupOpts="-f ${refFile} ${bamFile}";
fi

# https://www.edgebio.com/variant-calling-ion-torrent-data

if [ ${readType} = "IonTorrent" ] ; then
    pileupOpts="-d 10000 -L 10000 -Q 7 -h 50 -o 10 -e 17 -m 10 -f ${refFile} ${bamFile}";
fi

# for Variant annotation with annovar
# http://www.cureffi.org/2012/09/07/an-alternative-exome-sequencing-pipeline-using-bowtie2-and-samtools/

if [ -z ${regFile} ]; then
    #echo "samtools mpileup ${pileupOpts} | ~/scripts/mpileup2Proportion.pl -m 10 | perl -pe 's/ +/,/g'" > /dev/stderr
    samtools mpileup ${pileupOpts} | ~/scripts/mpileup2Proportion.pl -m 10 | perl -pe 's/ +/,/g'
else
    #echo "samtools mpileup -l ${regFile} ${pileupOpts} | ~/scripts/mpileup2Proportion.pl -m 10 | perl -pe 's/ +/,/g'" > /dev/stderr
    samtools mpileup -l ${regFile} ${pileupOpts} | ~/scripts/mpileup2Proportion.pl -m 10 | perl -pe 's/ +/,/g'
fi

# for x in IonXpress_0*.bam; do echo -n ${x} "..."; ~/scripts/bam2proportion.sh IonTorrent ${x} /data/all/david/diagnostic/design/IAD31669_Submitted.bed > proportion_Diagnostics_AmpliSeq_QUT_NGS11_98_099/$(basename ${x} .bam).pileupprop.csv; echo "done"; done
# for x in IonXpress_0*.bam; do echo -n ${x} "..."; ~/scripts/bam2proportion.sh IonTorrent ${x} > proportion_Shani_Mito_Run18_71_082/$(basename ${x} .bam).pileupprop.csv; echo "done"; done

exit 0
