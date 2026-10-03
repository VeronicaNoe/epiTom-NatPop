#!/bin/bash
SAMPLE="$(cat "$1")"
CONTIG_LIST="/mnt/ssd123/vibanez/19_ont-mapping/ab_run/polished_contigs_used_in_ragtag.ids"
OUT_DIR="/mnt/disk2/vibanez/otherAnalysis/19_pedigree-denovo-assembly/af_chr-split"
IN_BED="${OUT_DIR}/${SAMPLE}.prefiltered.bed"

awk -v outdir="${OUT_DIR}" -v sample="${SAMPLE}" '
BEGIN { FS=OFS="\t" }
NR==FNR {
    keep[$1]=1
    next
}
$1 in keep {
    out = outdir "/" sample "_" $1 ".bed"
    print > out
}
' "${CONTIG_LIST}" "${IN_BED}"
