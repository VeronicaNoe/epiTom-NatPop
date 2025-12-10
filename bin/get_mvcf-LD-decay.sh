#!/bin/bash
#plink1.9
MM="$( cat $1 | cut -d'_' -f1 )"
ANNO="$( cat $1 | cut -d'_' -f2 )"
GROUP="$( cat $1 | cut -d'_' -f3 )"
INDIR="/mnt/disk2/vibanez/06_get-meth-vcf/ab_epialleles"

if [[ "$GROUP" == "all" ]]; then
    OUTDIR="/mnt/disk2/vibanez/10_data-analysis/Fig2/ab_get-LD/bb_general"
    if [[ "$ANNO" == "general" ]]; then
        SAMPLE="${MM}"
    else
        SAMPLE="${MM}_${ANNO}"
    fi

else
    OUTDIR="/mnt/disk2/vibanez/10_data-analysis/Fig2/ab_get-LD/bc_by-groups"
    SAMPLE="${MM}_${ANNO}_${GROUP}"
fi

plink --vcf "$INDIR/${SAMPLE}.vcf.gz" \
         --threads 40 \
         --double-id \
         --allow-extra-chr \
         --maf 0.05 \
         --make-bed \
         --geno 0.2 \
         --mind 0.2 \
         -r2 gz \
         --ld-window-kb 2000 \
      	--ld-window-r2 0 \
         --out "${OUTDIR}/${SAMPLE}_LD"

#~/bin/ld_decay_calc.py -i $OUTDIR/$SAMPLE"_LD.ld.gz" -o $OUTDIR/$SAMPLE
python3 ~/bin/ld_decay_calc_py3.py --ld $OUTDIR/$SAMPLE"_LD.ld.gz" \
	--out ${OUTDIR}/${SAMPLE}"_LD-decay.tsv"
