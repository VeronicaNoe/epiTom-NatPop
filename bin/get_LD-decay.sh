#!/bin/bash
GROUP="$( cat $1 | cut -d'-' -f1,2)"
#!/bin/bash

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

plink1.9 --vcf "$INDIR/${SAMPLE}.vcf.gz" \
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

python ~/bin/ld_decay_calc.py -i $OUTPUT/$NAME".ld.gz" -o $OUTPUT/$NAME
