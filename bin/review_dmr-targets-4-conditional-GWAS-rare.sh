#!/usr/bin/env bash
SAMPLE="$(cat "$1")"
SIG_GWAS="/mnt/disk2/vibanez/10_data-analysis/Fig5/aa_GWAS-metabolome/bd_results/sig"
WDIR="/mnt/disk2/vibanez/otherAnalysis/09_conditional-metabolite-GWAS-rare-variants"
TDIR="$WDIR/aa_get-target-dmrs"
mkdir -p "$TDIR"

PTHRESH="1.048466e-07"
TOP_DMR_N=1

chr_to_num() {
  local c="$1"
  c="${c#chr}"
  c="${c#ch}"
  c="$(echo "$c" | sed 's/^0\+\([0-9]\)/\1/; s/^0$/0/')"
  echo "$c"
}

baseline_stats() {
  local ps_gz="$1" dmr="$2"
  zcat "$ps_gz" | awk -v d="$dmr" 'BEGIN{FS="\t"} $1==d {print $2"\t"$3"\t"$4; exit}'
}

TAG="${SAMPLE%.ps.gz}"
F="$SIG_GWAS/${TAG}.ps.gz"
bn="$TAG"
meta="${bn%%.leaf-metabolites_*}"
rest="${bn#*.leaf-metabolites_}"

IFS="_" read -r -a parts <<< "$rest"
mm="${parts[0]}"
index="${parts[-1]}"
if (( ${#parts[@]} > 2 )); then
  kinship_label="$(IFS="_"; echo "${parts[*]:1:${#parts[@]}-2}")"
else
  kinship_label="${parts[1]}"
fi

TARGETS="$TDIR/${SAMPLE}.targets.tsv"
echo -e "meta\tmm\tkinship_label\tindex\tbaseline_ps_gz\tlead_dmr\tlead_chr_raw\tlead_chr_num\tlead_pos\tbeta_base\tbeta_sd_base\tp_base" > "$TARGETS"

leads="$(zcat "$F" \
  | gawk -v pth="$PTHRESH" 'BEGIN{FS="[ \t]+"} ($4+0)<=pth {print $0}' \
  | sort -gk4,4 | head -n "$TOP_DMR_N" \
  | awk 'BEGIN{FS="[ \t]+"}{print $1}' )"

if [[ -z "$leads" ]]; then
  leads="$(zcat "$F" | sort -gk4,4 | head -n "$TOP_DMR_N" | awk 'BEGIN{FS="[ \t]+"}{print $1}')"
fi

while read -r dmr; do
  [[ -z "$dmr" ]] && continue
  chr_raw="${dmr%%:*}"
  pos_raw="$(echo "$dmr" | awk -F: '{print $2}')"
  chr_num="$(chr_to_num "$chr_raw")"

  stats="$(baseline_stats "$F" "$dmr" || true)"
  if [[ -z "$stats" ]]; then
    beta_base="NA"; beta_sd_base="NA"; p_base="NA"
  else
    beta_base="$(echo "$stats" | cut -f1)"
    beta_sd_base="$(echo "$stats" | cut -f2)"
    p_base="$(echo "$stats" | cut -f3)"
  fi

  echo -e "${meta}\t${mm}\t${kinship_label}\t${index}\t${F}\t${dmr}\t${chr_raw}\t${chr_num}\t${pos_raw}\t${beta_base}\t${beta_sd_base}\t${p_base}" >> "$TARGETS"
done <<< "$leads"
