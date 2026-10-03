#!/usr/bin/env bash
set -euo pipefail

SV_VCF="/mnt/disk2/vibanez/otherAnalysis/13_DMR-GWAS-SV/aa_markers/graphpanSV.sharedSamples.vcf.gz"
TFAM="/mnt/disk2/vibanez/10_data-analysis/Fig5/aa_GWAS-metabolome/ba_markers/DMR_general_leaf_LD.tfam"

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 svID [svID2 ...]" >&2
  exit 1
fi

tmpdir="$(mktemp -d)"
trap 'rm -rf "$tmpdir"' EXIT

awk '{print $2}' "$TFAM" > "$tmpdir/tfam.samples"

for sv_id in "$@"; do
  echo "=============================="
  echo "SV: $sv_id"

  bcftools view -i "ID=\"$sv_id\"" "$SV_VCF" | awk '!/^##/' > "$tmpdir/one_sv.tsv"

  if [[ ! -s "$tmpdir/one_sv.tsv" ]] || [[ $(wc -l < "$tmpdir/one_sv.tsv") -lt 2 ]]; then
    echo "Not found in VCF"
    continue
  fi

  awk '
  BEGIN{FS=OFS="\t"}
  FNR==NR{
    tfam[$1]=1
    next
  }
  NR==1{
    for(i=10;i<=NF;i++) samp[i]=$i
    next
  }
  NR==2{
    n_all=0; n_shared=0; n_nonmiss=0; n_miss=0
    for(i=10;i<=NF;i++){
      s=samp[i]
      n_all++
      if(!(s in tfam)) continue
      n_shared++

      split($i,a,":")
      g=a[1]

      if(g=="" || g ~ /\./){
        dose=-9
        n_miss++
      } else {
        m=split(g,al,/[\/|]/)
        dose=0
        bad=0
        for(j=1;j<=m;j++){
          if(al[j] !~ /^[0-9]+$/){ bad=1; break }
          if((al[j]+0)!=0) dose++
        }
        if(bad){
          dose=-9
          n_miss++
        } else {
          n_nonmiss++
        }
      }

      print s, g, dose
    }

    print "---- SUMMARY ----" > "/dev/stderr"
    print "all_vcf_samples\t" n_all > "/dev/stderr"
    print "shared_with_tfam\t" n_shared > "/dev/stderr"
    print "nonmissing_gt\t" n_nonmiss > "/dev/stderr"
    print "missing_gt\t" n_miss > "/dev/stderr"
  }' "$tmpdir/tfam.samples" "$tmpdir/one_sv.tsv" | sort
done
