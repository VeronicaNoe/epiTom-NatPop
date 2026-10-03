LIST="/mnt/disk2/vibanez/otherAnalysis/13_DMR-GWAS-SV/af_get-nonSig-GIF/00.0_list2move-from-sig2nonsig.tsv"
BASE="/mnt/disk2/vibanez/otherAnalysis/13_DMR-GWAS-SV/ad_results"

tail -n +2 "$LIST" | while IFS= read -r id; do
    chr="${id%%_*}"
    rest="${id#*_}"
    dmr="${rest%%_*}"

    src="$BASE/$dmr/$chr/sig"
    dst="$BASE/$dmr/$chr/nonSig-GIF"

    mkdir -p "$dst"
    [[ -d "$src" ]] || continue

    find "$src" -maxdepth 1 \( -type f -o -type l \) -name "${id}.*" -exec mv -n -t "$dst" {} +
done
