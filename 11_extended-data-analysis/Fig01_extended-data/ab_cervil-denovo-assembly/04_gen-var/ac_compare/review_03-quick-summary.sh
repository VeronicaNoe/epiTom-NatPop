echo "================ SV summary ================================="
awk -F'\t' 'NR>1{c[$8]++; n++} END{for(k in c) printf "%s\t%d\t%.4f\n", k, c[k], c[k]/n}' sv.classified.tsv | sort
echo "================ Small variants summary ====================="
awk -F'\t' 'NR>1{c[$10]++; n++} END{for(k in c) printf "%s\t%d\t%.4f\n", k, c[k], c[k]/n}' smallvars.classified.tsv | sort
