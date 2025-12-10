# for CG
cat $1 | awk '{OFS="\t"}{print $1,$2,$3,$9/$8}' |\
sed 's/,/./g' > /mnt/disk2/vibanez/02_methylkit/ae_cumulative-DMR/aa_data/$1
# C
cat $1 | awk '{OFS="\t"}{print $1,$2,$3,$6/$5}' |\
sed 's/,/./g' > /mnt/disk2/vibanez/02_methylkit/ae_cumulative-DMR/aa_data/$1
