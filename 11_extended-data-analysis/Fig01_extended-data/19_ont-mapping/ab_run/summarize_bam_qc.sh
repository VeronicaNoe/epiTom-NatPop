#!/usr/bin/env bash
set -euo pipefail

BASE="/mnt/disk2/vibanez/otherAnalysis/19_pedigree-denovo-assembly"
BAMDIR="/mnt/ssd123/vibanez/19_ont-mapping/ab_run/merged_by_sample/mods_bam"
OUTDIR="${BASE}/methylation_qc/bam_qc"

JOBS="${JOBS:-3}"
THREADS_PER_JOB="${THREADS_PER_JOB:-4}"

mkdir -p "${OUTDIR}"

summarize_bam_one() {
    bam="$1"
    base=$(basename "$bam")
    sample="${base%.mods.merged.sorted.bam}"
    generation="${sample%%-*}"

    covout="${OUTDIR}/${sample}.samtools_coverage.tsv"
    flagout="${OUTDIR}/${sample}.flagstat.txt"

    samtools flagstat -@ "${THREADS_PER_JOB}" "$bam" > "$flagout"
    samtools coverage "$bam" > "$covout"

    awk -v sample="$sample" \
        -v generation="$generation" \
        'BEGIN { FS=OFS="\t" }
        NR > 1 {
            len = $3 - $2 + 1
            total_len += len
            total_covbases += $5
            weighted_depth += $7 * len
            total_reads += $4
        }
        END {
            genome_coverage_pct = (total_len > 0 ? 100 * total_covbases / total_len : "NA")
            mean_depth = (total_len > 0 ? weighted_depth / total_len : "NA")
            print sample, generation, total_len, total_covbases, genome_coverage_pct, mean_depth, total_reads
        }' "$covout"
}

export -f summarize_bam_one
export OUTDIR THREADS_PER_JOB

{
    echo -e "sample\tgeneration\tgenome_size_bp\tcovered_bases\tgenome_coverage_pct\tmean_depth\ttotal_reads_contig_sum"

    find "${BAMDIR}" -name "*.mods.merged.sorted.bam" -print0 \
        | xargs -0 -n 1 -P "${JOBS}" bash -c 'summarize_bam_one "$0"'
} > "${OUTDIR}/bam_coverage_summary.tsv"

echo "Done:"
echo "${OUTDIR}/bam_coverage_summary.tsv"
