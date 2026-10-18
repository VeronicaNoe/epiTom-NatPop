#!/usr/bin/env bash
#SBATCH --job-name=cervil-syri
#SBATCH --cpus-per-task=24
#SBATCH --mem=128G
#SBATCH --time=2-00:00:00
#SBATCH --output=/work/project/ijpb/vibanez/14_cervil-project/07_cervil-denovo-assembly/06_report/logs/syri_%j.out
#SBATCH --error=/work/project/ijpb/vibanez/14_cervil-project/07_cervil-denovo-assembly/06_report/logs/syri_%j.err

set -euo pipefail

THREADS="${SLURM_CPUS_PER_TASK:-24}"

PROJECT="/work/project/ijpb/vibanez/14_cervil-project/07_cervil-denovo-assembly"
REPORT="${PROJECT}/06_report"
OUTDIR="${REPORT}/syri_plotsr_chr12"

REF0="/home/vibanez/work/01_mapping/01.0_ref/SL5.0.fasta"
QRY0="${PROJECT}/03_ragtag/ragtag.scaffold.fasta"

REF="${OUTDIR}/SL5.chr12.renamed.fasta"
QRY="${OUTDIR}/cervil.ragtag.chr12.renamed.fasta"

mkdir -p "${OUTDIR}" "${REPORT}/logs"
cd "${OUTDIR}"

source /etc/profile.d/modules.sh 2>/dev/null || true
module purge
module load devel/Miniconda/Miniconda3
module load bioinfo/Minimap/2-2.28
module load bioinfo/samtools/1.23
module load bioinfo/SyRI/1.6.3

echo "[check tools]" >&2
which minimap2 >&2
which samtools >&2
which syri >&2

echo "[1/5] prepare SL5 chr1-12 renamed as Chr1-Chr12" >&2
printf "1\n2\n3\n4\n5\n6\n7\n8\n9\n10\n11\n12\n" > sl5_chr12.list

samtools faidx "${REF0}" $(cat sl5_chr12.list) \
  | awk '
      /^>/ {
        gsub(/^>/, "", $0);
        print ">Chr"$0;
        next
      }
      {print}
    ' > "${REF}"

samtools faidx "${REF}"

echo "[2/5] prepare Cervil RagTag chr1-12 renamed as Chr1-Chr12" >&2
printf "1_RagTag\n2_RagTag\n3_RagTag\n4_RagTag\n5_RagTag\n6_RagTag\n7_RagTag\n8_RagTag\n9_RagTag\n10_RagTag\n11_RagTag\n12_RagTag\n" > ragtag_chr12.list

samtools faidx "${QRY0}" $(cat ragtag_chr12.list) \
  | awk '
      /^>/ {
        gsub(/^>/, "", $0);
        gsub(/_RagTag$/, "", $0);
        print ">Chr"$0;
        next
      }
      {print}
    ' > "${QRY}"

samtools faidx "${QRY}"

echo "[3/5] minimap2 alignment" >&2
minimap2 -ax asm5 --eqx -t "${THREADS}" "${REF}" "${QRY}" \
  | samtools sort -@ "${THREADS}" -o cervil_vs_SL5.sorted.bam -

samtools index -@ "${THREADS}" cervil_vs_SL5.sorted.bam

echo "[4/5] SyRI" >&2
syri \
  -c cervil_vs_SL5.sorted.bam \
  -r "${REF}" \
  -q "${QRY}" \
  -F B \
  --prefix cervil_vs_SL5. \
  > syri.stdout.log \
  2> syri.stderr.log

echo "[5/5] genomes file for plotsr" >&2
cat > genomes.txt <<EOF
#file	name	tags
${REF}	SL5	lw:1.5
${QRY}	Cervil	lw:1.5
EOF

echo "[done] SyRI finished" >&2
ls -lh *syri* genomes.txt >&2
