# Cervil pooled reference assembly report

## Assembly workflow
- ONT pooled reads
- hifiasm primary contig assembly
- Dorado polishing
- RagTag scaffolding against SL5.0
- QC with QUAST and BUSCO
- Annotation with Liftoff and EDTA

## Main assembly files
- Raw contig assembly: `/work/project/ijpb/vibanez/14_cervil-project/07_cervil-denovo-assembly/01_hifiasm/cervil.primary.fasta`
- Polished assembly: `/work/project/ijpb/vibanez/14_cervil-project/07_cervil-denovo-assembly/02_polish/cervil.polished.fasta`
- RagTag assembly: `/work/project/ijpb/vibanez/14_cervil-project/07_cervil-denovo-assembly/03_ragtag/ragtag.scaffold.fasta`

## Key tables
- `06_report/assembly_summary.tsv`
- `06_report/liftoff_summary.tsv`
- `06_report/edta_summary.tsv`

## Key figures
- `06_report/hifiasm_raw.vs_SL5.pdf`
- `06_report/dorado_polish.vs_SL5.pdf`
- `06_report/ragtag.vs_SL5.pdf`

## Interpretation points to complete
- Compare contiguity across raw, polished, and ragtag
- Evaluate BUSCO completeness
- Evaluate whether RagTag improves chromosome-scale representation without obvious structural artifacts
- Summarize gene transfer success from Liftoff
- Summarize repeat annotation from EDTA
- Decide final working assembly for downstream analyses
