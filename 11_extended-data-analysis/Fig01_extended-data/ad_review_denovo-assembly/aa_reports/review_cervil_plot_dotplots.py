#!/usr/bin/env python3
from pathlib import Path
import matplotlib.pyplot as plt

PROJECT = Path("/work/project/ijpb/vibanez/14_cervil-project/07_cervil-denovo-assembly")
DOTDIR = PROJECT / "06_report"
OUTDIR = PROJECT / "06_report"
OUTDIR.mkdir(parents=True, exist_ok=True)

for paf in sorted(DOTDIR.glob("*.vs_SL5.paf")):
    xs = []
    ys = []

    with open(paf) as f:
        for line in f:
            parts = line.rstrip().split("\t")
            if len(parts) < 12:
                continue

            qstart = int(parts[2])
            qend   = int(parts[3])
            tstart = int(parts[7])
            tend   = int(parts[8])

            xs.extend([tstart, tend, None])
            ys.extend([qstart, qend, None])

    plt.figure(figsize=(8, 8))
    plt.plot(xs, ys, linewidth=0.2)
    plt.xlabel("SL5 coordinate")
    plt.ylabel("Cervil assembly coordinate")
    plt.title(paf.stem.replace(".vs_SL5", " vs SL5"))
    plt.tight_layout()

    out = OUTDIR / f"{paf.stem}.pdf"
    plt.savefig(out)
    plt.close()
    print(f"Wrote: {out}")
