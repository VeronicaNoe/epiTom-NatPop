import sys

mode = sys.argv[1]   # "small" or "sv"
infile = sys.argv[2]
outfile = sys.argv[3]

def carrier(gt: str) -> int:
    gt = gt.strip()
    if gt in {".", "./.", ".|.", "0/0", "0|0"}:
        return 0
    alleles = gt.replace("|", "/").split("/")
    return int(any(a not in {"0", "."} for a in alleles))

with open(infile) as f, open(outfile, "w") as out:
    if mode == "small":
        out.write("chrom\tstart\tend\tpos\tref\talt\tvariant_type\tn_g0\tn_g5\tclass\tG0-C3\tG0-C6\tG0-C9\tG5-C1\tG5-C2\tG5-C3\n")
        for line in f:
            chrom, pos, ref, alt, g0c3, g0c6, g0c9, g5c1, g5c2, g5c3 = line.rstrip("\n").split("\t")
            pos = int(pos)
            alts = alt.split(",")
            span = max([len(ref)] + [len(a) for a in alts])
            start = pos - 1
            end = start + span
            vtype = "SNP" if all(len(ref) == 1 and len(a) == 1 for a in alts) else "INDEL"
            gts = [g0c3, g0c6, g0c9, g5c1, g5c2, g5c3]

            c = [carrier(x) for x in gts]
            n_g0 = sum(c[:3])
            n_g5 = sum(c[3:])

            if n_g0 > 0 and n_g5 == 0:
                cls = "fixed_G0" if n_g0 == 3 else "seg_G0"
            elif n_g5 > 0 and n_g0 == 0:
                cls = "fixed_G5" if n_g5 == 3 else "seg_G5"
            elif n_g0 == 3 and n_g5 == 3:
                cls = "shared_all"
            elif n_g0 > 0 and n_g5 > 0:
                cls = "mixed_both"
            else:
                cls = "none"

            out.write(f"{chrom}\t{start}\t{end}\t{pos}\t{ref}\t{alt}\t{vtype}\t{n_g0}\t{n_g5}\t{cls}\t" + "\t".join(gts) + "\n")

    else:
        out.write("chrom\tstart\tend\tpos\tsvtype\tn_g0\tn_g5\tclass\tG0-C3\tG0-C6\tG0-C9\tG5-C1\tG5-C2\tG5-C3\n")
        for line in f:
            chrom, pos, end, svtype, g0c3, g0c6, g0c9, g5c1, g5c2, g5c3 = line.rstrip("\n").split("\t")
            pos = int(pos)
            end = int(end) if end != "." else pos
            start = pos - 1
            gts = [g0c3, g0c6, g0c9, g5c1, g5c2, g5c3]

            c = [carrier(x) for x in gts]
            n_g0 = sum(c[:3])
            n_g5 = sum(c[3:])

            if n_g0 > 0 and n_g5 == 0:
                cls = "fixed_G0" if n_g0 == 3 else "seg_G0"
            elif n_g5 > 0 and n_g0 == 0:
                cls = "fixed_G5" if n_g5 == 3 else "seg_G5"
            elif n_g0 == 3 and n_g5 == 3:
                cls = "shared_all"
            elif n_g0 > 0 and n_g5 > 0:
                cls = "mixed_both"
            else:
                cls = "none"

            out.write(f"{chrom}\t{start}\t{end}\t{pos}\t{svtype}\t{n_g0}\t{n_g5}\t{cls}\t" + "\t".join(gts) + "\n")
