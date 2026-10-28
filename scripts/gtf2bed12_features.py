from collections import defaultdict
import re

feature = snakemake.params.feature
attr_re = re.compile(r'(\w+) "([^"]*)"')

groups = defaultdict(list)  # id -> [(chrom, start0, end, strand)]
with open(snakemake.input[0]) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != feature:
            continue
        attrs = dict(attr_re.findall(f[8]))
        # NCBI prokaryote GTFs often leave transcript_id empty
        gid = attrs.get("transcript_id") or attrs.get("gene_id")
        if gid:
            groups[gid].append((f[0], int(f[3]) - 1, int(f[4]), f[6]))

rows = []
for gid, blocks in groups.items():
    if len({b[0] for b in blocks}) > 1:
        continue
    chrom, strand = blocks[0][0], blocks[0][3]
    blocks = sorted((b[1], b[2]) for b in blocks)
    start, end = blocks[0][0], max(b[1] for b in blocks)
    sizes = ",".join(str(e - s) for s, e in blocks) + ","
    starts = ",".join(str(s - start) for s, _ in blocks) + ","
    rows.append((chrom, start, end, gid, "0", strand, start, end, "0", len(blocks), sizes, starts))

rows.sort(key=lambda r: (r[0], r[1]))
with open(snakemake.output.bed, "w") as out:
    for r in rows:
        out.write("\t".join(map(str, r)) + "\n")
