#!/usr/bin/env python3

import collections
import os
import sys

import pysam


def main() -> int:
    bam_path = sys.argv[1]
    classification_dir = sys.argv[2]
    out_path = sys.argv[3]

    subset_map = {
        "4qA_all": "4qA_all-reads-ID.txt",
        "4qA_complete": "4qA_complete-reads-ID.txt",
        "4qB_all": "4qB_all-reads-ID.txt",
        "4qB_complete": "4qB_complete-reads-ID.txt",
        "chimeric": "chimeric-reads-ID.txt",
        "chr10_all": "chr10_all-reads-ID.txt",
        "chr10_complete": "chr10_complete-reads-ID.txt",
    }

    ids_by_subset = {}
    union_ids = set()
    for subset_name, filename in subset_map.items():
        path = os.path.join(classification_dir, filename)
        ids = set()
        if os.path.exists(path):
            with open(path, "r", encoding="utf-8") as fh:
                ids = {line.strip() for line in fh if line.strip()}
        ids_by_subset[subset_name] = ids
        union_ids.update(ids)

    hp_by_read = {}
    has_any_hp = False
    if union_ids:
        bam = pysam.AlignmentFile(bam_path, "rb", check_sq=False)
        for aln in bam.fetch(until_eof=True):
            rid = aln.query_name
            if rid not in union_ids:
                continue
            hp = None
            if aln.has_tag("HP"):
                try:
                    hp = str(aln.get_tag("HP"))
                    has_any_hp = True
                except Exception:
                    hp = None
            is_primary = (not aln.is_secondary) and (not aln.is_supplementary)
            prev = hp_by_read.get(rid)
            if prev is None:
                hp_by_read[rid] = (is_primary, hp if hp is not None else "UNSET")
            elif is_primary and not prev[0]:
                hp_by_read[rid] = (True, hp if hp is not None else "UNSET")
        bam.close()

    with open(out_path, "w", encoding="utf-8") as out:
        out.write("subset_name\thp\tcount\ttotal_subset_reads\thas_hp_tags\n")
        for subset_name in sorted(subset_map):
            ids = ids_by_subset.get(subset_name, set())
            counts = collections.Counter()
            for rid in ids:
                hp = hp_by_read.get(rid, (False, "UNSET"))[1]
                counts[hp] += 1
            if not counts:
                counts["UNSET"] = 0
            total = len(ids)
            for hp, count in sorted(counts.items(), key=lambda x: x[0]):
                out.write(
                    f"{subset_name}\t{hp}\t{count}\t{total}\t{str(has_any_hp).lower()}\n"
                )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
