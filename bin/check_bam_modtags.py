#!/usr/bin/env python3

import sys

import pysam


def main() -> int:
    bam_path = sys.argv[1]
    bam = pysam.AlignmentFile(bam_path, "rb", check_sq=False)
    found = False
    for idx, aln in enumerate(bam.fetch(until_eof=True)):
        if aln.has_tag("MM") or aln.has_tag("Mm"):
            found = True
            break
        if idx >= 2000:
            break
    bam.close()
    return 0 if found else 1


if __name__ == "__main__":
    raise SystemExit(main())
