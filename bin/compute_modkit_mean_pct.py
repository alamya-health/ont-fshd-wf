#!/usr/bin/env python3

import gzip
import sys


def main() -> int:
    bed_path = sys.argv[1]
    weighted = 0.0
    coverage = 0.0
    with gzip.open(bed_path, "rt") as fh:
        for line in fh:
            if not line.strip() or line.startswith("track"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 11:
                continue
            try:
                cov = float(parts[9])
                pct = float(parts[10])
            except ValueError:
                continue
            weighted += cov * pct
            coverage += cov
    print(f"{(weighted / coverage) if coverage else 0.0:.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
