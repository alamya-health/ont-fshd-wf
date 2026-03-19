#!/usr/bin/env python3

import gzip
import re
import sys


def count_vcf_records(path: str) -> int:
    opener = gzip.open if path.endswith(".gz") else open
    count = 0
    with opener(path, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            count += 1
    return count


def parse_relevant(path: str):
    genes = set()
    clinvar_hits = 0
    records = 0
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            records += 1
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue
            info = cols[7]
            if "CLNDN=" in info and "Facioscapulohumeral".lower() in info.lower():
                clinvar_hits += 1
            match = re.search(r"ANN=([^;]+)", info)
            if match:
                for ann in match.group(1).split(","):
                    parts = ann.split("|")
                    if len(parts) > 3 and parts[3]:
                        genes.add(parts[3])
    return records, clinvar_hits, ",".join(sorted(genes)) if genes else "NA"


def main() -> int:
    _, clair3_vcf_gz, relevant_vcf_gz, hg38_sv_vcf, t2t_sv_vcf, summary_out = sys.argv[1:]

    total_clair3 = count_vcf_records(clair3_vcf_gz)
    relevant_count, clinvar_hits, genes = parse_relevant(relevant_vcf_gz)
    hg38_sv_count = count_vcf_records(hg38_sv_vcf)
    t2t_sv_count = count_vcf_records(t2t_sv_vcf)

    with open(summary_out, "w", encoding="utf-8") as fh:
        fh.write(
            "status\ttotal_clair3_variants\trelevant_variant_records\tclinvar_fshd_hits\thg38_sv_records\tt2t_sv_records\trelevant_genes\n"
        )
        fh.write(
            f"ready\t{total_clair3}\t{relevant_count}\t{clinvar_hits}\t{hg38_sv_count}\t{t2t_sv_count}\t{genes}\n"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
