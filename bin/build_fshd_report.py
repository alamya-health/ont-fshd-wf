#!/usr/bin/env python3

import csv
import gzip
import html
import os
import sys
from collections import Counter, defaultdict


REPEAT_BP = 3303


def clean_cell(value):
    if value is None:
        return ""
    if not isinstance(value, str):
        return value
    value = value.strip()
    if len(value) >= 2 and value[0] == '"' and value[-1] == '"':
        value = value[1:-1]
    return value.strip()


def read_semicolon_csv(path):
    rows = []
    if not path or not os.path.exists(path):
        return rows
    with open(path, "r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter=";")
        for row in reader:
            rows.append({clean_cell(k): clean_cell(v) for k, v in row.items()})
    return rows


def read_tsv(path):
    rows = []
    if not path or not os.path.exists(path):
        return rows
    with open(path, "r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            rows.append({clean_cell(k): clean_cell(v) for k, v in row.items()})
    return rows


def to_int(value, default=0):
    try:
        text = clean_cell(value)
        if text in ("", "NA", "nan", "None"):
            return default
        return int(float(text.replace(",", ".")))
    except Exception:
        return default


def to_float(value, default=0.0):
    try:
        text = clean_cell(value)
        if text in ("", "NA", "nan", "None"):
            return default
        return float(text.replace(",", "."))
    except Exception:
        return default


def truncate(text, width=80):
    text = clean_cell(text)
    return text if len(text) <= width else text[: width - 3] + "..."


def metric_card(label, value, note="", tone="default"):
    return (
        f"<div class='metric {tone}'>"
        f"<div class='metric-label'>{html.escape(label)}</div>"
        f"<div class='metric-value'>{html.escape(str(value))}</div>"
        f"<div class='metric-note'>{html.escape(note)}</div>"
        "</div>"
    )


def token_class(token):
    token = clean_cell(token)
    if token.startswith("c4"):
        return "repeat-short" if token.endswith("S") else ("repeat-long" if token.endswith("L") else "repeat")
    if token.startswith("c10"):
        return "offtarget"
    if "D4F104S1" in token or "CLUHP4" in token:
        return "anchor-prox"
    if "pLAM" in token:
        return "anchor-dist"
    if "DUX4_end" in token:
        return "anchor-dist"
    if token.startswith("end_"):
        return "soft-end"
    return "misc"


def display_token(token):
    token = clean_cell(token)
    if token == "CLUHP4-201_exon1":
        return "CLUHP4"
    if token == "4qA_D4F104S1":
        return "4qA D4F104S1"
    if token == "4qA_pLAM":
        return "4qA pLAM"
    if token == "DUX4_end":
        return "DUX4 end"
    return token.replace("_", " ")


def parse_repeat_tokens(seq):
    seq = clean_cell(seq)
    if not seq:
        return []
    return [token.strip() for token in seq.split(" - ") if token.strip()]


def has_proximal_anchor(tokens):
    return any("D4F104S1" in token or "CLUHP4" in token for token in tokens)


def has_distal_anchor(tokens):
    return any("pLAM" in token or "DUX4_end" in token for token in tokens)


def anchor_mode(tokens):
    prox = has_proximal_anchor(tokens)
    dist = has_distal_anchor(tokens)
    if prox and dist:
        return "both"
    if prox:
        return "proximal"
    if dist:
        return "distal"
    return "internal"


def contraction_assessment(complete_rows, threshold):
    complete_rus = [to_int(row.get("RU_count"), -1) for row in complete_rows if to_int(row.get("RU_count"), -1) > 0]
    if not complete_rows:
        return (
            "inconclusive",
            f"No complete 4qA-spanning read was classified, so a contraction at <= {threshold} repeat units cannot be called directly.",
        )
    if any(ru <= threshold for ru in complete_rus):
        support = sum(1 for ru in complete_rus if ru <= threshold)
        return (
            "supported",
            f"{support} complete 4qA read(s) measured at <= {threshold} repeat units.",
        )
    shortest = min(complete_rus) if complete_rus else "NA"
    return (
        "not_supported",
        f"Complete 4qA reads were present, but the shortest measured array was {shortest} repeat units.",
    )


def hist_svg(values, color, label):
    if not values:
        return "<p class='muted'>No repeat-unit measurements were available for this panel.</p>"
    counts = Counter(values)
    xs = sorted(counts)
    max_count = max(counts.values())
    bar_width = 44
    gap = 12
    width = max(460, len(xs) * (bar_width + gap) + 70)
    height = 230
    bars = []
    for idx, x in enumerate(xs):
        count = counts[x]
        bar_h = int((count / max_count) * 132)
        xpos = 44 + idx * (bar_width + gap)
        ypos = 176 - bar_h
        bars.append(
            f"<rect x='{xpos}' y='{ypos}' width='{bar_width}' height='{bar_h}' rx='10' fill='{color}' />"
        )
        bars.append(f"<text x='{xpos + bar_width / 2}' y='198' text-anchor='middle' class='axis'>{x}</text>")
        bars.append(
            f"<text x='{xpos + bar_width / 2}' y='{max(28, ypos - 8)}' text-anchor='middle' class='count'>{count}</text>"
        )
    return (
        f"<svg viewBox='0 0 {width} {height}' class='histogram' role='img' aria-label='{html.escape(label)}'>"
        f"<text x='28' y='26' class='title'>{html.escape(label)}</text>"
        f"<line x1='26' y1='176' x2='{width - 18}' y2='176' stroke='rgba(58,45,36,0.22)' stroke-width='2'/>"
        + "".join(bars)
        + "</svg>"
    )


def repeat_strip_html(tokens):
    if not tokens:
        return "<div class='repeat-strip muted'>No repeat architecture resolved.</div>"
    blocks = []
    for token in tokens:
        css = token_class(token)
        blocks.append(
            f"<span class='unit {css}' title='{html.escape(token)}'>{html.escape(display_token(token))}</span>"
        )
    return "<div class='repeat-strip'>" + "".join(blocks) + "</div>"


def table_html(headers, rows, klass=""):
    if not rows:
        return "<p class='muted'>No rows available.</p>"
    head = "".join(f"<th>{html.escape(h)}</th>" for h in headers)
    body = []
    for row in rows:
        body.append("<tr>" + "".join(f"<td>{cell}</td>" for cell in row) + "</tr>")
    return f"<table class='{klass}'><thead><tr>{head}</tr></thead><tbody>{''.join(body)}</tbody></table>"


def resolve_relative_artifact(root_dir, value):
    value = clean_cell(value)
    if not root_dir or not value or value in ("NA", "stub"):
        return None
    basename = os.path.basename(value)
    direct = os.path.join(root_dir, basename)
    if os.path.exists(direct):
        return direct
    nested = os.path.join(os.path.dirname(root_dir), value)
    if os.path.exists(nested):
        return nested
    full = os.path.join(root_dir, value)
    if os.path.exists(full):
        return full
    return None


def read_methyl_bed(path):
    rows = []
    if not path or not os.path.exists(path):
        return rows
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        for line in fh:
            if not line.strip() or line.startswith("track"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 11:
                continue
            rows.append(
                {
                    "chrom": parts[0],
                    "start": int(parts[1]),
                    "end": int(parts[2]),
                    "coverage": to_float(parts[9], 0.0),
                    "pct": to_float(parts[10], 0.0),
                }
            )
    return rows


def weighted_methyl_pct(rows):
    weighted = 0.0
    coverage = 0.0
    for row in rows:
        weighted += row["coverage"] * row["pct"]
        coverage += row["coverage"]
    return (weighted / coverage) if coverage else None


def summarize_windows(rows, start, end, window_bp):
    windows = []
    if not rows or end <= start or window_bp <= 0:
        return windows
    n = max(1, ((end - start) + window_bp - 1) // window_bp)
    sums = [0.0] * n
    covs = [0.0] * n
    for row in rows:
        idx = min(n - 1, max(0, (row["start"] - start) // window_bp))
        sums[idx] += row["coverage"] * row["pct"]
        covs[idx] += row["coverage"]
    for idx in range(n):
        win_start = start + idx * window_bp
        win_end = min(end, win_start + window_bp)
        pct = (sums[idx] / covs[idx]) if covs[idx] else None
        windows.append(
            {
                "label": idx + 1,
                "start": win_start,
                "end": win_end,
                "pct": pct,
                "coverage": covs[idx],
            }
        )
    return windows


def parse_bed_annotations(path):
    rows = []
    if not path or not os.path.exists(path):
        return rows
    with open(path, "r", encoding="utf-8") as fh:
        for raw in fh:
            if not raw.strip() or raw.startswith("#"):
                continue
            chrom, start, end, *rest = raw.rstrip("\n").split("\t")
            rows.append(
                {
                    "chrom": chrom,
                    "start": int(start),
                    "end": int(end),
                    "label": rest[0] if rest else f"{chrom}:{start}-{end}",
                }
            )
    return rows


def methyl_profile_svg(rows, annotations, title):
    if not rows:
        return "<p class='muted'>No locus methylation pileup was available.</p>"
    start = min(row["start"] for row in rows)
    end = max(row["end"] for row in rows)
    bins = summarize_windows(rows, start, end, max(1, (end - start) // 48))
    width = 1040
    height = 260
    left = 58
    right = 28
    top = 32
    bottom = 44
    plot_w = width - left - right
    plot_h = height - top - bottom
    pts = []
    bars = []
    for idx, bucket in enumerate(bins):
        frac = idx / max(1, len(bins) - 1)
        x = left + frac * plot_w
        pct = bucket["pct"] if bucket["pct"] is not None else 0.0
        y = top + plot_h - (pct / 100.0) * plot_h
        pts.append(f"{x:.1f},{y:.1f}")
        if idx < len(bins) - 1:
            nxt = left + ((idx + 1) / max(1, len(bins) - 1)) * plot_w
            bars.append(
                f"<rect x='{x:.1f}' y='{y:.1f}' width='{max(3.0, nxt - x)}' height='{top + plot_h - y:.1f}' fill='rgba(188,84,52,0.14)'/>"
            )
    anno_html = []
    for anno in annotations:
        if anno["chrom"] != rows[0]["chrom"]:
            continue
        if anno["end"] < start or anno["start"] > end:
            continue
        x1 = left + ((max(start, anno["start"]) - start) / max(1, end - start)) * plot_w
        x2 = left + ((min(end, anno["end"]) - start) / max(1, end - start)) * plot_w
        label = anno["label"]
        tone = "rgba(38,120,108,0.16)" if "D4Z4" in label or "DUX4" in label else "rgba(181,138,43,0.18)"
        anno_html.append(
            f"<rect x='{x1:.1f}' y='{top - 4}' width='{max(2.0, x2 - x1):.1f}' height='{plot_h + 8}' fill='{tone}'/>"
        )
        anno_html.append(
            f"<text x='{x1 + 4:.1f}' y='{top + 12}' class='anno'>{html.escape(label)}</text>"
        )
    y_ticks = []
    for pct in (0, 25, 50, 75, 100):
        y = top + plot_h - (pct / 100.0) * plot_h
        y_ticks.append(f"<line x1='{left}' y1='{y:.1f}' x2='{width - right}' y2='{y:.1f}' class='gridline'/>")
        y_ticks.append(f"<text x='{left - 10}' y='{y + 4:.1f}' text-anchor='end' class='axis'>{pct}%</text>")
    return (
        f"<svg viewBox='0 0 {width} {height}' class='methyl-svg' role='img' aria-label='{html.escape(title)}'>"
        f"<text x='{left}' y='22' class='title'>{html.escape(title)}</text>"
        + "".join(anno_html)
        + "".join(y_ticks)
        + "".join(bars)
        + f"<polyline fill='none' stroke='#bc5434' stroke-width='3' points='{' '.join(pts)}'/>"
        + f"<line x1='{left}' y1='{top + plot_h:.1f}' x2='{width - right}' y2='{top + plot_h:.1f}' class='axisline'/>"
        + f"<text x='{left}' y='{height - 10}' class='axis'>{rows[0]['chrom']}:{start:,}</text>"
        + f"<text x='{width - right}' y='{height - 10}' text-anchor='end' class='axis'>{rows[0]['chrom']}:{end:,}</text>"
        + "</svg>"
    )


def methyl_tile_grid(windows, title):
    if not windows:
        return "<p class='muted'>No repeat-sized methylation bins were available.</p>"
    cells = []
    for bucket in windows:
        pct = bucket["pct"]
        coverage = bucket["coverage"]
        if pct is None:
            tone = "#ece8df"
            text = "NA"
        else:
            green = int(230 - min(100, pct) * 1.0)
            red = int(242 - min(100, pct) * 0.35)
            blue = int(218 - min(100, pct) * 0.95)
            tone = f"rgb({red},{green},{blue})"
            text = f"{pct:.1f}%"
        cells.append(
            "<div class='methyl-tile'>"
            f"<div class='methyl-swatch' style='background:{tone}'></div>"
            f"<div class='methyl-label'>Bin {bucket['label']}</div>"
            f"<div class='methyl-value'>{html.escape(text)}</div>"
            f"<div class='methyl-note'>cov {coverage:.1f}</div>"
            "</div>"
        )
    return f"<div><h3>{html.escape(title)}</h3><div class='methyl-tiles'>{''.join(cells)}</div></div>"


def pick_representative_rows(rows, limit=18):
    def rank(row):
        tokens = parse_repeat_tokens(row.get("repeat_sequence"))
        return (
            0 if has_proximal_anchor(tokens) and has_distal_anchor(tokens) else 1,
            -to_int(row.get("RU_count"), -1),
            0 if clean_cell(row.get("PAS.type")) == "4qA_PAS" else 1,
            row.get("read.id", ""),
        )

    return sorted(rows, key=rank)[:limit]


def coverage_rows_from_tsv(path):
    rows = []
    if not path or not os.path.exists(path):
        return rows
    with open(path, "r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            rows.append({clean_cell(k): clean_cell(v) for k, v in row.items()})
    return rows


def blast_anchor_counts(blast_rows):
    counts = Counter()
    by_read = defaultdict(set)
    for row in blast_rows:
        read_id = clean_cell(row.get("read.id"))
        region = clean_cell(row.get("region"))
        if read_id and region:
            by_read[read_id].add(region)
    for read_regions in by_read.values():
        if "A_pLAM_4qA" in read_regions:
            counts["4qA pLAM"] += 1
        if "A_D4F104S1_4qA" in read_regions:
            counts["4qA D4F104S1"] += 1
        if "CLUHP4_201_exon1" in read_regions:
            counts["CLUHP4 exon1"] += 1
        if "DUX4_end" in read_regions:
            counts["DUX4 end"] += 1
    return counts


def main() -> int:
    (
        sample_id,
        classification_dir,
        subset_dir,
        flagstat_txt,
        coverage_tsv,
        hap_tsv,
        methylation_dir,
        methyl_tsv,
        variant_tsv,
        contraction_threshold,
        locus_bed,
        html_out,
        summary_out,
    ) = sys.argv[1:]
    contraction_threshold = int(contraction_threshold)

    category_files = {
        "4qA_all": "A_4qA_all-reads.csv",
        "4qA_complete": "A_4qA_complete-reads.csv",
        "4qB_all": "B_4qB_all-reads.csv",
        "4qB_complete": "B_4qB_complete-reads.csv",
        "chr10_all": "chr10_all-reads.csv",
        "chr10_complete": "chr10_complete-reads.csv",
        "chimeric": "chimeric_reads.csv",
        "chr4_undefined": "chr4_undefined_all-reads.csv",
    }
    data = {
        name: read_semicolon_csv(os.path.join(classification_dir, filename))
        for name, filename in category_files.items()
    }
    counts = {name: len(rows) for name, rows in data.items()}

    rows_4qa_all = data["4qA_all"]
    rows_4qa_complete = data["4qA_complete"]
    ru_4qa_all = [to_int(row.get("RU_count"), -1) for row in rows_4qa_all if to_int(row.get("RU_count"), -1) > 0]
    ru_4qa_complete = [to_int(row.get("RU_count"), -1) for row in rows_4qa_complete if to_int(row.get("RU_count"), -1) > 0]
    shortest_4qa = min(ru_4qa_complete) if ru_4qa_complete else None
    best_partial_4qa = max(ru_4qa_all) if ru_4qa_all else None

    contraction_status, contraction_note = contraction_assessment(rows_4qa_complete, contraction_threshold)
    contracted = contraction_status == "supported"

    if contraction_status == "supported":
        headline = f"Complete 4qA read support for an allele at or below {contraction_threshold} repeat units was observed."
    elif contraction_status == "not_supported":
        headline = f"Complete 4qA reads were present, but none supported a contraction at {contraction_threshold} repeat units or fewer."
    elif counts["4qA_all"] > 0:
        headline = "Partial 4qA-supporting reads were present, but no complete 4qA-spanning read was classified."
    else:
        headline = "No 4qA-supporting reads were classified in this sample."

    flagstat_headline = "NA"
    if os.path.exists(flagstat_txt):
        with open(flagstat_txt, "r", encoding="utf-8") as fh:
            flagstat_headline = clean_cell(fh.readline())

    coverage_rows = coverage_rows_from_tsv(coverage_tsv)
    coverage_focus = [
        row for row in coverage_rows if row.get("#chrom") == "chr4" and row.get("label") in ("D4Z4_4q35", "DUX4_gene_body_chr4", "pLAM_4qA", "DUX4_end_chr4")
    ]

    hap_rows = read_tsv(hap_tsv)
    has_hp_tags = any((row.get("has_hp_tags", "").lower() == "true") for row in hap_rows)
    dominant_hp_4qa = "NA"
    hp_counts = defaultdict(int)
    for row in hap_rows:
        if row.get("subset_name") == "4qA_all":
            hp_counts[row.get("hp", "NA")] += to_int(row.get("count"), 0)
    if hp_counts:
        dominant_hp_4qa = max(hp_counts.items(), key=lambda kv: kv[1])[0]

    methyl_rows = read_tsv(methyl_tsv)
    methyl_by_subset = {row.get("subset_name"): row for row in methyl_rows if row.get("subset_name")}
    methyl_4qa_all = methyl_by_subset.get("4qA_all", {})
    methyl_bed_path = resolve_relative_artifact(methylation_dir, methyl_4qa_all.get("pileup_bed_gz"))
    methyl_records = read_methyl_bed(methyl_bed_path)
    methyl_mean = weighted_methyl_pct(methyl_records)
    locus_annotations = parse_bed_annotations(locus_bed)
    repeat_windows = []
    if methyl_records:
        meth_start = min(row["start"] for row in methyl_records)
        meth_end = max(row["end"] for row in methyl_records)
        repeat_windows = summarize_windows(methyl_records, meth_start, meth_end, REPEAT_BP)

    variant_rows = read_tsv(variant_tsv)
    variant = variant_rows[0] if variant_rows else {
        "status": "unknown",
        "total_clair3_variants": "0",
        "relevant_variant_records": "0",
        "clinvar_fshd_hits": "0",
        "hg38_sv_records": "0",
        "t2t_sv_records": "0",
        "relevant_genes": "NA",
    }

    pas_counter = Counter(
        clean_cell(row.get("PAS.type", "NA")) or "NA" for row in rows_4qa_all
    )
    proximal_rows = [row for row in rows_4qa_all if has_proximal_anchor(parse_repeat_tokens(row.get("repeat_sequence")))]
    distal_rows = [row for row in rows_4qa_all if has_distal_anchor(parse_repeat_tokens(row.get("repeat_sequence")))]
    best_proximal_ru = max((to_int(row.get("RU_count"), -1) for row in proximal_rows), default=None)
    best_distal_ru = max((to_int(row.get("RU_count"), -1) for row in distal_rows), default=None)

    subset_manifest = os.path.join(subset_dir, f"{sample_id}.classified-subsets.manifest.tsv")
    subset_rows = read_tsv(subset_manifest)

    blast_4qa_rows = read_semicolon_csv(os.path.join(classification_dir, "blast_results", "4qA_reads_blast.csv"))
    anchor_counts = blast_anchor_counts(blast_4qa_rows)

    summary_fields = {
        "sample_id": sample_id,
        "contracted_4qA_observed": str(contracted).lower(),
        "shortest_complete_4qA_ru": shortest_4qa if shortest_4qa is not None else "NA",
        "complete_4qA_reads": counts["4qA_complete"],
        "all_4qA_reads": counts["4qA_all"],
        "complete_4qB_reads": counts["4qB_complete"],
        "all_4qB_reads": counts["4qB_all"],
        "complete_chr10_reads": counts["chr10_complete"],
        "all_chr10_reads": counts["chr10_all"],
        "chimeric_reads": counts["chimeric"],
        "has_hp_tags": str(has_hp_tags).lower(),
        "dominant_hp_4qA": dominant_hp_4qa,
        "methyl_4qA_all_mean_pct": methyl_by_subset.get("4qA_all", {}).get("mean_percent_modified", "NA"),
        "methyl_4qA_complete_mean_pct": methyl_by_subset.get("4qA_complete", {}).get("mean_percent_modified", "NA"),
        "methyl_chimeric_mean_pct": methyl_by_subset.get("chimeric", {}).get("mean_percent_modified", "NA"),
        "variant_status": variant.get("status", "unknown"),
        "variant_total_clair3_variants": variant.get("total_clair3_variants", "0"),
        "variant_relevant_records": variant.get("relevant_variant_records", "0"),
        "variant_clinvar_fshd_hits": variant.get("clinvar_fshd_hits", "0"),
        "variant_hg38_sv_records": variant.get("hg38_sv_records", "0"),
        "variant_t2t_sv_records": variant.get("t2t_sv_records", "0"),
        "variant_relevant_genes": variant.get("relevant_genes", "NA"),
        "contraction_assessment": contraction_status,
        "best_partial_4qA_ru": best_partial_4qa if best_partial_4qa is not None else "NA",
        "best_proximal_4qA_ru": best_proximal_ru if best_proximal_ru is not None else "NA",
        "best_distal_4qA_ru": best_distal_ru if best_distal_ru is not None else "NA",
        "locus_4qA_methyl_mean_pct": f"{methyl_mean:.3f}" if methyl_mean is not None else "NA",
    }

    with open(summary_out, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(summary_fields.keys()), delimiter="\t")
        writer.writeheader()
        writer.writerow(summary_fields)

    representative_rows = pick_representative_rows(rows_4qa_all, limit=18)
    representative_table = table_html(
        ["Read", "RU span", "Anchor", "PAS", "Repeat architecture"],
        [
            [
                html.escape(clean_cell(row.get("read.id"))),
                html.escape(clean_cell(row.get("RU_count")) or "NA"),
                html.escape(anchor_mode(parse_repeat_tokens(row.get("repeat_sequence")))),
                html.escape(clean_cell(row.get("PAS.type")) or "NA"),
                repeat_strip_html(parse_repeat_tokens(row.get("repeat_sequence"))),
            ]
            for row in representative_rows
        ],
    )

    anchor_table = table_html(
        ["Signal", "Reads"],
        [
            [html.escape(label), html.escape(str(anchor_counts.get(label, 0)))]
            for label in ("4qA D4F104S1", "CLUHP4 exon1", "4qA pLAM", "DUX4 end")
        ],
        "compact",
    )

    coverage_table = table_html(
        ["Label", "Reads", "Covered bp", "Breadth", "Mean depth"],
        [
            [
                html.escape(row.get("label", "")),
                html.escape(row.get("numreads", "")),
                html.escape(row.get("covbases", "")),
                html.escape(row.get("coverage", "")),
                html.escape(row.get("meandepth", "")),
            ]
            for row in coverage_focus
        ],
        "compact",
    )

    subset_table = table_html(
        ["Subset", "Read count", "Status"],
        [
            [
                html.escape(row.get("subset_name", "")),
                html.escape(row.get("read_count", "0")),
                html.escape(row.get("status", "NA")),
            ]
            for row in subset_rows
        ],
        "compact",
    )

    pas_summary = ", ".join(f"{k}: {v}" for k, v in sorted(pas_counter.items()) if k) or "No PAS assignments available."

    variant_genes = clean_cell(variant.get("relevant_genes", "NA"))
    variant_note = truncate(variant_genes, 180)

    html_body = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>{html.escape(sample_id)} FSHD VNTR report</title>
  <style>
    :root {{
      --ink: #261f1a;
      --paper: #f6efe5;
      --panel: rgba(255,255,255,0.80);
      --edge: rgba(69,54,42,0.12);
      --burnt: #bc5434;
      --burnt-soft: #f1c6b8;
      --sage: #2f7568;
      --gold: #b6882b;
      --slate: #6e635b;
      --shadow: 0 18px 44px rgba(80, 63, 49, 0.10);
    }}
    * {{ box-sizing: border-box; }}
    body {{
      margin: 0;
      font-family: "Georgia", "Iowan Old Style", serif;
      color: var(--ink);
      background:
        radial-gradient(circle at top left, rgba(188,84,52,0.20), transparent 28%),
        radial-gradient(circle at top right, rgba(47,117,104,0.18), transparent 22%),
        linear-gradient(180deg, #fcf7f1 0%, #f4ebdf 100%);
    }}
    .shell {{ max-width: 1280px; margin: 0 auto; padding: 34px 22px 48px; }}
    .hero, .section {{
      background: var(--panel);
      border: 1px solid var(--edge);
      border-radius: 28px;
      box-shadow: var(--shadow);
    }}
    .hero {{
      padding: 28px 30px 30px;
      position: relative;
      overflow: hidden;
    }}
    .hero::after {{
      content: "";
      position: absolute;
      width: 320px;
      height: 320px;
      border-radius: 50%;
      right: -120px;
      bottom: -170px;
      background: radial-gradient(circle, rgba(188,84,52,0.14), transparent 72%);
    }}
    h1 {{ margin: 0 0 10px; font-size: 2.6rem; line-height: 1.02; }}
    h2 {{ margin: 0 0 14px; font-size: 1.25rem; }}
    h3 {{ margin: 0 0 10px; font-size: 1rem; }}
    p {{ line-height: 1.58; }}
    .lede {{ max-width: 860px; font-size: 1.08rem; color: #3d3128; }}
    .badge-row {{ display: flex; flex-wrap: wrap; gap: 10px; margin-bottom: 14px; }}
    .badge {{
      display: inline-flex;
      align-items: center;
      padding: 8px 14px;
      border-radius: 999px;
      font-size: 0.92rem;
      border: 1px solid transparent;
    }}
    .badge.partial {{ background: rgba(182,136,43,0.14); color: #7d5a07; border-color: rgba(182,136,43,0.20); }}
    .badge.good {{ background: rgba(47,117,104,0.14); color: #1f6257; border-color: rgba(47,117,104,0.22); }}
    .badge.warn {{ background: rgba(188,84,52,0.12); color: #8f3518; border-color: rgba(188,84,52,0.18); }}
    .hero-note {{
      margin-top: 14px;
      border-left: 4px solid var(--burnt);
      background: rgba(188,84,52,0.08);
      padding: 14px 16px;
      border-radius: 14px;
    }}
    .metric-grid, .dual-grid {{
      display: grid;
      gap: 16px;
    }}
    .metric-grid {{
      margin-top: 22px;
      grid-template-columns: repeat(auto-fit, minmax(170px, 1fr));
    }}
    .dual-grid {{
      grid-template-columns: repeat(auto-fit, minmax(380px, 1fr));
    }}
    .metric {{
      background: rgba(255,255,255,0.82);
      border: 1px solid var(--edge);
      border-radius: 18px;
      padding: 16px 18px;
    }}
    .metric.good {{ background: linear-gradient(180deg, rgba(239,249,245,0.96), rgba(255,255,255,0.84)); }}
    .metric.warn {{ background: linear-gradient(180deg, rgba(253,244,239,0.96), rgba(255,255,255,0.84)); }}
    .metric.soft {{ background: linear-gradient(180deg, rgba(252,248,239,0.96), rgba(255,255,255,0.84)); }}
    .metric-label {{ font-size: 0.82rem; text-transform: uppercase; letter-spacing: 0.08em; color: var(--slate); }}
    .metric-value {{ font-size: 2rem; margin-top: 8px; font-weight: 700; }}
    .metric-note {{ margin-top: 6px; font-size: 0.92rem; color: var(--slate); }}
    .section {{ margin-top: 24px; padding: 24px; }}
    .muted {{ color: var(--slate); }}
    table {{
      width: 100%;
      border-collapse: collapse;
      margin-top: 14px;
      font-size: 0.94rem;
      background: rgba(255,255,255,0.64);
      border-radius: 16px;
      overflow: hidden;
    }}
    th, td {{
      padding: 10px 12px;
      border-bottom: 1px solid rgba(69,54,42,0.10);
      text-align: left;
      vertical-align: top;
    }}
    th {{
      font-size: 0.82rem;
      text-transform: uppercase;
      letter-spacing: 0.05em;
      color: var(--slate);
      background: rgba(190,165,132,0.10);
    }}
    tr:last-child td {{ border-bottom: 0; }}
    .compact td, .compact th {{ padding: 8px 10px; font-size: 0.90rem; }}
    .histogram, .methyl-svg {{ width: 100%; height: auto; border-radius: 18px; background: rgba(255,250,245,0.78); border: 1px solid rgba(69,54,42,0.10); }}
    .title {{ font-size: 16px; fill: #45362a; font-weight: 700; }}
    .axis {{ font-size: 11px; fill: #6e635b; }}
    .count {{ font-size: 12px; fill: #45362a; }}
    .gridline {{ stroke: rgba(69,54,42,0.10); stroke-width: 1; }}
    .axisline {{ stroke: rgba(69,54,42,0.26); stroke-width: 1.5; }}
    .anno {{ font-size: 10px; fill: #355f57; }}
    .repeat-strip {{
      display: flex;
      flex-wrap: wrap;
      gap: 6px;
      min-width: 280px;
    }}
    .unit {{
      display: inline-flex;
      align-items: center;
      justify-content: center;
      min-height: 32px;
      padding: 6px 10px;
      border-radius: 10px;
      font-size: 0.82rem;
      border: 1px solid rgba(69,54,42,0.10);
      white-space: nowrap;
    }}
    .repeat {{ background: #f2d3c9; }}
    .repeat-short {{ background: #f7e1ce; border-color: rgba(182,136,43,0.22); }}
    .repeat-long {{ background: #dbeeea; border-color: rgba(47,117,104,0.20); }}
    .anchor-prox {{ background: #dce9f7; }}
    .anchor-dist {{ background: #dff1d8; }}
    .offtarget {{ background: #ece7e0; color: #776e67; }}
    .soft-end {{ background: #f3efe8; color: #6e635b; }}
    .misc {{ background: #eee8fb; }}
    .legend {{
      display: flex;
      flex-wrap: wrap;
      gap: 10px 16px;
      margin-top: 12px;
      font-size: 0.9rem;
      color: var(--slate);
    }}
    .legend-chip {{
      display: inline-block;
      width: 14px;
      height: 14px;
      border-radius: 4px;
      vertical-align: text-bottom;
      margin-right: 6px;
      border: 1px solid rgba(69,54,42,0.08);
    }}
    .methyl-tiles {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(110px, 1fr));
      gap: 10px;
      margin-top: 12px;
    }}
    .methyl-tile {{
      background: rgba(255,255,255,0.74);
      border: 1px solid var(--edge);
      border-radius: 16px;
      padding: 10px;
    }}
    .methyl-swatch {{
      height: 10px;
      border-radius: 999px;
      margin-bottom: 8px;
    }}
    .methyl-label {{ font-size: 0.78rem; text-transform: uppercase; letter-spacing: 0.06em; color: var(--slate); }}
    .methyl-value {{ font-size: 1.1rem; margin-top: 6px; font-weight: 700; }}
    .methyl-note {{ font-size: 0.82rem; color: var(--slate); margin-top: 4px; }}
    .footnote {{
      margin-top: 12px;
      font-size: 0.90rem;
      color: var(--slate);
    }}
  </style>
</head>
<body>
  <div class="shell">
    <section class="hero">
      <div class="badge-row">
        <span class="badge partial">FSHD D4Z4 / DUX4 locus</span>
        <span class="badge {'good' if contracted else 'warn' if contraction_status == 'not_supported' else 'partial'}">{html.escape(contraction_status.replace('_', ' '))}</span>
        <span class="badge {'good' if counts['4qA_complete'] else 'partial'}">{counts['4qA_complete']} complete 4qA reads</span>
      </div>
      <h1>{html.escape(sample_id)}</h1>
      <p class="lede">{html.escape(headline)}</p>
      <div class="hero-note">
        <strong>Contraction assessment:</strong> {html.escape(contraction_note)}
      </div>
      <div class="metric-grid">
        {metric_card("4qA-supporting reads", counts["4qA_all"], "All permissive 4qA-classified reads", "good" if counts["4qA_all"] else "soft")}
        {metric_card("Best observed 4qA span", best_partial_4qa if best_partial_4qa is not None else "NA", "Largest partial or complete RU count on 4qA reads", "soft")}
        {metric_card("4qA PAS-positive reads", pas_counter.get("4qA_PAS", 0), "ATTAAA-compatible pLAM support", "good" if pas_counter.get("4qA_PAS", 0) else "soft")}
        {metric_card("4qA locus methylation", f"{methyl_mean:.1f}%" if methyl_mean is not None else methyl_4qa_all.get("mean_percent_modified", "NA"), "Weighted mean across modkit pileup positions", "good" if methyl_mean is not None else "soft")}
        {metric_card("Complete-read contraction", "yes" if contracted else "no" if contraction_status == "not_supported" else "inconclusive", f"Threshold: <= {contraction_threshold} RU", "good" if contracted else "warn" if contraction_status == "not_supported" else "soft")}
        {metric_card("Variant branch", variant.get("status", "unknown"), truncate(variant_note, 60), "soft")}
      </div>
      <p class="footnote"><strong>Locus BAM flagstat:</strong> {html.escape(flagstat_headline)}</p>
    </section>

    <section class="section">
      <h2>4qA VNTR Evidence</h2>
      <div class="dual-grid">
        <div>
          {hist_svg(ru_4qa_all, "#bc5434", "Observed repeat-unit spans on 4qA-supporting reads")}
          <p class="footnote">This sample has only partial 4qA reads. The histogram shows observed repeat spans on permissive reads, not full allele sizes.</p>
        </div>
        <div>
          <h3>Anchor support across 4qA-supporting reads</h3>
          {anchor_table}
          <p class="footnote"><strong>PAS summary:</strong> {html.escape(pas_summary)}</p>
        </div>
      </div>
      <div class="legend">
        <span><span class="legend-chip" style="background:#dce9f7"></span>proximal anchors</span>
        <span><span class="legend-chip" style="background:#f2d3c9"></span>c4 repeat units</span>
        <span><span class="legend-chip" style="background:#f7e1ce"></span>short distal repeat (`c4S`)</span>
        <span><span class="legend-chip" style="background:#dff1d8"></span>distal pLAM / DUX4 markers</span>
        <span><span class="legend-chip" style="background:#ece7e0"></span>chr10-like / off-target marker</span>
      </div>
      {representative_table}
    </section>

    <section class="section">
      <h2>Contraction at {contraction_threshold} Repeat Units</h2>
      <div class="metric-grid">
        {metric_card("Complete 4qA reads at <= threshold", sum(1 for row in rows_4qa_complete if 0 < to_int(row.get("RU_count"), -1) <= contraction_threshold), "Direct contraction support", "good" if contracted else "soft")}
        {metric_card("Best proximal-side span", best_proximal_ru if best_proximal_ru is not None else "NA", "Reads carrying CLUHP4 / D4F104S1 anchors", "soft")}
        {metric_card("Best distal-side span", best_distal_ru if best_distal_ru is not None else "NA", "Reads carrying pLAM / DUX4-end anchors", "soft")}
        {metric_card("Interpretation", "inconclusive" if contraction_status == "inconclusive" else ("supported" if contracted else "not supported"), "Complete 4qA-spanning reads are required for a direct call", "warn" if contraction_status == "inconclusive" else "good" if contracted else "soft")}
      </div>
      <p class="footnote">Partial 4qA reads can show how many repeat units were observed from one side of the locus, but they do not by themselves prove a full allele length. This sample currently sits in that partial-evidence category.</p>
    </section>

    <section class="section">
      <h2>Locus Methylation</h2>
      <div class="dual-grid">
        <div>
          {methyl_profile_svg(methyl_records, locus_annotations, "4qA locus methylation across observed chr4 positions")}
          <p class="footnote">This panel uses the modkit pileup emitted from the 4qA subset BAM. Shaded bands mark annotated chr4 D4Z4/DUX4 landmarks from the workflow BED.</p>
        </div>
        <div>
          {methyl_tile_grid(repeat_windows, "Repeat-sized methylation bins (~3.3 kb windows)")}
          <p class="footnote">These bins are reference-coordinate windows sized like a D4Z4 repeat unit. When only partial reads are present, they are best interpreted as an approximate locus methylation profile rather than exact allele-resolved repeat units.</p>
        </div>
      </div>
    </section>

    <section class="section">
      <h2>Support Tables</h2>
      <div class="dual-grid">
        <div>
          <h3>Classified subsets</h3>
          {subset_table}
        </div>
        <div>
          <h3>Relevant locus intervals</h3>
          {coverage_table}
        </div>
      </div>
    </section>

    <section class="section">
      <h2>Variant Context</h2>
      <p class="muted">The variant branch is preserved here as context, but this report is intentionally centered on the VNTR locus and permissive 4qA evidence.</p>
      {table_html(
          ["Status", "Clair3 calls", "Relevant SNV/indels", "ClinVar FSHD hits", "HG38 SVs", "T2T SVs", "Genes"],
          [[
              html.escape(variant.get("status", "unknown")),
              html.escape(variant.get("total_clair3_variants", "0")),
              html.escape(variant.get("relevant_variant_records", "0")),
              html.escape(variant.get("clinvar_fshd_hits", "0")),
              html.escape(variant.get("hg38_sv_records", "0")),
              html.escape(variant.get("t2t_sv_records", "0")),
              html.escape(truncate(variant_genes, 140)),
          ]],
          "compact",
      )}
    </section>
  </div>
</body>
</html>
"""

    with open(html_out, "w", encoding="utf-8") as fh:
        fh.write(html_body)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
