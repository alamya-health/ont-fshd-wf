#!/usr/bin/env python3

import csv
import gzip
import html
import math
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


def ru_values(rows):
    return [to_int(row.get("RU_count"), -1) for row in rows if to_int(row.get("RU_count"), -1) > 0]


def summarize_allele_candidate(label, all_rows, complete_rows, collective_length=None):
    ru_all = ru_values(all_rows)
    ru_complete = ru_values(complete_rows)
    best_partial = max(ru_all) if ru_all else None
    proximal_support = 0
    distal_support = 0
    for row in all_rows:
        tokens = parse_repeat_tokens(row.get("repeat_sequence"))
        if has_proximal_anchor(tokens):
            proximal_support += 1
        if has_distal_anchor(tokens):
            distal_support += 1

    if ru_complete:
        sorted_complete = sorted(ru_complete)
        midpoint = len(sorted_complete) // 2
        estimate_ru = sorted_complete[midpoint]
        spread = f"{sorted_complete[0]}-{sorted_complete[-1]}" if len(sorted_complete) > 1 else str(estimate_ru)
        return {
            "label": label,
            "estimate_ru": estimate_ru,
            "display_ru": str(estimate_ru),
            "method": "direct_complete_read",
            "confidence": "high",
            "support_reads": len(complete_rows),
            "best_partial_ru": best_partial,
            "proximal_support_reads": proximal_support,
            "distal_support_reads": distal_support,
            "note": f"{label} has {len(complete_rows)} complete spanning read(s); complete-read RU counts ranged {spread}.",
        }

    if collective_length and collective_length.get("estimate_ru") is not None:
        return {
            "label": label,
            "estimate_ru": collective_length["estimate_ru"],
            "display_ru": str(collective_length["estimate_ru"]),
            "method": collective_length["method"],
            "confidence": collective_length["confidence"],
            "support_reads": len(all_rows),
            "best_partial_ru": best_partial,
            "proximal_support_reads": proximal_support,
            "distal_support_reads": distal_support,
            "note": collective_length["note"],
        }

    if best_partial is not None:
        side_note = []
        if proximal_support:
            side_note.append(f"{proximal_support} proximal-side")
        if distal_support:
            side_note.append(f"{distal_support} distal-side")
        side_text = " and ".join(side_note) if side_note else "partial"
        return {
            "label": label,
            "estimate_ru": None,
            "display_ru": f">={best_partial}",
            "method": "partial_read_lower_bound",
            "confidence": "low",
            "support_reads": len(all_rows),
            "best_partial_ru": best_partial,
            "proximal_support_reads": proximal_support,
            "distal_support_reads": distal_support,
            "note": f"{label} is only supported by {side_text} reads, so the current result is a lower bound of at least {best_partial} repeat units rather than a full allele-length estimate.",
        }

    return {
        "label": label,
        "estimate_ru": None,
        "display_ru": "NA",
        "method": "unresolved",
        "confidence": "none",
        "support_reads": len(all_rows),
        "best_partial_ru": best_partial,
        "proximal_support_reads": proximal_support,
        "distal_support_reads": distal_support,
        "note": f"No usable {label} repeat-span evidence was available for allele-length inference.",
    }


def size_call_lower_bound(allele):
    display = clean_cell(allele.get("display_ru"))
    if not display or display == "NA":
        return None
    if display.startswith(">="):
        try:
            return int(display[2:])
        except Exception:
            return None
    try:
        return int(display)
    except Exception:
        return None


def normalize_size_display(value, fallback=None):
    text = clean_cell(value)
    if text in ("", "NA", "None", "none"):
        return fallback if fallback is not None else "NA"
    return text


def split_pooled_repeat_candidates(label_prefix, rows, proximal_tokens, distal_tokens):
    usable = []
    for row in rows:
        ru = to_int(row.get("RU_count"), -1)
        if ru <= 0:
            continue
        tokens = parse_repeat_tokens(row.get("repeat_sequence"))
        prox = any(tok in token for token in tokens for tok in proximal_tokens)
        dist = any(tok in token for token in tokens for tok in distal_tokens)
        side = 1 if prox and not dist else -1 if dist and not prox else 0
        usable.append({"row": row, "ru": ru, "side": side})

    if len(usable) < 4:
        return []

    seed_a = max(usable, key=lambda item: (item["ru"], item["side"]))
    seed_b = max(usable, key=lambda item: (item["ru"], -item["side"]))
    centroids = [(seed_a["ru"], seed_a["side"]), (seed_b["ru"], seed_b["side"])]
    groups = {0: usable, 1: []}

    for _ in range(12):
        groups = {0: [], 1: []}
        for item in usable:
            distances = [
                (item["ru"] - ru_c) ** 2 + (1.5 * (item["side"] - side_c)) ** 2
                for ru_c, side_c in centroids
            ]
            groups[0 if distances[0] <= distances[1] else 1].append(item)

        updated = []
        for idx in (0, 1):
            if groups[idx]:
                updated.append(
                    (
                        sum(item["ru"] for item in groups[idx]) / len(groups[idx]),
                        sum(item["side"] for item in groups[idx]) / len(groups[idx]),
                    )
                )
            else:
                updated.append(centroids[idx])
        if updated == centroids:
            break
        centroids = updated

    candidates = []
    for idx, items in sorted(groups.items(), key=lambda kv: max((x["ru"] for x in kv[1]), default=-1), reverse=True):
        if not items:
            continue
        max_ru = max(item["ru"] for item in items)
        prox_support = sum(1 for item in items if item["side"] > 0)
        dist_support = sum(1 for item in items if item["side"] < 0)
        candidates.append(
            {
                "label": f"{label_prefix} candidate {len(candidates) + 1}",
                "estimate_ru": None,
                "display_ru": f">={max_ru}",
                "method": "two_cluster_partial_ru",
                "confidence": "low",
                "support_reads": len(items),
                "best_partial_ru": max_ru,
                "proximal_support_reads": prox_support,
                "distal_support_reads": dist_support,
                "note": f"{label_prefix} pooled reads were heuristically split into two low-confidence clusters; this cluster supports a lower bound of at least {max_ru} repeat units from {len(items)} partial reads.",
            }
        )
    return candidates


def hist_svg(values, color, label):
    if not values:
        return "<p class='muted'>No repeat-unit measurements were available for this panel.</p>"
    counts = Counter(values)
    xs = sorted(counts)
    max_count = max(counts.values())
    bar_width = 44
    gap = 12
    width = max(460, len(xs) * (bar_width + gap) + 70)
    height = 270
    bars = []
    y_ticks = []
    for frac, tick in ((0.0, 0), (0.5, math.ceil(max_count / 2)), (1.0, max_count)):
        y = 196 - frac * 132
        y_ticks.append(f"<line x1='26' y1='{y:.1f}' x2='{width - 18}' y2='{y:.1f}' class='gridline'/>")
        y_ticks.append(f"<text x='18' y='{y + 4:.1f}' text-anchor='end' class='axis'>{tick}</text>")
    for idx, x in enumerate(xs):
        count = counts[x]
        bar_h = int((count / max_count) * 132)
        xpos = 44 + idx * (bar_width + gap)
        ypos = 196 - bar_h
        bars.append(
            f"<rect x='{xpos}' y='{ypos}' width='{bar_width}' height='{bar_h}' rx='10' fill='{color}' />"
        )
        bars.append(f"<text x='{xpos + bar_width / 2}' y='218' text-anchor='middle' class='axis'>{x}</text>")
        bars.append(
            f"<text x='{xpos + bar_width / 2}' y='{max(28, ypos - 8)}' text-anchor='middle' class='count'>{count}</text>"
        )
    return (
        f"<svg viewBox='0 0 {width} {height}' class='histogram' role='img' aria-label='{html.escape(label)}'>"
        f"<text x='28' y='26' class='title'>{html.escape(label)}</text>"
        + "".join(y_ticks)
        + f"<line x1='26' y1='196' x2='{width - 18}' y2='196' stroke='rgba(58,45,36,0.22)' stroke-width='2'/>"
        + "".join(bars)
        + f"<text x='{width / 2:.1f}' y='{height - 12}' text-anchor='middle' class='axis'>Observed repeat-unit span on classified reads (RU)</text>"
        + f"<text x='16' y='{height / 2:.1f}' text-anchor='middle' class='axis' transform='rotate(-90 16 {height / 2:.1f})'>Read count</text>"
        + f"<rect x='{width - 170}' y='34' width='14' height='14' rx='4' fill='{color}'/>"
        + f"<text x='{width - 150}' y='46' class='axis'>Histogram bars = read support</text>"
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


def methyl_records_for_subset(methyl_by_subset, methylation_dir, subset_name):
    subset_row = methyl_by_subset.get(subset_name, {})
    pileup_path = resolve_relative_artifact(methylation_dir, subset_row.get("pileup_bed_gz"))
    return read_methyl_bed(pileup_path)


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


def annotation_start(annotations, label):
    for row in annotations:
        if row.get("label") == label:
            return row["start"]
    return None


def repeat_unit_matrix(rows, repeat_anchor_start, repeat_bp, bins_per_unit=60, min_occupied_cols=8):
    if not rows or repeat_anchor_start is None:
        return []

    matrix = defaultdict(lambda: [{"weighted": 0.0, "coverage": 0.0} for _ in range(bins_per_unit)])
    occupied = defaultdict(set)
    total_cov = defaultdict(float)
    total_weighted = defaultdict(float)

    for row in rows:
        pos = row["start"]
        repeat_idx = math.floor((pos - repeat_anchor_start) / repeat_bp)
        offset = (pos - repeat_anchor_start) % repeat_bp
        bin_idx = min(bins_per_unit - 1, max(0, int((offset / repeat_bp) * bins_per_unit)))
        matrix[repeat_idx][bin_idx]["weighted"] += row["coverage"] * row["pct"]
        matrix[repeat_idx][bin_idx]["coverage"] += row["coverage"]
        occupied[repeat_idx].add(bin_idx)
        total_cov[repeat_idx] += row["coverage"]
        total_weighted[repeat_idx] += row["coverage"] * row["pct"]

    repeat_ids = sorted(idx for idx in matrix if len(occupied[idx]) >= min_occupied_cols and idx <= 0)
    if not repeat_ids:
        repeat_ids = sorted(idx for idx in matrix if len(occupied[idx]) >= min_occupied_cols)

    repeat_rows = []
    for order, repeat_idx in enumerate(repeat_ids, start=1):
        cells = []
        for cell in matrix[repeat_idx]:
            if cell["coverage"] > 0:
                cells.append(cell["weighted"] / cell["coverage"])
            else:
                cells.append(None)
        pct = (total_weighted[repeat_idx] / total_cov[repeat_idx]) if total_cov[repeat_idx] else None
        repeat_rows.append(
            {
                "repeat_idx": repeat_idx,
                "ordinal": order,
                "occupied_cols": len(occupied[repeat_idx]),
                "coverage": total_cov[repeat_idx],
                "pct": pct,
                "cells": cells,
                "is_distal_repeat": repeat_idx == 0,
            }
        )
    return repeat_rows


def with_repeat_labels(repeat_rows, prefix):
    labelled = []
    for row in repeat_rows:
        copy = dict(row)
        base = "distal" if row.get("is_distal_repeat") else f"RU {row.get('ordinal')}"
        copy["display_label"] = f"{prefix} {base}"
        labelled.append(copy)
    return labelled


def repeat_pileup_svg(repeat_rows, bins_per_unit, title):
    if not repeat_rows:
        return "<p class='muted'>No repeat-folded methylation pileup was available.</p>"

    left = 84
    top = 34
    cell_w = 11
    cell_h = 12
    gap = 2
    width = left + bins_per_unit * cell_w + 56
    height = top + len(repeat_rows) * (cell_h + gap) + 54

    def color_for_pct(pct):
        if pct is None:
            return "#efe8df"
        pct = max(0.0, min(100.0, pct))
        red = int(244 - pct * 0.28)
        green = int(231 - pct * 1.18)
        blue = int(219 - pct * 0.92)
        return f"rgb({red},{green},{blue})"

    blocks = []
    for row_idx, repeat_row in enumerate(repeat_rows):
        y = top + row_idx * (cell_h + gap)
        label = repeat_row.get("display_label")
        if not label:
            label = f"RU {repeat_row['ordinal']}"
            if repeat_row["is_distal_repeat"]:
                label += " distal"
        blocks.append(
            f"<text x='{left - 10}' y='{y + 9}' text-anchor='end' class='axis'>{html.escape(label)}</text>"
        )
        if repeat_row["is_distal_repeat"]:
            blocks.append(
                f"<rect x='{left - 4}' y='{y - 2}' width='{bins_per_unit * cell_w + 8}' height='{cell_h + 4}' fill='none' stroke='#2f7568' stroke-width='2' rx='6'/>"
            )
        for col_idx, pct in enumerate(repeat_row["cells"]):
            x = left + col_idx * cell_w
            blocks.append(
                f"<rect x='{x}' y='{y}' width='{cell_w - 1}' height='{cell_h}' rx='2' fill='{color_for_pct(pct)}'/>"
            )

    ticks = []
    for frac, label in ((0.0, "0 bp"), (0.33, "~1.1 kb"), (0.66, "~2.2 kb"), (1.0, "3.3 kb")):
        x = left + frac * (bins_per_unit - 1) * cell_w
        ticks.append(f"<line x1='{x:.1f}' y1='{top - 6}' x2='{x:.1f}' y2='{height - 26}' class='gridline'/>")
        ticks.append(f"<text x='{x:.1f}' y='{height - 8}' text-anchor='middle' class='axis'>{html.escape(label)}</text>")

    legend_x = left
    legend_y = height - 34
    legend = [
        f"<text x='{legend_x}' y='{legend_y - 10}' class='axis'>Mean CpG methylation legend</text>",
        f"<rect x='{legend_x}' y='{legend_y}' width='22' height='10' rx='3' fill='{color_for_pct(0)}'/>",
        f"<text x='{legend_x + 28}' y='{legend_y + 9}' class='axis'>0%</text>",
        f"<rect x='{legend_x + 66}' y='{legend_y}' width='22' height='10' rx='3' fill='{color_for_pct(50)}'/>",
        f"<text x='{legend_x + 94}' y='{legend_y + 9}' class='axis'>50%</text>",
        f"<rect x='{legend_x + 142}' y='{legend_y}' width='22' height='10' rx='3' fill='{color_for_pct(100)}'/>",
        f"<text x='{legend_x + 170}' y='{legend_y + 9}' class='axis'>100%</text>",
        f"<rect x='{legend_x + 238}' y='{legend_y - 1}' width='22' height='12' rx='4' fill='none' stroke='#2f7568' stroke-width='2'/>",
        f"<text x='{legend_x + 266}' y='{legend_y + 9}' class='axis'>Distal DUX4-containing repeat</text>",
        f"<text x='{legend_x}' y='{legend_y + 26}' class='axis'>X-axis bins tile one 3.3 kb D4Z4 repeat from proximal to distal sequence.</text>",
    ]

    return (
        f"<svg viewBox='0 0 {width} {height}' class='methyl-svg' role='img' aria-label='{html.escape(title)}'>"
        f"<text x='{left}' y='22' class='title'>{html.escape(title)}</text>"
        + "".join(ticks)
        + "".join(blocks)
        + "".join(legend)
        + "</svg>"
    )


def repeat_bin_strip_svg(repeat_rows, title):
    if not repeat_rows:
        return "<p class='muted'>No repeat-bin summary was available.</p>"

    width = max(520, len(repeat_rows) * 18 + 80)
    height = 224
    left = 42
    top = 34
    plot_h = 82
    bar_w = 14

    def color_for_pct(pct):
        if pct is None:
            return "#efe8df"
        pct = max(0.0, min(100.0, pct))
        red = int(244 - pct * 0.28)
        green = int(231 - pct * 1.18)
        blue = int(219 - pct * 0.92)
        return f"rgb({red},{green},{blue})"

    max_cov = max((row["coverage"] for row in repeat_rows), default=1.0) or 1.0
    parts = []
    y_ticks = []
    for frac, label in ((0.0, "0"), (0.5, f"{max_cov / 2:.0f}"), (1.0, f"{max_cov:.0f}")):
        y = top + plot_h - frac * plot_h
        y_ticks.append(f"<line x1='{left - 6}' y1='{y:.1f}' x2='{width - 18}' y2='{y:.1f}' class='gridline'/>")
        y_ticks.append(f"<text x='{left - 12}' y='{y + 4:.1f}' text-anchor='end' class='axis'>{label}</text>")
    for idx, row in enumerate(repeat_rows):
        x = left + idx * 18
        bar_h = (row["coverage"] / max_cov) * plot_h if max_cov else 0.0
        y = top + plot_h - bar_h
        parts.append(
            f"<rect x='{x}' y='{y:.1f}' width='{bar_w}' height='{bar_h:.1f}' rx='4' fill='{color_for_pct(row['pct'])}'/>"
        )
        if row["is_distal_repeat"]:
            parts.append(
                f"<rect x='{x - 2}' y='{top - 6}' width='{bar_w + 4}' height='{plot_h + 12}' fill='none' stroke='#2f7568' stroke-width='2' rx='6'/>"
            )
        if idx == 0 or idx == len(repeat_rows) - 1 or row["is_distal_repeat"]:
            label = "distal" if row["is_distal_repeat"] else str(row["ordinal"])
            parts.append(f"<text x='{x + bar_w / 2}' y='{top + plot_h + 18}' text-anchor='middle' class='axis'>{html.escape(label)}</text>")

    return (
        f"<svg viewBox='0 0 {width} {height}' class='methyl-svg' role='img' aria-label='{html.escape(title)}'>"
        f"<text x='{left}' y='22' class='title'>{html.escape(title)}</text>"
        + "".join(y_ticks)
        + f"<line x1='{left - 6}' y1='{top + plot_h:.1f}' x2='{width - 18}' y2='{top + plot_h:.1f}' class='axisline'/>"
        + "".join(parts)
        + f"<rect x='{left}' y='{height - 46}' width='16' height='12' rx='4' fill='{color_for_pct(80)}'/>"
        + f"<text x='{left + 22}' y='{height - 36}' class='axis'>Bar fill = mean CpG methylation</text>"
        + f"<rect x='{left + 200}' y='{height - 47}' width='16' height='14' rx='4' fill='none' stroke='#2f7568' stroke-width='2'/>"
        + f"<text x='{left + 222}' y='{height - 36}' class='axis'>Outlined bar = distal repeat</text>"
        + f"<text x='{left}' y='{height - 16}' class='axis'>X-axis orders observed repeat bins across the permissive chr4 locus. Bar height reflects summed modkit coverage.</text>"
        + f"<text x='16' y='{top + plot_h / 2:.1f}' text-anchor='middle' class='axis' transform='rotate(-90 16 {top + plot_h / 2:.1f})'>Coverage-weighted support</text>"
        + "</svg>"
    )


def repeat_methyl_track_svg(repeat_rows, title):
    if not repeat_rows:
        return "<p class='muted'>No repeat-level methylation track was available.</p>"

    left = 56
    top = 42
    box_w = 34
    gap = 8
    width = max(720, left * 2 + len(repeat_rows) * (box_w + gap))
    height = 222

    def color_for_pct(pct):
        if pct is None:
            return "#efe8df"
        pct = max(0.0, min(100.0, pct))
        red = int(244 - pct * 0.28)
        green = int(231 - pct * 1.18)
        blue = int(219 - pct * 0.92)
        return f"rgb({red},{green},{blue})"

    parts = [
        f"<text x='{left}' y='22' class='title'>{html.escape(title)}</text>",
        f"<text x='{left}' y='{height - 18}' class='axis'>Each box is one repeat bin observed in at least one supporting read set. Box color is mean CpG methylation and the green outline marks the distal repeat.</text>",
    ]

    for idx, row in enumerate(repeat_rows):
        x = left + idx * (box_w + gap)
        pct = row["pct"]
        label = row.get("display_label") or ("distal" if row["is_distal_repeat"] else f"RU {row['ordinal']}")
        coverage = row["coverage"]
        parts.append(
            f"<rect x='{x}' y='{top}' width='{box_w}' height='74' rx='10' fill='{color_for_pct(pct)}' stroke='rgba(69,54,42,0.12)' stroke-width='1'/>"
        )
        if row["is_distal_repeat"]:
            parts.append(
                f"<rect x='{x - 3}' y='{top - 3}' width='{box_w + 6}' height='80' rx='13' fill='none' stroke='#2f7568' stroke-width='2'/>"
            )
        parts.append(f"<text x='{x + box_w / 2}' y='{top + 20}' text-anchor='middle' class='track-label'>{html.escape(label)}</text>")
        parts.append(
            f"<text x='{x + box_w / 2}' y='{top + 42}' text-anchor='middle' class='axis'>{html.escape('NA' if pct is None else f'{pct:.1f}%')}</text>"
        )
        parts.append(
            f"<text x='{x + box_w / 2}' y='{top + 98}' text-anchor='middle' class='axis'>cov {coverage:.0f}</text>"
        )

    parts.extend(
        [
            f"<rect x='{left}' y='{height - 50}' width='16' height='12' rx='4' fill='{color_for_pct(80)}'/>",
            f"<text x='{left + 22}' y='{height - 40}' class='axis'>Color scale = mean CpG methylation</text>",
            f"<rect x='{left + 210}' y='{height - 51}' width='16' height='14' rx='4' fill='none' stroke='#2f7568' stroke-width='2'/>",
            f"<text x='{left + 232}' y='{height - 40}' class='axis'>Outlined box = distal repeat</text>",
        ]
    )

    return (
        f"<svg viewBox='0 0 {width} {height}' class='methyl-svg' role='img' aria-label='{html.escape(title)}'>"
        + "".join(parts)
        + "</svg>"
    )


def collective_tiling_assessment(observed_repeat_bins, threshold):
    if observed_repeat_bins is None:
        return (
            "unknown",
            "No repeat-resolved methylation tiling was available for a collective copy-number estimate.",
        )
    if observed_repeat_bins <= threshold:
        return (
            "possible_short",
            f"Observed 4qA-supporting methylation tiling covered about {observed_repeat_bins} repeat bins, which is compatible with a short array but is not a direct single-read call.",
        )
    return (
        "not_short",
        f"Observed 4qA-supporting methylation tiles covered about {observed_repeat_bins} repeat bins, which argues against a simple <= {threshold} repeat-unit contraction in the reads that were observed.",
    )


def estimate_repeat_units_collective(repeat_rows, has_proximal, has_distal):
    filtered = [row for row in repeat_rows if row["repeat_idx"] <= 0]
    if not filtered:
        return {
            "estimate_ru": None,
            "observed_bins": 0,
            "gap_bins": None,
            "method": "unavailable",
            "confidence": "none",
            "note": "No repeat-folded methylation bins were available for collective size estimation.",
        }

    min_idx = min(row["repeat_idx"] for row in filtered)
    max_idx = max(row["repeat_idx"] for row in filtered)
    span_ru = (max_idx - min_idx) + 1
    observed_bins = len(filtered)
    gap_bins = max(0, span_ru - observed_bins)

    if has_proximal and has_distal:
        confidence = "low" if gap_bins > 2 else "medium"
        note = (
            f"Permissive 4qA reads tile repeat bins {min_idx} through {max_idx} relative to the distal DUX4-containing repeat, "
            f"supporting an inferred collective array length of about {span_ru} repeat units. "
            f"{observed_bins} bins were directly observed and {gap_bins} bins are interpolated across gaps, so this is an inferred collective size estimate rather than a spanning-read call."
        )
        return {
            "estimate_ru": span_ru,
            "observed_bins": observed_bins,
            "gap_bins": gap_bins,
            "method": "collective_tiling_distal_anchor",
            "confidence": confidence,
            "note": note,
        }

    side = "proximal" if has_proximal else "distal" if has_distal else "internal"
    note = (
        f"Only {side} 4qA support was available across repeat-folded methylation bins, so the current data support a partial lower-bound view but not a balanced collective size estimate."
    )
    return {
        "estimate_ru": None,
        "observed_bins": observed_bins,
        "gap_bins": gap_bins,
        "method": f"{side}_only_partial_tiling",
        "confidence": "low",
        "note": note,
    }


def locus_annotation_style(label):
    clean = clean_cell(label)
    if any(token in clean for token in ("D4Z4", "DUX4")):
        return {
            "fill": "rgba(38,120,108,0.18)",
            "text": "#245d54",
            "legend": "D4Z4 / DUX4 interval",
        }
    if "pLAM" in clean:
        return {
            "fill": "rgba(188,84,52,0.16)",
            "text": "#8f3518",
            "legend": "pLAM / distal haplotype marker",
        }
    return {
        "fill": "rgba(181,138,43,0.18)",
        "text": "#7d5a07",
        "legend": "Anchor / marker interval",
    }


def pretty_annotation_label(label):
    return clean_cell(label).replace("_", " ")


def format_locus_tick(pos):
    return f"{pos / 1_000_000:.3f} Mb"


def methyl_profile_svg(rows, annotations, title):
    if not rows:
        return "<p class='muted'>No locus methylation pileup was available.</p>"
    start = min(row["start"] for row in rows)
    end = max(row["end"] for row in rows)
    bins = summarize_windows(rows, start, end, max(1, (end - start) // 48))
    width = 1280
    right = 34
    left = 78
    base_top = 38
    bottom = 92
    plot_h = 210
    chrom = rows[0]["chrom"]
    plot_w = width - left - right
    anno_specs = []
    for anno in annotations:
        if anno["chrom"] != chrom:
            continue
        if anno["end"] < start or anno["start"] > end:
            continue
        x1 = left + ((max(start, anno["start"]) - start) / max(1, end - start)) * plot_w
        x2 = left + ((min(end, anno["end"]) - start) / max(1, end - start)) * plot_w
        style = locus_annotation_style(anno["label"])
        label = pretty_annotation_label(anno["label"])
        label_width = max(72.0, min(180.0, len(label) * 5.8))
        anno_specs.append(
            {
                "x1": x1,
                "x2": x2,
                "mid": (x1 + x2) / 2.0,
                "label": label,
                "label_width": label_width,
                "fill": style["fill"],
                "text": style["text"],
                "legend": style["legend"],
            }
        )

    anno_specs.sort(key=lambda item: item["x1"])
    lane_ends = []
    for spec in anno_specs:
        center = min(width - right - spec["label_width"] / 2.0, max(left + spec["label_width"] / 2.0, spec["mid"]))
        label_left = center - spec["label_width"] / 2.0
        label_right = center + spec["label_width"] / 2.0
        lane = None
        for idx, last_right in enumerate(lane_ends):
            if label_left > last_right + 10:
                lane = idx
                lane_ends[idx] = label_right
                break
        if lane is None:
            lane = len(lane_ends)
            lane_ends.append(label_right)
        spec["label_x"] = center
        spec["lane"] = lane

    top = base_top + len(lane_ends) * 16
    height = top + plot_h + bottom
    left = 78
    plot_w = width - left - right
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
    legend_map = {}
    for spec in anno_specs:
        anno_html.append(
            f"<rect x='{spec['x1']:.1f}' y='{top - 4}' width='{max(2.0, spec['x2'] - spec['x1']):.1f}' height='{plot_h + 8}' fill='{spec['fill']}'/>"
        )
        label_y = base_top + spec["lane"] * 16
        anno_html.append(
            f"<line x1='{spec['mid']:.1f}' y1='{label_y + 3:.1f}' x2='{spec['mid']:.1f}' y2='{top - 6:.1f}' stroke='{spec['text']}' stroke-width='1.2' stroke-dasharray='2 2'/>"
        )
        anno_html.append(
            f"<text x='{spec['label_x']:.1f}' y='{label_y:.1f}' text-anchor='middle' class='anno'>{html.escape(spec['label'])}</text>"
        )
        legend_map[spec["legend"]] = spec["fill"]
    y_ticks = []
    for pct in (0, 25, 50, 75, 100):
        y = top + plot_h - (pct / 100.0) * plot_h
        y_ticks.append(f"<line x1='{left}' y1='{y:.1f}' x2='{width - right}' y2='{y:.1f}' class='gridline'/>")
        y_ticks.append(f"<text x='{left - 10}' y='{y + 4:.1f}' text-anchor='end' class='axis'>{pct}%</text>")
    x_ticks = []
    for frac, pos in ((0.0, start), (0.25, start + (end - start) * 0.25), (0.5, start + (end - start) * 0.5), (0.75, start + (end - start) * 0.75), (1.0, end)):
        x = left + frac * plot_w
        x_ticks.append(f"<line x1='{x:.1f}' y1='{top + plot_h:.1f}' x2='{x:.1f}' y2='{top + plot_h + 6:.1f}' class='axisline'/>")
        x_ticks.append(f"<text x='{x:.1f}' y='{top + plot_h + 20:.1f}' text-anchor='middle' class='axis'>{html.escape(format_locus_tick(pos))}</text>")

    legend_parts = [
        f"<line x1='{left}' y1='{height - 54}' x2='{left + 24}' y2='{height - 54}' stroke='#bc5434' stroke-width='3'/>",
        f"<text x='{left + 32}' y='{height - 50}' class='axis'>Mean CpG methylation profile</text>",
        f"<rect x='{left}' y='{height - 38}' width='16' height='12' rx='4' fill='rgba(188,84,52,0.14)'/>",
        f"<text x='{left + 24}' y='{height - 28}' class='axis'>Coverage-weighted methylation area</text>",
    ]
    legend_x = left + 280
    for idx, (legend_label, fill) in enumerate(sorted(legend_map.items())):
        legend_parts.append(f"<rect x='{legend_x}' y='{height - 38}' width='16' height='12' rx='4' fill='{fill}'/>")
        legend_parts.append(f"<text x='{legend_x + 24}' y='{height - 28}' class='axis'>{html.escape(legend_label)}</text>")
        legend_x += 230
    return (
        f"<svg viewBox='0 0 {width} {height}' class='methyl-svg' role='img' aria-label='{html.escape(title)}'>"
        f"<text x='{left}' y='22' class='title'>{html.escape(title)}</text>"
        + "".join(anno_html)
        + "".join(y_ticks)
        + "".join(x_ticks)
        + "".join(bars)
        + f"<polyline fill='none' stroke='#bc5434' stroke-width='3' points='{' '.join(pts)}'/>"
        + f"<line x1='{left}' y1='{top + plot_h:.1f}' x2='{width - right}' y2='{top + plot_h:.1f}' class='axisline'/>"
        + f"<text x='{width / 2:.1f}' y='{height - 8}' text-anchor='middle' class='axis'>Locus coordinate ({html.escape(chrom)}, T2T-CHM13)</text>"
        + f"<text x='18' y='{top + plot_h / 2:.1f}' text-anchor='middle' class='axis' transform='rotate(-90 18 {top + plot_h / 2:.1f})'>Mean CpG methylation</text>"
        + "".join(legend_parts)
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
        contraction_threshold,
        locus_bed,
        methylation_bed,
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
    rows_4qb_all = data["4qB_all"]
    rows_4qb_complete = data["4qB_complete"]
    rows_chr10_all = data["chr10_all"]
    rows_chr10_complete = data["chr10_complete"]
    ru_4qa_all = ru_values(rows_4qa_all)
    ru_4qa_complete = ru_values(rows_4qa_complete)
    ru_4qb_all = ru_values(rows_4qb_all)
    ru_4qb_complete = ru_values(rows_4qb_complete)
    ru_chr10_all = ru_values(rows_chr10_all)
    shortest_4qa = min(ru_4qa_complete) if ru_4qa_complete else None
    best_partial_4qa = max(ru_4qa_all) if ru_4qa_all else None

    contraction_status, contraction_note = contraction_assessment(rows_4qa_complete, contraction_threshold)
    contracted = contraction_status == "supported"

    flagstat_headline = "NA"
    if os.path.exists(flagstat_txt):
        with open(flagstat_txt, "r", encoding="utf-8") as fh:
            flagstat_headline = clean_cell(fh.readline())

    coverage_rows = coverage_rows_from_tsv(coverage_tsv)
    coverage_focus_chr4 = [
        row for row in coverage_rows if row.get("#chrom") == "chr4" and row.get("label") in ("D4Z4_4q35", "DUX4_gene_body_chr4", "pLAM_4qA", "DUX4_end_chr4")
    ]
    coverage_focus_chr10 = [
        row
        for row in coverage_rows
        if row.get("#chrom") == "chr10"
        and row.get("label") in ("DUX4_homologue_chr10", "CTRL_D4Z4_chr10", "CLUHP4_201_exon1_chr10", "D4F104S1_10q", "pLAM_10qA")
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
    methyl_annotations = parse_bed_annotations(methylation_bed)
    repeat_anchor_start = annotation_start(methyl_annotations, "distal_RU_gene_body")
    repeat_rows = repeat_unit_matrix(methyl_records, repeat_anchor_start, REPEAT_BP, bins_per_unit=60, min_occupied_cols=8)
    chr10_repeat_anchor_start = annotation_start(methyl_annotations, "distal_RU_gene_body_chr10")
    methyl_records_chr4_all = methyl_records_for_subset(methyl_by_subset, methylation_dir, "D4Z4_chr4_only") or methyl_records
    methyl_records_chr10_all = methyl_records_for_subset(methyl_by_subset, methylation_dir, "D4Z4_chr10_only")
    repeat_rows_chr4_all = repeat_unit_matrix(methyl_records_chr4_all, repeat_anchor_start, REPEAT_BP, bins_per_unit=60, min_occupied_cols=8)
    repeat_rows_chr10_all = repeat_unit_matrix(methyl_records_chr10_all, chr10_repeat_anchor_start, REPEAT_BP, bins_per_unit=60, min_occupied_cols=8)
    combined_repeat_rows = with_repeat_labels(repeat_rows_chr4_all, "4q") + with_repeat_labels(repeat_rows_chr10_all, "10q")
    observed_repeat_bins = len(repeat_rows) if repeat_rows else None
    collective_status, collective_note = collective_tiling_assessment(observed_repeat_bins, contraction_threshold)
    distal_repeat_pct = next((row["pct"] for row in repeat_rows if row["is_distal_repeat"]), None)

    if contraction_status == "supported":
        headline = f"Complete 4qA read support for an allele at or below {contraction_threshold} repeat units was observed."
    elif contraction_status == "not_supported":
        headline = f"Complete 4qA reads were present, but none supported a contraction at {contraction_threshold} repeat units or fewer."
    elif counts["4qA_all"] > 0:
        if collective_status == "not_short" and observed_repeat_bins is not None:
            headline = f"Partial 4qA-supporting reads were present without a complete spanning read, but collective 4qA tiling covered about {observed_repeat_bins} repeat bins and does not look contraction-like."
        else:
            headline = "Partial 4qA-supporting reads were present, but no complete 4qA-spanning read was classified."
    else:
        headline = "No 4qA-supporting reads were classified in this sample."

    pas_counter = Counter(
        clean_cell(row.get("PAS.type", "NA")) or "NA" for row in rows_4qa_all
    )
    proximal_rows = [row for row in rows_4qa_all if has_proximal_anchor(parse_repeat_tokens(row.get("repeat_sequence")))]
    distal_rows = [row for row in rows_4qa_all if has_distal_anchor(parse_repeat_tokens(row.get("repeat_sequence")))]
    best_proximal_ru = max((to_int(row.get("RU_count"), -1) for row in proximal_rows), default=None)
    best_distal_ru = max((to_int(row.get("RU_count"), -1) for row in distal_rows), default=None)
    collective_length = estimate_repeat_units_collective(repeat_rows, bool(proximal_rows), bool(distal_rows))
    allele_4qa = summarize_allele_candidate("4qA-permissive allele", rows_4qa_all, rows_4qa_complete, collective_length)
    allele_4qb = summarize_allele_candidate("4qB-permissive allele", rows_4qb_all, rows_4qb_complete)
    allele_chr10 = summarize_allele_candidate("chr10 pooled homologous signal", rows_chr10_all, rows_chr10_complete)
    chr10_candidates = split_pooled_repeat_candidates("chr10", rows_chr10_all, ("10qA_D4F104S1", "CLUHP4"), ("10qA_pLAM", "DUX4_end"))
    while len(chr10_candidates) < 2:
        chr10_candidates.append(
            {
                "label": f"chr10 candidate {len(chr10_candidates) + 1}",
                "estimate_ru": None,
                "display_ru": "NA",
                "method": "unresolved",
                "confidence": "none",
                "support_reads": 0,
                "best_partial_ru": None,
                "proximal_support_reads": 0,
                "distal_support_reads": 0,
                "note": "No separate chr10 allele-like cluster could be resolved from the pooled partial reads.",
            }
        )
    allele_4qa["display_ru"] = normalize_size_display(allele_4qa.get("display_ru"), f">={best_partial_4qa}" if best_partial_4qa is not None else None)
    allele_4qb["display_ru"] = normalize_size_display(allele_4qb.get("display_ru"), f">={max(ru_4qb_all)}" if ru_4qb_all else None)
    allele_chr10["display_ru"] = normalize_size_display(allele_chr10.get("display_ru"), f">={max(ru_chr10_all)}" if ru_chr10_all else None)
    chr10_candidates[0]["display_ru"] = normalize_size_display(chr10_candidates[0].get("display_ru"), f">={chr10_candidates[0].get('best_partial_ru')}" if chr10_candidates[0].get("best_partial_ru") is not None else None)
    chr10_candidates[1]["display_ru"] = normalize_size_display(chr10_candidates[1].get("display_ru"), f">={chr10_candidates[1].get('best_partial_ru')}" if chr10_candidates[1].get("best_partial_ru") is not None else None)
    chr4_total_lb = sum(v for v in (size_call_lower_bound(allele_4qa), size_call_lower_bound(allele_4qb)) if v is not None)
    chr10_total_lb = sum(v for v in (size_call_lower_bound(chr10_candidates[0]), size_call_lower_bound(chr10_candidates[1])) if v is not None)
    genome_total_lb = None if chr10_total_lb is None else chr4_total_lb + chr10_total_lb

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
        "contraction_assessment": contraction_status,
        "best_partial_4qA_ru": best_partial_4qa if best_partial_4qa is not None else "NA",
        "best_proximal_4qA_ru": best_proximal_ru if best_proximal_ru is not None else "NA",
        "best_distal_4qA_ru": best_distal_ru if best_distal_ru is not None else "NA",
        "locus_4qA_methyl_mean_pct": f"{methyl_mean:.3f}" if methyl_mean is not None else "NA",
        "observed_4qA_repeat_bins": observed_repeat_bins if observed_repeat_bins is not None else "NA",
        "collective_4qA_tiling_assessment": collective_status,
        "distal_repeat_methyl_mean_pct": f"{distal_repeat_pct:.3f}" if distal_repeat_pct is not None else "NA",
        "estimated_4qA_repeat_units": collective_length["estimate_ru"] if collective_length["estimate_ru"] is not None else "NA",
        "estimated_4qA_repeat_units_method": collective_length["method"],
        "estimated_4qA_repeat_units_confidence": collective_length["confidence"],
        "estimated_4qA_repeat_units_note": collective_length["note"],
        "chr4_allele_1_label": "4qA",
        "chr4_allele_1_size_call": allele_4qa["display_ru"],
        "chr4_allele_1_method": allele_4qa["method"],
        "chr4_allele_1_confidence": allele_4qa["confidence"],
        "chr4_allele_1_support_reads": allele_4qa["support_reads"],
        "chr4_allele_2_label": "4qB",
        "chr4_allele_2_size_call": allele_4qb["display_ru"],
        "chr4_allele_2_method": allele_4qb["method"],
        "chr4_allele_2_confidence": allele_4qb["confidence"],
        "chr4_allele_2_support_reads": allele_4qb["support_reads"],
        "chr10_pooled_size_call": allele_chr10["display_ru"],
        "chr10_pooled_method": allele_chr10["method"],
        "chr10_pooled_confidence": allele_chr10["confidence"],
        "chr10_pooled_support_reads": allele_chr10["support_reads"],
        "chr10_allele_1_size_call": chr10_candidates[0]["display_ru"],
        "chr10_allele_1_method": chr10_candidates[0]["method"],
        "chr10_allele_1_confidence": chr10_candidates[0]["confidence"],
        "chr10_allele_1_support_reads": chr10_candidates[0]["support_reads"],
        "chr10_allele_2_size_call": chr10_candidates[1]["display_ru"],
        "chr10_allele_2_method": chr10_candidates[1]["method"],
        "chr10_allele_2_confidence": chr10_candidates[1]["confidence"],
        "chr10_allele_2_support_reads": chr10_candidates[1]["support_reads"],
        "chr4_total_d4z4_units_lower_bound": chr4_total_lb if chr4_total_lb else "NA",
        "chr10_total_d4z4_units_lower_bound": chr10_total_lb if chr10_total_lb is not None else "NA",
        "genome_total_d4z4_units_lower_bound": genome_total_lb if genome_total_lb is not None else "NA",
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

    coverage_table_chr4 = table_html(
        ["Label", "Reads", "Covered bp", "Breadth", "Mean depth"],
        [
            [
                html.escape(row.get("label", "")),
                html.escape(row.get("numreads", "")),
                html.escape(row.get("covbases", "")),
                html.escape(row.get("coverage", "")),
                html.escape(row.get("meandepth", "")),
            ]
            for row in coverage_focus_chr4
        ],
        "compact",
    )

    coverage_table_chr10 = table_html(
        ["Label", "Reads", "Covered bp", "Breadth", "Mean depth"],
        [
            [
                html.escape(row.get("label", "")),
                html.escape(row.get("numreads", "")),
                html.escape(row.get("covbases", "")),
                html.escape(row.get("coverage", "")),
                html.escape(row.get("meandepth", "")),
            ]
            for row in coverage_focus_chr10
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
    .metric-value {{
      font-size: clamp(1.15rem, 2.1vw, 2rem);
      margin-top: 8px;
      font-weight: 700;
      line-height: 1.15;
      overflow-wrap: anywhere;
      word-break: break-word;
      hyphens: auto;
    }}
    .metric-note {{
      margin-top: 6px;
      font-size: 0.92rem;
      color: var(--slate);
      line-height: 1.35;
      overflow-wrap: anywhere;
      word-break: break-word;
      hyphens: auto;
    }}
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
    .track-label {{ font-size: 10px; fill: #45362a; font-weight: 700; }}
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
      <div class="hero-note" style="border-left-color: var(--sage); background: rgba(47,117,104,0.08); margin-top: 10px;">
        <strong>Collective repeat tiling:</strong> {html.escape(collective_note)}
      </div>
      <div class="hero-note" style="border-left-color: var(--gold); background: rgba(182,136,43,0.08); margin-top: 10px;">
        <strong>Estimated 4qA repeat length:</strong> {html.escape(collective_length["note"])}
      </div>
      <div class="hero-note" style="border-left-color: #7058a3; background: rgba(112,88,163,0.08); margin-top: 10px;">
        <strong>Diploid chr4 allele view:</strong> 4qA-permissive candidate {html.escape(allele_4qa["display_ru"])} RU ({html.escape(allele_4qa["method"])}), 4qB-permissive candidate {html.escape(allele_4qb["display_ru"])} RU ({html.escape(allele_4qb["method"])}). This view does not rely on haplotags.
      </div>
      <div class="hero-note" style="border-left-color: #4f6fa8; background: rgba(79,111,168,0.08); margin-top: 10px;">
        <strong>Genome-wide D4Z4 lower bound:</strong> chr4 candidates sum to at least {html.escape(str(chr4_total_lb) if chr4_total_lb else 'NA')} repeat units, chr10 heuristic candidates support at least {html.escape(chr10_candidates[0]["display_ru"])} and {html.escape(chr10_candidates[1]["display_ru"])} repeat units, for a conservative genome-wide lower bound of {html.escape(str(genome_total_lb) if genome_total_lb is not None else 'NA')} D4Z4 units.
      </div>
      <div class="metric-grid">
        {metric_card("4qA-supporting reads", counts["4qA_all"], "All permissive 4qA-classified reads", "good" if counts["4qA_all"] else "soft")}
        {metric_card("4qB-supporting reads", counts["4qB_all"], "All permissive 4qB-classified reads", "good" if counts["4qB_all"] else "soft")}
        {metric_card("chr10-supporting reads", counts["chr10_all"], "All chr10-homologous reads", "good" if counts["chr10_all"] else "soft")}
        {metric_card("Best observed read span", best_partial_4qa if best_partial_4qa is not None else "NA", "Largest partial or complete RU count on one 4qA read", "soft")}
        {metric_card("4qA allele size call", allele_4qa["display_ru"], f"{allele_4qa['confidence']} confidence {allele_4qa['method']}", "good" if allele_4qa["confidence"] in ("high", "medium") else "soft")}
        {metric_card("4qB allele size call", allele_4qb["display_ru"], f"{allele_4qb['confidence']} confidence {allele_4qb['method']}", "good" if allele_4qb["confidence"] in ("high", "medium") else "soft")}
        {metric_card("chr10 pooled size call", allele_chr10["display_ru"], f"{allele_chr10['confidence']} confidence {allele_chr10['method']}", "good" if allele_chr10["confidence"] in ("high", "medium") else "soft")}
        {metric_card("Collective tiled repeat bins", observed_repeat_bins if observed_repeat_bins is not None else "NA", "Observed repeat-sized bins across permissive 4qA methylation tiling", "good" if collective_status == "not_short" else "warn" if collective_status == "possible_short" else "soft")}
        {metric_card("4qA PAS-positive reads", pas_counter.get("4qA_PAS", 0), "ATTAAA-compatible pLAM support", "good" if pas_counter.get("4qA_PAS", 0) else "soft")}
        {metric_card("4qA locus methylation", f"{methyl_mean:.1f}%" if methyl_mean is not None else methyl_4qa_all.get("mean_percent_modified", "NA"), "Weighted mean across modkit pileup positions", "good" if methyl_mean is not None else "soft")}
        {metric_card("Distal repeat methylation", f"{distal_repeat_pct:.1f}%" if distal_repeat_pct is not None else "NA", "Folded methylation mean on the distal DUX4-containing repeat", "good" if distal_repeat_pct is not None else "soft")}
        {metric_card("Complete-read contraction", "yes" if contracted else "no" if contraction_status == "not_supported" else "inconclusive", f"Threshold: <= {contraction_threshold} RU", "good" if contracted else "warn" if contraction_status == "not_supported" else "soft")}
      </div>
      <p class="footnote"><strong>Locus BAM flagstat:</strong> {html.escape(flagstat_headline)}</p>
    </section>

    <section class="section">
      <h2>Diploid chr4 Allele Candidates</h2>
      {table_html(
          ["Candidate", "Distal haplotype", "Size call", "Method", "Confidence", "Support reads", "Interpretation"],
          [
              [
                  "Allele candidate 1",
                  "4qA-permissive",
                  html.escape(allele_4qa["display_ru"]),
                  html.escape(allele_4qa["method"]),
                  html.escape(allele_4qa["confidence"]),
                  html.escape(str(allele_4qa["support_reads"])),
                  html.escape(allele_4qa["note"]),
              ],
              [
                  "Allele candidate 2",
                  "4qB-permissive",
                  html.escape(allele_4qb["display_ru"]),
                  html.escape(allele_4qb["method"]),
                  html.escape(allele_4qb["confidence"]),
                  html.escape(str(allele_4qb["support_reads"])),
                  html.escape(allele_4qb["note"]),
              ],
          ],
      )}
      <p class="footnote">This diploid summary deliberately avoids haplotags. Instead it treats 4qA- and 4qB-permissive read groups as chr4 allele candidates and reports either a direct size, a collective estimate, or a conservative lower bound depending on the evidence available.</p>
    </section>

    <section class="section">
      <h2>Genome-Wide D4Z4 Burden</h2>
      {table_html(
          ["Component", "Size call", "Method", "Support reads", "Interpretation"],
          [
              [
                  "chr4 4qA-permissive candidate",
                  html.escape(allele_4qa["display_ru"]),
                  html.escape(allele_4qa["method"]),
                  html.escape(str(allele_4qa["support_reads"])),
                  html.escape(allele_4qa["note"]),
              ],
              [
                  "chr4 4qB-permissive candidate",
                  html.escape(allele_4qb["display_ru"]),
                  html.escape(allele_4qb["method"]),
                  html.escape(str(allele_4qb["support_reads"])),
                  html.escape(allele_4qb["note"]),
              ],
              [
                  "chr10 pooled homologous signal",
                  html.escape(allele_chr10["display_ru"]),
                  html.escape(allele_chr10["method"]),
                  html.escape(str(allele_chr10["support_reads"])),
                  html.escape(allele_chr10["note"]),
              ],
              [
                  "chr10 heuristic candidate 1",
                  html.escape(chr10_candidates[0]["display_ru"]),
                  html.escape(chr10_candidates[0]["method"]),
                  html.escape(str(chr10_candidates[0]["support_reads"])),
                  html.escape(chr10_candidates[0]["note"]),
              ],
              [
                  "chr10 heuristic candidate 2",
                  html.escape(chr10_candidates[1]["display_ru"]),
                  html.escape(chr10_candidates[1]["method"]),
                  html.escape(str(chr10_candidates[1]["support_reads"])),
                  html.escape(chr10_candidates[1]["note"]),
              ],
          ],
      )}
      <div class="metric-grid" style="margin-top:16px;">
        {metric_card("chr4 total lower bound", chr4_total_lb if chr4_total_lb else "NA", "Sum of the two chr4 allele-candidate size calls or lower bounds", "soft")}
        {metric_card("chr10 candidate 1", chr10_candidates[0]["display_ru"], "Heuristic split of pooled chr10 partial reads", "soft")}
        {metric_card("chr10 candidate 2", chr10_candidates[1]["display_ru"], "Heuristic split of pooled chr10 partial reads", "soft")}
        {metric_card("Genome total lower bound", genome_total_lb if genome_total_lb is not None else "NA", "Conservative lower bound from resolved chr4 candidates plus chr10 pooled signal", "soft")}
      </div>
      <p class="footnote">The T2T-CHM13 reference contains D4Z4-bearing loci on both chr4 and chr10, so reads can align and be classified there. The chr10 candidate split shown here is a low-confidence heuristic decomposition of pooled partial reads, not a Kivvi-style allele assembly.</p>
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
        {metric_card("Collective tiling view", "not contraction-like" if collective_status == "not_short" else "possibly short" if collective_status == "possible_short" else "NA", "Inference from all observed 4qA methylation tiles, not a single spanning read", "good" if collective_status == "not_short" else "warn" if collective_status == "possible_short" else "soft")}
        {metric_card("Inferred array length", collective_length["estimate_ru"] if collective_length["estimate_ru"] is not None else "NA", f"Method: {collective_length['method']}", "soft")}
      </div>
      <p class="footnote">Partial 4qA reads do not make a direct molecular call on their own, but the combined 4qA tiling can still show whether the observed read set behaves more like a long or short array.</p>
    </section>

    <section class="section">
      <h2>Repeat-Resolved Methylation</h2>
      <div class="dual-grid">
        <div>
          {repeat_pileup_svg(repeat_rows, 60, "Folded 3.3 kb D4Z4 methylation pileup")}
          <p class="footnote">Each row is one observed repeat-sized bin across the permissive 4qA read set. Columns represent position within a single 3.3 kb D4Z4 unit, and the cell color reflects mean methylation across modkit pileup calls folded onto that motif.</p>
        </div>
        <div>
          {repeat_bin_strip_svg(repeat_rows, "Observed repeat bins across the permissive 4qA locus")}
          <p class="footnote">This summary uses the distal DUX4-containing repeat as the 3.3 kb folding anchor. It does not require a single read to span the whole array; instead it summarizes how the observed 4qA-supporting reads tile across repeat-sized bins.</p>
        </div>
      </div>
      <div style="margin-top: 18px;">
        {repeat_methyl_track_svg(repeat_rows, "Methylation track across all repeat bins observed in 4qA-supporting reads")}
        <p class="footnote">This is the most literal repeat-bin methylation summary available from the current run: every box corresponds to one repeat-sized bin recovered from the `4qA_all` modkit pileup. It does not yet include 4qB or chr10 reads, because the current methylation branch only profiles `4qA_all`, `4qA_complete`, and `chimeric` subsets.</p>
      </div>
      <div class="dual-grid" style="margin-top: 18px;">
        <div>
          {repeat_pileup_svg(repeat_rows_chr4_all, 60, "All D4Z4 repeat bins observed at chr4")}
          <p class="footnote">This panel uses the broader `D4Z4_chr4_only` subset when available, and falls back to `4qA_all` on older runs.</p>
        </div>
        <div>
          {repeat_pileup_svg(repeat_rows_chr10_all, 60, "All D4Z4 repeat bins observed at chr10")}
          <p class="footnote">This panel requires methylation profiling of the `D4Z4_chr10_only` subset. Older runs that only profiled `4qA_all` will show this as unavailable.</p>
        </div>
      </div>
      <div style="margin-top: 18px;">
        {repeat_methyl_track_svg(combined_repeat_rows, "All observed D4Z4 repeat bins regardless of locus")}
        <p class="footnote">This combined track appends the chr4 and chr10 repeat bins into one view. It is intended as a broad locus-level methylation summary rather than an allele-resolved call.</p>
      </div>
    </section>

    <section class="section">
      <h2>4q and 10q Locus Coverage</h2>
      <div class="dual-grid">
        <div>
          {methyl_profile_svg(methyl_records_chr4_all or methyl_records, locus_annotations, "chr4 D4Z4 / DUX4 locus methylation")}
          <p class="footnote">The chr4 panel summarizes methylation across the permissive D4Z4 / DUX4 locus context. Shaded intervals correspond to annotated D4Z4, DUX4, pLAM, and anchor-marker regions from the locus BED.</p>
        </div>
        <div>
          {methyl_profile_svg(methyl_records_chr10_all, locus_annotations, "chr10 D4Z4 homologue locus methylation")}
          <p class="footnote">The chr10 panel uses the broader `D4Z4_chr10_only` subset when available. Older runs that only profiled 4qA-supporting reads will show this as unavailable.</p>
        </div>
      </div>
      <div class="dual-grid" style="margin-top: 18px;">
        <div>
          <h3>chr4 annotated interval summary</h3>
          {coverage_table_chr4}
        </div>
        <div>
          <h3>chr10 annotated interval summary</h3>
          {coverage_table_chr10}
        </div>
      </div>
      <p class="footnote">An IGV-style viewer is the next natural layer here: chr4 permissive locus, chr10 homologue locus, and a single-repeat motif window with the methylation pileup painted on top.</p>
    </section>

    <section class="section">
      <h2>Support Tables</h2>
      <div class="dual-grid">
        <div>
          <h3>Classified subsets</h3>
          {subset_table}
        </div>
        <div>
          <h3>Repeat-span support</h3>
          {table_html(
              ["Signal", "Value"],
              [
                  ["Direct complete-read contraction call", "supported" if contracted else "not supported" if contraction_status == "not_supported" else "inconclusive"],
                  ["4qA-permissive allele size call", allele_4qa["display_ru"]],
                  ["4qB-permissive allele size call", allele_4qb["display_ru"]],
                  ["chr10 pooled size call", allele_chr10["display_ru"]],
                  ["chr10 heuristic candidate 1", chr10_candidates[0]["display_ru"]],
                  ["chr10 heuristic candidate 2", chr10_candidates[1]["display_ru"]],
                  ["chr4 total lower bound", str(chr4_total_lb) if chr4_total_lb else "NA"],
                  ["Genome total lower bound", str(genome_total_lb) if genome_total_lb is not None else "NA"],
                  ["Observed repeat bins in collective 4qA tiling", str(observed_repeat_bins) if observed_repeat_bins is not None else "NA"],
                  ["Estimated collective 4qA length", str(collective_length["estimate_ru"]) if collective_length["estimate_ru"] is not None else "NA"],
                  ["Collective estimate confidence", collective_length["confidence"]],
                  ["Best proximal anchored RU span", str(best_proximal_ru) if best_proximal_ru is not None else "NA"],
                  ["Best distal anchored RU span", str(best_distal_ru) if best_distal_ru is not None else "NA"],
                  ["Weighted 4qA methylation mean", f"{methyl_mean:.1f}%" if methyl_mean is not None else "NA"],
              ],
              "compact",
          )}
        </div>
      </div>
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
