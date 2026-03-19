#!/usr/bin/env python3

import csv
import html
import os
import sys
from collections import Counter, defaultdict


def read_semicolon_csv(path):
    rows = []
    if not os.path.exists(path):
        return rows
    with open(path, "r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter=";")
        for row in reader:
            rows.append(
                {k.strip(): (v.strip() if isinstance(v, str) else v) for k, v in row.items()}
            )
    return rows


def read_tsv(path):
    rows = []
    if not os.path.exists(path):
        return rows
    with open(path, "r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            rows.append(row)
    return rows


def to_int(value, default=0):
    try:
        return int(float(str(value).replace(",", ".")))
    except Exception:
        return default


def truncate(text, width=70):
    text = (text or "").strip()
    return text if len(text) <= width else text[: width - 3] + "..."


def metric_card(label, value, note=""):
    return (
        f"<div class='metric'><div class='metric-label'>{html.escape(label)}</div>"
        f"<div class='metric-value'>{html.escape(str(value))}</div>"
        f"<div class='metric-note'>{html.escape(note)}</div></div>"
    )


def table_html(headers, rows, klass=""):
    if not rows:
        return "<p class='muted'>No rows available.</p>"
    head = "".join(f"<th>{html.escape(h)}</th>" for h in headers)
    body_rows = []
    for row in rows:
        body_rows.append("<tr>" + "".join(f"<td>{cell}</td>" for cell in row) + "</tr>")
    return f"<table class='{klass}'><thead><tr>{head}</tr></thead><tbody>{''.join(body_rows)}</tbody></table>"


def svg_hist(values, color, label):
    if not values:
        return "<p class='muted'>No repeat-unit measurements available.</p>"
    counts = Counter(values)
    xs = sorted(counts)
    max_count = max(counts.values())
    bar_width = 46
    gap = 14
    width = max(420, len(xs) * (bar_width + gap) + 60)
    height = 220
    bars = []
    for idx, x in enumerate(xs):
        count = counts[x]
        bar_h = int((count / max_count) * 130)
        xpos = 40 + idx * (bar_width + gap)
        ypos = 170 - bar_h
        bars.append(
            f"<rect x='{xpos}' y='{ypos}' width='{bar_width}' height='{bar_h}' rx='8' fill='{color}' opacity='0.88' />"
        )
        bars.append(
            f"<text x='{xpos + bar_width / 2}' y='190' text-anchor='middle' class='axis'>{x}</text>"
        )
        bars.append(
            f"<text x='{xpos + bar_width / 2}' y='{max(24, ypos - 8)}' text-anchor='middle' class='count'>{count}</text>"
        )
    return f"""<svg viewBox='0 0 {width} {height}' class='histogram' role='img' aria-label='{html.escape(label)}'>
      <line x1='24' y1='170' x2='{width - 18}' y2='170' stroke='rgba(80,56,36,0.28)' stroke-width='2'/>
      {''.join(bars)}
      <text x='28' y='28' class='title'>{html.escape(label)}</text>
    </svg>"""


def main() -> int:
    (
        sample_id,
        classification_dir,
        subset_dir,
        flagstat_txt,
        coverage_tsv,
        hap_tsv,
        methyl_tsv,
        variant_tsv,
        contraction_threshold,
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

    ru_4qa = [to_int(row.get("RU_count")) for row in data["4qA_complete"] if row.get("RU_count")]
    ru_4qb = [to_int(row.get("RU_count")) for row in data["4qB_complete"] if row.get("RU_count")]
    ru_10q = [to_int(row.get("RU_count")) for row in data["chr10_complete"] if row.get("RU_count")]

    shortest_4qa = min(ru_4qa) if ru_4qa else None
    contracted = shortest_4qa is not None and shortest_4qa <= contraction_threshold
    partial_only_4qa = counts["4qA_all"] > 0 and counts["4qA_complete"] == 0

    if contracted:
        headline = (
            f"Complete permissive 4qA reads with <= {contraction_threshold} repeat units were observed."
        )
    elif shortest_4qa is not None:
        headline = (
            f"Complete 4qA reads were observed, but the shortest measured array exceeded {contraction_threshold} repeat units."
        )
    elif partial_only_4qa:
        headline = "Partial 4qA evidence was observed, but no complete 4qA-spanning read was classified."
    else:
        headline = "No 4qA-classified reads were observed in the current sample."

    flagstat_headline = "NA"
    if os.path.exists(flagstat_txt):
        with open(flagstat_txt, "r", encoding="utf-8") as fh:
            flagstat_headline = fh.readline().strip()

    coverage_rows = []
    if os.path.exists(coverage_tsv):
        with open(coverage_tsv, "r", encoding="utf-8") as fh:
            for raw in fh:
                if not raw.strip() or raw.startswith("#"):
                    continue
                parts = raw.rstrip("\n").split("\t")
                if len(parts) >= 7:
                    coverage_rows.append(parts[:7])

    hap_rows = read_tsv(hap_tsv)
    has_hp_tags = any((row.get("has_hp_tags", "").lower() == "true") for row in hap_rows)
    dominant_hp_4qa = "NA"
    hp_counts = defaultdict(int)
    for row in hap_rows:
        if row.get("subset_name") == "4qA_all":
            hp = row.get("hp", "NA")
            hp_counts[hp] += to_int(row.get("count"))
    if hp_counts:
        dominant_hp_4qa = max(hp_counts.items(), key=lambda kv: kv[1])[0]

    methyl_rows = read_tsv(methyl_tsv)
    methyl_by_subset = {row["subset_name"]: row for row in methyl_rows if "subset_name" in row}
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
        row.get("PAS.type", "NA") for row in data["4qA_complete"] if row.get("PAS.type")
    )
    top_4qa_rows = sorted(
        data["4qA_complete"],
        key=lambda row: (to_int(row.get("RU_count"), 10**9), row.get("read.id", "")),
    )[:15]
    top_chimeric_rows = sorted(
        data["chimeric"],
        key=lambda row: (to_int(row.get("RU_count"), 10**9), row.get("read.id", "")),
    )[:12]

    subset_manifest = os.path.join(subset_dir, f"{sample_id}.classified-subsets.manifest.tsv")
    subset_rows = read_tsv(subset_manifest)

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
        "methyl_4qA_all_mean_pct": methyl_by_subset.get("4qA_all", {}).get(
            "mean_percent_modified", "NA"
        ),
        "methyl_4qA_complete_mean_pct": methyl_by_subset.get("4qA_complete", {}).get(
            "mean_percent_modified", "NA"
        ),
        "methyl_chimeric_mean_pct": methyl_by_subset.get("chimeric", {}).get(
            "mean_percent_modified", "NA"
        ),
        "variant_status": variant.get("status", "unknown"),
        "variant_total_clair3_variants": variant.get("total_clair3_variants", "0"),
        "variant_relevant_records": variant.get("relevant_variant_records", "0"),
        "variant_clinvar_fshd_hits": variant.get("clinvar_fshd_hits", "0"),
        "variant_hg38_sv_records": variant.get("hg38_sv_records", "0"),
        "variant_t2t_sv_records": variant.get("t2t_sv_records", "0"),
        "variant_relevant_genes": variant.get("relevant_genes", "NA"),
    }

    with open(summary_out, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(summary_fields.keys()), delimiter="\t")
        writer.writeheader()
        writer.writerow(summary_fields)

    top_4qa_table = table_html(
        ["Read", "RU count", "PAS", "Warning", "Repeat order"],
        [
            [
                html.escape(row.get("read.id", "")),
                html.escape(str(row.get("RU_count", ""))),
                html.escape(row.get("PAS.type", "NA")),
                html.escape(truncate(row.get("warning", "NA"), 80)),
                html.escape(truncate(row.get("repeat_sequence", "NA"), 90)),
            ]
            for row in top_4qa_rows
        ],
    )

    chimeric_table = table_html(
        ["Read", "RU count", "Front", "Back", "Warning"],
        [
            [
                html.escape(row.get("read.id", "")),
                html.escape(str(row.get("RU_count", ""))),
                html.escape(row.get("D4F104S1_type", "NA")),
                html.escape(row.get("pLAM_type", "NA")),
                html.escape(truncate(row.get("warning", "NA"), 90)),
            ]
            for row in top_chimeric_rows
        ],
    )

    coverage_table = table_html(
        ["Region", "Start", "End", "Reads", "Covered bp", "Breadth", "Mean depth"],
        [[html.escape(col) for col in row] for row in coverage_rows],
        "compact",
    )

    subset_table = table_html(
        ["Subset", "Read IDs", "Status"],
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

    hap_table = table_html(
        ["Subset", "HP", "Count", "Total subset reads"],
        [
            [
                html.escape(row.get("subset_name", "")),
                html.escape(row.get("hp", "")),
                html.escape(row.get("count", "")),
                html.escape(row.get("total_subset_reads", "")),
            ]
            for row in hap_rows
        ],
        "compact",
    )

    methyl_table = table_html(
        ["Subset", "Status", "Mean % modified", "Subset reads", "Donor reads"],
        [
            [
                html.escape(row.get("subset_name", "")),
                html.escape(row.get("status", "")),
                html.escape(row.get("mean_percent_modified", "NA")),
                html.escape(row.get("subset_read_count", "0")),
                html.escape(row.get("donor_read_count", "0")),
            ]
            for row in methyl_rows
        ],
        "compact",
    )

    variant_table = table_html(
        ["Status", "Clair3 calls", "Relevant SNV/indels", "ClinVar FSHD hits", "HG38 SVs", "T2T SVs", "Genes"],
        [[
            html.escape(variant.get("status", "unknown")),
            html.escape(variant.get("total_clair3_variants", "0")),
            html.escape(variant.get("relevant_variant_records", "0")),
            html.escape(variant.get("clinvar_fshd_hits", "0")),
            html.escape(variant.get("hg38_sv_records", "0")),
            html.escape(variant.get("t2t_sv_records", "0")),
            html.escape(variant.get("relevant_genes", "NA")),
        ]],
        "compact",
    )

    pas_note = ", ".join(f"{k}: {v}" for k, v in sorted(pas_counter.items())) if pas_counter else "No complete 4qA PAS assignments available."

    html_body = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>{html.escape(sample_id)} FSHD report</title>
  <style>
    :root {{
      --ink: #2c1f18;
      --paper: #fbf6ef;
      --card: rgba(255,255,255,0.75);
      --edge: rgba(93,62,39,0.18);
      --accent: #c65d32;
      --accent-soft: #efc7b3;
      --teal: #2f7f73;
      --gold: #b7852b;
      --muted: #6f6259;
      --shadow: 0 16px 40px rgba(108, 74, 48, 0.12);
    }}
    * {{ box-sizing: border-box; }}
    body {{
      margin: 0;
      font-family: "Georgia", "Iowan Old Style", serif;
      color: var(--ink);
      background:
        radial-gradient(circle at top left, rgba(226, 150, 113, 0.35), transparent 28%),
        radial-gradient(circle at top right, rgba(78, 143, 133, 0.22), transparent 22%),
        linear-gradient(180deg, #fdf8f3 0%, #f6efe6 100%);
    }}
    .shell {{ max-width: 1220px; margin: 0 auto; padding: 36px 24px 48px; }}
    .hero {{
      background: linear-gradient(135deg, rgba(255,255,255,0.82), rgba(250,235,220,0.82));
      border: 1px solid var(--edge);
      border-radius: 28px;
      padding: 30px;
      box-shadow: var(--shadow);
      position: relative;
      overflow: hidden;
    }}
    .hero::after {{
      content: "";
      position: absolute;
      inset: auto -100px -120px auto;
      width: 260px;
      height: 260px;
      border-radius: 50%;
      background: radial-gradient(circle, rgba(198,93,50,0.14), transparent 72%);
    }}
    h1 {{ margin: 0 0 10px; font-size: 2.35rem; line-height: 1.05; }}
    h2 {{ margin: 0 0 14px; font-size: 1.2rem; letter-spacing: 0.02em; }}
    p {{ line-height: 1.55; }}
    .lede {{ max-width: 840px; font-size: 1.08rem; color: #3c2c22; }}
    .badge {{
      display: inline-flex;
      align-items: center;
      gap: 8px;
      padding: 8px 14px;
      border-radius: 999px;
      font-size: 0.92rem;
      background: rgba(198, 93, 50, 0.12);
      color: #8d3b17;
      border: 1px solid rgba(198, 93, 50, 0.18);
    }}
    .grid {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
      gap: 16px;
      margin-top: 22px;
    }}
    .metric {{
      background: var(--card);
      border: 1px solid var(--edge);
      border-radius: 18px;
      padding: 16px 18px;
      box-shadow: var(--shadow);
    }}
    .metric-label {{ font-size: 0.9rem; text-transform: uppercase; letter-spacing: 0.08em; color: var(--muted); }}
    .metric-value {{ font-size: 2rem; margin-top: 8px; font-weight: 700; }}
    .metric-note {{ margin-top: 6px; color: var(--muted); font-size: 0.92rem; }}
    .section {{
      margin-top: 26px;
      background: rgba(255,255,255,0.78);
      border: 1px solid var(--edge);
      border-radius: 24px;
      padding: 24px;
      box-shadow: var(--shadow);
    }}
    .section-grid {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(320px, 1fr));
      gap: 22px;
    }}
    .muted {{ color: var(--muted); }}
    table {{
      width: 100%;
      border-collapse: collapse;
      margin-top: 14px;
      font-size: 0.95rem;
      background: rgba(255,255,255,0.65);
      border-radius: 16px;
      overflow: hidden;
    }}
    th, td {{
      padding: 10px 12px;
      border-bottom: 1px solid rgba(93,62,39,0.12);
      text-align: left;
      vertical-align: top;
    }}
    th {{ font-size: 0.84rem; text-transform: uppercase; letter-spacing: 0.06em; color: var(--muted); background: rgba(191,153,117,0.10); }}
    tr:last-child td {{ border-bottom: 0; }}
    .compact td, .compact th {{ padding: 8px 10px; font-size: 0.9rem; }}
    .histogram {{ width: 100%; height: auto; background: rgba(255,250,246,0.82); border-radius: 18px; border: 1px solid rgba(93,62,39,0.10); }}
    .histogram .title {{ font-size: 16px; fill: #4b3628; font-weight: 700; }}
    .histogram .axis {{ font-size: 12px; fill: #6f6259; }}
    .histogram .count {{ font-size: 12px; fill: #4b3628; }}
    .callout {{
      border-left: 5px solid var(--accent);
      background: rgba(198,93,50,0.08);
      padding: 14px 16px;
      border-radius: 14px;
      margin-top: 18px;
    }}
  </style>
</head>
<body>
  <div class="shell">
    <section class="hero">
      <div class="badge">FSHD D4Z4 / DUX4 profiling</div>
      <h1>{html.escape(sample_id)}</h1>
      <p class="lede">{html.escape(headline)}</p>
      <p class="muted"><strong>Locus BAM flagstat:</strong> {html.escape(flagstat_headline)}</p>
      <div class="grid">
        {metric_card("4qA complete reads", counts["4qA_complete"], f"Shortest RU: {shortest_4qa if shortest_4qa is not None else 'NA'}")}
        {metric_card("4qA all reads", counts["4qA_all"], "Permissive haplotype evidence")}
        {metric_card("4qB complete reads", counts["4qB_complete"], "Non-permissive chr4 haplotype")}
        {metric_card("chr10 complete reads", counts["chr10_complete"], "Homologous reference reads")}
        {metric_card("Chimeric reads", counts["chimeric"], "Potential rearrangements or hybrids")}
        {metric_card("Dominant 4qA HP", dominant_hp_4qa, "When HP tags are present")}
        {metric_card("Variant branch", variant.get("status", "unknown"), f"Relevant genes: {variant.get('relevant_genes', 'NA')}")}
      </div>
      <div class="callout">
        <strong>PAS summary:</strong> {html.escape(pas_note)}
      </div>
    </section>

    <section class="section">
      <h2>Repeat Unit Distributions</h2>
      <div class="section-grid">
        <div>{svg_hist(ru_4qa, "#c65d32", "Complete 4qA repeat-unit counts")}</div>
        <div>{svg_hist(ru_4qb, "#2f7f73", "Complete 4qB repeat-unit counts")}</div>
        <div>{svg_hist(ru_10q, "#b7852b", "Complete chr10 repeat-unit counts")}</div>
      </div>
    </section>

    <section class="section">
      <h2>Top 4qA Complete Reads</h2>
      {top_4qa_table}
    </section>

    <section class="section">
      <h2>Chimeric Read Highlights</h2>
      {chimeric_table}
    </section>

    <section class="section">
      <h2>Subset and Coverage Tracking</h2>
      <div class="section-grid">
        <div>
          <h3>Classified subsets</h3>
          {subset_table}
        </div>
        <div>
          <h3>Coverage over target intervals</h3>
          {coverage_table}
        </div>
      </div>
    </section>

    <section class="section">
      <h2>Haplotags and Methylation</h2>
      <div class="section-grid">
        <div>
          <h3>HP tag breakdown</h3>
          <p class="muted">HP tags detected: {html.escape(str(has_hp_tags).lower())}</p>
          {hap_table}
        </div>
        <div>
          <h3>Repaired methylation summary</h3>
          {methyl_table}
        </div>
      </div>
    </section>

    <section class="section">
      <h2>Variant Calling</h2>
      <p class="muted">HG38 realignment, Clair3 phasing, Whatshap haplotagging, Sniffles2 structural variant calling, and ClinVar/snpEff filtering for FSHD-relevant loci.</p>
      {variant_table}
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
