#!/usr/bin/env python3

import csv
import html
import sys


def main() -> int:
    tsv_out = sys.argv[1]
    html_out = sys.argv[2]
    paths = sys.argv[3:]

    rows = []
    for path in paths:
        with open(path, "r", encoding="utf-8") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            rows.extend(list(reader))

    fieldnames = list(rows[0].keys()) if rows else [
        "sample_id",
        "contracted_4qA_observed",
        "shortest_complete_4qA_ru",
        "complete_4qA_reads",
        "all_4qA_reads",
        "complete_4qB_reads",
        "all_4qB_reads",
        "complete_chr10_reads",
        "all_chr10_reads",
        "chimeric_reads",
        "has_hp_tags",
        "dominant_hp_4qA",
        "methyl_4qA_all_mean_pct",
        "methyl_4qA_complete_mean_pct",
        "methyl_chimeric_mean_pct",
    ]

    with open(tsv_out, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    contracted_count = sum(
        1 for row in rows if row.get("contracted_4qA_observed", "").lower() == "true"
    )

    table_rows = []
    for row in rows:
        table_rows.append(
            "<tr>"
            + "".join(
                f"<td>{html.escape(str(row.get(field, '')))}</td>" for field in fieldnames
            )
            + "</tr>"
        )

    html_doc = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>FSHD cohort summary</title>
  <style>
    body {{
      margin: 0;
      font-family: "Georgia", serif;
      background: linear-gradient(180deg, #faf4ee 0%, #f3eadf 100%);
      color: #31231c;
    }}
    .shell {{
      max-width: 1240px;
      margin: 0 auto;
      padding: 32px 24px 48px;
    }}
    .hero, .panel {{
      background: rgba(255,255,255,0.82);
      border: 1px solid rgba(93,62,39,0.16);
      border-radius: 24px;
      padding: 24px;
      box-shadow: 0 14px 32px rgba(88,63,45,0.08);
    }}
    .hero p {{ color: #6b5d54; }}
    table {{
      width: 100%;
      border-collapse: collapse;
      margin-top: 16px;
      font-size: 0.92rem;
    }}
    th, td {{
      padding: 9px 10px;
      border-bottom: 1px solid rgba(93,62,39,0.10);
      text-align: left;
    }}
    th {{
      font-size: 0.82rem;
      text-transform: uppercase;
      letter-spacing: 0.05em;
      color: #6b5d54;
    }}
  </style>
</head>
<body>
  <div class="shell">
    <section class="hero">
      <h1>FSHD cohort summary</h1>
      <p>Samples summarized: {len(rows)}. Samples with contracted complete 4qA signal: {contracted_count}.</p>
    </section>
    <section class="panel" style="margin-top: 24px;">
      <table>
        <thead>
          <tr>{''.join(f'<th>{html.escape(field)}</th>' for field in fieldnames)}</tr>
        </thead>
        <tbody>
          {''.join(table_rows)}
        </tbody>
      </table>
    </section>
  </div>
</body>
</html>
"""

    with open(html_out, "w", encoding="utf-8") as fh:
        fh.write(html_doc)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
