#!/usr/bin/env nextflow

process PAS_CHECK_FSHD {

  tag "${sample_id}"

  publishDir "${params.output_dir}/pas-check-fshd/${sample_id}", mode: 'copy', overwrite: true

  input:
    tuple val(sample_id), path(locus_bam), path(locus_bai)

  output:
    tuple val(sample_id), path("${sample_id}.PAS.txt")

  script:
    """
    set -euo pipefail

    if [[ ! -f "${locus_bam}.bai" ]]; then
      ln -sf "${locus_bai}" "${locus_bam}.bai"
    fi

    python3 - "${locus_bam}" "${sample_id}.PAS.txt" <<'PY'
import sys
import pysam

bam_path = sys.argv[1]
out_path = sys.argv[2]

regions = [
    ("chr4", 193543619, 6, "ATTAAA", "ATCAAA"),
    ("chr10", 134726542, 6, "ATTAAA", "ATCAAA"),
]

rows = []
bam = pysam.AlignmentFile(bam_path, "rb")
for chrom, start_pos, length, q4a_seq, q10a_seq in regions:
    start0 = start_pos - 1
    end0 = start0 + length
    for read in bam.fetch(chrom, start0, end0):
        if read.is_secondary or read.is_supplementary or read.is_unmapped:
            continue
        aligned_pairs = dict(read.get_aligned_pairs(matches_only=False))
        reverse_lookup = {ref_pos: query_pos for query_pos, ref_pos in aligned_pairs.items() if ref_pos is not None}
        bases = []
        for i in range(length):
            ref_pos = start0 + i
            query_pos = reverse_lookup.get(ref_pos)
            if query_pos is None or read.query_sequence is None:
                bases.append('-')
            else:
                bases.append(read.query_sequence[query_pos])
        read_seq = ''.join(bases)
        if read_seq == q4a_seq:
            status = "4qA_PAS"
        elif read_seq == q10a_seq:
            status = "10qA_PAS"
        else:
            status = "PAS_disrupted"
        rows.append((read.query_name, read_seq, status))
bam.close()

with open(out_path, "w", encoding="utf-8") as fh:
    fh.write("read.id\\tPAS.seq\\tPAS.type\\n")
    seen = set()
    for row in rows:
        if row[0] in seen:
            continue
        seen.add(row[0])
        fh.write("\\t".join(row) + "\\n")
PY
    """
}
