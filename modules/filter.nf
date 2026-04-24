/*
 * FILTER_SAMPLES
 *
 * Reads QUAST transposed_report.tsv and applies three QC thresholds:
 *   1. Total assembly length  >= params.min_genome_size
 *   2. Contig count           <= params.max_contigs
 *   3. N50                    >= params.min_n50
 *
 * Samples that pass all three emit on the `passed` channel.
 * Samples that fail emit on the `failed` channel with a human-readable
 * fail_reason.txt explaining every threshold that was breached.
 *
 * Input tuple:  ( sample_id, quast_dir, fasta )
 * Output:
 *   passed  -> ( sample_id, fasta_copy )   — only passing samples
 *   failed  -> ( sample_id, fail_reason.txt ) — only failing samples
 */

process FILTER_SAMPLES {

    tag "$sample"

    // Minimal Python 3.10 environment — no heavy deps needed
    conda "conda-forge::python=3.10"

    // Keep QC logs alongside the main results
    publishDir "${params.outdir}/qc/filter", mode: 'copy', saveAs: { fn ->
        fn == "fail_reason.txt" ? "${sample}_fail_reason.txt" : null
    }

    input:
    tuple val(sample), path(quast_dir), path(fasta)

    output:
    tuple val(sample), path("passed_fasta.fasta"), emit: passed,  optional: true
    tuple val(sample), path("fail_reason.txt"),    emit: failed,  optional: true

    script:
    // All threshold values are injected from params at pipeline compile time
    """
    python3 << 'EOF'
import sys
import os
import shutil
import csv

# ── Thresholds (injected by Nextflow from params) ────────────────────────────
MIN_SIZE    = int("${params.min_genome_size}")
MAX_CONTIGS = int("${params.max_contigs}")
MIN_N50     = int("${params.min_n50}")
SAMPLE      = "${sample}"
QUAST_TSV   = "${quast_dir}/transposed_report.tsv"
FASTA_IN    = "${fasta}"

# ── Parse QUAST transposed_report.tsv ────────────────────────────────────────
# transposed_report.tsv has metrics as columns, one assembly per row.
# Column names vary slightly between QUAST versions so we normalise them.
if not os.path.isfile(QUAST_TSV):
    sys.exit(f"[FILTER_SAMPLES] ERROR: {QUAST_TSV} not found for sample {SAMPLE}")

with open(QUAST_TSV, "r") as fh:
    reader = csv.DictReader(fh, delimiter="\\t")
    rows = list(reader)

if not rows:
    sys.exit(f"[FILTER_SAMPLES] ERROR: {QUAST_TSV} is empty for sample {SAMPLE}")

row = rows[0]   # one row per assembly in transposed format

# ── Column name normalisation ─────────────────────────────────────────────────
# QUAST sometimes writes "# contigs (>= 0 bp)" instead of "# contigs".
# We look for partial matches so version differences don't break the pipeline.
def find_col(row, candidates):
    for key in row:
        for c in candidates:
            if c.lower() in key.lower():
                return key
    return None

col_size    = find_col(row, ["Total length"])
col_contigs = find_col(row, ["# contigs"])
col_n50     = find_col(row, ["N50"])

missing = [name for name, col in
           [("Total length", col_size),
            ("# contigs",    col_contigs),
            ("N50",          col_n50)]
           if col is None]

if missing:
    available = list(row.keys())
    sys.exit(
        f"[FILTER_SAMPLES] ERROR: Could not find columns {missing} in {QUAST_TSV}.\\n"
        f"Available columns: {available}"
    )

# ── Extract values ────────────────────────────────────────────────────────────
try:
    genome_size = int(row[col_size].replace(",", ""))
    contigs     = int(row[col_contigs].replace(",", ""))
    n50         = int(row[col_n50].replace(",", ""))
except ValueError as e:
    sys.exit(f"[FILTER_SAMPLES] ERROR: Could not parse numeric value — {e}")

# ── Apply thresholds ──────────────────────────────────────────────────────────
fail_reasons = []

if genome_size < MIN_SIZE:
    fail_reasons.append(
        f"Assembly size {genome_size:,} bp is below minimum {MIN_SIZE:,} bp"
    )

if contigs > MAX_CONTIGS:
    fail_reasons.append(
        f"Contig count {contigs} exceeds maximum {MAX_CONTIGS} "
        f"(assembly too fragmented)"
    )

if n50 < MIN_N50:
    fail_reasons.append(
        f"N50 {n50:,} bp is below minimum {MIN_N50:,} bp "
        f"(contigs too short for reliable gene detection)"
    )

# ── Emit result ───────────────────────────────────────────────────────────────
if not fail_reasons:
    # Copy — not symlink — so Nextflow can safely stage the file across
    # work directories and on network/cloud executors (S3, GCS, etc.)
    shutil.copy(FASTA_IN, "passed_fasta.fasta")
    print(f"[FILTER_SAMPLES] PASS  {SAMPLE} | "
          f"size={genome_size:,}  contigs={contigs}  N50={n50:,}")
else:
    with open("fail_reason.txt", "w") as fh:
        fh.write(f"Sample: {SAMPLE}\\n")
        fh.write(f"Assembly size : {genome_size:,} bp\\n")
        fh.write(f"Contig count  : {contigs}\\n")
        fh.write(f"N50           : {n50:,} bp\\n")
        fh.write("\\nFail reasons:\\n")
        for i, reason in enumerate(fail_reasons, 1):
            fh.write(f"  {i}. {reason}\\n")
    print(f"[FILTER_SAMPLES] FAIL  {SAMPLE} — " + "; ".join(fail_reasons),
          file=sys.stderr)

EOF
    """
}
