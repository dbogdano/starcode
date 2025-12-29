#!/usr/bin/env bash
set -euo pipefail

# Benchmark Starcode counts-input mode vs raw FASTQ on macOS/Linux
# Creates a temporary dataset from misc/sample_counts.tsv with adjustable SCALE
# Reports time and memory (max RSS) for both runs.

ROOT_DIR="$(cd "$(dirname "$0")"/.. && pwd)"
cd "$ROOT_DIR"

# Ensure starcode binary exists
if [[ ! -x ./starcode ]]; then
  echo "Building starcode..."
  make -s
fi

SCALE=${SCALE:-50000}   # counts per unique (default 50k)
DIST=${DIST:-1}
THREADS=${THREADS:-1}
TMPDIR="$(mktemp -d)"
trap 'rm -rf "$TMPDIR"' EXIT

TSV_ORIG="misc/sample_counts.tsv"
TSV="${TMPDIR}/counts.tsv"
FASTQ="${TMPDIR}/reads.fastq"

# Prepare TSV with actual numeric counts
awk -v S="$SCALE" 'BEGIN{OFS="\t"} {gsub(/SCALE/, S); print $0}' "$TSV_ORIG" > "$TSV"

# Generate FASTQ by expanding counts per sequence
# For each TSV line seq\tcount, emit count FASTQ records with dummy header and quality
awk 'BEGIN{FS="\t"; OFS="\n"}
{
  seq=$1; count=$2; q="";
  for (i=1;i<=length(seq);i++) q=q"I";  # Q40
  for (i=1;i<=count;i++) {
    print "@read_" NR "_" i, seq, "+", q;
  }
}' "$TSV" > "$FASTQ"

# Pick a portable time command
TIME_CMD="/usr/bin/time"
TIME_ARGS="-l"  # macOS
MEM_KEY_MAC="maximum resident set size"
MEM_KEY_LNX="Maximum resident set size"
if [[ ! -x $TIME_CMD ]]; then
  if command -v gtime >/dev/null 2>&1; then
    TIME_CMD="gtime"
    TIME_ARGS="-v"
  else
    # Fallback: use shell built-in time without -v; RSS may be unavailable
    TIME_CMD="/usr/bin/time"
    TIME_ARGS="-l"
  fi
fi

run_and_measure() {
  local label="$1"; shift
  local logfile="${TMPDIR}/${label}.time"
  # shellcheck disable=SC2086
  # shellcheck disable=SC2086
  $TIME_CMD $TIME_ARGS "$@" 2>"$logfile" || true
  local mem=""
  mem=$(grep -E "${MEM_KEY_MAC}|${MEM_KEY_LNX}" "$logfile" | sed -E 's/.*: *([0-9]+).*/\1/' | tail -n1)
  echo "--- ${label} time output (excerpt) ---"
  grep -E "${MEM_KEY_MAC}|${MEM_KEY_LNX}|real|user|sys" "$logfile" || true
  echo "$label: mem=${mem:-unknown}"
}

echo "Running counts-input mode (uniques only)..."
run_and_measure counts "./starcode" --counts-input -i "$TSV" -o /dev/null --threads "$THREADS" --dist "$DIST"

echo "Running raw FASTQ mode (all reads)..."
run_and_measure raw "./starcode" -i "$FASTQ" -o /dev/null --threads "$THREADS" --dist "$DIST"

echo "\nSummary:"
if [[ -f "${TMPDIR}/counts.time" && -f "${TMPDIR}/raw.time" ]]; then
  counts_mem=$(sed -nE 's/.*(maximum resident set size|Maximum resident set size).*: *([0-9]+).*/\2/p' "${TMPDIR}/counts.time" | tail -n1)
  raw_mem=$(sed -nE 's/.*(maximum resident set size|Maximum resident set size).*: *([0-9]+).*/\2/p' "${TMPDIR}/raw.time" | tail -n1)
  echo "  counts-input max RSS: ${counts_mem:-unknown}"
  echo "  raw FASTQ   max RSS: ${raw_mem:-unknown}"
fi

echo "\nFiles:"
echo "  TSV:   $TSV"
echo "  FASTQ: $FASTQ"
