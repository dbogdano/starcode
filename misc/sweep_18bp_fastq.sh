#!/usr/bin/env bash
set -Eeuo pipefail

# Parameter sweep for barcode clustering on FASTQ input.
# - Subsamples first N reads from FASTQ/FASTQ.GZ
# - Converts to counts TSV (SEQ\tCOUNT)
# - Runs starcode across algorithm/parameter combinations
# - Writes one summary TSV with cluster count + timing + max RSS

ROOT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
STARCODE_BIN="${STARCODE_BIN:-$ROOT_DIR/starcode}"

log() {
  printf "[%s] %s\n" "$(date +%H:%M:%S)" "$*"
}

usage() {
  cat <<'EOF'
Usage:
  misc/sweep_18bp_fastq.sh -i INPUT.fastq[.gz] [options]

Required:
  -i FILE              Input FASTQ or FASTQ.GZ

Options:
  -n READS             Reads to subsample from head (default: 5000000)
  -t THREADS           Starcode threads (default: 1)
  -a PATTERN           Optional --allow pattern (e.g. NNHNNYRNNNNYRNNHNN)
  -o DIR               Output directory (default: sweep_YYYYmmdd_HHMMSS)
  --dists "LIST"       Space-separated distances (default: "1 2")
  --ratios "LIST"      Space-separated cluster ratios for MP (default: "2 3 5")
  --sort-mem MEM       sort memory cap (default: 8G)
  --print-clusters     Include --print-clusters in all starcode runs
  --keep-temp          Keep temporary files
  -h, --help           Show this help

Output:
  summary.tsv columns:
    algo\tdist\tratio\tclusters\tmax_rss\treal_s\toutput_file

Algorithms swept:
  mp, sphere, cc, cc_stream
EOF
}

INPUT=""
READS=5000000
THREADS=1
ALLOW_PATTERN=""
OUTDIR=""
DISTS="1 2"
RATIOS="2 3 5"
SORT_MEM="8G"
PRINT_CLUSTERS=0
KEEP_TEMP=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i) INPUT="$2"; shift 2 ;;
    -n) READS="$2"; shift 2 ;;
    -t) THREADS="$2"; shift 2 ;;
    -a) ALLOW_PATTERN="$2"; shift 2 ;;
    -o) OUTDIR="$2"; shift 2 ;;
    --dists) DISTS="$2"; shift 2 ;;
    --ratios) RATIOS="$2"; shift 2 ;;
    --sort-mem) SORT_MEM="$2"; shift 2 ;;
    --print-clusters) PRINT_CLUSTERS=1; shift ;;
    --keep-temp) KEEP_TEMP=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 2 ;;
  esac
done

if [[ -z "$INPUT" ]]; then
  echo "Error: -i INPUT.fastq[.gz] is required" >&2
  usage
  exit 2
fi

if [[ ! -f "$INPUT" ]]; then
  echo "Error: input not found: $INPUT" >&2
  exit 2
fi

if [[ ! -x "$STARCODE_BIN" ]]; then
  echo "starcode binary not found at $STARCODE_BIN; attempting build..."
  (cd "$ROOT_DIR" && make -s)
fi

if [[ -z "$OUTDIR" ]]; then
  OUTDIR="$ROOT_DIR/sweep_$(date +%Y%m%d_%H%M%S)"
fi
mkdir -p "$OUTDIR"

RUN_LOG="$OUTDIR/run.log"
ERR_LOG="$OUTDIR/error.log"

on_error() {
  local line="$1"
  local code="$2"
  {
    echo "ERROR: sweep failed"
    echo "line: $line"
    echo "exit_code: $code"
    echo "last_command: ${BASH_COMMAND:-unknown}"
  } >> "$ERR_LOG"
  echo "Sweep failed. See: $ERR_LOG" >&2
}
trap 'on_error ${LINENO} $?' ERR

TMPDIR_SWEEP="$(mktemp -d "${TMPDIR:-/tmp}/starcode_sweep.XXXXXX")"
if [[ "$KEEP_TEMP" -eq 1 ]]; then
  echo "Keeping temp dir: $TMPDIR_SWEEP"
else
  trap 'rm -rf "$TMPDIR_SWEEP"' EXIT
fi

SUB_FASTQ="$TMPDIR_SWEEP/subsample.fastq"
COUNTS_TSV="$TMPDIR_SWEEP/subsample.counts.tsv"
SUMMARY="$OUTDIR/summary.tsv"

# Write summary header early so the output directory is never empty.
printf "algo\tdist\tratio\tclusters\tmax_rss\treal_s\toutput_file\n" > "$SUMMARY"

log "Preparing subsample ($READS reads) ..."
if [[ "$INPUT" == *.gz ]]; then
  gzip -dc "$INPUT" | awk -v max="$READS" '{ print; if (NR % 4 == 0) {r++; if (r >= max) exit} }' > "$SUB_FASTQ"
else
  awk -v max="$READS" '{ print; if (NR % 4 == 0) {r++; if (r >= max) exit} }' "$INPUT" > "$SUB_FASTQ"
fi

if [[ ! -s "$SUB_FASTQ" ]]; then
  echo "Subsample FASTQ is empty. Check input format/path." >> "$ERR_LOG"
  exit 1
fi

log "Converting subsample FASTQ -> counts TSV ..."
awk 'NR % 4 == 2' "$SUB_FASTQ" \
  | LC_ALL=C sort -S "$SORT_MEM" -T "$TMPDIR_SWEEP" \
  | uniq -c \
  | awk '{print $2"\t"$1}' > "$COUNTS_TSV"

if [[ ! -s "$COUNTS_TSV" ]]; then
  echo "Counts TSV is empty after conversion." >> "$ERR_LOG"
  exit 1
fi

TIME_CMD="/usr/bin/time"
TIME_ARGS="-l"
if ! command -v /usr/bin/time >/dev/null 2>&1; then
  if command -v gtime >/dev/null 2>&1; then
    TIME_CMD="gtime"
    TIME_ARGS="-v"
  fi
fi

extract_rss() {
  local log="$1"
  local rss
  rss=$(sed -nE 's/.*(maximum resident set size|Maximum resident set size).*: *([0-9]+).*/\2/p' "$log" | tail -n1)
  if [[ -z "$rss" ]]; then
    echo "NA"
  else
    echo "$rss"
  fi
}

extract_real_s() {
  local log="$1"
  local r
  r=$(sed -nE 's/^\s*([0-9]+\.?[0-9]*)\s+real$/\1/p' "$log" | tail -n1)
  if [[ -z "$r" ]]; then
    r="NA"
  fi
  echo "$r"
}

run_case() {
  local algo="$1"
  local dist="$2"
  local ratio="$3"
  shift 3
  local out_file="$OUTDIR/${algo}_d${dist}_r${ratio}.out"
  local log_file="$OUTDIR/${algo}_d${dist}_r${ratio}.time"

  local cmd=("$STARCODE_BIN" -q --counts-input -i "$COUNTS_TSV" -o "$out_file" --dist "$dist" --threads "$THREADS")

  if [[ -n "$ALLOW_PATTERN" ]]; then
    cmd+=(--allow "$ALLOW_PATTERN")
  fi
  if [[ "$PRINT_CLUSTERS" -eq 1 ]]; then
    cmd+=(--print-clusters)
  fi

  cmd+=("$@")

  # shellcheck disable=SC2086
  "$TIME_CMD" $TIME_ARGS "${cmd[@]}" 2> "$log_file"

  local clusters
  clusters=$(wc -l < "$out_file" | tr -d ' ')
  local rss
  rss=$(extract_rss "$log_file")
  local real_s
  real_s=$(extract_real_s "$log_file")

  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$algo" "$dist" "$ratio" "$clusters" "$rss" "$real_s" "$out_file" >> "$SUMMARY"
}

log "Running sweep ..."
for d in $DISTS; do
  for r in $RATIOS; do
    run_case mp "$d" "$r" --cluster-ratio "$r"
  done
  run_case sphere "$d" NA --sphere
  run_case cc "$d" NA --connected-comp
  run_case cc_stream "$d" NA --connected-comp --stream-clusters
done

log "Done. Summary: $SUMMARY"
log "Tip: sort by clusters/time to find stable fast settings."
