#!/usr/bin/env bash
set -euo pipefail

usage() { echo "Usage: $0 -i <indir> [-o <outdir>] -p <jobs>" >&2; exit 1; }

indir="" outdir="" jobs=""
while getopts ":i:o:p:h" opt; do
  case "$opt" in
    i) indir="$OPTARG" ;;
    o) outdir="$OPTARG" ;;
    p) jobs="$OPTARG" ;;
    h) usage ;;
    *) usage ;;
  esac
done

[[ -n "$indir" && -n "$jobs" ]] || usage
[[ -d "$indir" ]] || { echo "ERROR: not a directory: $indir" >&2; exit 1; }
[[ "$jobs" =~ ^[0-9]+$ ]] || { echo "ERROR: -p must be an integer" >&2; exit 1; }
(( jobs >= 1 )) || { echo "ERROR: -p must be >= 1" >&2; exit 1; }

outdir="${outdir:-${indir%/}/.fastq_read_stats}"
mkdir -p "$outdir"

decomp() {
  local f="$1"
  if [[ "$f" == *.gz ]]; then
    if command -v pigz >/dev/null 2>&1; then pigz -dc "$f"; else gzip -dc "$f"; fi
  else
    cat "$f"
  fi
}

process_one() {
  local f="$1"
  local outdir="$2"

  local bn base out
  bn="$(basename "$f")"
  base="${bn%.gz}"
  base="${base%.fastq}"
  base="${base%.fq}"
  out="$outdir/${base}.meanQ.tsv"

  {
    printf "read_id\tlength\tmean_Q\n"
    decomp "$f" | awk '
      # Robust FASTQ parser (handles wrapped sequence/quality lines)
      # Computes mean Q from Phred+33 qualities.
      function meanq(qual,   i,c,sum,L) {
        L = length(qual)
        if (L == 0) return "NA"
        sum = 0
        for (i=1; i<=L; i++) {
          c = substr(qual, i, 1)
          sum += (index(ORD, c) - 1) - 33
        }
        return sum / L
      }
      BEGIN {
        # ASCII lookup: index(ORD, char)-1 = ASCII code
        ORD = ""
        for (k=0; k<128; k++) ORD = ORD sprintf("%c", k)
        state = 0
      }
      state==0 {
        # Expect header
        if ($0 ~ /^@/) {
          id = $0
          sub(/^@/,"",id); sub(/ .*/,"",id)
          seq = ""
          qual = ""
          state = 1
        }
        next
      }
      state==1 {
        # Read sequence lines until +
        if ($0 ~ /^\+/) { state = 2; next }
        seq = seq $0
        next
      }
      state==2 {
        # Read quality lines until length matches sequence
        qual = qual $0
        if (length(qual) >= length(seq)) {
          L = length(seq)
          mq = meanq(substr(qual,1,L))
          printf "%s\t%d\t%.4f\n", id, L, mq
          state = 0
        }
        next
      }
    '
  } > "$out"
}
export -f process_one decomp

# Recursive file discovery
mapfile -d '' files < <(find "$indir" -type f \( -name "*.fq" -o -name "*.fastq" -o -name "*.fq.gz" -o -name "*.fastq.gz" \) -print0)

((${#files[@]} > 0)) || { echo "No FASTQ files found under: $indir" >&2; exit 1; }

# Parallel per input file
if command -v parallel >/dev/null 2>&1; then
  printf "%s\0" "${files[@]}" \
    | parallel -0 -j "$jobs" --halt now,fail=1 --line-buffer \
        'process_one "{}" "'"$outdir"'"'
else
  printf "%s\0" "${files[@]}" \
    | xargs -0 -n 1 -P "$jobs" -I {} bash -c 'process_one "$1" "$2"' _ {} "$outdir"
fi

echo "Done. Outputs in: $outdir" >&2
