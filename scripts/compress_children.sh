#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# Usage:
#   ./compress_children.sh [-i parent_dir] [-m mode] [-d direction] [-j jobs] [-s suffix]
#
# Description:
#   Compresses or extracts the immediate children of a directory using
#   parallel processing. Files are compressed with pigz and directories are
#   archived as .tar.gz. Each file or directory is processed independently.
#
# Options:
#   -i   Parent directory whose immediate children will be processed
#        Default: current directory (.)
#
#   -m   Mode: what to operate on
#        dirs   → directories only
#        files  → files only
#        both   → directories and files
#        Default: both
#
#   -d   Direction: operation to perform
#        compress → create .gz or .tar.gz archives
#        extract  → extract .gz and .tar.gz archives
#        Default: compress
#
#   -j   Number of parallel jobs (GNU parallel)
#        Default: number of CPU cores
#
#   -s   Optional file suffix filter when mode=files
#        Example: .fasta .fastq .bam
#
# Examples:
#
#   Compress everything in the current directory
#     ./compress_children.sh
#
#   Compress all children of a directory
#     ./compress_children.sh -i /data
#
#   Compress only directories
#     ./compress_children.sh -i /data -m dirs
#
#   Compress only .fasta files
#     ./compress_children.sh -i /data -m files -s .fasta
#
#   Extract all compressed archives
#     ./compress_children.sh -i /data -d extract
#
# Requirements:
#   - GNU parallel
#   - pigz
#   - tar
#
# Notes:
#   - Only immediate children of the parent directory are processed
#   - Each item is handled independently and in parallel
#   - Originals are removed after successful compression or extraction
# =============================================================================

usage() {
  cat <<EOF
Usage: $0 [-i parent_dir] [-m mode] [-d direction] [-j jobs] [-s suffix]

Defaults:
  parent_dir : .
  mode       : both
  direction  : compress
  jobs       : number of CPU cores

Run without arguments to process the current directory.
EOF
  exit 1
}

# defaults
parent_dir="."
mode="both"
direction="compress"
jobs="$(nproc)"
suffix=""

while getopts "i:m:d:j:s:h" opt; do
  case "$opt" in
    i) parent_dir="$OPTARG" ;;
    m) mode="$OPTARG" ;;
    d) direction="$OPTARG" ;;
    j) jobs="$OPTARG" ;;
    s) suffix="$OPTARG" ;;
    h) usage ;;
    *) usage ;;
  esac
done

command -v parallel >/dev/null || { echo "GNU parallel required" >&2; exit 1; }
command -v pigz >/dev/null || { echo "pigz required" >&2; exit 1; }

case "$mode" in
  dirs|files|both) ;;
  *) echo "Invalid mode: $mode" >&2; exit 1 ;;
esac

case "$direction" in
  compress|extract) ;;
  *) echo "Invalid direction: $direction" >&2; exit 1 ;;
esac

if [ ! -d "$parent_dir" ]; then
  echo "Parent directory does not exist: $parent_dir" >&2
  exit 1
fi

process_one() {
  local path base parent out
  path="$(realpath "$1")"
  base="$(basename "$path")"
  parent="$(dirname "$path")"

  if [ "$direction" = "compress" ]; then
    if [ -d "$path" ]; then
      out="${path}.tar.gz"
      [ -e "$out" ] && return 0
      tar -C "$parent" -cf - "$base" | pigz -c > "$out" && rm -rf "$path"

    elif [ -f "$path" ]; then
      case "$base" in
        *.gz|*.tar.gz) return 0 ;;
      esac
      pigz -c "$path" > "${path}.gz" && rm -f "$path"
    fi

  else
    case "$base" in
      *.tar.gz)
        tar -C "$parent" -xzf "$path" && rm -f "$path"
        ;;
      *.gz)
        out="${path%.gz}"
        pigz -dc "$path" > "$out" && rm -f "$path"
        ;;
    esac
  fi
}

export -f process_one
export direction

if [ "$mode" = "dirs" ]; then
  find "$parent_dir" -mindepth 1 -maxdepth 1 -type d -print0 \
    | parallel -0 -j "$jobs" --no-run-if-empty process_one {}

elif [ "$mode" = "files" ]; then
  if [ -n "$suffix" ]; then
    find "$parent_dir" -mindepth 1 -maxdepth 1 -type f -name "*$suffix" -print0 \
      | parallel -0 -j "$jobs" --no-run-if-empty process_one {}
  else
    find "$parent_dir" -mindepth 1 -maxdepth 1 -type f -print0 \
      | parallel -0 -j "$jobs" --no-run-if-empty process_one {}
  fi

else
  find "$parent_dir" -mindepth 1 -maxdepth 1 -print0 \
    | parallel -0 -j "$jobs" --no-run-if-empty process_one {}
fi