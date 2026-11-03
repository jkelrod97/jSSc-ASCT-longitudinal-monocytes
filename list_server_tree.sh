#!/usr/bin/env bash
# List the folder structure of the project on the server so we can decide what
# the pipeline needs. Read-only: it only runs find/du and writes one text file.
#
#   bash list_server_tree.sh [ROOT]      (default ROOT: /home/jelrod/scl-CITE-SEQ)
#
# Output: ~/sample_swap_logs/server_tree.txt

set -euo pipefail
R="${1:-/home/jelrod/scl-CITE-SEQ}"
OUT=~/sample_swap_logs/server_tree.txt
mkdir -p ~/sample_swap_logs
# Don't descend into cellranger internals / per-sample outputs (huge, not needed here)
PRUNE=( \( -name SC_MULTI_CS -o -name outs -o -name per_sample_outs -o -name '*.__swaptmp_*' \) )

{
  echo "# ROOT: $R   ($(date -Iseconds), $(hostname))"
  echo "# ---- directories"
  find "$R" "${PRUNE[@]}" -prune -print -o -type d -print
  echo "# ---- files (outside cellranger internals, max depth 6)"
  find "$R" -maxdepth 6 "${PRUNE[@]}" -prune -o -type f -printf '%s\t%p\n'
  echo "# ---- symlinks (outside cellranger internals)"
  find "$R" -name SC_MULTI_CS -prune -o -type l -printf '%p -> %l\n'
  echo "# ---- sizes of first two levels"
  du -sh "$R"/*/ "$R"/*/*/ 2>/dev/null || true
} | sed "s#$R/##g" > "$OUT"

echo "Wrote $OUT ($(wc -l < "$OUT") lines)"
