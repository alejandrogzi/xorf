#!/usr/bin/env bash
set -euo pipefail

profile=${1:-test}
mode=${2:-run}
engine=${TEST_ENGINE:-docker}
root=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)

case "$profile" in
  test) profiles=("$profile") ;;
  all) profiles=(test) ;;
  *) echo "usage: $0 {test|all} [run|verify]" >&2; exit 2 ;;
esac

if [[ "$mode" != verify ]]; then
  nextflow_bin=${NEXTFLOW_BIN:-$(command -v nextflow || true)}

  if [[ -z "$nextflow_bin" ]]; then
    echo "Nextflow is not on PATH; run 'mamba activate nextflow' or set NEXTFLOW_BIN" >&2
    exit 2
  fi

  for selected in "${profiles[@]}"; do
    "$nextflow_bin" run "$root/src" -profile "$selected,$engine" -ansi-log false
  done
fi

python "$root/assets/scripts/verify.py" "$profile"
