#!/usr/bin/env bash
# Runs every check for the paper and writes checks/output.txt. Needs Python 3 only.
set -euo pipefail
cd "$(dirname "$0")"
{
  echo "== Lemma 2.1: walk enumerators and the determinant formula"
  python3 path_lemma.py
  echo
  echo "== Section 3: shear coordinates by the S/Z rule equal d_(i-1) - d_i"
  python3 shear_rule.py 10
  echo
  echo "== Main theorem: mutation equals the weighted ideal sum, and all roots are simple and negative"
  python3 disk_polygons.py 9
  echo
  echo "== Section 4: the examples of Kang, Lee and Lim"
  python3 kll_examples.py
  echo
  echo "all checks passed"
} | tee output.txt
