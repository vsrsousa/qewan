#!/usr/bin/env bash
set -euo pipefail

# Simple sequential runner for the test flow. Adjust executable names and paths
# to pseudopotentials and outdir as needed.

OUTDIR=./out
mkdir -p "$OUTDIR"

echo "1) SCF"
pw.x -in scf.in > scf.out 2>&1

echo "2) BANDS"
pw.x -in bands.in > bands.out 2>&1

echo "3) PROJWFC"
projwfc.x -in projwfc.in > projwfc.out 2>&1

echo "4) NSCF"
pw.x -in nscf.in > nscf.out 2>&1

echo "5) PW2WANNIER90"
pw2wannier90.x -in pw2wan.in > pw2wan.out 2>&1

echo "6) WANNIER90"
wannier90 Fe_bcc > wannier90.out 2>&1 || true

echo "Flow finished (check *.out files)."
