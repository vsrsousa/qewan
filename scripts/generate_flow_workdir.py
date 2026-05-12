#!/usr/bin/env python3
"""Generate a full QE+Wannier input flow under <workdir>/<seedname>/<step>/

This script uses `qewan.io` generators and places outputs in the
conventional subfolders: `run_scf`, `run_nscf`, `run_bands`, `run_wan`,
`run_pw2wannier`, `run_projwfc` under the provided `workdir/seedname`.

Usage:
  python scripts/generate_flow_workdir.py path/to.cif --workdir tests/flows --seedname GdCu --pseudo-dir /abs/path/to/pseudos
"""
import argparse
import os
from qewan.io import read_cif, atoms_to_pw_input, write_pw_input, generate_wannier_win, generate_pw2wannier_input, generate_projwfc_input


def ensure_dir(p):
    os.makedirs(p, exist_ok=True)


def write_pw(atoms, calculation, outdir, **kwargs):
    ensure_dir(outdir)
    inp = atoms_to_pw_input(atoms, calculation=calculation, **kwargs)
    path = os.path.join(outdir, 'pw.{}.in'.format('scf' if calculation == 'scf' else calculation if calculation != 'nscf' else 'nscf'))
    write_pw_input(path, inp)
    print(f'Wrote {calculation} input: {path}')
    return path


def main():
    p = argparse.ArgumentParser()
    p.add_argument('cif', help='Path to CIF file')
    p.add_argument('--workdir', required=True, help='Top-level workdir where <workdir>/<seedname> will be created')
    p.add_argument('--seedname', required=True, help='Seedname used as folder basename')
    p.add_argument('--pseudo-dir', default='./pseudos', help='Pseudo directory string to write in inputs (passed verbatim to generators)')
    p.add_argument('--nspin', type=int, default=1)
    p.add_argument('--kpoints', default='6 6 6')
    args = p.parse_args()

    atoms = read_cif(args.cif)
    base = os.path.join(args.workdir, args.seedname)

    # SCF
    scf_dir = os.path.join(base, 'run_scf')
    ensure_dir(scf_dir)
    write_pw(atoms, 'scf', scf_dir, kpoints=args.kpoints, pseudos=None, pseudo_dir=args.pseudo_dir, nspin=args.nspin)

    # NSCF
    nscf_dir = os.path.join(base, 'run_nscf')
    ensure_dir(nscf_dir)
    write_pw(atoms, 'nscf', nscf_dir, kpoints=args.kpoints, pseudos=None, pseudo_dir=args.pseudo_dir, nspin=args.nspin)

    # Bands
    bands_dir = os.path.join(base, 'run_bands')
    ensure_dir(bands_dir)
    write_pw(atoms, 'bands', bands_dir, kpoints=args.kpoints, pseudos=None, pseudo_dir=args.pseudo_dir, nspin=args.nspin)

    # Wannier .win
    wan_dir = os.path.join(base, 'run_wan')
    ensure_dir(wan_dir)
    win = generate_wannier_win(atoms, wan_dir, seedname=args.seedname)
    print(f'Wrote wannier90 .win: {win}')

    # pw2wannier
    pw2_dir = os.path.join(base, 'run_pw2wannier')
    ensure_dir(pw2_dir)
    p2 = generate_pw2wannier_input(pw2_dir, prefix=args.seedname, nscf_outdir=nscf_dir, seedname=args.seedname)
    print(f'Wrote pw2wannier input: {p2}')

    # projwfc
    proj_dir = os.path.join(base, 'run_projwfc')
    ensure_dir(proj_dir)
    proj = generate_projwfc_input(proj_dir, prefix=args.seedname, outdir_pw=nscf_dir)
    print(f'Wrote projwfc input: {proj}')


if __name__ == '__main__':
    main()
