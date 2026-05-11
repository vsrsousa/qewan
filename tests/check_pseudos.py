#!/usr/bin/env python3
from qewan.io import read_cif, atoms_to_pw_input
from pathlib import Path
import sys

cif = 'tests/cifs/SrMnO3.cif'
pseudo_dir = '/home/vinicius/scratch/data/pseudos/SSSP_1.3.0_PBE_efficiency'

atoms = read_cif(cif)

inp = atoms_to_pw_input(atoms, calculation='nscf', kpoints='6 6 6', pseudos=None, pseudo_dir=pseudo_dir, nspin=1)

lines = inp.splitlines()
# print pseudo_dir line and ATOMIC_SPECIES block
for i, ln in enumerate(lines):
    if 'pseudo_dir' in ln:
        print(ln)
    if ln.strip() == 'ATOMIC_SPECIES':
        print('\nATOMIC_SPECIES')
        # print following lines until ATOMIC_POSITIONS
        for j in range(i+1, len(lines)):
            if lines[j].strip().startswith('ATOMIC_POSITIONS'):
                break
            if lines[j].strip():
                print(lines[j])
        break

# also print a short snippet of the full generated input for manual inspection
print('\n--- INPUT SNIPPET ---\n')
print('\n'.join(lines[:60]))
