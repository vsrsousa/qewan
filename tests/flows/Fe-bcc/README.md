Fe-bcc flow inputs

This folder contains example input files to run a test flow for bcc Fe:

- scf.in        : self-consistent calculation
- bands.in      : band structure calculation (along a short high-symmetry path)
- projwfc.in    : projection analysis input for `projwfc.x`
- nscf.in       : non-self-consistent calculation on dense k-grid
- pw2wan.in     : input for `pw2wannier90.x` (writes the necessary data for Wannier90)
- Fe_bcc.win    : minimal `wannier90.win` file for generating Wannier functions
- run_flow.sh   : simple sequential runner (adjust executables/pseudo paths)

Notes / placeholders
- Pseudopotential: replace `Fe.pseudo.UPF` in the inputs with an appropriate pseudopotential file present in `pseudo/`.
- `outdir` is set to `./out` by default; ensure this directory exists or change the path.
- The inputs are basic templates meant to be adapted to your pseudopotentials, functional, and convergence needs.

Usage example:
```bash
bash run_flow.sh
```
