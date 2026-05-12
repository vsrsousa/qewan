# qewan

qewan is a small Python CLI to generate Quantum ESPRESSO (`pw.x`) and Wannier90 (`.win`) input files from a CIF file, following the project's conventions.

Features
- Generate SCF, bands, NSCF, Wannier90 `.win` files, and helper inputs for `pw2wannier` and `projwfc`.
- `--workdir` + `--seedname` layout: outputs are organized as `workdir/seedname/<scf,bands,nscf,wan,projwfc>` by default.
- Spin-aware outputs: when `--nspin 2` the tool writes per-spin files with `_up` and `_dn` suffixes and includes `spin_component` where applicable.
- Wannier90 `.win` files: include `unit_cell_cart` (Å), `begin atoms_frac`, `mp_grid = N N N` when known, and a `begin kpoints` block with coordinates (weights are stored in the NSCF `kpoints_list.txt`).
- `pw2wannier` inputs: written using the `&INPUTPP` namelist. If a `pw.nscf.in` inside the provided NSCF folder contains a literal `outdir = '...'` in its `&CONTROL`, that string is used verbatim for `outdir` in the `pw2wannier` input; otherwise the absolute NSCF path is used.
- Pseudopotential handling: when absolute UPF paths are provided, `pseudo_dir` is set to the common directory and `ATOMIC_SPECIES` lists only UPF basenames (no full paths).
- `ecutrho` inference: derived from `ecutwfc` and the pseudo type detected in the UPF header (PAW/USPP → multiplier 8, NC → multiplier 4).

Quick install

```bash
pip install -r requirements.txt
python -m pip install -e .
```

**Ambiente de desenvolvimento**

- **Ativar venv antes de rodar comandos Python**: execute `source .venv/bin/activate` na raiz do repositório. Isto garante que o intérprete e dependências corretas sejam usados.


Basic examples

Note: the project prefers `--workdir` and `--seedname` to produce a reproducible folder layout. If you previously used `--outdir`, it will still work but `--workdir` is recommended.

```bash
# input CIF
CIF=tests/cifs/FeO.cif

# generate SCF and NSCF (explicit MP grid)
python -m qewan.cli run-scf $CIF --workdir demo_run --seedname feo --kpoints "6 6 6"
python -m qewan.cli run-nscf $CIF --workdir demo_run --seedname feo --kpoints "6 6 6" --save-kmesh

# generate Wannier inputs (the command will read demo_run/feo/nscf/kpoints_list.txt if available)
python -m qewan.cli run-wan $CIF --workdir demo_run --seedname feo --nspin 1

# generate pw2wannier inputs (these will be written inside the wan folder)
python -m qewan.cli run-pw2wannier $CIF --workdir demo_run --seedname feo --nscf-outdir demo_run/feo/nscf --nspin 2
```

Pseudopotentials

- Use `--pseudos` to map elements to UPF filenames (comma-separated `Element:File.UPF`).
- If absolute UPF paths are provided, `qewan` will set `pseudo_dir` to the common directory containing the UPF files; `ATOMIC_SPECIES` always uses basenames only.

Wannier90 and k-points

- If the NSCF step saved an explicit k-point list (`kpoints_list.txt`) or MP grid (`kmesh.txt`), the `.win` generator will include `mp_grid = N N N` and a `begin kpoints` block with the explicit coordinates. The weights remain in `kpoints_list.txt`.
- You can pass `--bands-file` to `run-wan` to add a `bands_file = <path>` hint in the `.win` header.

Files of interest

- `qewan/io.py` — core input generators.
- `qewan/cli.py` — Click CLI entrypoints.

Notes

- The CLI exposes flags for `--write-mmn/--no-write-mmn`, `--write-amn`, `--write-unk`, and `--auto-projections/--no-auto-projections` that control `pw2wannier` and `.win` output options.
- When `--nspin 2` is used, per-spin `.win` and `pw2wannier` inputs are produced with `_up`/`_dn` suffixes and the appropriate `spin_component` fields.

Contributions & further work

- If you want the README restructured (API reference, extended examples, CI snippets or unit tests), I can open a PR with a longer document.

---

