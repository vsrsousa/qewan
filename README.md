# qewan

qewan is a small Python CLI/package to convert a CIF into a set of inputs for Quantum ESPRESSO (pw.x) and Wannier90 (.win), following project conventions.

Key behaviors
- pw2wannier input namelist: the generator writes `&INPUTPP` (not `&pw2wannier`).
- pw2wannier `outdir` resolution: if a `pw.nscf.in` file exists inside the provided NSCF folder, the literal `outdir = '...'` value from its `&CONTROL` block is used verbatim. Otherwise the absolute NSCF path is written.
- Spin support: use `--nspin`/`nspin` options; when `nspin=2` the CLI writes per-spin outputs with `_up` and `_dn` suffixes and includes `spin_component = 'up'/'down'` when appropriate.
- Wannier `.win`: `mp_grid = N N N` (when known) is written immediately above `begin kpoints` and `begin kpoints` contains coordinates only (no weights). The `.win` also contains `unit_cell_cart` in angstrom and `begin atoms_frac` when fractional coordinates are used.

Recent behavior notes
- For `calculation='bands'` the package now generates a k-path in the Brillouin zone by default using `seekpath`. `seekpath` is required for bands generation (listed in `requirements.txt`) and the generator will raise a clear error if it's missing.
- Pseudopotential handling: when the user provides absolute UPF paths via `--pseudos`, the generator will set `pseudo_dir` in the produced `pw.x` inputs to the common directory containing those UPF files (so you do not need to copy pseudos into the outdir). `ATOMIC_SPECIES` always lists only the UPF basenames (no full paths).
- `ecutrho` is inferred automatically from `ecutwfc` and the UPF header pseudo type: PAW/USPP use multiplier 8, norm-conserving use multiplier 4. The detector reads the UPF header (not filenames) to determine pseudo type.


Quick example (run the demo test performed by the maintainer):

```bash
# CIF to test
CIF=/Users/vinicius/Documents/cifs/GdNi4Si.cif

# run full flow into demo_run/gdni4si
python -m qewan.cli run-scf $CIF --outdir demo_run/gdni4si/scf --kpoints "6 6 6"
python -m qewan.cli run-bands $CIF --outdir demo_run/gdni4si/bands --kpath --npoints 60
python -m qewan.cli run-projwfc $CIF --outdir demo_run/gdni4si/projwfc --nscf-outdir demo_run/gdni4si/nscf
python -m qewan.cli run-nscf $CIF --outdir demo_run/gdni4si/nscf --kpoints "4 4 4"
python -m qewan.cli run-wan $CIF --outdir demo_run/gdni4si/wan --nscf-dir demo_run/gdni4si/nscf --seedname qewan_test
python -m qewan.cli run-pw2wannier $CIF --outdir demo_run/gdni4si/pw2wann --wan-dir demo_run/gdni4si/wan --nscf-outdir demo_run/gdni4si/nscf --seedname qewan_test
```

Pseudopotential handling
- You can control the `pseudo_dir` string written into the QE `&CONTROL` block using `--pseudo-dir`.
- You can provide an explicit per-element mapping using `--pseudos`, for example:

```bash
python -m qewan.cli run-scf /path/to/str.cif --outdir scf --pseudo-dir ./mypseudos --pseudos "Gd:Gd.UPF,Co:Co.UPF"
```

If `--pseudos` is omitted, the generator defaults to `<Element>.UPF` for each species and writes the `pseudo_dir` string you supplied (default `./pseudos`).

Notes
- If you want `pw2wannier` to use an exact `outdir` string (for example `./tmp`) make sure `pw.nscf.in` in the NSCF folder contains that `outdir` assignment.
- The CLI exposes flags for `--write-mmn/--no-write-mmn`, `--write-amn`, `--write-unk`, and `--auto-projections/--no-auto-projections`.

Files of interest
- `qewan/io.py` — core generators.
- `qewan/cli.py` — Click-based CLI.

If you want me to add a unit test or CI snippet covering the `&INPUTPP` header and `outdir` extraction, tell me and I'll add it.

New CLI flags (added)
- `run-wan` now accepts `--wannier-plot/--no-wannier-plot` (default: `--no-wannier-plot`).
- `run-wan` now accepts `--write-xyz/--no-write-xyz` (default: `--write-xyz`).
- `run-wan` now accepts `--write-hr/--no-write-hr` (default: `--write-hr`).

These flags control the corresponding entries in the generated `.win` file (`wannier_plot`, `write_xyz`, `write_hr`) and default to the values used by the maintainer during testing.

Per-spin behavior
- When running with `--nspin 2`, `run-wan` generates separate `.win` files for each spin channel (suffix `_up` and `_dn`).
- `run-pw2wannier` when used with `--nspin 2` creates `pw2wannier_up.in` and `pw2wannier_dn.in` colocated with the `.win` files and includes `spin_component = 'up'/'down'`.

Quick example showing the new flags:

```bash
# generate per-spin .win with custom output flags
python -m qewan.cli run-wan /path/to/structure.cif --outdir demo_run/wan --nscf-dir demo_run/nscf --seedname myseed --nspin 2 --wannier-plot --no-write-xyz --no-write-hr
```

README updates: these notes summarize recent changes; if you'd like the README restructured (API, examples, CLI reference), I can open a PR with a longer doc.
# qewan

Ferramenta para gerar arquivos de entrada do Quantum ESPRESSO (pw.x) e Wannier90 a partir de um arquivo CIF. Suporta geração de SCF, bands, nscf e arquivos de entrada do Wannier90, além de templates de submissão para cluster (SLURM).

Instalação rápida:

```bash
pip install -r requirements.txt
python -m pip install -e .
```

Exemplo de uso (CLI):

```bash
qewan run_scf structure.cif --outdir run_scf
qewan run_bands structure.cif --outdir run_bands
qewan run_nscf structure.cif --outdir run_nscf
qewan run_wan structure.cif --outdir run_wan
```

Configurar pseudopotenciais: edite o dicionário `pseudos` nos templates ou passe via argumentos.

Padrões de nomes de arquivos gerados:

- `projwfc` defaults:
	- `filproj` → `<prefix>.bands.dat.proj` (ex.: `qewan.bands.dat.proj` quando `prefix=qewan`)
	- `fillowdin` → `lowdin.txt`

- NSCF outputs when using an explicit MP grid (`Nx Ny Nz`):
	- `kpoints_list.txt` → contains lines `kx ky kz weight` (weights normalized to 1/num_points)
	- `kmesh.txt` → contains the original MP string, e.g. `8 8 8` (used by the .win generator)

Wannier (.win) notes:

- When an explicit NSCF k-point list or an `mp_grid` is available, the generated `.win` will include `mp_grid = N N N` immediately above the `begin kpoints` block.
- The `.win` k-point block contains coordinates only (no weights) — weights must be read from the NSCF `kpoints_list.txt`.
- You can pass a `--bands-file` to `run_wan` and it will write a `bands_file = <path>` hint in the `.win` header to point to your bands output.

Example end-to-end sequence (recommended):

```bash
python -m qewan.cli run_scf structure.cif --outdir demo_run/scf --kpoints "6 6 6"
python -m qewan.cli run_bands structure.cif --outdir demo_run/bands --kpath --npoints 60
python -m qewan.cli run_projwfc structure.cif --outdir demo_run/projwfc --prefix qewan
python -m qewan.cli run_nscf structure.cif --outdir demo_run/nscf --kpoints "8 8 8"
python -m qewan.cli run_wan structure.cif --outdir demo_run/wann --nscf-dir demo_run/nscf --bands-file demo_run/bands/pw.bands.out
```

Se quiser, posso ajustar os nomes de saída ou o formato dos arquivos para compatibilidade com um cluster específico — diga como você prefere os prefixes/paths.

Coordenadas e célula convencional
--------------------------------

- `ATOMIC_POSITIONS` em todos os inputs `pw.x` (SCF/BANDS/NSCF) são escritas em coordenadas fracionárias (`crystal`). Isso é independente do tipo de célula extraída pelo ASE: as posições sempre aparecem como fracionárias nos arquivos gerados.
- Por padrão a CLI tenta reduzir a estrutura para a célula primitiva via `spglib` para evitar átomos equivalentes repetidos na célula (isso reduz o número de posições escritas). Se você preferir usar a célula convencional gerada pelo ASE (mais átomos), passe a flag `--conventional-cell` nos comandos que aceitam CIF.
- Se `spglib` não estiver instalado a redução é ignorada; instale com `pip install spglib` para ativar a redução automática.
