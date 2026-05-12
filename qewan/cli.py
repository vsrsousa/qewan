import click
import os
from qewan.io import read_cif, atoms_to_pw_input, write_pw_input, generate_wannier_win, generate_projwfc_input, reduce_atoms_by_symmetry
from qewan.io import generate_pw2wannier_input
from qewan.cluster import write_slurm_script, submit_slurm
import subprocess


def _canonical_wan_dir(requested_outdir: str, explicit_wan: str | None = None) -> str:
    """Return a canonical directory to store Wannier-related files.

    Rules:
    - if `explicit_wan` is provided, return it
    - if the requested_outdir basename contains 'pw2wann' or 'pw2wannier',
      replace that token with 'wan' (e.g. 'qewan_pw2wann' -> 'qewan_wan')
    - if requested_outdir already ends with '_wan', return it
    - otherwise append '_wan' to the basename
    """
    if explicit_wan:
        return explicit_wan
    parent = os.path.dirname(requested_outdir) or '.'
    base = os.path.basename(requested_outdir)
    lb = base.lower()
    if 'pw2wann' in lb:
        # preserve original case around the token where possible
        new_base = base.replace('pw2wann', 'wan').replace('pw2wannier', 'wan')
        return os.path.join(parent, new_base)
    if base.endswith('_wan') or base.lower().endswith('wan'):
        return requested_outdir
    return os.path.join(parent, base + '_wan')


@click.group()
def cli():
    """qewan CLI: generate QE + Wannier inputs from CIF"""
    pass


@cli.command()
@click.argument('cif', type=click.Path(exists=True))
@click.option('--workdir', default='.', help='Top-level workdir to place results under <workdir>/<seedname>/{scf,nscf,bands,projwfc,wan}')
@click.option('--seedname', default='qewan', help='Seedname used to create subfolder under --workdir (default: qewan)')
@click.option('--kpoints', default='6 6 6')
@click.option('--pseudo-dir', default='./pseudos', help='Pseudo directory string to write in &CONTROL (will be placed verbatim)')
@click.option('--pseudos', default=None, help="Comma-separated element:pseudo pairs, e.g. 'Gd:Gd.UPF,Co:Co.UPF'. If omitted, defaults to <element>.UPF")
@click.option('--conventional-cell/--no-conventional-cell', default=False, help='Use CIF conventional cell as parsed by ASE (default: reduce by symmetry)')
@click.option('--nspin', default=1, type=int, help='Set `nspin` in &SYSTEM (1 or 2)')
@click.option('--starting-magnetization', default=None, help='Comma-separated species:mag pairs, e.g. "Fe:0.5,Co:0.2"')
@click.option('--occupations', default='smearing', help="Set occupations in &SYSTEM (e.g. 'smearing')")
@click.option('--smearing', default='mv', help="Set smearing type in &SYSTEM (e.g. 'mv')")
@click.option('--degauss', default=0.02, type=float, help='Set degauss value in &SYSTEM (e.g. 0.02)')
@click.option('--ecutrho', default=None, type=float, help='Set ecutrho in &SYSTEM (e.g. 320.0)')
@click.option('--ecutwfc', default=40.0, type=float, help='Set ecutwfc in &SYSTEM (e.g. 40.0)')
def run_scf(cif, workdir, seedname, kpoints, pseudo_dir, pseudos, conventional_cell, nspin, starting_magnetization, occupations, smearing, degauss, ecutrho, ecutwfc):
    atoms = read_cif(cif)
    # reduce by symmetry unless user requests conventional cell
    if not conventional_cell:
        try:
            reduced, _ = reduce_atoms_by_symmetry(atoms, to_primitive=True)
            if reduced is not None:
                atoms = reduced
                click.echo(f'Using symmetry-reduced primitive cell with {len(atoms)} atoms')
        except Exception:
            pass
    # parse pseudos mapping (element:pseudo) and starting magnetization string passed by click
    orig_pseudos = None
    if pseudos:
        orig_pseudos = {}
        for pair in pseudos.split(','):
            if ':' in pair:
                k, v = pair.split(':', 1)
                orig_pseudos[k.strip()] = v.strip()
    # prepare mapping to pass to atoms_to_pw_input: use basenames when desired
    pseudos_map = None
    if orig_pseudos:
        # pass full paths to the generator so it can inspect UPF headers;
        # atoms_to_pw_input will write only basenames in ATOMIC_SPECIES
        pseudos_map = {k: v for k, v in orig_pseudos.items()}
    smap = None
    if starting_magnetization:
        smap = {}
        for pair in starting_magnetization.split(','):
            if ':' in pair:
                k, v = pair.split(':', 1)
                try:
                    smap[k.strip()] = float(v)
                except Exception:
                    pass
    # If user supplied explicit pseudo paths, copy them into outdir/pseudo_dir
    # BEFORE generating the input so the generator can inspect the actual files
    try:
        import shutil
        if orig_pseudos:
            target = os.path.join(outdir, pseudo_dir)
            os.makedirs(target, exist_ok=True)
            for elm, src in orig_pseudos.items():
                if os.path.exists(src):
                    dst = os.path.join(target, os.path.basename(src))
                    shutil.copy(src, dst)
                    # pass the copied full path to the generator so it reads the copied file
                    if pseudos_map is None:
                        pseudos_map = {}
                    pseudos_map[elm] = dst
                    click.echo(f'Copied pseudo {src} -> {dst}')
                else:
                    # keep whatever was provided (basename or path)
                    if pseudos_map is None:
                        pseudos_map = {}
                    pseudos_map[elm] = src
    except Exception:
        pass

    inp = atoms_to_pw_input(atoms, calculation='scf', kpoints=kpoints, pseudos=pseudos_map, pseudo_dir=pseudo_dir, nspin=nspin, starting_magnetization=smap, occupations=occupations, smearing=smearing, degauss=degauss, ecutrho=ecutrho, ecutwfc=ecutwfc)
    # compute canonical outdir: <workdir>/<seedname>/scf
    outdir = os.path.join(workdir or '.', seedname or 'qewan', 'scf')
    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, 'pw.scf.in')
    write_pw_input(path, inp)
    click.echo(f'Wrote SCF input: {path}')


@cli.command()
@click.argument('cif', type=click.Path(exists=True))
@click.option('--workdir', default='.', help='Top-level workdir to place results under <workdir>/<seedname>/{scf,nscf,bands,projwfc,wan}')
@click.option('--seedname', default='qewan', help='Seedname used to create subfolder under --workdir (default: qewan)')
@click.option('--kpoints', default='6 6 6')
@click.option('--kpath/--no-kpath', default=False, help='Use seekpath to generate band path')
@click.option('--npoints', default=40, help='Approx points along total path')
@click.option('--slurm/--no-slurm', default=False, help='Write SLURM submission script')
@click.option('--submit/--no-submit', default=False, help='Submit the SLURM job after creating script')
@click.option('--ntasks', default=16)
@click.option('--time', default='04:00:00')
@click.option('--conventional-cell/--no-conventional-cell', default=False, help='Use CIF conventional cell as parsed by ASE (default: reduce by symmetry)')
@click.option('--nspin', default=1, type=int, help='Set `nspin` in &SYSTEM (1 or 2)')
@click.option('--starting-magnetization', default=None, help='Comma-separated species:mag pairs, e.g. "Fe:0.5,Co:0.2"')
@click.option('--occupations', default='smearing', help="Set occupations in &SYSTEM (e.g. 'smearing')")
@click.option('--smearing', default='mv', help="Set smearing type in &SYSTEM (e.g. 'mv')")
@click.option('--degauss', default=0.02, type=float, help='Set degauss value in &SYSTEM (e.g. 0.02)')
@click.option('--ecutrho', default=None, type=float, help='Set ecutrho in &SYSTEM (e.g. 320.0)')
@click.option('--pseudo-dir', default='./pseudos', help='Pseudo directory string to write in &CONTROL (will be placed verbatim)')
@click.option('--pseudos', default=None, help="Comma-separated element:pseudo pairs, e.g. 'Gd:Gd.UPF,Co:Co.UPF'. If omitted, defaults to <element>.UPF")
@click.option('--ecutwfc', default=40.0, type=float, help='Set ecutwfc in &SYSTEM (e.g. 40.0)')
def run_bands(cif, workdir, seedname, kpoints, kpath, npoints, slurm, submit, ntasks, time, conventional_cell, nspin, starting_magnetization, occupations, smearing, degauss, ecutrho, pseudo_dir, pseudos, ecutwfc):
    atoms = read_cif(cif)
    if not conventional_cell:
        try:
            reduced, _ = reduce_atoms_by_symmetry(atoms, to_primitive=True)
            if reduced is not None:
                atoms = reduced
                click.echo(f'Using symmetry-reduced primitive cell with {len(atoms)} atoms')
        except Exception:
            pass
    orig_pseudos = None
    if 'pseudos' in locals() and pseudos is not None:
        orig_pseudos = {}
        for pair in pseudos.split(','):
            if ':' in pair:
                k, v = pair.split(':', 1)
                orig_pseudos[k.strip()] = v.strip()
    # prepare mapping to pass to atoms_to_pw_input: use basenames when desired
    pseudos_map = None
    if orig_pseudos:
        pseudos_map = {k: v for k, v in orig_pseudos.items()}

    smap = None
    if starting_magnetization:
        smap = {}
        for pair in starting_magnetization.split(','):
            if ':' in pair:
                k, v = pair.split(':', 1)
                try:
                    smap[k.strip()] = float(v)
                except Exception:
                    pass
    # copy pseudos into outdir/pseudo_dir before generating input so detector reads them
    try:
        import shutil
        if orig_pseudos:
            target = os.path.join(outdir, pseudo_dir)
            os.makedirs(target, exist_ok=True)
            for elm, src in orig_pseudos.items():
                if os.path.exists(src):
                    dst = os.path.join(target, os.path.basename(src))
                    shutil.copy(src, dst)
                    if pseudos_map is None:
                        pseudos_map = {}
                    pseudos_map[elm] = dst
                    click.echo(f'Copied pseudo {src} -> {dst}')
                else:
                    if pseudos_map is None:
                        pseudos_map = {}
                    pseudos_map[elm] = src
    except Exception:
        pass

    inp = atoms_to_pw_input(atoms, calculation='bands', kpoints=kpoints, kpath=kpath, npoints=npoints, pseudos=pseudos_map, pseudo_dir=pseudo_dir, nspin=nspin, starting_magnetization=smap, occupations=occupations, smearing=smearing, degauss=degauss, ecutrho=ecutrho, ecutwfc=ecutwfc)
    outdir = os.path.join(workdir or '.', seedname or 'qewan', 'bands')
    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, 'pw.bands.in')
    write_pw_input(path, inp)
    click.echo(f'Wrote bands input: {path}')
    if slurm:
        cmd = f"mpirun -np $SLURM_NTASKS pw.x -inp {os.path.basename(path)} > {os.path.basename(path)}.out"
        script = write_slurm_script(outdir, cmd, ntasks=ntasks, time=time, job_name='qewan_bands')
        click.echo(f'Wrote SLURM script: {script}')
        if submit:
            res = submit_slurm(script)
            click.echo(f'sbatch output: {res.stdout.strip()} {res.stderr.strip()}')


@cli.command()
@click.argument('cif', type=click.Path(exists=True))
@click.option('--workdir', default='.', help='Top-level workdir to place results under <workdir>/<seedname>/{scf,nscf,bands,projwfc,wan}')
@click.option('--seedname', default='qewan', help='Seedname used to create subfolder under --workdir (default: qewan)')
@click.option('--prefix', default='qewan')
@click.option('--outdir-pw', default='./tmp')
@click.option('--lsym/--no-lsym', default=False)
@click.option('--lwrite-overlaps/--no-lwrite-overlaps', default=False)
@click.option('--filproj', default=None)
@click.option('--fillowdin', default='lowdin.txt')
def run_projwfc(cif, workdir, seedname, prefix, outdir_pw, lsym, lwrite_overlaps, filproj, fillowdin):
    # cif argument kept for interface consistency (not used)
    outdir = os.path.join(workdir or '.', seedname or 'qewan', 'projwfc')
    os.makedirs(outdir, exist_ok=True)
    path = generate_projwfc_input(outdir, prefix=prefix, outdir_pw=outdir_pw, lsym=lsym, lwrite_overlaps=lwrite_overlaps, filproj=filproj, fillowdin=fillowdin)
    click.echo(f'Wrote projwfc input: {path}')


@cli.command()
@click.argument('cif', type=click.Path(exists=True))
@click.option('--prefix', default='qewan')
@click.option('--seedname', default='qewan', help='Seedname to use in pw2wannier input (overrides prefix)')
@click.option('--workdir', default='.', help='Top-level workdir to place results under <workdir>/<seedname>/{scf,nscf,bands,projwfc,wan}')
@click.option('--nspin', default=1, type=int, help='Number of spin channels (1 or 2). When 2, generate per-spin inputs')
@click.option('--nscf-outdir', default='./tmp', help='Path to nscf outdir where pw.x wrote data')
@click.option('--write-mmn/--no-write-mmn', default=True, help='Set write_mmn in pw2wannier input')
@click.option('--write-amn/--no-write-amn', default=True, help='Set write_amn in pw2wannier input')
@click.option('--write-unk/--no-write-unk', default=False, help='Set write_unk in pw2wannier input')
@click.option('--auto-projections/--no-auto-projections', default=True, help='Include atom_proj in pw2wannier input when using automatic projections')
@click.option('--wan-dir', default=None, type=click.Path(exists=False), help='Directory where Wannier .win files live; pw2wannier inputs/outputs will be written here')
def run_pw2wannier(cif, prefix, seedname, workdir, nscf_outdir, write_mmn, write_amn, write_unk, auto_projections, nspin, wan_dir):
    """Generate a pw2wannier90 input file pointing to the NSCF outdir."""
    # If spin-polarized, generate separate inputs per spin channel in subdirs
    seed = seedname or prefix
    # determine directory where pw2wannier input should be written
    default_wan = os.path.join(workdir or '.', seedname or 'qewan', 'wan')
    write_dir = wan_dir or default_wan
    os.makedirs(write_dir, exist_ok=True)
    if int(nspin) == 2:
        for ch in ('up', 'dn'):
            # generator expects 'down' string (not 'dn') for spin_component
            spin_comp = 'up' if ch == 'up' else 'down'
            p = generate_pw2wannier_input(
                write_dir,
                prefix=prefix,
                nscf_outdir=nscf_outdir,
                seedname=f"{seed}_{ch}",
                write_mmn=write_mmn,
                write_amn=write_amn,
                write_unk=write_unk,
                auto_projections=auto_projections,
                spin_component=spin_comp,
            )
            # rename to per-spin filename to avoid overwrite and avoid using seedname in filename
            dst = os.path.join(write_dir, f'pw2wannier_{ch}.in')
            try:
                if os.path.exists(dst):
                    os.remove(dst)
                os.rename(os.path.join(write_dir, 'pw2wannier.in'), dst)
                click.echo(f'Wrote pw2wannier input: {dst}')
            except Exception:
                # fallback: report original path
                click.echo(f'Wrote pw2wannier input: {p}')
        else:
            path = generate_pw2wannier_input(
                write_dir,
                prefix=prefix,
                nscf_outdir=nscf_outdir,
                seedname=seed,
                write_mmn=write_mmn,
                write_amn=write_amn,
                write_unk=write_unk,
                auto_projections=auto_projections,
            )
            click.echo(f'Wrote pw2wannier input: {path}')


@cli.command()
@click.argument('cif', type=click.Path(exists=True))
@click.option('--kpoints', default='6 6 6')
@click.option('--save-kmesh/--no-save-kmesh', default=True, help='Save the nscf k-mesh to kmesh.txt in outdir')
@click.option('--conventional-cell/--no-conventional-cell', default=False, help='Use CIF conventional cell as parsed by ASE (default: reduce by symmetry)')
@click.option('--workdir', default='.', help='Top-level workdir to place results under <workdir>/<seedname>/{scf,nscf,bands,projwfc,wan}')
@click.option('--seedname', default='qewan', help='Seedname used to create subfolder under --workdir (default: qewan)')
@click.option('--nspin', default=1, type=int, help='Set `nspin` in &SYSTEM (1 or 2)')
@click.option('--starting-magnetization', default=None, help='Comma-separated species:mag pairs, e.g. "Fe:0.5,Co:0.2"')
@click.option('--occupations', default='smearing', help="Set occupations in &SYSTEM (e.g. 'smearing')")
@click.option('--smearing', default='mv', help="Set smearing type in &SYSTEM (e.g. 'mv')")
@click.option('--degauss', default=0.02, type=float, help='Set degauss value in &SYSTEM (e.g. 0.02)')
@click.option('--ecutrho', default=None, type=float, help='Set ecutrho in &SYSTEM (e.g. 320.0)')
@click.option('--ecutwfc', default=40.0, type=float, help='Set ecutwfc in &SYSTEM (e.g. 40.0)')
@click.option('--pseudo-dir', default='./pseudos', help='Pseudo directory string to write in &CONTROL (will be placed verbatim)')
@click.option('--pseudos', default=None, help="Comma-separated element:pseudo pairs, e.g. 'Gd:Gd.UPF,Co:Co.UPF'. If omitted, defaults to <element>.UPF")
def run_nscf(cif, kpoints, save_kmesh, conventional_cell, nspin, starting_magnetization, occupations, smearing, degauss, ecutrho, pseudo_dir, pseudos, ecutwfc, workdir='.', seedname='qewan'):
    atoms = read_cif(cif)
    if not conventional_cell:
        try:
            reduced, _ = reduce_atoms_by_symmetry(atoms, to_primitive=True)
            if reduced is not None:
                atoms = reduced
                click.echo(f'Using symmetry-reduced primitive cell with {len(atoms)} atoms')
        except Exception:
            pass
    # If kpoints is a triple like 'Nx Ny Nz', generate explicit list and write it
    if kpoints and len(kpoints.split()) == 3:
        from qewan.io import generate_mp_kpoints_list
        pts, block = generate_mp_kpoints_list(atoms, kpoints)
        orig_pseudos = None
        if 'pseudos' in locals() and pseudos is not None:
            orig_pseudos = {}
            for pair in pseudos.split(','):
                if ':' in pair:
                    k, v = pair.split(':', 1)
                    orig_pseudos[k.strip()] = v.strip()
        pseudos_map = None
        if orig_pseudos:
            pseudos_map = {k: v for k, v in orig_pseudos.items()}

        smap = None
        if starting_magnetization:
            smap = {}
            for pair in starting_magnetization.split(','):
                if ':' in pair:
                    k, v = pair.split(':', 1)
                    try:
                        smap[k.strip()] = float(v)
                    except Exception:
                        pass
        # copy pseudos into outdir/pseudo_dir before generating input so detector reads them
        try:
            import shutil
            if orig_pseudos:
                target = os.path.join(outdir, pseudo_dir)
                os.makedirs(target, exist_ok=True)
                for elm, src in orig_pseudos.items():
                    if os.path.exists(src):
                        dst = os.path.join(target, os.path.basename(src))
                        shutil.copy(src, dst)
                        if pseudos_map is None:
                            pseudos_map = {}
                        pseudos_map[elm] = dst
                        click.echo(f'Copied pseudo {src} -> {dst}')
                    else:
                        if pseudos_map is None:
                            pseudos_map = {}
                        pseudos_map[elm] = src
        except Exception:
            pass

        outdir = os.path.join(workdir or '.', seedname or 'qewan', 'nscf')
        inp = atoms_to_pw_input(atoms, calculation='nscf', kpoints=block, pseudos=pseudos_map, pseudo_dir=pseudo_dir, nspin=nspin, starting_magnetization=smap, occupations=occupations, smearing=smearing, degauss=degauss, ecutrho=ecutrho, ecutwfc=ecutwfc)
        path = os.path.join(outdir, 'pw.nscf.in')
        os.makedirs(outdir, exist_ok=True)
        write_pw_input(path, inp)
        click.echo(f'Wrote nscf input (explicit k-list): {path}')
        if save_kmesh:
            os.makedirs(outdir, exist_ok=True)
            kmesh_path = os.path.join(outdir, 'kpoints_list.txt')
            total = len(pts)
            weight = 1.0 / float(total) if total > 0 else 1.0
            with open(kmesh_path, 'w') as f:
                for p in pts:
                    f.write(f"{p[0]:.8f} {p[1]:.8f} {p[2]:.8f} {weight:.8f}\n")
            # also save the original mp grid string so callers can read mp_grid
            kmesh_info = os.path.join(outdir, 'kmesh.txt')
            with open(kmesh_info, 'w') as f:
                f.write(kpoints.strip() + "\n")
            click.echo(f'Wrote explicit k-points list: {kmesh_path}')
    else:
        smap = None
        if starting_magnetization:
            smap = {}
            for pair in starting_magnetization.split(','):
                if ':' in pair:
                    k, v = pair.split(':', 1)
                    try:
                        smap[k.strip()] = float(v)
                    except Exception:
                        pass
        inp = atoms_to_pw_input(atoms, calculation='nscf', kpoints=kpoints, pseudos=pseudos_map, pseudo_dir=pseudo_dir, nspin=nspin, starting_magnetization=smap, occupations=occupations, smearing=smearing, degauss=degauss, ecutrho=ecutrho, ecutwfc=ecutwfc)
        outdir = os.path.join(workdir or '.', seedname or 'qewan', 'nscf')
        os.makedirs(outdir, exist_ok=True)
        path = os.path.join(outdir, 'pw.nscf.in')
        write_pw_input(path, inp)
        click.echo(f'Wrote nscf input: {path}')
        # copy pseudos
        try:
            import shutil
            if pseudos_map:
                target = os.path.join(outdir, pseudo_dir)
                os.makedirs(target, exist_ok=True)
                for elm, src in pseudos_map.items():
                    if os.path.exists(src):
                        dst = os.path.join(target, os.path.basename(src))
                        shutil.copy(src, dst)
                        click.echo(f'Copied pseudo {src} -> {dst}')
        except Exception:
            pass
        if save_kmesh and kpoints:
            os.makedirs(outdir, exist_ok=True)
            kmesh_path = os.path.join(outdir, 'kmesh.txt')
            with open(kmesh_path, 'w') as f:
                f.write(kpoints.strip() + "\n")
            click.echo(f'Wrote k-mesh info: {kmesh_path}')


@cli.command()
@click.argument('cif', type=click.Path(exists=True))
@click.option('--workdir', default='.', help='Top-level workdir to place results under <workdir>/<seedname>/{scf,nscf,bands,projwfc,wan}')
@click.option('--num-wann', default=10)
@click.option('--num-bands', default=None, type=int)
@click.option('--bands-plot/--no-bands-plot', default=True)
@click.option('--dis-win-max', default=70.0, type=float)
@click.option('--dis-froz-max', default=30.0, type=float)
@click.option('--dis-num-iter', default=500, type=int)
@click.option('--num-iter', default=400, type=int)
@click.option('--dis-mix-ratio', default=1.0, type=float)
@click.option('--projections', default=None, help='Projections block, e.g. "Fe:s;p;d" or leave empty for random')
@click.option('--atoms-frac/--atoms-cart', default=True)
@click.option('--kpoint-path/--no-kpoint-path', default=True, help='Include high-symmetry path block generated by seekpath')
@click.option('--nscf-dir', type=click.Path(exists=False), default=None, help='Path to nscf output directory to read the kpoints_list')
@click.option('--kmesh', default=None, help="Explicit kmesh 'Nx Ny Nz' to use for mp_grid in the .win (overrides --nscf-dir)")
@click.option('--bands-file', default=None, help='Path to bands output file to reference from the .win')
@click.option('--seedname', default='qewan', help='Seedname (basename) for the .win and pw2wannier seedname')
@click.option('--auto-projections/--no-auto-projections', default=True, help='Enable automatic projection selection in .win')
@click.option('--include-projections/--no-include-projections', default=True, help='Include explicit projections block in the .win (ignored when auto-projections enabled)')
@click.option('--conventional-cell/--no-conventional-cell', default=False, help='Use CIF conventional cell as parsed by ASE (default: reduce by symmetry)')
@click.option('--nspin', default=1, type=int, help='Number of spin channels (1 or 2). When 2, generate per-spin .win files')
@click.option('--wannier-plot/--no-wannier-plot', default=False, help='Enable wannier_plot in .win')
@click.option('--write-xyz/--no-write-xyz', default=True, help='Enable write_xyz in .win')
@click.option('--write-hr/--no-write-hr', default=True, help='Enable write_hr in .win')
def run_wan(cif, workdir, num_wann, num_bands, bands_plot, dis_win_max, dis_froz_max, dis_num_iter, num_iter, dis_mix_ratio, projections, atoms_frac, kpoint_path, nscf_dir, kmesh, bands_file, seedname, auto_projections, include_projections, conventional_cell, nspin, wannier_plot, write_xyz, write_hr):
    atoms = read_cif(cif)
    if not conventional_cell:
        try:
            reduced, _ = reduce_atoms_by_symmetry(atoms, to_primitive=True)
            if reduced is not None:
                atoms = reduced
                click.echo(f'Using symmetry-reduced primitive cell with {len(atoms)} atoms')
        except Exception:
            pass
    use_kmesh = None
    if kmesh:
        use_kmesh = kmesh
    elif nscf_dir:
        kmesh_path = os.path.join(nscf_dir, 'kmesh.txt')
        if os.path.exists(kmesh_path):
            with open(kmesh_path, 'r') as f:
                use_kmesh = f.read().strip()
    # If nscf produced an explicit kpoints_list, prefer to read it
    kpoints_list = None
    if nscf_dir:
        kp_list_path = os.path.join(nscf_dir, 'kpoints_list.txt')
        if os.path.exists(kp_list_path):
            with open(kp_list_path, 'r') as f:
                kpoints_list = []
                for line in f.readlines():
                    if not line.strip():
                        continue
                    parts = line.split()
                    # accept lines with 3 or 4 columns: kx ky kz [weight]
                    coords = tuple(map(float, parts[:3]))
                    kpoints_list.append(coords)

    # compute canonical outdir: <workdir>/<seedname>/wan and create it
    outdir = os.path.join(workdir or '.', seedname or 'qewan', 'wan')
    os.makedirs(outdir, exist_ok=True)

    # When spin-polarized, create separate .win files per spin channel in subdirs
    if int(nspin) == 2:
        for ch in ('up', 'dn'):
            p = generate_wannier_win(
                atoms,
                outdir,
                num_wann=num_wann,
                num_bands=num_bands,
                bands_plot=bands_plot,
                dis_win_max=dis_win_max,
                dis_froz_max=dis_froz_max,
                dis_num_iter=dis_num_iter,
                num_iter=num_iter,
                mp_grid=use_kmesh,
                dis_mix_ratio=dis_mix_ratio,
                projections=projections,
                atoms_frac=atoms_frac,
                kpoint_path=kpoint_path,
                kpoints_list=kpoints_list,
                bands_file=bands_file,
                    seedname=f"{seedname}_{ch}",
                    auto_projections=auto_projections,
                    include_projections=(include_projections and (not auto_projections)),
                    dis_froz_proj=True,
                    dis_proj_max=0.95,
                        dis_proj_min=0.02,
                        wannier_plot=wannier_plot,
                        write_xyz=write_xyz,
                        write_hr=write_hr,
            )
            click.echo(f'Wrote wannier90 .win: {p}')
    else:
        path = generate_wannier_win(
        atoms,
        outdir,
        num_wann=num_wann,
        num_bands=num_bands,
        bands_plot=bands_plot,
        dis_win_max=dis_win_max,
        dis_froz_max=dis_froz_max,
        dis_num_iter=dis_num_iter,
        num_iter=num_iter,
        mp_grid=use_kmesh,
        dis_mix_ratio=dis_mix_ratio,
        projections=projections,
        atoms_frac=atoms_frac,
        kpoint_path=kpoint_path,
        kpoints_list=kpoints_list,
        bands_file=bands_file,
        seedname=seedname,
        auto_projections=auto_projections,
        include_projections=(include_projections and (not auto_projections)),
        dis_froz_proj=True,
        dis_proj_max=0.95,
        dis_proj_min=0.02,
        wannier_plot=wannier_plot,
        write_xyz=write_xyz,
        write_hr=write_hr,
    )
    


if __name__ == '__main__':
    cli()
 
