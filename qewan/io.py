from ase.io import read
from ase import Atoms
from typing import Dict, Optional, Tuple, List
import os
import math
import re
try:
    from seekpath import get_path
except Exception:
    get_path = None


def read_cif(path: str) -> Atoms:
    """Read a CIF file and return an ASE Atoms object."""
    return read(path)


def detect_pseudo_type_from_upf(path: str, head_bytes: int = 65536) -> dict:
    """Parse UPF header (PP_HEADER) and infer pseudo type.

        Returns a dict with keys:
            - pseudo_type_attr: value of pseudo_type attribute if present (str or None)
            - attrs: dict of all attributes found on the PP_HEADER tag (keys lowercased)
            - inferred: one of 'paw', 'us', 'nc', or 'unknown'
            - snippet: short preview of the PP_HEADER content
    """
    result = {"pseudo_type_attr": None, "attrs": {}, "inferred": "unknown", "snippet": ""}
    try:
        with open(path, 'r', errors='ignore') as f:
            head = f.read(head_bytes)
    except Exception:
        return result

    # Find a PP_HEADER opening tag (self-closing or with closing tag). Capture attributes string.
    m = re.search(r'(?is)<pp_header\b([^>]*)/?>', head)
    attrs = {}
    if m:
        attr_str = m.group(1)
        # find attr="value" pairs
        for k, v in re.findall(r"(\w[\w-]*)\s*=\s*\"([^\"]*)\"", attr_str):
            attrs[k.lower()] = v.strip()
        result['attrs'] = attrs
        if 'pseudo_type' in attrs:
            result['pseudo_type_attr'] = attrs['pseudo_type']
        # snippet: first 300 chars of the attr string
        result['snippet'] = (attr_str.strip().replace('\n', '\\n'))[:300]

    # If there is a full <PP_HEADER>...</PP_HEADER> block, inspect its body for short tokens
    m_block = re.search(r'(?is)<pp_header\b[^>]*>(.*?)</pp_header>', head)
    if m_block:
        pp_body = m_block.group(1)
        # prefer showing a snippet of the body if attributes absent
        if not result.get('snippet'):
            result['snippet'] = pp_body.strip().replace('\n', '\\n')[:300]
        # inspect lines for clear tokens like 'US' or 'PAW' or 'Ultrasoft'
        lines = [ln.strip() for ln in pp_body.splitlines() if ln.strip()]
        for ln in lines[:10]:
            lnl = ln.lower()
            # line starting with US (common in many UPF headers)
            if re.match(r'(?i)^us\b', ln):
                result['inferred'] = 'uspp'
                if not result.get('pseudo_type_attr'):
                    result['pseudo_type_attr'] = 'USPP'
                return result
            if 'ultrasoft' in lnl or 'uspp' in lnl:
                result['inferred'] = 'uspp'
                if not result.get('pseudo_type_attr'):
                    result['pseudo_type_attr'] = 'USPP'
                return result
            if re.match(r'(?i)^paw\b', ln) or 'paw' in lnl:
                result['inferred'] = 'paw'
                return result
            if 'norm' in lnl or 'norm-conserv' in lnl or 'normconserv' in lnl or re.search(r'(?i)oncv', lnl):
                result['inferred'] = 'nc'
                return result

    # If pseudo_type attribute present, use it
    ptype = result.get('pseudo_type_attr')
    if ptype:
        pl = ptype.lower()
        if 'paw' in pl:
            result['inferred'] = 'paw'
            return result
        if 'uspp' in pl or 'ultrasoft' in pl:
            result['inferred'] = 'uspp'
            # prefer to expose a normalized pseudo_type_attr when not present
            if not result.get('pseudo_type_attr'):
                result['pseudo_type_attr'] = 'USPP'
            return result
        if 'norm' in pl or 'nc' in pl or 'oncv' in pl:
            result['inferred'] = 'nc'
            return result

    # Check explicit boolean attributes
    if attrs.get('is_paw', '').upper() == 'T':
        result['inferred'] = 'paw'
        return result
    if attrs.get('is_ultrasoft', '').upper() == 'T':
        result['inferred'] = 'uspp'
        if not result.get('pseudo_type_attr'):
            result['pseudo_type_attr'] = 'USPP'
        return result

    # Fallback: look for child tag <pseudo_type> inside the file
    m2 = re.search(r'(?i)<pseudo_type>([^<]+)</pseudo_type>', head)
    if m2:
        val = m2.group(1).strip().lower()
        if 'paw' in val:
            result['inferred'] = 'paw'
            result['pseudo_type_attr'] = m2.group(1).strip()
            return result
        if 'uspp' in val or 'ultrasoft' in val:
            result['inferred'] = 'uspp'
            result['pseudo_type_attr'] = m2.group(1).strip() or 'USPP'
            return result
        if 'norm' in val or 'nc' in val or 'oncv' in val:
            result['inferred'] = 'nc'
            result['pseudo_type_attr'] = m2.group(1).strip()
            return result

    # Last resort: keyword search in header
    low = head.lower()
    if 'paw' in low:
        result['inferred'] = 'paw'
        return result
    if 'ultrasoft' in low or 'uspp' in low:
        result['inferred'] = 'uspp'
        if not result.get('pseudo_type_attr'):
            result['pseudo_type_attr'] = 'USPP'
        return result
    if 'norm-conserv' in low or 'normconserv' in low or re.search(r'(?i)oncv', low):
        result['inferred'] = 'nc'
        return result

    return result


def _squeeze_blank_lines(s: str) -> str:
    """Replace 3+ consecutive newlines (or newlines with spaces) with two newlines."""
    return re.sub(r"\n\s*\n\s*\n+", "\n\n", s)


def reduce_atoms_by_symmetry(atoms: Atoms, symprec: float = 1e-5, to_primitive: bool = True):
    """Return a reduced/conventional cell using spglib symmetry analysis.

    - If `to_primitive` is True, attempt to find the primitive cell via
      `spglib.find_primitive`. Returns (reduced_atoms, mapping) where mapping
      is None for primitive result.
    - If `to_primitive` is False, returns a reduced list of symmetry-unique
      atoms (one representative per equivalence class) plus the
      `equivalent_atoms` mapping from the dataset.

    If `spglib` is not installed, raises ImportError with installation hint.
    """
    try:
        import spglib
    except Exception:
        raise ImportError("spglib is required for symmetry reduction. Install with `pip install spglib`.")

    lattice = atoms.get_cell()
    positions = atoms.get_scaled_positions()
    numbers = atoms.get_atomic_numbers()
    cell = (lattice, positions, numbers)

    if to_primitive:
        prim = spglib.find_primitive(cell, symprec=symprec)
        if prim is None:
            # cannot find primitive; return original
            return atoms.copy(), None
        prim_lattice, prim_positions, prim_numbers = prim
        prim_atoms = Atoms(numbers=prim_numbers, cell=prim_lattice, scaled_positions=prim_positions)
        return prim_atoms, None
    else:
        dataset = spglib.get_symmetry_dataset(cell, symprec=symprec)
        eq_atoms = dataset.get('equivalent_atoms', None)
        if eq_atoms is None:
            return atoms.copy(), None
        # pick one representative per equivalence class (first occurrence)
        reps = {}
        for idx, rep in enumerate(eq_atoms):
            reps.setdefault(rep, idx)
        unique_indices = [reps[k] for k in sorted(reps.keys())]
        red_numbers = [numbers[i] for i in unique_indices]
        red_positions = [positions[i] for i in unique_indices]
        red_atoms = Atoms(numbers=red_numbers, cell=lattice, scaled_positions=red_positions)
        return red_atoms, eq_atoms


def atoms_to_pw_input(atoms: Atoms, calculation: str = "scf", kpoints: Optional[str] = None, pseudos: Optional[Dict[str, str]] = None, pseudo_dir: str = './pseudos', kpath: bool = False, npoints: int = 40, nspin: int = 1, starting_magnetization: Optional[Dict[str, float]] = None, occupations: Optional[str] = None, smearing: Optional[str] = None, degauss: Optional[float] = None, ecutrho: Optional[float] = None, ecutwfc: Optional[float] = None) -> str:
    """Generate a minimal pw.x input string from ASE Atoms.

    - calculation: one of 'scf', 'bands', 'nscf'
    - kpoints: either a string 'Nx Ny Nz' or a full K_POINTS block starting with 'K_POINTS'
    - pseudos: dict mapping element -> pseudo filename
    """
    cell = atoms.get_cell()
    lattice_str = "\n".join(["  %16.10f %16.10f %16.10f" % tuple(cell[i]) for i in range(3)])

    # use fractional (scaled) coordinates for ATOMIC_POSITIONS
    atomic_positions = []
    try:
        scaled_positions = atoms.get_scaled_positions()
    except Exception:
        # fallback to converting cartesian to scaled
        from numpy import array, linalg
        C = array(atoms.get_cell()).T
        invC = linalg.inv(C)
        scaled_positions = (array(atoms.get_positions()) @ invC).tolist()
    for sym, pos in zip(atoms.get_chemical_symbols(), scaled_positions):
        atomic_positions.append("%s %16.10f %16.10f %16.10f" % (sym, pos[0], pos[1], pos[2]))

    # Preserve the order of first appearance in ATOMIC_POSITIONS
    symbols = atoms.get_chemical_symbols()
    species = []
    for s in symbols:
        if s not in species:
            species.append(s)
    if pseudos is None:
        pseudos = {s: f"{s}.UPF" for s in species}

    # Ensure we write only basenames in ATOMIC_SPECIES (not absolute paths)
    def _pseudo_basename(p):
        if not p or not isinstance(p, str):
            return p
        try:
            # if it looks like a path, use basename
            if os.path.sep in p or p.startswith('./') or p.startswith('..'):
                return os.path.basename(p)
        except Exception:
            pass
        return p

    pseudos_str = "\n".join([f"{s} {_pseudo_basename(pseudos.get(s, s + '.UPF'))}" for s in species])

    # If pseudos were provided as absolute paths and those files exist,
    # set `pseudo_dir` to the directory containing them so the generated
    # `pseudo_dir` in the pw input points to the real UPF location when
    # we are *not* copying pseudos into the outdir.
    try:
        abs_paths = [v for v in pseudos.values() if isinstance(v, str) and os.path.isabs(v)]
        if abs_paths and any(os.path.exists(p) for p in abs_paths):
            pseudo_dir = os.path.commonpath([os.path.dirname(p) for p in abs_paths])
    except Exception:
        pass

    # K-points selection logic.
    # For `calculation=='bands'` we require `seekpath` to be available and
    # generate a k-path by default. If seekpath is missing, raise an error
    # (seekpath is listed as a required dependency in requirements.txt).
    if calculation == 'bands' and get_path is None:
        raise RuntimeError("seekpath is required to generate band k-paths; install the dependency listed in requirements.txt")

    # Default for `calculation=='bands'`: generate a k-path via seekpath when
    # available. Respect an explicit `K_POINTS` block if provided by the user.
    if calculation == 'bands' and get_path is not None:
        if kpoints is not None:
            kp_strip = kpoints.strip()
            if kp_strip.upper().startswith('K_POINTS'):
                kpoints_block = kpoints if kpoints.endswith('\n') else kpoints + '\n'
            else:
                kpoints_block = generate_kpath(atoms, npoints=npoints)
        else:
            kpoints_block = generate_kpath(atoms, npoints=npoints)
    else:
        # allow passing a full K_POINTS block (starts with 'K_POINTS')
        if kpoints is None:
            kpoints = "6 6 6"
        kp_strip = kpoints.strip()
        if kp_strip.upper().startswith('K_POINTS'):
            # assume user passed a full block
            kpoints_block = kpoints if kpoints.endswith('\n') else kpoints + '\n'
        else:
            kpoints_block = f"K_POINTS automatic\n{kpoints} 0 0 0\n"

    # Add extra flags for nscf runs (disable symmetry/ inversion)
    extra_system_flags = ""
    if calculation == 'nscf':
        extra_system_flags = "    nosym = .true.\n    noinv = .true.\n"

    # spin polarization options
    spin_lines = [f"    nspin = {int(nspin)}\n"]
    if int(nspin) > 1:
        if starting_magnetization is None:
            starting_magnetization = {}
        # starting_magnetization(i) uses 1-based species indexing
        for i, s in enumerate(species, start=1):
            val = starting_magnetization.get(s, 0.0)
            spin_lines.append(f"    starting_magnetization({i}) = {float(val)}\n")
    spin_block = "".join(spin_lines)
    # smearing / occupations options
    smearing_lines = []
    if occupations is not None:
        smearing_lines.append(f"    occupations = '{occupations}'\n")
    if smearing is not None:
        smearing_lines.append(f"    smearing = '{smearing}'\n")
    if degauss is not None:
        try:
            dval = float(degauss)
            smearing_lines.append(f"    degauss = {dval}\n")
        except Exception:
            pass
    smearing_block = "".join(smearing_lines)
    ecutrho_block = ""
    # determine ecutwfc (use provided or default 40)
    try:
        ecutwfc_val = float(ecutwfc) if ecutwfc is not None else 40.0
    except Exception:
        ecutwfc_val = 40.0
    # if user provided ecutrho, use it; otherwise infer from pseudos and ecutwfc
    if ecutrho is not None:
        try:
            ecr = float(ecutrho)
            ecutrho_block = f"    ecutrho = {ecr}\n"
        except Exception:
            ecutrho_block = ""
    else:
        # infer multiplier: prefer 8 if any pseudo is PAW or US/ultrasoft, else 4
        try:
            # Use the dedicated UPF header parser to determine pseudo type
            has_uspp_or_paw = False
            has_nc = False
            if pseudos:
                for v in pseudos.values():
                    if not v or not isinstance(v, str):
                        continue
                    candidates = []
                    # absolute path
                    if os.path.isabs(v) and os.path.exists(v):
                        candidates.append(v)
                    # check in provided pseudo_dir
                    pd = pseudo_dir or './pseudos'
                    pjoin = os.path.join(pd, v)
                    if os.path.exists(pjoin):
                        candidates.append(pjoin)
                    # basename in pseudo_dir
                    b = os.path.basename(v)
                    pjoin2 = os.path.join(pd, b)
                    if os.path.exists(pjoin2):
                        candidates.append(pjoin2)
                    # relative or direct path
                    if os.path.exists(v):
                        candidates.append(v)

                    for cpath in candidates:
                        try:
                            res = detect_pseudo_type_from_upf(cpath)
                        except Exception:
                            res = {'inferred': 'unknown'}
                        inf = res.get('inferred', 'unknown')
                        if inf in ('paw', 'uspp'):
                            has_uspp_or_paw = True
                            break
                        if inf == 'nc':
                            has_nc = True
                    if has_uspp_or_paw:
                        break

            mult = 8 if has_uspp_or_paw else 4
            ecr = float(mult * ecutwfc_val)
            ecutrho_block = f"    ecutrho = {ecr}\n"
        except Exception:
            ecutrho_block = ""

    # Build &SYSTEM block in requested order: ecutwfc, ecutrho, occupations/smearing/degauss, extra flags, spin
    ecutwfc_line = f"    ecutwfc = {float(ecutwfc_val)}\n"
    system_block = ecutwfc_line + ecutrho_block + smearing_block + extra_system_flags + spin_block

    input_template = f"""
&CONTROL
    calculation = '{calculation}'
    prefix = 'qewan'
    outdir = './tmp'
    pseudo_dir = '{pseudo_dir}'
/
&SYSTEM
    ibrav = 0
    nat = {len(atoms)}
    ntyp = {len(species)}
{system_block}/ 
&ELECTRONS
    conv_thr = 1.0d-8
/

CELL_PARAMETERS angstrom
{lattice_str}

ATOMIC_SPECIES
{pseudos_str}

ATOMIC_POSITIONS crystal
{chr(10).join(atomic_positions)}

{kpoints_block}
"""
    return _squeeze_blank_lines(input_template)


def generate_kpath(atoms: Atoms, npoints: int = 40) -> str:
    """Generate a K_POINTS crystal_b block following the high-symmetry path using seekpath.

    Falls back to an automatic grid if seekpath unavailable or fails.
    """
    if get_path is None:
        return f"K_POINTS automatic\n6 6 6 0 0 0\n"

    cell = [list(vec) for vec in atoms.get_cell()]
    try:
        positions_frac = atoms.get_scaled_positions().tolist()
    except Exception:
        from numpy import array, linalg
        positions = atoms.get_positions()
        C = array(cell).T
        invC = linalg.inv(C)
        positions_frac = (array(positions) @ invC).tolist()
    numbers = atoms.get_atomic_numbers().tolist()
    structure = (cell, positions_frac, numbers)

    try:
        path_data = get_path(structure)
    except Exception:
        return f"K_POINTS automatic\n6 6 6 0 0 0\n"

    point_coords = path_data.get('point_coords', {})
    path = path_data.get('path', {})
    if isinstance(path, dict):
        ordered_segments = list(path.values())
    elif isinstance(path, list):
        ordered_segments = path
    else:
        ordered_segments = []

    klist = []
    total_segments = len(ordered_segments)
    if total_segments == 0:
        return f"K_POINTS automatic\n6 6 6 0 0 0\n"

    points_per_seg = max(2, int(math.ceil(npoints / total_segments)))
    for seg in ordered_segments:
        for i in range(len(seg) - 1):
            a = seg[i]
            b = seg[i + 1]
            ca = point_coords.get(a)
            cb = point_coords.get(b)
            if ca is None or cb is None:
                continue
            for t in range(points_per_seg):
                frac = t / float(points_per_seg)
                kx = ca[0] * (1 - frac) + cb[0] * frac
                ky = ca[1] * (1 - frac) + cb[1] * frac
                kz = ca[2] * (1 - frac) + cb[2] * frac
                klist.append((kx, ky, kz))
    last_seg = ordered_segments[-1]
    if last_seg:
        last_label = last_seg[-1]
        if last_label in point_coords:
            klist.append(tuple(point_coords[last_label]))

    uniq = []
    for k in klist:
        if not uniq or any(abs(k[i] - uniq[-1][i]) > 1e-6 for i in range(3)):
            uniq.append(k)

    lines = [f"K_POINTS crystal_b", str(len(uniq))]
    for k in uniq:
        lines.append(f"{k[0]:.8f} {k[1]:.8f} {k[2]:.8f} 1")
    return "\n".join(lines) + "\n"


def generate_mp_kpoints_list(atoms: Atoms, kmesh: str) -> Tuple[List[Tuple[float, float, float]], str]:
    """Generate an explicit Monkhorst-Pack k-point list (fractional coords) from a string 'Nx Ny Nz'.

    Returns (pts, block) where pts is list of tuples and block is PW K_POINTS block with normalized weights.
    """
    parts = kmesh.split()
    if len(parts) != 3:
        raise ValueError("kmesh must be 'Nx Ny Nz'")
    try:
        nx, ny, nz = [int(p) for p in parts]
    except Exception:
        raise ValueError("kmesh values must be integers")

    pts = []
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                kx = i / float(nx)
                ky = j / float(ny)
                kz = k / float(nz)
                pts.append((kx, ky, kz))

    total = len(pts)
    weight = 1.0 / float(total) if total > 0 else 1.0
    lines = ["K_POINTS crystal", str(len(pts))]
    for p in pts:
        lines.append(f"{p[0]:.8f} {p[1]:.8f} {p[2]:.8f} {weight:.8f}")
    block = "\n".join(lines) + "\n"
    return pts, block


def generate_projwfc_input(outdir: str, prefix: str = 'qewan', outdir_pw: str = './tmp', lsym: bool = False, lwrite_overlaps: bool = False, filproj: Optional[str] = None, fillowdin: str = 'lowdin.txt') -> str:
    """Write a simple projwfc.x input file with common options.

    Defaults: filproj -> '{prefix}.bands.dat.proj', fillowdin -> 'lowdin.txt'
    """
    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, 'projwfc.in')
    lsym_str = '.true.' if lsym else '.false.'
    lov_str = '.true.' if lwrite_overlaps else '.false.'
    filproj = filproj or f"{prefix}.bands.dat.proj"
    content = f"""&PROJWFC
    prefix = '{prefix}'
    outdir = '{outdir_pw}'
    lsym = {lsym_str}
    lwrite_overlaps = {lov_str}
    filproj = '{filproj}'
    fillowdin = '{fillowdin}'
/
"""
    content = _squeeze_blank_lines(content)
    with open(path, 'w') as f:
        f.write(content)
    return path


def generate_pw2wannier_input(
    outdir: str,
    prefix: str = 'qewan',
    nscf_outdir: str = './tmp',
    seedname: Optional[str] = None,
    write_mmn: bool = True,
    write_amn: bool = True,
    write_unk: bool = False,
    auto_projections: bool = True,
    spin_component: Optional[str] = None,
) -> str:
    """Write a minimal pw2wannier90 input file that points to the NSCF `outdir`.

    If a `pw.nscf.in` file exists inside the provided `nscf_outdir` directory,
    the `outdir` value from its `&CONTROL` block is used verbatim. Otherwise
    the absolute path to `nscf_outdir` is written.
    """
    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, 'pw2wannier.in')

    outdir_value = None
    try:
        if os.path.isdir(nscf_outdir):
            candidate = os.path.join(nscf_outdir, 'pw.nscf.in')
            if os.path.exists(candidate):
                with open(candidate, 'r') as f:
                    txt = f.read()
                m = re.search(r"outdir\s*=\s*['\"]([^'\"]+)['\"]", txt)
                if m:
                    outdir_value = m.group(1)
        if outdir_value is None:
            outdir_value = os.path.abspath(nscf_outdir)
    except Exception:
        outdir_value = os.path.abspath(nscf_outdir)

    mmn = '.true.' if write_mmn else '.false.'
    amn = '.true.' if write_amn else '.false.'
    unk = '.true.' if write_unk else '.false.'
    seed = seedname or prefix
    atom_proj_line = "  atom_proj = .true.\n" if auto_projections else ""
    spin_comp_line = (
        f"  spin_component = '{spin_component}'\n"
        if spin_component in ('up', 'down')
        else ""
    )

    content = f"""&INPUTPP
  outdir = '{outdir_value}'
  prefix = '{prefix}'
  seedname = '{seed}'
{atom_proj_line}{spin_comp_line}  write_mmn = {mmn}
  write_amn = {amn}
  write_unk = {unk}
/
"""

    content = _squeeze_blank_lines(content)
    with open(path, 'w') as f:
        f.write(content)
    return path


def write_pw_input(path: str, content: str):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as f:
        f.write(content)


def generate_wannier_win(
    atoms: Atoms,
    outdir: str,
    num_wann: int = 10,
    num_bands: Optional[int] = None,
    bands_plot: bool = True,
    dis_win_max: float = 70.0,
    dis_froz_max: float = 30.0,
    dis_num_iter: int = 500,
    num_iter: int = 400,
    mp_grid: Optional[str] = None,
    dis_mix_ratio: float = 1.0,
    projections: Optional[str] = None,
    include_projections: bool = True,
    atoms_frac: bool = True,
    kpoint_path: bool = False,
    kpoints_list: Optional[list] = None,
    bands_file: Optional[str] = None,
    seedname: str = 'qewan',
    auto_projections: bool = True,
    dis_froz_proj: bool = True,
    dis_proj_max: float = 0.95,
    dis_proj_min: float = 0.02,
    wannier_plot: bool = False,
    write_xyz: bool = True,
    write_hr: bool = True,
) -> str:
    """Create a Wannier90 .win using a richer template.

    Behavior:
    - If `kpoints_list` provided, write `mp_grid = N N N` (if mp_grid known) immediately above `begin kpoints` and write coords-only kpoints (no weights).
    - If `mp_grid` provided but no explicit list, generate coords-only list from mp_grid and write mp_grid above begin kpoints.
    - `bands_file` is optional hint written as `bands_file = <path>` in header.
    """
    num_bands_line = f"num_bands         =   {num_bands}" if num_bands is not None else ""
    bands_plot_line = "bands_plot        = true" if bands_plot else "bands_plot        = false"
    bands_file_line = f"bands_file        = {bands_file}" if bands_file else ""
    auto_proj_line = "auto_projections = .true." if auto_projections else "auto_projections = .false."
    dis_froz_proj_line = "dis_froz_proj = .true." if dis_froz_proj else "dis_froz_proj = .false."
    dis_proj_max_line = f"dis_proj_max =   {dis_proj_max:.2f}"
    dis_proj_min_line = f"dis_proj_min =   {dis_proj_min:.2f}"
    dis_win_max_line = f"dis_win_max       = {dis_win_max:.1f}d0"
    dis_froz_max_line = f"dis_froz_max      = {dis_froz_max:.1f}d0"
    dis_num_iter_line = f"dis_num_iter      = {dis_num_iter}"
    # default wannier output flags (can be overridden by args)
    wannier_plot_line = "wannier_plot        = .true." if wannier_plot else "wannier_plot        = .false."
    write_xyz_line = "write_xyz           = .true." if write_xyz else "write_xyz           = .false."
    write_hr_line = "write_hr            = .true." if write_hr else "write_hr            = .false."
    # select disentanglement block: projection-based or window-based
    if auto_projections:
        disent_block = "\n".join([auto_proj_line, dis_froz_proj_line, dis_proj_max_line, dis_proj_min_line])
    else:
        disent_block = "\n".join([dis_win_max_line, dis_froz_max_line, dis_num_iter_line])
    projections_block = "random" if not projections else projections

    kpoints_block = ""
    if kpoints_list is not None and len(kpoints_list) > 0:
        lines: List[str] = []
        if mp_grid:
            lines.append(f"mp_grid           = {mp_grid}")
        lines.append("begin kpoints")
        for p in kpoints_list:
            lines.append(f"  {p[0]:.8f}  {p[1]:.8f}  {p[2]:.8f}")
        lines.append("end kpoints")
        kpoints_block = "\n".join(lines) + "\n"
    elif mp_grid:
        try:
            pts, _block = generate_mp_kpoints_list(atoms, mp_grid)
        except Exception:
            pts = []
        if pts:
            lines = [f"mp_grid           = {mp_grid}", "begin kpoints"]
            for p in pts:
                lines.append(f"  {p[0]:.8f}  {p[1]:.8f}  {p[2]:.8f}")
            lines.append("end kpoints")
            kpoints_block = "\n".join(lines) + "\n"

    kpoint_path_block = ""
    if kpoint_path and get_path is not None:
        try:
            cell = [list(vec) for vec in atoms.get_cell()]
            positions_frac = atoms.get_scaled_positions().tolist()
            numbers = atoms.get_atomic_numbers().tolist()
            structure = (cell, positions_frac, numbers)
            pdata = get_path(structure)
            point_coords = pdata.get('point_coords', {})
            path = pdata.get('path', {})
            segments = []
            if isinstance(path, dict):
                for seg in path.values():
                    segments.extend(seg)
            elif isinstance(path, list):
                for seg in path:
                    segments.extend(seg)
            kp_lines = ["begin kpoint_path"]
            for i in range(len(segments) - 1):
                a = segments[i]
                b = segments[i + 1]
                if a == b:
                    continue
                ca = point_coords.get(a)
                cb = point_coords.get(b)
                if ca is None or cb is None:
                    continue
                # skip degenerate pairs (identical coordinates)
                if all(abs(ca[j] - cb[j]) < 1e-8 for j in range(3)):
                    continue
                kp_lines.append(f"{a} {ca[0]:.4f} {ca[1]:.4f} {ca[2]:.4f}\t {b} {cb[0]:.4f} {cb[1]:.4f} {cb[2]:.4f}")
            kp_lines.append("end kpoint_path")
            kpoint_path_block = "\n".join(kp_lines) + "\n"
        except Exception:
            kpoint_path_block = ""

    # write cell in angstroms for unit_cell_cart
    cell = atoms.get_cell()
    cell_lines = "\n".join([" ".join([f"{x:.5f}" for x in vec]) for vec in cell])

    if atoms_frac:
        atoms_lines = []
        fracs = atoms.get_scaled_positions()
        for s, p in zip(atoms.get_chemical_symbols(), fracs):
            atoms_lines.append(f"{s}  {p[0]:.3f}  {p[1]:.3f}  {p[2]:.3f}")
        atoms_block = "begin atoms_frac\n" + "\n".join(atoms_lines) + "\nend atoms_frac\n"
    else:
        atoms_lines = []
        carts = atoms.get_positions()
        for s, p in zip(atoms.get_chemical_symbols(), carts):
            atoms_lines.append(f"{s} {p[0]} {p[1]} {p[2]}")
        atoms_block = "begin atoms_cart\n" + "\n".join(atoms_lines) + "\nend atoms_cart\n"

    # only include projections block if requested and not using auto_projections
    # avoid backslashes inside f-string expressions by preparing projections section
    if include_projections:
        proj_section = 'begin projections\n' + projections_block + '\nend projections\n'
    else:
        proj_section = ''

    win = f"""
{num_bands_line}
num_wann          =   {num_wann}

{bands_plot_line}
{wannier_plot_line}
{write_xyz_line}
{write_hr_line}
{bands_file_line}
{disent_block}

num_iter          = {num_iter}

dis_mix_ratio = {dis_mix_ratio}

begin unit_cell_cart
angstrom
{cell_lines}
end unit_cell_cart

{atoms_block}

{proj_section}

{kpoint_path_block}

{kpoints_block}
"""

    outpath = os.path.join(outdir, f'{seedname}.win')
    os.makedirs(outdir, exist_ok=True)
    with open(outpath, 'w') as f:
        f.write(_squeeze_blank_lines(win))
    return outpath
