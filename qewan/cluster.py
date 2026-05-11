import os
import subprocess
from typing import Optional

def write_slurm_script(outdir: str, command: str, ntasks: int = 16, time: str = '04:00:00', partition: str = 'compute', job_name: str = 'qewan', modules: Optional[str] = None) -> str:
    os.makedirs(outdir, exist_ok=True)
    script_path = os.path.join(outdir, 'submit.slurm')
    modules_lines = ''
    if modules:
        modules_lines = f'module load {modules}\n'

    script = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --ntasks={ntasks}
#SBATCH --time={time}
#SBATCH --partition={partition}
#SBATCH --output=%x.%j.out

{modules_lines}
cd {outdir}

echo "Running: {command}"
{command}
"""
    with open(script_path, 'w') as f:
        f.write(script)
    os.chmod(script_path, 0o755)
    return script_path


def submit_slurm(script_path: str) -> subprocess.CompletedProcess:
    return subprocess.run(['sbatch', script_path], check=False, capture_output=True, text=True)
