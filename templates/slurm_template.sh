
#!/bin/bash
# Generic SLURM template (fill variables or use qewan.cluster.write_slurm_script)
#SBATCH --job-name={{JOB_NAME}}
#SBATCH --ntasks={{NTASKS}}
#SBATCH --time={{TIME}}
#SBATCH --partition={{PARTITION}}
#SBATCH --output=%x.%j.out

module load espresso

# Run PWscf
mpirun -np $SLURM_NTASKS pw.x -inp pw.scf.in > pw.scf.out

# Run nscf/bands if present (example)
# mpirun -np $SLURM_NTASKS pw.x -inp pw.nscf.in > pw.nscf.out

# Run wannier90 steps (example)
# wannier90.x qewan.win > qewan.wout
