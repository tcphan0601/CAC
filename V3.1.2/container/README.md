# CAC Container Management (Singularity)

This directory contains instructions and definition files to build and run the Concurrent Atomistic-Continuum (CAC) application version 3.1.2 using containers.

## 1. Building the Images
```bash
cd V3.1.2
sudo singularity build cac.sif ./container/cac.def
```

## 2. Local Usage
Use --bind option to give singularity access to folders containing necessary data for your simulations
```bash
mpirun -n 8 singularity run --bind /data1 cac.sif -in input.in
```

## 3. HPC Usage (SLURM)
To run CAC on a cluster, use the Singularity image (.sif). This allows the container to leverage the host's MPI libraries for high-performance communication.
A sample slurm script is provided in container/slurm_script. You will need to change check which module to load for singularity or if it is already included (which singularity).
