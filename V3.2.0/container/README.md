# CAC Container Management (Singularity)

This directory contains instructions and definition files to build and run the Concurrent Atomistic-Continuum (CAC) application using containers.

## 1. Install singularity
This is for installing singularity on RHEL (or any RHEL-based distro). For other distributions, check their own instruction
```bash
sudo dnf install singularity-ce
```

## 1. Building the Images
```bash
cd V3.2.0
sudo singularity build cac.sif ./container/cac.def
```

## 2. Local Usage
Use --bind option to give singularity access to folders containing necessary data for your simulations

```bash
mpirun -n 8 singularity run --bind [list-of-paths-to-folders-containing-your-simulation-data] [path-to-your-cac.sif] -in input.in
mpirun -n 8 singularity run --bind /data1   /home/tphan4/cac.sif -in input.in
```

## 3. HPC Usage (SLURM)
To run CAC on a cluster, use the Singularity image (.sif). This allows the container to leverage the host's MPI libraries for high-performance communication.
A sample slurm script is provided in container/slurm_script. You will need to change check which module to load for singularity (module load singularitypro for Expanse)  or if it is already included (which singularity on Anvil).
