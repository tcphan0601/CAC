# CAC Container Management (Singularity)

This directory contains instructions and definition files to build and run the Concurrent Atomistic-Continuum (CAC) application using containers.

## 1. Install singularity
This is for installing singularity on RHEL (or any RHEL-based distro). For other distributions, check their own instruction
```bash
sudo dnf install singularity-ce
```

## 1. Building the Images
You can download the image from here: https://drive.google.com/drive/folders/1QwPFWlkIBCKFQR0iyiMpP0tp7-PNjPDx?usp=drive_link
Or build your own image on local machine. The image file can be copied to and used on any machine and no need to be rebuild.

```bash
cd V3.2.0
sudo singularity build cac.sif ./container/cac.def
```

## 2. Local Usage
Use --bind option to give singularity access to folders containing necessary data for your simulations, such as data files, input script, potential files, output directories, if needed. By default, the current directory that it is called on is automatically binded. Once a folder is binded, all of it's sub directories are binded.

```bash
mpirun -n 8 singularity run --bind [list-of-paths-to-folders-containing-your-simulation-data] [path-to-your-cac.sif] -in input.in
mpirun -n 8 singularity run --bind /data1   /home/tphan4/cac.sif -in input.in
```

## 3. HPC Usage (SLURM)
To run CAC on a cluster, use the Singularity image (.sif). This allows the container to leverage the host's MPI libraries for high-performance communication.
A sample slurm script is provided in container/slurm_script. You will need to change check which module to load for singularity (module load singularitypro for Expanse)  or if it is already included (which singularity on Anvil).
