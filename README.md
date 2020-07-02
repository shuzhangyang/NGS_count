# NGS_count
Count sgRNAs in fastq files

The MPI version (count_spacers_mpi.py, a distributed version) is acompanied by a slurm script (runmpi.slurm, run in bash by sbatch), which scale up the job across a cluster of multiple nodes. It is suitable for large number of fastq files.
