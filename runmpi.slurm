#!/usr/bin/bash
#SBATCH --job-name=MPIshu
#SBATCH -p super  	             # partition requested: super, 256GB, 32GB, GPU
##SBATCH -N 2	   	               # nodes requested
#SBATCH -n 20    	               # ntasks requested; 1+num_proc in py file
##SBATCH --ntasks-per-node=5
##SBATCH --cpus-per-task=6
#SBATCH --mem=350G	             # memory in Mb 381952 for GPU
#SBATCH -t 01:00:00	             # time requested: 10-00:00:00
#SBATCH -o soutfile.txt  	       # send stdout to outfile
#SBATCH -e serrfile.txt  	       # send stderr to errfile

# Clear the environment from any previously loaded modules
#module purge > /dev/null 2>&1

#activate python3.7 env (ml37) that has mpi4py installed
source activate ml37

# run the mpi program
mpirun -np $SLURM_NTASKS python -m mpi4py.futures count_spacers_mpi.py

# submit job by running slurm script on head node after cd to workng dir
# sbatch runmpi.slurm
