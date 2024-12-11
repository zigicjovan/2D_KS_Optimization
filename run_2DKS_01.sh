#!/bin/bash
#SBATCH --account=def-bprotas           # RAC account for Bartek, to be used for priority jobs
#SBATCH --mail-user=zigicj@mcmaster.ca # email address
#SBATCH --mail-type=ALL                 # Send email regarding all events with job
#SBATCH --nodes=2		        # number of nodes (32 cpus per node)
#SBATCH --ntasks-per-node=32	        # number of MPI processes
#SBATCH --mem=0M		        # total memory required; default unit is megabytes
#SBATCH --time=00-02:00		        # time (DD-HH:MM) or (HH:MM:SS)
#SBATCH --output=output_DNS_T10_dt1e-5_X1.01Y1.01_lin_macheps.out		# output file
#SBATCH --job-name=DNS01            	# keep track of jobs

module load StdEnv/2020 fftw-mpi/3.3.8 netcdf-fortran-mpi/4.5.2
rm -rf DNS_T10_dt1e-5_X1.01Y1.01_lin_macheps
mkdir DNS_T10_dt1e-5_X1.01Y1.01_lin_macheps
make clean
make type=DNS
srun ./prog			        # run the job
seff $SLURM_JOBID		        # short summary of cpu and memory efficiency, printed to output (for more details use: sacct -j <jobid>)
