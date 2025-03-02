#!/bin/bash
#SBATCH --account=def-bprotas           # RAC account for Bartek, to be used for priority jobs
#SBATCH --mail-user=zigicj@mcmaster.ca # email address
#SBATCH --mail-type=ALL                 # Send email regarding all events with job
#SBATCH --nodes=1		        # number of nodes (32 cpus per node)
#SBATCH --ntasks-per-node=32	        # number of MPI processes
#SBATCH --mem=0M		        # total memory required; default unit is megabytes
#SBATCH --time=00-00:15		        # time (DD-HH:MM) or (HH:MM:SS)
#SBATCH --output=output_DNS_T40_dt1e-3_X1.2Y1.2_lin_machepssinL.out		# output file
#SBATCH --job-name=DNST40L12            	# keep track of jobs

module load StdEnv/2020 fftw-mpi/3.3.8 netcdf-fortran-mpi/4.5.2
rm -rf /home/zigicj/scratch/2D_KS_Optimization-1/DNS_T40_dt1e-3_X1.2Y1.2_lin_machepssinL/
mkdir /home/zigicj/scratch/2D_KS_Optimization-1/DNS_T40_dt1e-3_X1.2Y1.2_lin_machepssinL/
make clean
make type=DNS
srun ./prog			        # run the job
seff $SLURM_JOBID		        # short summary of cpu and memory efficiency, printed to output (for more details use: sacct -j <jobid>)
