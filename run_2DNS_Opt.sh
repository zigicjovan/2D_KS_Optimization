#!/bin/bash
#SBATCH --account=rrg-bprotas           # RAC account for Bartek, to be used for priority jobs
#SBATCH --mail-user=mathap1@mcmaster.ca # email address
#SBATCH --mail-type=ALL                 # Send email regarding all events with job
#SBATCH --nodes=1		        # number of nodes (32 cpus per node)
#SBATCH --ntasks-per-node=32	        # number of MPI processes
#SBATCH --mem=8192M		        # total memory required; default unit is megabytes
#SBATCH --time=01-00:00		        # time (DD-HH:MM) or (HH:MM:SS)
#SBATCH --output=ICOpt-%J.out		# output file
#SBATCH --job-name=OptIC            	# keep track of jobs

module restore 2DNS_modules

srun ./prog			        # run the job
seff $SLURM_JOBID		        # short summary of cpu and memory efficiency, printed to output (for more details use: sacct -j <jobid>)
