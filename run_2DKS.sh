#!/bin/bash
#SBATCH --account=def-bprotas                               # RAC account for Bartek, to be used for priority jobs
#SBATCH --mail-user=zigicj@mcmaster.ca                      # email address
#SBATCH --mail-type=ALL                                     # Send email regarding all events with job
#SBATCH --nodes=1		                                    # number of nodes (32 cpus per node)
#SBATCH --ntasks-per-node=32	                            # number of MPI processes
#SBATCH --mem=0M		                                    # total memory required; default unit is megabytes
#SBATCH --time=00-00:30		                                # time (DD-HH:MM) or (HH:MM:SS)
#SBATCH --output=output_DNS_T10_N128_dt1e-4_X3.20Y3.20_sineL.out	# output file
#SBATCH --job-name=DNST10X3.20Y3.20sineL            	        # keep track of jobs

module load StdEnv/2020 fftw-mpi/3.3.8 netcdf-fortran-mpi/4.5.2
rm -rf /home/zigicj/projects/def-bprotas/zigicj/2D_KS_Optimization/DNS_T10_N128_dt1e-4_X3.20Y3.20_sineL/
mkdir /home/zigicj/projects/def-bprotas/zigicj/2D_KS_Optimization/DNS_T10_N128_dt1e-4_X3.20Y3.20_sineL/
make clean
make type=DNS
srun ./prog			    # run the job
seff $SLURM_JOBID		# short summary of cpu and memory efficiency, printed to output (for more details use: sacct -j <jobid>)
cp -r . /home/zigicj/projects/def-bprotas/zigicj/2D_KS_Optimization/DNS_T10_N128_dt1e-4_X3.20Y3.20_sineL/