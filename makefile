# Makefile to easily update the compilation of a program (prog)
# --------
#
# by Alain Veilleux, 4 August 1993
#    Modified to run 2D Navier-Stokes DNS by Pritpal Matharu: 2020/07/17
#    Modified to run 2D Kuramoto-Sivashinsky DNS by Jovan Zigic: 2024/11/22
#
# GOAL AND FUNCTIONING OF THIS SCRIPT:
#    Script in the form of a "Makefile" allowing to update a program containing
#    multiple separated routines on the disk. This script is not executed by itself,
#    but is instead read and interpreted by the "make" command. When it is called,
#    the "make" command verifies the dates of the various components your program is
#    built from. Only routines that were modified after the last compilation of the
#    program are recompiled in object form (files ending in .o). Recompiled .o files
#    are subsequently linked together to form an updated version of the final program.
#
# TO ADAPT THIS SCRIPT TO YOUR PROGRAM:
#    Modify the contents of the variables hereunder. Comments will guide you how and
#    where.
#
# USING "make" ON THE UNIX COMMAND LINE:
#    1- Type "make" to update the whole program.
#    2- Type "make RoutineName" to only update the RoutineName routine.
#    3- Type "make type=sim_type" to define the object files to create for the simulation.
#       sim_type can be Opt, DNS, or Kappa (by default Opt, for optimization scheme).
#
#================================================================================================
#
# Before running make command, must load the modules!!!
# Load Modules : StdEnv/2020 fftw-mpi/3.3.8 netcdf-fortran-mpi/4.5.2
# Modules were loaded using:    module load StdEnv/2020 fftw-mpi/3.3.8 netcdf-fortran-mpi/4.5.2
# Then saved as "2DNS_modules": module save 2DKS_modules
# We now reload these modules : module restore 2DKS_modules
#
# To send the job (prog) to the scheduler, appropriately modify the shell file "run_2DNS_****.sh"
# then run: sbatch run_2DKS.sh
#
#================================================================================================
#
# The instructions outlined in this makefile are equivalent to manually running the following commands:
# module load fftw-mpi/3.3.8
# module load netcdf-fortran-mpi/4.5.2
# export FFTWDIR=/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/MPI/intel2020/openmpi4/fftw-mpi/3.3.8/
# export nCDF_DIR=/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/MPI/intel2020/openmpi4/netcdf-fortran-mpi/4.5.2
# export MKLROOT=/cvmfs/soft.computecanada.ca/easybuild/software/2017/Core/imkl/11.3.4.258/mkl
# mpif90 -lmpi -o prog global_variables.f90 nse_initialize.f90 data_ops.f90 fftwfunction.f90 function_ops.f90 solvers.f90 dnsNS2D_main.f90 -lm -I${FFTWDIR}/include -L${FFTWDIR}/lib -lfftw3_mpi -lfftw3 -I${nCDF_DIR}/include -L${nCDF_DIR}/lib -lnetcdff -lnetcdf
#
#================================================================================================
# Instructions to run program (Jovan Zigic, 2024/11/18):
# On startup:
# 	check work_pathname and scratch_pathname in global_variables.f90
# 	mkdir bin_files (if necessary)
# 	module load StdEnv/2020 fftw-mpi/3.3.8 netcdf-fortran-mpi/4.5.2
# For quick tests:
#   mkdir DNS_T10_dt1e-3_X1.0Y1.0
# 	make clean
# 	make type=DNS
# 	sbatch --account=def-bprotas --nodes=2 --ntasks-per-node=32 --mem=0M --time=0-01:00  ./prog > output_DNS_T400_dt1e-5_X1.0Y1.0_lin_macheps &
# To check job status:
# 	seff [job id] 
# 	sq
# Tp copy all files to a testing subdirectory:
# 	mkdir test1
#   cd test1
# 	scp /home/zigicj/projects/def-bprotas/zigicj/2D_KS_Optimization/* .
# To cancel:
# 	scancel [job id] 
# To add time:
# 	scontrol update jobid=[job_id] TimeLimit=0-02:00
# To remove a subdirectory:
# 	rm -rf [subdirectory] 
# To debug [not on VScode, only on terminal]:
# 	[install XQuartz] 
# 	ssh -Y zigicj@graham.alliancecan.ca
# 	[enter CCDB password and 2-factor] 
# 	module load StdEnv/2020 fftw-mpi/3.3.8 netcdf-fortran-mpi/4.5.2 ddt-cpu
# 	[turn on -g to use debugger (see below)]
# 	make clean
# 	make type=DNS
# 	salloc --x11 --time=0-3:00 --mem-per-cpu=4G --ntasks=4 -A def-bprotas
# 	ddt ./prog
#================================================================================================


#====================  Definition of variables  =====================
# Remark : variables are sometimes called "macros" in Makefiles.

# Compiler to use (MPI FORTRAN)
CompilerName= mpif90

# Libraries
BASIC_LIB  = -lm -lmpi                # Linking flags
FFTW_DIR   =                          # Include FFTW3 Library
FFTW_LIB   = -lfftw3_mpi -lfftw3 -lm  # Load FFTW3 Library
nCDF_DIR   =                          # Include nCDF Library
nCDF_LIB   = -lnetcdff -lnetcdf       # Load nCDF Library
#FFTW_DIR   = -I/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/MPI/intel2020/openmpi4/fftw-mpi/3.3.8/include/                            # Include FFTW3 Library
#FFTW_LIB   = -L/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/MPI/intel2020/openmpi4/fftw-mpi/3.3.8/lib -lfftw3_mpi -lfftw3 -lm         # Load FFTW3 Library
#nCDF_DIR   = -I/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/MPI/intel2020/openmpi4/netcdf-fortran-mpi/4.5.2/include                   # Include nCDF Library
#nCDF_LIB   = -L/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/MPI/intel2020/openmpi4/netcdf-fortran-mpi/4.5.2/lib -lnetcdff -lnetcdf    # Load nCDF Library
#FFTW_DIR   = -I/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx512/MPI/intel2020/openmpi4/fftw-mpi/3.3.8/include/                            # Include FFTW3 Library
#FFTW_LIB   = -L/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx512/MPI/intel2020/openmpi4/fftw-mpi/3.3.8/lib -lfftw3_mpi -lfftw3 -lm         # Load FFTW3 Library
#nCDF_DIR   = -I/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx512/MPI/intel2020/openmpi4/netcdf-fortran-mpi/4.5.2/include                   # Include nCDF Library
#nCDF_LIB   = -L/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx512/MPI/intel2020/openmpi4/netcdf-fortran-mpi/4.5.2/lib -lnetcdff -lnetcdf    # Load nCDF Library
##MKLROOT    = -I/cvmfs/soft.computecanada.ca/easybuild/software/2019/Core/imkl/2020.1.217/mkl
#MKLFFTW = -L/home/mathap1/mklfftw3/ -lfftw3x_cdft_lp64

# Compilation options: the below options are usually used to compile FORTRAN
#                      code. You can assign other values than those suggested
#                      in the "CompilationOptions" variables.
#CompilationOptions= -O3
# Remove the below "#" to activate compilation in debug mode
#CompilationOptions= -g
# Remove the below "#" to use "gprof", which indicates the computation time in
#    each subroutine
#CompilationOptions= -O3 -pg

# List of routines to compile: here we list all object files that are needed.
# Put a "\" at the end of each line that if you want to continue the list of
#    routines on the following line.

# Different routines, depending of type of simulation
# DNS
#ObjDNS= global_variables.o fftwfunction.o nse_initialize.o data_ops.o function_ops.o solvers.o dnsNS2D_main.o
ObjDNS= global_variables.o fftwfunction.o KS_initialize.o data_ops.o function_ops.o solvers.o dns2DKS_main.o
# Kappa Test
ObjKappa= global_variables.o fftwfunction.o nse_initialize.o data_ops.o function_ops.o solvers.o optimization.o KappaNS2D_main.o
# Optimization
ObjOpt= global_variables.o fftwfunction.o nse_initialize.o data_ops.o function_ops.o solvers.o optimization.o OptNS2D_main.o
#ObjectFiles= global_variables.o fftwfunction.o nse_initialize.o data_ops.o function_ops.o solvers.o optimization.o OptNS2D_main.o

# Determine if override at commandline was given
type= DNS 			# Default is to make routine for optimization scheme
ifeq ($(type), DNS)		# Make routine for DNS
        ObjectFiles= $(ObjDNS)
else ifeq ($(type), Kappa)	# Make routine for Kappa Test
        ObjectFiles= $(ObjKappa)
else				# Make routine for optimzation
    	ObjectFiles= $(ObjOpt)
endif

# Name of the final executable
ProgramOut= prog
#=====  End of variable definitions =====

# To send the job (prog) to the scheduler, appropriately modify the shell file "run_2DNS_****.sh"
# then run: sbatch run_2DNS_****.sh

#===============  There is nothing to change below this line  =============


# Defines a rule: how to build an object file (ending in ".o")
#                 from a source file (ending in ".f90")
# note: "$<" symbols will be replaced by the name of the file that is compiled
# Compiling Fortran files:
%.o: %.f90
	$(CompilerName) $(CompilationOptions) -c $(nCDF_DIR) $(FFTW_DIR) $<

# Dependencies of the main executable on the object files (".o") it is built from.
# The dependency of object files on source files (".f") is implied by the above
# implicit rules.
$(ProgramOut): $(ObjectFiles)
	$(CompilerName) $(CompilationOptions) -o $(ProgramOut) $(ObjectFiles) $(BASIC_LIB) $(FFTW_LIB) $(nCDF_LIB)
#       $(CompilerName) $(CompilationOptions) -o $(ProgramOut) $(ObjectFiles) $(BASIC_LIB) $(FFTW_LIB) $(nCDF_LIB) $(MKLFFTW)

clean:
	rm -f *.o *~ *.mod
