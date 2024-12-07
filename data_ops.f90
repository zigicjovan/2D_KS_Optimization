!===================================================================================================
! MODULE CONTAINS ROUTINES REPRESENTING OPERATIONS FOR WRITING OR READING FILES
!
! Author: Jovan Zigic (inherited from Pritpal Matharu)                                          
! McMaster University                                 
! Date: 2024/12/06                                  
!
! CONTAINS:
! (*) save_bin               - Saves bin file of vorticity field
! (*) read_bin               - Reads bin file of vorticity field
! (*) save_NS_Opt            - Saves optimal IC and corresponding diagnostic quantities
! (*) save_NS_DNS            - Saves vorticity field and diagnostic quantities from DNS
! (*) read_IC                - Reads initial condition from MATLAB
! (*) read_NS_Opt            - Reads initial condition from previous optimization run
! (*) read_BS_Opt            - Reads initial condition from for bootstrapping
! (*) save_NS_vorticity      - Saves the 2D vorticity field
! (*) save_NS_fwd            - Saves vorticity field and diagnostic (Old verison of save_NS_DNS)
! (*) save_field_R2toR1_ncdf - Saves 2D field uses netcdf (Old verison for save_NS_fwd)
! (*) save_field_R1toR1_ncdf - Saves 1D field uses netcdf (Old verison for save_NS_fwd)
! (*) ncdf_error_handle      - Prints error code from netcdf
!===================================================================================================
MODULE data_ops
  IMPLICIT NONE ! Prevent using implicit typing rules for entire module

  CONTAINS
    !==================================================================
    ! *** Save bin file ***
    ! Input: mydata  - field in physical space
    !==================================================================
    SUBROUTINE save_bin(mydata)
      ! Load variables
      USE global_variables, ONLY: pr, rank, n_nse, RESOL, nx_dim, visc, Lx, Ly, local_Ny, scratch_pathname, IC_type, normconstr, Grad_type, endTime, ell, Time_iter, Nchar, lchar, tchar, viscchar, Statinfo
      USE mpi                 ! Use MPI module (binding works well with fftw libraries)
      ! Initialize variables
      COMPLEX(pr), DIMENSION(:,:), INTENT(IN)  :: mydata       ! Vorticity field, to be saved
      INTEGER                                  :: file_handle  ! File handle for writing the file
      INTEGER(KIND=MPI_OFFSET_KIND)            :: my_offset    ! Offset for process to write in correct location
      CHARACTER(6)                             :: indexchar    ! Time iteration as character
      CHARACTER(2)                             :: timechar     ! Final time as character
      CHARACTER(200)                           :: filename     ! Filename for writing the file

      ! Iteration number as a character
      WRITE(indexchar, '(i6.6)') Time_iter
      WRITE(timechar, '(i2.2)')  int(endTime*1.0e1)
      ! Filename path for saving in the scratch folder, for current timestep
      filename = TRIM(scratch_pathname)//"Vort_fwdNS_"//IC_type//"_Norm"//normconstr//"_Grad"//Grad_type//"_N"//Nchar//"_NU"//TRIM(ADJUSTL(viscchar))//"_L"//TRIM(ADJUSTL(lchar))//"_T"//TRIM(ADJUSTL(tchar))//"_"//indexchar//".bin"
      ! Delete any existing file
      IF (rank==0) THEN
        CALL MPI_File_delete(filename, MPI_INFO_NULL, Statinfo)
      END IF
      ! Open the file for writing
      CALL MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, file_handle, Statinfo)
      ! Calculate the processors offset in the file
!      my_offset = rank*ny_dim*pr
      my_offset = 2*rank*nx_dim*pr
      ! Move the pointer for writing to the correct place for writing
      CALL MPI_File_seek(file_handle, my_offset, MPI_SEEK_SET, Statinfo)
      ! Write using the individual file pointer
!      CALL MPI_File_write(file_handle, mydata, ny_dim, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, Statinfo)
      CALL MPI_File_write(file_handle, mydata, 2*nx_dim, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, Statinfo)
      ! Close file
      CALL MPI_File_close(file_handle, Statinfo)
    END SUBROUTINE save_bin

    !==================================================================
    ! *** Read bin file ***
    ! Output: mydata  - field in physical space
    !==================================================================
    SUBROUTINE read_bin(mydata)
      ! Load variables
      USE global_variables, ONLY: pr, rank, n_nse, RESOL, nx_dim, visc, Lx, Ly, local_Ny, scratch_pathname, IC_type, normconstr, Grad_type, endTime, ell, Time_iter, Nchar, lchar, tchar, viscchar, Statinfo
      USE mpi                 ! Use MPI module (binding works well with fftw libraries)
      ! Initialize variables
      COMPLEX(pr), DIMENSION(:,:), INTENT(OUT) :: mydata       ! Vorticity field, to be read
      INTEGER                                  :: file_handle  ! File handle for writing the file
      INTEGER(KIND=MPI_OFFSET_KIND)            :: my_offset    ! Offset for process to write in correct location
      CHARACTER(6)                             :: indexchar    ! Time iteration as character
      CHARACTER(2)                             :: timechar     ! Final time as character
      CHARACTER(200)                           :: filename     ! Filename for writing the file

      ! Iteration number as a character
      WRITE(indexchar, '(i6.6)') Time_iter
      WRITE(timechar, '(i2.2)')  int(endTime*1.0e1)
      ! Filename path for saving in the scratch folder, for current timestep
      filename = TRIM(scratch_pathname)//"Vort_fwdNS_"//IC_type//"_Norm"//normconstr//"_Grad"//Grad_type//"_N"//Nchar//"_NU"//TRIM(ADJUSTL(viscchar))//"_L"//TRIM(ADJUSTL(lchar))//"_T"//TRIM(ADJUSTL(tchar))//"_"//indexchar//".bin"
      ! Open the file for reading
      CALL MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, file_handle, Statinfo)
      ! Calculate the processors offset in the file
!      my_offset = rank*ny_dim*pr
      my_offset = 2*rank*nx_dim*pr
      ! Move the pointer for reading to the correct place for writing
      CALL MPI_File_seek(file_handle, my_offset, MPI_SEEK_SET, Statinfo)
      ! Read using the individual file pointer
      CALL MPI_File_read(file_handle, mydata, 2*nx_dim, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, Statinfo)
      ! Close file
      CALL MPI_File_close(file_handle, Statinfo)
      ! Delete file after reading (to ensure that files don't accumulate)
      IF (rank==0) THEN
        CALL MPI_File_delete(filename, MPI_INFO_NULL, Statinfo)
      END IF
    END SUBROUTINE read_bin

    !==================================================================
    ! *** Save values from optimization scheme ***
    ! Input:   w0    - initial guess vorticity field in physical space
    !           w    - optimal vorticity field in physical space
    !         Pal    - palinstrophy vector
    !         Ens    - enstrophy vector
    !    KinEnerg    - kinetic energy vector
    !           t    - time vector
    !           J    - cost functional vector
    !==================================================================
    SUBROUTINE save_NS_Opt(w0, w, Pal, Ens, KinEnerg, t, J)
      ! Load variables
      USE global_variables    ! Declares variables and defines parameters of the system
      USE netcdf              ! Use netcdf for saving files
      USE mpi                 ! Use MPI module (binding works well with fftw libraries)
      ! Initialize variables
!      REAL(pr), DIMENSION(1:n_nse(1),1:local_Ny), INTENT(IN) :: w         ! Vorticity field, to be saved
      REAL(pr), DIMENSION(:,:), INTENT(IN) :: w0        ! Starting vorticity field, to be saved
      REAL(pr), DIMENSION(:,:), INTENT(IN) :: w         ! Vorticity field, to be saved
      REAL(pr), DIMENSION(:),   INTENT(IN) :: Pal       ! Palinstrophy vector, to be saved
      REAL(pr), DIMENSION(:),   INTENT(IN) :: Ens       ! Enstrophy vector, to be saved
      REAL(pr), DIMENSION(:),   INTENT(IN) :: KinEnerg  ! Kinetic energy vector, to be saved
      REAL(pr), DIMENSION(:),   INTENT(IN) :: t         ! Time vector, to be saved
      REAL(pr), DIMENSION(:),   INTENT(IN) :: J         ! Cost functional, to be saved
      CHARACTER(200)                       :: filename                                                       ! Filename for writing the file
      INTEGER, DIMENSION(1:2)              :: starts, counts, dimids                                         ! Positioning for saving
      INTEGER                              :: ncout, ncid                                                    ! Temporary integers for saving
      INTEGER                              :: varid, varid0, varidP, varidE, varidK, varidt, varidJ, varidit ! Temporary integers for saving
      INTEGER                              :: x_dimid, y_dimid, t_dimid, j_dimid, p_dimid                    ! Temporary dimensions for saving
      INTEGER                              :: ii, t_len, j_len                                               ! Temporary integer for looping and reference lengths

      ! Filename path for saving file
      filename = TRIM(work_pathname)//"Optimization_"//IC_type//"_Norm"//normconstr//"_Grad"//Grad_type//"_N"//Nchar//"_NU"//TRIM(ADJUSTL(viscchar))//"_L"//TRIM(ADJUSTL(lchar))//"_T"//TRIM(ADJUSTL(tchar))//".nc"
      CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)

      IF (rank==0) THEN
        t_len = size(t)
        j_len = size(J)
        ncout = nf90_create(filename, NF90_CLOBBER, ncid)
        IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
        PRINT *, " Created."
      END IF
      CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)

      ! Initialize file for writing
      IF (rank==0) THEN
        PRINT *, " Before."
        ncout = nf90_def_dim(ncid, "x", n_nse(1), x_dimid)
        PRINT *, " After."
        IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
        PRINT *, " x defined."
        ncout = nf90_def_dim(ncid, "y", NF90_UNLIMITED, y_dimid)
        IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
        PRINT *, " y defined."
        ncout = nf90_def_dim(ncid, 't', t_len, t_dimid)
        IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
        PRINT *, " t defined."
        ncout = nf90_def_dim(ncid, 'p', j_len, j_dimid)
        IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
        PRINT *, " p defined."
!        ncout = nf90_def_dim(ncid, 'iter', 1, p_dimid)
        dimids = (/ x_dimid, y_dimid /)
        PRINT *, " dimids defined."

        ncout = nf90_def_var(ncid, TRIM("w"), NF90_DOUBLE, dimids, varid)
        IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
        PRINT *, " w defined."
        ncout = nf90_def_var(ncid, TRIM("w0"), NF90_DOUBLE, dimids, varid0)
        IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
        PRINT *, " w0 defined."

        ncout = nf90_def_var(ncid, TRIM("Palin"), NF90_DOUBLE, (/ t_dimid /), varidP)
        IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
        PRINT *, " Palin defined."
        ncout = nf90_def_var(ncid, TRIM("Enst"), NF90_DOUBLE, (/ t_dimid /), varidE)
        IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
        PRINT *, " Enst defined."
        ncout = nf90_def_var(ncid, TRIM("Kin"), NF90_DOUBLE, (/ t_dimid /), varidK)
        IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
        PRINT *, " Kin defined."

        ncout = nf90_def_var(ncid, TRIM("tvec"), NF90_DOUBLE, (/ t_dimid /), varidt)
        IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
        PRINT *, " tvec defined."
        ncout = nf90_def_var(ncid, TRIM("Jvec"), NF90_DOUBLE, (/ j_dimid /), varidJ)
        IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
        PRINT *, " Jvec defined."
!        ncout = nf90_def_var(ncid, TRIM("p"), NF90_DOUBLE, (/ p_dimid /), varidp)
        ncout = nf90_enddef(ncid)
        PRINT *, " end of define."

        IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
        ncout = nf90_put_var(ncid, varidP, Pal)
        IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
        PRINT *, " put Palin."
        ncout = nf90_put_var(ncid, varidE, Ens)
        IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
        PRINT *, " put Ens."
        ncout = nf90_put_var(ncid, varidK, KinEnerg)
        IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
        PRINT *, " put KinEnerg."
        ncout = nf90_put_var(ncid, varidt, t)
        IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
        PRINT *, " put t."
        ncout = nf90_put_var(ncid, varidJ, J)
        IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
        PRINT *, " put J."
!        ncout = nf90_put_var(ncid, varidp, iter)

        ncout = nf90_close(ncid)
        IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
        PRINT *, " Closed."
      END IF
      CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)

      ! Positioning
      starts = (/ 1, rank*local_Ny+1 /)
      counts = (/ n_nse(1), local_Ny /)

      !!--------------------------
      !! START netCDF ROUTINES
      !!--------------------------
      ! Optimal field
      DO ii=0,np-1
        IF (rank==ii) THEN
          ncout = nf90_open(filename, NF90_WRITE, ncid)
          IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
          ncout = nf90_inq_varid(ncid, TRIM("w"), varid)
          IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
          ncout = nf90_put_var(ncid, varid, w, start = starts, count = counts)
          IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
          ncout = nf90_close(ncid)
          IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
        END IF
        CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)
      END DO
      ! Initial guess
      DO ii=0,np-1
        IF (rank==ii) THEN
          ncout = nf90_open(filename, NF90_WRITE, ncid)
          IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
          ncout = nf90_inq_varid(ncid, TRIM("w0"), varid0)
          IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
          ncout = nf90_put_var(ncid, varid0, w0, start = starts, count = counts)
          IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
          ncout = nf90_close(ncid)
          IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
        END IF
        CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)
      END DO
      CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)

    END SUBROUTINE save_NS_Opt

    !==================================================================
    ! *** Save values from DNS run in one file ***
    ! Input:   w0    - initial condition vorticity field in physical space
    !         Pal    - palinstrophy vector
    !         Ens    - enstrophy vector
    !    KinEnerg    - kinetic energy vector
    !           t    - time vector
    !        ipL2    - L^2 inner product vector
    !        ipH1    - H^1 inner product vector
    !        ipH2    - H^2 inner product vector
    !        ipHn1   - H^(-1) inner product vector
    !==================================================================
    SUBROUTINE save_NS_DNS(w0, Pal, Ens, KinEnerg, t, ipL2, ipH1, ipH2, ipHn1)
      ! Load variables
      USE global_variables    ! Declares variables and defines parameters of the system
      USE netcdf              ! Use netcdf for saving files
      USE mpi                 ! Use MPI module (binding works well with fftw libraries)
      ! Initialize variables
!      REAL(pr), DIMENSION(1:n_nse(1),1:local_Ny), INTENT(IN) :: w         ! Vorticity field, to be saved
      REAL(pr), DIMENSION(:,:), INTENT(IN) :: w0        ! Starting vorticity field, to be saved
      REAL(pr), DIMENSION(:),   INTENT(IN) :: Pal       ! Palinstrophy vector, to be saved
      REAL(pr), DIMENSION(:),   INTENT(IN) :: Ens       ! Enstrophy vector, to be saved
      REAL(pr), DIMENSION(:),   INTENT(IN) :: KinEnerg  ! Kinetic energy vector, to be saved
      REAL(pr), DIMENSION(:),   INTENT(IN) :: t         ! Time vector, to be saved
      REAL(pr), DIMENSION(:),   INTENT(IN) :: ipL2      ! inner product vector, to be saved
      REAL(pr), DIMENSION(:),   INTENT(IN) :: ipH1      ! inner product vector, to be saved
      REAL(pr), DIMENSION(:),   INTENT(IN) :: ipH2      ! inner product vector, to be saved
      REAL(pr), DIMENSION(:),   INTENT(IN) :: ipHn1     ! inner product vector, to be saved
      CHARACTER(200)                       :: filename                                        ! Filename for writing the file
      INTEGER, DIMENSION(1:2)              :: starts, counts, dimids                          ! Positioning for saving
      INTEGER                              :: ncout, ncid                                     ! Temporary integers for saving
      INTEGER                              :: varid, varid0, varidP, varidE, varidK, varidt   ! Temporary integers for saving
      INTEGER                              :: varidLONE, varidHONE, varidHTWO, varidHNONE     ! Temporary integers for saving
      INTEGER                              :: x_dimid, y_dimid, t_dimid, j_dimid, p_dimid     ! Temporary dimensions for saving
      INTEGER                              :: ii, t_len                                       ! Temporary integer for looping and reference length

      ! Filename path for saving file
      filename = TRIM(work_pathname)//"DNS_"//IC_type//"_Norm"//normconstr//"_Grad"//Grad_type//"_N"//Nchar//"_L"//TRIM(ADJUSTL(lchar))//"_T"//TRIM(ADJUSTL(tchar))//".nc"

      t_len = size(t)
      ! Initialize file for writing
      IF (rank==0) THEN
        ncout = nf90_create(filename, NF90_CLOBBER, ncid=ncid)
        IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
        ncout = nf90_def_dim(ncid, "x", n_nse(1), x_dimid)
        IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
        ncout = nf90_def_dim(ncid, "y", NF90_UNLIMITED, y_dimid)
        IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
        ncout = nf90_def_dim(ncid, 't', t_len, t_dimid)
        IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
        dimids = (/ x_dimid, y_dimid /)

        ncout = nf90_def_var(ncid, TRIM("w0"), NF90_DOUBLE, dimids, varid0)
        IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)

        ncout = nf90_def_var(ncid, TRIM("Palin"), NF90_DOUBLE, (/ t_dimid /), varidP)
        IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
        ncout = nf90_def_var(ncid, TRIM("Enst"), NF90_DOUBLE, (/ t_dimid /), varidE)
        IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
        ncout = nf90_def_var(ncid, TRIM("Kin"), NF90_DOUBLE, (/ t_dimid /), varidK)
        IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
        ncout = nf90_def_var(ncid, TRIM("InnerProduct_L2"), NF90_DOUBLE, (/ t_dimid /), varidLONE)
        IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
        ncout = nf90_def_var(ncid, TRIM("InnerProduct_H1"), NF90_DOUBLE, (/ t_dimid /), varidHONE)
        IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
        ncout = nf90_def_var(ncid, TRIM("InnerProduct_H2"), NF90_DOUBLE, (/ t_dimid /), varidHTWO)
        IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
        ncout = nf90_def_var(ncid, TRIM("InnerProduct_Hn1"), NF90_DOUBLE, (/ t_dimid /), varidHNONE)
        IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)

        ncout = nf90_def_var(ncid, TRIM("tvec"), NF90_DOUBLE, (/ t_dimid /), varidt)
        IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
        ncout = nf90_enddef(ncid)

        IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
        ncout = nf90_put_var(ncid, varidP, Pal)
        IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
        ncout = nf90_put_var(ncid, varidE, Ens)
        IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
        ncout = nf90_put_var(ncid, varidK, KinEnerg)
        IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
        ncout = nf90_put_var(ncid, varidLONE, ipL2)
        IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
        ncout = nf90_put_var(ncid, varidHONE, ipH1)
        IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
        ncout = nf90_put_var(ncid, varidHTWO, ipH2)
        IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
        ncout = nf90_put_var(ncid, varidHNONE, ipHn1)
        IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
        ncout = nf90_put_var(ncid, varidt, t)
        IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)

        ncout = nf90_close(ncid)
        IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
      END IF
      CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)

      ! Positioning
      starts = (/ 1, rank*local_Ny+1 /)
      counts = (/ n_nse(1), local_Ny /)

      !!--------------------------
      !! START netCDF ROUTINES
      !!--------------------------
      ! Initial condition
      DO ii=0,np-1
        IF (rank==ii) THEN
          ncout = nf90_open(filename, NF90_WRITE, ncid)
          IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
          ncout = nf90_inq_varid(ncid, TRIM("w0"), varid0)
          IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
          ncout = nf90_put_var(ncid, varid0, w0, start = starts, count = counts)
          IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
          ncout = nf90_close(ncid)
          IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
        END IF
        CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)
      END DO
      CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)

    END SUBROUTINE save_NS_DNS

    !==================================================================
    ! *** Read initial condition from MATLAB ***
    ! Input:   myIC  - type of initial guess from optimization
    !         size_x - size of data in x direction
    !         size_y - size of data in y direction
    ! Output: mydata - vorticity field, in physical space
    !==================================================================
    SUBROUTINE read_IC(myIC, size_x, size_y, mydata)
      ! Load variables
      USE global_variables, ONLY: pr, rank, n_nse, local_Ny, RESOL
      ! Initialize variables
      REAL(pr), DIMENSION(:,:), INTENT(OUT) :: mydata         ! Vorticity in physical space
      CHARACTER(len=*),          INTENT(IN) :: myIC           ! Description string
      INTEGER,                   INTENT(IN) :: size_x, size_y ! Size of data in x and y
      INTEGER                               :: i,j,k          ! Integers for looping
      CHARACTER(4)                          :: Nchar          ! Resolution as character
      CHARACTER(200)                        :: filename       ! Filename for writing the file


      ! Write resolution as a character
      WRITE(Nchar, '(i4.4)') RESOL
      ! Filename path for reading in project folder, for random initial condition
      filename = TRIM("/project/def-bprotas/mathap1/randIC/NS_"//myIC//"IC_N")//Nchar//".dat"

      ! Open file of random initial condition from MATLAB
      OPEN (20, FILE = filename, STATUS = 'OLD')

      ! Read data into appropriate location of the initial vorticity matrix
      ! ***This method of reading the data is inefficient, and there are more efficient ways to do this!
      DO k = 0, rank
        DO j=1,size_y ! Loop through columns to read
          ! Read data from file into appropriate location of vorticity
          ! Note for rank > 0, the data is written over until appropriate spot is reached
          ! This aspect makes the read inefficient!
          READ(20, *) ( mydata(i,j), i=1,size_x) ! Read in the x data (entire column) at once
        END DO
      END DO
      ! Close file
      CLOSE(20)
    END SUBROUTINE read_IC

    !==================================================================
    ! *** Read optimal initial condition ***
    ! Input:   myIC  - type of initial guess from optimization
    ! Output:   w    - optimal vorticity field in physical space
    !==================================================================
    SUBROUTINE read_NS_Opt( myIC, w )
      ! Load variables
      USE global_variables    ! Declares variables and defines parameters of the system
      USE netcdf              ! Use netcdf for saving files
      USE mpi                 ! Use MPI module (binding works well with fftw libraries)
      ! Initialize variables
!      REAL(pr), DIMENSION(1:n_nse(1),1:local_Ny), INTENT(IN) :: w         ! Vorticity field, to be saved
      CHARACTER(len=*),          INTENT(IN) :: myIC           ! Description string
      REAL(pr), DIMENSION(:,:), INTENT(OUT) :: w              ! Optimal initial condition, in physical space
      CHARACTER(200)                        :: filename       ! Filename for writing the file
      INTEGER, DIMENSION(1:2)               :: starts, counts ! Positioning for saving
      INTEGER                               :: ncout, ncid    ! Temporary integers for saving
      INTEGER                               :: varid          ! Temporary integers for saving
      INTEGER                               :: ii             ! Temporary integer for looping

      ! Filename path for reading file
      filename = TRIM(work_pathname)//"Optimization_"//myIC//"_Norm"//normconstr//"_Grad"//Grad_type//"_N"//Nchar//"_NU"//TRIM(ADJUSTL(viscchar))//"_L"//TRIM(ADJUSTL(lchar))//"_T"//TRIM(ADJUSTL(tchar))//".nc"

      CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)

      ! Positioning
      starts = (/ 1, rank*local_Ny+1 /)
      counts = (/ n_nse(1), local_Ny /)

      !!--------------------------
      !! START netCDF ROUTINES
      !!--------------------------
      ! Optimal field
      DO ii=0,np-1
        IF (rank==ii) THEN
          ncout = nf90_open(filename, NF90_NOWRITE, ncid)
          IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
          ncout = nf90_inq_varid(ncid, TRIM("w"), varid)
          IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
          ncout = nf90_get_var(ncid, varid, w, start = starts, count = counts)
          IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
          ncout = nf90_close(ncid)
          IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
        END IF
        CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)
      END DO

    END SUBROUTINE read_NS_Opt

    !==================================================================
    ! *** Read optimal initial condition for bootstrapping ***
    ! Input:   myIC  - type of initial guess from optimization
    ! Output:   w    - optimal vorticity field in physical space
    !==================================================================
    SUBROUTINE read_BS_Opt( myIC, w )
      ! Load variables
      USE global_variables!, ONLY: pr, n_nse, visc, local_Ny, work_pathname, IC_type
      USE netcdf              ! Use netcdf for saving files
      USE mpi                 ! Use MPI module (binding works well with fftw libraries)
      ! Initialize variables
!      REAL(pr), DIMENSION(1:n_nse(1),1:local_Ny), INTENT(IN) :: w         ! Vorticity field, to be saved
      CHARACTER(len=*),          INTENT(IN) :: myIC           ! Description string
      REAL(pr), DIMENSION(:,:), INTENT(OUT) :: w              ! Optimal initial condition, in physical space
      CHARACTER(200)                        :: filename       ! Filename for writing the file
      CHARACTER(4)                          :: NcharP         ! Resolution as character
      CHARACTER(8)                          :: lcharP         ! Sobolev parameter as character
      CHARACTER(10)                         :: visccharP      ! Viscosity as character
      CHARACTER(10)                         :: tcharP         ! Final time as character
      INTEGER, DIMENSION(1:2)               :: starts, counts ! Positioning for saving
      INTEGER                               :: ncout, ncid    ! Temporary integers for saving
      INTEGER                               :: varid          ! Temporary integers for saving
      INTEGER                               :: ii             ! Temporary integer for looping

      ! Previous resolution as a character
      WRITE(NcharP, '(i4.4)') RESOLP
      ! Previous viscosity as a character
      WRITE(visccharP, '(ES10.2)') viscP
      ! Previous Sobolev parameter as a character
      WRITE(lcharP, '(ES8.1)') ellP
      ! Previous final time as a character
      WRITE(tcharP, '(ES10.2)') endTimeP
      ! Filename path for saving file
      filename = TRIM(work_pathname)//"Optimization_"//myIC//"_Norm"//normconstr//"_Grad"//Grad_type//"_N"//NcharP//"_NU"//TRIM(ADJUSTL(visccharP))//"_L"//TRIM(ADJUSTL(lcharP))//"_T"//TRIM(ADJUSTL(tcharP))//".nc"

      CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)

      ! Positioning
      starts = (/ 1, rank*local_Ny+1 /)
      counts = (/ n_nse(1), local_Ny /)

      !!--------------------------
      !! START netCDF ROUTINES
      !!--------------------------
      ! Optimal field
      DO ii=0,np-1
        IF (rank==ii) THEN
          ncout = nf90_open(filename, NF90_NOWRITE, ncid)
          IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
          ncout = nf90_inq_varid(ncid, TRIM("w"), varid)
          IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
          ncout = nf90_get_var(ncid, varid, w, start = starts, count = counts)
          IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
          ncout = nf90_close(ncid)
          IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
        END IF
        CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)
      END DO

    END SUBROUTINE read_BS_Opt

    !==================================================================
    ! *** Save vorticity ***
    ! Input:     w   - vorticity field in physical space
    !        myindex - current time iteration
    !        mytype  - descriptor
    !==================================================================
    SUBROUTINE save_NS_vorticity(w, myindex, mytype)
      ! Load variables
      USE global_variables, ONLY: pr, n_nse, RESOL, visc, Lx, Ly, local_Ny, work_pathname, IC_type
      ! Initialize variables
!      REAL(pr), DIMENSION(1:n_nse(1),1:local_Ny), INTENT(IN) :: w         ! Vorticity field, to be saved
      REAL(pr), DIMENSION(:,:),                   INTENT(IN) :: w         ! Vorticity field, to be saved
      CHARACTER(len=*),                           INTENT(IN) :: mytype    ! Description string
      INTEGER,                                    INTENT(IN) :: myindex   ! Current time iteration
      CHARACTER(4)                                           :: Nchar     ! Resolution as character
      CHARACTER(6)                                           :: indexchar ! Time iteration as character
      !CHARACTER(13)                                          :: viscchar  ! Viscosity as character
      CHARACTER(200)                                         :: filename  ! Filename for writing the file

      ! Resolution as a character
      WRITE(Nchar, '(i4.4)') RESOL
      ! Iteration number as a character
      WRITE(indexchar, '(i6.6)') myindex
      ! Viscosity as a character
      !WRITE(viscchar, '(i13.13)') int(visc*(1.0e12)) ! 2DNS
      ! Filename path for saving in the scratch folder, for current timestep
      filename = TRIM(work_pathname)//"Vorticity_"//IC_type//"_N"//Nchar//"_"//mytype//"_"//indexchar//".nc"
      ! Save the vorticity as a netCDF file
      CALL save_field_R2toR1_ncdf(w, "w", filename)
    END SUBROUTINE save_NS_vorticity

    !==================================================================
    ! *** Save values from forward solve ***
    ! Input:   w     - vorticity field in physical space
    !       KinEnerg - kinetic energy vector
    !         Ens    - enstrophy vector
    !           t    - time vector
    !        myindex - current time iteration
    !==================================================================
    SUBROUTINE save_NS_fwd(w, KinEnerg, Ens, Pal, t, myindex)
      ! Load variables
      USE global_variables, ONLY: pr, n_nse, RESOL, visc, Lx, Ly, local_Ny, work_pathname, IC_type
      ! Initialize variables
!      REAL(pr), DIMENSION(1:n_nse(1),1:local_Ny), INTENT(IN) :: w         ! Vorticity field, to be saved
      REAL(pr), DIMENSION(:,:),                   INTENT(IN) :: w         ! Vorticity field, to be saved
      REAL(pr), DIMENSION(:),                     INTENT(IN) :: KinEnerg  ! Kinetic energy vector, to be saved
      REAL(pr), DIMENSION(:),                     INTENT(IN) :: Ens       ! Enstrophy vector, to be saved
      REAL(pr), DIMENSION(:),                     INTENT(IN) :: Pal       ! Palinstrophy vector, to be saved
      REAL(pr), DIMENSION(:),                     INTENT(IN) :: t         ! Time vector, to be saved
      INTEGER,                                    INTENT(IN) :: myindex   ! Current time iteration
      CHARACTER(4)                                           :: Nchar     ! Resolution as character
      CHARACTER(6)                                           :: indexchar ! Time iteration as character
      CHARACTER(13)                                          :: viscchar  ! Viscosity as character
      CHARACTER(200)                                         :: filename  ! Filename for writing the file

      ! Resolution as a character
      WRITE(Nchar, '(i4.4)') RESOL
      ! Iteration number as a character
      WRITE(indexchar, '(i6.6)') myindex
      ! Viscosity as a character
      WRITE(viscchar, '(i13.13)') int(visc*(1.0e12))
      ! Filename path for saving vorticity
      filename = TRIM(work_pathname)//"Vorticity_"//IC_type//"_N"//Nchar//"_NU"//viscchar//"_FWDFinal_"//indexchar//".nc"
      CALL save_field_R2toR1_ncdf(w, "w", filename) ! Save vorticity field
      ! Filename path for saving enstrophy
      filename = TRIM(work_pathname)//"Vorticity_"//IC_type//"_N"//Nchar//"_NU"//viscchar//"_Kin_"//indexchar//".nc"
      CALL save_field_R1toR1_ncdf(KinEnerg, "Kin", filename) ! Save Kinetic Energy
      ! Filename path for saving enstrophy
      filename = TRIM(work_pathname)//"Vorticity_"//IC_type//"_N"//Nchar//"_NU"//viscchar//"_Enst_"//indexchar//".nc"
      CALL save_field_R1toR1_ncdf(Ens, "Enst", filename) ! Save Enstrophy
      ! Filename path for saving enstrophy
      filename = TRIM(work_pathname)//"Vorticity_"//IC_type//"_N"//Nchar//"_NU"//viscchar//"_Palin_"//indexchar//".nc"
      CALL save_field_R1toR1_ncdf(Pal, "Palin", filename) ! Save Palinstrophy
      ! Filename path for saving time vector
      filename = TRIM(work_pathname)//"Vorticity_"//IC_type//"_N"//Nchar//"_NU"//viscchar//"_t_"//indexchar//".nc"
      CALL save_field_R1toR1_ncdf(t, "tvec", filename)  ! Save time vector
    END SUBROUTINE save_NS_fwd

    !==================================================================
    ! *** Save 2D scalar field in netCDF format ***
    ! Input:   myfield  - 2D field (vorticity) in physical space
    !        field_name - name of the field to save
    !        file_name  - string for saving file
    !==================================================================
    SUBROUTINE save_field_R2toR1_ncdf(myfield, field_name, file_name)
      ! Load variables
      USE global_variables    ! Declares variables and defines parameters of the system
      ! Load subroutines
      USE netcdf              ! Use netcdf for saving files
      USE mpi                 ! Use MPI module (binding works well with fftw libraries)
      ! Initialize variables
!      REAL(pr), DIMENSION(1:n_nse(1),1:local_Ny), INTENT(IN) :: myfield                       ! 2D scalar field, in physical space
      REAL(pr), DIMENSION(:,:), INTENT(IN)                   :: myfield                       ! 2D scalar field, in physical space
      CHARACTER(len=*)                                       :: field_name                    ! Name of field
      CHARACTER(len=*)                                       :: file_name                     ! Filename
      INTEGER, DIMENSION(1:2)                                :: starts, counts, dimids        ! Positioning for saving
      INTEGER                                                :: ncout, ncid, varid            ! Temporary integers for saving
      INTEGER                                                :: x_dimid, y_dimid              ! Temporary dimensions for saving
      INTEGER                                                :: ii                            ! Temporary integer for looping

      ! Initialize file for writing
      IF (rank==0) THEN
        ncout = nf90_create(file_name, NF90_CLOBBER, ncid=ncid)
        IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
        ncout = nf90_def_dim(ncid, "x", n_nse(1), x_dimid)
        IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
        ncout = nf90_def_dim(ncid, "y", NF90_UNLIMITED, y_dimid)
        IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
        dimids = (/ x_dimid, y_dimid /)

        ncout = nf90_def_var(ncid, TRIM(field_name), NF90_DOUBLE, dimids, varid)
        IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)

        ncout = nf90_enddef(ncid)
        IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
        ncout = nf90_close(ncid)
        IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
      END IF
      CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)

      ! Positioning
      starts = (/ 1, rank*local_Ny+1 /)
      counts = (/ n_nse(1), local_Ny /)

      !!--------------------------
      !! START netCDF ROUTINES
      !!--------------------------
      DO ii=0,np-1
        IF (rank==ii) THEN
          ncout = nf90_open(file_name, NF90_WRITE, ncid)
          IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)

          ncout = nf90_inq_varid(ncid, TRIM(field_name), varid)
          IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
          ncout = nf90_put_var(ncid, varid, myfield, start = starts, count = counts)
          IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)

          ncout = nf90_close(ncid)
          IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
        END IF
        CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)
      END DO
      CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)
    END SUBROUTINE save_field_R2toR1_ncdf

    !==================================================================
    ! *** Save 1D scalar field in netCDF format ***
    ! Input:   myfield  - 1D field (vector, function of time)
    !        field_name - name of the field to save
    !        file_name  - string for saving file
    !==================================================================
    SUBROUTINE save_field_R1toR1_ncdf(myfield, field_name, file_name)
      ! Load variables
      USE global_variables    ! Declares variables and defines parameters of the system
      ! Load subroutines
      USE netcdf              ! Use netcdf for saving files
      USE mpi                 ! Use MPI module (binding works well with fftw libraries)
      ! Initialize variables
      REAL(pr), DIMENSION(:), INTENT(IN) :: myfield            ! 1D scalar field
      CHARACTER(len=*)                   :: field_name         ! Name of field
      CHARACTER(len=*)                   :: file_name          ! Filename
      INTEGER                            :: ncout, ncid, varid ! Temporary integers for saving
      INTEGER                            :: t_dimid            ! Temporary dimensions for saving
      INTEGER                            :: ii                 ! Temporary integer for size

      ! Assign size
      ii = size(myfield)
      ! Write to file with only one processor
      IF (rank == 0) THEN
        ncout = nf90_create(file_name, NF90_CLOBBER, ncid)
        IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
        ncout = nf90_def_dim(ncid, 't', ii, t_dimid)
        IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
        ncout = nf90_def_var(ncid, field_name, NF90_DOUBLE, (/ t_dimid /), varid)
        IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
        ncout = nf90_enddef(ncid)
        IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
        ncout = nf90_put_var(ncid, varid, myfield)
        IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
        ncout = nf90_close(ncid)
        IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
      END IF
      CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)
    END SUBROUTINE save_field_R1toR1_ncdf

    !==================================================================
    ! *** NETCDF Error Handle Routine ***
    ! Input: nerror - Error code
    !==================================================================
    SUBROUTINE ncdf_error_handle(nerror)
      USE netcdf ! Use netcdf for saving files
      ! Initialize variables
      INTEGER, INTENT(IN) :: nerror       ! Error code
      CHARACTER(80)       :: error_string ! Error code as a string

      PRINT *, " Printing Error. "
      ! Translate error code to a string
      error_string = NF90_STRERROR(nerror)
      ! Print error code
      PRINT *, " Error reading netCDF file. "//error_string
    END SUBROUTINE ncdf_error_handle


END MODULE
