!###############################################################################
!#
!# Useful netCDF functions
!#
!# ELR 18/02/2014
!#
!###############################################################################

MODULE netcdf_mod

  USE netcdf

  IMPLICIT none

  PUBLIC nc_err_check, err_check

CONTAINS

!###############################################################################
!# netCDF error check (from DBC's in regrid_scrip) 
!###############################################################################
  SUBROUTINE nc_err_check( status,message,callFrom )

! Stop if an error is found by a netCDF call.

    IMPLICIT NONE
    INTEGER, INTENT(in) :: status  !  error value from netCDF
    CHARACTER (len=*), OPTIONAL, INTENT(in) :: message
    CHARACTER (len=*), OPTIONAL, INTENT(in) :: callFrom
    CHARACTER(len=50), PARAMETER :: procName = 'nc_err_check'  !  name of this procedure

    IF ( status /= nf90_noerr ) THEN
      WRITE (*,*) TRIM( nf90_strerror(status) )
      IF ( PRESENT(message) ) WRITE(*,*) TRIM(message)
      IF ( PRESENT(callFrom) ) WRITE(*,*) 'Called from ',TRIM(callFrom)
      WRITE(*,*)'Stopping in ',TRIM(procName)
      STOP 2 
    END IF

  END SUBROUTINE nc_err_check

!###############################################################################
!# error check (from DBC's in regrid_scrip) 
!###############################################################################
  SUBROUTINE err_check( status,message,callFrom )

! Stop if an error is found 

    IMPLICIT NONE
    INTEGER, INTENT(in) :: status  !  error value 
    CHARACTER (len=*), OPTIONAL, INTENT(in) :: message
    CHARACTER (len=*), OPTIONAL, INTENT(in) :: callFrom
    CHARACTER(len=50), PARAMETER :: procName = 'err_check'  !  name of this procedure

    IF ( status /= 0 ) THEN
      IF ( PRESENT(message) ) WRITE(*,*) TRIM(message)
      IF ( PRESENT(callFrom) ) WRITE(*,*) 'Called from ',TRIM(callFrom)
      WRITE(*,*)'Stopping in ',TRIM(procName)
      STOP 3
    END IF

  END SUBROUTINE err_check

!###############################################################################
!# create netcdf file of appropriate formate (copied from 
!# ~/Data-Area/projects/PAGODA/simple_model/fortran/netcdf_mod.f90
!# then edited massively)
!###############################################################################

  SUBROUTINE create_nc_file( ndimMax,npoints,nx,ny,nz,nvar,missingValue,       &
&         daily,monthly,charYear,fileName,title,            &
&         varNdim,varNlev,latitude,longitude,varDims,                &
&         varLongName,varName,varUnits,ncID,timeVarID,varID,regLatLon,useT )

! Create a new netCDF dataset, clobbering any existing file.
! Variables are expected to be REAL.
! The names of the dimensions etc are set for WaterMIP runs and might well have
! to be changed for other experiments.

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! Scalar arguments with intent(in)
!-------------------------------------------------------------------------------
  INTEGER, INTENT(in) :: ndimMax    !  maximum allowed number of dimensions
  INTEGER, INTENT(in) :: npoints    !  number of defined points (only used if compress=T?)
  INTEGER, INTENT(in) :: nx         !  size of grid in x
  INTEGER, INTENT(in) :: ny         !  size of grid in y
  INTEGER, INTENT(in) :: nz         !  size of grid in z
  INTEGER, INTENT(in) :: nvar       !  number of variables to be created (not co-ord vars)
  REAL, INTENT(in) :: missingValue  !  missing data value
  LOGICAL, INTENT(in) :: daily      !  TRUE means data are daily
  LOGICAL, INTENT(in) :: monthly    !  TRUE means data are monthly
  LOGICAL, INTENT(in) :: regLatLon  !  TRUE if grid is a lat-lon grid
!                            This is not currently checked! Could be input.
  LOGICAL, INTENT(in) :: useT       !  TRUE if we have time-varying data
  CHARACTER(len=*), INTENT(in) :: charYear  !  year (yyyy)
  CHARACTER(len=*), INTENT(in) :: fileName  !  name of dataset to be created
  CHARACTER(len=*), INTENT(in) :: title     !  used for title attribute

!-------------------------------------------------------------------------------
! Array arguments with intent(in)
!-------------------------------------------------------------------------------
  INTEGER, INTENT(in) :: varNdim(nvar)  !  number of dimensions in each variable
  INTEGER, INTENT(in) :: varNlev(nvar)  !  number of z levels for each variable
  REAL, INTENT(in) :: latitude(nx,ny)     !  latitude of each point
  REAL, INTENT(in) :: longitude(nx,ny)     !  longitude of each point
  CHARACTER(len=*), INTENT(in) :: varDims(nvar,ndimMax)  !  dimensions of each variable
  CHARACTER(len=*), INTENT(in) :: varLongName(nvar)      !  long name of each variable
  CHARACTER(len=*), INTENT(in) :: varName(nvar)          !  long name of each variable
  CHARACTER(len=*), INTENT(in) :: varUnits(nvar)         !  units of each variable

!-------------------------------------------------------------------------------
! Scalar arguments with intent(out)
!-------------------------------------------------------------------------------
  INTEGER, INTENT(out) :: ncID       !  netCDF ID of newly-created dataset
  INTEGER, INTENT(out) :: timeVarID  !  netCDF ID of time variable

!-------------------------------------------------------------------------------
! Array arguments with intent(out)
!-------------------------------------------------------------------------------
  INTEGER, INTENT(out) :: varID(nvar)  !  netCDF ID of each variable

!-------------------------------------------------------------------------------
! Optional array arguments with intent(in).
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Local scalar parameters.
!-------------------------------------------------------------------------------
  CHARACTER(len=20), PARAMETER :: procName = 'create_nc_file'  !  name of this procedure

!-------------------------------------------------------------------------------
! Local scalar variables.
!-------------------------------------------------------------------------------
  INTEGER :: idim,ivar,ix,iy  !  work
  INTEGER :: indexID    !  netCDF ID for index (e.g. land) variable
  INTEGER :: latID      !  netCDF ID for latitude variable
  INTEGER :: lonID      !  netCDF ID for longitude variable
  INTEGER :: ndim       !  work
  INTEGER :: nzMax      !  maximum number of levels in any one variable
  INTEGER :: indexDimID  !  netCDF ID for index dimension
  INTEGER :: tDimID     !  netCDF ID for time dimension
  INTEGER :: xDimID     !  netCDF ID for x dimension
  INTEGER :: yDimID     !  netCDF ID for y dimension
  INTEGER :: zDimID     !  netCDF ID for z dimension

  LOGICAL :: useZ       !  TRUE if we have a multi-level variable

  CHARACTER(len=100) :: cwork  !  workspace

!-------------------------------------------------------------------------------
! Local array variables.
!-------------------------------------------------------------------------------
  INTEGER :: dimIDs(ndimMax)  !  IDs of dimensions needed for a particular variable
!                                These are currently x,y,z,t.

!-------------------------------------------------------------------------------
! Check arguments.
!-------------------------------------------------------------------------------
  IF ( daily .AND. monthly ) THEN
    WRITE(*,*)'ERROR in ',TRIM(procName)
    WRITE(*,*)'At most one of daily and monthly can be TRUE.'
    STOP 2
  END IF

!-------------------------------------------------------------------------------
! Decide what dimensions we need.
!-------------------------------------------------------------------------------
  useZ = .FALSE.
  nzMax = 1

! Currently assuming no z levels.
! Look for variables with more than one level.
! Loop over all output variables, looking for those in this file.
  DO ivar=1,nvar
    IF ( varNlev(ivar) > nzMax ) THEN
      useZ = .TRUE.
      nzMax = varNlev(ivar)
    END IF
  END DO

!-------------------------------------------------------------------------------
! Get dataset ID.
!-------------------------------------------------------------------------------
  CALL nc_err_check( nf90_create( fileName,nf90_clobber,ncID ),                &
&                   'nf90_create '//TRIM(fileName),procName)

!-------------------------------------------------------------------------------
! Create an unlimited dimension.
!-------------------------------------------------------------------------------
  IF (useT) CALL nc_err_check( nf90_def_dim( ncID,'time',nf90_unlimited,tDimID ),        &
&                   'nf90_def_dim time',procName )

!-------------------------------------------------------------------------------
! Create other dimensions.
! Currently assuming a (regular) lat-lon grid.
!-------------------------------------------------------------------------------
  CALL nc_err_check( nf90_def_dim( ncID,'x',nx,xDimID ),                   &
&                   'nf90_def_dim lon',procName )

  CALL nc_err_check( nf90_def_dim( ncID,'y',ny,yDimID ),                   &
&                   'nf90_def_dim lat',procName )

  IF ( useZ ) CALL nc_err_check( nf90_def_dim( ncID,'z',nzMax,zDimID ),        &
&                               'nf90_def_dim z',procName )

!-------------------------------------------------------------------------------
! Create dimension variables and attributes.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Time.
!-------------------------------------------------------------------------------

! Get ID for variable.
  IF (useT) THEN
    CALL nc_err_check( nf90_def_var( ncID,'time',nf90_float,(/tDimID/),timeVarID ),                       &
&                   'nf90_def_var time',procName )

! Write attributes.
    IF ( monthly ) THEN
      cwork = 'months since ' // charYear // '-01-01'
    ELSE IF ( daily ) THEN
!   Assuming daily file starts first of month.
      cwork = 'days since ' // charYear // '-01-01 00:00:00'
    ELSE
      WRITE(*,*)'No attribute code for this time interval.'
      WRITE(*,*)'Stopping in ',TRIM(procName)
      STOP 2
    END IF
    CALL nc_err_check( nf90_put_att( ncID,timeVarID,'units',TRIM(cwork) ),       &
  &                   'nf90_put_att time',procName )

    CALL nc_err_check( nf90_put_att( ncID,timeVarID,'title','Time' ),            &
  &                   'nf90_put_att time',procName )

    CALL nc_err_check( nf90_put_att( ncID,timeVarID,'long_name','Time axis' ),   &
  &                   'nf90_put_att time',procName )
  END IF

!-------------------------------------------------------------------------------
! Longitude.
!-------------------------------------------------------------------------------
! Get ID for variable.
  IF ( regLatLon ) THEN
    CALL nc_err_check( nf90_def_var( ncID,'lon',nf90_float,(/xDimID/),lonID ), &
&                   'nf90_def_var lon',procName )
  ELSE
    CALL nc_err_check( nf90_def_var( ncID,'lon',nf90_float,                    &
&                    (/xDimID,yDimID/),lonID ), &
&                   'nf90_def_var lon',procName )
  END IF

! Write attributes.
  CALL nc_err_check( nf90_put_att( ncID,lonID,'units','degrees_east' ),      &
&                   'nf90_put_att lon',procName )

  CALL nc_err_check( nf90_put_att( ncID,lonID,'valid_min',MINVAL(longitude(:,1)) ),&
&                   'nf90_put_att lon',procName )

  CALL nc_err_check( nf90_put_att( ncID,lonID,'valid_max',MAXVAL(longitude(:,1)) ),&
&                   'nf90_put_att lon',procName )

  CALL nc_err_check( nf90_put_att( ncID,lonID,'long_name','Longitude' ),     &
&                   'nf90_put_att lon',procName )

!-------------------------------------------------------------------------------
! Latitude.
!-------------------------------------------------------------------------------
! Get ID for variable.
  IF ( regLatLon ) THEN
    CALL nc_err_check( nf90_def_var( ncID,'lat',nf90_float,(/yDimID/),latID ), &
&                   'nf90_def_var lat',procName )
  ELSE
    CALL nc_err_check( nf90_def_var( ncID,'lat',nf90_float,                    &
&                    (/xDimID,yDimID/),latID ), &
&                   'nf90_def_var lat',procName )
  END IF

! Write attributes.
  CALL nc_err_check( nf90_put_att( ncID,latID,'units','degrees_north' ),       &
&                   'nf90_put_att lat',procName )

  CALL nc_err_check( nf90_put_att( ncID,latID,'valid_min',MINVAL(latitude(1,:)) ),  &
&                   'nf90_put_att lat',procName )

  CALL nc_err_check( nf90_put_att( ncID,latID,'valid_max',MAXVAL(latitude(1,:)) ),  &
&                   'nf90_put_att lat',procName )

  CALL nc_err_check( nf90_put_att( ncID,latID,'long_name','Latitude' ),        &
&                   'nf90_put_att lat',procName )

!-------------------------------------------------------------------------------
! Global attributes. Hardwired!
!-------------------------------------------------------------------------------
  CALL nc_err_check( nf90_put_att( ncID,nf90_global,'title',TRIM(title) ),                       &
&                   'nf90_put_att global',procName )  

!-------------------------------------------------------------------------------
! Create data variables and attributes.
!-------------------------------------------------------------------------------

! Loop over all output variables.
  DO ivar=1,nvar

!   Work out what dimensions this variable needs.
!   Note that dimension names are hardwired here.
!   Minimum is x,t.
    dimIDs(:) = 0   !  for clarity
    ndim = 0
    DO idim=1,varNDim(ivar)
      IF ( varDims(ivar,idim) == 'x' ) THEN
        ndim = ndim + 1
        dimIDs(ndim) = xDimID
      ELSE IF ( varDims(ivar,idim) == 'y' ) THEN
        ndim = ndim + 1
        dimIDs(ndim) = yDimID
      ELSE IF ( varDims(ivar,idim) == 'z' ) THEN
        ndim = ndim + 1
        dimIDs(ndim) = zDimID
      ELSE IF ( varDims(ivar,idim) == 't' ) THEN
        ndim = ndim + 1
        dimIDs(ndim) = tDimID
      END IF
    END DO  !  idim

!   Get ID for variable.
    CALL nc_err_check( nf90_def_var( ncID,TRIM(varName(ivar)),nf90_float     &
&                       ,dimIDs(1:ndim),varID(ivar) ),                      &
&                     'nf90_def_var '//TRIM((varName(ivar))),procName )

!   Write attributes.
    CALL nc_err_check( nf90_put_att( ncID,varID(ivar),'units',TRIM(varUnits(ivar))), &
&                     'nf90_put_att '//TRIM((varName(ivar))),procName )

    CALL nc_err_check( nf90_put_att( ncID,varID(ivar),'missing_value',missingValue ), &
&                     'nf90_put_att '//TRIM((varName(ivar))),procName )

    CALL nc_err_check( nf90_put_att( ncID,varID(ivar),'long_name',TRIM(varLongName(ivar)) ), &
&                     'nf90_put_att '//TRIM((varName(ivar))),procName )

  END DO

!-------------------------------------------------------------------------------
! Finished definitions, leave define mode.
!-------------------------------------------------------------------------------
  CALL nc_err_check( nf90_enddef( ncID ), 'nf90_enddef',procName )

!-------------------------------------------------------------------------------
! Write lat and lon variables.
!-------------------------------------------------------------------------------
  IF ( regLatLon ) THEN
    CALL nc_err_check( nf90_put_var( ncID,lonID,longitude(:,1) ),              &
&                     'nf90_put_var longitude',procName )

    CALL nc_err_check( nf90_put_var( ncID,latID,latitude(1,:) ),               &
&                     'nf90_put_var latitude',procName )
  ELSE
    CALL nc_err_check( nf90_put_var( ncID,lonID,longitude(:,:) ),              &
&                     'nf90_put_var longitude',procName )


    CALL nc_err_check( nf90_put_var( ncID,latID,latitude(:,:) ),               &
&                     'nf90_put_var latitude',procName )
  END IF

  END SUBROUTINE create_nc_file

!###############################################################################
!# Write data with no time dimension
!# ~/Data-Area/projects/PAGODA/simple_model/fortran/netcdf_mod.f90
!###############################################################################

  SUBROUTINE nc_write_data( nx,ny,nz,ncID,varID,valVector )

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! Scalar arguments with intent(in)
!-------------------------------------------------------------------------------
  INTEGER, INTENT(in) :: nx,ny,nz !  size of output grid
  INTEGER, INTENT(in) :: ncID  !  netCDF dataset ID
  INTEGER, INTENT(in) :: varID !  netCDF variable ID

!-------------------------------------------------------------------------------
! Array arguments with intent(in)
!-------------------------------------------------------------------------------
  REAL, INTENT(in) :: valVector(nx,ny,nz) 

!-------------------------------------------------------------------------------
! Local scalar parameters.
!-------------------------------------------------------------------------------
  CHARACTER(len=20), PARAMETER :: procName = 'nc_write_data'  !  name of this procedure

!-------------------------------------------------------------------------------
!   Reshape to a grid and write to an (x,y,z) variable.
!-------------------------------------------------------------------------------
  CALL nc_err_check(                                                           &
&       nf90_put_var( ncID,varID,valVector,start=(/1,1,1/) ), &
&      'nf90_put_var',procName )
    

  END SUBROUTINE nc_write_data


!###############################################################################
!# close the file
!# ~/Data-Area/projects/PAGODA/simple_model/fortran/netcdf_mod.f90
!###############################################################################

  SUBROUTINE nc_close_input(ncIDin)

  IMPLICIT NONE

  INTEGER, INTENT(inout) :: ncIDin

  CHARACTER(len=20):: procName='nc_close'
!-------------------------------------------------------------------------------
! Close input file.
!-------------------------------------------------------------------------------
  CALL nc_err_check( nf90_close(ncIDin), 'nf90_close',procName)
  ncIDin = -1   !  reset so we know that no dataset is open
    
  END SUBROUTINE nc_close_input



END MODULE netcdf_mod
