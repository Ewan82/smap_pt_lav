!###############################################################################
!#
!# Read input files
!#
!# ELR 18/02/2014
!#
!###############################################################################

MODULE read_input_mod

  IMPLICIT none

  PUBLIC :: read_map, nx, ny, mapVar, MU_GLOBAL, SEQ, ISSOIL, T_SAND,   &
            T_SILT, T_CLAY, T_BULK_DENSITY, S_SAND, S_SILT, S_CLAY,           &
            S_BULK_DENSITY, read_params,              &
            VAN_G_N_IN , VAN_G_ALPHA_IN , VAN_G_THETAR_IN , VAN_G_THETAS_IN , &
            VAN_G_KSat_IN, latVar, lonVar
  


  SAVE

  INTEGER :: nx, ny
  INTEGER :: ncID   ! netCDF file id  

  REAL, ALLOCATABLE :: mapVar(:,:)
  REAL, ALLOCATABLE :: latVar(:,:)
  REAL, ALLOCATABLE :: lonVar(:,:)

  INTEGER, ALLOCATABLE :: MU_GLOBAL(:) , SEQ(:) , ISSOIL(:)
  REAL, ALLOCATABLE :: T_SAND(:) , T_SILT(:) , T_CLAY(:) , T_OC(:),           &
                       T_BULK_DENSITY(:),                                     &
                       S_SAND(:) , S_SILT(:) , S_CLAY(:) , S_OC(:),           &
                       S_BULK_DENSITY(:)
  
  REAL, ALLOCATABLE :: VAN_G_N_IN(:,:) , VAN_G_ALPHA_IN(:,:) ,                &
                       VAN_G_THETAR_IN(:,:) , VAN_G_THETAS_IN(:,:) ,          &
                       VAN_G_KSat_IN(:,:) 

  CHARACTER(len=10), ALLOCATABLE :: SU_SYM90(:) 
  
CONTAINS

!###############################################################################
!# Read netCDF map file
!###############################################################################
  SUBROUTINE read_map

    USE netcdf_mod
    USE init_mod, ONLY : inMapFileName, inMapVarName, l_LUT

    IMPLICIT none

    INTEGER :: dimID  ! netCDF dimension id  
    INTEGER :: varID  ! netCDF variable id  

    CHARACTER(len=8) :: procName='read_map'

    CHARACTER(len=3) :: inLatVarName = 'lat'
    CHARACTER(len=3) :: inLonVarName = 'lon'

    PRINT *, "Reading soil map"
!-------------------------------------------------------------------------------
! Open file
!-------------------------------------------------------------------------------
    WRITE (*,*) " Opening ", TRIM(inMapFileName)

    CALL nc_err_check( nf90_open( inMapFileName,nf90_nowrite,ncID ),          &
&                    'nf90_open '//TRIM(inMapFileName),procName )

!-------------------------------------------------------------------------------
! Get dimension lengths
!-------------------------------------------------------------------------------
    CALL nc_err_check( nf90_inq_dimid( ncID,'x',dimID ),          &
&                   'nf90_inq_dimid x',procName )
    CALL nc_err_check( nf90_inquire_dimension( ncID,dimID,len=nx),   &
&                   'nf90_inquire_dimension x',procName )

    CALL nc_err_check( nf90_inq_dimid( ncID,'y',dimID ),          &
&                   'nf90_inq_dimid y',procName )
    CALL nc_err_check( nf90_inquire_dimension( ncID,dimID,len=ny),   &
&                   'nf90_inquire_dimension y',procName )


!-------------------------------------------------------------------------------
! Allocate array(s)
!-------------------------------------------------------------------------------
    ALLOCATE(mapVar(nx,ny))

!-------------------------------------------------------------------------------
! Read data
!-------------------------------------------------------------------------------
    CALL nc_err_check( nf90_inq_varid(ncID,inMapVarName,varID),               &
&                   'nf90_inq_varid '//inMapVarName,procName )
    CALL nc_err_check( nf90_get_var( ncID,varID,mapVar(:,:)),                  &
&                   'nf90_get_var '//inMapVarName,procName )

!-------------------------------------------------------------------------------
! Read lat and lon if we need it
!-------------------------------------------------------------------------------
    IF (.NOT. l_LUT) THEN
      ALLOCATE(latVar(nx,ny))
      CALL nc_err_check( nf90_inq_varid(ncID,inLatVarName,varID),             &
&                   'nf90_inq_varid '//inLatVarName,procName )
    CALL nc_err_check( nf90_get_var( ncID,varID,latVar(:,:)),                 &
&                   'nf90_get_var '//inLatVarName,procName )

      ALLOCATE(lonVar(nx,ny))
      CALL nc_err_check( nf90_inq_varid(ncID,inLonVarName,varID),             &
&                   'nf90_inq_varid '//inLonVarName,procName )
    CALL nc_err_check( nf90_get_var( ncID,varID,LonVar(:,:)),                 &
&                   'nf90_get_var '//inLonVarName,procName )

    END IF


!-------------------------------------------------------------------------------
! Close file
!-------------------------------------------------------------------------------
    CALL nc_err_check( nf90_close( ncID ),'nf90_close',procName )

  END SUBROUTINE read_map

!###############################################################################
!# Read parameter table
!###############################################################################
  SUBROUTINE read_params

    USE init_mod, ONLY : inParamFileName, nMapUnits

    IMPLICIT none

    INTEGER :: ID , MU_SOURCE2 , SU_CODE74 , SU_CODE85 , SU_CODE90 ,          &
               T_TEXTURE , DRAINAGE , AWC_CLASS , PHASE1 , PHASE2 , ROOTS ,   &
               IL , SWR , ADD_PROP , REF_DEPTH 

    REAL :: SHARE , T_GRAVEL , T_REF_BULK_DENSITY ,                           &
            T_PH_H2O , T_CEC_CLAY , T_CEC_SOIL , T_BS , T_TEB , T_CACO3 ,     &
            T_CASO4 , T_ESP , T_ECE , S_GRAVEL , S_REF_BULK_DENSITY ,         &
            S_PH_H2O , S_CEC_CLAY , S_CEC_SOIL ,             &
            S_BS , S_TEB , S_CACO3 , S_CASO4 , S_ESP , S_ECE

    CHARACTER(len=10) :: MU_SOURCE1 , SU_SYM74 , SU_SYM85 ,        &
                         T_USDA_TEX_CLASS , S_USDA_TEX_CLASS

    INTEGER :: inUnitP

    ! Loop variables
    INTEGER :: i

    PRINT *, "Reading soil data"
!-------------------------------------------------------------------------------
! Allocate variables
!-------------------------------------------------------------------------------

    ! Variables that we actually use get an array
    ALLOCATE(MU_GLOBAL(nMapUnits)) 
    ALLOCATE(SEQ(nMapUnits)) 
    ALLOCATE(SU_SYM90(nMapUnits)) 
    ALLOCATE(ISSOIL(nMapUnits))
    ALLOCATE(T_SAND(nMapUnits)) 
    ALLOCATE(T_BULK_DENSITY(nMapUnits)) 
    ALLOCATE(T_SILT(nMapUnits)) 
    ALLOCATE(T_CLAY(nMapUnits)) 
    ALLOCATE(T_OC(nMapUnits)) 
    ALLOCATE(S_SAND(nMapUnits)) 
    ALLOCATE(S_BULK_DENSITY(nMapUnits)) 
    ALLOCATE(S_SILT(nMapUnits)) 
    ALLOCATE(S_CLAY(nMapUnits)) 
    ALLOCATE(S_OC(nMapUnits)) 

!-------------------------------------------------------------------------------
! Open file
!-------------------------------------------------------------------------------
    WRITE (*,*) " Opening ", TRIM(inParamFileName)
    OPEN( inUnitP, file=TRIM(inParamFileName), action='read', status='old')

    ! Skip the first line
    READ(inUnitP, *)

    DO i=1,nMapUnits
      READ(inUnitP , *) ID , MU_GLOBAL(i) , MU_SOURCE1 , MU_SOURCE2 ,         &
                       ISSOIL(i) , SHARE , SEQ(i) , SU_SYM74 , SU_CODE74 ,    &
                       SU_SYM85 , SU_CODE85 , SU_SYM90(i) , SU_CODE90 ,       &
                       T_TEXTURE , DRAINAGE , REF_DEPTH , AWC_CLASS ,         &
                       PHASE1 , PHASE2 , ROOTS , IL , SWR , ADD_PROP ,        &
                       T_GRAVEL , T_SAND(i) , T_SILT(i) , T_CLAY(i) ,         &
                       T_USDA_TEX_CLASS , T_REF_BULK_DENSITY ,                &
                       T_BULK_DENSITY(i) , T_OC(i) , T_PH_H2O , T_CEC_CLAY ,  &
                       T_CEC_SOIL , T_BS , T_TEB , T_CACO3 , T_CASO4 ,        &
                       T_ESP , T_ECE , S_GRAVEL , S_SAND(i) , S_SILT(i) ,     &
                       S_CLAY(i) , S_USDA_TEX_CLASS , S_REF_BULK_DENSITY ,    &
                       S_BULK_DENSITY(i) , S_OC(i) , S_PH_H2O , S_CEC_CLAY ,  &
                       S_CEC_SOIL , S_BS , S_TEB , S_CACO3 , S_CASO4 ,        &
                       S_ESP , S_ECE
    END DO

    CLOSE(inUnitP)
  ENDSUBROUTINE read_params

!###############################################################################
!# Read parameter table
!###############################################################################
  SUBROUTINE read_vg_params

    USE init_mod, ONLY : inVGFileName, nSoilClass, nSoilInLay

    USE netcdf_mod, ONLY : err_check


    IMPLICIT none

    CHARACTER(len=14) :: procName='read_vg_params'

    INTEGER :: inUnitV

    ! Loop variables
    INTEGER :: i, j
    INTEGER :: soilClass

    PRINT *, "Reading van Genuchten look-up-table"
    IF (nSoilInLay<1 .OR. nSoilInLay>2)                                       &
      CALL err_check(1,'ERROR: Can have one or two input soil layers',procName)

!-------------------------------------------------------------------------------
! Allocate variables
!-------------------------------------------------------------------------------

    ! Variables that we actually use get an array
    ALLOCATE(VAN_G_N_IN(nSoilClass,nSoilInLay)) 
    ALLOCATE(VAN_G_ALPHA_IN(nSoilClass,nSoilInLay)) 
    ALLOCATE(VAN_G_THETAR_IN(nSoilClass,nSoilInLay)) 
    ALLOCATE(VAN_G_THETAS_IN(nSoilClass,nSoilInLay)) 
    ALLOCATE(VAN_G_KSat_IN(nSoilClass,nSoilInLay)) 

!-------------------------------------------------------------------------------
! Open file
!-------------------------------------------------------------------------------
    WRITE (*,*) " Opening ", TRIM(inVGFileName)
    OPEN( inUnitV, file=TRIM(inVGFileName), action='read', status='old')

!-------------------------------------------------------------------------------
! Read values
!-------------------------------------------------------------------------------
    DO i=1,nSoilClass
      DO j=1,nSoilInLay
        READ(inUnitV , *) soilClass, VAN_G_KSat_IN(i,j), VAN_G_THETAS_IN(i,j),&
                         VAN_G_THETAR_IN(i,j), VAN_G_ALPHA_IN(i,j),           &
                         VAN_G_N_IN(i,j)
        IF (soilClass/=i)                                                     &
          CALL err_check(1,'ERROR: Wrong number of soil layers',procName)
      END DO
    END DO

    CLOSE(inUnitV)
  ENDSUBROUTINE read_vg_params

END MODULE read_input_mod


