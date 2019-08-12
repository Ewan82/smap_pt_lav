!###############################################################################
!#
!# Write the output LUT or ancil
!#
!# ELR 19/02/2014
!#
!###############################################################################
!#
!# Added writing of netCDF ancil 
!#
!# ELR 10/03/2014
!#
!###############################################################################

MODULE write_mod

  IMPLICIT none

  PUBLIC write_LUT, write_nc_ancil , XWILT_GRID , XCRIT_GRID , XSATN_GRID ,    &
         ONE_OVER_NMINUSONE_GRID , ONE_OVER_ALPHA_GRID , SATHH_GRID , B_GRID , &
         KSAT_GRID , HCON_GRID , HCAP_GRID

  REAL, ALLOCATABLE :: XWILT_GRID(:,:,:)               ! wilting point
  REAL, ALLOCATABLE :: XCRIT_GRID(:,:,:)               ! critical point
  REAL, ALLOCATABLE :: XSATN_GRID(:,:,:)               ! saturation point
  REAL, ALLOCATABLE :: ONE_OVER_NMINUSONE_GRID(:,:,:)  ! 1/(n-1) [is B in CAP]
  REAL, ALLOCATABLE :: ONE_OVER_ALPHA_GRID(:,:,:)      ! 1/alpha [SUCTION]
  REAL, ALLOCATABLE :: SATHH_GRID(:,:,:)               ! Sat soil water press(m)
  REAL, ALLOCATABLE :: B_GRID(:,:,:)                   ! B&C exponent
  REAL, ALLOCATABLE :: KSAT_GRID(:,:,:)                ! saturated soil cond
  REAL, ALLOCATABLE :: HCON_GRID(:,:,:)                ! thermal conducivity
  REAL, ALLOCATABLE :: HCAP_GRID(:,:,:)                ! heat capacity


CONTAINS

!###############################################################################
!# Write LUT
!###############################################################################
  SUBROUTINE write_LUT

    USE init_mod, ONLY : nSoilInLay, outFileName, psi_c, psi_w, nSoilTop,      &
                         nSoilSub, nSoil, soilDepth, l_vg, hcapMethod,         &
                         hconMethod, inMapFileName
    USE read_input_mod, ONLY : VAN_G_N_IN , VAN_G_ALPHA_IN , VAN_G_THETAR_IN , &
                               VAN_G_THETAS_IN , VAN_G_KSat_IN
                               
    USE soil_mod, ONLY : nMuUnique, XWILT_OUT, XCRIT_OUT, XSATN_OUT,           &
                         ONE_OVER_NMINUSONE_OUT, KSAT_OUT, ONE_OVER_ALPHA_OUT ,&
                         MU_GLOBAL_UN, SU_SYM90_UN, T_SOIL_TEX, S_SOIL_TEX,    &
                         HCAP_OUT, HCON_OUT, SATHH_OUT, B_OUT

    IMPLICIT none

    INTEGER :: isoil, itop
    INTEGER :: outUnit
    INTEGER :: indxTop, indxSub
    INTEGER :: iout
    CHARACTER(len=10) :: fmtstr

    PRINT *, "Writing data to look up table"
    ! Topsoil always first
    indxTop=1

    ! If only one soil layer in, subsoil is same as topsoil, otherwise it's 
    ! second
    indxSub=nSoilInLay

    PRINT *, " Opening ", TRIM(outFileName), " for writing"

    ! Open file
    OPEN( outUnit, file=TRIM(outFileName), action='write' )

    ! Headers
    IF (l_vg) THEN
      WRITE( outUnit, * ) " # HWSD soil series van Genuchten parameters for JULES"
    ELSE
      WRITE( outUnit, * ) " # HWSD soil series Brooks & Corey parameters for JULES"
    END IF

    WRITE( outUnit, * ) " # Created using soil textures from "//TRIM(inMapFileName)//"." 
    IF (hcapMethod==1) THEN
      WRITE( outUnit, * ) " # hcap method 1."
    ELSE IF (hcapMethod==2) THEN
      WRITE( outUnit, * ) " # hcap method 2 (Johansen)."
    END IF

    IF (hconMethod==1) THEN
      WRITE( outUnit, * ) " # hcon method 1 (Farouki)." 
    ELSE IF (hconMethod==2) THEN
      WRITE( outUnit, * ) " # hcon method 2 (Peters Lidard/Johansen)." 
    ELSE IF (hconMethod==3) THEN
      WRITE( outUnit, * ) " # hcon method 3 (Lu)." 
    END IF


    WRITE( outUnit, * ) "################################################################################"
    WRITE( outUnit, * ) "# Format of this file:"
    WRITE( outUnit, * ) "#   Number of soil layers"
    WRITE( outUnit, * ) "#   Soil layer thicknesses (m); top layer first"
    WRITE( outUnit, * ) "#   Number of soil types"
    WRITE( outUnit, * ) "#   The entry for each soil type consists of:"
    WRITE( outUnit, * ) "#     Soil number followed by its name"
    WRITE( outUnit, * ) "#     For each layer (top layer first):"
    IF (l_vg) THEN
      WRITE( outUnit, * ) "#       1/alpha - van Genuchten parameter (m-1)"
      WRITE( outUnit, * ) "#       1/(n-1) - van Genuchten parameter"
    ELSE
      WRITE( outUnit, * ) "#       Hydraulic conductivity at saturation (m)"
      WRITE( outUnit, * ) "#       b - Cosby / Brooks & Corey exponent"
    END IF
    WRITE( outUnit, * ) "#       heat capacity of dry soil (J m-3 K-1)"
    WRITE( outUnit, * ) "#       thermal conductivity of dry soil (W m-1 K-1)"
    WRITE( outUnit, * ) "#       hydraulic conductivity at saturation (kg m-2 s-1)"
    WRITE( outUnit, '(" "AF8.2A)' ) "#       volumetric water content at critical point (defined as ", psi_c, " kPa)"
    WRITE( outUnit, '(" "AF8.2A)' ) "#       volumetric water content at saturation  (defined as 0 kPa)"
    WRITE( outUnit, * ) "#       volumetric water content at wilting point (defined as ", psi_w, " kPa)"
    IF (l_vg) THEN
      WRITE( outUnit, * ) "#    NOTE: the above water contents are defined relative to the residual water"
      WRITE( outUnit, * ) "#    content. "
    END IF

    WRITE( outUnit, * ) "################################################################################"

    ! Run specific values
    WRITE( outUnit, '(I6)' ) nSoil
    WRITE( fmtstr, '("("I2"(F8.4))")' ) nSoil
    WRITE( outUnit, fmtstr ) soilDepth
    WRITE( outUnit, '(I6)' ) nMuUnique

    ! Loop over soils
    DO isoil=1,nMuUnique
    
      ! Soil ID
      WRITE( outUnit, '(I6","A)' ) MU_GLOBAL_UN(isoil), SU_SYM90_UN(isoil) 

      ! Topsoils
      DO itop=1,nSoilTop
        IF (l_vg) THEN
          iout=T_SOIL_TEX(isoil)
          WRITE( outUnit, '(F7.3","F7.4","E10.3","5(F8.5","))' )              &
                 ONE_OVER_ALPHA_OUT(iout,indxTop),                  &
                 ONE_OVER_NMINUSONE_OUT(iout,indxTop),              &
                 HCAP_OUT(iout,indxTop),                            &
                 HCON_OUT(iout,indxTop),                            &
                 KSAT_OUT(iout,indxTop),                            &
                 XCRIT_OUT(iout,indxTop),                           &
                 XSATN_OUT(iout,indxTop),                           &
                 XWILT_OUT(iout,indxTop)
        ELSE
          iout=isoil
          WRITE( outUnit, '(F7.3","F7.4","E10.3","5(F8.5","))' )              &
                 SATHH_OUT(iout,indxTop),                          &
                 B_OUT(iout,indxTop),                               &
                 HCAP_OUT(iout,indxTop),                            &
                 HCON_OUT(iout,indxTop),                            &
                 KSAT_OUT(iout,indxTop),                            &
                 XCRIT_OUT(iout,indxTop),                           &
                 XSATN_OUT(iout,indxTop),                           &
                 XWILT_OUT(iout,indxTop)
        END IF
      END DO

      ! Subsoils
      DO itop=1,nSoilSub
        IF (l_vg) THEN
          iout=S_SOIL_TEX(isoil)
          WRITE( outUnit, '(F7.3","F7.4","E10.3","5(F8.5","))' )    &
                 ONE_OVER_ALPHA_OUT(iout,indxSub),                  &
                 ONE_OVER_NMINUSONE_OUT(iout,indxSub),              &
                 HCAP_OUT(iout,indxSub),                            &
                 HCON_OUT(iout,indxSub),                            &
                 KSAT_OUT(iout,indxSub),                            &
                 XCRIT_OUT(iout,indxSub),                           &
                 XSATN_OUT(iout,indxSub),                           &
                 XWILT_OUT(iout,indxSub)
        ELSE
          iout=isoil
          WRITE( outUnit, '(F7.3","F7.4","E10.3","5(F8.5","))' )    &
                 SATHH_OUT(iout,indxSub),                          &
                 B_OUT(iout,indxSub),                               &
                 HCAP_OUT(iout,indxSub),                            &
                 HCON_OUT(iout,indxSub),                            &
                 KSAT_OUT(iout,indxSub),                            &
                 XCRIT_OUT(iout,indxSub),                           &
                 XSATN_OUT(iout,indxSub),                           &
                 XWILT_OUT(iout,indxSub)
        END IF



      END DO

    END DO
    CLOSE(outUnit)
  END SUBROUTINE write_LUT

!###############################################################################
!# Write netCDF ancillary file
!###############################################################################
  SUBROUTINE write_nc_ancil

    USE init_mod, ONLY : nSoil, nSoilInLay, outFileName, l_vg, missingValR,  &
                         hcapMethod, hconMethod, inMapFileName


    USE read_input_mod, ONLY : nx, ny, latVar, lonVar

    USE netcdf_mod, ONLY : create_nc_file, nc_close_input, nc_write_data

    IMPLICIT none

    INTEGER :: indxTop, indxSub
    INTEGER :: outUnit
    INTEGER :: timeVarID

    INTEGER, PARAMETER :: nVar=8
    INTEGER, PARAMETER :: nDimMax=3

    INTEGER, PARAMETER :: varNdim(nVar) = 3
    INTEGER :: varNlev(nVar)
    INTEGER :: varID(nVar)
    INTEGER :: i

    ! Not compressing the grid
    LOGICAL, PARAMETER :: useT=.FALSE.
    LOGICAL, PARAMETER :: regLatLon=.FALSE.

    ! title
    CHARACTER(len=200) :: title
    CHARACTER(len=20) :: varDims(nvar,nDimMax)
    CHARACTER(len=16) :: varName(nvar)
    CHARACTER(len=42) :: varLongName(nvar)
    CHARACTER(len=10) :: varUnits(nvar)

    PRINT *, "Writing data to netCDF ancillary file"
    ! Topsoil always first
    indxTop=1

    ! If only one soil layer in, subsoil is same as topsoil, otherwise it's 
    ! second
    indxSub=nSoilInLay

    varNlev(:)=nSoil

    ! Title for the file
    IF (l_vg) THEN
      title="HWSD soil series van Genuchten parameters for JULES."
    ELSE
      title="HWSD soil series Brooks & Corey parameters for JULES."
    END IF
    title=TRIM(title)//" Created using soil textures from "//TRIM(inMapFileName)//"." 
    IF (hcapMethod==1) THEN
      title=TRIM(title)//" hcap method 1."
    ELSE IF (hcapMethod==2) THEN
      title=TRIM(title)//" hcap method 2 (Johansen)."
    END IF

    IF (hconMethod==1) THEN
      title=TRIM(title)//" hcon method 1 (Farouki)." 
    ELSE IF (hconMethod==2) THEN
      title=TRIM(title)//" hcon method 2 (Peters Lidard/Johansen)." 
    ELSE IF (hconMethod==3) THEN
      title=TRIM(title)//" hcon method 3 (Lu)." 
    END IF

    PRINT *, " Opening ", TRIM(outFileName), " for writing"

    ! Variable dimensions (all the same)
    DO i=1,nVar
      varDims(i,:)=(/'x','y','z'/)
    END DO

    ! Variable names and long names
    ! Note that for vG we still keep the B&C names

    ! Not any more - now we use the appropriate names by replacing them if
    ! we're using van Genuchten

    varName = (/'sathh           ', 'b               ', 'hcap            ', &
                'hcon            ', 'satcon          ', 'vcrit           ', &
                'vsat            ', 'vwilt           '/) 

    varLongName = (/'Hydraulic conductivity at saturation      ', &
                    'b - Cosby / Brooks & Corey exponent       ', &
                    'Heat capacity of dry soil                 ', &
                    'Thermal conductivity of dry soil          ', &
                    'Hydraulic conductivity at saturation      ', &
                    'Volumetric water content at critical point', &
                    'Volumetric water content at saturation    ', &
                    'Volumetric water content at wilting point ' /)

    IF (l_vg) THEN
      varName(1)='oneoveralpha    '
      varName(2)='oneovernminusone'
      varLongName(1)='Reciprocal of van G alpha parameter       '
      varLongName(2)='Reciprocal of van G n parameter minus one '
    END IF

    varUnits = (/ 'm         ', '-         ', 'J m-3 K-1 ', 'W m-1 K-1 ', & 
                  'kg m-2 s-1', 'kPa       ', 'kPa       ', 'kPa       ' /)
    
    CALL create_nc_file( nDimMax , 1 , nx , ny , nSoil , nVar ,          &
                         missingValR , .FALSE. , .FALSE. ,      &
                         '' , outFileName , TRIM(title),           &
                         varNdim , varNlev , latVar, lonVar,varDims, &
                         varLongName , varName , varUnits , outUnit ,        &
                         timeVarID , varID , regLatLon , useT )

    CALL put_on_grid

! Write the various variables
      IF (l_vg) THEN
        CALL nc_write_data( nx, ny, nSoil, outUnit, varID(1), ONE_OVER_ALPHA_GRID)
        CALL nc_write_data( nx, ny, nSoil, outUnit, varID(2), ONE_OVER_NMINUSONE_GRID)
      ELSE
        CALL nc_write_data( nx, ny, nSoil, outUnit, varID(1), SATHH_GRID )
        CALL nc_write_data( nx, ny, nSoil, outUnit, varID(2), B_GRID )
      END IF

      CALL nc_write_data( nx, ny, nSoil, outUnit, varID(3), HCAP_GRID)
      CALL nc_write_data( nx, ny, nSoil, outUnit, varID(4), HCON_GRID)
      CALL nc_write_data( nx, ny, nSoil, outUnit, varID(5), KSAT_GRID)
      CALL nc_write_data( nx, ny, nSoil, outUnit, varID(6), XCRIT_GRID)
      CALL nc_write_data( nx, ny, nSoil, outUnit, varID(7), XSATN_GRID)
      CALL nc_write_data( nx, ny, nSoil, outUnit, varID(8), XWILT_GRID)

    CALL nc_close_input(outUnit)

  END SUBROUTINE write_nc_ancil

!###############################################################################
!# Put the params on a full grid if necessary
!###############################################################################
  SUBROUTINE put_on_grid

    USE read_input_mod, ONLY : mapVar, nx, ny

    USE init_mod, ONLY : nSoil, nSoilInLay, nSoilTop, nSoilSub, l_vg, &
                         missingValR

    USE soil_mod, ONLY : ONE_OVER_ALPHA_OUT, ONE_OVER_NMINUSONE_OUT, &
                         SATHH_OUT, B_OUT, HCAP_OUT, HCON_OUT, KSAT_OUT, &
                         XCRIT_OUT, XSATN_OUT, XWILT_OUT, nMuUnique, &
                         muUnique, T_SOIL_TEX, S_SOIL_TEX

    IMPLICIT none

    INTEGER, ALLOCATABLE :: indxs(:,:)
    INTEGER :: isoil, ix, iy
    INTEGER :: indxTop, indxSub
    INTEGER :: iout

    INTEGER :: missingVal = -9999

    REAL, PARAMETER :: thresh=0.01

    PRINT *, "Putting on grid"

    ALLOCATE( XWILT_GRID(nx,ny,nSoil))               
    ALLOCATE( XCRIT_GRID(nx,ny,nSoil))               
    ALLOCATE( XSATN_GRID(nx,ny,nSoil))               
    IF (l_vg) THEN
      ALLOCATE( ONE_OVER_NMINUSONE_GRID(nx,ny,nSoil))  
      ALLOCATE( ONE_OVER_ALPHA_GRID(nx,ny,nSoil))      
    ELSE
      ALLOCATE( SATHH_GRID(nx,ny,nSoil))               
      ALLOCATE( B_GRID(nx,ny,nSoil))                  
    END IF
    ALLOCATE( KSAT_GRID(nx,ny,nSoil))               
    ALLOCATE( HCON_GRID(nx,ny,nSoil))              
    ALLOCATE( HCAP_GRID(nx,ny,nSoil))              
    
    ALLOCATE( indxs(nx,ny) )

    ! Initialise to missing values
    XWILT_GRID(:,:,:)=missingValR
    XCRIT_GRID(:,:,:)=missingValR
    XSATN_GRID(:,:,:)=missingValR
    IF (l_vg) THEN
      ONE_OVER_NMINUSONE_GRID(:,:,:)=missingValR
      ONE_OVER_ALPHA_GRID(:,:,:)=missingValR
    ELSE
      SATHH_GRID(:,:,:)=missingValR
      B_GRID(:,:,:)=missingValR
    END IF
    KSAT_GRID(:,:,:)=missingValR
    HCON_GRID(:,:,:)=missingValR
    HCAP_GRID(:,:,:)=missingValR
    
    indxs(:,:)=missingVal 

    ! Topsoil always first
    indxTop=1

    ! If only one soil layer in, subsoil is same as topsoil, otherwise it's 
    ! second
    indxSub=nSoilInLay

    ! Get the indices to go from the individual soil types back up to the grid
    DO isoil=1,nMuUnique
      WHERE(ABS(mapVar(:,:)-muUnique(isoil))<thresh) indxs=isoil
    END DO

    IF (l_vg) THEN
      ! van Genuchten only 
      DO ix=1,nx
        DO iy=1,ny

          IF (indxs(ix,iy)/=missingVal) THEN
            ! Topsoil 
            iout=T_SOIL_TEX(indxs(ix,iy))

            ! Select by soil class
            ONE_OVER_ALPHA_GRID(ix,iy,1:nSoilTop) = ONE_OVER_ALPHA_OUT(iout,indxTop)
            ONE_OVER_NMINUSONE_GRID(ix,iy,1:nSoilTop) = ONE_OVER_NMINUSONE_OUT(iout,indxTop)
            KSAT_GRID(ix,iy,1:nSoilTop) = KSAT_OUT(iout,indxTop) 
            XCRIT_GRID(ix,iy,1:nSoilTop) = XCRIT_OUT(iout,indxTop) 
            XSATN_GRID(ix,iy,1:nSoilTop) = XSATN_OUT(iout,indxTop) 
            XWILT_GRID(ix,iy,1:nSoilTop) = XWILT_OUT(iout,indxTop) 

            ! Select by soil type 
            HCAP_GRID(ix,iy,1:nSoilTop) = HCAP_OUT(indxs(ix,iy),indxTop) 
            HCON_GRID(ix,iy,1:nSoilTop) = HCON_OUT(indxs(ix,iy),indxTop) 

            ! subsoil
            iout=S_SOIL_TEX(indxs(ix,iy))

            ! Select by soil class
            ONE_OVER_ALPHA_GRID(ix,iy,nSoilTop+1:nSoil) = ONE_OVER_ALPHA_OUT(iout,indxSub)
            ONE_OVER_NMINUSONE_GRID(ix,iy,nSoilTop+1:nSoil) = ONE_OVER_NMINUSONE_OUT(iout,indxSub)
            KSAT_GRID(ix,iy,nSoilTop+1:nSoil) = KSAT_OUT(iout,indxSub) 
            XCRIT_GRID(ix,iy,nSoilTop+1:nSoil) = XCRIT_OUT(iout,indxSub) 
            XSATN_GRID(ix,iy,nSoilTop+1:nSoil) = XSATN_OUT(iout,indxSub) 
            XWILT_GRID(ix,iy,nSoilTop+1:nSoil) = XWILT_OUT(iout,indxSub) 

            ! Select by soil type
            HCAP_GRID(ix,iy,nSoilTop+1:nSoil) = HCAP_OUT(indxs(ix,iy),indxSub) 
            HCON_GRID(ix,iy,nSoilTop+1:nSoil) = HCON_OUT(indxs(ix,iy),indxSub) 

          END IF

        END DO
      END DO
    ELSE
      ! B&C only
      DO ix=1,nx
        DO iy=1,ny

          IF (indxs(ix,iy)/=missingVal) THEN

            iout=indxs(ix,iy)

            ! Topsoil 
            ! Select by soil type
            SATHH_GRID(ix,iy,1:nSoilTop) = SATHH_OUT(iout,indxTop)
            B_GRID(ix,iy,1:nSoilTop) = B_OUT(iout,indxTop)
            HCAP_GRID(ix,iy,1:nSoilTop) = HCAP_OUT(iout,indxTop) 
            HCON_GRID(ix,iy,1:nSoilTop) = HCON_OUT(iout,indxTop) 
            KSAT_GRID(ix,iy,1:nSoilTop) = KSAT_OUT(iout,indxTop) 
            XCRIT_GRID(ix,iy,1:nSoilTop) = XCRIT_OUT(iout,indxTop) 
            XSATN_GRID(ix,iy,1:nSoilTop) = XSATN_OUT(iout,indxTop) 
            XWILT_GRID(ix,iy,1:nSoilTop) = XWILT_OUT(iout,indxTop) 

            ! subsoil
            ! Select by soil type
            SATHH_GRID(ix,iy,nSoilTop+1:nSoil) = SATHH_OUT(iout,indxSub)
            B_GRID(ix,iy,nSoilTop+1:nSoil) = B_OUT(iout,indxSub)
            HCAP_GRID(ix,iy,nSoilTop+1:nSoil) = HCAP_OUT(iout,indxSub) 
            HCON_GRID(ix,iy,nSoilTop+1:nSoil) = HCON_OUT(iout,indxSub) 
            KSAT_GRID(ix,iy,nSoilTop+1:nSoil) = KSAT_OUT(iout,indxSub) 
            XCRIT_GRID(ix,iy,nSoilTop+1:nSoil) = XCRIT_OUT(iout,indxSub) 
            XSATN_GRID(ix,iy,nSoilTop+1:nSoil) = XSATN_OUT(iout,indxSub) 
            XWILT_GRID(ix,iy,nSoilTop+1:nSoil) = XWILT_OUT(iout,indxSub) 
          END IF
        END DO
      END DO

    END IF
    


  END SUBROUTINE put_on_grid


END MODULE write_mod
