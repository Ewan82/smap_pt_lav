!###############################################################################
!#
!# Clean up allocated arrays
!#
!# ELR 19/02/2014
!#
!###############################################################################

MODULE clean_mod

  IMPLICIT none

  PUBLIC clean_up

CONTAINS

!###############################################################################
!# Clean up
!###############################################################################
  SUBROUTINE clean_up

    USE netcdf
    USE netcdf_mod

    USE init_mod, ONLY : soilDepth, l_vg, l_lut

    USE read_input_mod, ONLY : mapVar, MU_GLOBAL, SEQ, ISSOIL, T_SAND,        &
                               T_SILT, T_CLAY, T_OC, S_SAND, S_SILT, S_CLAY,  &
                               S_OC, VAN_G_N_IN ,                       &
                               VAN_G_ALPHA_IN , VAN_G_THETAR_IN ,             &
                               VAN_G_THETAS_IN , VAN_G_KSat_IN


    USE soil_mod, ONLY : muUnique, muUnIndx, paramIndx, MU_GLOBAL_UN,         &
                         SU_SYM90_UN, SEQ_UN, &
                         ISSOIL_UN, T_SAND_UN, T_SILT_UN, T_CLAY_UN, T_OC_UN, &
                         S_SAND_UN, S_SILT_UN, S_CLAY_UN, S_OC_UN, T_SOIL_TEX,&
                         S_SOIL_TEX, XWILT_OUT, XCRIT_OUT, XSATN_OUT,         &
                         ONE_OVER_NMINUSONE_OUT, KSAT_OUT, ONE_OVER_ALPHA_OUT,&
                         B_OUT, SATHH_OUT

    USE write_mod, ONLY : XWILT_GRID, XCRIT_GRID, XSATN_GRID,&
                         ONE_OVER_NMINUSONE_GRID, KSAT_GRID,                  &
                         ONE_OVER_ALPHA_GRID, B_GRID, SATHH_GRID

!-------------------------------------------------------------------------------

    IMPLICIT none

    CHARACTER(len=8) :: procName='clean_up'

    PRINT *, "Cleaning up"
!-------------------------------------------------------------------------------
! Deallocate variables
!-------------------------------------------------------------------------------

    ! init_mod
    DEALLOCATE(soilDepth)

    ! read_input_mod
    DEALLOCATE (mapVar)

    DEALLOCATE(MU_GLOBAL) 
    DEALLOCATE(SEQ) 
    DEALLOCATE(ISSOIL)
    DEALLOCATE(T_SAND) 
    DEALLOCATE(T_SILT) 
    DEALLOCATE(T_CLAY) 
    DEALLOCATE(T_OC) 
    DEALLOCATE(S_SAND) 
    DEALLOCATE(S_SILT) 
    DEALLOCATE(S_CLAY) 
    DEALLOCATE(S_OC) 
    IF (l_vg) THEN
      DEALLOCATE(VAN_G_N_IN) 
      DEALLOCATE(VAN_G_ALPHA_IN) 
      DEALLOCATE(VAN_G_THETAR_IN) 
      DEALLOCATE(VAN_G_THETAS_IN) 
      DEALLOCATE(VAN_G_KSat_IN) 
    END IF

    ! soil_mod
    DEALLOCATE(muUnique)
    DEALLOCATE(muUnIndx)
    DEALLOCATE(paramIndx)
    DEALLOCATE(MU_GLOBAL_UN)
    DEALLOCATE(SU_SYM90_UN)
    DEALLOCATE(SEQ_UN)
    DEALLOCATE(ISSOIL_UN)
    DEALLOCATE(T_SAND_UN)
    DEALLOCATE(T_SILT_UN)
    DEALLOCATE(T_CLAY_UN)
    DEALLOCATE(T_OC_UN)
    DEALLOCATE(S_SAND_UN)
    DEALLOCATE(S_SILT_UN)
    DEALLOCATE(S_CLAY_UN)
    DEALLOCATE(S_OC_UN)
    IF (l_vg) THEN
      DEALLOCATE(T_SOIL_TEX)
      DEALLOCATE(S_SOIL_TEX)
    END IF
    DEALLOCATE(XWILT_OUT)
    DEALLOCATE(XCRIT_OUT)
    DEALLOCATE(XSATN_OUT)
    DEALLOCATE(KSAT_OUT)
    IF (l_vg) THEN
      DEALLOCATE(ONE_OVER_NMINUSONE_OUT)
      DEALLOCATE(ONE_OVER_ALPHA_OUT)
    ELSE
      DEALLOCATE(B_OUT)
      DEALLOCATE(SATHH_OUT)
    END IF

    IF ( .NOT. l_lut ) THEN 
      DEALLOCATE(XWILT_GRID)
      DEALLOCATE(XCRIT_GRID)
      DEALLOCATE(XSATN_GRID)
      DEALLOCATE(KSAT_GRID)
      IF (l_vg) THEN
        DEALLOCATE(ONE_OVER_NMINUSONE_GRID)
        DEALLOCATE(ONE_OVER_ALPHA_GRID)
      ELSE
        DEALLOCATE(B_GRID)
        DEALLOCATE(SATHH_GRID)
      END IF
    END IF




  END SUBROUTINE clean_up


END MODULE clean_mod
