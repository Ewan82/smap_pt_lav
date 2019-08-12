!###############################################################################
!# 
!# Code to create a JULES-ready LUT from a netCDF soil map file and a text 
!# param file, using HWSD data
!#
!# ELR 18/02/2014
!#
!###############################################################################

PROGRAM main

  USE init_mod, ONLY : init, l_vg

  USE read_input_mod, ONLY : read_map, read_params, read_vg_params
  
  USE clean_mod, ONLY : clean_up

  USE soil_mod, ONLY : get_unique, get_params, get_param_indices,             &
                       get_soil_texture, get_vg, get_bc, get_therm

  USE write_mod, ONLY : write_LUT

  IMPLICIT none

  PRINT *, "Creating LUT"

  CALL init

  CALL read_map

  CALL read_params 

  IF (l_vg) CALL read_vg_params 

  CALL get_unique

  CALL get_param_indices

  CALL get_params

  IF (l_vg) THEN
    CALL get_soil_texture 
    CALL get_vg
  ELSE
    CALL get_bc
  END IF

  CALL get_therm

  CALL write_LUT

  CALL clean_up

  PRINT *, "Finished"

END PROGRAM main
