!###############################################################################
!#
!# Initialisation, including reading input file
!#
!# ELR 17/02/2014
!#
!###############################################################################

MODULE init_mod

  IMPLICIT none

  PUBLIC init, inMapFileName, inParamFileName, outFileName, inMapVarName,     &
         nMapUnits, hcapMethod, hconMethod, l_vg, inVGFileName, nSoilClass,   &
         nSoilInLay, psi_w, psi_c, nSoilTop, nSoilSub, nSoil, soilDepth,      &
         l_LUT, missingValR

  SAVE

  CHARACTER(len=200) ::  inMapFileName
  CHARACTER(len=50) ::  inMapVarName
  CHARACTER(len=200) ::  inParamFileName
  CHARACTER(len=200) ::  inVGFileName
  CHARACTER(len=200) ::  outFileName

  INTEGER :: nMapUnits
  INTEGER :: hcapMEthod, hconMethod
  INTEGER :: nSoilClass
  INTEGER :: nSoilInLay
  INTEGER :: nSoilTop, nSoilSub, nSoil

  LOGICAL :: l_vg
  LOGICAL :: l_loge
  LOGICAL :: l_LUT

  REAL :: psi_w, psi_c

  REAL :: missingValR

  REAL, ALLOCATABLE :: soilDepth(:)

CONTAINS

!###############################################################################

  SUBROUTINE init

  IMPLICIT none

  WRITE (*,*) "Reading from stdin - if nothing is happening I suggest redirecting a file"

  READ(*,*) inMapFileName
  READ(*,*) inMapVarName
  READ(*,*) missingValR
  READ(*,*) inParamFileName
  READ(*,*) nMapUnits
  READ(*,*) hcapMethod, hconMethod
  READ(*,*) l_vg
  IF (l_vg) THEN
    ! Read vg LUT file name
    READ(*,*) inVGFileName
    ! Skip loge flag
    READ(*,*)
    l_loge=.FALSE.
  ELSE
    ! Skip vg LUT file name
    READ(*,*)
    inVGFileName=''
    ! Read loge flag
    READ(*,*) l_loge
  END IF
  READ(*,*) nSoilClass,nSoilInLay
  READ(*,*) l_LUT
  READ(*,*) outFileName
  READ(*,*) psi_w, psi_c
  READ(*,*) nSoilTop, nSoilSub
  nSoil=nSoilTop+nSoilSub
  ALLOCATE(soilDepth(nSoil))
  READ(*,*) soilDepth

  END SUBROUTINE init

END MODULE init_mod
