!###############################################################################
!#
!# Subroutines for sorting out the soil data
!#
!# ELR 18/02/2014
!#
!###############################################################################

MODULE soil_mod

  IMPLICIT none

  PUBLIC get_unique, muUnique, muUnIndx, nMuUnique, paramIndx, get_params,    &
         get_param_indices, get_soil_texture, get_vg , XWILT_OUT , XCRIT_OUT ,&
         XSATN_OUT , ONE_OVER_NMINUSONE_OUT , KSAT_OUT , ONE_OVER_ALPHA_OUT,  &
         SATHH_OUT, B_OUT, HCON_OUT, HCAP_OUT, MU_GLOBAL_UN, SU_SYM90_UN,     &
         T_SOIL_TEX, S_SOIL_TEX, get_therm

  SAVE

  INTEGER, ALLOCATABLE :: muUnique(:)
  INTEGER, ALLOCATABLE :: muUnIndx(:,:)
  INTEGER, ALLOCATABLE :: paramIndx(:)

  INTEGER :: missingVal = -9999
  !REAL :: missingValR = -9999.0
  
  INTEGER :: nMuUnique

  INTEGER, ALLOCATABLE :: MU_GLOBAL_UN(:), SEQ_UN(:), ISSOIL_UN(:)
  REAL, ALLOCATABLE :: T_SAND_UN(:), T_SILT_UN(:), T_CLAY_UN(:), T_OC_UN(:),  &
                       T_BULK_DENSITY_UN(:),                                  &
                       S_SAND_UN(:), S_SILT_UN(:), S_CLAY_UN(:), S_OC_UN(:),  &
                       S_BULK_DENSITY_UN(:)

  CHARACTER(len=10), ALLOCATABLE :: SU_SYM90_UN(:)

  ! Van Genuchten variables calculated from the input params
  REAL, ALLOCATABLE :: XWILT_OUT(:,:)               ! wilting point
  REAL, ALLOCATABLE :: XCRIT_OUT(:,:)               ! critical point
  REAL, ALLOCATABLE :: XSATN_OUT(:,:)               ! saturation point
  REAL, ALLOCATABLE :: XSATN_TMP(:)               ! saturation point on all 
                                                    ! soil types
  REAL, ALLOCATABLE :: ONE_OVER_NMINUSONE_OUT(:,:)  ! 1/(n-1) [is B in CAP]
  REAL, ALLOCATABLE :: ONE_OVER_ALPHA_OUT(:,:)      ! 1/alpha [SUCTION]
  REAL, ALLOCATABLE :: SATHH_OUT(:,:)               ! Sat soil water pressure(m)
  REAL, ALLOCATABLE :: B_OUT(:,:)                   ! B&C exponent
  REAL, ALLOCATABLE :: KSAT_OUT(:,:)                ! saturated soil cond
  REAL, ALLOCATABLE :: HCON_OUT(:,:)                ! thermal conducivity
  REAL, ALLOCATABLE :: HCAP_OUT(:,:)                ! heat capacity

  INTEGER, ALLOCATABLE :: T_SOIL_TEX(:), S_SOIL_TEX(:)

  ! Parameters
  REAL, PARAMETER :: G = 9.80665             ! van_genuchten.F
  REAL, PARAMETER :: RHO_WATER = 1000.0      ! van_genuchten.F
  REAL, PARAMETER :: ROCK_DENSITY_C = 2700.0 ! soil_control.F

  ! Silt values copied from sand (apparently temporary)
  ! Clay and sand chosen to reproduce the dry thermal conductivity and capacity 
  ! values quoted in table 4.1 of "The Frozen Earth"
  REAL, PARAMETER :: CAP_CLAY = 2.373E6   ! Volumetric heat capacity of clay,
  REAL, PARAMETER :: CAP_SAND = 2.133E6   ! sand and  
  REAL, PARAMETER :: CAP_SILT = 2.133E6   ! silt (J/m3/K).
  REAL, PARAMETER :: CON_CLAY = 1.16      ! Thermal conductivity of clay,
  REAL, PARAMETER :: CON_SAND = 1.57      ! sand and 
  REAL, PARAMETER :: CON_SILT = 1.57      ! silt minerals (W/m/K).
  REAL, PARAMETER :: CON_AIR  = 0.025     ! Heat conductivity of air (W/m/K)


  ! JOHANSEN (1975) proposes an average value for soil
  REAL, PARAMETER :: JOHANSEN_HCAP=1.942E6

  ! Threshold for equality
  REAL, PARAMETER :: mThresh=1e-6

CONTAINS

!###############################################################################
!# Get the unique variables
!###############################################################################
  SUBROUTINE get_unique

    USE read_input_mod, ONLY : mapVar, nx, ny

    USE init_mod, ONLY : nMapUnits

    IMPLICIT none

    INTEGER :: i,j

    PRINT *, "Finding the soil types required for look up table"

    ALLOCATE(muUnique(nx*ny))
    muUnique=missingVal

    ALLOCATE(muUnIndx(nx*ny,2))
    muUnIndx=missingVal


    nMuUnique=0

    DO i=1,nx
      DO j=1,ny
        IF ( .NOT. ANY(muUnique==mapVar(i,j)) )  THEN
          nMuUnique=nMuUnique+1
          muUnique(nMuUnique)=mapVar(i,j)
          muUnIndx(nMuUnique,1)=i
          muUnIndx(nMuUnique,2)=j
        END IF
      END DO
    END DO

  END SUBROUTINE get_unique

!###############################################################################
!# Get the indices needed from param file
!###############################################################################
  SUBROUTINE get_param_indices

    USE read_input_mod, ONLY : MU_GLOBAL, SEQ
    USE init_mod, ONLY : nMapUnits

    IMPLICIT none

    INTEGER :: i
    INTEGER :: indx

    LOGICAL :: found

    PRINT *, "Finding indices of the required soil types"
    ALLOCATE(paramIndx(nMuUnique))

    DO i=1,nMuUnique

      indx=1
      DO

        ! Find where in the parameter file we have the dominant soil of the
        ! required mapping unit
        IF ( MU_GLOBAL(indx)==muUnique(i) .AND. SEQ(indx)==1 ) THEN
          EXIT
        ELSE
          indx=indx+1
        END IF

        paramIndx(i)=indx

      END DO

    END DO

  END SUBROUTINE get_param_indices

!###############################################################################
!# Get the params that we need
!###############################################################################
  SUBROUTINE get_params

    USE netcdf_mod, ONLY : err_check
    USE init_mod, ONLY : missingValR
    USE read_input_mod, ONLY : MU_GLOBAL, SU_SYM90, SEQ, ISSOIL,              &
                               T_SAND, T_SILT, T_CLAY, T_OC, T_BULK_DENSITY,  &
                               S_SAND, S_SILT, S_CLAY, S_OC, S_BULK_DENSITY

    IMPLICIT none

    INTEGER :: i
    CHARACTER(len=10) :: procName='get_params'

    PRINT *, "Getting data for the required soil types"

    ALLOCATE(MU_GLOBAL_UN(nMuUnique))
    ALLOCATE(SU_SYM90_UN(nMuUnique))
    ALLOCATE(SEQ_UN(nMuUnique))
    ALLOCATE(ISSOIL_UN(nMuUnique))
    ALLOCATE(T_SAND_UN(nMuUnique))
    ALLOCATE(T_SILT_UN(nMuUnique))
    ALLOCATE(T_CLAY_UN(nMuUnique))
    ALLOCATE(T_OC_UN(nMuUnique))
    ALLOCATE(T_BULK_DENSITY_UN(nMuUnique))
    ALLOCATE(S_SAND_UN(nMuUnique))
    ALLOCATE(S_SILT_UN(nMuUnique))
    ALLOCATE(S_CLAY_UN(nMuUnique))
    ALLOCATE(S_OC_UN(nMuUnique))
    ALLOCATE(S_BULK_DENSITY_UN(nMuUnique))

    MU_GLOBAL_UN(:)=MU_GLOBAL(paramIndx(:))
    SU_SYM90_UN(:)=SU_SYM90(paramIndx(:))
    SEQ_UN(:)=SEQ(paramIndx(:))
    ISSOIL_UN(:)=ISSOIL(paramIndx(:))

    ! Get unique values
    T_CLAY_UN(:)=T_CLAY(paramIndx(:))
    T_SAND_UN(:)=T_SAND(paramIndx(:))
    T_SILT_UN(:)=T_SILT(paramIndx(:))
    T_OC_UN(:)=T_OC(paramIndx(:))
    T_BULK_DENSITY_UN(:)=T_BULK_DENSITY(paramIndx(:))
    S_CLAY_UN(:)=S_CLAY(paramIndx(:))
    S_SAND_UN(:)=S_SAND(paramIndx(:))
    S_SILT_UN(:)=S_SILT(paramIndx(:))
    S_OC_UN(:)=S_OC(paramIndx(:))
    S_BULK_DENSITY_UN(:)=S_BULK_DENSITY(paramIndx(:))

    ! Convert from percentages to fractions
    WHERE (ABS((T_CLAY_UN-missingValR)/missingValR)>=mThresh) T_CLAY_UN = T_CLAY_UN/100.0
    WHERE (ABS((T_SAND_UN-missingValR)/missingValR)>=mThresh) T_SAND_UN = T_SAND_UN/100.0
    WHERE (ABS((T_SILT_UN-missingValR)/missingValR)>=mThresh) T_SILT_UN = T_SILT_UN/100.0
    WHERE (ABS((T_OC_UN-missingValR)/missingValR)>=mThresh) T_OC_UN = T_OC_UN/100.0
    WHERE (ABS((S_CLAY_UN-missingValR)/missingValR)>=mThresh) S_CLAY_UN = S_CLAY_UN/100.0
    WHERE (ABS((S_SAND_UN-missingValR)/missingValR)>=mThresh) S_SAND_UN = S_SAND_UN/100.0
    WHERE (ABS((S_SILT_UN-missingValR)/missingValR)>=mThresh) S_SILT_UN = S_SILT_UN/100.0
    WHERE (ABS((S_OC_UN-missingValR)/missingValR)>=mThresh) S_OC_UN = S_OC_UN/100.0

    ! Convert from kg/dm3 to kg/m3
    WHERE (ABS((T_BULK_DENSITY_UN-missingValR)/missingValR)>=mThresh)                                    &
                              T_BULK_DENSITY_UN = T_BULK_DENSITY_UN*1000.0
    WHERE (ABS((S_BULK_DENSITY_UN-missingValR)/missingValR)>=mThresh)                                    &
                              S_BULK_DENSITY_UN = S_BULK_DENSITY_UN*1000.0

    ! Error checks
    ! These checks are on integers, so can use exact equality
    IF (ANY(SEQ_UN(:)/=1)) &
      CALL err_check(1,'ERROR: found some sub-dominant seq values',procName)

    IF (ANY(ISSOIL_UN(:)/=1)) &
      CALL err_check(1,'ERROR: found some non-soil points',procName)

    ! Now back to thresholds again
    !IF (ANY(ABS(T_SAND_UN(:)-missingValR)/missingValR<mThresh)) &
      !CALL err_check(1,'ERROR: missing topsoil sand',procName)

    !IF (ANY(ABS(T_SILT_UN(:)-missingValR)/missingValR<mThresh)) &
      !CALL err_check(1,'ERROR: missing topsoil silt',procName)

    !IF (ANY(ABS(T_CLAY_UN(:)-missingValR)/missingValR<mThresh)) &
      !CALL err_check(1,'ERROR: missing topsoil clay',procName)

    !IF (ANY(ABS(T_OC_UN(:)-missingValR)/missingValR<mThresh)) &
      !CALL err_check(1,'ERROR: missing topsoil organic carbon',procName)

    ! Replace missing sub-soil values with top soil
    WHERE(ABS((S_SAND_UN(:)-missingValR)/missingValR)<mThresh) S_SAND_UN = T_SAND_UN 
    WHERE(ABS((S_CLAY_UN(:)-missingValR)/missingValR)<mThresh) S_CLAY_UN = T_CLAY_UN 
    WHERE(ABS((S_SILT_UN(:)-missingValR)/missingValR)<mThresh) S_SILT_UN = T_SILT_UN 
    WHERE(ABS((S_OC_UN(:)-missingValR)/missingValR)<mThresh) S_OC_UN = T_OC_UN 
      
  END SUBROUTINE get_params

!###############################################################################
!# Soil texture
!# Heavily modified bit of code from hwsd_soil.F and AncilMod_igbp_soil.F90
!# in CAP 6.1
!###############################################################################
  SUBROUTINE get_soil_texture
    
    IMPLICIT none

    REAL, PARAMETER :: CLAY_H=0.6,                                            &
                       CLAY_M=0.35,                                           &
                       CLAY_L=0.18,                                           &
                       SILT_H=0.5,                                            &
                       SAND_H=0.65,                                           &
                       OC_H=10.0,                                             &
                       OC_M=3.0

    INTEGER, PARAMETER :: S_CR=1,                                             &
                          S_MD=2,                                             &
                          S_MF=3,                                             &
                          S_FI=4,                                             &
                          S_VF=5,                                             &
                          S_OR=6


    PRINT *, "Finding soil texture classes"

    ALLOCATE(T_SOIL_TEX(nMuUnique))
    T_SOIL_TEX(:)=missingVal
    ALLOCATE(S_SOIL_TEX(nMuUnique))
    S_SOIL_TEX(:)=missingVal

    
    CALL soil_texture(T_CLAY_UN(:) , T_SAND_UN(:) , T_SILT_UN(:) ,            &
                      T_OC_UN(:) , T_SOIL_TEX(:))
    
    CALL soil_texture(S_CLAY_UN(:) , S_SAND_UN(:) , S_SILT_UN(:) ,            &
                      S_OC_UN(:) , S_SOIL_TEX(:))

  END SUBROUTINE get_soil_texture

!###############################################################################
!# Soil texture function
!# Heavily modified bit of code from hwsd_soil.F and AncilMod_igbp_soil.F90
!# in CAP 6.1
!###############################################################################
  SUBROUTINE soil_texture(F_CLAY,F_SAND,F_SILT,F_OC,SOIL_TEX)
    
    USE netcdf_mod, ONLY : err_check

    IMPLICIT none

    REAL, INTENT(in)  :: F_CLAY(:), F_SAND(:), F_SILT(:), F_OC(:)
    INTEGER, INTENT(out) :: SOIL_TEX(:)
    CHARACTER(len=12) :: procName='soil_texture'

    REAL, PARAMETER :: CLAY_H=0.6,                                            &
                       CLAY_M=0.35,                                           &
                       CLAY_L=0.18,                                           &
                       SILT_H=0.5,                                            &
                       SAND_H=0.65,                                           &
                       OC_H=0.10,                                             &
                       OC_M=0.03
    ! In the CAP, the organic carbon fractions were actually percentages(!)
    ! Have set to fractions here for consistency

    INTEGER, PARAMETER :: S_CR=1,                                             &
                          S_MD=2,                                             &
                          S_MF=3,                                             &
                          S_FI=4,                                             &
                          S_VF=5,                                             &
                          S_OR=6

    PRINT *, "Calculating derived van Genuchten properties"

    WHERE (F_CLAY .ge. CLAY_H) 
      SOIL_TEX=S_VF ! very fine soil
    ELSEWHERE (F_CLAY .ge. CLAY_M )  
      SOIL_TEX=S_FI ! fine soil
    END WHERE

    WHERE ((F_CLAY.lt.CLAY_M).and.(F_CLAY.ge.0.0) .and.           & 
          ((F_SAND.lt.SAND_H) .or.(F_CLAY.gt.CLAY_L)))            &
                                               SOIL_TEX=S_MD ! medium soil

    WHERE ((F_SILT.ge.SILT_H).and.(F_CLAY.lt.CLAY_M))             &
                                               SOIL_TEX=S_MF ! medium-fine soil

    WHERE ((F_SAND.ge.SAND_H).and.(F_CLAY.le.CLAY_L))                         &
                                               SOIL_TEX=S_CR ! coarse soil

    WHERE (F_OC .ge. OC_H)                     SOIL_TEX=S_OR ! organic soil

    WHERE ((F_OC.ge.OC_M).and.(SOIL_TEX.eq.1)) SOIL_TEX=S_MD ! medium soil

    IF (ANY(SOIL_TEX==missingVal)) THEN
      CALL err_check(1,'ERROR: missing soil textures after allocation',procName)
    END IF
  END SUBROUTINE soil_texture

!###############################################################################
!# Convert suction (kPa) to head (m)
!###############################################################################
  REAL FUNCTION head(suction)

    IMPLICIT none

    REAL :: suction

    head = -1000.0 * suction / (G * RHO_WATER)

    RETURN

  END FUNCTION head

!###############################################################################
!# Get the Brooks&Cory/Cosby params for each soil type
!###############################################################################
  SUBROUTINE get_bc

    USE init_mod, ONLY : nSoilInLay, psi_w, psi_c, l_loge, missingValR

    IMPLICIT none

    INTEGER :: i

    ! Allocate variables
    ALLOCATE(XWILT_OUT(nMuUnique,nSoilInLay))
    ALLOCATE(XCRIT_OUT(nMuUnique,nSoilInLay))
    ALLOCATE(XSATN_OUT(nMuUnique,nSoilInLay))
    ALLOCATE(SATHH_OUT(nMuUnique,nSoilInLay))
    ALLOCATE(KSAT_OUT(nMuUnique,nSoilInLay))
    ALLOCATE(B_OUT(nMuUnique,nSoilInLay))

    ! First get topsoil
    CALL COSBY( nMuUnique, T_CLAY_UN(:) , T_SILT_UN(:) , T_SAND_UN(:) ,       &
                XWILT_OUT(:,1) , XCRIT_OUT(:,1) , XSATN_OUT(:,1) ,            &
                KSAT_OUT(:,1) , B_OUT(:,1) , SATHH_OUT(:,1) , l_loge ,        &
                PSI_W , PSI_C , missingValR )

    ! Then subsoil if necessary
    IF (nSoilInLay==2) THEN
      CALL COSBY( nMuUnique, S_CLAY_UN(:) , S_SILT_UN(:) , S_SAND_UN(:) ,     &
                  XWILT_OUT(:,2) , XCRIT_OUT(:,2) , XSATN_OUT(:,2) ,          &
                  KSAT_OUT(:,2) , B_OUT(:,2) , SATHH_OUT(:,2) , l_loge ,      &
                  PSI_W , PSI_C , missingValR )
    END IF


  END SUBROUTINE get_bc

!###############################################################################
! Cosby (aka Brooks and Corey) parameters
!###############################################################################
  SUBROUTINE COSBY( NPNTS , F_CLAY , F_SILT , F_SAND , V_WILT , V_CRIT ,      &
                    V_SAT , SATCON , B , SATHH , LLOG_E , PSI_W , PSI_C ,     &
                    missing )

    USE netcdf_mod, ONLY : err_check

    IMPLICIT none

    INTEGER, INTENT(in) :: NPNTS    ! IN Number of gridpoints.

    LOGICAL, INTENT(in) :: LLOG_E   ! IN T to use log base e, false to use 
                                    ! log base 10 in Cosby equations

    REAL, INTENT(in) :: F_CLAY(:)   ! IN Fraction of clay.
    REAL, INTENT(in) :: F_SILT(:)   ! IN Fraction of silt.
    REAL, INTENT(in) :: F_SAND(:)   ! IN Fraction of sand.
    REAL, INTENT(in) :: PSI_C       ! matrix potential for crit pt
    REAL, INTENT(in) :: PSI_W       ! matrix potential for wilt pt
    REAL, INTENT(in) :: missing     ! missing value

    REAL, INTENT(out) :: B(:)       ! OUT Clapp-Hornberger exponent.
    REAL, INTENT(out) :: SATCON(:)  ! OUT Saturated hydraulic conductivity
                                    !     (kg/m2/s).
    REAL, INTENT(out) :: SATHH(:)   ! OUT Saturated soil water pressure (m).
    REAL, INTENT(out) :: V_CRIT(:)  ! OUT Volumetric soil moisture
                                    !     concentration below which
                                    !     transpiration is water limited
                                    !     (m3 H2O/m3 soil).
    REAL, INTENT(out) :: V_SAT(:)   ! OUT Volumetric soil moisture
                                    !     concentration at saturation
                                    !     (m3 H2O/m3 soil).
    REAL, INTENT(out) :: V_WILT(:)  ! OUT Volumetric soil moisture
                                    !     concentration below which
                                    !     transpiration ceases
                                    !     (m3 H2O/m3 soil).

    ! Local
    REAL :: HEAD_C, HEAD_W

    ! Loop counter
    INTEGER :: I
  
    ! Initialize
    V_WILT(:)=missing
    V_CRIT(:)=missing
    V_SAT(:)=missing
    B(:)=missing
    SATCON(:)=missing
    SATHH(:)=missing

! Convert matrix water potentials to heads
    HEAD_C=head(PSI_C)
    HEAD_W=head(PSI_W)

    ! Loop over soils
    DO I=1,NPNTS
!-----------------------------------------------------------------------
! Calculate the primary hydraulic parameters from the multiple
! regression relationships of Cosby et al. (1984)
! NB: Cosby are ambiguous in their paper as to what log base
! is used. We have always assumed base e but most other users
! assume base 10 which is probably the correct one to use.
! Therefore, an option, LLOG_E, has been included so that the
! the user can choose which to use.
!-----------------------------------------------------------------------
      IF ( ABS((F_CLAY(I)-missing)/missing)>=mThresh .AND. &
           ABS((F_SAND(I)-missing)/missing)>=mThresh .AND. &
           ABS((F_SILT(I)-missing)/missing)>=mThresh ) THEN

        B(I)=3.10+15.70*F_CLAY(I)-0.3*F_SAND(I)
        IF(LLOG_E) THEN
          SATHH(I)=0.01*EXP(2.17-0.63*F_CLAY(I)-1.58*F_SAND(I))
          SATCON(I)=EXP(-5.55-0.64*F_CLAY(I)+1.26*F_SAND(I))
        ELSE
          SATHH(I)=0.01*10.0**(2.17-0.63*F_CLAY(I)-1.58*F_SAND(I))
          SATCON(I)=10.0**(-2.75-0.64*F_CLAY(I)+1.26*F_SAND(I))
        ENDIF
        V_SAT(I)=0.505-0.037*F_CLAY(I)-0.142*F_SAND(I)

!-----------------------------------------------------------------------
! Calculate the volumetric soil moisture concentration at the
! wilting point
!-----------------------------------------------------------------------
        V_WILT(I)=V_SAT(I)*(SATHH(I)/HEAD_W)**(1.0/B(I))

!-----------------------------------------------------------------------
! Calculate the volumetric soil moisture concentration at the
! critical point
!-----------------------------------------------------------------------
        V_CRIT(I)=V_SAT(I)*(SATHH(I)/HEAD_C)**(1.0/B(I))

      END IF
    END DO

  END SUBROUTINE COSBY


!###############################################################################
!# Get the van Genuchten params for each soil class and layer (top/sub)
!###############################################################################
  SUBROUTINE get_vg

    USE init_mod, ONLY : nSoilClass, nSoilInLay, psi_w, psi_c, missingValR

    USE read_input_mod, ONLY : VAN_G_N_IN , VAN_G_ALPHA_IN , VAN_G_THETAR_IN ,&
                                VAN_G_THETAS_IN , VAN_G_KSat_IN 
  

    IMPLICIT none

    INTEGER :: i

    ! Allocate variables
    ALLOCATE(XWILT_OUT(nSoilClass,nSoilInLay))
    ALLOCATE(XCRIT_OUT(nSoilClass,nSoilInLay))
    ALLOCATE(XSATN_OUT(nSoilClass,nSoilInLay))
    ALLOCATE(ONE_OVER_NMINUSONE_OUT(nSoilClass,nSoilInLay))
    ALLOCATE(KSAT_OUT(nSoilClass,nSoilInLay))
    ALLOCATE(ONE_OVER_ALPHA_OUT(nSoilClass,nSoilInLay))

    ! First do calculations
    DO i=1,nSoilInLay
      CALL van_genuchten(nSoilClass, VAN_G_N_IN(:,i) , VAN_G_ALPHA_IN(:,i) ,  &
                         VAN_G_THETAR_IN(:,i) , VAN_G_THETAS_IN(:,i) ,        &
                         VAN_G_KSat_IN(:,i) , XWILT_OUT(:,i) ,                &
                         XCRIT_OUT(:,i) , XSATN_OUT(:,i) ,                    &
                         ONE_OVER_NMINUSONE_OUT(:,i) , KSAT_OUT(:,i) ,   &
                         ONE_OVER_ALPHA_OUT(:,i) , psi_w , psi_c, missingValR)

    END DO


  END SUBROUTINE get_vg

!###############################################################################
!# Van Genuchten calculation
!# Heavily modified version of van_genuchten.F from the CAP 6.1
!###############################################################################
  SUBROUTINE van_genuchten(N_SOIL, VAN_G_N , VAN_G_ALPHA , VAN_G_THETAR ,     &
                           VAN_G_THETAS , VAN_G_KSat , XWILT , XCRIT ,        &
                           XSATN , ONE_OVER_NMINUSONE , KSAT ,                &
                           ONE_OVER_ALPHA , PSI_W , PSI_C, missing )

    USE netcdf_mod, ONLY : err_check

    IMPLICIT none

    ! Input scalars
    INTEGER, INTENT(in) :: N_SOIL       !IN Size of arrays

    ! Input arrays
    REAL, INTENT(in) :: VAN_G_N(:)      !IN Van Genuchten N parameter
    REAL, INTENT(in) :: VAN_G_ALPHA(:)  !IN Van Genuchten alpha parameter 
    REAL, INTENT(in) :: VAN_G_THETAR(:) !IN Van Genuchten Theta_R
    REAL, INTENT(in) :: VAN_G_THETAS(:) !IN Van Genuchten Theta_S
    REAL, INTENT(in) :: VAN_G_KSat(:)   !IN Van Genuchten KSat

    ! Input scalars
    REAL, INTENT(in) :: PSI_C           !IN matrix water potential for crit pt
    REAL, INTENT(in) :: PSI_W           !IN matrix water potential for wilt pt
    REAL, INTENT(in) :: missing     ! missing value

    ! Output arrays
    REAL, INTENT(out) :: XWILT(:)               !OUT wilting point
    REAL, INTENT(out) :: XCRIT(:)               !OUT critical point
    REAL, INTENT(out) :: XSATN(:)               !OUT saturation point
    REAL, INTENT(out) :: ONE_OVER_NMINUSONE(:)  !OUT 1/(n-1) [is B in CAP]
    REAL, INTENT(out) :: KSAT(:)                !OUT saturated soil conductivity
    REAL, INTENT(out) :: ONE_OVER_ALPHA(:)      !OUT 1/alpha [is SUCTION in CAP]

    ! Work scalars
    REAL :: HEAD_C, HEAD_W
    INTEGER :: ITYPE
    REAL :: M, THETA_CRIT, THETA_WILT, THETA_SMR


    ! Convert matrix water potentials to heads
    HEAD_C=head(PSI_C)
    HEAD_W=head(PSI_W)

    ! Initialize
    XWILT(:)=missing
    XCRIT(:)=missing
    XSATN(:)=missing
    ONE_OVER_NMINUSONE(:)=missing
    KSAT(:)=missing
    ONE_OVER_ALPHA(:)=missing

    DO ITYPE=1,N_SOIL

    ! Can't have n=0 
      IF(VAN_G_N(ITYPE).GT.0.0) THEN

        M = 1. - 1./VAN_G_N(ITYPE)

! The following are actually THETA_CRIT-THETA_R and THETA_WILT-THETA_R
        THETA_CRIT = (1.+(VAN_G_ALPHA(ITYPE)*HEAD_C) **VAN_G_N(ITYPE))**(-M) * &
                     (VAN_G_THETAS(ITYPE)-VAN_G_THETAR(ITYPE)) 

        THETA_WILT = (1.+(VAN_G_ALPHA(ITYPE)*HEAD_W) **VAN_G_N(ITYPE))**(-M) * &
                     (VAN_G_THETAS(ITYPE)-VAN_G_THETAR(ITYPE))

! This is  (THETA_SATN-THETA_R)
        THETA_SMR = (VAN_G_THETAS(ITYPE)-VAN_G_THETAR(ITYPE))

! Write values to output arrays
        XWILT(ITYPE)=THETA_WILT
        XCRIT(ITYPE)=THETA_CRIT
        XSATN(ITYPE)=THETA_SMR
        ONE_OVER_NMINUSONE(ITYPE)=1.0/(VAN_G_N(ITYPE)-1.0)
        KSAT(ITYPE)=VAN_G_Ksat(ITYPE) 
        ONE_OVER_ALPHA(ITYPE)=1.0/(VAN_G_ALPHA(ITYPE))  

      END IF

    END DO


  END SUBROUTINE van_genuchten

!###############################################################################
!# Soil thermal properties
!###############################################################################
  SUBROUTINE get_therm

    USE init_mod, ONLY : nSoilInLay, hcapMethod, hconMethod, missingValR, l_vg

    IMPLICIT none

    PRINT *, "Calculating thermal properties"

    ALLOCATE(HCON_OUT(nMuUnique,nSoilInLay))
    ALLOCATE(HCAP_OUT(nMuUnique,nSoilInLay))

    ALLOCATE(XSATN_TMP(nMuUnique))

    ! Get the top soil
    IF (l_vg) THEN
      XSATN_TMP(:)=XSATN_OUT(T_SOIL_TEX(:),1)
    ELSE
      XSATN_TMP(:)=XSATN_OUT(:,1)
    END IF

    CALL THERM_COND( nMuUnique, T_CLAY_UN(:) , T_SILT_UN(:) , T_SAND_UN(:) , &
                     XSATN_TMP , T_BULK_DENSITY_UN(:) , hconMethod ,    &
                     HCON_OUT(:,1) , ROCK_DENSITY_C , missingValR )

    CALL HEAT_CAPACITY(T_CLAY_UN(:) , T_SILT_UN(:) , T_SAND_UN(:) ,         &
                       hcapMethod, XSATN_TMP , HCAP_OUT(:,1) ,         &
                       missingValR )

    ! Get the sub soil if necessary
    IF (l_vg) THEN
      XSATN_TMP(:)=XSATN_OUT(S_SOIL_TEX(:),2)
    ELSE
      XSATN_TMP(:)=XSATN_OUT(:,2)
    END IF
    IF (nSoilInLay==2) THEN
      CALL THERM_COND( nMuUnique, S_CLAY_UN(:) , S_SILT_UN(:) , S_SAND_UN(:) , &
                       XSATN_TMP , S_BULK_DENSITY_UN(:) , hconMethod ,    &
                       HCON_OUT(:,2) , ROCK_DENSITY_C , missingValR )

      CALL HEAT_CAPACITY(S_CLAY_UN(:) , S_SILT_UN(:) , S_SAND_UN(:) ,          &
                         hcapMethod, XSATN_TMP , HCAP_OUT(:,2) ,          &
                         missingValR )

    END IF

    DEALLOCATE(XSATN_TMP)

  END SUBROUTINE get_therm

!###############################################################################
!# Thermal conductivity
!# Modified version of therm_cond.F from CAP 6.1
!# Three methods are currently.
!# METHOD 1 - Farouki (1981) method, usually used in conjunction with
!#            Cosby hydrological parameters in Clapp Hornberger hydraulics
!#
!# METHOD 2 - Peters Lidard (Johnasen) method as function of porosity
!#            and soil bulk density
!#
!# METHOD 3 - Lu method as simple regression function of porosity
!###############################################################################
  SUBROUTINE THERM_COND( NPOINTS , F_CLAY , F_SILT , F_SAND , CHI_S ,         &
                         SOIL_BULK_DENSITY , METHOD , HCON , ROCK_DENSITY,    &
                         missing )

    IMPLICIT none

    INTEGER, INTENT(in) :: NPOINTS            !IN number of points
    INTEGER, INTENT(in) :: METHOD             !IN method to use

    REAL, INTENT(in) :: F_CLAY(:) !IN fraction of clay
    REAL, INTENT(in) :: F_SILT(:) !IN fraction of silt
    REAL, INTENT(in) :: F_SAND(:) !IN fraction of sand
    REAL, INTENT(in) :: CHI_S(:)  !IN volumetric soil moisture at saturation     
    REAL, INTENT(in) :: SOIL_BULK_DENSITY(:) !IN soil bulk density
    REAL, INTENT(in) :: ROCK_DENSITY  !IN rock density
    REAL, INTENT(in) :: missing       !IN missing value 

    REAL, INTENT(out) :: HCON(:) !OUT thermal conductivity

    SELECT CASE (METHOD)

      CASE(1)
        CALL FAROUKI_THERM(NPOINTS,F_CLAY,F_SILT,F_SAND,HCON,CHI_S,missing)
 
       CASE(2)
         CALL PETERS_LIDARD_THERM(NPOINTS,CHI_S,HCON,SOIL_BULK_DENSITY,        &
                                  ROCK_DENSITY,missing)

       CASE(3)
         CALL LU_THERM(NPOINTS,CHI_S,HCON,missing)

     END SELECT

  END SUBROUTINE THERM_COND

!###############################################################################
!# Heat capacity
!# Modified version of heat_capacity.F from CAP 6.1
!###############################################################################
  SUBROUTINE HEAT_CAPACITY(F_CLAY , F_SILT , F_SAND , METHOD , CHI_S , HCAP,  &
                           missing )

    IMPLICIT none

    ! Input
    INTEGER, INTENT(in) :: METHOD                 !IN method to use

    REAL, INTENT(in) :: F_CLAY(:) !IN fraction of clay
    REAL, INTENT(in) :: F_SILT(:) !IN fraction of silt
    REAL, INTENT(in) :: F_SAND(:) !IN fraction of sand
    REAL, INTENT(in) :: CHI_S(:)  !IN volumetric soil moisture at saturation
    REAL, INTENT(in) :: missing       !IN missing value 

    ! Output
    REAL, INTENT(out) :: HCAP(:)   !OUT calculated heat capacity

    ! Local
    INTEGER :: I ! do loop variable

    REAL :: AVE_SOIL_HEATCAP  ! average heat capacity of soil
    REAL :: HEATCAP_CLAY      ! heat capacity for clay   
    REAL :: HEATCAP_SILT      ! heat capacity for silt   
    REAL :: HEATCAP_SAND      ! heat capacity for sand   

    IF(METHOD.EQ.1) THEN
! Value calculated according to relative fractions of clay/silt/sand
! For now just use the default values from CAP 6.1

      HEATCAP_CLAY=CAP_CLAY
      HEATCAP_SILT=CAP_SILT
      HEATCAP_SAND=CAP_SAND

      WHERE (ABS((F_CLAY(:)-missing)/missing)>=mThresh .AND. &
             ABS((F_SAND(:)-missing)/missing)>=mThresh .AND. &
             ABS((F_SILT(:)-missing)/missing)>=mThresh .AND. &
             ABS((CHI_S(:)-missing)/missing)>=mThresh )     &
        HCAP(:) = (1.0-CHI_S(:)) * (F_CLAY(:)*HEATCAP_CLAY +                  &
                                    F_SAND(:)*HEATCAP_SAND +                  &
                                    F_SILT(:)*HEATCAP_SILT)

    ELSE

      AVE_SOIL_HEATCAP=JOHANSEN_HCAP

      WHERE (ABS((CHI_S(:)-missing)/missing)>=mThresh ) &
        HCAP(:)=(1-CHI_S(:))*AVE_SOIL_HEATCAP

    END IF



  END SUBROUTINE HEAT_CAPACITY

!###############################################################################
!# Thermal conductivity
!# Modified version of farouki_therm.F from CAP 6.1
!###############################################################################
  SUBROUTINE FAROUKI_THERM( NPOINTS , F_CLAY , F_SILT , F_SAND , THERM_COND , &
                            CHI_S , missing )

    IMPLICIT none

    INTEGER, INTENT(in) :: NPOINTS                !IN number of points

    REAL, INTENT(in) :: F_CLAY(:)        !IN fraction of clay
    REAL, INTENT(in) :: F_SILT(:)        !IN fraction of silt
    REAL, INTENT(in) :: F_SAND(:)        !IN fraction of sand
    REAL, INTENT(in) :: CHI_S(:)         !IN volumetric saturation point
    REAL, INTENT(in) :: missing     ! missing value

    REAL, INTENT(out) :: THERM_COND(:)    !OUT dry thermal conductivity

    INTEGER :: i

    ! Initialise 
    THERM_COND(:)=missing

    DO I=1,NPOINTS

      IF(CHI_S(I).GT.0 .AND. & 
         ABS((F_CLAY(I)-missing)/missing)>=mThresh .AND. &
         ABS((F_SAND(I)-missing)/missing)>=mThresh .AND. &
         ABS((F_SILT(I)-missing)/missing)>=mThresh .AND. &
         ABS((CHI_S(I)-missing)/missing)>=mThresh ) THEN

        THERM_COND(I)=(CON_AIR**CHI_S(I))                                     &
                     *(CON_CLAY**(F_CLAY(I)*(1-CHI_S(I))))                    &
                     *(CON_SAND**(F_SAND(I)*(1-CHI_S(I))))                    &
                     *(CON_SILT**(F_SILT(I)*(1-CHI_S(I))))

      ENDIF
      
    ENDDO


  END SUBROUTINE FAROUKI_THERM

!###############################################################################
!# Thermal conductivity
!# Modified version of peters_lidard_therm.F from CAP 6.1
!###############################################################################
  SUBROUTINE PETERS_LIDARD_THERM( POINTS , CHI_S , THERM_COND ,               &
                                  SOIL_BULK_DENSITY , ROCK_DENSITY , missing )

    IMPLICIT none

    ! Input variables
    INTEGER, INTENT(in) :: POINTS       !IN number of points

    REAL, INTENT(in) :: CHI_S(:) !IN volumetric concentration at saturation pnt
    REAL, INTENT(in) :: SOIL_BULK_DENSITY(:) !IN soil bulk density
    REAL, INTENT(in) :: ROCK_DENSITY  !IN rock density
    REAL, INTENT(in) :: missing     ! missing value

    ! Output variables
    REAL, INTENT(out) :: THERM_COND(:) !OUT thermal conductivity

    ! Local variables
    INTEGER :: IPOINT       ! do loop variable

    REAL :: GAMMA_DS     ! dry soil density

    ! Initialise 
    THERM_COND(:)=missing

    DO IPOINT=1,POINTS

      IF(CHI_S(IPOINT).GT.0 .AND. & 
         ABS((CHI_S(IPOINT)-missing)/missing)>=mThresh .AND. &
         ABS((SOIL_BULK_DENSITY(IPOINT)-missing)/missing)>=mThresh ) THEN


        GAMMA_DS=SOIL_BULK_DENSITY(IPOINT)

        THERM_COND(IPOINT)=(0.135*GAMMA_DS+64.7)/ (ROCK_DENSITY-0.947*GAMMA_DS)

      ENDIF


    END DO

  END SUBROUTINE 

!###############################################################################
!# Thermal conductivity
!# Modified version of lu_therm.F from CAP 6.1
!###############################################################################
  SUBROUTINE LU_THERM( POINTS , CHI_S , THERM_COND , missing )

    IMPLICIT none

    INTEGER, INTENT(in) :: POINTS       !IN number of points

    REAL, INTENT(in) :: CHI_S(POINTS) !IN volumetric concentration at satn
    REAL, INTENT(in) :: missing     ! missing value
    REAL, INTENT(out) :: THERM_COND(POINTS) !OUT thermal conductivity

    INTEGER :: IPOINT

    ! Initialise 
    THERM_COND(:)=missing

    DO IPOINT=1,POINTS

      IF(CHI_S(IPOINT).GT.0 .AND. & 
         ABS((CHI_S(IPOINT)-missing)/missing)>=mThresh ) THEN
        THERM_COND(IPOINT)=-0.56*CHI_S(IPOINT)+0.51
      ENDIF

    ENDDO

  END SUBROUTINE LU_THERM

!###############################################################################
!# Fin
!###############################################################################
END MODULE soil_mod
