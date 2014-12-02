!Crown Copyright 2014 AWE.
!
! This file is part of TeaLeaf.
!
! TeaLeaf is free software: you can redistribute it and/or modify it under 
! the terms of the GNU General Public License as published by the 
! Free Software Foundation, either version 3 of the License, or (at your option) 
! any later version.
!
! TeaLeaf is distributed in the hope that it will be useful, but 
! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
! FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more 
! details.
!
! You should have received a copy of the GNU General Public License along with 
! TeaLeaf. If not, see http://www.gnu.org/licenses/.

!>  @brief Fortran kernel to update the external halo cells in a chunk.
!>  @author David Beckingsale, Wayne Gaudin
!>  @details Updates halo cells for the required fields at the required depth
!>  for any halo cells that lie on an external boundary. The location and type
!>  of data governs how this is carried out. External boundaries are always
!>  reflective.

MODULE update_halo_kernel_module

CONTAINS

  SUBROUTINE update_halo_kernel(x_min,x_max,y_min,y_max,z_min,z_max,                &
                        chunk_neighbours,                                           &
                        density,                                                   &
                        energy0,                                                    &
                        energy1,                                                    &
                        u,                                                          &
                        p,                                                          &
                        sd,                                                         &
                        fields,                                                     &
                        depth                                                       )
  IMPLICIT NONE

  INTEGER :: x_min,x_max,y_min,y_max,z_min,z_max
  INTEGER, DIMENSION(6) :: chunk_neighbours
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+2) :: density,energy0
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+2) :: energy1
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+2) :: u, p, sd

  ! These need to be kept consistent with the data module to avoid use statement
  INTEGER,      PARAMETER :: CHUNK_LEFT   =1    &
                            ,CHUNK_RIGHT  =2    &
                            ,CHUNK_BOTTOM =3    &
                            ,CHUNK_TOP    =4    &
                            ,CHUNK_BACK   =5    &
                            ,CHUNK_FRONT  =6    &
                            ,EXTERNAL_FACE=-1

  INTEGER,      PARAMETER :: FIELD_DENSITY    = 1         &
                            ,FIELD_ENERGY0    = 2         &
                            ,FIELD_ENERGY1    = 3         &
                            ,FIELD_U          = 4         &
                            ,FIELD_P          = 5         &
                            ,FIELD_SD         = 6         &
                            ,NUM_FIELDS       = 6

  INTEGER :: fields(NUM_FIELDS),depth

  INTEGER :: j,k,l

!$OMP PARALLEL PRIVATE(j,k)

  ! Update values in external halo cells based on depth and fields requested
  ! Even though half of these loops look the wrong way around, it should be noted
  ! that depth is either 1 or 2 so that it is more efficient to always thread
  ! loop along the mesh edge.
  IF(fields(FIELD_DENSITY).EQ.1) THEN
    IF(chunk_neighbours(CHUNK_BOTTOM).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO l=z_min-depth,z_max+depth
        DO j=x_min-depth,x_max+depth
          DO k=1,depth
            density(j,1-k,l)=density(j,0+k,l)
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF(chunk_neighbours(CHUNK_TOP).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO l=z_min-depth,z_max+depth
        DO j=x_min-depth,x_max+depth
          DO k=1,depth
            density(j,y_max+k,l)=density(j,y_max+1-k,l)
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF(chunk_neighbours(CHUNK_LEFT).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO l=z_min-depth,z_max+depth
        DO k=y_min-depth,y_max+depth
          DO j=1,depth
            density(1-j,k,l)=density(0+j,k,l)
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF(chunk_neighbours(CHUNK_RIGHT).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO l=z_min-depth,z_max+depth
        DO k=y_min-depth,y_max+depth
          DO j=1,depth
            density(x_max+j,k,l)=density(x_max+1-j,k,l)
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF(chunk_neighbours(CHUNK_BACK).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO k=y_min-depth,y_max+depth
        DO j=x_min-depth,x_max+depth
          DO l=1,depth
            density(j,k,1-l)=density(j,k,0+l)
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF(chunk_neighbours(CHUNK_FRONT).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO k=y_min-depth,y_max+depth
        DO j=x_min-depth,x_max+depth
          DO l=1,depth
            density(j,k,z_max+l)=density(j,k,z_max+1-l)
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
  ENDIF

  IF(fields(FIELD_ENERGY0).EQ.1) THEN
    IF(chunk_neighbours(CHUNK_BOTTOM).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO l=z_min-depth,z_max+depth
        DO j=x_min-depth,x_max+depth
          DO k=1,depth
            energy0(j,1-k,l)=energy0(j,0+k,l)
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF(chunk_neighbours(CHUNK_TOP).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO l=z_min-depth,z_max+depth
        DO j=x_min-depth,x_max+depth
          DO k=1,depth
            energy0(j,y_max+k,l)=energy0(j,y_max+1-k,l)
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF(chunk_neighbours(CHUNK_LEFT).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO l=z_min-depth,z_max+depth
        DO k=y_min-depth,y_max+depth
          DO j=1,depth
            energy0(1-j,k,l)=energy0(0+j,k,l)
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF(chunk_neighbours(CHUNK_RIGHT).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO l=z_min-depth,z_max+depth
        DO k=y_min-depth,y_max+depth
          DO j=1,depth
            energy0(x_max+j,k,l)=energy0(x_max+1-j,k,l)
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF(chunk_neighbours(CHUNK_BACK).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO k=y_min-depth,y_max+depth
        DO j=x_min-depth,x_max+depth
          DO l=1,depth
            energy0(j,k,1-l)=energy0(j,k,0+l)
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF(chunk_neighbours(CHUNK_FRONT).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO k=y_min-depth,y_max+depth
        DO j=x_min-depth,x_max+depth
          DO l=1,depth
            energy0(j,k,z_max+l)=energy0(j,k,z_max+1-l)
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
  ENDIF

  IF(fields(FIELD_ENERGY1).EQ.1) THEN
    IF(chunk_neighbours(CHUNK_BOTTOM).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO l=z_min-depth,z_max+depth
        DO j=x_min-depth,x_max+depth
          DO k=1,depth
            energy1(j,1-k,l)=energy1(j,0+k,l)
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF(chunk_neighbours(CHUNK_TOP).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO l=z_min-depth,z_max+depth
        DO j=x_min-depth,x_max+depth
          DO k=1,depth
            energy1(j,y_max+k,l)=energy1(j,y_max+1-k,l)
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF(chunk_neighbours(CHUNK_LEFT).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO l=z_min-depth,z_max+depth
        DO k=y_min-depth,y_max+depth
          DO j=1,depth
            energy1(1-j,k,l)=energy1(0+j,k,l)
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF(chunk_neighbours(CHUNK_RIGHT).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO l=z_min-depth,z_max+depth
        DO k=y_min-depth,y_max+depth
          DO j=1,depth
            energy1(x_max+j,k,l)=energy1(x_max+1-j,k,l)
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF(chunk_neighbours(CHUNK_BACK).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO k=y_min-depth,y_max+depth
        DO j=x_min-depth,x_max+depth
          DO l=1,depth
            energy1(j,k,1-l)=energy1(j,k,0+l)
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF(chunk_neighbours(CHUNK_FRONT).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO k=y_min-depth,y_max+depth
        DO j=x_min-depth,x_max+depth
          DO l=1,depth
            energy1(j,k,z_max+l)=energy1(j,k,z_max+1-l)
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
  ENDIF

  IF(fields(FIELD_p).EQ.1) THEN
    IF(chunk_neighbours(CHUNK_BOTTOM).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO l=z_min-depth,z_max+depth
        DO j=x_min-depth,x_max+depth
          DO k=1,depth
            p(j,1-k,l)=p(j,0+k,l)
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF(chunk_neighbours(CHUNK_TOP).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO l=z_min-depth,z_max+depth
        DO j=x_min-depth,x_max+depth
          DO k=1,depth
            p(j,y_max+k,l)=p(j,y_max+1-k,l)
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF(chunk_neighbours(CHUNK_LEFT).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO l=z_min-depth,z_max+depth
        DO k=y_min-depth,y_max+depth
          DO j=1,depth
            p(1-j,k,l)=p(0+j,k,l)
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF(chunk_neighbours(CHUNK_RIGHT).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO l=z_min-depth,z_max+depth
        DO k=y_min-depth,y_max+depth
          DO j=1,depth
            p(x_max+j,k,l)=p(x_max+1-j,k,l)
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF(chunk_neighbours(CHUNK_BACK).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO k=y_min-depth,y_max+depth
        DO j=x_min-depth,x_max+depth
          DO l=1,depth
            p(j,k,1-l)=p(j,k,0+l)
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF(chunk_neighbours(CHUNK_FRONT).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO k=y_min-depth,y_max+depth
        DO j=x_min-depth,x_max+depth
          DO l=1,depth
            p(j,k,z_max+l)=p(j,k,z_max+1-l)
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
  ENDIF
  IF(fields(FIELD_u).EQ.1) THEN
    IF(chunk_neighbours(CHUNK_BOTTOM).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO l=z_min-depth,z_max+depth
        DO j=x_min-depth,x_max+depth
          DO k=1,depth
            u(j,1-k,l)=u(j,0+k,l)
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF(chunk_neighbours(CHUNK_TOP).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO l=z_min-depth,z_max+depth
        DO j=x_min-depth,x_max+depth
          DO k=1,depth
            u(j,y_max+k,l)=u(j,y_max+1-k,l)
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF(chunk_neighbours(CHUNK_LEFT).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO l=z_min-depth,z_max+depth
        DO k=y_min-depth,y_max+depth
          DO j=1,depth
            u(1-j,k,l)=u(0+j,k,l)
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF(chunk_neighbours(CHUNK_RIGHT).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO l=z_min-depth,z_max+depth
        DO k=y_min-depth,y_max+depth
          DO j=1,depth
            u(x_max+j,k,l)=u(x_max+1-j,k,l)
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF(chunk_neighbours(CHUNK_BACK).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO k=y_min-depth,y_max+depth
        DO j=x_min-depth,x_max+depth
          DO l=1,depth
            u(j,k,1-l)=u(j,k,0+l)
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF(chunk_neighbours(CHUNK_FRONT).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO k=y_min-depth,y_max+depth
        DO j=x_min-depth,x_max+depth
          DO l=1,depth
            u(j,k,z_max+l)=u(j,k,z_max+1-l)
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
  ENDIF

  IF(fields(FIELD_sd).EQ.1) THEN
    IF(chunk_neighbours(CHUNK_BOTTOM).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO l=z_min-depth,z_max+depth
        DO j=x_min-depth,x_max+depth
          DO k=1,depth
            sd(j,1-k,l)=sd(j,0+k,l)
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF(chunk_neighbours(CHUNK_TOP).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO l=z_min-depth,z_max+depth
        DO j=x_min-depth,x_max+depth
          DO k=1,depth
            sd(j,y_max+k,l)=sd(j,y_max+1-k,l)
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF(chunk_neighbours(CHUNK_LEFT).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO l=z_min-depth,z_max+depth
        DO k=y_min-depth,y_max+depth
          DO j=1,depth
            sd(1-j,k,l)=sd(0+j,k,l)
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF(chunk_neighbours(CHUNK_RIGHT).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO l=z_min-depth,z_max+depth
        DO k=y_min-depth,y_max+depth
          DO j=1,depth
            sd(x_max+j,k,l)=sd(x_max+1-j,k,l)
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF(chunk_neighbours(CHUNK_BACK).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO k=y_min-depth,y_max+depth
        DO j=x_min-depth,x_max+depth
          DO l=1,depth
            sd(j,k,1-l)=sd(j,k,0+l)
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF(chunk_neighbours(CHUNK_FRONT).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO k=y_min-depth,y_max+depth
        DO j=x_min-depth,x_max+depth
          DO l=1,depth
            sd(j,k,z_max+l)=sd(j,k,z_max+1-l)
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
  ENDIF

!$OMP END PARALLEL

END SUBROUTINE update_halo_kernel

END  MODULE update_halo_kernel_module
