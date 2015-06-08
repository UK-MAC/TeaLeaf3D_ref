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

  ! These need to be kept consistent with the data module to avoid use statement
  INTEGER,private,PARAMETER :: CHUNK_LEFT   =1    &
                            ,CHUNK_RIGHT  =2    &
                            ,CHUNK_BOTTOM =3    &
                            ,CHUNK_TOP    =4    &
                            ,CHUNK_BACK   =5    &
                            ,CHUNK_FRONT  =6    &
                            ,EXTERNAL_FACE=-1

  INTEGER,private,PARAMETER :: FIELD_DENSITY    = 1         &
                            ,FIELD_ENERGY0    = 2         &
                            ,FIELD_ENERGY1    = 3         &
                            ,FIELD_U          = 4         &
                            ,FIELD_P          = 5         &
                            ,FIELD_SD         = 6         &
                            ,FIELD_R          = 7         &
                            ,NUM_FIELDS       = 7

CONTAINS

  SUBROUTINE update_halo_kernel(x_min,x_max,y_min,y_max,z_min, z_max,halo_exchange_depth, &
                        chunk_neighbours,                                           &
                        density,                                                    &
                        energy0,                                                    &
                        energy1,                                                    &
                        u,                                                          &
                        p,                                                          &
                        sd,                                                         &
                        fields,                                                     &
                        reflective_boundary,                                        &
                        depth                                                       )
  IMPLICIT NONE

  INTEGER :: x_min,x_max,y_min,y_max,z_min,z_max,halo_exchange_depth
  LOGICAL :: reflective_boundary
  INTEGER, DIMENSION(6) :: chunk_neighbours
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max,z_min-halo_exchange_depth:z_max+halo_exchange_depth) :: density,energy0,energy1, u, sd, p

  INTEGER :: fields(NUM_FIELDS),depth

!$OMP PARALLEL

  ! Update values in external halo cells based on depth and fields requested
  ! Even though half of these loops look the wrong way around, it should be noted
  ! that depth is either 1 or 2 so that it is more efficient to always thread
  ! loop along the mesh edge.

  IF(fields(FIELD_DENSITY).EQ.1) THEN
    CALL update_halo_cell(x_min, x_max, y_min, y_max,z_min, z_max,   &
      halo_exchange_depth, chunk_neighbours, density, depth)
  ENDIF

  IF(fields(FIELD_ENERGY0).EQ.1) THEN
    CALL update_halo_cell(x_min, x_max, y_min, y_max,z_min, z_max,   &
      halo_exchange_depth, chunk_neighbours, energy0, depth)
  ENDIF

  IF(fields(FIELD_ENERGY1).EQ.1) THEN
    CALL update_halo_cell(x_min, x_max, y_min, y_max,z_min, z_max,   &
      halo_exchange_depth, chunk_neighbours, energy1, depth)
  ENDIF

  IF (reflective_boundary .EQV. .TRUE.) THEN
    IF(fields(FIELD_U).EQ.1) THEN
      CALL update_halo_cell(x_min, x_max, y_min, y_max,z_min, z_max,   &
        halo_exchange_depth, chunk_neighbours, u, depth)
    ENDIF

    IF(fields(FIELD_p).EQ.1) THEN
      CALL update_halo_cell(x_min, x_max, y_min, y_max,z_min, z_max,   &
        halo_exchange_depth, chunk_neighbours, p, depth)
    ENDIF

    IF(fields(FIELD_sd).EQ.1) THEN
      CALL update_halo_cell(x_min, x_max, y_min, y_max,z_min, z_max,   &
        halo_exchange_depth, chunk_neighbours, sd, depth)
    ENDIF
  ENDIF

!$OMP END PARALLEL

END SUBROUTINE update_halo_kernel

SUBROUTINE update_halo_cell(x_min,x_max,y_min,y_max,z_min, z_max, halo_exchange_depth,    &
                        chunk_neighbours,               &
                        mesh,                           &
                        depth                           )
  IMPLICIT NONE

  INTEGER :: x_min,x_max,y_min,y_max,z_min, z_max, halo_exchange_depth
  INTEGER, DIMENSION(6) :: chunk_neighbours
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max,z_min-halo_exchange_depth:z_max+halo_exchange_depth) :: mesh

  INTEGER :: depth

  INTEGER :: j,k,l

  IF(chunk_neighbours(CHUNK_LEFT).EQ.EXTERNAL_FACE) THEN
!$OMP DO COLLAPSE(2)
    DO l=z_min-depth,z_max+depth
      DO k=y_min-depth,y_max+depth
        DO j=1,depth
          mesh(1-j,k,l)=mesh(0+j,k,l)
        ENDDO
      ENDDO
    ENDDO
!$OMP END DO NOWAIT
  ENDIF
  IF(chunk_neighbours(CHUNK_RIGHT).EQ.EXTERNAL_FACE) THEN
!$OMP DO COLLAPSE(2)
    DO l=z_min-depth,z_max+depth
      DO k=y_min-depth,y_max+depth
        DO j=1,depth
          mesh(x_max+j,k,l)=mesh(x_max+1-j,k,l)
        ENDDO
      ENDDO
    ENDDO
!$OMP END DO NOWAIT
  ENDIF

  ! Don't need barrier if depth is only 1
  IF(depth .gt. 1) then
!$OMP BARRIER
  ENDIF

  IF(chunk_neighbours(CHUNK_BOTTOM).EQ.EXTERNAL_FACE) THEN
!$OMP DO COLLAPSE(2)
    DO l=z_min-depth,z_max+depth
      DO k=1,depth
        DO j=x_min-depth,x_max+depth
          mesh(j,1-k,l)=mesh(j,0+k,l)
        ENDDO
      ENDDO
    ENDDO
!$OMP END DO NOWAIT
  ENDIF
  IF(chunk_neighbours(CHUNK_TOP).EQ.EXTERNAL_FACE) THEN
!$OMP DO COLLAPSE(2)
    DO l=z_min-depth,z_max+depth
      DO k=1,depth
        DO j=x_min-depth,x_max+depth
          mesh(j,y_max+k,l)=mesh(j,y_max+1-k,l)
        ENDDO
      ENDDO
    ENDDO
!$OMP END DO NOWAIT
  ENDIF

  ! Again, only need barrier if depth is 1
  IF(depth .gt. 1) then
!$OMP BARRIER
  ENDIF

  IF(chunk_neighbours(CHUNK_BACK).EQ.EXTERNAL_FACE) THEN
!$OMP DO COLLAPSE(2)
    DO l=1,depth
      DO k=y_min-depth,y_max+depth
        DO j=x_min-depth,x_max+depth
          mesh(j,k,1-l)=mesh(j,k,0+l)
        ENDDO
      ENDDO
    ENDDO
!$OMP END DO NOWAIT
  ENDIF
  IF(chunk_neighbours(CHUNK_FRONT).EQ.EXTERNAL_FACE) THEN
!$OMP DO COLLAPSE(2)
    DO l=1,depth
      DO k=y_min-depth,y_max+depth
        DO j=x_min-depth,x_max+depth
          mesh(j,k,z_max+l)=mesh(j,k,z_max+1-l)
        ENDDO
      ENDDO
    ENDDO
!$OMP END DO NOWAIT
  ENDIF

END SUBROUTINE update_halo_cell

END MODULE update_halo_kernel_module

