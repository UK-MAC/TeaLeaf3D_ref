!cROWn Copyright 2014 AWE.
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

!>  @brief Fortran heat conduction kernel
!>  @author Michael Boulton, Wayne Gaudin
!>  @details Implicitly calculates the change in temperature using CG method

MODULE tea_leaf_kernel_cg_module

IMPLICIT NONE

CONTAINS

SUBROUTINE tea_leaf_kernel_init_cg_fortran(x_min,             &
                           x_max,             &
                           y_min,             &
                           y_max,             &
                           z_min,             &
                           z_max,             &
                           density,           &
                           energy,            &
                           u,                 &
                           p,           &
                           r,           &
                           Mi,          &
                           w,           &
                           z,           &
                           Kx,          &
                           Ky,          &
                           Kz,          &
                           rx,          &
                           ry,          &
                           rz,          &
                           rro,         &
                           coef,        &
                           preconditioner_on)

  IMPLICIT NONE

  LOGICAL :: preconditioner_on
  INTEGER(KIND=4):: x_min,x_max,y_min,y_max, z_min, z_max
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+2) :: density
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+2) :: energy
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+2) :: u
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+2) :: p
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+2) :: r
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+2) :: Mi
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+2) :: w
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+2) :: z
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+2) :: Kx, Ky, Kz

  INTEGER(KIND=4) :: coef
  INTEGER(KIND=4) :: j,k,n,l

  REAL(kind=8) :: rro
  REAL(KIND=8) ::  rx, ry, rz

   INTEGER         ::            CONDUCTIVITY        = 1 &
                                ,RECIP_CONDUCTIVITY  = 2

  rro = 0.0_8
  p = 0.0_8
  r = 0.0_8

!$OMP PARALLEL
!$OMP DO 
  DO l=z_min-2,z_max+2
    DO k=y_min-2, y_max+2
      DO j=x_min-2, x_max+2
        u(j, k, l) = energy(j, k, l)*density(j, k, l)
      ENDDO
    ENDDO
  ENDDO
!$OMP END DO

  IF(coef .EQ. RECIP_CONDUCTIVITY) THEN
!$OMP DO 
    ! use w as temp val
  DO l=z_min-1,z_max+1
    DO k=y_min-1,y_max+1
      DO j=x_min-1,x_max+1
         w(j, k, l)=1.0_8/density(j, k, l)
      ENDDO
    ENDDO
  ENDDO
!$OMP END DO
  ELSE IF(coef .EQ. CONDUCTIVITY) THEN
!$OMP DO
  DO l=z_min-1,z_max+1
    DO k=y_min-1,y_max+1
      DO j=x_min-1,x_max+1
         w(j, k, l)=density(j, k, l)
      ENDDO
    ENDDO
  ENDDO
!$OMP END DO
  ENDIF

!$OMP DO
  DO l=z_min,z_max+1
   DO k=y_min,y_max+1
     DO j=x_min,x_max+1
          Kx(j, k, l)=(w(j-1,k  , l  ) + w(j, k, l))/(2.0_8*w(j-1,k  , l  )*w(j, k, l))
          Ky(j, k, l)=(w(j  ,k-1, l  ) + w(j, k, l))/(2.0_8*w(j  ,k-1, l  )*w(j, k, l))
          Kz(j, k, l)=(w(j  ,k  , l-1) + w(j, k, l))/(2.0_8*w(j  ,k  , l-1)*w(j, k, l))
     ENDDO
   ENDDO
  ENDDO
!$OMP END DO

!$OMP DO
  DO l=z_min,z_max
    DO k=y_min,y_max
        DO j=x_min,x_max
            w(j, k, l) = (1.0_8                                      &
                + rx*(Kx(j+1, k, l) + Kx(j, k, l))  &
                + ry*(Ky(j, k+1, l) + Ky(j, k, l))                      &
                + rz*(Kz(j, k, l+1) + Kz(j, k, l)))*u(j, k, l)             &
                - rx*(Kx(j+1, k, l)*u(j+1, k, l) + Kx(j, k, l)*u(j-1, k, l))    &
                - ry*(Ky(j, k+1, l)*u(j, k+1, l) + Ky(j, k, l)*u(j, k-1, l))  &
                - rz*(Kz(j, k, l+1)*u(j, k, l+1) + Kz(j, k, l)*u(j, k, l-1))

            r(j, k, l) = u(j, k, l) - w(j, k, l)
        ENDDO
    ENDDO
  ENDDO

  IF (preconditioner_on) then
!$OMP DO REDUCTION(+:rro)
  DO l=z_min,z_max
    DO k=y_min,y_max
        DO j=x_min,x_max
            ! inverse diagonal used as preconditioner
            Mi(j, k, l) = (1.0_8                                     &
                + rx*(Kx(j+1, k, l) + Kx(j, k, l))  &
                + ry*(Ky(j, k+1, l) + Ky(j, k, l))                      &
                + rz*(Kz(j, k, l+1) + Kz(j, k, l)))
            Mi(j, k, l) = 1.0_8/Mi(j, k, l)

            z(j, k, l) = Mi(j, k, l)*r(j, k, l)
            p(j, k, l) = z(j, k, l)

            rro = rro + r(j, k, l)*p(j, k, l);
        ENDDO
    ENDDO
  ENDDO
!$OMP END DO
  ELSE
!$OMP DO REDUCTION(+:rro)
  DO l=z_min,z_max
    DO k=y_min,y_max
        DO j=x_min,x_max
            p(j, k, l) = r(j, k, l)

            rro = rro + r(j, k, l)*p(j, k, l);
        ENDDO
    ENDDO
  ENDDO
!$OMP END DO
  ENDIF
!$OMP END PARALLEL

END SUBROUTINE tea_leaf_kernel_init_cg_fortran

SUBROUTINE tea_leaf_kernel_solve_cg_fortran_calc_w(x_min,             &
                           x_max,             &
                           y_min,             &
                           y_max,             &
                           z_min,             &
                           z_max,             &
                           p,            &
                           w,     &
                           Kx,  &
                           Ky,            &
                           Kz,            &
                           rx, &
                           ry, &
                           rz, &
                           pw)

  IMPLICIT NONE

  INTEGER(KIND=4):: x_min,x_max,y_min,y_max, z_min, z_max
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+2) :: p
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+2) :: w
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+2) :: Kx, Ky, Kz

    REAL(KIND=8) ::  rx, ry, rz

    INTEGER(KIND=4) :: j,k,n, l
    REAL(kind=8) :: pw

    pw = 0.0_08

!$OMP PARALLEL
!$OMP DO REDUCTION(+:pw)
  DO l=z_min,z_max
    DO k=y_min,y_max
        DO j=x_min,x_max
            w(j, k, l) = (1.0_8                                      &
                + rx*(Kx(j+1, k, l) + Kx(j, k, l))  &
                + ry*(Ky(j, k+1, l) + Ky(j, k, l))                      &
                + rz*(Kz(j, k, l+1) + Kz(j, k, l)))*p(j, k, l)             &
                - rx*(Kx(j+1, k, l)*p(j+1, k, l) + Kx(j, k, l)*p(j-1, k, l))    &
                - ry*(Ky(j, k+1, l)*p(j, k+1, l) + Ky(j, k, l)*p(j, k-1, l))  &
                - rz*(Kz(j, k, l+1)*p(j, k, l+1) + Kz(j, k, l)*p(j, k, l-1))

            pw = pw + w(j, k, l)*p(j, k, l)
        ENDDO
    ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE tea_leaf_kernel_solve_cg_fortran_calc_w

SUBROUTINE tea_leaf_kernel_solve_cg_fortran_calc_ur(x_min,             &
                           x_max,             &
                           y_min,             &
                           y_max,             &
                           z_min,             &
                           z_max,             &
                           u,                &
                           p,            &
                           r,            &
                           Mi,                &
                           w,     &
                           z,     &
                           alpha, &
                           rrn, &
                           preconditioner_on)

  IMPLICIT NONE

  LOGICAL :: preconditioner_on
  INTEGER(KIND=4):: x_min,x_max,y_min,y_max, z_min, z_max
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+2) :: u
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+2) :: p
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+2) :: r
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+2) :: Mi
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+2) :: w
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+2) :: z

    INTEGER(KIND=4) :: j,k,n,l
    REAL(kind=8) :: alpha, rrn

    rrn = 0.0_08

!$OMP PARALLEL
  IF (preconditioner_on) THEN
!$OMP DO REDUCTION(+:rrn)
  DO l=z_min,z_max
    DO k=y_min,y_max
        DO j=x_min,x_max
            u(j, k, l) = u(j, k, l) + alpha*p(j, k, l)
            r(j, k, l) = r(j, k, l) - alpha*w(j, k, l)
            z(j, k, l) = Mi(j, k, l)*r(j, k, l)
            rrn = rrn + r(j, k, l)*z(j, k, l)
        ENDDO
    ENDDO
  ENDDO
!$OMP END DO
  ELSE
!$OMP DO REDUCTION(+:rrn)
  DO l=z_min,z_max
    DO k=y_min,y_max
        DO j=x_min,x_max
            u(j, k, l) = u(j, k, l) + alpha*p(j, k, l)
            r(j, k, l) = r(j, k, l) - alpha*w(j, k, l)
            rrn = rrn + r(j, k, l)*r(j, k, l)
        ENDDO
    ENDDO
  ENDDO
!$OMP END DO
  ENDIF
!$OMP END PARALLEL

END SUBROUTINE tea_leaf_kernel_solve_cg_fortran_calc_ur

SUBROUTINE tea_leaf_kernel_solve_cg_fortran_calc_p(x_min,             &
                           x_max,             &
                           y_min,             &
                           y_max,             &
                           z_min,             &
                           z_max,             &
                           p,            &
                           r,            &
                           z,     &
                           beta, &
                           preconditioner_on)

  IMPLICIT NONE

  LOGICAL :: preconditioner_on
  INTEGER(KIND=4):: x_min,x_max,y_min,y_max, z_min, z_max
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+2) :: p, r, z

    REAL(kind=8) :: error

    INTEGER(KIND=4) :: j,k,n,l
    REAL(kind=8) :: beta

!$OMP PARALLEL
  IF (preconditioner_on) THEN
!$OMP DO
  DO l=z_min,z_max
    DO k=y_min,y_max
        DO j=x_min,x_max
            p(j, k, l) = z(j, k, l) + beta*p(j, k, l)
        ENDDO
    ENDDO
  ENDDO
!$OMP END DO
  ELSE
!$OMP DO
  DO l=z_min,z_max
    DO k=y_min,y_max
        DO j=x_min,x_max
            p(j, k, l) = r(j, k, l) + beta*p(j, k, l)
        ENDDO
    ENDDO
  ENDDO
!$OMP END DO
  ENDIF
!$OMP END PARALLEL

END SUBROUTINE tea_leaf_kernel_solve_cg_fortran_calc_p

END MODULE

