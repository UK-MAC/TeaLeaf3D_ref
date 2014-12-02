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
!>  @author David Beckingsale, Wayne Gaudin
!>  @details Implicitly calculates the change in temperature using a Jacobi iteration

MODULE tea_leaf_kernel_module

CONTAINS

SUBROUTINE tea_leaf_kernel_init(x_min,             &
                           x_max,             &
                           y_min,             &
                           y_max,             &
                           z_min,             &
                           z_max,             &
                           density,           &
                           energy,            &
                           u0,                &
                           u1,                &
                           un,                &
                           Kx,                &
                           Ky,                &
                           Kz,                &
                           coef)

  IMPLICIT NONE

   INTEGER         ::            CONDUCTIVITY        = 1 &
                                ,RECIP_CONDUCTIVITY  = 2

  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,z_min,z_max
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+2) :: density, energy, u0, Kx, Ky, Kz, u1, un

  INTEGER(KIND=4) :: coef

  INTEGER(KIND=4) :: j,k,n,l


! CALC DIFFUSION COEFFICIENT
!$OMP PARALLEL
  IF(coef .EQ. RECIP_CONDUCTIVITY) THEN
!$OMP DO 
  DO l=z_min-1,z_max+2
    DO k=y_min-1,y_max+2
      DO j=x_min-1,x_max+2
         un(j, k, l)=1.0_8/density(j, k, l)
      ENDDO
    ENDDO
  ENDDO
!$OMP END DO
  ELSE IF(coef .EQ. CONDUCTIVITY) THEN
!$OMP DO
  DO l=z_min-1,z_max+2
    DO k=y_min-1,y_max+2
      DO j=x_min-1,x_max+2
         un(j, k, l)=density(j, k, l)
      ENDDO
    ENDDO
  ENDDO
!$OMP END DO
  ENDIF

!$OMP DO
  DO l=z_min,z_max+1
  DO k=y_min,y_max+1
    DO j=x_min,x_max+1
         Kx(j, k, l)=(un(j-1,k  ,l)+un(j, k, l))/(2.0_8*un(j-1,k  ,l)*un(j, k, l))
         Ky(j, k, l)=(un(j  ,k-1,l)+un(j, k, l))/(2.0_8*un(j,  k-1,l)*un(j, k, l))
         Kz(j, k, l)=(un(j, k, l-1)+un(j, k, l))/(2.0_8*un(j,  k,l-1)*un(j, k, l))
    ENDDO
  ENDDO
  ENDDO
!$OMP END DO

!$OMP DO 
  DO l=z_min-1,z_max+1
  DO k=y_min-1, y_max+1
    DO j=x_min-1, x_max+1
      u0(j, k, l) = energy(j, k, l) * density(j, k, l)
      u1(j, k, l) = u0(j, k, l)
    ENDDO
  ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL


END SUBROUTINE tea_leaf_kernel_init

SUBROUTINE tea_leaf_kernel_solve(x_min,       &
                           x_max,             &
                           y_min,             &
                           y_max,             &
                           z_min,             &
                           z_max,             &
                           rx,                &
                           ry,                &
                           rz,                &
                           Kx,                &
                           Ky,                &
                           Kz,                &
                           error,           &
                           u0,                &
                           u1,                &
                           un)


  IMPLICIT NONE

  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,z_min,z_max
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+2) :: u0, un, u1, Kx, Ky, Kz

  REAL(KIND=8) :: ry,rx, rz,error

  INTEGER(KIND=4) :: j,k,l

  error = 0.0_8

!$OMP PARALLEL
!$OMP DO
  DO l=z_min-1,z_max+1
    DO k=y_min-1, y_max+1
      DO j=x_min-1, x_max+1
        un(j, k, l) = u1(j, k, l)
      ENDDO
    ENDDO
  ENDDO
!$OMP END DO

!$OMP DO REDUCTION(+:error)
  DO l=z_min,z_max
    DO k=y_min, y_max
      DO j=x_min, x_max
        u1(j, k, l) = (u0(j, k, l) + rx*(Kx(j+1,k  ,l)*un(j+1,k  ,l) + Kx(j, k, l)*un(j-1,k  ,l)) &
                                  + ry*(Ky(j  ,k+1,l)*un(j  ,k+1,l) + Ky(j, k, l)*un(j  ,k-1,l)) &
                                  + rz*(Kz(j  ,k,l+1)*un(j  ,k,l+1) + Kz(j, k, l)*un(j  ,k,l-1))) &
                             /(1.0_8 &
                                + rx*(Kx(j, k, l)+Kx(j+1,k,l)) &
                                + ry*(Ky(j, k, l)+Ky(j,k+1,l)) &
                                + rz*(Kz(j, k, l)+Kz(j,k,l+1)))

        error = error +  ABS(u1(j, k, l)-un(j, k, l))
      ENDDO
    ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE tea_leaf_kernel_solve

! Finalise routine is used by both implementations
SUBROUTINE tea_leaf_kernel_finalise(x_min,    &
                           x_max,             &
                           y_min,             &
                           y_max,             &
                           z_min,             &
                           z_max,             &
                           energy,            &
                           density,           &
                           u)

  IMPLICIT NONE

  INTEGER(KIND=4):: x_min,x_max,y_min,y_max, z_min, z_max
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+2) :: u
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+2) :: energy
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+2) :: density

  INTEGER(KIND=4) :: j,k,l

!$OMP PARALLEL DO 
  DO l=z_min, z_max
  DO k=y_min, y_max
    DO j=x_min, x_max
      energy(j, k, l) = u(j, k, l) / density(j, k, l)
    ENDDO
  ENDDO
  ENDDO
!$OMP END PARALLEL DO

END SUBROUTINE tea_leaf_kernel_finalise

SUBROUTINE tea_leaf_calc_residual(x_min,       &
                           x_max,             &
                           y_min,             &
                           y_max,             &
                           z_min,             &
                           z_max,             &
                           u ,                &
                           u0,                &
                           r,                &
                           Kx,                &
                           Ky,                &
                           Kz,                &
                           rx, ry, rz)

  IMPLICIT NONE

  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,z_min,z_max
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+2) :: u0, u, r
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+2) :: Kx, ky, kz

  REAL(KIND=8) :: smvp, rx, ry, rz

  INTEGER(KIND=4) :: j,k,l

!$OMP PARALLEL
!$OMP DO PRIVATE(smvp)
  DO l=z_min,z_max
    DO k=y_min,y_max
        DO j=x_min,x_max
            smvp = (1.0_8                                      &
                + rx*(Kx(j+1, k, l) + Kx(j, k, l)) &
                + ry*(Ky(j, k+1, l) + Ky(j, k, l))                      &
                + rz*(Kz(j, k, l+1) + Kz(j, k, l)))*u(j, k, l)             &
                - rx*(Kx(j+1, k, l)*u(j+1, k, l) + Kx(j, k, l)*u(j-1, k, l)) &
                - ry*(Ky(j, k+1, l)*u(j, k+1, l) + Ky(j, k, l)*u(j, k-1, l))  &
                - rz*(Kz(j, k, l+1)*u(j, k, l+1) + Kz(j, k, l)*u(j, k, l-1))
            r(j, k, l) = u0(j, k, l) - smvp
        ENDDO
    ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE tea_leaf_calc_residual

END MODULE tea_leaf_kernel_module

