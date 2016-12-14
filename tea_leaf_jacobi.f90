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

MODULE tea_leaf_kernel_jacobi_module

CONTAINS

SUBROUTINE tea_leaf_kernel_jacobi_solve(x_min, &
                           x_max,              &
                           y_min,              &
                           y_max,              &
                           z_min,              &
                           z_max,              &
                           halo_exchange_depth,&
                           rx,                 &
                           ry,                 &
                           rz,                 &
                           Kx,                 &
                           Ky,                 &
                           Kz,                 &
                           error,              &
                           u0,                 &
                           u1,                 &
                           un)


  IMPLICIT NONE

  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,z_min,z_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,&
                          y_min-halo_exchange_depth:y_max+halo_exchange_depth,&
                          z_min-halo_exchange_depth:z_max+halo_exchange_depth) :: u0, un, u1, Kx, Ky, Kz

  REAL(KIND=8) :: ry,rx, rz,error

  INTEGER(KIND=4) :: j,k,l

  error = 0.0_8

!$OMP PARALLEL
!$OMP DO
  DO l=z_min-1, z_max+1
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
                                /(1.0_8 + rx*(Kx(j, k, l)+Kx(j+1,k,l)) + ry*(Ky(j, k, l)+Ky(j,k+1,l)) &
                                + rz*(Kz(j, k, l)+Kz(j,k,l+1)))

        error = error +  ABS(u1(j, k, l)-un(j, k, l))
      ENDDO
    ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE tea_leaf_kernel_jacobi_solve

END MODULE tea_leaf_kernel_jacobi_module

