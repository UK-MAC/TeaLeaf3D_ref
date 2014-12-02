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

!>  @brief Fortran field summary kernel
!>  @author David Beckingsale, Wayne Gaudin
!>  @details The total mass, internal energy, temperature are calculated

MODULE field_summary_kernel_module

CONTAINS

SUBROUTINE field_summary_kernel(x_min,x_max,y_min,y_max,z_min,z_max,&
                                volume,                  &
                                density,                 &
                                energy0,                 &
                                u,                       &
                                vol,mass,ie,temp         )

  IMPLICIT NONE

  INTEGER      :: x_min,x_max,y_min,y_max,z_min,z_max
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+2) :: volume,u
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+2) :: density,energy0
  REAL(KIND=8) :: vol,mass,ie,temp

  INTEGER      :: j,k,l
  REAL(KIND=8) :: cell_vol,cell_mass

  vol=0.0
  mass=0.0
  ie=0.0
  temp=0.0

!$OMP PARALLEL
!$OMP DO PRIVATE(cell_vol,cell_mass,j,k) REDUCTION(+ : vol,mass,ie, temp)
  DO l=z_min,z_max
    DO k=y_min,y_max
      DO j=x_min,x_max
        cell_vol=volume(j,k,l)
        cell_mass=cell_vol*density(j,k,l)
        vol=vol+cell_vol
        mass=mass+cell_mass
        ie=ie+cell_mass*energy0(j,k,l)
        temp = temp + cell_mass*u(j,k,l)
      ENDDO
    ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE field_summary_kernel

END MODULE field_summary_kernel_module
