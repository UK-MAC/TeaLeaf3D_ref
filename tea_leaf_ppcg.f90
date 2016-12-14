!Crown Copyright 2016 AWE.
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
!>  @author Michael Boulton, Wayne Gaudin, Douglas Shanks
!>  @details Implicitly calculates the change in temperature using polynomially preconditioned CG

MODULE tea_leaf_kernel_ppcg_module

USE tea_leaf_kernel_common_module
USE tea_leaf_kernel_cheby_module
USE definitions_module, only: tl_ppcg_active ! added for ppcg init

IMPLICIT NONE

CONTAINS
 
SUBROUTINE tea_leaf_kernel_ppcg_init_sd(x_min,              &
                                        x_max,              &
                                        y_min,              &
                                        y_max,              &
                                        z_min,              &
                                        z_max,              &
                                        halo_exchange_depth,&
                                        r,                  &
                                        kx,                 &
                                        ky,                 &
                                        kz,                 &
                                        sd,                 &
                                        z,                  &
                                        cp,                 &
                                        bfp,                &
                                        Mi,                 &
                                        rx, ry, rz,         &
                                        theta,              &
                                        rtemp,              &
                                        utemp,              &
                                        preconditioner_type)

  IMPLICIT NONE

  INTEGER :: preconditioner_type
  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,z_min,z_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,&
                          y_min-halo_exchange_depth:y_max+halo_exchange_depth,&
                          z_min-halo_exchange_depth:z_max+halo_exchange_depth) :: r, sd, kx, ky, z, mi, Kz, rtemp, utemp
  REAL(KIND=8), DIMENSION(x_min:x_max,y_min:y_max,z_min:z_max) :: cp, bfp
  REAL(KIND=8) :: theta, rx, ry, rz

  INTEGER :: j,k,l

!$OMP PARALLEL

  IF (preconditioner_type .NE. TL_PREC_NONE) THEN
    ! Use a block Jacobi preconditioner
    IF (preconditioner_type .EQ. TL_PREC_JAC_BLOCK) THEN
      CALL tea_block_solve(x_min, x_max, y_min, y_max, z_min, z_max, halo_exchange_depth, r, z, cp,&
                           bfp, Kx, Ky, Kz, rx, ry, rz)
    ! Use a point diagonal Jacobi preconditioner                       
    ELSE IF (preconditioner_type .EQ. TL_PREC_JAC_DIAG) THEN
      CALL tea_diag_solve(x_min, x_max, y_min, y_max, z_min, z_max, halo_exchange_depth, r, z, Mi)
    ENDIF

!$OMP DO
  DO l=z_min,z_max
    DO k=y_min,y_max
      DO j=x_min,x_max
        sd(j, k, l) = z(j, k, l)/theta
         rtemp(j, k, l) = r(j, k, l)
        utemp(j, k, l) = sd(j, k, l)
      ENDDO
    ENDDO
  ENDDO
!$OMP END DO NOWAIT
  ELSE
!$OMP DO
  DO l=z_min,z_max
    DO k=y_min,y_max
      DO j=x_min,x_max
        sd(j, k, l) = r(j, k, l)/theta
        rtemp(j, k, l) = r(j, k, l)
        utemp(j, k, l) = sd(j, k, l)
      ENDDO
    ENDDO
  ENDDO
!$OMP END DO NOWAIT
  ENDIF
!$OMP END PARALLEL

END SUBROUTINE tea_leaf_kernel_ppcg_init_sd

SUBROUTINE tea_leaf_kernel_ppcg_inner(x_min,              &
                                      x_max,              &
                                      y_min,              &
                                      y_max,              &
                                      z_min,              &
                                      z_max,              &
                                      halo_exchange_depth,&
                                      x_min_bound,        &
                                      x_max_bound,        &
                                      y_min_bound,        &
                                      y_max_bound,        &
                                      z_min_bound,        &
                                      z_max_bound,        &
                                      alpha,              &
                                      beta,               &
                                      rx, ry, rz,         &
                                      inner_step,         &  
                                      u,                  &
                                      r,                  &
                                      rtemp,              &
                                      utemp,              &
                                      Kx,                 &
                                      Ky,                 &
                                      Kz,                 &
                                      sd,                 &
                                      z,                  &
                                      cp,                 &
                                      bfp,                &
                                      Mi,                 &
                                      preconditioner_type)
  IMPLICIT NONE

  INTEGER :: preconditioner_type
  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,z_min,z_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,&
                          y_min-halo_exchange_depth:y_max+halo_exchange_depth,&
                          z_min-halo_exchange_depth:z_max+halo_exchange_depth) :: u, r, Kx, Ky, sd, Kz, z, Mi, rtemp, utemp
  REAL(KIND=8), DIMENSION(x_min:x_max,y_min:y_max,z_min:z_max) :: cp, bfp
  INTEGER(KIND=4) :: j,k,l
  REAL(KIND=8), DIMENSION(:) :: alpha, beta
  REAL(KIND=8) :: smvp, rx, ry, rz

  INTEGER(KIND=4) :: x_min_bound, x_max_bound, y_min_bound, y_max_bound, z_min_bound, z_max_bound, inner_step

!$OMP PARALLEL PRIVATE(smvp)
!$OMP DO
  DO l=z_min_bound,z_max_bound
    DO k=y_min_bound,y_max_bound
        DO j=x_min_bound,x_max_bound
            smvp = (1.0_8                                                     &
                + rx*(Kx(j+1, k, l) + Kx(j, k, l))                            &
                + ry*(Ky(j, k+1, l) + Ky(j, k, l))                            &
                + rz*(Kz(j, k, l+1) + Kz(j, k, l)))*sd(j, k, l)               &
                - rx*(Kx(j+1, k, l)*sd(j+1, k, l) + Kx(j, k, l)*sd(j-1, k, l))&
                - ry*(Ky(j, k+1, l)*sd(j, k+1, l) + Ky(j, k, l)*sd(j, k-1, l))&
                - rz*(Kz(j, k, l+1)*sd(j, k, l+1) + Kz(j, k, l)*sd(j, k, l-1))

            rtemp(j, k, l) = rtemp(j, k, l) - smvp

        ENDDO
    ENDDO
  ENDDO
!$OMP END DO

    IF (preconditioner_type .NE. TL_PREC_NONE) THEN
      ! Use a block Jacobi preconditioner
      IF (preconditioner_type .EQ. TL_PREC_JAC_BLOCK) THEN
        CALL tea_block_solve(x_min, x_max, y_min, y_max, z_min, z_max, halo_exchange_depth, rtemp, z, cp,&
                             bfp, Kx, Ky, Kz, rx, ry, rz)
      ! Use a point diagonal Jacobi preconditioner                       
      ELSE IF (preconditioner_type .EQ. TL_PREC_JAC_DIAG) THEN
        CALL tea_diag_solve(x_min, x_max, y_min, y_max, z_min, z_max, halo_exchange_depth, rtemp, z, Mi)
      ENDIF

!$OMP DO
  DO l=z_min_bound,z_max_bound
    DO k=y_min_bound,y_max_bound
        DO j=x_min_bound,x_max_bound
        sd(j, k, l) = alpha(inner_step)*sd(j, k, l) + beta(inner_step)*z(j, k, l)
        utemp(j, k, l) = utemp(j, k, l) + sd(j, k, l)
      ENDDO
    ENDDO
  ENDDO
!$OMP END DO NOWAIT
  ELSE
!$OMP DO
  DO l=z_min_bound,z_max_bound
    DO k=y_min_bound,y_max_bound
        DO j=x_min_bound,x_max_bound
            sd(j, k, l) = alpha(inner_step)*sd(j, k, l) + beta(inner_step)*rtemp(j, k, l)
            utemp(j, k, l) = utemp(j, k, l) + sd(j, k, l)
        ENDDO
    ENDDO
  ENDDO
!$OMP END DO NOWAIT
  ENDIF
!$OMP END PARALLEL

END SUBROUTINE

SUBROUTINE tea_leaf_ppcg_calc_zrnorm_kernel(x_min,              &
                                            x_max,              &
                                            y_min,              &
                                            y_max,              &
                                            z_min,              &
                                            z_max,              &
                                            halo_exchange_depth,&
                                            z,                  &
                                            r,                  &
                                            preconditioner_type,&
                                            norm)

  IMPLICIT NONE

  INTEGER :: preconditioner_type
  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,z_min,z_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,&
                          y_min-halo_exchange_depth:y_max+halo_exchange_depth,&
                          z_min-halo_exchange_depth:z_max+halo_exchange_depth) :: r, z
  REAL(KIND=8) :: norm
  integer :: j, k,l

  norm = 0.0_8

!$OMP PARALLEL
  IF (preconditioner_type .NE. TL_PREC_NONE .OR. tl_ppcg_active) THEN
!$OMP DO REDUCTION(+:norm)
  DO l=z_min,z_max
    DO k=y_min,y_max
        DO j=x_min,x_max
            norm = norm + r(j, k, l)*z(j, k, l)
        ENDDO
    ENDDO
  ENDDO
!$OMP END DO

  ELSE
!$OMP DO REDUCTION(+:norm)
  DO l=z_min,z_max
    DO k=y_min,y_max
        DO j=x_min,x_max
            norm = norm + r(j, k, l)*r(j, k, l)
        ENDDO
    ENDDO
  ENDDO
!$OMP END DO
  ENDIF
!$OMP END PARALLEL
end SUBROUTINE tea_leaf_ppcg_calc_zrnorm_kernel

SUBROUTINE tea_calc_ls_coefs(ch_alphas, ch_betas, eigmin, eigmax, &
                             theta, ppcg_inner_steps)

  INTEGER :: ppcg_inner_steps
  REAL(KIND=8), DIMENSION(ppcg_inner_steps) :: ch_alphas, ch_betas
  REAL(KIND=8) :: eigmin, eigmax, theta

  ! TODO
  CALL tea_calc_ch_coefs(ch_alphas, ch_betas, eigmin, eigmax, &
                         theta, ppcg_inner_steps)

END SUBROUTINE

! NEW ROUTINES

! Calculate the residual, using FCG(1) to minimise rounding error

SUBROUTINE tea_leaf_kernel_ppcg_calc_rrn(x_min,              &
                                         x_max,              &
                                         y_min,              &
                                         y_max,              &
                                         z_min,              &
                                         z_max,              &
                                         halo_exchange_depth,&
                                         r,                  &
                                         rstore,             &
                                         z,                  &
                                         rrn)

  IMPLICIT NONE

  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,z_min,z_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,&
                          y_min-halo_exchange_depth:y_max+halo_exchange_depth,&
                          z_min-halo_exchange_depth:z_max+halo_exchange_depth) :: r, rstore, z

  INTEGER :: j,k,l
  REAL(KIND=8) ::  rrn

  rrn = 0.0_8

!$OMP PARALLEL
!$OMP DO REDUCTION(+:rrn)
  DO l=z_min,z_max
   DO k=y_min,y_max
    DO j=x_min,x_max
        rrn = rrn + ( r(j, k, l) - rstore(j, k, l) )*z(j, k, l)      
    ENDDO
   ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE

! Store the previous residual

SUBROUTINE tea_leaf_kernel_ppcg_store_r(x_min,              &
                                        x_max,              &
                                        y_min,              &
                                        y_max,              &
                                        z_min,              &
                                        z_max,              &
                                        halo_exchange_depth,&
                                        r,                  &
                                        rstore)

  IMPLICIT NONE

  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,z_min,z_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,&
                          y_min-halo_exchange_depth:y_max+halo_exchange_depth,&
                          z_min-halo_exchange_depth:z_max+halo_exchange_depth) :: r, rstore

  INTEGER :: j,k,l

!$OMP PARALLEL
!$OMP DO 
  DO l=z_min,z_max
   DO k=y_min,y_max
    DO j=x_min,x_max
        rstore(j, k, l) = r(j, k, l)
    ENDDO
   ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE
! Routine to update z

SUBROUTINE tea_leaf_kernel_ppcg_update_z(x_min,             &
                                        x_max,              &
                                        y_min,              &
                                        y_max,              &
                                        z_min,              &
                                        z_max,              &
                                        halo_exchange_depth,&
                                        z,                  &
                                        utemp)

  IMPLICIT NONE

  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,z_min,z_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,&
                          y_min-halo_exchange_depth:y_max+halo_exchange_depth,&
                          z_min-halo_exchange_depth:z_max+halo_exchange_depth) :: z, utemp

  INTEGER :: j,k,l

!$OMP PARALLEL
!$OMP DO 
  DO l=z_min,z_max
   DO k=y_min,y_max
    DO j=x_min,x_max
        z(j, k, l) = utemp(j, k, l)
    ENDDO
   ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE
!PPCG initialisation routine

SUBROUTINE tea_leaf_kernel_ppcg_init(x_min,              &
                                     x_max,              &
                                     y_min,              &
                                     y_max,              &
                                     z_min,              &
                                     z_max,              &
                                     halo_exchange_depth,&
                                     p,                  &
                                     r,                  &
                                     z,                  &
                                     rx, ry, rz,         &
                                     Kx,                 &
                                     Ky,                 &
                                     Kz,                 &
                                     Mi,                 &
                                     cp,                 &
                                     bfp,                &
                                     rro,                &
                                     preconditioner_type,&
                                     step)

  IMPLICIT NONE
  
  INTEGER :: preconditioner_type
  INTEGER :: step
  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,z_min,z_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,&
                          y_min-halo_exchange_depth:y_max+halo_exchange_depth,&
                          z_min-halo_exchange_depth:z_max+halo_exchange_depth) :: p, r, z, Mi, Kx, Ky, Kz
  REAL(KIND=8), DIMENSION(x_min:x_max,y_min:y_max,z_min:z_max) :: cp, bfp
  INTEGER :: j,k,l
  REAL(KIND=8) ::  rro,rx, ry, rz 
  
!$OMP PARALLEL
  IF (step == 1 .OR. step == 3) rro = 0.0_8

  IF (step == 1) THEN

!$OMP DO
  DO l=z_min,z_max
    DO k=y_min,y_max
      DO j=x_min,x_max
           p(j, k, l) = 0.0_8
           z(j, k, l) = 0.0_8
      ENDDO
    ENDDO
  ENDDO    
!$OMP END DO
  
  ELSEIF (step == 3) THEN
   
!$OMP DO
  DO l=z_min,z_max
    DO k=y_min,y_max
      DO j=x_min,x_max
           p(j, k, l) = 0.0_8
      ENDDO
    ENDDO
  ENDDO
!$OMP END DO
  
  ENDIF 
  
  IF (preconditioner_type .NE. TL_PREC_NONE .OR. (tl_ppcg_active .AND. step == 3)) THEN
  
    ! We don't apply the preconditioner on the final application of the polynomial acceleration   
    IF (step == 1 .OR. step == 2) THEN
      ! Use a block Jacobi preconditioner
      IF (preconditioner_type .EQ. TL_PREC_JAC_BLOCK) THEN
        CALL tea_block_solve(x_min, x_max, y_min, y_max, z_min, z_max, halo_exchange_depth, r, z, cp,&
                             bfp, Kx, Ky, Kz, rx, ry, rz)
      ! Use a point diagonal Jacobi preconditioner                       
      ELSE IF (preconditioner_type .EQ. TL_PREC_JAC_DIAG) THEN
        CALL tea_diag_solve(x_min, x_max, y_min, y_max, z_min, z_max, halo_exchange_depth, r, z, Mi)
      ENDIF     
    ENDIF
    
    IF ( step == 1 .OR. step ==3 ) THEN  
!$OMP DO
    DO l=z_min,z_max
      DO k=y_min,y_max
          DO j=x_min,x_max
              p(j, k, l) = z(j, k, l)
          ENDDO
      ENDDO
    ENDDO      
!$OMP END DO NOWAIT
    ENDIF    
  ELSE 
    IF (step == 1) THEN  
!$OMP DO
    DO l=z_min,z_max
      DO k=y_min,y_max
          DO j=x_min,x_max
               p(j, k, l) = r(j, k, l)
          ENDDO
      ENDDO
    ENDDO      
!$OMP END DO NOWAIT
    ENDIF
  ENDIF  
  IF (step == 1 .OR. step == 3) THEN
!$OMP DO REDUCTION(+:rro)
    DO l=z_min,z_max
      DO k=y_min,y_max
        DO j=x_min,x_max
          rro = rro + r(j, k, l)*p(j, k, l)
        ENDDO
      ENDDO
    ENDDO      
!$OMP END DO NOWAIT
  ENDIF 
!$OMP END PARALLEL 

END SUBROUTINE

END MODULE


