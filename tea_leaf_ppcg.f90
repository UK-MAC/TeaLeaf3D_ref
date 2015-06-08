MODULE tea_leaf_kernel_ppcg_module

USE tea_leaf_kernel_common_module
USE tea_leaf_kernel_cheby_module

IMPLICIT NONE

CONTAINS

SUBROUTINE tea_leaf_kernel_ppcg_init_sd(x_min,             &
                                        x_max,             &
                                        y_min,             &
                                        y_max,             &
                                        z_min,             &
                                        z_max,             &
                                        halo_exchange_depth,             &
                                        r,                 &
                                        kx,                 &
                                        ky,                 &
                                        kz,                 &
                                        sd,                &
                                        z,                &
                                        cp,                &
                                        bfp,                &
                                        Mi,                &
                                        rx, ry, rz,         &
                                        theta,             &
                                        preconditioner_type)

  IMPLICIT NONE

  INTEGER :: preconditioner_type
  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,z_min,z_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth,z_min-halo_exchange_depth:z_max+halo_exchange_depth) :: r, sd, kx, ky, z, mi, Kz
  REAL(KIND=8), DIMENSION(x_min:x_max,y_min:y_max,z_min:z_max) :: cp, bfp
  REAL(KIND=8) :: theta, theta_r, rx, ry, rz

  INTEGER :: j,k,l

  theta_r = 1.0_8/theta

!$OMP PARALLEL

  IF (preconditioner_type .NE. TL_PREC_NONE) THEN

    IF (preconditioner_type .EQ. TL_PREC_JAC_BLOCK) THEN
      CALL tea_block_solve(x_min, x_max, y_min, y_max, z_min, z_max,    &
        halo_exchange_depth, r, z, cp, bfp, Kx, Ky, Kz, rx, ry, rz)
    ELSE IF (preconditioner_type .EQ. TL_PREC_JAC_DIAG) THEN
      CALL tea_diag_solve(x_min, x_max, y_min, y_max, z_min, z_max, &
        halo_exchange_depth, r, z, Mi, Kx, Ky, Kz, rx, ry, rz)
    ENDIF

!$OMP DO
  DO l=z_min,z_max
    DO k=y_min,y_max
      DO j=x_min,x_max
        sd(j, k, l) = z(j, k, l)*theta_r
      ENDDO
    ENDDO
  ENDDO
!$OMP END DO NOWAIT
  ELSE
!$OMP DO
  DO l=z_min,z_max
    DO k=y_min,y_max
      DO j=x_min,x_max
        sd(j, k, l) = r(j, k, l)*theta_r
      ENDDO
    ENDDO
  ENDDO
!$OMP END DO NOWAIT
  ENDIF
!$OMP END PARALLEL

END SUBROUTINE tea_leaf_kernel_ppcg_init_sd

SUBROUTINE tea_leaf_kernel_ppcg_inner(x_min,             &
                                      x_max,             &
                                      y_min,             &
                                      y_max,             &
                                      z_min,             &
                                      z_max,             &
                                      halo_exchange_depth,             &
                                      outer_step,              &
                                      alpha,             &
                                      beta,              &
                                      rx, ry, rz,         &
                                      ppcg_cur_step, tl_ppcg_inner_steps,   &
                                      u,                 &
                                      r,                 &
                                      Kx,                &
                                      Ky,                &
                                      Kz,                &
                                      sd,                &
                                      z,                &
                                      cp,                &
                                      bfp,                &
                                      Mi,                &
                                      preconditioner_type)
  IMPLICIT NONE

  INTEGER :: preconditioner_type
  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,z_min,z_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth,z_min-halo_exchange_depth:z_max+halo_exchange_depth) :: u, r, Kx, Ky, sd, kz, z, Mi
  REAL(KIND=8), DIMENSION(x_min:x_max,y_min:y_max,z_min:z_max) :: cp, bfp
  INTEGER(KIND=4) :: j,k,l, outer_step
  REAL(KIND=8), DIMENSION(:) :: alpha, beta
  REAL(KIND=8) :: smvp, rx, ry, rz

  INTEGER(KIND=4) :: bounds_extra, ppcg_cur_step, tl_ppcg_inner_steps, inner_step

  inner_step = outer_step

  DO bounds_extra = halo_exchange_depth-1, 0, -1
!$OMP PARALLEL PRIVATE(smvp)
!$OMP DO
  DO l=z_min-bounds_extra,z_max+bounds_extra
    DO k=y_min-bounds_extra,y_max+bounds_extra
        DO j=x_min-bounds_extra,x_max+bounds_extra
            smvp = (1.0_8                                      &
                + rx*(Kx(j+1, k, l) + Kx(j, k, l))  &
                + ry*(Ky(j, k+1, l) + Ky(j, k, l))                      &
                + rz*(Kz(j, k, l+1) + Kz(j, k, l)))*sd(j, k, l)             &
                - rx*(Kx(j+1, k, l)*sd(j+1, k, l) + Kx(j, k, l)*sd(j-1, k, l))    &
                - ry*(Ky(j, k+1, l)*sd(j, k+1, l) + Ky(j, k, l)*sd(j, k-1, l))  &
                - rz*(Kz(j, k, l+1)*sd(j, k, l+1) + Kz(j, k, l)*sd(j, k, l-1))

            r(j, k, l) = r(j, k, l) - smvp

            u(j, k, l) = u(j, k, l) + sd(j, k, l)
        ENDDO
    ENDDO
  ENDDO
!$OMP END DO

  IF (preconditioner_type .NE. TL_PREC_NONE) THEN

    IF (preconditioner_type .EQ. TL_PREC_JAC_BLOCK) THEN
      CALL tea_block_solve(x_min, x_max, y_min, y_max, z_min, z_max,    &
        halo_exchange_depth, r, z, cp, bfp, Kx, Ky, Kz, rx, ry, rz)
    ELSE IF (preconditioner_type .EQ. TL_PREC_JAC_DIAG) THEN
      CALL tea_diag_solve(x_min, x_max, y_min, y_max, z_min, z_max, &
        halo_exchange_depth, r, z, Mi, Kx, Ky, Kz, rx, ry, rz)
    ENDIF

!$OMP DO
  DO l=z_min-bounds_extra,z_max+bounds_extra
    DO k=y_min-bounds_extra,y_max+bounds_extra
      DO j=x_min-bounds_extra,x_max+bounds_extra
        sd(j, k, l) = alpha(inner_step)*sd(j, k, l) + beta(inner_step)*z(j, k, l)
      ENDDO
    ENDDO
  ENDDO
!$OMP END DO NOWAIT
  ELSE
!$OMP DO
  DO l=z_min-bounds_extra,z_max+bounds_extra
    DO k=y_min-bounds_extra,y_max+bounds_extra
        DO j=x_min-bounds_extra,x_max+bounds_extra
            sd(j, k, l) = alpha(inner_step)*sd(j, k, l) + beta(inner_step)*r(j, k, l)
        ENDDO
    ENDDO
  ENDDO
!$OMP END DO NOWAIT
  ENDIF
!$OMP END PARALLEL

  inner_step = inner_step + 1
  IF (inner_step .ge. tl_ppcg_inner_steps) EXIT

  END DO

END SUBROUTINE

SUBROUTINE tea_leaf_ppcg_calc_zrnorm_kernel(x_min, &
                          x_max,             &
                          y_min,             &
                          y_max,             &
                          z_min,             &
                          z_max,             &
                          halo_exchange_depth,             &
                          z, r,               &
                          preconditioner_type,    &
                          norm)

  IMPLICIT NONE

  INTEGER :: preconditioner_type
  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,z_min,z_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max,z_min-halo_exchange_depth:z_max+halo_exchange_depth) :: r, z
  REAL(KIND=8) :: norm
  integer :: j, k,l

  norm = 0.0_8

!$OMP PARALLEL
  IF (preconditioner_type .NE. TL_PREC_NONE) THEN
!$OMP DO REDUCTION(+:norm)
  DO l=z_min,z_max
    DO k=y_min,y_max
        DO j=x_min,x_max
            norm = norm + z(j, k, l)*r(j, k, l)
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

END MODULE


