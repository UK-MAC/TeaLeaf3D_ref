MODULE tea_leaf_kernel_ppcg_module

USE tea_leaf_kernel_module
USE tea_leaf_kernel_cheby_module

IMPLICIT NONE

CONTAINS

SUBROUTINE tea_leaf_kernel_ppcg_init_sd(x_min,             &
                           x_max,             &
                           y_min,             &
                           y_max,             &
                           z_min,             &
                           z_max,             &
                           r,            &
                           sd, &
                           theta)

  IMPLICIT NONE

  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,z_min,z_max
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+2) :: r, sd
  REAL(KIND=8) :: theta

  INTEGER :: j,k,l

!$OMP PARALLEL
!$OMP DO
  DO l=z_min,z_max
    DO k=y_min,y_max
        DO j=x_min,x_max
            sd(j, k, l) = r(j, k, l)/theta
        ENDDO
    ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE tea_leaf_kernel_ppcg_init_sd

SUBROUTINE tea_leaf_kernel_ppcg_inner(x_min,             &
                           x_max,             &
                           y_min,             &
                           y_max,             &
                           z_min,             &
                           z_max,             &
                           step, &
                           alpha, &
                           beta, &
                           rx, ry, rz, &
                           u, &
                           r, &
                           Kx, &
                           Ky, &
                           Kz, &
                           sd)
  IMPLICIT NONE

  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,z_min,z_max
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+2) :: u, r, Kx, Ky, sd, kz
  INTEGER(KIND=4) :: j,k,l, step
  REAL(KIND=8), DIMENSION(:) :: alpha, beta
  REAL(KIND=8) :: smvp, rx, ry, rz

!$OMP PARALLEL
!$OMP DO
  DO l=z_min,z_max
    DO k=y_min,y_max
        DO j=x_min,x_max
            smvp = (1.0_8                                      &
                + rx*(Kx(j+1, k, l) + Kx(j, k, l)) &
                + ry*(Ky(j, k+1, l) + Ky(j, k, l))                      &
                + rz*(Kz(j, k, l+1) + Kz(j, k, l)))*sd(j, k, l)             &
                - rx*(Kx(j+1, k, l)*sd(j+1, k, l) + Kx(j, k, l)*sd(j-1, k, l)) &
                - ry*(Ky(j, k+1, l)*sd(j, k+1, l) + Ky(j, k, l)*sd(j, k-1, l))  &
                - rz*(Kz(j, k, l+1)*sd(j, k, l+1) + Kz(j, k, l)*sd(j, k, l-1))
            r(j, k, l) = r(j, k, l) - smvp
            u(j, k, l) = u(j, k, l) + sd(j, k, l)
        ENDDO
    ENDDO
  ENDDO
!$OMP END DO
!$OMP DO
  DO l=z_min,z_max
    DO k=y_min,y_max
        DO j=x_min,x_max
            sd(j, k, l) = alpha(step)*sd(j, k, l) + beta(step)*r(j, k, l)
        ENDDO
    ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

end SUBROUTINE

SUBROUTINE tea_leaf_kernel_ppcg_init_p(x_min,             &
                           x_max,             &
                           y_min,             &
                           y_max,             &
                           z_min,             &
                           z_max,             &
                           p,                &
                           r,            &
                           rro)

  IMPLICIT NONE

  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,z_min,z_max
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+2) :: p, r

  INTEGER :: j,k,l
  REAL(KIND=8) ::  rro

  rro = 0.0_8

!$OMP PARALLEL
!$OMP DO REDUCTION(+:rro)
  DO l=z_min,z_max
  DO k=y_min,y_max
    DO j=x_min,x_max
      p(j, k, l) = r(j, k, l)
      rro = rro + p(j, k, l)*r(j, k, l)
    ENDDO
  ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE

SUBROUTINE tea_calc_ls_coefs(ch_alphas, ch_betas, eigmin, eigmax, &
                             theta, ppcg_inner_steps)

  INTEGER :: n, ppcg_inner_steps
  REAL(KIND=8), DIMENSION(ppcg_inner_steps) :: ch_alphas, ch_betas
  REAL(KIND=8) :: eigmin, eigmax

  REAL(KIND=8) :: theta, delta, sigma, rho_old, rho_new, cur_alpha, cur_beta

  ! TODO
  CALL tea_calc_ch_coefs(ch_alphas, ch_betas, eigmin, eigmax, &
                         theta, ppcg_inner_steps)

END SUBROUTINE

END MODULE


