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

!>  @brief Driver for the heat conduction kernel
!>  @author David Beckingsale, Wayne Gaudin
!>  @details Invokes the user specified kernel for the heat conduction

MODULE tea_leaf_module

  USE report_module
  USE data_module
  USE tea_leaf_kernel_common_module
  USE tea_leaf_kernel_jacobi_module
  USE tea_leaf_kernel_cg_module
  USE tea_leaf_kernel_ppcg_module
  USE tea_leaf_kernel_cheby_module
  USE update_halo_module

  IMPLICIT NONE

CONTAINS

SUBROUTINE tea_leaf()

  IMPLICIT NONE

!$ INTEGER :: OMP_GET_THREAD_NUM
  INTEGER :: c, n
  REAL(KIND=8) :: ry,rx,rz,error,exact_error,initial_residual

  INTEGER :: fields(NUM_FIELDS)

  REAL(KIND=8) :: timer,halo_time,solve_time,init_time,reset_time,dot_product_time

  ! For CG solver
  REAL(KIND=8) :: rro, pw, rrn, alpha, beta

  ! For chebyshev solver
  REAL(KIND=8), DIMENSION(max_iters) :: cg_alphas, cg_betas
  REAL(KIND=8), DIMENSION(max_iters) :: ch_alphas, ch_betas
  REAL(KIND=8),SAVE :: eigmin, eigmax, theta, cn
  INTEGER :: est_itc, cheby_calc_steps, max_cheby_iters, info, ppcg_inner_iters
  LOGICAL :: ch_switch_check
  LOGICAL, SAVE :: first=.TRUE.

  INTEGER :: cg_calc_steps

  REAL(KIND=8) :: cg_time, ch_time, total_solve_time, ch_per_it, cg_per_it, iteration_time

  cg_time = 0.0_8
  ch_time = 0.0_8
  cg_calc_steps = 0
  ppcg_inner_iters = 0
  ch_switch_check = .false.

  total_solve_time = 0.0_8
  init_time = 0.0_8
  halo_time = 0.0_8
  solve_time = 0.0_8
  dot_product_time = 0.0_8

  IF(coefficient .NE. RECIP_CONDUCTIVITY .AND. coefficient .NE. conductivity) THEN
    CALL report_error('tea_leaf', 'unknown coefficient option')
  ENDIF

  error = 1e10
  cheby_calc_steps = 0
  cg_calc_steps = 0

  initial_residual = 0.0_8

  total_solve_time = timer()

  DO c=1,chunks_per_task

    IF(chunks(c)%task.EQ.parallel%task) THEN

      ! INIT
      IF (profiler_on) init_time=timer()

      fields=0
      fields(FIELD_ENERGY1) = 1
      fields(FIELD_DENSITY) = 1

      IF (profiler_on) halo_time=timer()
      CALL update_halo(fields,halo_exchange_depth)
      IF (profiler_on) init_time = init_time + (timer()-halo_time)

      IF (use_fortran_kernels) THEN
        rx = dt/(chunks(c)%field%celldx(chunks(c)%field%x_min)**2)
        ry = dt/(chunks(c)%field%celldy(chunks(c)%field%y_min)**2)

        CALL tea_leaf_kernel_init_common(chunks(c)%field%x_min, &
            chunks(c)%field%x_max,                                  &
            chunks(c)%field%y_min,                                  &
            chunks(c)%field%y_max,                                  &
            chunks(c)%field%z_min,                                  &
            chunks(c)%field%z_max,                                  &
            halo_exchange_depth,                                  &
            chunks(c)%chunk_neighbours,                             &
            reflective_boundary,                                    &
            chunks(c)%field%density,                                &
            chunks(c)%field%energy1,                                &
            chunks(c)%field%u,                                      &
            chunks(c)%field%u0,                                      &
            chunks(c)%field%vector_r,                               &
            chunks(c)%field%vector_w,                               &
            chunks(c)%field%vector_Kx,                              &
            chunks(c)%field%vector_Ky,                              &
            chunks(c)%field%vector_Kz,                              &
            chunks(c)%field%tri_cp,   &
            chunks(c)%field%tri_bfp,    &
            chunks(c)%field%vector_Mi,    &
            rx, ry, rz, tl_preconditioner_type, coefficient)
      ENDIF

      fields=0
      fields(FIELD_U) = 1

      IF (profiler_on) halo_time=timer()
      CALL update_halo(fields,1)
      IF (profiler_on) init_time = init_time + (timer()-halo_time)

      IF(use_fortran_kernels) THEN
        CALL tea_leaf_calc_residual(chunks(c)%field%x_min,&
            chunks(c)%field%x_max,                       &
            chunks(c)%field%y_min,                       &
            chunks(c)%field%y_max,                       &
            chunks(c)%field%z_min,                       &
            chunks(c)%field%z_max,                       &
            halo_exchange_depth,                       &
            chunks(c)%field%u,                           &
            chunks(c)%field%u0,                 &
            chunks(c)%field%vector_r,                 &
            chunks(c)%field%vector_Kx,                 &
            chunks(c)%field%vector_Ky,                 &
            chunks(c)%field%vector_Kz,                 &
            rx, ry, rz)
        CALL tea_leaf_calc_2norm_kernel(chunks(c)%field%x_min,        &
            chunks(c)%field%x_max,                       &
            chunks(c)%field%y_min,                       &
            chunks(c)%field%y_max,                       &
            chunks(c)%field%z_min,                       &
            chunks(c)%field%z_max,                       &
            halo_exchange_depth,                       &
            chunks(c)%field%vector_r,                 &
            initial_residual)
      ENDIF

      IF (profiler_on) dot_product_time=timer()
      CALL tea_allsum(initial_residual)
      IF (profiler_on) profiler%dot_product= profiler%dot_product+ (timer() - dot_product_time)
      IF (profiler_on) init_time = init_time + (timer()-dot_product_time)

      initial_residual=SQRT(initial_residual)
      IF(parallel%boss.AND.verbose_on) THEN
!$      IF(OMP_GET_THREAD_NUM().EQ.0) THEN
          WRITE(g_out,*)"Initial residual ",initial_residual
!$      ENDIF
      ENDIF

      IF(tl_use_cg .OR. tl_use_chebyshev .OR. tl_use_ppcg) THEN
        ! All 3 of these solvers use the CG kernels
        IF(use_fortran_kernels) THEN
          CALL tea_leaf_kernel_init_cg_fortran(chunks(c)%field%x_min, &
              chunks(c)%field%x_max,                       &
              chunks(c)%field%y_min,                       &
              chunks(c)%field%y_max,                       &
              chunks(c)%field%z_min,                       &
              chunks(c)%field%z_max,                       &
              halo_exchange_depth,                       &
              chunks(c)%field%vector_p,                 &
              chunks(c)%field%vector_r,                 &
              chunks(c)%field%vector_Mi,                 &
              chunks(c)%field%vector_z,                 &
              chunks(c)%field%vector_Kx,                 &
              chunks(c)%field%vector_Ky,                 &
              chunks(c)%field%vector_Kz,                 &
              chunks(c)%field%tri_cp,   &
              chunks(c)%field%tri_bfp,    &
              rx, ry, rz, rro, tl_preconditioner_type)
        ENDIF

        ! and globally sum rro
        IF (profiler_on) dot_product_time=timer()
        CALL tea_allsum(rro)
        IF (profiler_on) profiler%dot_product= profiler%dot_product+ (timer() - dot_product_time)
        IF (profiler_on) init_time = init_time + (timer()-dot_product_time)

        ! need to update p when using CG due to matrix/vector multiplication
        fields=0
        fields(FIELD_U) = 1
        fields(FIELD_P) = 1

        IF (profiler_on) halo_time=timer()
        CALL update_halo(fields,1)
        IF (profiler_on) init_time=init_time+(timer()-halo_time)

        fields=0
        fields(FIELD_P) = 1
      ELSEIF (tl_use_jacobi) THEN
        fields=0
        fields(FIELD_U) = 1
      ENDIF

      IF (profiler_on) profiler%tea_init = profiler%tea_init + (timer() - init_time)

      IF (profiler_on) solve_time = timer()

      DO n=1,max_iters

        iteration_time = timer()

        IF (ch_switch_check .eqv. .false.) THEN
          IF (tl_ch_cg_errswitch) THEN
              ! either the ABS(error) has got below tolerance, or it's already going - minimum 20 steps to converge eigenvalues
              ch_switch_check = (cheby_calc_steps .GT. 0) .OR. (ABS(error) .LE. tl_ch_cg_epslim) .AND. (n .GE. 20)
          ELSE
              ! enough steps have passed and ABS(error) < 1, otherwise it's nowhere near converging on eigenvalues
              ch_switch_check = (n .GE. tl_ch_cg_presteps) .AND. (ABS(error) .le. 1.0_8)
          ENDIF
        ENDIF

        IF ((tl_use_chebyshev .OR. tl_use_ppcg) .AND. ch_switch_check) THEN
          ! on the first chebyshev steps, find the eigenvalues, coefficients,
          ! and expected number of iterations
          IF (cheby_calc_steps .EQ. 0) THEN
            ! maximum number of iterations in chebyshev solver
            max_cheby_iters = max_iters - n + 2
            rro = error

            IF(first) THEN
              ! calculate eigenvalues
              CALL tea_calc_eigenvalues(cg_alphas, cg_betas, eigmin, eigmax, &
                  max_iters, n-1, info)
              first=.FALSE.
              IF (info .NE. 0) CALL report_error('tea_leaf', 'Error in calculating eigenvalues')
              eigmin = eigmin * 0.95
              eigmax = eigmax * 1.05
            ENDIF

            IF (tl_use_chebyshev) THEN
              ! calculate chebyshev coefficients
              CALL tea_calc_ch_coefs(ch_alphas, ch_betas, eigmin, eigmax, &
                  theta, max_cheby_iters)

              ! don't need to update p any more
              fields = 0
              fields(FIELD_U) = 1
            ELSE IF (tl_use_ppcg) THEN
              ! To avoid some irritating bounds checking in ppcg inner iterations
              max_cheby_iters = max_cheby_iters + halo_exchange_depth

              ! currently also calculate chebyshev coefficients
              CALL tea_calc_ls_coefs(ch_alphas, ch_betas, eigmin, eigmax, &
                  theta, tl_ppcg_inner_steps)
            ENDIF

            cn = eigmax/eigmin

            IF (parallel%boss) THEN
!$            IF(OMP_GET_THREAD_NUM().EQ.0) THEN
100 FORMAT("Eigen min",e14.6," Eigen max",e14.6," Condition number",e14.6," Error",e14.6)
                WRITE(g_out,'(a,i3,a,e15.7)') "Switching after ",n," CG its, error ",rro
                WRITE(g_out, 100) eigmin,eigmax,cn,error
                WRITE(0,'(a,i3,a,e15.7)') "Switching after ",n," CG its, error ",rro
                WRITE(0, 100) eigmin,eigmax,cn,error
!$            ENDIF
            ENDIF
          ENDIF

          IF (tl_use_chebyshev) THEN
              IF (cheby_calc_steps .EQ. 0) THEN
                CALL tea_leaf_cheby_first_step(c, ch_alphas, ch_betas, fields, &
                    error, rx, ry, rz, theta, cn, max_cheby_iters, est_itc, solve_time)

                cheby_calc_steps = 1
              ELSE
                  IF(use_fortran_kernels) THEN
                      CALL tea_leaf_kernel_cheby_iterate(chunks(c)%field%x_min,&
                          chunks(c)%field%x_max,                       &
                          chunks(c)%field%y_min,                       &
                          chunks(c)%field%y_max,                       &
                          chunks(c)%field%z_min,                       &
                          chunks(c)%field%z_max,                       &
                          halo_exchange_depth,                       &
                          chunks(c)%field%u,                           &
                          chunks(c)%field%u0,                          &
                          chunks(c)%field%vector_p,                 &
                          chunks(c)%field%vector_r,                 &
                          chunks(c)%field%vector_Mi,                 &
                          chunks(c)%field%vector_w,                 &
                          chunks(c)%field%vector_z,                 &
                          chunks(c)%field%vector_Kx,                 &
                          chunks(c)%field%vector_Ky,                 &
                          chunks(c)%field%vector_Kz,                 &
                          chunks(c)%field%tri_cp,   &
                          chunks(c)%field%tri_bfp,    &
                          ch_alphas, ch_betas, max_cheby_iters,        &
                          rx, ry, rz, cheby_calc_steps, tl_preconditioner_type)
                  ENDIF

                  ! after estimated number of iterations has passed, calc resid.
                  ! Leaving 10 iterations between each global reduction won't affect
                  ! total time spent much if at all (number of steps spent in
                  ! chebyshev is typically O(300+)) but will greatly reduce global
                  ! synchronisations needed
                  IF ((n .GE. est_itc) .AND. (MOD(n, 10) .eq. 0)) THEN
                    IF(use_fortran_kernels) THEN
                      CALL tea_leaf_calc_2norm_kernel(chunks(c)%field%x_min,        &
                            chunks(c)%field%x_max,                       &
                            chunks(c)%field%y_min,                       &
                            chunks(c)%field%y_max,                       &
                            chunks(c)%field%z_min,                       &
                            chunks(c)%field%z_max,                       &
                            halo_exchange_depth,                       &
                            chunks(c)%field%vector_r,                 &
                            error)
                    ENDIF

                    IF (profiler_on) dot_product_time=timer()
                    CALL tea_allsum(error)
                    IF (profiler_on) profiler%dot_product= profiler%dot_product+ (timer() - dot_product_time)
                    IF (profiler_on) solve_time = solve_time + (timer()-dot_product_time)
                  ENDIF
              ENDIF
          ELSE IF (tl_use_ppcg) THEN
            IF (cheby_calc_steps .EQ. 0) THEN
              cheby_calc_steps = 1

              fields=0
              fields(FIELD_U) = 1

              IF (profiler_on) halo_time=timer()
              CALL update_halo(fields,1)
              IF (profiler_on) solve_time = solve_time + (timer()-halo_time)

              IF(use_fortran_kernels) THEN
                CALL tea_leaf_calc_residual(chunks(c)%field%x_min,&
                    chunks(c)%field%x_max,                       &
                    chunks(c)%field%y_min,                       &
                    chunks(c)%field%y_max,                       &
                    chunks(c)%field%z_min,                       &
                    chunks(c)%field%z_max,                       &
                    halo_exchange_depth,                       &
                    chunks(c)%field%u,                           &
                    chunks(c)%field%u0,                 &
                    chunks(c)%field%vector_r,                 &
                    chunks(c)%field%vector_Kx,                 &
                    chunks(c)%field%vector_Ky,                 &
                    chunks(c)%field%vector_Kz,                 &
                    rx, ry, rz)
              ENDIF
            ENDIF

            IF(use_fortran_kernels) THEN
              CALL tea_leaf_kernel_solve_cg_fortran_calc_w(chunks(c)%field%x_min,&
                  chunks(c)%field%x_max,                       &
                  chunks(c)%field%y_min,                       &
                  chunks(c)%field%y_max,                       &
                  chunks(c)%field%z_min,                       &
                  chunks(c)%field%z_max,                       &
                  halo_exchange_depth,                       &
                  chunks(c)%field%vector_p,                 &
                  chunks(c)%field%vector_w,                 &
                  chunks(c)%field%vector_Kx,                 &
                  chunks(c)%field%vector_Ky,                 &
                  chunks(c)%field%vector_Kz,                 &
                  rx, ry, rz, pw)
            ENDIF

            IF (profiler_on) dot_product_time=timer()
            CALL tea_allsum(pw)
            IF (profiler_on) profiler%dot_product= profiler%dot_product+ (timer() - dot_product_time)
            IF (profiler_on) solve_time = solve_time + (timer()-dot_product_time)

            alpha = rro/pw

            IF(use_fortran_kernels) THEN
              CALL tea_leaf_kernel_solve_cg_fortran_calc_ur(chunks(c)%field%x_min,&
                  chunks(c)%field%x_max,                                          &
                  chunks(c)%field%y_min,                                          &
                  chunks(c)%field%y_max,                                          &
                  chunks(c)%field%z_min,                                          &
                  chunks(c)%field%z_max,                                          &
                  halo_exchange_depth,                                          &
                  chunks(c)%field%u,                                              &
                  chunks(c)%field%vector_p,                                       &
                  chunks(c)%field%vector_r,                                       &
                  chunks(c)%field%vector_Mi,                                      &
                  chunks(c)%field%vector_w,                                       &
                  chunks(c)%field%vector_z,                                       &
                  chunks(c)%field%tri_cp,   &
                  chunks(c)%field%tri_bfp,    &
                  chunks(c)%field%vector_Kx,                              &
                  chunks(c)%field%vector_Ky,                              &
                  chunks(c)%field%vector_Kz,                              &
                  rx, ry, rz, &
                  alpha, rrn, tl_preconditioner_type)
            ENDIF

            ! not using rrn, so don't do a tea_allsum

            CALL tea_leaf_run_ppcg_inner_steps(ch_alphas, ch_betas, theta, &
                rx, ry, rz, tl_ppcg_inner_steps, c, solve_time)
            ppcg_inner_iters = ppcg_inner_iters + tl_ppcg_inner_steps

            IF(use_fortran_kernels) THEN
              CALL tea_leaf_ppcg_calc_zrnorm_kernel(chunks(c)%field%x_min,        &
                    chunks(c)%field%x_max,                       &
                    chunks(c)%field%y_min,                       &
                    chunks(c)%field%y_max,                       &
                    chunks(c)%field%z_min,                       &
                    chunks(c)%field%z_max,                       &
                    halo_exchange_depth,                       &
                    chunks(c)%field%vector_z,                 &
                    chunks(c)%field%vector_r,                 &
                    tl_preconditioner_type, rrn)
            ENDIF

            IF (profiler_on) dot_product_time=timer()
            CALL tea_allsum(rrn)
            IF (profiler_on) profiler%dot_product= profiler%dot_product+ (timer() - dot_product_time)
            IF (profiler_on) solve_time = solve_time + (timer()-dot_product_time)

            beta = rrn/rro

            IF(use_fortran_kernels) THEN
              CALL tea_leaf_kernel_solve_cg_fortran_calc_p(chunks(c)%field%x_min,&
                  chunks(c)%field%x_max,                       &
                  chunks(c)%field%y_min,                       &
                  chunks(c)%field%y_max,                       &
                  chunks(c)%field%z_min,                       &
                  chunks(c)%field%z_max,                       &
                  halo_exchange_depth,                       &
                  chunks(c)%field%vector_p,                 &
                  chunks(c)%field%vector_r,                 &
                  chunks(c)%field%vector_z,                 &
                  beta, tl_preconditioner_type)
            ENDIF

            error = rrn
            rro = rrn
          ENDIF

          cheby_calc_steps = cheby_calc_steps + 1
        ELSEIF(tl_use_cg .OR. tl_use_chebyshev .OR. tl_use_ppcg) THEN
          fields(FIELD_P) = 1
          cg_calc_steps = cg_calc_steps + 1

          pw = 0.0_08

          IF(use_fortran_kernels) THEN
            CALL tea_leaf_kernel_solve_cg_fortran_calc_w(chunks(c)%field%x_min,&
                chunks(c)%field%x_max,                       &
                chunks(c)%field%y_min,                       &
                chunks(c)%field%y_max,                       &
                chunks(c)%field%z_min,                       &
                chunks(c)%field%z_max,                       &
                halo_exchange_depth,                       &
                chunks(c)%field%vector_p,                 &
                chunks(c)%field%vector_w,                 &
                chunks(c)%field%vector_Kx,                 &
                chunks(c)%field%vector_Ky,                 &
                chunks(c)%field%vector_Kz,                 &
                rx, ry, rz, pw)
          ENDIF

          IF (profiler_on) dot_product_time=timer()
          CALL tea_allsum(pw)
          IF (profiler_on) profiler%dot_product= profiler%dot_product+ (timer() - dot_product_time)
          IF (profiler_on) solve_time = solve_time + (timer()-dot_product_time)

          alpha = rro/pw
          cg_alphas(n) = alpha

          rrn = 0.0_8

          IF(use_fortran_kernels) THEN
            CALL tea_leaf_kernel_solve_cg_fortran_calc_ur(chunks(c)%field%x_min,&
                chunks(c)%field%x_max,                       &
                chunks(c)%field%y_min,                       &
                chunks(c)%field%y_max,                       &
                chunks(c)%field%z_min,                       &
                chunks(c)%field%z_max,                       &
                halo_exchange_depth,                       &
                chunks(c)%field%u,                           &
                chunks(c)%field%vector_p,                 &
                chunks(c)%field%vector_r,                 &
                chunks(c)%field%vector_Mi,                 &
                chunks(c)%field%vector_w,                 &
                chunks(c)%field%vector_z,                 &
                chunks(c)%field%tri_cp,   &
                chunks(c)%field%tri_bfp,    &
                chunks(c)%field%vector_Kx,                              &
                chunks(c)%field%vector_Ky,                              &
                chunks(c)%field%vector_Kz,                              &
                rx, ry, rz, &
                alpha, rrn, tl_preconditioner_type)
          ENDIF

          IF (profiler_on) dot_product_time=timer()
          CALL tea_allsum(rrn)
          IF (profiler_on) profiler%dot_product= profiler%dot_product+ (timer() - dot_product_time)
          IF (profiler_on) solve_time = solve_time + (timer()-dot_product_time)

          beta = rrn/rro
          cg_betas(n) = beta

          IF(use_fortran_kernels) THEN
            CALL tea_leaf_kernel_solve_cg_fortran_calc_p(chunks(c)%field%x_min,&
                chunks(c)%field%x_max,                       &
                chunks(c)%field%y_min,                       &
                chunks(c)%field%y_max,                       &
                chunks(c)%field%z_min,                       &
                chunks(c)%field%z_max,                       &
                halo_exchange_depth,                       &
                chunks(c)%field%vector_p,                 &
                chunks(c)%field%vector_r,                 &
                chunks(c)%field%vector_z,                 &
                beta, tl_preconditioner_type)
          ENDIF

          error = rrn
          rro = rrn
        ELSEIF(tl_use_jacobi) THEN
          IF(use_fortran_kernels) THEN
            CALL tea_leaf_kernel_jacobi_solve(chunks(c)%field%x_min,&
                chunks(c)%field%x_max,                       &
                chunks(c)%field%y_min,                       &
                chunks(c)%field%y_max,                       &
                chunks(c)%field%z_min,                       &
                chunks(c)%field%z_max,                       &
                halo_exchange_depth,                       &
                rx,                                          &
                ry,                                          &
                rz,                                          &
                chunks(c)%field%vector_Kx,                 &
                chunks(c)%field%vector_Ky,                 &
                chunks(c)%field%vector_Kz,                 &
                error, &
                chunks(c)%field%u0,                          &
                chunks(c)%field%u,                           &
                chunks(c)%field%vector_r)
          ENDIF

          ! error for jacobi is calculated recursively and is not very accurate,
          ! so do this every so often to see whether it has actually converged
          IF (MOD(n, 50) .EQ. 0) THEN
            IF (profiler_on) halo_time = timer()
            CALL update_halo(fields,1)
            IF (profiler_on) solve_time = solve_time + (timer()-halo_time)

            IF(use_fortran_kernels) THEN
              CALL tea_leaf_calc_residual(chunks(c)%field%x_min,&
                  chunks(c)%field%x_max,                       &
                  chunks(c)%field%y_min,                       &
                  chunks(c)%field%y_max,                       &
                  chunks(c)%field%z_min,                       &
                  chunks(c)%field%z_max,                       &
                  halo_exchange_depth,                       &
                  chunks(c)%field%u,                           &
                  chunks(c)%field%u0,                 &
                  chunks(c)%field%vector_r,                 &
                  chunks(c)%field%vector_Kx,                 &
                  chunks(c)%field%vector_Ky,                 &
                  chunks(c)%field%vector_Kz,                 &
                  rx, ry, rz)
              CALL tea_leaf_calc_2norm_kernel(chunks(c)%field%x_min,        &
                  chunks(c)%field%x_max,                       &
                  chunks(c)%field%y_min,                       &
                  chunks(c)%field%y_max,                       &
                  chunks(c)%field%z_min,                       &
                  chunks(c)%field%z_max,                       &
                  halo_exchange_depth,                       &
                  chunks(c)%field%vector_r,                 &
                  error)
            ENDIF
          ENDIF

          IF (profiler_on) dot_product_time=timer()
          CALL tea_allsum(error)
          IF (profiler_on) profiler%dot_product= profiler%dot_product+ (timer() - dot_product_time)
          IF (profiler_on) solve_time = solve_time + (timer()-dot_product_time)
        ENDIF

        ! updates u and possibly p
        IF (profiler_on) halo_time = timer()
        CALL update_halo(fields,1)
        IF (profiler_on) solve_time = solve_time + (timer()-halo_time)

        IF (profiler_on) THEN
          IF (tl_use_chebyshev .AND. ch_switch_check) THEN
            ch_time=ch_time+(timer()-iteration_time)
          ELSE
            cg_time=cg_time+(timer()-iteration_time)
          ENDIF
        ENDIF

        error=SQRT(error)
        IF(parallel%boss.AND.verbose_on) THEN
!$        IF(OMP_GET_THREAD_NUM().EQ.0) THEN
            WRITE(g_out,*)"Residual ",error
!$        ENDIF
        ENDIF

        IF (abs(error) .LT. eps*initial_residual) EXIT

      ENDDO

      IF (tl_check_result) THEN
        fields = 0
        fields(FIELD_U) = 1

        IF (profiler_on) halo_time = timer()
        CALL update_halo(fields,1)
        IF (profiler_on) solve_time = solve_time + (timer()-halo_time)

        IF(use_fortran_kernels) THEN
          CALL tea_leaf_calc_residual(chunks(c)%field%x_min,&
              chunks(c)%field%x_max,                       &
              chunks(c)%field%y_min,                       &
              chunks(c)%field%y_max,                       &
              chunks(c)%field%z_min,                       &
              chunks(c)%field%z_max,                       &
              halo_exchange_depth,                       &
              chunks(c)%field%u,                           &
              chunks(c)%field%u0,                 &
              chunks(c)%field%vector_r,                 &
              chunks(c)%field%vector_Kx,                 &
              chunks(c)%field%vector_Ky,                 &
              chunks(c)%field%vector_Kz,                 &
              rx, ry, rz)
          CALL tea_leaf_calc_2norm_kernel(chunks(c)%field%x_min,        &
              chunks(c)%field%x_max,                       &
              chunks(c)%field%y_min,                       &
              chunks(c)%field%y_max,                       &
              chunks(c)%field%z_min,                       &
              chunks(c)%field%z_max,                       &
              halo_exchange_depth,                       &
              chunks(c)%field%vector_r,                 &
              exact_error)
        ENDIF

        IF (profiler_on) dot_product_time=timer()
        CALL tea_allsum(exact_error)
        IF (profiler_on) profiler%dot_product= profiler%dot_product+ (timer() - dot_product_time)
        IF (profiler_on) solve_time = solve_time + (timer()-dot_product_time)

        exact_error = SQRT(exact_error)
      ENDIF

      IF (profiler_on) profiler%tea_solve = profiler%tea_solve + (timer() - solve_time)

      IF (parallel%boss) THEN
!$      IF(OMP_GET_THREAD_NUM().EQ.0) THEN

102 FORMAT('Conduction error ',e14.7)
          WRITE(g_out,102) error/initial_residual
          WRITE(0,102) error/initial_residual

          IF (tl_check_result) THEN
101 FORMAT('EXACT error calculated as', e14.7)
            WRITE(0, 101) exact_error/initial_residual
            WRITE(g_out, 101) exact_error/initial_residual
          ENDIF

          WRITE(g_out,"('Iteration count ',i8)") n-1
          WRITE(0,"('Iteration count ', i8)") n-1
          IF(tl_use_ppcg) THEN
103 FORMAT('PPCG Iteration count', i8, ' (Total ',i8,')')
            WRITE(g_out,103) ppcg_inner_iters, ppcg_inner_iters + n-1
            WRITE(0,103) ppcg_inner_iters, ppcg_inner_iters + n-1
          ENDIF
!$      ENDIF
      ENDIF

      ! RESET
      IF (profiler_on) reset_time=timer()

      IF(use_fortran_kernels) THEN
          CALL tea_leaf_kernel_finalise(chunks(c)%field%x_min, &
              chunks(c)%field%x_max,                           &
              chunks(c)%field%y_min,                           &
              chunks(c)%field%y_max,                           &
              chunks(c)%field%z_min,                           &
              chunks(c)%field%z_max,                           &
              halo_exchange_depth,                           &
              chunks(c)%field%energy1,                         &
              chunks(c)%field%density,                        &
              chunks(c)%field%u)
      ENDIF

      fields=0
      fields(FIELD_ENERGY1) = 1

      IF (profiler_on) halo_time=timer()
      CALL update_halo(fields,1)
      IF (profiler_on) reset_time = reset_time + (timer()-halo_time)

      IF (profiler_on) profiler%tea_reset = profiler%tea_reset + (timer() - reset_time)

    ENDIF

  ENDDO

  IF (profiler_on .AND. parallel%boss) THEN
    total_solve_time = (timer() - total_solve_time)
    WRITE(0, "(a16,f16.10,a7,i7,a16,f16.10)") "Solve Time",total_solve_time,"Its",n,"Time Per It",total_solve_time/n
    WRITE(g_out, "(a16,f16.10,a7,i7,a16,f16.10)") "Solve Time",total_solve_time,"Its",n,"Time Per It",total_solve_time/n
  ENDIF

  IF (profiler_on .AND. tl_use_chebyshev) THEN
    CALL tea_sum(ch_time)
    CALL tea_sum(cg_time)
    IF (parallel%boss) THEN
      cg_per_it = MERGE((cg_time/cg_calc_steps)/parallel%max_task, 0.0_8, cg_calc_steps .GT. 0)
      ch_per_it = MERGE((ch_time/cheby_calc_steps)/parallel%max_task, 0.0_8, cheby_calc_steps .GT. 0)

      WRITE(0, "(a3, a16, a7, a16, a7)") "", "Time", "Its", "Per it", "Ratio"
      WRITE(0, "(a3, f16.10, i7, f16.10, f7.2)") &
          "CG", cg_time + 0.0_8, cg_calc_steps, cg_per_it, 1.0_8
      WRITE(0, "(a3, f16.10, i7, f16.10, f7.2)") "CH", ch_time + 0.0_8, cheby_calc_steps, &
          ch_per_it, MERGE(ch_per_it/cg_per_it, 0.0_8, cheby_calc_steps .GT. 0)
      WRITE(0, "('Chebyshev actually took ', i6, ' (' i6, ' off guess)')") &
          cheby_calc_steps, cheby_calc_steps-est_itc

      WRITE(g_out, "(a3, a16, a7, a16, a7)") "", "Time", "Its", "Per it", "Ratio"
      WRITE(g_out, "(a3, f16.10, i7, f16.10, f7.2)") &
          "CG", cg_time + 0.0_8, cg_calc_steps, cg_per_it, 1.0_8
      WRITE(g_out, "(a3, f16.10, i7, f16.10, f7.2)") "CH", ch_time + 0.0_8, cheby_calc_steps, &
          ch_per_it, MERGE(ch_per_it/cg_per_it, 0.0_8, cheby_calc_steps .GT. 0)
      WRITE(g_out, "('Chebyshev actually took ', i6, ' (' i6, ' off guess)')") &
          cheby_calc_steps, cheby_calc_steps-est_itc
    ENDIF
  ENDIF

END SUBROUTINE tea_leaf

SUBROUTINE tea_leaf_run_ppcg_inner_steps(ch_alphas, ch_betas, theta, &
    rx, ry, rz, tl_ppcg_inner_steps, c, solve_time)

  INTEGER :: fields(NUM_FIELDS)
  INTEGER :: c, tl_ppcg_inner_steps, ppcg_cur_step
  REAL(KIND=8) :: rx, ry, rz, theta
  REAL(KIND=8) :: halo_time, timer, solve_time
  REAL(KIND=8), DIMENSION(max_iters) :: ch_alphas, ch_betas

  INTEGER(KIND=4) :: x_min_bound, x_max_bound, y_min_bound, y_max_bound, z_min_bound, z_max_bound, inner_step, bounds_extra

  fields = 0
  fields(FIELD_U) = 1

  IF (profiler_on) halo_time=timer()
  CALL update_halo(fields,1)
  IF (profiler_on) solve_time = solve_time + (timer() - halo_time)

  IF(use_fortran_kernels) THEN
     CALL tea_leaf_calc_residual(chunks(c)%field%x_min,&
         chunks(c)%field%x_max,                        &
         chunks(c)%field%y_min,                        &
         chunks(c)%field%y_max,                        &
        chunks(c)%field%z_min,                       &
        chunks(c)%field%z_max,                       &
        halo_exchange_depth,                       &
        chunks(c)%field%u,                           &
        chunks(c)%field%u0,                 &
        chunks(c)%field%vector_r,                 &
        chunks(c)%field%vector_Kx,                 &
        chunks(c)%field%vector_Ky,                 &
        chunks(c)%field%vector_Kz,                 &
        rx, ry, rz)
    CALL tea_leaf_kernel_ppcg_init_sd(chunks(c)%field%x_min,&
        chunks(c)%field%x_max,                              &
        chunks(c)%field%y_min,                              &
        chunks(c)%field%y_max,                              &
        chunks(c)%field%z_min,                              &
        chunks(c)%field%z_max,                              &
        halo_exchange_depth,                              &
        chunks(c)%field%vector_r,                           &
        chunks(c)%field%vector_Kx,                        &
        chunks(c)%field%vector_Ky,                        &
        chunks(c)%field%vector_Kz,                        &
        chunks(c)%field%vector_sd,                          &
        chunks(c)%field%vector_z,                          &
        chunks(c)%field%tri_cp,                          &
        chunks(c)%field%tri_bfp,                          &
        chunks(c)%field%vector_Mi,                          &
        rx, ry, rz,                      &
        theta, tl_preconditioner_type)
  ENDIF

  ! inner steps
  DO ppcg_cur_step=1,tl_ppcg_inner_steps,halo_exchange_depth

    fields = 0
    fields(FIELD_SD) = 1
    fields(FIELD_R) = 1

    IF (profiler_on) halo_time = timer()
    CALL update_halo(fields,halo_exchange_depth)
    IF (profiler_on) solve_time = solve_time + (timer()-halo_time)

    inner_step = ppcg_cur_step

    DO bounds_extra = halo_exchange_depth-1, 0, -1
      IF (chunks(c)%chunk_neighbours(CHUNK_LEFT).EQ.EXTERNAL_FACE) THEN
        x_min_bound = chunks(c)%field%x_min
      ELSE
        x_min_bound = chunks(c)%field%x_min - bounds_extra
      ENDIF

      IF (chunks(c)%chunk_neighbours(CHUNK_RIGHT).EQ.EXTERNAL_FACE) THEN
        x_max_bound = chunks(c)%field%x_max
      ELSE
        x_max_bound = chunks(c)%field%x_max + bounds_extra
      ENDIF

      IF (chunks(c)%chunk_neighbours(CHUNK_BOTTOM).EQ.EXTERNAL_FACE) THEN
        y_min_bound = chunks(c)%field%y_min
      ELSE
        y_min_bound = chunks(c)%field%y_min - bounds_extra
      ENDIF

      IF (chunks(c)%chunk_neighbours(CHUNK_TOP).EQ.EXTERNAL_FACE) THEN
        y_max_bound = chunks(c)%field%y_max
      ELSE
        y_max_bound = chunks(c)%field%y_max + bounds_extra
      ENDIF

      IF (chunks(c)%chunk_neighbours(CHUNK_BACK).EQ.EXTERNAL_FACE) THEN
        z_min_bound = chunks(c)%field%z_min
      ELSE
        z_min_bound = chunks(c)%field%z_min - bounds_extra
      ENDIF

      IF (chunks(c)%chunk_neighbours(CHUNK_FRONT).EQ.EXTERNAL_FACE) THEN
        z_max_bound = chunks(c)%field%z_max
      ELSE
        z_max_bound = chunks(c)%field%z_max + bounds_extra
      ENDIF

      IF(use_fortran_kernels) THEN
        CALL tea_leaf_kernel_ppcg_inner(chunks(c)%field%x_min,&
            chunks(c)%field%x_max,                            &
            chunks(c)%field%y_min,                            &
            chunks(c)%field%y_max,                            &
            chunks(c)%field%z_min,                            &
            chunks(c)%field%z_max,                            &
            halo_exchange_depth,                            &
            x_min_bound,                                    &
            x_max_bound,                                    &
            y_min_bound,                                    &
            y_max_bound,                                    &
            z_min_bound,                                    &
            z_max_bound,                                    &
            ch_alphas, ch_betas,                              &
            rx, ry, rz,                                       &
            inner_step,                                     &
            chunks(c)%field%u,                                &
            chunks(c)%field%vector_r,                         &
            chunks(c)%field%vector_Kx,                        &
            chunks(c)%field%vector_Ky,                        &
            chunks(c)%field%vector_Kz,                        &
            chunks(c)%field%vector_sd,                        &
            chunks(c)%field%vector_z,                          &
            chunks(c)%field%tri_cp,                          &
            chunks(c)%field%tri_bfp,                          &
            chunks(c)%field%vector_Mi,                          &
            tl_preconditioner_type)
      ENDIF

      fields = 0
      fields(FIELD_SD) = 1

      IF (profiler_on) halo_time = timer()
      CALL update_boundary(fields, bounds_extra)
      IF (profiler_on) solve_time = solve_time + (timer()-halo_time)

      inner_step = inner_step + 1
      IF (inner_step .gt. tl_ppcg_inner_steps) EXIT
    ENDDO
  ENDDO

  fields = 0
  fields(FIELD_P) = 1

END SUBROUTINE tea_leaf_run_ppcg_inner_steps

SUBROUTINE tea_leaf_cheby_first_step(c, ch_alphas, ch_betas, fields, &
    error, rx, ry, rz, theta, cn, max_cheby_iters, est_itc, solve_time)

  IMPLICIT NONE

  integer :: c, est_itc, max_cheby_iters
  integer, dimension(:) :: fields
  REAL(KIND=8) :: it_alpha, cn, gamm, bb, error, rx, ry, rz, theta
  REAL(KIND=8), DIMENSION(:) :: ch_alphas, ch_betas
  REAL(KIND=8) :: halo_time, timer, dot_product_time, solve_time

  ! calculate 2 norm of u0
  IF(use_fortran_kernels) THEN
    CALL tea_leaf_calc_2norm_kernel(chunks(c)%field%x_min,        &
          chunks(c)%field%x_max,                       &
          chunks(c)%field%y_min,                       &
          chunks(c)%field%y_max,                       &
          chunks(c)%field%z_min,                       &
          chunks(c)%field%z_max,                       &
          halo_exchange_depth,                       &
          chunks(c)%field%u0,                 &
          bb)
  ENDIF

  IF (profiler_on) dot_product_time=timer()
  CALL tea_allsum(bb)
  IF (profiler_on) profiler%dot_product= profiler%dot_product+ (timer() - dot_product_time)
  IF (profiler_on) solve_time = solve_time + (timer()-dot_product_time)

  ! initialise 'p' array
  IF(use_fortran_kernels) THEN
    CALL tea_leaf_kernel_cheby_init(chunks(c)%field%x_min,&
          chunks(c)%field%x_max,                          &
          chunks(c)%field%y_min,                          &
          chunks(c)%field%y_max,                          &
          chunks(c)%field%z_min,                          &
          chunks(c)%field%z_max,                          &
          halo_exchange_depth,                          &
          chunks(c)%field%u,                              &
          chunks(c)%field%u0,                             &
          chunks(c)%field%vector_p,                       &
          chunks(c)%field%vector_r,                       &
          chunks(c)%field%vector_Mi,                      &
          chunks(c)%field%vector_w,                       &
          chunks(c)%field%vector_z,                       &
          chunks(c)%field%vector_Kx,                      &
          chunks(c)%field%vector_Ky,                      &
          chunks(c)%field%vector_Kz,                      &
          chunks(c)%field%tri_cp,   &
          chunks(c)%field%tri_bfp,    &
          rx, ry, rz, theta, tl_preconditioner_type)
  ENDIF

  IF (profiler_on) halo_time = timer()
  CALL update_halo(fields,1)
  IF (profiler_on) solve_time = solve_time + (timer()-halo_time)

  IF(use_fortran_kernels) THEN
      CALL tea_leaf_kernel_cheby_iterate(chunks(c)%field%x_min,&
          chunks(c)%field%x_max,                               &
          chunks(c)%field%y_min,                               &
          chunks(c)%field%y_max,                               &
          chunks(c)%field%z_min,                       &
          chunks(c)%field%z_max,                       &
          halo_exchange_depth,                       &
          chunks(c)%field%u,                           &
          chunks(c)%field%u0,                          &
          chunks(c)%field%vector_p,                 &
          chunks(c)%field%vector_r,                 &
          chunks(c)%field%vector_Mi,                 &
          chunks(c)%field%vector_w,                 &
          chunks(c)%field%vector_z,                 &
          chunks(c)%field%vector_Kx,                 &
          chunks(c)%field%vector_Ky,                 &
          chunks(c)%field%vector_Kz,                 &
          chunks(c)%field%tri_cp,   &
          chunks(c)%field%tri_bfp,    &
          ch_alphas, ch_betas, max_cheby_iters,        &
          rx, ry, rz, 1, tl_preconditioner_type)
  ENDIF

  IF(use_fortran_kernels) THEN
    CALL tea_leaf_calc_2norm_kernel(chunks(c)%field%x_min,        &
          chunks(c)%field%x_max,                       &
          chunks(c)%field%y_min,                       &
          chunks(c)%field%y_max,                       &
          chunks(c)%field%z_min,                       &
          chunks(c)%field%z_max,                       &
          halo_exchange_depth,                       &
          chunks(c)%field%vector_r,                 &
          error)
  ENDIF

  IF (profiler_on) dot_product_time=timer()
  CALL tea_allsum(error)
  IF (profiler_on) profiler%dot_product= profiler%dot_product+ (timer() - dot_product_time)
  IF (profiler_on) solve_time = solve_time + (timer()-dot_product_time)

  it_alpha = EPSILON(1.0_8)*bb/(4.0_8*error)
  gamm = (SQRT(cn) - 1.0_8)/(SQRT(cn) + 1.0_8)
  est_itc = NINT(LOG(it_alpha)/(2.0_8*LOG(gamm)))

  IF (parallel%boss) THEN
      WRITE(g_out,'(a11)')"est itc"
      WRITE(g_out,'(11i11)')est_itc
      WRITE(0,'(a11)')"est itc"
      WRITE(0,'(11i11)')est_itc
  ENDIF

END SUBROUTINE tea_leaf_cheby_first_step

END MODULE tea_leaf_module
