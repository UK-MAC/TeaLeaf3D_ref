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

!>  @brief Communication Utilities
!>  @author David Beckingsale, Wayne Gaudin
!>  @details Contains all utilities required to run TeaLeaf in a distributed
!>  environment, including initialisation, mesh decompostion, reductions and
!>  halo exchange using explicit buffers.
!>
!>  Note the halo exchange is currently coded as simply as possible and no
!>  optimisations have been implemented, such as post receives before sends or packing
!>  buffers with multiple data fields. This is intentional so the effect of these
!>  optimisations can be measured on large systems, as and when they are added.
!>
!>  Even without these modifications TeaLeaf weak scales well on moderately sized
!>  systems of the order of 10K cores.

MODULE tea_module

  USE data_module
  USE definitions_module
  USE pack_module
  !USE MPI

  IMPLICIT NONE

  include "mpif.h"

CONTAINS

SUBROUTINE tea_barrier

  INTEGER :: err

  CALL MPI_BARRIER(mpi_cart_comm,err)

END SUBROUTINE tea_barrier

SUBROUTINE tea_abort

  INTEGER :: ierr,err

  CALL MPI_ABORT(mpi_cart_comm,ierr,err)

END SUBROUTINE tea_abort

SUBROUTINE tea_finalize

  INTEGER :: err

  CLOSE(g_out)
  CALL FLUSH(0)
  CALL FLUSH(6)
  CALL FLUSH(g_out)
  CALL MPI_FINALIZE(err)

END SUBROUTINE tea_finalize

SUBROUTINE tea_init_comms

  IMPLICIT NONE

  INTEGER :: err,rank,size

  rank=0
  size=1

  CALL MPI_INIT(err)

  CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,err)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,size,err)

  parallel%parallel=.TRUE.
  parallel%task=rank

  IF(rank.EQ.0) THEN
    parallel%boss=.TRUE.
  ENDIF

  parallel%boss_task=0
  parallel%max_task=size

END SUBROUTINE tea_init_comms

SUBROUTINE tea_get_num_chunks(count)

  IMPLICIT NONE

  INTEGER :: count

! Should be changed so there can be more than one chunk per mpi task

  count=parallel%max_task

END SUBROUTINE tea_get_num_chunks

SUBROUTINE tea_decompose(x_cells,y_cells,z_cells,left,right,bottom,top,back,front)

  ! This decomposes the mesh into a number of chunks.
  ! The number of chunks may be a multiple of the number of mpi tasks
  ! Picks split with minimal surface area to volume ratio

  IMPLICIT NONE

  INTEGER :: x_cells,y_cells,z_cells
  INTEGER, dimension(1:) :: left,right,top,bottom,back,front
  INTEGER :: c,delta_x,delta_y,delta_z

  INTEGER  :: chunk_x,chunk_y,chunk_z,mod_x,mod_y,mod_z

  INTEGER  :: err

  ! 1 chunk
  c = 1

  ! Get destinations/sources
  CALL mpi_cart_shift(mpi_cart_comm, 0, 1,      &
    chunks(c)%chunk_neighbours(chunk_bottom),   &
    chunks(c)%chunk_neighbours(chunk_top),      &
    err)
  CALL mpi_cart_shift(mpi_cart_comm, 1, 1,      &
    chunks(c)%chunk_neighbours(chunk_left),     &
    chunks(c)%chunk_neighbours(chunk_right),    &
    err)
  CALL mpi_cart_shift(mpi_cart_comm, 2, 1,      &
    chunks(c)%chunk_neighbours(chunk_back),     &
    chunks(c)%chunk_neighbours(chunk_front),    &
    err)

  WHERE (chunks(c)%chunk_neighbours .EQ. mpi_proc_null)
    chunks(c)%chunk_neighbours = external_face
  END WHERE

  chunk_y = mpi_dims(1)
  chunk_x = mpi_dims(2)
  chunk_z = mpi_dims(3)

  delta_x=x_cells/chunk_x
  delta_y=y_cells/chunk_y
  delta_z=z_cells/chunk_z
  mod_x=MOD(x_cells,chunk_x)
  mod_y=MOD(y_cells,chunk_y)
  mod_z=MOD(z_cells,chunk_z)

  left(c) = mpi_coords(2)*delta_x + 1
  if (mpi_coords(2) .le. mod_x) then
    left(c) = left(c) + mpi_coords(2)
  else
    left(c) = left(c) + mod_x
  endif
  right(c) = left(c)+delta_x - 1
  if (mpi_coords(2) .lt. mod_x) then
    right(c) = right(c) + 1
  endif

  bottom(c) = mpi_coords(1)*delta_y + 1
  if (mpi_coords(1) .le. mod_y) then
    bottom(c) = bottom(c) + mpi_coords(1)
  else
    bottom(c) = bottom(c) + mod_y
  endif
  top(c) = bottom(c)+delta_y - 1
  if (mpi_coords(1) .lt. mod_y) then
    top(c) = top(c) + 1
  endif

  back(c) = mpi_coords(3)*delta_z + 1
  if (mpi_coords(3) .le. mod_z) then
    back(c) = back(c) + mpi_coords(3)
  else
    back(c) = back(c) + mod_z
  endif
  front(c) = back(c)+delta_z - 1
  if (mpi_coords(3) .lt. mod_z) then
    front(c) = front(c) + 1
  endif

  IF(parallel%boss)THEN
    WRITE(g_out,*)
    WRITE(g_out,*)"Decomposing the mesh into ",chunk_x," by ",chunk_y," by ",chunk_z," chunks"
    WRITE(g_out,*)
  ENDIF

END SUBROUTINE tea_decompose

SUBROUTINE tea_allocate_buffers(chunk)

  IMPLICIT NONE

  INTEGER      :: chunk

  INTEGER      :: bt_size, lr_size, fb_size
  INTEGER,parameter :: num_buffered=6

  INTEGER           :: allocate_extra_size

  allocate_extra_size = max(2, halo_allocate_depth)

  lr_size = num_buffered*(chunks(chunk)%field%y_max + 2*allocate_extra_size)*(chunks(chunk)%field%z_max + 2*allocate_extra_size)*halo_allocate_depth
  bt_size = num_buffered*(chunks(chunk)%field%x_max + 2*allocate_extra_size)*(chunks(chunk)%field%z_max + 2*allocate_extra_size)*halo_allocate_depth
  fb_size = num_buffered*(chunks(chunk)%field%x_max + 2*allocate_extra_size)*(chunks(chunk)%field%y_max + 2*allocate_extra_size)*halo_allocate_depth

  ! Unallocated buffers for external boundaries caused issues on some systems so they are now
  !  all allocated
  IF(parallel%task.EQ.chunks(chunk)%task)THEN
    !IF(chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) THEN
      ALLOCATE(chunks(chunk)%left_snd_buffer(lr_size))
      ALLOCATE(chunks(chunk)%left_rcv_buffer(lr_size))
    !ENDIF
    !IF(chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) THEN
      ALLOCATE(chunks(chunk)%right_snd_buffer(lr_size))
      ALLOCATE(chunks(chunk)%right_rcv_buffer(lr_size))
    !ENDIF
    !IF(chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) THEN
      ALLOCATE(chunks(chunk)%bottom_snd_buffer(bt_size))
      ALLOCATE(chunks(chunk)%bottom_rcv_buffer(bt_size))
    !ENDIF
    !IF(chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) THEN
      ALLOCATE(chunks(chunk)%top_snd_buffer(bt_size))
      ALLOCATE(chunks(chunk)%top_rcv_buffer(bt_size))
    !ENDIF
    !IF(chunks(chunk)%chunk_neighbours(chunk_back).NE.external_face) THEN
      ALLOCATE(chunks(chunk)%back_snd_buffer(fb_size))
      ALLOCATE(chunks(chunk)%back_rcv_buffer(fb_size))
    !ENDIF
    !IF(chunks(chunk)%chunk_neighbours(chunk_front).NE.external_face) THEN
      ALLOCATE(chunks(chunk)%front_snd_buffer(fb_size))
      ALLOCATE(chunks(chunk)%front_rcv_buffer(fb_size))
    !ENDIF
  ENDIF

  chunks(chunk)%left_snd_buffer = 0
  chunks(chunk)%right_snd_buffer = 0
  chunks(chunk)%bottom_snd_buffer = 0
  chunks(chunk)%top_snd_buffer = 0
  chunks(chunk)%back_snd_buffer = 0
  chunks(chunk)%front_snd_buffer = 0

  chunks(chunk)%left_rcv_buffer = 0
  chunks(chunk)%right_rcv_buffer = 0
  chunks(chunk)%bottom_rcv_buffer = 0
  chunks(chunk)%top_rcv_buffer = 0
  chunks(chunk)%back_rcv_buffer = 0
  chunks(chunk)%front_rcv_buffer = 0

END SUBROUTINE tea_allocate_buffers

SUBROUTINE tea_exchange(fields,depth)

  IMPLICIT NONE

    INTEGER      :: fields(NUM_FIELDS),depth, chunk,err,field
    INTEGER,DIMENSION(NUM_FIELDS)      :: left_right_offset,bottom_top_offset,back_front_offset
    INTEGER      :: end_pack_index_left_right, end_pack_index_bottom_top,end_pack_index_back_front
    INTEGER      :: message_count_lr, message_count_ud, message_count_fb
    INTEGER      :: exchange_size_lr, exchange_size_ud, exchange_size_fb
    INTEGER, dimension(4)                 :: request_lr, request_ud, request_fb
    INTEGER, dimension(MPI_STATUS_SIZE,4) :: status_lr, status_ud, status_fb
    LOGICAL :: test_complete

    ! Assuming 1 patch per task, this will be changed
    chunk = 1

    IF (ALL(chunks(chunk)%chunk_neighbours .eq. external_face)) return

    exchange_size_lr = depth*(chunks(chunk)%field%y_max+2*depth)*(chunks(chunk)%field%z_max+2*depth)
    exchange_size_ud = depth*(chunks(chunk)%field%x_max+2*depth)*(chunks(chunk)%field%z_max+2*depth)
    exchange_size_fb = depth*(chunks(chunk)%field%x_max+2*depth)*(chunks(chunk)%field%y_max+2*depth)

    request_lr = 0
    message_count_lr = 0
    request_ud = 0
    message_count_ud = 0
    request_fb = 0
    message_count_fb = 0

    end_pack_index_left_right=0
    end_pack_index_bottom_top=0
    end_pack_index_back_front=0
    left_right_offset = 0
    bottom_top_offset = 0
    back_front_offset = 0
    DO field=1,NUM_FIELDS
      IF(fields(field).EQ.1) THEN
        left_right_offset(field)=end_pack_index_left_right
        bottom_top_offset(field)=end_pack_index_bottom_top
        back_front_offset(field)=end_pack_index_back_front
        end_pack_index_left_right=end_pack_index_left_right + exchange_size_lr
        end_pack_index_bottom_top=end_pack_index_bottom_top + exchange_size_ud
        end_pack_index_back_front=end_pack_index_back_front + exchange_size_fb
      ENDIF
    ENDDO

    IF(chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) THEN
      ! do left exchanges
        CALL tea_pack_left(chunk, fields, depth, left_right_offset)

      !send and recv messagse to the left
      CALL tea_send_recv_message_left(chunks(chunk)%left_snd_buffer,                      &
                                         chunks(chunk)%left_rcv_buffer,                      &
                                         chunk,end_pack_index_left_right,                    &
                                         1, 2,                                               &
                                         request_lr(message_count_lr+1), request_lr(message_count_lr+2))
      message_count_lr = message_count_lr + 2
    ENDIF

    IF(chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) THEN
      ! do right exchanges
        CALL tea_pack_right(chunk, fields, depth, left_right_offset)

      !send message to the right
      CALL tea_send_recv_message_right(chunks(chunk)%right_snd_buffer,                     &
                                          chunks(chunk)%right_rcv_buffer,                     &
                                          chunk,end_pack_index_left_right,                    &
                                          2, 1,                                               &
                                          request_lr(message_count_lr+1), request_lr(message_count_lr+2))
      message_count_lr = message_count_lr + 2
    ENDIF

!    IF (depth .eq. 1) THEN
!      test_complete = .false.
!      ! don't have to transfer now
!      CALL MPI_TESTALL(message_count_lr, request_lr, test_complete, status_lr, err)
!    ELSE
!      test_complete = .true.
      !make a call to wait / sync
      CALL MPI_WAITALL(message_count_lr,request_lr,status_lr,err)
!    ENDIF

    !unpack in left direction
    IF(chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) THEN
      else
        CALL tea_unpack_left(fields, chunk, depth,                      &
                                chunks(chunk)%left_rcv_buffer,             &
                                left_right_offset)
    ENDIF

    !unpack in right direction
    IF(chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) THEN
        CALL tea_unpack_right(fields, chunk, depth,                     &
                                 chunks(chunk)%right_rcv_buffer,           &
                                 left_right_offset)
    ENDIF


    IF(chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) THEN
      ! do bottom exchanges
        CALL tea_pack_bottom(chunk, fields, depth, bottom_top_offset)

      !send message downwards
      CALL tea_send_recv_message_bottom(chunks(chunk)%bottom_snd_buffer,                     &
                                           chunks(chunk)%bottom_rcv_buffer,                     &
                                           chunk,end_pack_index_bottom_top,                     &
                                           3, 4,                                                &
                                           request_ud(message_count_ud+1), request_ud(message_count_ud+2))
      message_count_ud = message_count_ud + 2
    ENDIF

    IF(chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) THEN
      ! do top exchanges
        CALL tea_pack_top(chunk, fields, depth, bottom_top_offset)

      !send message upwards
      CALL tea_send_recv_message_top(chunks(chunk)%top_snd_buffer,                           &
                                        chunks(chunk)%top_rcv_buffer,                           &
                                        chunk,end_pack_index_bottom_top,                        &
                                        4, 3,                                                   &
                                        request_ud(message_count_ud+1), request_ud(message_count_ud+2))
      message_count_ud = message_count_ud + 2
    ENDIF

    !need to make a call to wait / sync
    CALL MPI_WAITALL(message_count_ud,request_ud,status_ud,err)

    !unpack in top direction
    IF( chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face ) THEN
        CALL tea_unpack_top(fields, chunk, depth,                       &
                               chunks(chunk)%top_rcv_buffer,               &
                               bottom_top_offset)
    ENDIF

    !unpack in bottom direction
    IF(chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) THEN
        CALL tea_unpack_bottom(fields, chunk, depth,                   &
                                 chunks(chunk)%bottom_rcv_buffer,         &
                                 bottom_top_offset)
    ENDIF

    IF(chunks(chunk)%chunk_neighbours(chunk_back).NE.external_face) THEN
      ! do back exchanges
        CALL tea_pack_back(chunk, fields, depth, back_front_offset)

      !send message downwards
      CALL tea_send_recv_message_back(chunks(chunk)%back_snd_buffer,                        &
                                           chunks(chunk)%back_rcv_buffer,                      &
                                           chunk,end_pack_index_back_front,                    &
                                           5, 6,                                               &
                                           request_fb(message_count_fb+1), request_fb(message_count_fb+2))
      message_count_fb = message_count_fb + 2
    ENDIF

    IF(chunks(chunk)%chunk_neighbours(chunk_front).NE.external_face) THEN
      ! do front exchanges
        CALL tea_pack_front(chunk, fields, depth, back_front_offset)

      !send message upwards
      CALL tea_send_recv_message_front(chunks(chunk)%front_snd_buffer,                       &
                                        chunks(chunk)%front_rcv_buffer,                         &
                                        chunk,end_pack_index_back_front,                        &
                                        6, 5,                                                   &
                                        request_fb(message_count_fb+1), request_fb(message_count_fb+2))
      message_count_fb = message_count_fb + 2
    ENDIF

    !need to make a call to wait / sync
    CALL MPI_WAITALL(message_count_fb,request_fb,status_fb,err)

    !unpack in front direction
    IF( chunks(chunk)%chunk_neighbours(chunk_front).NE.external_face ) THEN
        CALL tea_unpack_front(fields, chunk, depth,                       &
                               chunks(chunk)%front_rcv_buffer,               &
                               back_front_offset)
    ENDIF

    !unpack in back direction
    IF(chunks(chunk)%chunk_neighbours(chunk_back).NE.external_face) THEN
        CALL tea_unpack_back(fields, chunk, depth,                   &
                                 chunks(chunk)%back_rcv_buffer,         &
                                 back_front_offset)
    ENDIF

END SUBROUTINE tea_exchange

SUBROUTINE tea_send_recv_message_left(left_snd_buffer, left_rcv_buffer,      &
                                         chunk, total_size,                     &
                                         tag_send, tag_recv,                    &
                                         req_send, req_recv)

  REAL(KIND=8)    :: left_snd_buffer(:), left_rcv_buffer(:)
  INTEGER         :: left_task
  INTEGER         :: chunk
  INTEGER         :: total_size, tag_send, tag_recv, err
  INTEGER         :: req_send, req_recv

  left_task =chunks(chunk)%chunk_neighbours(chunk_left)

  CALL MPI_ISEND(left_snd_buffer,total_size,MPI_DOUBLE_PRECISION,left_task,tag_send &
                ,mpi_cart_comm,req_send,err)

  CALL MPI_IRECV(left_rcv_buffer,total_size,MPI_DOUBLE_PRECISION,left_task,tag_recv &
                ,mpi_cart_comm,req_recv,err)

END SUBROUTINE tea_send_recv_message_left

SUBROUTINE tea_send_recv_message_right(right_snd_buffer, right_rcv_buffer,   &
                                          chunk, total_size,                    &
                                          tag_send, tag_recv,                   &
                                          req_send, req_recv)

  IMPLICIT NONE

  REAL(KIND=8) :: right_snd_buffer(:), right_rcv_buffer(:)
  INTEGER      :: right_task
  INTEGER      :: chunk
  INTEGER      :: total_size, tag_send, tag_recv, err
  INTEGER      :: req_send, req_recv

  right_task=chunks(chunk)%chunk_neighbours(chunk_right)

  CALL MPI_ISEND(right_snd_buffer,total_size,MPI_DOUBLE_PRECISION,right_task,tag_send, &
                 mpi_cart_comm,req_send,err)

  CALL MPI_IRECV(right_rcv_buffer,total_size,MPI_DOUBLE_PRECISION,right_task,tag_recv, &
                 mpi_cart_comm,req_recv,err)

END SUBROUTINE tea_send_recv_message_right

SUBROUTINE tea_send_recv_message_top(top_snd_buffer, top_rcv_buffer,     &
                                        chunk, total_size,                  &
                                        tag_send, tag_recv,                 &
                                        req_send, req_recv)

    IMPLICIT NONE

    REAL(KIND=8) :: top_snd_buffer(:), top_rcv_buffer(:)
    INTEGER      :: top_task
    INTEGER      :: chunk
    INTEGER      :: total_size, tag_send, tag_recv, err
    INTEGER      :: req_send, req_recv

    top_task=chunks(chunk)%chunk_neighbours(chunk_top)

    CALL MPI_ISEND(top_snd_buffer,total_size,MPI_DOUBLE_PRECISION,top_task,tag_send, &
                   mpi_cart_comm,req_send,err)

    CALL MPI_IRECV(top_rcv_buffer,total_size,MPI_DOUBLE_PRECISION,top_task,tag_recv, &
                   mpi_cart_comm,req_recv,err)

END SUBROUTINE tea_send_recv_message_top

SUBROUTINE tea_send_recv_message_bottom(bottom_snd_buffer, bottom_rcv_buffer,        &
                                           chunk, total_size,                           &
                                           tag_send, tag_recv,                          &
                                           req_send, req_recv)

  IMPLICIT NONE

  REAL(KIND=8) :: bottom_snd_buffer(:), bottom_rcv_buffer(:)
  INTEGER      :: bottom_task
  INTEGER      :: chunk
  INTEGER      :: total_size, tag_send, tag_recv, err
  INTEGER      :: req_send, req_recv

  bottom_task=chunks(chunk)%chunk_neighbours(chunk_bottom)

  CALL MPI_ISEND(bottom_snd_buffer,total_size,MPI_DOUBLE_PRECISION,bottom_task,tag_send &
                ,mpi_cart_comm,req_send,err)

  CALL MPI_IRECV(bottom_rcv_buffer,total_size,MPI_DOUBLE_PRECISION,bottom_task,tag_recv &
                ,mpi_cart_comm,req_recv,err)

END SUBROUTINE tea_send_recv_message_bottom

SUBROUTINE tea_send_recv_message_back(back_snd_buffer, back_rcv_buffer,     &
                                         chunk, total_size,                  &
                                         tag_send, tag_recv,                 &
                                         req_send, req_recv)

  IMPLICIT NONE

  REAL(KIND=8) :: back_snd_buffer(:), back_rcv_buffer(:)
  INTEGER      :: back_task
  INTEGER      :: chunk
  INTEGER      :: total_size, tag_send, tag_recv, err
  INTEGER      :: req_send, req_recv

  back_task=chunks(chunk)%chunk_neighbours(chunk_back)

  CALL MPI_ISEND(back_snd_buffer,total_size,MPI_DOUBLE_PRECISION,back_task,tag_send, &
                 mpi_cart_comm,req_send,err)

  CALL MPI_IRECV(back_rcv_buffer,total_size,MPI_DOUBLE_PRECISION,back_task,tag_recv, &
                 mpi_cart_comm,req_recv,err)

END SUBROUTINE tea_send_recv_message_back

SUBROUTINE tea_send_recv_message_front(front_snd_buffer, front_rcv_buffer,        &
                                          chunk, total_size,                         &
                                          tag_send, tag_recv,                        &
                                          req_send, req_recv)

  IMPLICIT NONE

  REAL(KIND=8) :: front_snd_buffer(:), front_rcv_buffer(:)
  INTEGER      :: front_task
  INTEGER      :: chunk
  INTEGER      :: total_size, tag_send, tag_recv, err
  INTEGER      :: req_send, req_recv

  front_task=chunks(chunk)%chunk_neighbours(chunk_front)

  CALL MPI_ISEND(front_snd_buffer,total_size,MPI_DOUBLE_PRECISION,front_task,tag_send &
                ,mpi_cart_comm,req_send,err)

  CALL MPI_IRECV(front_rcv_buffer,total_size,MPI_DOUBLE_PRECISION,front_task,tag_recv &
                ,mpi_cart_comm,req_recv,err)

END SUBROUTINE tea_send_recv_message_front

SUBROUTINE tea_sum(value)

  ! Only sums to the master

  IMPLICIT NONE

  REAL(KIND=8) :: value

  REAL(KIND=8) :: total

  INTEGER :: err

  total=value

  CALL MPI_REDUCE(value,total,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,mpi_cart_comm,err)

  value=total

END SUBROUTINE tea_sum

SUBROUTINE tea_allsum(value)

  ! Global reduction for CG solver

  IMPLICIT NONE

  REAL(KIND=8) :: value

  REAL(KIND=8) :: total

  INTEGER :: err

  total=value

  CALL MPI_ALLREDUCE(value,total,1,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_cart_comm,err)

  value=total

END SUBROUTINE tea_allsum

SUBROUTINE tea_min(value)

  IMPLICIT NONE

  REAL(KIND=8) :: value

  REAL(KIND=8) :: minimum

  INTEGER :: err

  minimum=value

  CALL MPI_ALLREDUCE(value,minimum,1,MPI_DOUBLE_PRECISION,MPI_MIN,mpi_cart_comm,err)

  value=minimum

END SUBROUTINE tea_min

SUBROUTINE tea_max(value)

  IMPLICIT NONE

  REAL(KIND=8) :: value

  REAL(KIND=8) :: maximum

  INTEGER :: err

  maximum=value

  CALL MPI_ALLREDUCE(value,maximum,1,MPI_DOUBLE_PRECISION,MPI_MAX,mpi_cart_comm,err)

  value=maximum

END SUBROUTINE tea_max

SUBROUTINE tea_allgather(value,values)

  IMPLICIT NONE

  REAL(KIND=8) :: value

  REAL(KIND=8) :: values(parallel%max_task)

  INTEGER :: err

  values(1)=value ! Just to ensure it will work in serial

  CALL MPI_ALLGATHER(value,1,MPI_DOUBLE_PRECISION,values,1,MPI_DOUBLE_PRECISION,mpi_cart_comm,err)

END SUBROUTINE tea_allgather

SUBROUTINE tea_check_error(error)

  IMPLICIT NONE

  INTEGER :: error

  INTEGER :: maximum

  INTEGER :: err

  maximum=error

  CALL MPI_ALLREDUCE(error,maximum,1,MPI_INTEGER,MPI_MAX,mpi_cart_comm,err)

  error=maximum

END SUBROUTINE tea_check_error


END MODULE tea_module
