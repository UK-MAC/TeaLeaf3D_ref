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

!>  @brief Fortran mpi buffer packing kernel
!>  @author Wayne Gaudin
!>  @details Packs/unpacks mpi send and receive buffers

MODULE pack_kernel_module

CONTAINS

SUBROUTINE tea_pack_message_left(x_min,x_max,y_min,y_max,z_min,z_max,halo_exchange_depth,field,                 &
                                left_snd_buffer,                                           &
                                depth,x_inc, y_inc, z_inc,   &
                                buffer_offset)

  IMPLICIT NONE

  INTEGER      :: depth,x_min,x_max,y_min,y_max,z_min,z_max,halo_exchange_depth
  INTEGER      :: j,k,l,x_inc,y_inc,z_inc,index,buffer_offset

  REAL(KIND=8) :: field(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth,z_min-halo_exchange_depth:z_max+halo_exchange_depth)
  REAL(KIND=8) :: left_snd_buffer(:)

  ! Pack

!$OMP DO
  DO l=z_min-depth,z_max+z_inc+depth
    DO k=y_min-depth,y_max+y_inc+depth
      DO j=1,depth
        index=buffer_offset + &
            j + &
            (k+depth-1)*depth + &
            (l+depth-1)*(y_max+y_inc+2*depth)*depth
        left_snd_buffer(index)=field(x_min+x_inc-1+j,k,l)
      ENDDO
    ENDDO
  ENDDO
!$OMP END DO NOWAIT

END SUBROUTINE tea_pack_message_left

SUBROUTINE tea_unpack_message_left(x_min,x_max,y_min,y_max,z_min,z_max,halo_exchange_depth,field,                 &
                                left_rcv_buffer,                                           &
                                depth,x_inc, y_inc, z_inc,   &
                                buffer_offset)

  IMPLICIT NONE

  INTEGER      :: depth,x_min,x_max,y_min,y_max,z_min,z_max,halo_exchange_depth
  INTEGER      :: j,k,l,x_inc,y_inc,z_inc,index,buffer_offset

  REAL(KIND=8) :: field(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth,z_min-halo_exchange_depth:z_max+halo_exchange_depth)
  REAL(KIND=8) :: left_rcv_buffer(:)

  ! Unpack

!$OMP DO
  DO l=z_min-depth,z_max+z_inc+depth
    DO k=y_min-depth,y_max+y_inc+depth
      DO j=1,depth
        index=buffer_offset + (j+(k+depth-1)*depth) + ((l+depth-1)*(y_max+y_inc+2*depth)*depth)
        field(x_min-j,k,l)=left_rcv_buffer(index)
      ENDDO
    ENDDO
  ENDDO
!$OMP END DO NOWAIT

END SUBROUTINE tea_unpack_message_left

SUBROUTINE tea_pack_message_right(x_min,x_max,y_min,y_max,z_min,z_max,halo_exchange_depth,field,                 &
                                right_snd_buffer,                                           &
                                depth,x_inc, y_inc, z_inc,   &
                                buffer_offset)

  IMPLICIT NONE

  INTEGER      :: depth,x_min,x_max,y_min,y_max,z_min,z_max,halo_exchange_depth
  INTEGER      :: j,k,l,x_inc,y_inc,z_inc,index,buffer_offset

  REAL(KIND=8) :: field(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth,z_min-halo_exchange_depth:z_max+halo_exchange_depth)
  REAL(KIND=8) :: right_snd_buffer(:)

  ! Pack

!$OMP DO
  DO l=z_min-depth,z_max+z_inc+depth
    DO k=y_min-depth,y_max+y_inc+depth
      DO j=1,depth
        index=buffer_offset + (j+(k+depth-1)*depth) + ((l+depth-1)*(y_max+y_inc+2*depth)*depth)
        right_snd_buffer(index)=field(x_max+1-j,k,l)
      ENDDO
    ENDDO
  ENDDO
!$OMP END DO NOWAIT

END SUBROUTINE tea_pack_message_right

SUBROUTINE tea_unpack_message_right(x_min,x_max,y_min,y_max,z_min,z_max,halo_exchange_depth,field,                 &
                                right_rcv_buffer,                                           &
                                depth,x_inc, y_inc, z_inc,   &
                                buffer_offset)

  IMPLICIT NONE

  INTEGER      :: depth,x_min,x_max,y_min,y_max,z_min,z_max,halo_exchange_depth
  INTEGER      :: j,k,l,x_inc,y_inc,z_inc,index,buffer_offset

  REAL(KIND=8) :: field(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth,z_min-halo_exchange_depth:z_max+halo_exchange_depth)
  REAL(KIND=8) :: right_rcv_buffer(:)

  ! Unpack

!$OMP DO
  DO l=z_min-depth,z_max+z_inc+depth
    DO k=y_min-depth,y_max+y_inc+depth
      DO j=1,depth
        index=buffer_offset + (j+(k+depth-1)*depth) + ((l+depth-1)*(y_max+y_inc+2*depth)*depth)
        field(x_max+x_inc+j,k,l)=right_rcv_buffer(index)
      ENDDO
    ENDDO
  ENDDO
!$OMP END DO NOWAIT

END SUBROUTINE tea_unpack_message_right

SUBROUTINE tea_pack_message_top(x_min,x_max,y_min,y_max,z_min,z_max,halo_exchange_depth,field,                 &
                                top_snd_buffer,                                           &
                                depth,x_inc, y_inc, z_inc,   &
                                buffer_offset)

  IMPLICIT NONE

  INTEGER      :: depth,x_min,x_max,y_min,y_max,z_min,z_max,halo_exchange_depth
  INTEGER      :: j,k,l,x_inc,y_inc,z_inc,index,buffer_offset

  REAL(KIND=8) :: field(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth,z_min-halo_exchange_depth:z_max+halo_exchange_depth)
  REAL(KIND=8) :: top_snd_buffer(:)

  ! Pack

  DO k=1,depth
!$OMP DO
    DO l=z_min-depth,z_max+z_inc+depth
      DO j=x_min-depth,x_max+x_inc+depth
        index= buffer_offset + &
            (j+depth) + &
            (l-1+depth)*(x_max+x_inc+2*depth) + &
            (k-1)*((x_max+x_inc+2*depth)*(z_max+z_inc+2*depth))
        top_snd_buffer(index)=field(j,y_max+1-k,l)
      ENDDO
    ENDDO
!$OMP END DO NOWAIT
  ENDDO

END SUBROUTINE tea_pack_message_top

SUBROUTINE tea_unpack_message_top(x_min,x_max,y_min,y_max,z_min,z_max,halo_exchange_depth,field,                 &
                                top_rcv_buffer,                                           &
                                depth,x_inc, y_inc, z_inc,   &
                                buffer_offset)

  IMPLICIT NONE

  INTEGER      :: depth,x_min,x_max,y_min,y_max,z_min,z_max,halo_exchange_depth
  INTEGER      :: j,k,l,x_inc,y_inc,z_inc,index,buffer_offset

  REAL(KIND=8) :: field(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth,z_min-halo_exchange_depth:z_max+halo_exchange_depth)
  REAL(KIND=8) :: top_rcv_buffer(:)

  ! Unpack

  DO k=1,depth
!$OMP DO
    DO l=z_min-depth,z_max+z_inc+depth
      DO j=x_min-depth,x_max+x_inc+depth
        index= buffer_offset + (j+depth) + (l-1+depth)*(x_max+x_inc+2*depth) + (k-1)*((x_max+x_inc+2*depth)*(z_max+z_inc+2*depth))
        field(j,y_max+y_inc+k,l)=top_rcv_buffer(index)
      ENDDO
    ENDDO
!$OMP END DO NOWAIT
  ENDDO

END SUBROUTINE tea_unpack_message_top

SUBROUTINE tea_pack_message_bottom(x_min,x_max,y_min,y_max,z_min,z_max,halo_exchange_depth,field,                 &
                                bottom_snd_buffer,                                           &
                                depth,x_inc, y_inc, z_inc,   &
                                buffer_offset)

  IMPLICIT NONE

  INTEGER      :: depth,x_min,x_max,y_min,y_max,z_min,z_max,halo_exchange_depth
  INTEGER      :: j,k,l,x_inc,y_inc,z_inc,index,buffer_offset

  REAL(KIND=8) :: field(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth,z_min-halo_exchange_depth:z_max+halo_exchange_depth)
  REAL(KIND=8) :: bottom_snd_buffer(:)

  ! Pack

  DO k=1,depth
!$OMP DO
    DO l=z_min-depth,z_max+z_inc+depth
      DO j=x_min-depth,x_max+x_inc+depth
        index= buffer_offset + (j+depth) + (l-1+depth)*(x_max+x_inc+2*depth) + (k-1)*((x_max+x_inc+2*depth)*(z_max+z_inc+2*depth))
        bottom_snd_buffer(index)=field(j,y_min+y_inc-1+k,l)
      ENDDO
    ENDDO
!$OMP END DO NOWAIT
  ENDDO

END SUBROUTINE tea_pack_message_bottom

SUBROUTINE tea_unpack_message_bottom(x_min,x_max,y_min,y_max,z_min,z_max,halo_exchange_depth,field,                 &
                                bottom_rcv_buffer,                                           &
                                depth,x_inc, y_inc, z_inc,   &
                                buffer_offset)

  IMPLICIT NONE

  INTEGER      :: depth,x_min,x_max,y_min,y_max,z_min,z_max,halo_exchange_depth
  INTEGER      :: j,k,l,x_inc,y_inc,z_inc,index,buffer_offset

  REAL(KIND=8) :: field(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth,z_min-halo_exchange_depth:z_max+halo_exchange_depth)
  REAL(KIND=8) :: bottom_rcv_buffer(:)

  ! Unpack

  DO k=1,depth
!$OMP DO
    DO l=z_min-depth,z_max+z_inc+depth
      DO j=x_min-depth,x_max+x_inc+depth
        index= buffer_offset + (j+depth) + (l-1+depth)*(x_max+x_inc+2*depth) + (k-1)*((x_max+x_inc+2*depth)*(z_max+z_inc+2*depth))
        field(j,y_min-k,l)=bottom_rcv_buffer(index)
      ENDDO
    ENDDO
!$OMP END DO NOWAIT
  ENDDO

END SUBROUTINE tea_unpack_message_bottom

SUBROUTINE tea_pack_message_back(x_min,x_max,y_min,y_max,z_min,z_max,halo_exchange_depth,field,                 &
                                back_snd_buffer,                                           &
                                depth,x_inc, y_inc, z_inc,   &
                                buffer_offset)

  IMPLICIT NONE

  INTEGER      :: depth,x_min,x_max,y_min,y_max,z_min,z_max,halo_exchange_depth
  INTEGER      :: j,k,l,x_inc,y_inc,z_inc,index,buffer_offset

  REAL(KIND=8) :: field(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth,z_min-halo_exchange_depth:z_max+halo_exchange_depth)
  REAL(KIND=8) :: back_snd_buffer(:)

  ! Pack

  DO l=1,depth
!$OMP DO
    DO k=y_min-depth,y_max+y_inc+depth
      DO j=x_min-depth,x_max+x_inc+depth
        index= buffer_offset + &
            (j+depth) + &
            (k-1+depth)*(x_max+x_inc+2*depth) + &
            (l-1)*((x_max+x_inc+2*depth)*(y_max+y_inc+2*depth))
        back_snd_buffer(index)=field(j,k,z_min+z_inc-1+l)
      ENDDO
    ENDDO
!$OMP END DO NOWAIT
  ENDDO

END SUBROUTINE tea_pack_message_back

SUBROUTINE tea_unpack_message_back(x_min,x_max,y_min,y_max,z_min,z_max,halo_exchange_depth,field,                 &
                                back_rcv_buffer,                                           &
                                depth,x_inc, y_inc, z_inc,   &
                                buffer_offset)

  IMPLICIT NONE

  INTEGER      :: depth,x_min,x_max,y_min,y_max,z_min,z_max,halo_exchange_depth
  INTEGER      :: j,k,l,x_inc,y_inc,z_inc,index,buffer_offset

  REAL(KIND=8) :: field(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth,z_min-halo_exchange_depth:z_max+halo_exchange_depth)
  REAL(KIND=8) :: back_rcv_buffer(:)

  ! Unpack

  DO l=1,depth
!$OMP DO
    DO k=y_min-depth,y_max+y_inc+depth
      DO j=x_min-depth,x_max+x_inc+depth
        index= buffer_offset + (j+depth) + (k-1+depth)*(x_max+x_inc+2*depth) + (l-1)*((x_max+x_inc+2*depth)*(y_max+y_inc+2*depth))
        field(j,k,z_min-l)=back_rcv_buffer(index)
      ENDDO
    ENDDO
!$OMP END DO NOWAIT
  ENDDO

END SUBROUTINE tea_unpack_message_back

SUBROUTINE tea_pack_message_front(x_min,x_max,y_min,y_max,z_min,z_max,halo_exchange_depth,field,                 &
                                front_snd_buffer,                                           &
                                depth,x_inc, y_inc, z_inc,   &
                                buffer_offset)

  IMPLICIT NONE

  INTEGER      :: depth,x_min,x_max,y_min,y_max,z_min,z_max,halo_exchange_depth
  INTEGER      :: j,k,l,x_inc,y_inc,z_inc,index,buffer_offset

  REAL(KIND=8) :: field(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth,z_min-halo_exchange_depth:z_max+halo_exchange_depth)
  REAL(KIND=8) :: front_snd_buffer(:)

  ! Pack

  DO l=1,depth
!$OMP DO
    DO k=y_min-depth,y_max+y_inc+depth
      DO j=x_min-depth,x_max+x_inc+depth
        index= buffer_offset + (j+depth) + (k-1+depth)*(x_max+x_inc+2*depth) + (l-1)*((x_max+x_inc+2*depth)*(y_max+y_inc+2*depth))
        front_snd_buffer(index)=field(j,k,z_max+1-l)
      ENDDO
    ENDDO
!$OMP END DO NOWAIT
  ENDDO

END SUBROUTINE tea_pack_message_front

SUBROUTINE tea_unpack_message_front(x_min,x_max,y_min,y_max,z_min,z_max,halo_exchange_depth,field,                 &
                                front_rcv_buffer,                                           &
                                depth,x_inc, y_inc, z_inc,   &
                                buffer_offset)

  IMPLICIT NONE

  INTEGER      :: depth,x_min,x_max,y_min,y_max,z_min,z_max,halo_exchange_depth
  INTEGER      :: j,k,l,x_inc,y_inc,z_inc,index,buffer_offset

  REAL(KIND=8) :: field(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth,z_min-halo_exchange_depth:z_max+halo_exchange_depth)
  REAL(KIND=8) :: front_rcv_buffer(:)

  ! Unpack

  DO l=1,depth
!$OMP DO
    DO k=y_min-depth,y_max+y_inc+depth
      DO j=x_min-depth,x_max+x_inc+depth
        index= buffer_offset + (j+depth) + (k-1+depth)*(x_max+x_inc+2*depth) + ((l-1)*(x_max+x_inc+2*depth)*(y_max+y_inc+2*depth))
        field(j,k,z_max+z_inc+l)=front_rcv_buffer(index)
      ENDDO
    ENDDO
!$OMP END DO NOWAIT
  ENDDO

END SUBROUTINE tea_unpack_message_front

END MODULE pack_kernel_module
