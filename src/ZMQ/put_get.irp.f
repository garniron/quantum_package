integer function zmq_put_dvector(zmq_to_qp_run_socket, worker_id, name, x, size_x)
  use f77_zmq
  implicit none
  BEGIN_DOC
! Put a float vector on the qp_run scheduler
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer, intent(in)            :: worker_id
  character*(*)                  :: name
  integer, intent(in)            :: size_x
  double precision, intent(in)   :: x(size_x)
  integer                        :: rc
  character*(256)                :: msg

  zmq_put_dvector = 0

  write(msg,'(A,1X,I8,1X,A200)') 'put_data '//trim(zmq_state), worker_id, name
  rc = f77_zmq_send(zmq_to_qp_run_socket,trim(msg),len(trim(msg)),ZMQ_SNDMORE)
  if (rc /= len(trim(msg))) then
    zmq_put_dvector = -1
    return 
  endif

  rc = f77_zmq_send(zmq_to_qp_run_socket,x,size_x*8,0)
  if (rc /= size_x*8) then
    zmq_put_dvector = -1
    return
  endif

  rc = f77_zmq_recv(zmq_to_qp_run_socket,msg,len(msg),0)
  if (msg(1:rc) /= 'put_data_reply ok') then
    zmq_put_dvector = -1
    return
  endif

end


integer function zmq_get_dvector(zmq_to_qp_run_socket, worker_id, name, x, size_x)
  use f77_zmq
  implicit none
  BEGIN_DOC
! Get a float vector from the qp_run scheduler
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer, intent(in)            :: worker_id
  integer, intent(in)            :: size_x
  character*(*), intent(in)      :: name
  double precision, intent(out)  :: x(size_x)
  integer                        :: rc
  character*(256)                :: msg

  PROVIDE zmq_state
  ! Success
  zmq_get_dvector = 0

  if (mpi_master) then
    write(msg,'(A,1X,I8,1X,A200)') 'get_data '//trim(zmq_state), worker_id, name
    rc = f77_zmq_send(zmq_to_qp_run_socket,trim(msg),len(trim(msg)),0)
    if (rc /= len(trim(msg))) then
      zmq_get_dvector = -1
      print *,  irp_here, 'rc /= len(trim(msg))', rc, len(trim(msg))
      go to 10
    endif

    rc = f77_zmq_recv(zmq_to_qp_run_socket,msg,len(msg),0)
    if (msg(1:14) /= 'get_data_reply') then
      print *,  irp_here, 'msg(1:14) /= get_data_reply', msg(1:14)
      zmq_get_dvector = -1
      go to 10
    endif

    rc = f77_zmq_recv(zmq_to_qp_run_socket,x,size_x*8,0)
    if (rc /= size_x*8) then
      print *,  irp_here, 'rc /= size_x*8', rc, size_x*8
      zmq_get_dvector = -1
      go to 10
    endif
  endif

  10 continue

  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  IRP_IF MPI
    integer :: ierr
    include 'mpif.h'
    call MPI_BCAST (zmq_get_dvector, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      print *,  irp_here//': Unable to broadcast zmq_get_dvector'
      stop -1
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call broadcast_chunks_double(x, size_x)
  IRP_ENDIF

end



integer function zmq_put_ivector(zmq_to_qp_run_socket, worker_id, name, x, size_x)
  use f77_zmq
  implicit none
  BEGIN_DOC
! Put a vector of integers on the qp_run scheduler
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer, intent(in)            :: worker_id
  character*(*)                  :: name
  integer, intent(in)            :: size_x
  integer, intent(in)            :: x(size_x)
  integer                        :: rc
  character*(256)                :: msg

  zmq_put_ivector = 0

  write(msg,'(A,1X,I8,1X,A200)') 'put_data '//trim(zmq_state), worker_id, name
  rc = f77_zmq_send(zmq_to_qp_run_socket,trim(msg),len(trim(msg)),ZMQ_SNDMORE)
  if (rc /= len(trim(msg))) then
    zmq_put_ivector = -1
    return 
  endif

  rc = f77_zmq_send(zmq_to_qp_run_socket,x,size_x*4,0)
  if (rc /= size_x*4) then
    zmq_put_ivector = -1
    return
  endif

  rc = f77_zmq_recv(zmq_to_qp_run_socket,msg,len(msg),0)
  if (msg(1:rc) /= 'put_data_reply ok') then
    zmq_put_ivector = -1
    return
  endif

end


integer function zmq_get_ivector(zmq_to_qp_run_socket, worker_id, name, x, size_x)
  use f77_zmq
  implicit none
  BEGIN_DOC
! Get a vector of integers from the qp_run scheduler
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer, intent(in)            :: worker_id
  integer, intent(in)            :: size_x
  character*(*), intent(in)      :: name
  integer, intent(out)           :: x(size_x)
  integer                        :: rc
  character*(256)                :: msg

  PROVIDE zmq_state
  ! Success
  zmq_get_ivector = 0

  if (mpi_master) then
    write(msg,'(A,1X,I8,1X,A200)') 'get_data '//trim(zmq_state), worker_id, name
    rc = f77_zmq_send(zmq_to_qp_run_socket,trim(msg),len(trim(msg)),0)
    if (rc /= len(trim(msg))) then
      zmq_get_ivector = -1
      go to 10
    endif

    rc = f77_zmq_recv(zmq_to_qp_run_socket,msg,len(msg),0)
    if (msg(1:14) /= 'get_data_reply') then
      zmq_get_ivector = -1
      go to 10
    endif

    rc = f77_zmq_recv(zmq_to_qp_run_socket,x,size_x*4,0)
    if (rc /= size_x*4) then
      zmq_get_ivector = -1
      go to 10
    endif
  endif

  10 continue

  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  IRP_IF MPI
    integer :: ierr
    include 'mpif.h'
    call MPI_BCAST (zmq_get_ivector, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      print *,  irp_here//': Unable to broadcast zmq_get_ivector'
      stop -1
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call broadcast_chunks_integer(x, size_x)
  IRP_ENDIF

end


integer function zmq_put8_dvector(zmq_to_qp_run_socket, worker_id, name, x, size_x)
  use f77_zmq
  implicit none
  BEGIN_DOC
! Put a float vector on the qp_run scheduler
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer, intent(in)            :: worker_id
  character*(*)                  :: name
  integer*8, intent(in)          :: size_x
  double precision, intent(in)   :: x(size_x)
  integer*8                      :: rc
  character*(256)                :: msg

  zmq_put8_dvector = 0

  write(msg,'(A,1X,I8,1X,A200)') 'put_data '//trim(zmq_state), worker_id, name
  rc = f77_zmq_send(zmq_to_qp_run_socket,trim(msg),len(trim(msg)),ZMQ_SNDMORE)
  if (rc /= len(trim(msg))) then
    zmq_put8_dvector = -1
    return 
  endif

  rc = f77_zmq_send(zmq_to_qp_run_socket,x,size_x*8,0)
  if (rc /= size_x*8) then
    zmq_put8_dvector = -1
    return
  endif

  rc = f77_zmq_recv(zmq_to_qp_run_socket,msg,len(msg),0)
  if (msg(1:rc) /= 'put_data_reply ok') then
    zmq_put8_dvector = -1
    return
  endif

end


integer function zmq_get8_dvector(zmq_to_qp_run_socket, worker_id, name, x, size_x)
  use f77_zmq
  implicit none
  BEGIN_DOC
! Get a float vector from the qp_run scheduler
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer, intent(in)            :: worker_id
  integer*8, intent(in)          :: size_x
  character*(*), intent(in)      :: name
  double precision, intent(out)  :: x(size_x)
  integer*8                      :: rc
  character*(256)                :: msg

  PROVIDE zmq_state
  ! Success
  zmq_get8_dvector = 0

  if (mpi_master) then
    write(msg,'(A,1X,I8,1X,A200)') 'get_data '//trim(zmq_state), worker_id, name
    rc = f77_zmq_send(zmq_to_qp_run_socket,trim(msg),len(trim(msg)),0)
    if (rc /= len(trim(msg))) then
      zmq_get8_dvector = -1
      print *,  irp_here, 'rc /= len(trim(msg))', rc, len(trim(msg))
      go to 10
    endif

    rc = f77_zmq_recv(zmq_to_qp_run_socket,msg,len(msg),0)
    if (msg(1:14) /= 'get_data_reply') then
      print *,  irp_here, 'msg(1:14) /= get_data_reply', msg(1:14)
      zmq_get8_dvector = -1
      go to 10
    endif

    rc = f77_zmq_recv(zmq_to_qp_run_socket,x,size_x*8,0)
    if (rc /= size_x*8) then
      print *,  irp_here, 'rc /= size_x*8', rc, size_x*8
      zmq_get8_dvector = -1
      go to 10
    endif
  endif

  10 continue

  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  IRP_IF MPI
    integer :: ierr
    include 'mpif.h'
    call MPI_BCAST (zmq_get8_dvector, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      print *,  irp_here//': Unable to broadcast zmq_get8_dvector'
      stop -1
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call broadcast_chunks_double(x, size_x)
  IRP_ENDIF

end



integer function zmq_put8_ivector(zmq_to_qp_run_socket, worker_id, name, x, size_x)
  use f77_zmq
  implicit none
  BEGIN_DOC
! Put a vector of integers on the qp_run scheduler
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer, intent(in)            :: worker_id
  character*(*)                  :: name
  integer*8, intent(in)          :: size_x
  integer, intent(in)            :: x(size_x)
  integer*8                      :: rc
  character*(256)                :: msg

  zmq_put8_ivector = 0

  write(msg,'(A,1X,I8,1X,A200)') 'put_data '//trim(zmq_state), worker_id, name
  rc = f77_zmq_send(zmq_to_qp_run_socket,trim(msg),len(trim(msg)),ZMQ_SNDMORE)
  if (rc /= len(trim(msg))) then
    zmq_put8_ivector = -1
    return 
  endif

  rc = f77_zmq_send(zmq_to_qp_run_socket,x,size_x*4,0)
  if (rc /= size_x*4) then
    zmq_put8_ivector = -1
    return
  endif

  rc = f77_zmq_recv(zmq_to_qp_run_socket,msg,len(msg),0)
  if (msg(1:rc) /= 'put_data_reply ok') then
    zmq_put8_ivector = -1
    return
  endif

end


integer function zmq_get8_ivector(zmq_to_qp_run_socket, worker_id, name, x, size_x)
  use f77_zmq
  implicit none
  BEGIN_DOC
! Get a vector of integers from the qp_run scheduler
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer, intent(in)            :: worker_id
  integer*8, intent(in)            :: size_x
  character*(*), intent(in)      :: name
  integer, intent(out)           :: x(size_x)
  integer*8                      :: rc
  character*(256)                :: msg

  PROVIDE zmq_state
  ! Success
  zmq_get8_ivector = 0

  if (mpi_master) then
    write(msg,'(A,1X,I8,1X,A200)') 'get_data '//trim(zmq_state), worker_id, name
    rc = f77_zmq_send(zmq_to_qp_run_socket,trim(msg),len(trim(msg)),0)
    if (rc /= len(trim(msg))) then
      zmq_get8_ivector = -1
      go to 10
    endif

    rc = f77_zmq_recv(zmq_to_qp_run_socket,msg,len(msg),0)
    if (msg(1:14) /= 'get_data_reply') then
      zmq_get8_ivector = -1
      go to 10
    endif

    rc = f77_zmq_recv(zmq_to_qp_run_socket,x,size_x*4,0)
    if (rc /= size_x*4) then
      zmq_get8_ivector = -1
      go to 10
    endif
  endif

  10 continue

  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  IRP_IF MPI
    integer :: ierr
    include 'mpif.h'
    call MPI_BCAST (zmq_get8_ivector, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      print *,  irp_here//': Unable to broadcast zmq_get8_ivector'
      stop -1
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call broadcast_chunks_integer(x, size_x)
  IRP_ENDIF

end



integer function zmq_put_int(zmq_to_qp_run_socket, worker_id, name, x)
  use f77_zmq
  implicit none
  BEGIN_DOC
! Put a vector of integers on the qp_run scheduler
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer, intent(in)            :: worker_id
  character*(*)                  :: name
  integer, intent(in)            :: x
  integer                        :: rc
  character*(256)                :: msg

  zmq_put_int = 0

  write(msg,'(A,1X,I8,1X,A200)') 'put_data '//trim(zmq_state), worker_id, name
  rc = f77_zmq_send(zmq_to_qp_run_socket,trim(msg),len(trim(msg)),ZMQ_SNDMORE)
  if (rc /= len(trim(msg))) then
    zmq_put_int = -1
    return 
  endif

  rc = f77_zmq_send(zmq_to_qp_run_socket,x,4,0)
  if (rc /= 4) then
    zmq_put_int = -1
    return
  endif

  rc = f77_zmq_recv(zmq_to_qp_run_socket,msg,len(msg),0)
  if (msg(1:rc) /= 'put_data_reply ok') then
    zmq_put_int = -1
    return
  endif

end

integer function zmq_get_int(zmq_to_qp_run_socket, worker_id, name, x)
  use f77_zmq
  implicit none
  BEGIN_DOC
! Get a vector of integers from the qp_run scheduler
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer, intent(in)            :: worker_id
  character*(*), intent(in)      :: name
  integer, intent(out)           :: x
  integer                        :: rc
  character*(256)                :: msg

  PROVIDE zmq_state
  ! Success
  zmq_get_int = 0

  if (mpi_master) then
    write(msg,'(A,1X,I8,1X,A200)') 'get_data '//trim(zmq_state), worker_id, name
    rc = f77_zmq_send(zmq_to_qp_run_socket,trim(msg),len(trim(msg)),0)
    if (rc /= len(trim(msg))) then
      zmq_get_int = -1
      go to 10
    endif

    rc = f77_zmq_recv(zmq_to_qp_run_socket,msg,len(msg),0)
    if (msg(1:14) /= 'get_data_reply') then
      zmq_get_int = -1
      go to 10
    endif

    rc = f77_zmq_recv(zmq_to_qp_run_socket,x,4,0)
    if (rc /= 4) then
      zmq_get_int = -1
      go to 10
    endif
  endif

  10 continue

end

