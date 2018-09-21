integer function zmq_put_psi(zmq_to_qp_run_socket,worker_id)
  use f77_zmq
  implicit none
  BEGIN_DOC
! Put the wave function on the qp_run scheduler
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer, intent(in)            :: worker_id
  character*(256)                :: msg

  integer, external              :: zmq_put_N_states
  integer, external              :: zmq_put_N_det
  integer, external              :: zmq_put_psi_det_size
  integer, external              :: zmq_put_psi_det
  integer, external              :: zmq_put_psi_coef

  zmq_put_psi = 0
  if (zmq_put_N_states(zmq_to_qp_run_socket, worker_id) == -1) then
    zmq_put_psi = -1
    return
  endif
  if (zmq_put_N_det(zmq_to_qp_run_socket, worker_id) == -1) then
    zmq_put_psi = -1
    return
  endif
  if (zmq_put_psi_det_size(zmq_to_qp_run_socket, worker_id) == -1) then
    zmq_put_psi = -1
    return
  endif
  if (zmq_put_psi_det(zmq_to_qp_run_socket, worker_id) == -1) then
    zmq_put_psi = -1
    return
  endif
  if (zmq_put_psi_coef(zmq_to_qp_run_socket, worker_id) == -1) then
    zmq_put_psi = -1
    return
  endif

end



BEGIN_TEMPLATE 

integer function zmq_put_$X(zmq_to_qp_run_socket,worker_id)
  use f77_zmq
  implicit none
  BEGIN_DOC
! Put $X on the qp_run scheduler
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer, intent(in)            :: worker_id
  integer                        :: rc
  character*(256)                :: msg

  zmq_put_$X = 0

  write(msg,'(A,1X,I8,1X,A200)') 'put_data '//trim(zmq_state), worker_id, '$X'
  rc = f77_zmq_send(zmq_to_qp_run_socket,trim(msg),len(trim(msg)),ZMQ_SNDMORE)
  if (rc /= len(trim(msg))) then
    zmq_put_$X = -1
    return
  endif

  rc = f77_zmq_send(zmq_to_qp_run_socket,$X,4,0)
  if (rc /= 4) then
    zmq_put_$X = -1
    return
  endif

  rc = f77_zmq_recv(zmq_to_qp_run_socket,msg,len(msg),0)
  if (msg(1:rc) /= 'put_data_reply ok') then
    zmq_put_$X = -1
    return
  endif

end

integer function zmq_get_$X(zmq_to_qp_run_socket, worker_id)
  use f77_zmq
  implicit none
  BEGIN_DOC
! Get $X from the qp_run scheduler
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer, intent(in)            :: worker_id
  integer                        :: rc
  character*(256)                :: msg

  PROVIDE zmq_state
  zmq_get_$X = 0
  if (mpi_master) then
    write(msg,'(A,1X,I8,1X,A200)') 'get_data '//trim(zmq_state), worker_id, '$X'
    rc = f77_zmq_send(zmq_to_qp_run_socket,trim(msg),len(trim(msg)),0)
    if (rc /= len(trim(msg))) then
      zmq_get_$X = -1
      go to 10
    endif

    rc = f77_zmq_recv(zmq_to_qp_run_socket,msg,len(msg),0)
    if (msg(1:14) /= 'get_data_reply') then
      zmq_get_$X = -1
      go to 10
    endif

    rc = f77_zmq_recv(zmq_to_qp_run_socket,$X,4,0)
    if (rc /= 4) then
      zmq_get_$X = -1
      go to 10
    endif

  endif

 10 continue
  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr

    call MPI_BCAST (zmq_get_$X, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to broadcast zmq_get_psi_det'
    endif
    call MPI_BCAST ($X, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to broadcast zmq_get_psi_det'
    endif
  IRP_ENDIF

end

SUBST [ X ]

N_states ;;
N_det ;;
psi_det_size ;;

END_TEMPLATE

integer*8 function zmq_put_psi_det(zmq_to_qp_run_socket,worker_id)
  use f77_zmq
  implicit none
  BEGIN_DOC
! Put psi_det on the qp_run scheduler
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer, intent(in)            :: worker_id
  integer*8                      :: rc8
  character*(256)                :: msg

  integer*8 :: zmq_put_i8matrix
  integer   :: ni, nj

  if (size(psi_det,kind=8) <= 8388608_8) then
    ni = size(psi_det,kind=4)
    nj = 1
  else
    ni = 8388608_8
    nj = int(size(psi_det,kind=8)/8388608_8,4) + 1
  endif
  rc8 = zmq_put_i8matrix(zmq_to_qp_run_socket, 1, 'psi_det', psi_det, ni, nj, size(psi_det,kind=8))
  zmq_put_psi_det = rc8
end

integer*8 function zmq_put_psi_coef(zmq_to_qp_run_socket,worker_id)
  use f77_zmq
  implicit none
  BEGIN_DOC
! Put psi_coef on the qp_run scheduler
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer, intent(in)            :: worker_id
  integer*8                      :: rc8
  character*(256)                :: msg

  zmq_put_psi_coef = 0

  integer*8 :: zmq_put_dmatrix
  integer   :: ni, nj

  if (size(psi_coef,kind=8) <= 8388608_8) then
    ni = size(psi_coef,kind=4)
    nj = 1
  else
    ni = 8388608
    nj = int(size(psi_coef,kind=8)/8388608_8,4) + 1
  endif
  rc8 = zmq_put_dmatrix(zmq_to_qp_run_socket, 1, 'psi_coef', psi_coef, ni, nj, size(psi_coef,kind=8) )
  zmq_put_psi_coef = rc8
end

integer*8 function zmq_get_psi_det(zmq_to_qp_run_socket,worker_id)
  use f77_zmq
  implicit none
  BEGIN_DOC
! Get psi_det on the qp_run scheduler
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer, intent(in)            :: worker_id
  integer*8                      :: rc8
  character*(256)                :: msg

  integer*8 :: zmq_get_i8matrix
  integer   :: ni, nj

  if (size(psi_det,kind=8) <= 8388608_8) then
    ni = size(psi_det,kind=4)
    nj = 1
  else
    ni = 8388608
    nj = int(size(psi_det,kind=8)/8388608_8,4) + 1
  endif
  rc8 = zmq_get_i8matrix(zmq_to_qp_run_socket, 1, 'psi_det', psi_det, ni, nj, size(psi_det,kind=8))
  zmq_get_psi_det = rc8
end

integer*8 function zmq_get_psi_coef(zmq_to_qp_run_socket,worker_id)
  use f77_zmq
  implicit none
  BEGIN_DOC
! get psi_coef on the qp_run scheduler
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer, intent(in)            :: worker_id
  integer*8                      :: rc8
  character*(256)                :: msg

  zmq_get_psi_coef = 0_8

  integer*8 :: zmq_get_dmatrix
  integer   :: ni, nj

  if (size(psi_coef,kind=8) <= 8388608_8) then
    ni = size(psi_coef,kind=4)
    nj = 1
  else
    ni = 8388608
    nj = int(size(psi_coef,kind=8)/8388608_8,4) + 1
  endif
  rc8 = zmq_get_dmatrix(zmq_to_qp_run_socket, 1, 'psi_coef', psi_coef, ni, nj, size(psi_coef,kind=8) )
  zmq_get_psi_coef = rc8
end

!---------------------------------------------------------------------------


integer function zmq_get_psi(zmq_to_qp_run_socket, worker_id)
  use f77_zmq
  implicit none
  BEGIN_DOC
! Get the wave function from the qp_run scheduler
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer, intent(in)            :: worker_id

  integer, external              :: zmq_get_N_states
  integer, external              :: zmq_get_N_det
  integer, external              :: zmq_get_psi_det_size
  integer*8, external            :: zmq_get_psi_det
  integer*8, external            :: zmq_get_psi_coef

  zmq_get_psi = 0

  if (zmq_get_N_states(zmq_to_qp_run_socket, worker_id) == -1) then
    zmq_get_psi = -1
    return
  endif
  if (zmq_get_N_det(zmq_to_qp_run_socket, worker_id) == -1) then
    zmq_get_psi = -1
    return
  endif
  if (zmq_get_psi_det_size(zmq_to_qp_run_socket, worker_id) == -1) then
    zmq_get_psi = -1
    return
  endif
  
  if (size(psi_det,kind=8) /= N_int*2_8*psi_det_size*bit_kind) then
    deallocate(psi_det)
    allocate(psi_det(N_int,2,psi_det_size))
  endif

  if (size(psi_coef,kind=8) /= psi_det_size*N_states) then
    deallocate(psi_coef)
    allocate(psi_coef(psi_det_size,N_states))
  endif

  if (zmq_get_psi_det(zmq_to_qp_run_socket, worker_id) == -1_8) then
    zmq_get_psi = -1
    return
  endif
  if (zmq_get_psi_coef(zmq_to_qp_run_socket, worker_id) == -1_8) then
    zmq_get_psi = -1
    return
  endif
  SOFT_TOUCH psi_det psi_coef psi_det_size N_det N_states

end




