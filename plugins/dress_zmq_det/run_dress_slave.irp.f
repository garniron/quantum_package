use bitmasks

BEGIN_PROVIDER [ integer, fragment_count ]
  implicit none
  BEGIN_DOC
    ! Number of fragments for the deterministic part
  END_DOC
  fragment_count = 1
END_PROVIDER


subroutine run_dress_slave(thread,iproce,energy)
  use f77_zmq
  use omp_lib
  implicit none

  double precision, intent(in)    :: energy(N_states_diag)
  integer,  intent(in)            :: thread, iproce
  integer                        :: rc, i, subset, i_generator

  integer                        :: worker_id, task_id, ctask, ltask
  character*(5120)                :: task

  integer(ZMQ_PTR),external      :: new_zmq_to_qp_run_socket
  integer(ZMQ_PTR)               :: zmq_to_qp_run_socket

  integer(ZMQ_PTR), external     :: new_zmq_push_socket
  integer(ZMQ_PTR)               :: zmq_socket_push

  logical :: done

  double precision,allocatable :: dress_detail(:)
  integer :: ind
  
  double precision,allocatable :: delta_ij_loc(:,:,:)
  integer :: h,p,n,i_state
  logical :: ok

  integer, allocatable :: int_buf(:)
  double precision, allocatable :: double_buf(:)
  integer(bit_kind), allocatable :: det_buf(:,:,:)
  integer :: N_buf(3)
  logical :: last
  !integer, external :: omp_get_thread_num 
  double precision, allocatable :: delta_det(:,:,:,:), cp(:,:,:,:)
  integer :: toothMwen
  logical :: fracted
  double precision :: fac
   

  if(iproce /= 0) stop "RUN DRESS SLAVE is OMP"
  
  allocate(delta_det(N_states, N_det, 0:comb_teeth+1, 2))
  allocate(cp(N_states, N_det, N_cp, 2))
  delta_det = 0d9
  cp = 0d0

  
  task(:) = CHAR(0)

  
   
  integer :: anproc, iproc, cur_cp, done_for(0:N_cp)
  integer, allocatable :: tasks(:)
  integer :: lastCp(Nproc)
  integer :: lastSent, lastSendable
  logical :: send
  integer(kind=OMP_LOCK_KIND) :: lck_det(0:comb_teeth+1)
  integer(kind=OMP_LOCK_KIND) :: lck_sto(0:N_cp+1)
  
  do i=0,N_cp+1
    call omp_init_lock(lck_sto(i))
  end do
  do i=0,comb_teeth+1
    call omp_init_lock(lck_det(i))
  end do 
  
  lastCp = 0
  lastSent = 0
  send = .false.
  done_for = 0
  anproc = min(nproc, N_det_generators)
  double precision :: hij, sij
  !call i_h_j_s2(psi_det(1,1,1),psi_det(1,1,2),N_int,hij, sij)
  
  hij = E0_denominator(1)  !PROVIDE BEFORE OMP PARALLEL
  
  !$OMP PARALLEL NUM_THREADS(nproc) DEFAULT(SHARED) &
  !$OMP PRIVATE(int_buf, double_buf, det_buf, delta_ij_loc, task, task_id) &
  !$OMP PRIVATE(lastSendable, toothMwen, fracted, fac) &
  !$OMP PRIVATE(i, cur_cp, send, i_generator, subset, iproc, N_buf) &
  !$OMP PRIVATE(zmq_to_qp_run_socket, zmq_socket_push, worker_id)
  task(:) = " "
  zmq_to_qp_run_socket = new_zmq_to_qp_run_socket()
  zmq_socket_push      = new_zmq_push_socket(thread)
  call connect_to_taskserver(zmq_to_qp_run_socket,worker_id,thread)
  if(worker_id == -1) then
    call end_zmq_to_qp_run_socket(zmq_to_qp_run_socket)
    call end_zmq_push_socket(zmq_socket_push,thread)
    stop "WORKER -1"
  end if
  
  
  iproc = omp_get_thread_num()+1
  allocate(int_buf(N_dress_int_buffer)) 
  allocate(double_buf(N_dress_double_buffer)) 
  allocate(det_buf(N_int, 2, N_dress_det_buffer)) 
  allocate(delta_ij_loc(N_states,N_det,2)) 
  do
    call get_task_from_taskserver(zmq_to_qp_run_socket,worker_id, task_id, task)
    task = task//"  0"
    if(task_id == 0) exit
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(task_id /= 0) then
      !print *, "TASKID", task_id
      !print *, "<=",task,"=>"
      read (task,*) subset, i_generator

      delta_ij_loc(:,:,:) = 0d0
      call generator_start(i_generator, iproc)
      call alpha_callback(delta_ij_loc, i_generator, subset, iproc)
      call generator_done(i_generator, int_buf, double_buf, det_buf, N_buf, iproc)

      N_buf = (/0,1,0/)

      call push_dress_results(zmq_socket_push, i_generator, delta_ij_loc, int_buf, double_buf, det_buf, N_buf, task_id)
      call task_done_to_taskserver(zmq_to_qp_run_socket,worker_id,task_id)
    end if

    if(task_id == 0) exit
  end do
  
  call disconnect_from_taskserver(zmq_to_qp_run_socket,worker_id)
  call end_zmq_to_qp_run_socket(zmq_to_qp_run_socket)
  call end_zmq_push_socket(zmq_socket_push,thread)
  !$OMP END PARALLEL

  do i=0,N_cp+1
    call omp_destroy_lock(lck_sto(i))
  end do
  do i=0,comb_teeth+1
    call omp_destroy_lock(lck_det(i))
  end do
end subroutine




subroutine push_dress_results(zmq_socket_push, ind, delta_loc, int_buf, double_buf, det_buf, N_bufi, task_id)
  use f77_zmq
  implicit none
  
  integer, parameter :: sendt = 4
  integer(ZMQ_PTR), intent(in)   :: zmq_socket_push
  double precision, intent(inout)   :: delta_loc(N_states, N_det, 2)
  real(kind=4), allocatable   :: delta_loc4(:,:,:)
  double precision, intent(in) :: double_buf(*)
  integer, intent(in) :: int_buf(*)
  integer(bit_kind), intent(in) :: det_buf(N_int, 2, *)
  integer, intent(in) :: N_bufi(3)
  integer :: N_buf(3)
  integer, intent(in) :: ind, task_id
  integer :: rc, i, j, k, l
  double precision :: contrib(N_states)
  real(sendt), allocatable :: r4buf(:,:,:)
  
  !if(ind < 10) then
  !  print *, "DRESS", ind, delta_loc(1,:5,1)
  !end if

    rc = f77_zmq_send( zmq_socket_push, ind, 4, ZMQ_SNDMORE)
  if(rc /= 4) stop "push"
  

    allocate(r4buf(N_states, N_det, 2))
    do i=1,2
    do j=1,N_det
    do k=1,N_states
       r4buf(k,j,i) = real(delta_loc(k,j,i), sendt)
    end do
    end do
    end do

    rc = f77_zmq_send( zmq_socket_push, r4buf(1,1,1), sendt*N_states*N_det, ZMQ_SNDMORE)
    if(rc /= sendt*N_states*N_det) stop "push"

    rc = f77_zmq_send( zmq_socket_push, r4buf(1,1,2), sendt*N_states*N_det, ZMQ_SNDMORE)
    if(rc /= sendt*N_states*N_det) stop "push"
   


    contrib = 0d0
    do i=1,N_det
      contrib(:) += delta_loc(:,i, 1) * psi_coef(i, :)
    end do

    rc = f77_zmq_send( zmq_socket_push, contrib, 8*N_states, ZMQ_SNDMORE)
    if(rc /= 8*N_states) stop "push"
    
    N_buf = N_bufi
    !N_buf = (/0,1,0/)

    rc = f77_zmq_send( zmq_socket_push, N_buf, 4*3, ZMQ_SNDMORE)
    if(rc /= 4*3) stop "push5" 
    
    if(N_buf(1) > N_dress_int_buffer) stop "run_dress_slave N_buf bad size?"
    if(N_buf(2) > N_dress_double_buffer) stop "run_dress_slave N_buf bad size?"
    if(N_buf(3) > N_dress_det_buffer) stop "run_dress_slave N_buf bad size?"

    
    if(N_buf(1) > 0) then
      rc = f77_zmq_send( zmq_socket_push, int_buf, 4*N_buf(1), ZMQ_SNDMORE)
      if(rc /= 4*N_buf(1)) stop "push6"
    end if
    
    if(N_buf(2) > 0) then
      rc = f77_zmq_send( zmq_socket_push, double_buf, 8*N_buf(2), ZMQ_SNDMORE)
      if(rc /= 8*N_buf(2)) stop "push8"
    end if

    if(N_buf(3) > 0) then
      rc = f77_zmq_send( zmq_socket_push, det_buf, 2*N_int*bit_kind*N_buf(3), ZMQ_SNDMORE)
      if(rc /= 2*N_int*bit_kind*N_buf(3)) stop "push10"
    end if

    rc = f77_zmq_send( zmq_socket_push, task_id, 4, 0)
    if(rc /= 4) stop "push11"

! Activate is zmq_socket_push is a REQ
IRP_IF ZMQ_PUSH
IRP_ELSE
  character*(2) :: ok
  rc = f77_zmq_recv( zmq_socket_push, ok, 2, 0)
IRP_ENDIF

end subroutine


BEGIN_PROVIDER [ real(4), real4buf, (N_states, N_det, 2) ]
  
END_PROVIDER


subroutine pull_dress_results(zmq_socket_pull, ind, delta_loc, int_buf, double_buf, det_buf, N_buf, task_id, contrib)
  use f77_zmq
  implicit none
  integer, parameter :: sendt = 4
  integer(ZMQ_PTR), intent(in)   :: zmq_socket_pull
  double precision, intent(inout) :: delta_loc(N_states, N_det, 2)
  double precision, intent(out) :: double_buf(*), contrib(N_states)
  integer, intent(out) :: int_buf(*)
  integer(bit_kind), intent(out) :: det_buf(N_int, 2, *)
  integer, intent(out) :: ind
  integer, intent(out) :: task_id
  integer :: rc, i, j, k
  integer, intent(out) :: N_buf(3)

  rc = f77_zmq_recv( zmq_socket_pull, ind, 4, 0)
  if(rc /= 4) stop "pulla"
 
    
    rc = f77_zmq_recv( zmq_socket_pull, real4buf(1,1,1), N_states*sendt*N_det, 0)
    if(rc /= sendt*N_states*N_det) stop "pullc"
  
    rc = f77_zmq_recv( zmq_socket_pull, real4buf(1,1,2), N_states*sendt*N_det, 0)
    if(rc /= sendt*N_states*N_det) stop "pulld"
    
    do i=1,2
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j,k)
    do j=1,N_det
    do k=1,N_states
       delta_loc(k,j,i) = real(real4buf(k,j,i), 8)
    end do
    end do
    end do
    
    
    
    rc = f77_zmq_recv( zmq_socket_pull, contrib, 8*N_states, 0)
    if(rc /= 8*N_states) stop "pullc"
  
    rc = f77_zmq_recv( zmq_socket_pull, N_buf, 4*3, 0)
    if(rc /= 4*3) stop "pull" 
    if(N_buf(1) > N_dress_int_buffer) stop "run_dress_slave N_buf bad size?"
    if(N_buf(2) > N_dress_double_buffer) stop "run_dress_slave N_buf bad size?"
    if(N_buf(3) > N_dress_det_buffer) stop "run_dress_slave N_buf bad size?"

    
    if(N_buf(1) > 0) then
      rc = f77_zmq_recv( zmq_socket_pull, int_buf, 4*N_buf(1), 0)
      if(rc /= 4*N_buf(1)) stop "pull1"
    end if
    
    if(N_buf(2) > 0) then
      rc = f77_zmq_recv( zmq_socket_pull, double_buf, 8*N_buf(2), 0)
      if(rc /= 8*N_buf(2)) stop "pull2"
    end if
    
    if(N_buf(3) > 0) then
      rc = f77_zmq_recv( zmq_socket_pull, det_buf, 2*N_int*bit_kind*N_buf(3), 0)
      if(rc /= 2*N_int*bit_kind*N_buf(3)) stop "pull3"
    end if

    rc = f77_zmq_recv( zmq_socket_pull, task_id, 4, 0)
    if(rc /= 4) stop "pull4"
! Activate is zmq_socket_pull is a REP
IRP_IF ZMQ_PUSH
IRP_ELSE
  rc = f77_zmq_send( zmq_socket_pull, 'ok', 2, 0)
IRP_ENDIF

end subroutine
 
 
            
