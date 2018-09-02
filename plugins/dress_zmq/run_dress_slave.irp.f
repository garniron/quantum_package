use bitmasks


subroutine run_dress_slave(thread,iproce,energy)
  use f77_zmq
  use omp_lib
  implicit none

  double precision, intent(in)    :: energy(N_states_diag)
  integer,  intent(in)            :: thread, iproce
  integer                        :: rc, i, subset, i_generator

  integer                        :: worker_id, ctask, ltask
  character*(5120)                :: task

  integer(ZMQ_PTR),external      :: new_zmq_to_qp_run_socket
  integer(ZMQ_PTR)               :: zmq_to_qp_run_socket

  integer(ZMQ_PTR), external     :: new_zmq_push_socket
  integer(ZMQ_PTR)               :: zmq_socket_push
  
  double precision,allocatable :: breve_delta_m(:,:,:)
  integer :: i_state,m,l,t,p,sum_f
  !integer, external :: omp_get_thread_num 
  double precision, allocatable :: delta_det(:,:,:,:), cp(:,:,:,:), edI(:)
  double precision, allocatable :: edI_task(:)
  integer, allocatable :: edI_index(:), edI_taskID(:)
  integer :: n_tasks
  
  integer :: iproc
  integer, allocatable :: f(:)
  integer :: cp_sent, cp_done
  integer :: cp_max(Nproc)
  integer :: will_send, task_id, purge_task_id
  integer(kind=OMP_LOCK_KIND) :: lck_det(0:pt2_N_teeth+1)
  integer(kind=OMP_LOCK_KIND) :: lck_sto(0:dress_N_cp+1), sending
  double precision :: fac
  double precision :: ending(1)
  integer, external :: zmq_get_dvector
      
  if(iproce /= 0) stop "RUN DRESS SLAVE is OMP"
  
  allocate(delta_det(N_states, N_det, 0:pt2_N_teeth+1, 2))
  allocate(cp(N_states, N_det, dress_N_cp, 2))
  allocate(edI(N_det_generators), f(N_det_generators))
  allocate(edI_index(N_det_generators), edI_task(N_det_generators))
  edI = 0d0
  f = 0
  delta_det = 0d0

  task(:) = CHAR(0)
  
  call omp_init_lock(sending)
  do i=0,dress_N_cp+1
    call omp_init_lock(lck_sto(i))
  end do
  do i=0,pt2_N_teeth+1
    call omp_init_lock(lck_det(i))
  end do 
  
  cp_done = 0
  cp_sent = 0
  will_send = 0

  double precision :: hij, sij, tmp
  logical :: purge
  purge_task_id = 0
  hij = E0_denominator(1)  !PROVIDE BEFORE OMP PARALLEL
  ending(1) = dble(dress_N_cp+1)
  !$OMP PARALLEL DEFAULT(SHARED) &
  !$OMP PRIVATE(breve_delta_m, task, task_id) &
  !$OMP PRIVATE(tmp,fac,m,l,t,sum_f,n_tasks) &
  !$OMP PRIVATE(i,p,will_send, i_generator, subset, iproc) &
  !$OMP PRIVATE(zmq_to_qp_run_socket, zmq_socket_push, worker_id)
  zmq_to_qp_run_socket = new_zmq_to_qp_run_socket()
  zmq_socket_push      = new_zmq_push_socket(thread)
  call connect_to_taskserver(zmq_to_qp_run_socket,worker_id,thread)
  if(worker_id == -1) then
    call end_zmq_to_qp_run_socket(zmq_to_qp_run_socket)
    call end_zmq_push_socket(zmq_socket_push,thread)
    stop "WORKER -1"
  end if
  
  iproc = omp_get_thread_num()+1
  allocate(breve_delta_m(N_states,N_det,2))
  
  
  do while(cp_done > cp_sent .or. m /= dress_N_cp+1)
    call get_task_from_taskserver(zmq_to_qp_run_socket,worker_id, task_id, task)
    task = task//"  0"
    if(task_id /= 0) then
      read (task,*) subset, i_generator
      m = dress_P(i_generator)
    else
      m = dress_N_cp + 1
      i= zmq_get_dvector(zmq_to_qp_run_socket, worker_id, "ending", ending, 1)
    end if
      
    will_send = 0
    
    !$OMP CRITICAL
    cp_max(iproc) = m
    cp_done = minval(cp_max)-1
    if(cp_done > cp_sent) then
      will_send = cp_sent + 1
      cp_sent = will_send
    end if
    if(purge_task_id == 0) then
      purge_task_id = task_id
      task_id = 0
    end if
    !$OMP END CRITICAL
    
    if(will_send /= 0 .and. will_send <= int(ending(1))) then
      breve_delta_m = 0d0
      
      do l=will_send, 1,-1
        breve_delta_m(:,:,1) += cp(:,:,l,1)
        breve_delta_m(:,:,2) += cp(:,:,l,2)
      end do

      breve_delta_m(:,:,:) = breve_delta_m(:,:,:) / dress_M_m(will_send)
      
      do t=dress_dot_t(will_send)-1,0,-1
        breve_delta_m(:,:,1) = breve_delta_m(:,:,1) + delta_det(:,:,t,1)
        breve_delta_m(:,:,2) = breve_delta_m(:,:,2) + delta_det(:,:,t,2)
      end do
      
      call omp_set_lock(sending)
      n_tasks = 0
      sum_f = 0
      do i=1,N_det_generators
        if(dress_P(i) <= will_send) sum_f = sum_f + f(i)
        if(dress_P(i) == will_send .and. f(i) /= 0) then
          n_tasks += 1
          edI_task(n_tasks) = edI(i)
          edI_index(n_tasks) = i
        end if
      end do
      call push_dress_results(zmq_socket_push, will_send, sum_f, edI_task, edI_index, breve_delta_m, 0, n_tasks)
      call omp_unset_lock(sending)
    end if
    
    if(m /= dress_N_cp+1) then    
      !UPDATE i_generator

      breve_delta_m(:,:,:) = 0d0
      call generator_start(i_generator, iproc)

      call alpha_callback(breve_delta_m, i_generator, subset, pt2_F(i_generator)*0 + 1, iproc)
      
      t = dress_T(i_generator)
    
      call omp_set_lock(lck_det(t))
      delta_det(:,:,t, 1) += breve_delta_m(:,:,1)
      delta_det(:,:,t, 2) += breve_delta_m(:,:,2)
      call omp_unset_lock(lck_det(t))
    
      do p=1,dress_N_cp
        if(dress_e(i_generator, p) /= 0d0) then
          fac = dress_e(i_generator, p)
          call omp_set_lock(lck_sto(p))
          cp(:,:,p,1) += breve_delta_m(:,:,1) * fac
          cp(:,:,p,2) += breve_delta_m(:,:,2) * fac
          call omp_unset_lock(lck_sto(p))
        end if
      end do
      
      tmp = 0d0
      do i=N_det,1,-1
        tmp += psi_coef(i, dress_stoch_istate)*breve_delta_m(dress_stoch_istate, i, 1)
      end do
      !$OMP ATOMIC
      edI(i_generator) += tmp
      !$OMP ATOMIC
      f(i_generator) += 1
      !push bidon
      if(task_id /= 0) then
        call push_dress_results(zmq_socket_push, 0, 0, edI_task, edI_index, breve_delta_m, task_id, 1)
      end if
    end if
  end do
  !$OMP BARRIER
  !$OMP SINGLE
    if(purge_task_id /= 0) then
      do while(int(ending(1)) == dress_N_cp+1)
        call sleep(1)
        i= zmq_get_dvector(zmq_to_qp_run_socket, worker_id, "ending", ending, 1)
      end do

      will_send = int(ending(1))
      breve_delta_m = 0d0
      
      do l=will_send, 1,-1
        breve_delta_m(:,:,1) += cp(:,:,l,1)
        breve_delta_m(:,:,2) += cp(:,:,l,2)
      end do

      breve_delta_m(:,:,:) = breve_delta_m(:,:,:) / dress_M_m(will_send)
      
      do t=dress_dot_t(will_send)-1,0,-1
        breve_delta_m(:,:,1) = breve_delta_m(:,:,1) + delta_det(:,:,t,1)
        breve_delta_m(:,:,2) = breve_delta_m(:,:,2) + delta_det(:,:,t,2)
      end do
      
      sum_f = 0
      do i=1,N_det_generators
        if(dress_P(i) <= will_send) sum_f = sum_f + f(i)
      end do
      call push_dress_results(zmq_socket_push, -will_send, sum_f, edI_task, edI_index, breve_delta_m, purge_task_id, 1)
    end if
   
  !$OMP END SINGLE
  call disconnect_from_taskserver(zmq_to_qp_run_socket,worker_id)
  call end_zmq_to_qp_run_socket(zmq_to_qp_run_socket)
  call end_zmq_push_socket(zmq_socket_push,thread)
  !$OMP END PARALLEL
  do i=0,dress_N_cp+1
    call omp_destroy_lock(lck_sto(i))
  end do
  do i=0,pt2_N_teeth+1
    call omp_destroy_lock(lck_det(i))
  end do
end subroutine



subroutine push_dress_results(zmq_socket_push, m_task, f, edI_task, edI_index, breve_delta_m, task_id, n_tasks)
  use f77_zmq
  implicit none
  integer(ZMQ_PTR), intent(in)   :: zmq_socket_push
  integer, intent(in) :: m_task, f, edI_index(n_tasks)
  double precision, intent(in) :: breve_delta_m(N_states, N_det, 2), edI_task(n_tasks)
  integer, intent(in) :: task_id, n_tasks
  integer :: rc, i, j, k
    rc = f77_zmq_send( zmq_socket_push, n_tasks, 4, ZMQ_SNDMORE)
    if(rc /= 4) stop "push1"
    
    rc = f77_zmq_send( zmq_socket_push, task_id, 4, ZMQ_SNDMORE)
    if(rc /= 4) stop "push2"
    
    rc = f77_zmq_send( zmq_socket_push, m_task, 4, ZMQ_SNDMORE)
    if(rc /= 4) stop "push3"

    rc = f77_zmq_send( zmq_socket_push, f, 4, ZMQ_SNDMORE)
    if(rc /= 4) stop "push4"
    
    rc = f77_zmq_send( zmq_socket_push, edI_task, 8*n_tasks, ZMQ_SNDMORE)
    if(rc /= 8*n_tasks) stop "push5"

    rc = f77_zmq_send( zmq_socket_push, edI_index, 4*n_tasks, ZMQ_SNDMORE)
    if(rc /= 4*n_tasks) stop "push6"
    
    if(m_task < 0) then
      rc = f77_zmq_send( zmq_socket_push, breve_delta_m, 8*N_det*N_states*2, 0)
      if(rc /= 8*N_det*N_states*2) stop "push6"
    else
      rc = f77_zmq_send( zmq_socket_push, breve_delta_m, 8, 0)
      if(rc /= 8) stop "push6"
    end if
! Activate is zmq_socket_pull is a REP
IRP_IF ZMQ_PUSH
IRP_ELSE
  character*(2) :: ok
  rc = f77_zmq_recv( zmq_socket_push, ok, 2, 0)
IRP_ENDIF
end subroutine




subroutine pull_dress_results(zmq_socket_pull, m_task, f, edI_task, edI_index, breve_delta_m, task_id, n_tasks)
  use f77_zmq
  implicit none
  integer(ZMQ_PTR), intent(in)   :: zmq_socket_pull
  integer, intent(out) :: m_task, f, edI_index(N_det_generators)
  double precision, intent(out) :: breve_delta_m(N_states, N_det, 2), edI_task(N_det_generators)
  integer, intent(out) :: task_id, n_tasks
  integer :: rc, i, j, k

    rc = f77_zmq_recv( zmq_socket_pull, n_tasks, 4, 0)
    if(rc /= 4) stop "pullc"
    
    rc = f77_zmq_recv( zmq_socket_pull, task_id, 4, 0)
    if(rc /= 4) stop "pull4"
    
    rc = f77_zmq_recv( zmq_socket_pull, m_task, 4, 0)
    if(rc /= 4) stop "pullc"

    rc = f77_zmq_recv( zmq_socket_pull, f, 4, 0)
    if(rc /= 4) stop "pullc"

    rc = f77_zmq_recv( zmq_socket_pull, edI_task, 8*n_tasks, 0)
    if(rc /= 8*n_tasks) stop "pullc"

    rc = f77_zmq_recv( zmq_socket_pull, edI_index, 4*n_tasks, 0)
    if(rc /= 4*n_tasks) stop "pullc"
    if(m_task < 0) then
      rc = f77_zmq_recv( zmq_socket_pull, breve_delta_m, 8*N_det*N_states*2, 0)
      if(rc /= 8*N_det*N_states*2) stop "pullc"
    else
      rc = f77_zmq_recv( zmq_socket_pull, breve_delta_m, 8, 0)
      if(rc /= 8) stop "pullc"
    end if
! Activate is zmq_socket_pull is a REP
IRP_IF ZMQ_PUSH
IRP_ELSE
  rc = f77_zmq_send( zmq_socket_pull, 'ok', 2, 0)
IRP_ENDIF

end subroutine
 
 
            
