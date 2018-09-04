BEGIN_PROVIDER [ integer, pt2_stoch_istate ]
 implicit none
 BEGIN_DOC
 ! State for stochatsic PT2 
 END_DOC
 pt2_stoch_istate = 1
END_PROVIDER


 BEGIN_PROVIDER [ integer, pt2_N_teeth ]
&BEGIN_PROVIDER [ integer, pt2_minDetInFirstTeeth ]
&BEGIN_PROVIDER [ integer, pt2_n_tasks_max ]
&BEGIN_PROVIDER [ integer, pt2_F, (N_det_generators) ]
  implicit none
  logical, external :: testTeethBuilding
  pt2_F(:) = 1
  !pt2_F(:N_det_generators/1000*0+50) = 1
  pt2_n_tasks_max = 16 ! N_det_generators/100 + 1
  
  if(N_det_generators < 1024) then
    pt2_minDetInFirstTeeth = 1
    pt2_N_teeth = 1
  else
    do pt2_N_teeth=32,1,-1
      pt2_minDetInFirstTeeth = min(5, N_det_generators)
      if(testTeethBuilding(pt2_minDetInFirstTeeth, pt2_N_teeth)) exit
    end do
  end if
END_PROVIDER


logical function testTeethBuilding(minF, N)
  implicit none
  integer, intent(in) :: minF, N
  integer :: n0, i
  double precision :: u0, Wt, r
  
  double precision, allocatable :: tilde_w(:), tilde_cW(:)
  integer, external :: dress_find_sample

  allocate(tilde_w(N_det_generators), tilde_cW(0:N_det_generators))
  
  tilde_cW(0) = 0d0
  do i=1,N_det_generators
    tilde_w(i)  = psi_coef_generators(i,pt2_stoch_istate)**2
    tilde_cW(i) = tilde_cW(i-1) + tilde_w(i)
  enddo
  tilde_cW(N_det_generators) = 1d0
 
  n0 = 0
  do
    u0 = tilde_cW(n0)
    r = tilde_cW(n0 + minF)
    Wt = (1d0 - u0) / dble(N)
    if(Wt >= r - u0) then
       testTeethBuilding = .true.
       return
    end if
    n0 += 1
    if(N_det_generators - n0 < minF * N) then
      testTeethBuilding = .false.
      return
    end if
  end do
  stop "exited testTeethBuilding"
end function



subroutine ZMQ_pt2(E, pt2,relative_error, absolute_error, error)
  use f77_zmq
  use selection_types
  
  implicit none
  
  character(len=64000)           :: task
  integer(ZMQ_PTR)               :: zmq_to_qp_run_socket, zmq_socket_pull
  integer, external              :: omp_get_thread_num
  double precision, intent(in)   :: relative_error, absolute_error, E(N_states)
  double precision, intent(out)  :: pt2(N_states),error(N_states)
  
  
  integer                        :: i, j, k
  
  double precision, external     :: omp_get_wtime
  double precision               :: state_average_weight_save(N_states), w(N_states)
  integer(ZMQ_PTR), external     :: new_zmq_to_qp_run_socket
  
  if (N_det < max(10,N_states)) then
    pt2=0.d0
    call ZMQ_selection(0, pt2)
    error(:) = 0.d0
  else
    
    state_average_weight_save(:) = state_average_weight(:)
    do pt2_stoch_istate=1,N_states
      SOFT_TOUCH pt2_stoch_istate
      state_average_weight(:) = 0.d0
      state_average_weight(pt2_stoch_istate) = 1.d0
      TOUCH state_average_weight 

      provide nproc pt2_F mo_bielec_integrals_in_map mo_mono_elec_integral pt2_w psi_selectors 
      
      print *, '========== ================= ================= ================='
      print *, ' Samples        Energy         Stat. Error         Seconds      '
      print *, '========== ================= ================= ================='
      
      call new_parallel_job(zmq_to_qp_run_socket, zmq_socket_pull, 'pt2')

      integer, external              :: zmq_put_psi
      integer, external              :: zmq_put_N_det_generators
      integer, external              :: zmq_put_N_det_selectors
      integer, external              :: zmq_put_dvector
      if (zmq_put_psi(zmq_to_qp_run_socket,1) == -1) then
        stop 'Unable to put psi on ZMQ server'
      endif
      if (zmq_put_N_det_generators(zmq_to_qp_run_socket, 1) == -1) then
        stop 'Unable to put N_det_generators on ZMQ server'
      endif
      if (zmq_put_N_det_selectors(zmq_to_qp_run_socket, 1) == -1) then
        stop 'Unable to put N_det_selectors on ZMQ server'
      endif
      if (zmq_put_dvector(zmq_to_qp_run_socket,1,'energy',pt2_e0_denominator,size(pt2_e0_denominator)) == -1) then
        stop 'Unable to put energy on ZMQ server'
      endif
      
      integer, external :: add_task_to_taskserver
      
      
      do i=1,N_det_generators
        do j=1,pt2_F(i) !!!!!!!!!!!! 
          write(task(1:20),'(I9,1X,I9''|'')') j, pt2_J(i)
          if (add_task_to_taskserver(zmq_to_qp_run_socket,trim(task(1:20))) == -1) then
            stop 'Unable to add task to task server'
          endif
        end do
      end do
      
      integer, external :: zmq_set_running
      if (zmq_set_running(zmq_to_qp_run_socket) == -1) then
        print *,  irp_here, ': Failed in zmq_set_running'
      endif

      
      !$OMP PARALLEL DEFAULT(shared) NUM_THREADS(nproc+1)            &
          !$OMP  PRIVATE(i)
      i = omp_get_thread_num()
      if (i==0) then
        call pt2_collector(zmq_socket_pull, E(pt2_stoch_istate),relative_error, absolute_error, w, error)
        pt2(pt2_stoch_istate) = w(pt2_stoch_istate)
      else
        call pt2_slave_inproc(i)
      endif
      !$OMP END PARALLEL
      call end_parallel_job(zmq_to_qp_run_socket, zmq_socket_pull, 'pt2')
      
      print *, '========== ================= ================= ================='
      
    enddo
    FREE pt2_stoch_istate 
    state_average_weight(:) = state_average_weight_save(:)
    TOUCH state_average_weight
  endif
  do k=N_det+1,N_states
    pt2(k) = 0.d0
  enddo

end subroutine


subroutine pt2_slave_inproc(i)
  implicit none
  integer, intent(in)            :: i

  call run_pt2_slave(1,i,pt2_e0_denominator)
end


subroutine pt2_collector(zmq_socket_pull, E, relative_error, absolute_error, pt2, error)
  use f77_zmq
  use selection_types
  use bitmasks
  implicit none

  
  integer(ZMQ_PTR), intent(in)   :: zmq_socket_pull
  double precision, intent(in) :: relative_error, absolute_error, E
  double precision, intent(out)  :: pt2(N_states), error(N_states)


  double precision, allocatable      :: eI(:,:), eI_task(:,:), S(:), S2(:)
  integer(ZMQ_PTR),external      :: new_zmq_to_qp_run_socket
  integer(ZMQ_PTR)               :: zmq_to_qp_run_socket
  integer, external :: zmq_delete_tasks
  integer, external :: pt2_find_sample

  integer :: more, n, i, p, c, t, n_tasks, U
  integer, allocatable :: task_id(:)
  integer, allocatable :: index(:)
  
  double precision, external :: omp_get_wtime
  double precision :: v, x, avg, eqt, E0
  double precision :: time, time0
  
  integer, allocatable :: f(:)
  logical, allocatable :: d(:) 

  zmq_to_qp_run_socket = new_zmq_to_qp_run_socket()
  allocate(task_id(pt2_n_tasks_max), index(pt2_n_tasks_max), f(N_det_generators))
  allocate(d(N_det_generators+1))
  allocate(eI(N_states, N_det_generators), eI_task(N_states, pt2_n_tasks_max))
  allocate(S(pt2_N_teeth+1), S2(pt2_N_teeth+1))
   
  S(:) = 0d0
  S2(:) = 0d0
  n = 1
  t = 0
  U = 0
  eI(:,:) = 0d0
  f(:) = pt2_F(:)
  d(:) = .false.
  n_tasks = 0
  E0 = E
  more = 1
  time0 = omp_get_wtime()

  do while (n <= N_det_generators)
    if(f(pt2_J(n)) == 0) then
      d(pt2_J(n)) = .true.
      do while(d(U+1))
        U += 1
      end do
      do while(t <= pt2_N_teeth)
        if(U >= pt2_n_0(t+1)) then
          t=t+1
          E0 = E
          do i=pt2_n_0(t),1,-1
            E0 += eI(pt2_stoch_istate, i)
          end do
        else
          exit
        end if
      end do

      c = pt2_R(n)
      if(c /= 0) then
        x = 0d0
        do p=pt2_N_teeth, 1, -1
          v = pt2_u_0 + pt2_W_T * (pt2_u(c) + dble(p-1))
          i = pt2_find_sample(v, pt2_cW)
          x += eI(pt2_stoch_istate, i) * pt2_W_T / pt2_w(i)
          S(p) += x
          S2(p) += x**2
        end do
        avg = S(t) / dble(c)
        eqt = (S2(t) / c) - (S(t)/c)**2
        eqt = sqrt(eqt / dble(c-1))
        pt2(pt2_stoch_istate) = E0-E+avg
        error(pt2_stoch_istate) = eqt
        time = omp_get_wtime()
        print '(G10.3, 2X, F16.10, 2X, G16.3, 2X, F16.4, A20)', c, avg+E0, eqt, time-time0, ''
      end if
      n += 1
    else if(more == 0) then
      exit
    else
      call pull_pt2_results(zmq_socket_pull, index, eI_task, task_id, n_tasks)
      if (zmq_delete_tasks(zmq_to_qp_run_socket,zmq_socket_pull,task_id,n_tasks,more) == -1) then
          stop 'Unable to delete tasks'
      endif
      do i=1,n_tasks
        eI(:, index(i)) += eI_task(:, i)
        f(index(i)) -= 1
      end do
    end if
  end do
  call end_zmq_to_qp_run_socket(zmq_to_qp_run_socket)
end subroutine


integer function pt2_find_sample(v, w)
  implicit none
  double precision, intent(in) :: v, w(0:N_det_generators)
  integer :: i,l,r

  l = 0
  r = N_det_generators

  do while(r-l > 1)
    i = (r+l) / 2
    if(w(i) < v) then
      l = i
    else
      r = i
    end if
  end do

  pt2_find_sample = r
end function


 BEGIN_PROVIDER[ integer, pt2_J, (N_det_generators)]
&BEGIN_PROVIDER[ double precision, pt2_u, (N_det_generators)]
&BEGIN_PROVIDER[ integer, pt2_R, (N_det_generators)]
  implicit none
  integer :: N_c, N_j, U, t, i
  double precision :: v
  logical, allocatable :: d(:)
  integer, external :: pt2_find_sample
  
  allocate(d(N_det_generators))
  
  pt2_R(:) = 0
  N_c = 0
  N_j = pt2_n_0(1)
  d(:) = .false.

  do i=1,N_j
      d(i) = .true.
      pt2_J(i) = i
  end do
  call random_seed(put=(/3211,64,6566,321,65,321,654,65,321,6321,654,65,321,621,654,65,321,65,654,65,321,65/)) 
  call RANDOM_NUMBER(pt2_u)
  call RANDOM_NUMBER(pt2_u)
  

  
  U = 0
  
  do while(N_j < N_det_generators)
    !ADD_COMB
    N_c += 1
    do t=0, pt2_N_teeth-1
      v = pt2_u_0 + pt2_W_T * (dble(t) + pt2_u(N_c))
      i = pt2_find_sample(v, pt2_cW)
      if(.not. d(i)) then
        N_j += 1
        pt2_J(N_j) = i
        d(i) = .true.
      end if
    end do
    
    pt2_R(N_j) = N_c
    
    !FILL_TOOTH
    do while(U < N_det_generators)
      U += 1
      if(.not. d(U)) then
        N_j += 1
        pt2_J(N_j) = U
        d(U) = .true.
        exit;
      end if
    end do
  enddo
  if(N_det_generators > 1) then
    pt2_R(N_det_generators-1) = 0
    pt2_R(N_det_generators) = N_c
  end if
END_PROVIDER


 BEGIN_PROVIDER [ double precision,     pt2_w, (N_det_generators) ] 
&BEGIN_PROVIDER [ double precision,     pt2_cW, (0:N_det_generators) ] 
&BEGIN_PROVIDER [ double precision,     pt2_W_T ]
&BEGIN_PROVIDER [ double precision,     pt2_u_0 ]
&BEGIN_PROVIDER [ integer,              pt2_n_0, (pt2_N_teeth+1) ]
  implicit none
  integer :: i, t
  double precision, allocatable :: tilde_w(:), tilde_cW(:)
  double precision :: r, tooth_width
  integer, external :: pt2_find_sample

  allocate(tilde_w(N_det_generators), tilde_cW(0:N_det_generators))
  
  tilde_cW(0) = 0d0
  
  do i=1,N_det_generators
    tilde_w(i)  = psi_coef_generators(i,pt2_stoch_istate)**2
    tilde_cW(i) = tilde_cW(i-1) + tilde_w(i)
  enddo
  
  pt2_n_0(1) = 0
  do
    pt2_u_0 = tilde_cW(pt2_n_0(1))
    r = tilde_cW(pt2_n_0(1) + pt2_minDetInFirstTeeth)
    pt2_W_T = (1d0 - pt2_u_0) / dble(pt2_N_teeth)
    if(pt2_W_T >= r - pt2_u_0) then
      exit
    end if
    pt2_n_0(1) += 1
    if(N_det_generators - pt2_n_0(1) < pt2_minDetInFirstTeeth * pt2_N_teeth) then
      stop "teeth building failed"
    end if
  end do
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
  do t=2, pt2_N_teeth
    r = pt2_u_0 + pt2_W_T * dble(t-1)
    pt2_n_0(t) = pt2_find_sample(r, tilde_cW)
  end do
  pt2_n_0(pt2_N_teeth+1) = N_det_generators
    
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  pt2_w(:pt2_n_0(1)) = tilde_w(:pt2_n_0(1))
  do t=1, pt2_N_teeth
    tooth_width = tilde_cW(pt2_n_0(t+1)) - tilde_cW(pt2_n_0(t))
    do i=pt2_n_0(t)+1, pt2_n_0(t+1)
      pt2_w(i) = tilde_w(i) * pt2_W_T / tooth_width
    end do
  end do
  
  pt2_cW(0) = 0d0
  do i=1,N_det_generators
    pt2_cW(i) = pt2_cW(i-1) + pt2_w(i)      
  end do
  pt2_n_0(pt2_N_teeth+1) = N_det_generators
END_PROVIDER





