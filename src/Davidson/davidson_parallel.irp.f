
!brought to you by garniroy inc.

use bitmasks
use f77_zmq

subroutine davidson_process(blockb, blocke, N, idx, vt, st)
  
  implicit none
  

  integer             , intent(in)        :: blockb, blocke
  integer             , intent(inout)     :: N
  integer             , intent(inout)     :: idx(dav_size)
  double precision    , intent(inout)     :: vt(N_states_diag, dav_size)
  double precision    , intent(inout)     :: st(N_states_diag, dav_size)
  
  integer :: i, j, sh, sh2, exa, ext, org_i, org_j, istate, ni, endi
  integer(bit_kind) :: sorted_i(N_int)
  double precision :: s2, hij
  logical, allocatable :: wrotten(:)
  
  allocate(wrotten(dav_size))
  wrotten = .false.

  do sh = blockb, blocke
  do sh2=1,sh
    exa = 0
    do ni=1,N_int
      exa = exa + popcnt(xor(version_(ni,sh,1), version_(ni,sh2,1)))
    end do
    if(exa > 2) then
      cycle
    end if
    
    do i=shortcut_(sh,1),shortcut_(sh+1,1)-1
      org_i = sort_idx_(i,1)
      if(sh==sh2) then
        endi = i-1
      else
        endi = shortcut_(sh2+1,1)-1
      end if
      do ni=1,N_int
        sorted_i(ni) = sorted_(ni,i,1)
      enddo
      
      do j=shortcut_(sh2,1),endi
        org_j = sort_idx_(j,1)
        ext = exa
        do ni=1,N_int
          ext = ext + popcnt(xor(sorted_i(ni), sorted_(ni,j,1)))
        end do
        if(ext <= 4) then
          call i_h_j (dav_det(1,1,org_j),dav_det(1,1,org_i),n_int,hij)
          call get_s2(dav_det(1,1,org_j),dav_det(1,1,org_i),n_int,s2) 
          if(.not. wrotten(org_i)) then
            wrotten(org_i) = .true.
            vt (:,org_i) = 0d0
            st (:,org_i) = 0d0
          end if
          if(.not. wrotten(org_j)) then
            wrotten(org_j) = .true.
            vt (:,org_j) = 0d0
            st (:,org_j) = 0d0
          end if
          do istate=1,N_states_diag
            vt (istate,org_i) += hij*dav_ut(istate,org_j)
            st (istate,org_i) += s2*dav_ut(istate,org_j)
            vt (istate,org_j) += hij*dav_ut(istate,org_i)
            st (istate,org_j) += s2*dav_ut(istate,org_i)
          enddo
        endif
      enddo
    enddo
  enddo
  enddo

  N = 0
  do i=1, dav_size
    if(wrotten(i)) then
      N = N+1
      do istate=1,N_states_diag
        vt (istate,N) = vt (istate,i)
        st (istate,N) = st (istate,i)
        idx(N) = i
      enddo
    end if
  end do
end subroutine




subroutine davidson_collect(blockb, blocke, N, idx, vt, st , v0t, s0t)
  implicit none


  integer             , intent(in)        :: blockb, blocke
  integer             , intent(in)        :: N
  integer             , intent(in)        :: idx(N)
  double precision    , intent(in)        :: vt(N_states_diag, N)
  double precision    , intent(in)        :: st(N_states_diag, N)
  double precision    , intent(inout)     :: v0t(N_states_diag,dav_size)
  double precision    , intent(inout)     :: s0t(N_states_diag,dav_size)

  integer :: i, j, k

  !DIR$ IVDEP
  do i=1,N
    k = idx(i)
    !DIR$ IVDEP
    do j=1,N_states_diag
      v0t(j,k) = v0t(j,k) + vt(j,i)
      s0t(j,k) = s0t(j,k) + st(j,i)
    enddo
  end do
end subroutine


subroutine davidson_init(zmq_to_qp_run_socket)
  use f77_zmq
  implicit none
  
  integer(ZMQ_PTR), intent(out) :: zmq_to_qp_run_socket
  
  touch dav_size dav_det dav_ut
  call new_parallel_job(zmq_to_qp_run_socket,"davidson")
end subroutine



subroutine davidson_add_task(zmq_to_qp_run_socket, blockb, blocke)
  use f77_zmq
  implicit none
  
  integer(ZMQ_PTR)    ,intent(in)    :: zmq_to_qp_run_socket
  integer             ,intent(in)    :: blockb, blocke
  character*(512)                    :: task 
  
  
  write(task,*) blockb, blocke
  call add_task_to_taskserver(zmq_to_qp_run_socket, task)
end subroutine



subroutine davidson_slave_inproc(i)
  implicit none
  integer, intent(in)            :: i

  call davidson_run_slave(1,i)
end


subroutine davidson_slave_tcp(i)
  implicit none
  integer, intent(in)            :: i
  
  call davidson_run_slave(0,i)
end



subroutine davidson_run_slave(thread,iproc)
  use f77_zmq
  implicit none

  integer,  intent(in)           :: thread, iproc

  integer                        :: worker_id, task_id, blockb, blocke
  character*(512)                :: task

  integer(ZMQ_PTR),external      :: new_zmq_to_qp_run_socket
  integer(ZMQ_PTR)               :: zmq_to_qp_run_socket

  integer(ZMQ_PTR), external     :: new_zmq_push_socket
  integer(ZMQ_PTR)               :: zmq_socket_push

  

  zmq_to_qp_run_socket = new_zmq_to_qp_run_socket()
  zmq_socket_push      = new_zmq_push_socket(thread)
  call connect_to_taskserver(zmq_to_qp_run_socket,worker_id,thread)
  if(worker_id == -1) then
    print *, "WORKER -1"
    call end_zmq_to_qp_run_socket(zmq_to_qp_run_socket)
    call end_zmq_push_socket(zmq_socket_push,thread)
    return
  end if
  
  call davidson_slave_work(zmq_to_qp_run_socket, zmq_socket_push, worker_id)
  call disconnect_from_taskserver(zmq_to_qp_run_socket,zmq_socket_push,worker_id)
  call end_zmq_to_qp_run_socket(zmq_to_qp_run_socket)
  call end_zmq_push_socket(zmq_socket_push,thread)
end subroutine



subroutine davidson_slave_work(zmq_to_qp_run_socket, zmq_socket_push, worker_id)
  use f77_zmq
  implicit none
  
  integer(ZMQ_PTR),intent(in)   :: zmq_to_qp_run_socket
  integer(ZMQ_PTR),intent(in)   :: zmq_socket_push
  integer,intent(in)            :: worker_id
  integer :: task_id
  character*(512) :: task
  

  integer                  :: blockb, blocke
  integer                  :: N
  integer             ,  allocatable      :: idx(:)
  double precision    ,  allocatable      :: vt(:,:)
  double precision    ,  allocatable      :: st(:,:)
  
  
  allocate(idx(dav_size)) 
  allocate(vt(N_states_diag, dav_size)) 
  allocate(st(N_states_diag, dav_size)) 
  
  
  do
    call get_task_from_taskserver(zmq_to_qp_run_socket,worker_id, task_id, task)
    if(task_id == 0) exit
    read (task,*) blockb, blocke
    
    call davidson_process(blockb, blocke, N, idx, vt, st)
    
    call task_done_to_taskserver(zmq_to_qp_run_socket,worker_id,task_id)
    call davidson_push_results(zmq_socket_push, blockb, blocke, N, idx, vt, st, task_id)
  end do
end subroutine



subroutine davidson_push_results(zmq_socket_push, blockb, blocke, N, idx, vt, st, task_id)
  use f77_zmq
  implicit none
  
  integer(ZMQ_PTR)    ,intent(in)    :: zmq_socket_push
  integer             ,intent(in)    :: task_id

  integer             ,intent(in)    :: blockb, blocke
  integer             ,intent(in)    :: N
  integer             ,intent(in)    :: idx(N)
  double precision    ,intent(in)    :: vt(N_states_diag, N)
  double precision    ,intent(in)    :: st(N_states_diag, N)
  integer                            :: rc

  rc = f77_zmq_send( zmq_socket_push, blockb, 4, ZMQ_SNDMORE)
  if(rc /= 4) stop "davidson_push_results failed to push blockb"

  rc = f77_zmq_send( zmq_socket_push, blocke, 4, ZMQ_SNDMORE)
  if(rc /= 4) stop "davidson_push_results failed to push blocke"
  
  rc = f77_zmq_send( zmq_socket_push, N, 4, ZMQ_SNDMORE)
  if(rc /= 4) stop "davidson_push_results failed to push N"

  rc = f77_zmq_send( zmq_socket_push, idx, 4*N, ZMQ_SNDMORE)
  if(rc /= 4*N) stop "davidson_push_results failed to push idx"

  rc = f77_zmq_send( zmq_socket_push, vt, 8*N_states_diag* N, ZMQ_SNDMORE)
  if(rc /= 8*N_states_diag* N) stop "davidson_push_results failed to push vt"

  rc = f77_zmq_send( zmq_socket_push, st, 8*N_states_diag* N, ZMQ_SNDMORE)
  if(rc /= 8*N_states_diag* N) stop "davidson_push_results failed to push st"

  rc = f77_zmq_send( zmq_socket_push, task_id, 4, 0)
  if(rc /= 4) stop "davidson_push_results failed to push task_id"
end subroutine



subroutine davidson_pull_results(zmq_socket_pull, blockb, blocke, N, idx, vt, st, task_id)
  use f77_zmq
  implicit none
  
  integer(ZMQ_PTR)    ,intent(in)     :: zmq_socket_pull
  integer             ,intent(out)    :: task_id
  integer             ,intent(out)    :: blockb, blocke
  integer             ,intent(out)    :: N
  integer             ,intent(out)    :: idx(dav_size)
  double precision    ,intent(out)    :: vt(N_states_diag, dav_size)
  double precision    ,intent(out)    :: st(N_states_diag, dav_size)

  integer                            :: rc

  rc = f77_zmq_recv( zmq_socket_pull, blockb, 4, 0)
  if(rc /= 4) stop "davidson_push_results failed to pull blockb"
  
  rc = f77_zmq_recv( zmq_socket_pull, blocke, 4, 0)
  if(rc /= 4) stop "davidson_push_results failed to pull blocke"
  
  rc = f77_zmq_recv( zmq_socket_pull, N, 4, 0)
  if(rc /= 4) stop "davidson_push_results failed to pull N"
  
  rc = f77_zmq_recv( zmq_socket_pull, idx, 4*N, 0)
  if(rc /= 4*N) stop "davidson_push_results failed to pull idx"

  rc = f77_zmq_recv( zmq_socket_pull, vt, 8*N_states_diag* N, 0)
  if(rc /= 8*N_states_diag* N) stop "davidson_push_results failed to pull vt"

  rc = f77_zmq_recv( zmq_socket_pull, st, 8*N_states_diag* N, 0)
  if(rc /= 8*N_states_diag* N) stop "davidson_push_results failed to pull st"

  rc = f77_zmq_recv( zmq_socket_pull, task_id, 4, 0)
  if(rc /= 4) stop "davidson_pull_results failed to pull task_id"
end subroutine



subroutine davidson_collector(zmq_to_qp_run_socket, zmq_socket_pull , v0, s0)
  use f77_zmq
  implicit none

  integer(ZMQ_PTR), intent(in)     :: zmq_to_qp_run_socket
  integer(ZMQ_PTR), intent(in)     :: zmq_socket_pull

  double precision    ,intent(inout)    :: v0(dav_size, N_states_diag)
  double precision    ,intent(inout)    :: s0(dav_size, N_states_diag)

  integer                          :: more, task_id
  

  integer                  :: blockb, blocke
  integer                  :: N
  integer             , allocatable      :: idx(:)
  double precision    , allocatable      :: vt(:,:), v0t(:,:), s0t(:,:)
  double precision    , allocatable      :: st(:,:)
  
  allocate(idx(dav_size)) 
  allocate(vt(N_states_diag, dav_size)) 
  allocate(st(N_states_diag, dav_size)) 
  allocate(v0t(N_states_diag, dav_size)) 
  allocate(s0t(N_states_diag, dav_size)) 
  
  v0t = 00.d0
  s0t = 00.d0

  more = 1
  
  do while (more == 1)
    call davidson_pull_results(zmq_socket_pull, blockb, blocke, N, idx, vt, st, task_id)
    !DIR$ FORCEINLINE
    call davidson_collect(blockb, blocke, N, idx, vt, st , v0t, s0t)
    call zmq_delete_task(zmq_to_qp_run_socket,zmq_socket_pull,task_id,more)
  end do
  deallocate(idx,vt,st)

  call dtranspose(v0t,size(v0t,1), v0, size(v0,1), N_states_diag, dav_size)
  call dtranspose(s0t,size(s0t,1), s0, size(s0,1), N_states_diag, dav_size)
  deallocate(v0t,s0t)
end subroutine


subroutine davidson_run(zmq_to_qp_run_socket , v0, s0)
  use f77_zmq
  implicit none
  
  integer(ZMQ_PTR), intent(in) :: zmq_to_qp_run_socket
  integer(ZMQ_PTR),external      :: new_zmq_to_qp_run_socket
  integer(ZMQ_PTR)             :: zmq_collector
  integer(ZMQ_PTR), external     :: new_zmq_pull_socket
  integer(ZMQ_PTR)               :: zmq_socket_pull
  
  integer :: i
  integer, external              :: omp_get_thread_num

  double precision    , intent(inout)     :: v0(dav_size, N_states_diag)
  double precision    , intent(inout)     :: s0(dav_size, N_states_diag)
  
  call zmq_set_running(zmq_to_qp_run_socket)
  
  zmq_collector = new_zmq_to_qp_run_socket()
  zmq_socket_pull = new_zmq_pull_socket()
  i = omp_get_thread_num()
  
  
  PROVIDE nproc
  
  !$OMP PARALLEL DEFAULT(shared) private(i) num_threads(nproc+2)
    i = omp_get_thread_num()
    if (i==0) then
      call davidson_collector(zmq_collector, zmq_socket_pull , v0, s0)
      call end_zmq_to_qp_run_socket(zmq_collector)
      call end_zmq_pull_socket(zmq_socket_pull)
      call davidson_miniserver_end()
    else if(i==1) then
      call davidson_miniserver_run()
    else
      call davidson_slave_inproc(i)
    endif
  !$OMP END PARALLEL
  call end_parallel_job(zmq_to_qp_run_socket, 'davidson')
end subroutine



subroutine davidson_miniserver_run()
  use f77_zmq
  implicit none
  integer(ZMQ_PTR)        responder
  character*(64)          address
  character(len=:), allocatable :: buffer
  integer                 rc
  
  allocate (character(len=20) :: buffer)
  address = 'tcp://*:11223'
  
  responder = f77_zmq_socket(zmq_context, ZMQ_REP)
  rc        = f77_zmq_bind(responder,address)
  
  do
    rc = f77_zmq_recv(responder, buffer, 5, 0)
    if (buffer(1:rc) /= 'end') then
      rc = f77_zmq_send (responder, dav_size, 4, ZMQ_SNDMORE)
      rc = f77_zmq_send (responder, dav_det, 16*N_int*dav_size, ZMQ_SNDMORE)
      rc = f77_zmq_send (responder, dav_ut, 8*dav_size*N_states_diag, 0)
    else
      rc = f77_zmq_send (responder, "end", 3, 0)
      exit
    endif
  enddo

  rc = f77_zmq_close(responder)
end subroutine


subroutine davidson_miniserver_end()
  implicit none
  use f77_zmq
  
  integer(ZMQ_PTR)        requester
  character*(64)          address
  integer                 rc
  character*(64)          buf
  
  address = trim(qp_run_address)//':11223'
  requester = f77_zmq_socket(zmq_context, ZMQ_REQ)
  rc        = f77_zmq_connect(requester,address)

  rc = f77_zmq_send(requester, "end", 3, 0)
  rc = f77_zmq_recv(requester, buf, 3, 0)
  rc = f77_zmq_close(requester)
end subroutine

  
subroutine davidson_miniserver_get()
  implicit none
  use f77_zmq
  
  integer(ZMQ_PTR)        requester
  character*(64)          address
  character*(20)          buffer
  integer                 rc
  
  address = trim(qp_run_address)//':11223'
  
  requester = f77_zmq_socket(zmq_context, ZMQ_REQ)
  rc        = f77_zmq_connect(requester,address)

  rc = f77_zmq_send(requester, "Hello", 5, 0)
  rc = f77_zmq_recv(requester, dav_size, 4, 0)
  TOUCH dav_size
  rc = f77_zmq_recv(requester, dav_det, 16*N_int*dav_size, 0)
  rc = f77_zmq_recv(requester, dav_ut, 8*dav_size*N_states_diag, 0)
  TOUCH dav_det dav_ut

  
end subroutine



BEGIN_PROVIDER [ integer(bit_kind), dav_det, (N_int, 2, dav_size) ]
END_PROVIDER


BEGIN_PROVIDER [ double precision, dav_ut, (N_states_diag, dav_size) ]
END_PROVIDER


BEGIN_PROVIDER [ integer, dav_size ]
END_PROVIDER


 BEGIN_PROVIDER [ integer, shortcut_, (0:dav_size+1, 2) ]
&BEGIN_PROVIDER [ integer(bit_kind), version_, (N_int, dav_size, 2) ]
&BEGIN_PROVIDER [ integer(bit_kind), sorted_, (N_int, dav_size, 2) ]
&BEGIN_PROVIDER [ integer, sort_idx_, (dav_size, 2) ]
  implicit none
  call sort_dets_ab_v(dav_det, sorted_(1,1,1), sort_idx_(1,1), shortcut_(0,1), version_(1,1,1), dav_size, N_int)
  call sort_dets_ba_v(dav_det, sorted_(1,1,2), sort_idx_(1,2), shortcut_(0,2), version_(1,1,2), dav_size, N_int)
END_PROVIDER