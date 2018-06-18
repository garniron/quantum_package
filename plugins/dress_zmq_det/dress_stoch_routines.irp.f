BEGIN_PROVIDER [ integer, fragment_first ]
  implicit none
  fragment_first = first_det_of_teeth(1)
END_PROVIDER


subroutine ZMQ_dress(E, dress, delta_out, delta_s2_out, relative_error, lndet)
  use f77_zmq
  
  implicit none
  
  integer, intent(in) :: lndet
  character(len=64000)           :: task
  character(len=3200)            :: temp 
  integer(ZMQ_PTR)               :: zmq_to_qp_run_socket, zmq_socket_pull
  integer, external              :: omp_get_thread_num
  double precision, intent(in)   :: E(N_states), relative_error
  double precision, intent(out)  :: dress(N_states)
  double precision, intent(out)  :: delta_out(N_states, N_det)
  double precision, intent(out)  :: delta_s2_out(N_states, N_det)
  
  double precision, allocatable  :: delta(:,:)
  double precision, allocatable  :: delta_s2(:,:)
  
  integer                        :: i, j, k, Ncp
  
  double precision, external     :: omp_get_wtime
  double precision               :: time
  integer, external              :: add_task_to_taskserver
  double precision               :: state_average_weight_save(N_states)
  task(:) = CHAR(0)
  temp(:) = CHAR(0)
  allocate(delta(N_states,N_det), delta_s2(N_states, N_det))
  state_average_weight_save(:) = state_average_weight(:)
  do dress_stoch_istate=1,N_states
    SOFT_TOUCH dress_stoch_istate
    state_average_weight(:) = 0.d0
    state_average_weight(dress_stoch_istate) = 1.d0
    TOUCH state_average_weight
    
    !provide psi_coef_generators
    provide nproc fragment_first fragment_count mo_bielec_integrals_in_map mo_mono_elec_integral dress_weight psi_selectors
    !print *, dress_e0_denominator
    
    print *, '========== ================= ================= ================='
    print *, ' Samples        Energy         Stat. Error         Seconds      '
    print *, '========== ================= ================= ================='
   
    call new_parallel_job(zmq_to_qp_run_socket,zmq_socket_pull, 'dress')
    
    integer, external              :: zmq_put_psi
    integer, external              :: zmq_put_N_det_generators
    integer, external              :: zmq_put_N_det_selectors
    integer, external              :: zmq_put_dvector
    integer, external              :: zmq_set_running

    if (zmq_put_psi(zmq_to_qp_run_socket,1) == -1) then
      stop 'Unable to put psi on ZMQ server'
    endif
    if (zmq_put_N_det_generators(zmq_to_qp_run_socket, 1) == -1) then
      stop 'Unable to put N_det_generators on ZMQ server'
    endif
    if (zmq_put_N_det_selectors(zmq_to_qp_run_socket, 1) == -1) then
      stop 'Unable to put N_det_selectors on ZMQ server'
    endif
    if (zmq_put_dvector(zmq_to_qp_run_socket,1,'energy',dress_e0_denominator,size(dress_e0_denominator)) == -1) then
      stop 'Unable to put energy on ZMQ server'
    endif
    if (zmq_put_dvector(zmq_to_qp_run_socket,1,"state_average_weight",state_average_weight,N_states) == -1) then
      stop 'Unable to put state_average_weight on ZMQ server'
    endif
    if (zmq_put_dvector(zmq_to_qp_run_socket,1,"dress_stoch_istate",real(dress_stoch_istate,8),1) == -1) then
      stop 'Unable to put dress_stoch_istate on ZMQ server'
    endif


    integer(ZMQ_PTR), external     :: new_zmq_to_qp_run_socket
    integer                        :: ipos, sz
    integer                        :: block(1), block_i, cur_tooth_reduce, ntas
    logical                        :: flushme
    block = 0
    block_i = 0
    cur_tooth_reduce = 0
    ipos=1
    ntas = 0
    do i=1,N_dress_jobs+1
      flushme = (i==N_dress_jobs+1 .or. block_i == size(block) .or. block_i >=cur_tooth_reduce )
      if(.not. flushme) flushme = (tooth_reduce(dress_jobs(i)) == 0 .or. tooth_reduce(dress_jobs(i)) /= cur_tooth_reduce)
      
      if(flushme .and. block_i > 0) then
        if(block(1) > fragment_first) then
          ntas += 1
          write(temp, '(I9,1X,60(I9,1X))') 0, block(:block_i)
          sz = len(trim(temp))+1
          temp(sz:sz) = '|'
          !write(task(ipos:ipos+20),'(I9,1X,I9,''|'')') 0, dress_jobs(i)
          write(task(ipos:ipos+sz), *) temp(:sz)
          !ipos += 20
          ipos += sz+1
          if (ipos > 63000 .or. i==N_dress_jobs+1) then
            if (add_task_to_taskserver(zmq_to_qp_run_socket,trim(task(1:ipos))) == -1) then
              stop 'Unable to add task to task server'
            endif
          
            ipos=1
          endif
        else
          if(block_i /= 1) stop "reduced fragmented dets"
          do j=1,fragment_count
            ntas += 1
            write(task(ipos:ipos+20),'(I9,1X,I9,''|'')') j, block(1)
            ipos += 20
            if (ipos > 63000 .or. i==N_dress_jobs+1) then
              ntas += 1
              if (add_task_to_taskserver(zmq_to_qp_run_socket,trim(task(1:ipos))) == -1) then
                stop 'Unable to add task to task server'
              endif
              ipos=1
            endif
          end do
        end if
        block_i = 0
        block = 0
      end if
      
      if(i /= N_dress_jobs+1) then
        cur_tooth_reduce = tooth_reduce(dress_jobs(i))
        block_i += 1
        block(block_i) = dress_jobs(i)
      end if
    end do
    if (zmq_set_running(zmq_to_qp_run_socket) == -1) then
      print *,  irp_here, ': Failed in zmq_set_running'
    endif
    
    !!$OMP PARALLEL DEFAULT(shared) NUM_THREADS(nproc)              &
    !    !$OMP  PRIVATE(i)
    !i = omp_get_thread_num()
    !if (i==0) then
      call dress_collector(zmq_socket_pull,E, relative_error, delta, delta_s2, dress,&
         dress_stoch_istate)
    !else
    !  call dress_slave_inproc(i)
    !endif
    !!$OMP END PARALLEL
    delta_out(dress_stoch_istate,1:N_det) = delta(dress_stoch_istate,1:N_det)
    delta_s2_out(dress_stoch_istate,1:N_det) = delta_s2(dress_stoch_istate,1:N_det)
    call end_parallel_job(zmq_to_qp_run_socket, zmq_socket_pull, 'dress')
    
    print *, '========== ================= ================= ================='
  enddo
  FREE dress_stoch_istate
  state_average_weight(:) = state_average_weight_save(:)
  TOUCH state_average_weight
  deallocate(delta,delta_s2)
  
end subroutine


subroutine dress_slave_inproc(i)
  implicit none
  integer, intent(in)            :: i
  
  call run_dress_slave(1,i,dress_e0_denominator)
end



subroutine dress_collector(zmq_socket_pull, E, relative_error, delta, delta_s2, dress, istate)
  use f77_zmq
  use bitmasks
  implicit none

  
  integer(ZMQ_PTR), intent(in)   :: zmq_socket_pull
  integer, intent(in)            :: istate

  double precision, intent(in)   :: relative_error, E(N_states)
  double precision, intent(out)  :: dress(N_states)
  double precision, allocatable  :: cp(:,:,:,:)

  double precision, intent(out)  :: delta(N_states, N_det)
  double precision, intent(out)  :: delta_s2(N_states, N_det)
  double precision, allocatable  :: delta_loc(:,:,:)
  double precision, allocatable  :: dress_detail(:,:)
  double precision               :: dress_mwen(N_states)
  integer(ZMQ_PTR),external      :: new_zmq_to_qp_run_socket
  integer(ZMQ_PTR)               :: zmq_to_qp_run_socket

  integer(ZMQ_PTR), external     :: new_zmq_pull_socket

  integer :: more
  integer :: i, j, k, i_state, N
  integer :: task_id, i_gen
  double precision, save :: time0 = -1.d0
  double precision :: time
  double precision, external :: omp_get_wtime
  integer :: last_cp
  integer :: delta_loc_cur, is, N_buf(3)
  integer, allocatable :: int_buf(:), agreg_for_cp(:)
  double precision, allocatable :: double_buf(:)
  integer(bit_kind), allocatable :: det_buf(:,:,:)
  integer, external :: zmq_delete_tasks
  print *, "MODULE HIDEUX DRESS_ZMQ_DET"

  last_cp = 10000000
  allocate(agreg_for_cp(N_cp))
  agreg_for_cp = 0
  allocate(int_buf(N_dress_int_buffer), double_buf(N_dress_double_buffer), det_buf(N_int,2,N_dress_det_buffer))
  delta_loc_cur = 1

  delta = 0d0
  delta_s2 = 0d0
  allocate(cp(N_states, N_det, N_cp, 2), dress_detail(N_states, N_det))
  allocate(delta_loc(N_states, N_det, 2))
  dress_detail = 0d0 ! -1000d0
  cp = 0d0
  character*(512) :: task

  zmq_to_qp_run_socket = new_zmq_to_qp_run_socket()
  more = 1
  if (time0 < 0.d0) then
      call wall_time(time0)
  endif
  logical :: loop, floop

  floop = .true.
  loop = .true.

  pullLoop : do while (loop)
    call pull_dress_results(zmq_socket_pull, i_gen, delta_loc, int_buf, double_buf, det_buf, N_buf, task_id, dress_mwen)
    if(floop) then
      call wall_time(time)
      print *, "first_pull", time-time0
      time0 = time
      floop = .false.
    end if
    
    
    !if(cur_cp == -1) then
      call dress_pulled(i_gen, int_buf, double_buf, det_buf, N_buf) 
      if (zmq_delete_tasks(zmq_to_qp_run_socket,zmq_socket_pull,task_id,1,more) == -1) then
        stop 'Unable to delete tasks'
      endif
      if(more == 0)  loop = .false. !stop 'loop = .false.' !!!!!!!!!!!!!!!!
      dress_detail(:, i_gen) = dress_mwen(:)
      dress(istate) = E(istate)+sum(dress_detail(istate, :))


      !print *, "DETAIL", ind, dress_mwen
    !else if(cur_cp > 0) then
      !if(ind == 0) cycle
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,N_det
        delta(:,i) += delta_loc(:,i,1)
        !cp(:,i,cur_cp,1) += delta_loc(:,i,1)
      end do

      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,N_det
        delta_s2(:,i) += delta_loc(:,i,2)
        !cp(:,i,cur_cp,2) += delta_loc(:,i,2)
      end do
      !$OMP END PARALLEL DO
    !end if
  end do pullLoop

  dress(istate) = E(istate)+sum(dress_detail(istate, :))
  print *, "DRESSING ENERGY", dress(istate)
  !dress(istate) = E(istate)+E0+avg
  call end_zmq_to_qp_run_socket(zmq_to_qp_run_socket)
end subroutine


integer function dress_find(v, w, sze, imin, imax)
  implicit none
  integer, intent(in) :: sze, imin, imax
  double precision, intent(in) :: v, w(sze)
  integer :: i,l,h
  integer, parameter :: block=64

  l = imin
  h = imax-1

  do while(h-l >= block)
    i = ishft(h+l,-1)
    if(w(i+1) > v) then
      h = i-1
    else
      l = i+1
    end if
  end do
  !DIR$ LOOP COUNT (64)
  do dress_find=l,h
    if(w(dress_find) >= v) then
      exit
    end if
  end do
end function


 BEGIN_PROVIDER [ integer, gen_per_cp ]
&BEGIN_PROVIDER [ integer, comb_teeth ]
&BEGIN_PROVIDER [ integer, N_cps_max ]
  implicit none
  BEGIN_DOC
! N_cps_max : max number of checkpoints
!
! gen_per_cp : number of generators per checkpoint
  END_DOC
  comb_teeth = 64
  N_cps_max = 16
  gen_per_cp = (N_det_generators / N_cps_max) + 1
END_PROVIDER


 BEGIN_PROVIDER [ integer, N_cp ]
&BEGIN_PROVIDER [ double precision, cps_N, (N_cps_max) ]
&BEGIN_PROVIDER [ integer, cp_first_tooth, (N_cps_max) ]
&BEGIN_PROVIDER [ integer, done_cp_at, (N_det_generators) ]
&BEGIN_PROVIDER [ integer, done_cp_at_det, (N_det_generators) ]
&BEGIN_PROVIDER [ integer, needed_by_cp, (0:N_cps_max) ]
&BEGIN_PROVIDER [ double precision, cps, (N_det_generators, N_cps_max) ]
&BEGIN_PROVIDER [ integer, N_dress_jobs ]
&BEGIN_PROVIDER [ integer, dress_jobs, (N_det_generators) ]
&BEGIN_PROVIDER [ double precision, comb, (N_det_generators) ]
&BEGIN_PROVIDER [ integer, tooth_reduce, (N_det_generators) ]
  implicit none
  logical, allocatable         :: computed(:), comp_filler(:)
  integer                        :: i, j, last_full, dets(comb_teeth)
  
  integer                        :: k, l, cur_cp, under_det(comb_teeth+1)
  integer :: cp_limit(N_cps_max)
  integer, allocatable :: iorder(:), first_cp(:)
  integer, allocatable :: filler(:)
  integer :: nfiller, lfiller, cfiller
  logical :: fracted
   
  integer :: first_suspect
  provide psi_coef_generators
  first_suspect = 1

  allocate(filler(n_det_generators))
  allocate(iorder(N_det_generators), first_cp(N_cps_max+1))
  allocate(computed(N_det_generators))
  allocate(comp_filler(N_det_generators))
  first_cp = 1
  cps = 0d0
  cur_cp = 1
  done_cp_at = 0
  done_cp_at_det = 0
  needed_by_cp = 0
  comp_filler = .false.
  computed = .false.
  cps_N = 1d0
  tooth_reduce = 0
  
  integer :: fragsize
  fragsize = N_det_generators / ((N_cps_max-1+1)*(N_cps_max-1+2)/2)

  do i=1,N_cps_max
    cp_limit(i) = fragsize * i * (i+1) / 2
  end do
  cp_limit(N_cps_max) = N_det*2

  N_dress_jobs = first_det_of_comb - 1
  do i=1, N_dress_jobs
    dress_jobs(i) = i
    computed(i) = .true.
  end do
  
  l=first_det_of_comb
  call random_seed(put=(/321,654,65,321,65,321,654,65,321,6321,654,65,321,621,654,65,321,65,654,65,321,65/))
  call RANDOM_NUMBER(comb)
  lfiller = 1
  nfiller = 1
  do i=1,N_det_generators
    !print *, i, N_dress_jobs
    comb(i) = comb(i) * comb_step
    !DIR$ FORCEINLINE
    call add_comb(comb(i), computed, cps(1, cur_cp), N_dress_jobs, dress_jobs)
    
    !if(N_dress_jobs / gen_per_cp > (cur_cp-1) .or. N_dress_jobs == N_det_generators) then
    if(N_dress_jobs > cp_limit(cur_cp) .or. N_dress_jobs == N_det_generators) then
      first_cp(cur_cp+1) = N_dress_jobs
      done_cp_at(N_dress_jobs) = cur_cp
      cps_N(cur_cp) = dfloat(i)
      if(N_dress_jobs /= N_det_generators) then
        cps(:, cur_cp+1) = cps(:,  cur_cp)
        cur_cp += 1
      end if
      
      if (N_dress_jobs == N_det_generators) then
        exit
      end if
    end if

    !!!!!!!!!!!!!!!!!!!!!!!!
    if(.TRUE.) then 
    do l=first_suspect,N_det_generators
      if((.not. computed(l))) then
        N_dress_jobs+=1
        dress_jobs(N_dress_jobs) = l
        computed(l) = .true.
        first_suspect = l
        exit
      end if
    end do
    
    if (N_dress_jobs == N_det_generators) exit

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ELSE
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do l=first_suspect,N_det_generators
      if((.not. computed(l)) .and. (.not. comp_filler(l))) exit
    end do
    first_suspect = l
    if(l > N_det_generators) cycle
    
    cfiller = tooth_of_det(l)-1
    if(cfiller > lfiller) then
      do j=1,nfiller-1
        if(.not. computed(filler(j))) then
          k=N_dress_jobs+1
          dress_jobs(k) = filler(j)
          N_dress_jobs = k
        end if  
        computed(filler(j)) = .true.
      end do
      nfiller = 2
      filler(1) = l
      lfiller = cfiller
    else
      filler(nfiller) = l
      nfiller += 1
    end if
    comp_filler(l) = .True.
    end if
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  enddo

  
  do j=1,nfiller-1
     if(.not. computed(filler(j)))then
       k=N_dress_jobs+1
       dress_jobs(k) = filler(j)
       N_dress_jobs = k
     end if
     computed(filler(j)) = .true.
  end do


  N_cp = cur_cp
  
  if(N_dress_jobs /= N_det_generators .or. N_cp > N_cps_max) then
    print *, N_dress_jobs, N_det_generators, N_cp, N_cps_max
    stop "error in jobs creation"
  end if

  cur_cp = 0
  do i=1,N_dress_jobs
    if(done_cp_at(i) /= 0) cur_cp = done_cp_at(i)
    done_cp_at(i) = cur_cp
    done_cp_at_det(dress_jobs(i)) = cur_cp
    needed_by_cp(cur_cp) += 1
  end do
  

  under_det = 0
  cp_first_tooth = 0
  do i=1,N_dress_jobs
    do j=comb_teeth+1,1,-1
      if(dress_jobs(i) <= first_det_of_teeth(j)) then
        under_det(j) = under_det(j) + 1
        if(under_det(j) == first_det_of_teeth(j))then
          do l=done_cp_at(i)+1, N_cp
            cps(:first_det_of_teeth(j)-1, l) = 0d0
            cp_first_tooth(l) = j
          end do
          cps(first_det_of_teeth(j), done_cp_at(i)+1) = &
              cps(first_det_of_teeth(j), done_cp_at(i)+1) * fractage(j)
        end if
      else
        exit
      end if
    end do
  end do
  cp_first_tooth(N_cp) = comb_teeth+1
  
  do i=1,N_det_generators
    do j=N_cp,2,-1
      cps(i,j) -= cps(i,j-1)
    end do
  end do

  iorder = -1

  cps(:, N_cp) = 0d0 

  iloop : do i=fragment_first+1,N_det_generators
    k = tooth_of_det(i)
    if(k == 0) cycle
    if (i == first_det_of_teeth(k)) cycle
    
    do j=1,N_cp
     if(cps(i, j) /= 0d0) cycle iloop
    end do
    
    tooth_reduce(i) = k
  end do iloop
  
  do i=1,N_det_generators
    if(tooth_reduce(dress_jobs(i)) == 0) dress_jobs(i) = dress_jobs(i)+N_det*2
  end do

  do i=1,N_cp-1
    call isort(dress_jobs(first_cp(i)+1),iorder,first_cp(i+1)-first_cp(i)-1)
  end do
 
 do i=1,N_det_generators
   if(dress_jobs(i) > N_det) dress_jobs(i) = dress_jobs(i) - N_det*2
 end do
END_PROVIDER


subroutine get_comb_val(stato, detail, cur_cp, val, istate)
  implicit none
  integer, intent(in)           :: cur_cp, istate
  integer                       :: first
  double precision, intent(in)  :: stato, detail(N_states, N_det_generators)
  double precision, intent(out) :: val
  double precision :: curs
  integer :: j, k
  integer, external :: dress_find

  curs = 1d0 - stato
  val = 0d0
  first = cp_first_tooth(cur_cp) 

  do j = comb_teeth, first, -1
    !DIR$ FORCEINLINE
    k = dress_find(curs, dress_cweight,size(dress_cweight), first_det_of_teeth(j), first_det_of_teeth(j+1))
    if(k == first_det_of_teeth(first)) then
     val += detail(istate, k) * dress_weight_inv(k) * comb_step * fractage(first)
    else
     val += detail(istate, k) * dress_weight_inv(k) * comb_step
    end if

    curs -= comb_step
  end do
end subroutine


subroutine get_comb(stato, dets)
  implicit none
  double precision, intent(in) :: stato
  integer, intent(out) :: dets(comb_teeth)
  double precision :: curs
  integer :: j
  integer, external :: dress_find

  curs = 1d0 - stato
  do j = comb_teeth, 1, -1
    !DIR$ FORCEINLINE
    dets(j) = dress_find(curs, dress_cweight,size(dress_cweight), first_det_of_teeth(j), first_det_of_teeth(j+1))
    curs -= comb_step
  end do
end subroutine


subroutine add_comb(com, computed, cp, N, tbc)
  implicit none
  double precision, intent(in) :: com
  integer, intent(inout) :: N
  double precision, intent(inout) :: cp(N_det)
  logical, intent(inout) :: computed(N_det_generators)
  integer, intent(inout) :: tbc(N_det_generators)
  integer :: i, k, l, dets(comb_teeth)
  
  !DIR$ FORCEINLINE
  call get_comb(com, dets)
  k=N+1
  do i = 1, comb_teeth
    l = dets(i)
    cp(l) += 1d0
    if(.not.(computed(l))) then
      tbc(k) = l
      k = k+1
      computed(l) = .true.
    end if
  end do
  N = k-1
end subroutine


BEGIN_PROVIDER [ integer, dress_stoch_istate ]
  implicit none
  dress_stoch_istate = 1
END_PROVIDER

 
 BEGIN_PROVIDER [ double precision, dress_weight, (N_det_generators) ]
&BEGIN_PROVIDER [ double precision, dress_weight_inv, (N_det_generators) ]
&BEGIN_PROVIDER [ double precision, dress_cweight, (N_det_generators) ]
&BEGIN_PROVIDER [ double precision, dress_cweight_cache, (N_det_generators) ]
&BEGIN_PROVIDER [ double precision, fractage, (comb_teeth) ]
&BEGIN_PROVIDER [ double precision, comb_step ]
&BEGIN_PROVIDER [ integer, first_det_of_teeth, (comb_teeth+1) ]
&BEGIN_PROVIDER [ integer, first_det_of_comb ]
&BEGIN_PROVIDER [ integer, tooth_of_det, (N_det_generators) ]
  implicit none
  integer :: i
  double precision :: norm_left, stato
  integer, external :: dress_find  

  dress_weight(1)  = psi_coef_generators(1,dress_stoch_istate)**2
  dress_cweight(1) = psi_coef_generators(1,dress_stoch_istate)**2
  
  do i=1,N_det_generators
    dress_weight(i) = psi_coef_generators(i,dress_stoch_istate)**2
  enddo

  ! Important to loop backwards for numerical precision
  dress_cweight(N_det_generators) = dress_weight(N_det_generators)
  do i=N_det_generators-1,1,-1
    dress_cweight(i) = dress_weight(i) + dress_cweight(i+1) 
  end do
  
  do i=1,N_det_generators
    dress_weight(i)  = dress_weight(i) / dress_cweight(1)
    dress_cweight(i) = dress_cweight(i) / dress_cweight(1)
  enddo

  do i=1,N_det_generators-1
    dress_cweight(i) = 1.d0 - dress_cweight(i+1) 
  end do
  dress_cweight(N_det_generators) = 1.d0
  
  norm_left = 1d0
  
  comb_step = 1d0/dfloat(comb_teeth)
  !print *, "comb_step", comb_step
  first_det_of_comb = 1
  do i=1,N_det_generators ! min(100,N_det_generators)
    first_det_of_comb = i
    if(dress_weight(i)/norm_left < comb_step) then
      exit
    end if
    norm_left -= dress_weight(i)
  end do
  first_det_of_comb = max(2,first_det_of_comb)
  call write_int(6, first_det_of_comb-1, 'Size of deterministic set')
  

  comb_step =  (1d0 - dress_cweight(first_det_of_comb-1)) * comb_step
  
  stato = 1d0 - comb_step
  iloc = N_det_generators
  do i=comb_teeth, 1, -1
    integer :: iloc
    iloc = dress_find(stato, dress_cweight, N_det_generators, 1, iloc)
    first_det_of_teeth(i) = iloc
    fractage(i) = (dress_cweight(iloc) - stato) / dress_weight(iloc)
    stato -= comb_step
  end do
  first_det_of_teeth(comb_teeth+1) = N_det_generators + 1
  first_det_of_teeth(1) = first_det_of_comb
  

  if(first_det_of_teeth(1) /= first_det_of_comb) then
     print *, 'Error in ', irp_here
     stop "comb provider"
  endif
  
  do i=1,N_det_generators
    dress_weight_inv(i) = 1.d0/dress_weight(i)
  enddo

  tooth_of_det(:first_det_of_teeth(1)-1) = 0
  do i=1,comb_teeth
    tooth_of_det(first_det_of_teeth(i):first_det_of_teeth(i+1)-1) = i
  end do
END_PROVIDER







