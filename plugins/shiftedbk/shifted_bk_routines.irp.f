use selection_types


 BEGIN_PROVIDER [ double precision, global_sum_alpha2, (N_states) ]
&BEGIN_PROVIDER [ double precision, slave_sum_alpha2, (N_states, Nproc) ]
  global_sum_alpha2 = 0d0
  slave_sum_alpha2 = 0d0 
END_PROVIDER


 BEGIN_PROVIDER [ double precision, fock_diag_tmp_, (2,mo_tot_num+1,Nproc) ]
&BEGIN_PROVIDER [ integer, n_det_add ]
&BEGIN_PROVIDER [ double precision, a_h_i, (N_det, Nproc) ]
&BEGIN_PROVIDER [ double precision, a_s2_i, (N_det, Nproc) ]
&BEGIN_PROVIDER [ type(selection_buffer), sb, (Nproc) ]
&BEGIN_PROVIDER [ type(selection_buffer), global_sb ]
&BEGIN_PROVIDER [ type(selection_buffer), mini_sb ]
&BEGIN_PROVIDER [ double precision, N_det_increase_factor ]
   implicit none
   fock_diag_tmp_(:,:,:) = 0.d0
   integer                        :: i
   
   N_det_increase_factor = dble(N_states)
   
   
   n_det_add = max(1, int(float(N_det) * N_det_increase_factor))
   call create_selection_buffer(n_det_add, n_det_add*2, global_sb)
   call create_selection_buffer(n_det_add, n_det_add*2, mini_sb)
   do i=1,Nproc
     call create_selection_buffer(n_det_add, n_det_add*2, sb(i))
   end do
   a_h_i = 0d0
   a_s2_i = 0d0
END_PROVIDER


 BEGIN_PROVIDER [ integer, N_dress_int_buffer ]
&BEGIN_PROVIDER [ integer, N_dress_double_buffer ]
&BEGIN_PROVIDER [ integer, N_dress_det_buffer ]
  implicit none
  N_dress_int_buffer = 1
  N_dress_double_buffer = n_det_add+N_states
  N_dress_det_buffer = n_det_add
END_PROVIDER


subroutine generator_done(i_gen, int_buf, double_buf, det_buf, N_buf, iproc)
  implicit none
  integer, intent(in) :: i_gen, iproc
  integer, intent(out) :: int_buf(N_dress_int_buffer), N_buf(3)
  double precision, intent(out) :: double_buf(N_dress_double_buffer)
  integer(bit_kind), intent(out) :: det_buf(N_int, 2, N_dress_det_buffer)
  integer :: i
  int_buf(:) = 0
  
    call sort_selection_buffer(sb(iproc))
   
    if(sb(iproc)%cur > 0) then
      !$OMP CRITICAL
      call merge_selection_buffers(sb(iproc), mini_sb)
      !call sort_selection_buffer(mini_sb)
      do i=1,Nproc
        mini_sb%mini = min(sb(i)%mini, mini_sb%mini)
      end do
      do i=1,Nproc
        sb(i)%mini = mini_sb%mini
      end do
      !$OMP END CRITICAL
    end if
      call truncate_to_mini(sb(iproc))
      det_buf(:,:,:sb(iproc)%cur) = sb(iproc)%det(:,:,:sb(iproc)%cur)
      double_buf(:sb(iproc)%cur) = sb(iproc)%val(:sb(iproc)%cur)
      double_buf(sb(iproc)%cur+1:sb(iproc)%cur+N_states) = slave_sum_alpha2(:,iproc)
      N_buf(1) = 1
      N_buf(2) = sb(iproc)%cur+N_states
      N_buf(3) = sb(iproc)%cur
    
    sb(iproc)%cur = 0
    slave_sum_alpha2(:,iproc) = 0d0
end subroutine


subroutine generator_start(i_gen, iproc)
  implicit none
  integer, intent(in) :: i_gen, iproc
  integer :: i
  
  call build_fock_tmp(fock_diag_tmp_(1,1,iproc),psi_det_generators(1,1,i_gen),N_int)
end subroutine


subroutine dress_pulled(ind, int_buf, double_buf, det_buf, N_buf)
  use bitmasks
  implicit none
  
  integer, intent(in) :: ind, N_buf(3)
  integer, intent(in) :: int_buf(*)
  double precision, intent(in) :: double_buf(*)
  integer(bit_kind), intent(in) :: det_buf(N_int,2,*)
  integer :: i
  
  do i=1,N_buf(3)
    call add_to_selection_buffer(global_sb, det_buf(1,1,i), double_buf(i))
  end do
  if(N_buf(3) + N_states /= N_buf(2)) stop "buf size"
  !$OMP CRITICAL
  global_sum_alpha2(:) += double_buf(N_buf(3)+1:N_buf(2))
  !$OMP END CRITICAL
end subroutine


subroutine delta_ij_done()
  use bitmasks
  implicit none
  integer :: i, old_det_gen
  integer(bit_kind), allocatable :: old_generators(:,:,:)
  
  allocate(old_generators(N_int, 2, N_det_generators))
  old_generators(:,:,:) = psi_det_generators(:,:,:N_det_generators)
  old_det_gen = N_det_generators
 

  ! Add buffer only when the last state is computed
  call unique_selection_buffer(global_sb)
  call sort_selection_buffer(global_sb)
  call fill_H_apply_buffer_no_selection(global_sb%cur,global_sb%det,N_int,0) 
  call copy_H_apply_buffer_to_wf()
  if (s2_eig.or.(N_states > 1) ) then
    call make_s2_eigenfunction
  endif
  call undress_with_alpha(old_generators, old_det_gen, psi_det(1,1,N_det_delta_ij+1), N_det-N_det_delta_ij)
  call save_wavefunction

end subroutine


subroutine undress_with_alpha(old_generators, old_det_gen, alpha, n_alpha)
  use bitmasks
  implicit none

  integer(bit_kind), intent(in)   :: alpha(N_int,2,n_alpha)
  integer, intent(in)             :: n_alpha
  integer, allocatable :: minilist(:)
  integer(bit_kind), allocatable :: det_minilist(:,:,:)
  double precision, allocatable :: delta_ij_loc(:,:,:,:)
  integer :: exc(0:2,2,2), h1, h2, p1, p2, s1, s2
  integer :: i, j, k, ex, n_minilist, iproc, degree
  double precision :: haa, contrib, phase, c_alpha(N_states,Nproc), s_c_alpha(N_states)
  logical :: ok
  integer, external :: omp_get_thread_num
 
  integer,intent(in) :: old_det_gen
  integer(bit_kind), intent(in) :: old_generators(N_int, 2, old_det_gen)
  
  allocate(minilist(N_det_delta_ij), det_minilist(N_int, 2, N_det_delta_ij), delta_ij_loc(N_states, N_det_delta_ij, 2, Nproc))
  
  c_alpha = 0d0
  delta_ij_loc = 0d0
  
  !$OMP PARALLEL DO DEFAULT(SHARED) SCHEDULE(STATIC) PRIVATE(i, j, iproc, n_minilist, ex) &
  !$OMP PRIVATE(det_minilist, minilist, haa, contrib, s_c_alpha)  &
  !$OMP PRIVATE(exc, h1, h2, p1, p2, s1, s2, phase, degree, ok)
  do i=n_alpha,1,-1
    iproc = omp_get_thread_num()+1
    if(mod(i,10000) == 0) print *, "UNDRESSING", i, "/", n_alpha, iproc
    n_minilist = 0
    ok = .false.
    
    do j=1, old_det_gen
      call get_excitation_degree(alpha(1,1,i), old_generators(1,1,j), ex, N_int)
      if(ex <= 2) then
        call get_excitation(old_generators(1,1,j), alpha(1,1,i), exc,degree,phase,N_int)
        call decode_exc(exc,degree,h1,p1,h2,p2,s1,s2)
        ok = (mo_class(h1)(1:1) == 'A' .or. mo_class(h1)(1:1) == 'I') .and. &
           (mo_class(p1)(1:1) == 'A' .or. mo_class(p1)(1:1) == 'V') 
        if(ok .and. degree == 2) then
             ok = (mo_class(h2)(1:1) == 'A' .or. mo_class(h2)(1:1) == 'I') .and. &
             (mo_class(p2)(1:1) == 'A' .or. mo_class(p2)(1:1) == 'V') 
        end if
        if(ok) exit
      end if
    end do
    
    if(.not. ok) cycle

    do j=1, N_det_delta_ij
      call get_excitation_degree(alpha(1,1,i), psi_det(1,1,j), ex, N_int)
      if(ex <= 2) then
        n_minilist += 1
        det_minilist(:,:,n_minilist) = psi_det(:,:,j)
        minilist(n_minilist) = j
      end if
    end do
    call i_h_j(alpha(1,1,i), alpha(1,1,i), N_int, haa)
    call dress_with_alpha_(N_states, N_det_delta_ij, N_int, delta_ij_loc(1,1,1,iproc),  &
        minilist, det_minilist, n_minilist, alpha(1,1,i), haa, contrib, s_c_alpha, iproc)
    
    c_alpha(:,iproc) += s_c_alpha(:)**2
  end do
  !$OMP END PARALLEL DO
  
  do i=2,Nproc
    delta_ij_loc(:,:,:,1) += delta_ij_loc(:,:,:,i)
    c_alpha(:,1) += c_alpha(:,i)
  end do
  

  delta_ij_tmp(:,:,1) -=  delta_ij_loc(:,:,1,1)
  delta_ij_tmp(:,:,2) -=  delta_ij_loc(:,:,2,1) 
  
  !print *, "SUM ALPHA2  PRE", global_sum_alpha2
  !global_sum_alpha2(:) -= c_alpha(:,1)
  print *, "SUM C_ALPHA^2 =", global_sum_alpha2(:)
  !print *, "*** DRESSINS DIVIDED BY 1+SUM C_ALPHA^2 ***"
  !do i=1,N_states
  !  delta_ij_tmp(i,:,:) = delta_ij_tmp(i,:,:) / (1d0 + global_sum_alpha2(i))
  !end do
  global_sum_alpha2 = 0d0 
end subroutine


subroutine dress_with_alpha_(Nstates,Ndet,Nint,delta_ij_loc,minilist, det_minilist, n_minilist, alpha, haa, contrib, c_alpha, iproc)
  use bitmasks
  implicit none
  BEGIN_DOC
  !delta_ij_loc(:,:,1) : dressing column for H
  !delta_ij_loc(:,:,2) : dressing column for S2
  !minilist : indices of determinants connected to alpha ( in psi_det )
  !n_minilist : size of minilist
  !alpha : alpha determinant
  END_DOC
  integer, intent(in)             :: Nint, Ndet, Nstates, n_minilist, iproc
  integer(bit_kind), intent(in)   :: alpha(Nint,2), det_minilist(Nint, 2, n_minilist)
  integer,intent(in)              :: minilist(n_minilist)
  double precision, intent(inout) :: delta_ij_loc(Nstates,N_det,2)
  double precision, intent(out) :: contrib, c_alpha(N_states)
  double precision,intent(in) :: haa
  double precision :: hij, sij
  double precision, external :: diag_H_mat_elem_fock
  integer                        :: i,j,k,l,m, l_sd
  double precision :: hdress, sdress
  double precision :: de, a_h_psi(Nstates)
  

  a_h_psi = 0d0
  

  do l_sd=1,n_minilist
    call i_h_j_s2(alpha,det_minilist(1,1,l_sd),N_int,hij, sij)
    a_h_i(l_sd, iproc) = hij
    a_s2_i(l_sd, iproc) = sij
    do i=1,Nstates
      a_h_psi(i) += hij * psi_coef(minilist(l_sd), i)
    end do
  end do

  contrib = 0d0

  do i=1,Nstates
    de = dress_E0_denominator(i) - haa
    if(DABS(de) < 1D-5) cycle

    c_alpha(i) = a_h_psi(i) / de
    contrib = min(contrib, c_alpha(i) * a_h_psi(i))
    
    do l_sd=1,n_minilist
      hdress = c_alpha(i) * a_h_i(l_sd, iproc)
      sdress = c_alpha(i) * a_s2_i(l_sd, iproc)
      !if(c_alpha(i) * a_s2_i(l_sd, iproc) > 1d-1) then
      !  call debug_det(det_minilist(1,1,l_sd), N_int)
      !  call debug_det(alpha,N_int)
      !end if
      delta_ij_loc(i, minilist(l_sd), 1) += hdress
      delta_ij_loc(i, minilist(l_sd), 2) += sdress
    end do
  end do
end subroutine


subroutine dress_with_alpha_buffer(Nstates,Ndet,Nint,delta_ij_loc, i_gen, minilist, det_minilist, n_minilist, alpha, iproc)
  use bitmasks
  implicit none
  BEGIN_DOC
  !delta_ij_loc(:,:,1) : dressing column for H
  !delta_ij_loc(:,:,2) : dressing column for S2
  !i_gen : generator index in psi_det_generators
  !minilist : indices of determinants connected to alpha ( in psi_det )
  !n_minilist : size of minilist
  !alpha : alpha determinant
  END_DOC
  integer, intent(in)             :: Nint, Ndet, Nstates, n_minilist, iproc, i_gen
  integer(bit_kind), intent(in)   :: alpha(Nint,2), det_minilist(Nint, 2, n_minilist)
  integer,intent(in)              :: minilist(n_minilist)
  double precision, intent(inout) :: delta_ij_loc(Nstates,N_det,2)
  double precision, external :: diag_H_mat_elem_fock
  double precision :: haa, contrib, c_alpha(N_states)


  
  haa =  diag_H_mat_elem_fock(psi_det_generators(1,1,i_gen),alpha,fock_diag_tmp_(1,1,iproc),N_int)

  call dress_with_alpha_(Nstates, Ndet, Nint, delta_ij_loc, minilist, det_minilist, n_minilist, alpha, haa, contrib, c_alpha, iproc)
  slave_sum_alpha2(:,iproc) += c_alpha(:)**2
  if(contrib < sb(iproc)%mini) then
    call add_to_selection_buffer(sb(iproc), alpha, contrib)
  end if
end subroutine


BEGIN_PROVIDER [ logical, initialize_E0_denominator ]
    implicit none
    BEGIN_DOC
    ! If true, initialize pt2_E0_denominator
    END_DOC
    initialize_E0_denominator = .True.
END_PROVIDER



