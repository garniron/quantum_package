use bitmasks 

subroutine generator_start(i_gen, iproc, interesting)
  implicit none
  integer, intent(in) :: i_gen, iproc
  logical, intent(inout) :: interesting
  integer :: i
  logical, external :: deteq
  PROVIDE dij
  interesting = .true.
  do i=1,N_det_ref
    if(deteq(psi_det_generators(1,1,i_gen), psi_ref(1,1,i), N_int)) then
      interesting = .false.
      exit
    end if
  end do
end subroutine

 BEGIN_PROVIDER [ double precision, hij_cache_, (N_det,Nproc) ]
&BEGIN_PROVIDER [ double precision, sij_cache_, (N_det,Nproc) ]
&BEGIN_PROVIDER [ double precision, dIa_hla_, (N_states,N_det,Nproc) ]
&BEGIN_PROVIDER [ double precision, dIa_sla_, (N_states,N_det,Nproc) ]
&BEGIN_PROVIDER [ integer(bit_kind), sorted_mini, (N_int,2,N_det,Nproc) ]
&BEGIN_PROVIDER [ integer, excs_ , (0:2,2,2,N_det,Nproc) ]
&BEGIN_PROVIDER [ integer, idx_buf , (N_det, Nproc) ]
&BEGIN_PROVIDER [ double precision, phases_, (N_det, Nproc) ]
BEGIN_DOC
 ! temporay arrays for dress_with_alpha_buffer. Avoids reallocation.
END_DOC
END_PROVIDER


 BEGIN_PROVIDER [ integer(bit_kind), psi_ref_detsorted, (N_int,2,N_det_ref) ]
&BEGIN_PROVIDER [ integer, psi_ref_detsorted_idx, (N_det_ref) ]
  implicit none

  psi_ref_detsorted = psi_ref(:,:,:N_det_ref)
  call sort_det(psi_ref_detsorted, psi_ref_detsorted_idx, N_det_ref, n_int)

END_PROVIDER


subroutine dress_with_alpha_buffer(Nstates, Ndet,Nint,delta_ij_loc, i_gen, minilist, det_minilist, n_minilist, alpha, iproc)
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
  integer(bit_kind)               :: dettmp(Nint,2), tmp
  double precision, intent(inout) :: delta_ij_loc(Nstates,N_det,2)
  double precision :: hij, sij
  double precision, external :: diag_H_mat_elem_fock
  double precision :: c_alpha(N_states)
  double precision :: hdress, sdress
  integer :: i, l_sd, j, k, i_I, s, ni
  logical :: ok
  double precision               :: phase, phase2
  integer                        :: degree, exc(0:2,2,2)
  integer(8), save :: diamond = 0
  if(n_minilist == 1) return
  !check if not linked to reference
  do i=1,n_minilist
    if(idx_non_ref_rev(minilist(i)) == 0) then
      return
    end if
  end do
  
  sorted_mini(:,:,:n_minilist,iproc) = det_minilist(:,:,:)
  call sort_det(sorted_mini(1,1,1,iproc), idx_buf(1,iproc), n_minilist, nint)
  
  c_alpha = 0d0
  
  do i=1,n_minilist
    !call get_excitation_degree(alpha, psi_ref(1,1,i_I), degree, nint)
    !if(degree > 4) cycle
    do s=1,2
    do ni=1,nint
      dettmp(ni,s) = alpha(ni,s)-sorted_mini(ni,s,i,iproc)
    end do
    end do
    i_I=1
    j=i+1

    diamondloop : do while(i_I <= N_det_ref .and. j <= n_minilist)
      
      do s=1,2
      do ni=nint,1,-1
        if(sorted_mini(ni,s,j,iproc) - psi_ref_detsorted(ni,s,i_I) > dettmp(ni,s)) then
          i_I += 1
          cycle diamondloop
        else if(sorted_mini(ni,s,j,iproc) - psi_ref_detsorted(ni,s,i_I) < dettmp(ni,s)) then
          j += 1
          cycle diamondloop
        end if
      end do
      end do
    
      !check potential diamond found

      do s=1,2
      do ni=1,nint
        tmp = ieor(sorted_mini(ni,s,i,iproc), sorted_mini(ni,s,j,iproc))
        tmp = ieor(tmp, psi_ref_detsorted(ni,s,i_I))
        tmp = ieor(tmp, alpha(ni,s))
        if(tmp /= 0_8) then
          !print *, "fake diamond spotted"
          !i_I += 1
          j += 1
          cycle diamondloop
        end if
      end do 
      end do
      !diamond += 1
      !if(mod(diamond,100000) == 1) print *, "diam", diamond
      !diamond found
      if(det_minilist(1,1,idx_buf(j,iproc)) /= sorted_mini(1,1,j,iproc)) stop "STOOPE"
      call get_excitation(psi_ref_detsorted(1,1,i_I),det_minilist(1,1,idx_buf(j,iproc)),exc,degree,phase,Nint)
      call get_excitation(alpha,det_minilist(1,1,idx_buf(i,iproc)),exc,degree,phase2,Nint)

      do s=1,Nstates
        c_alpha(s) += psi_ref_coef(psi_ref_detsorted_idx(i_I), s) * dij(psi_ref_detsorted_idx(i_I), idx_non_ref_rev(minilist(idx_buf(i,iproc))), s) &
            * dij(psi_ref_detsorted_idx(i_I), idx_non_ref_rev(minilist(idx_buf(j,iproc))), s) * phase * phase2
      end do
      !i_I += 1
      j += 1
    end do diamondloop
  end do

  if(maxval(c_alpha) == 0d0 .and. minval(c_alpha) == 0d0) return
  
  do i=1,n_minilist
    call i_h_j_s2(alpha,det_minilist(1,1,i),N_int,hij, sij)
    do s=1,Nstates
      hdress = c_alpha(s) * hij
      sdress = c_alpha(s) * sij
      delta_ij_loc(s, minilist(i), 1) += hdress
      delta_ij_loc(s, minilist(i), 2) += sdress
    end do
  end do
end subroutine



subroutine dress_with_alpha_buffer_neu(Nstates,Ndet,Nint,delta_ij_loc, i_gen, minilist, det_minilist, n_minilist, alpha, iproc)
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
  integer(bit_kind)               :: dettmp(Nint,2)
  double precision, intent(inout) :: delta_ij_loc(Nstates,N_det,2)
  double precision :: hij, sij
  double precision, external :: diag_H_mat_elem_fock
  double precision :: c_alpha(N_states)
  double precision :: hdress, sdress
  integer :: i, l_sd, j, k, i_I, s, ni
  logical :: ok
  double precision               :: phase, phase2
  integer                        :: degree, exc(0:2,2,2)
  integer(8), save :: diamond = 0
  if(n_minilist == 1) return
  !check if not linked to reference
  do i=1,n_minilist
    if(idx_non_ref_rev(minilist(i)) == 0) then
      return
    end if
  end do
  
  c_alpha = 0d0
  
  do i_I=1,N_det_ref
    call get_excitation_degree(alpha, psi_ref(1,1,i_I), degree, nint)
    if(degree > 4) cycle

    do i=1,n_minilist
    diamondloop : do j=i+1,n_minilist
      do s=1,2
      do ni=1,nint
        dettmp(ni,s) = ieor(det_minilist(ni,s,i), det_minilist(ni,s,j))
        dettmp(ni,s) = ieor(dettmp(ni,s), psi_ref(ni,s,i_I))
        dettmp(ni,s) = ieor(dettmp(ni,s), alpha(ni,s))
        if(dettmp(ni,s) /= 0_8) cycle diamondloop
      end do 
      end do
      !diamond found
      diamond += 1
      if(mod(diamond,10000) == 1) print *, "diam", diamond

      call get_excitation(psi_ref(1,1,i_I),det_minilist(1,1,j),exc,degree,phase,Nint)
      call get_excitation(alpha,det_minilist(1,1,i),exc,degree,phase2,Nint)

      do s=1,Nstates
        c_alpha(s) += psi_ref_coef(i_I, s) * dij(i_I, idx_non_ref_rev(minilist(i)), s) &
           * dij(i_I, idx_non_ref_rev(minilist(j)), s) * phase * phase2
      end do
    end do diamondloop
    end do
  end do
  
  if(maxval(c_alpha) == 0d0 .and. minval(c_alpha) == 0d0) return
  
  do i=1,n_minilist
    call i_h_j_s2(alpha,det_minilist(1,1,i),N_int,hij, sij)
    do s=1,Nstates
      hdress = c_alpha(s) * hij
      sdress = c_alpha(s) * sij
      delta_ij_loc(s, minilist(i), 1) += hdress
      delta_ij_loc(s, minilist(i), 2) += sdress
    end do
  end do
end subroutine



subroutine dress_with_alpha_buffer_old(Nstates,Ndet,Nint,delta_ij_loc, i_gen, minilist, det_minilist, n_minilist, alpha, iproc)
 use bitmasks
 implicit none
  BEGIN_DOC
  !delta_ij_loc(:,:,1) : dressing column for H
  !delta_ij_loc(:,:,2) : dressing column for S2
  !minilist : indices of determinants connected to alpha ( in psi_det_sorted )
  !n_minilist : size of minilist
  !alpha : alpha determinant
  END_DOC
  integer(bit_kind), intent(in)   :: alpha(Nint,2), det_minilist(Nint, 2, n_minilist)
  integer,intent(in)              :: minilist(n_minilist), n_minilist, iproc, i_gen, Nstates, Ndet, Nint
  double precision, intent(inout) :: delta_ij_loc(Nstates,Ndet,2)


  integer                        :: i,j,k,l,m
  integer                        :: degree1, degree2, degree

  double precision               :: hIk, hla, hIl, sla, dIk(Nstates), dka(Nstates), dIa(Nstates), hka
  double precision               :: phase, phase2
  integer                        :: exc(0:2,2,2)
  integer                        :: h1,h2,p1,p2,s1,s2
  integer(bit_kind)              :: tmp_det(Nint,2), ctrl
  integer                        :: i_state, k_sd, l_sd, m_sd, ll_sd, i_I
  double precision :: Delta_E_inv(Nstates)
  double precision :: sdress, hdress
  logical :: ok, ok2
  integer :: canbediamond

  PROVIDE mo_class dij N_int N_states elec_num n_act_orb 
  
  if(n_minilist == 1) return
  
  do i=1,n_minilist
    if(idx_non_ref_rev(minilist(i)) == 0) return
  end do

  if (perturbative_triples) then
    PROVIDE one_anhil fock_virt_total fock_core_inactive_total one_creat
  endif

  canbediamond = 0
  do l_sd=1,n_minilist
    call get_excitation(det_minilist(1,1,l_sd),alpha,exc,degree1,phase,Nint)
    call decode_exc(exc,degree1,h1,p1,h2,p2,s1,s2)
   
    ok = (mo_class(h1)(1:1) == 'A' .or. mo_class(h1)(1:1) == 'I') .and. &
         (mo_class(p1)(1:1) == 'A' .or. mo_class(p1)(1:1) == 'V') 
    if(ok .and. degree1 == 2) then
            ok = (mo_class(h2)(1:1) == 'A' .or. mo_class(h2)(1:1) == 'I') .and. &
                 (mo_class(p2)(1:1) == 'A' .or. mo_class(p2)(1:1) == 'V') 
    end if
    
    if(ok) then
      canbediamond += 1
      excs_(:,:,:,l_sd,iproc) = exc(:,:,:)
      phases_(l_sd, iproc) = phase
    else
      phases_(l_sd, iproc) = 0d0
    end if
    call i_h_j_s2(alpha,det_minilist(1,1,l_sd),Nint,hij_cache_(l_sd,iproc), sij_cache_(l_sd,iproc))
  enddo
  if(canbediamond <= 1) return

  do i_I=1,N_det_ref
    call get_excitation_degree(alpha,psi_ref(1,1,i_I),degree1,Nint)
    if (degree1 > 4) then
      cycle
    endif
    
    do i_state=1,Nstates
      dIa(i_state) = 0.d0
    enddo

    do k_sd=1,n_minilist
      if(phases_(k_sd,iproc) == 0d0) cycle
      call get_excitation_degree(psi_ref(1,1,i_I),det_minilist(1,1,k_sd),degree,Nint)
      if (degree > 2) then
        cycle
      endif
      
      phase = phases_(k_sd, iproc)
      exc(:,:,:) = excs_(:,:,:,k_sd,iproc)
      degree2 = exc(0,1,1) + exc(0,1,2)
      call apply_excitation(psi_ref(1,1,i_I), exc, tmp_det, ok, Nint)
      if((.not. ok) .and. (.not. perturbative_triples)) cycle

      do i_state=1,Nstates
        dka(i_state) = 0.d0
      enddo
      
      ok2 = .false.
      !do i_state=1,Nstates
      !  !if(dka(i_state) == 0) cycle
      !  dIk(i_state) = dij(i_I, idx_non_ref_rev(minilist(k_sd)), i_state)
      !  if(dIk(i_state) /= 0d0) then
      !    ok2 = .true.
      !  endif
      !enddo
      !if(.not. ok2) cycle

      if (ok) then
        phase2 = 0d0
        do l_sd=k_sd+1,n_minilist
          if(phases_(l_sd, iproc) == 0d0) cycle
          call get_excitation_degree(tmp_det,det_minilist(1,1,l_sd),degree,Nint)
          if (degree == 0) then
            do i_state=1,Nstates
              dIk(i_state) = dij(i_I, idx_non_ref_rev(minilist(k_sd)), i_state)
              if(dIk(i_state) /= 0d0) then
                if(phase2 == 0d0) call get_excitation(psi_ref(1,1,i_I),det_minilist(1,1,l_sd),exc,degree,phase2,Nint)
                dka(i_state) = dij(i_I, idx_non_ref_rev(minilist(l_sd)), i_state) * phase * phase2
              end if
            end do

            exit

          endif
        enddo
      else if (perturbative_triples) then
        hka = hij_cache_(k_sd,iproc)
        if (dabs(hka) > 1.d-12) then
          call get_delta_e_dyall_general_mp(psi_ref(1,1,i_I),alpha,Delta_E_inv)

          do i_state=1,Nstates
            ASSERT (Delta_E_inv(i_state) < 0.d0)
            dka(i_state) = hka / Delta_E_inv(i_state)
          enddo
        endif
      endif
      
           
      if (perturbative_triples.and. (degree2 == 1) ) then
          if(sum(popcnt(tmp_det(:,1))) /= elec_alpha_num) stop "STOP 1"
  if(sum(popcnt(tmp_det(:,2))) /= elec_beta_num) stop "STOP 2"
  if(sum(popcnt(tmp_det(:,1))) /= elec_alpha_num) stop "STOP 3"
  if(sum(popcnt(tmp_det(:,2))) /= elec_beta_num) stop "STOP 4"
  

          call i_h_j(psi_ref(1,1,i_I),tmp_det,Nint,hka)
          hka = hij_cache_(k_sd,iproc) - hka
          if (dabs(hka) > 1.d-12) then
            call get_delta_e_dyall_general_mp(psi_ref(1,1,i_I),alpha,Delta_E_inv)
            do i_state=1,Nstates
              ASSERT (Delta_E_inv(i_state) < 0.d0)
              dka(i_state) = hka / Delta_E_inv(i_state)
            enddo
          endif
      endif
      do i_state=1,Nstates
        dIa(i_state) = dIa(i_state) + dIk(i_state) * dka(i_state)
      enddo
    enddo
    
    ok2 = .false.
    do i_state=1,Nstates
      if(dIa(i_state) /= 0d0) ok2 = .true.
    enddo
    if(.not. ok2) cycle

    do l_sd=1,n_minilist
      k_sd = minilist(l_sd)
      hla = hij_cache_(l_sd,iproc)
      sla = sij_cache_(l_sd,iproc)
      do i_state=1,Nstates
        hdress =  dIa(i_state) * hla * psi_ref_coef(i_I,i_state)
        sdress =  dIa(i_state) * sla * psi_ref_coef(i_I,i_state)
        !!!$OMP ATOMIC
        delta_ij_loc(i_state,k_sd,1) += hdress
        !!!$OMP ATOMIC
        delta_ij_loc(i_state,k_sd,2) += sdress
      enddo
    enddo
  enddo
end subroutine





!! TESTS MINILIST
subroutine test_minilist(minilist, n_minilist, alpha)
  use bitmasks
  implicit none
  integer, intent(in)            :: n_minilist
  integer(bit_kind),intent(in)  :: alpha(N_int, 2)
  integer, intent(in) :: minilist(n_minilist)
  integer :: a, i, deg
  integer :: refc(N_det), testc(N_det)
  
  refc = 0
  testc = 0
  do i=1,N_det
    call get_excitation_degree(psi_det(1,1,i), alpha, deg, N_int) 
    if(deg <= 2) refc(i) = refc(i) + 1
  end do
  do i=1,n_minilist
    call get_excitation_degree(psi_det(1,1,minilist(i)), alpha, deg, N_int) 
    if(deg <= 2) then
      testc(minilist(i)) += 1
    else
      stop "NON LINKED IN MINILIST"
    end if
  end do
  
  do i=1,N_det
    if(refc(i) /= testc(i)) then
      print *, "MINILIST FAIL ", sum(refc), sum(testc), n_minilist
      exit
    end if
  end do
end subroutine


