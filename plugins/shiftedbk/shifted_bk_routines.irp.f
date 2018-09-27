 BEGIN_PROVIDER [ double precision, fock_diag_tmp_, (2,mo_tot_num+1,Nproc) ]
&BEGIN_PROVIDER [ double precision, a_h_i, (N_det, Nproc) ]
&BEGIN_PROVIDER [ double precision, a_s2_i, (N_det, Nproc) ]
   implicit none
   integer                        :: i
   
   fock_diag_tmp_(:,:,:) = 0.d0
   a_h_i = 0d0
   a_s2_i = 0d0
END_PROVIDER


subroutine generator_start(i_gen, iproc, interesting)
  implicit none
  integer, intent(in) :: i_gen, iproc
  logical, intent(inout) :: interesting
  integer :: i
  interesting = .true.
  call build_fock_tmp(fock_diag_tmp_(1,1,iproc),psi_det_generators(1,1,i_gen),N_int)
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
  double precision :: hij, sij
  double precision, external :: diag_H_mat_elem_fock
  double precision :: haa, contrib, c_alpha(N_states)
  double precision :: de, a_h_psi(N_states)
  double precision :: hdress, sdress
  integer :: i, l_sd

  
  haa =  diag_H_mat_elem_fock(psi_det_generators(1,1,i_gen),alpha,fock_diag_tmp_(1,1,iproc),N_int)

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
      delta_ij_loc(i, minilist(l_sd), 1) += hdress
      delta_ij_loc(i, minilist(l_sd), 2) += sdress
    end do
  end do

end subroutine


BEGIN_PROVIDER [ logical, initialize_E0_denominator ]
 implicit none
 BEGIN_DOC
 ! If true, initialize pt2_E0_denominator
 END_DOC
 initialize_E0_denominator = .True.
END_PROVIDER


