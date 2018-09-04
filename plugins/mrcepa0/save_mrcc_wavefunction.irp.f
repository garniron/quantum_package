program save_mrcc_wf
  implicit none
  
  threshold_generators = 1.d0
  threshold_selectors = 1.d0
  PROVIDE N_int psi_det
  TOUCH threshold_generators threshold_selectors

  mrmode=5
  read_wf = .True.
  SOFT_TOUCH read_wf mrmode
  call generate_all_alpha_beta_det_products
  
  call run1
  call run2
end

subroutine run1
  implicit none

  integer :: k
  double precision :: c_alpha(N_states)
  call set_generators_bitmasks_as_holes_and_particles

  call get_cc_coef(psi_det(1,1,1), c_alpha)
  !$OMP PARALLEL DO DEFAULT(SHARED) &
  !$OMP PRIVATE(k,c_alpha) SCHEDULE(static,64)
  do k=1,N_det
!    if (maxval(abs(psi_coef(k,1:N_states))) == 0.d0) then
      if (iand(k,1023) == 0) then
        print *,  k, '/', N_det
      endif
      call get_cc_coef(psi_det(1,1,k), c_alpha)
      psi_coef(k,1:N_states) = c_alpha(1:N_states)
!    endif
  enddo
  !$OMP END PARALLEL DO
  SOFT_TOUCH psi_coef
end

subroutine run2
  implicit none

  integer :: k
  double precision :: c_alpha(N_states)

  psi_det(1:N_int,1:2,1:N_det) = psi_det_sorted(1:N_int,1:2,1:N_det)
  psi_coef(1:N_det,1:N_states) = psi_coef_sorted(1:N_det,1:N_states)
!  do k=N_det,1,-1
!    if (maxval(abs(psi_coef(k,1:N_states))) > 0.d0) then
!      exit
!    endif
!  enddo
!  N_det = k
  SOFT_TOUCH N_det psi_coef psi_det
  call save_wavefunction
end

