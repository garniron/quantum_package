program var_pt2_ratio_run
  implicit none
  integer                        :: i,k

  
  double precision, allocatable  :: pt2(:), norm_pert(:), H_pert_diag(:)
  integer                        :: N_st, degree
  N_st = N_states
  allocate (pt2(N_st), norm_pert(N_st),H_pert_diag(N_st))
  character*(64)                 :: perturbation
  
  double precision, allocatable  :: psi_det_save(:,:,:), psi_coef_save(:,:)

  double precision :: E_fci, E_var, ratio, E_ref, selection_criterion_save
  integer :: Nmin, Nmax

  pt2 = 1.d0
  diag_algorithm = "Lapack"

  ratio = 0.d0
  Nmin=1
  do while (ratio < var_pt2_ratio)
    call H_apply_FCI(pt2, norm_pert, H_pert_diag,  N_st)
    psi_det = psi_det_sorted
    psi_coef = psi_coef_sorted
    soft_touch N_det psi_det psi_coef
    call diagonalize_CI
    ratio = (CI_energy(1) - HF_energy) / (CI_energy(1)+pt2(1) - HF_energy)
  enddo

  threshold_selectors = 1.d0
  threshold_generators = 0.999d0
  selection_criterion_save = selection_criterion
  call diagonalize_CI
  call H_apply_FCI_PT2(pt2, norm_pert, H_pert_diag,  N_st)
  E_ref = CI_energy(1) + pt2(1)
  threshold_selectors = 0.999d0
  threshold_generators = 0.99d0

  Nmax=N_det
  Nmin=1
  do while (Nmax-Nmin > 1)
    ratio = (CI_energy(1) - HF_energy) / (E_ref - HF_energy)

    if (ratio < var_pt2_ratio) then
      Nmin = N_det
      Nmax = max(Nmax,Nmin+10)
      ! Select new determinants
      call H_apply_FCI(pt2, norm_pert, H_pert_diag,  N_st)
      selection_criterion = selection_criterion_save
      SOFT_TOUCH selection_criterion selection_criterion_min selection_criterion_factor
    else
      Nmax = N_det
      N_det = Nmin + (Nmax-Nmin)/2
    endif

    psi_det = psi_det_sorted
    psi_coef = psi_coef_sorted
    soft_touch N_det psi_det psi_coef
    call diagonalize_CI
    call save_wavefunction
    print *,  'Det min, Det max: ', Nmin, Nmax
    print *,  'Ratio           : ', ratio, '  ~  ', var_pt2_ratio
    print *,  'HF_energy = ', HF_energy
    print *,  'Est FCI   = ', E_ref
    print *,  'N_det     = ', N_det
    print *,  'E         = ', CI_energy(1)
    call ezfio_set_full_ci_energy(CI_energy)
  enddo
  deallocate(pt2,norm_pert)
end
