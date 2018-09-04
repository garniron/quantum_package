program shifted_bk
  implicit none
  integer                        :: i,j,k
  double precision, allocatable  :: pt2(:)
  integer                        :: degree
  integer                        :: n_det_before
  double precision               :: threshold_davidson_in
  
  allocate (pt2(N_states))

  double precision               :: hf_energy_ref
  logical                        :: has
  double precision               :: relative_error, absolute_error
  integer                        :: N_states_p
  character*(512)                :: fmt

  PROVIDE psi_bilinear_matrix_columns_loc psi_det_alpha_unique psi_det_beta_unique
  PROVIDE psi_bilinear_matrix_rows psi_det_sorted_order psi_bilinear_matrix_order
  PROVIDE psi_bilinear_matrix_transp_rows_loc psi_bilinear_matrix_transp_columns
  PROVIDE psi_bilinear_matrix_transp_order


  pt2 = -huge(1.e0)
  threshold_davidson_in = threshold_davidson
  threshold_davidson = threshold_davidson_in * 100.d0
  SOFT_TOUCH threshold_davidson

  call diagonalize_CI_dressed
  call save_wavefunction
  
  call ezfio_has_hartree_fock_energy(has)
  if (has) then
    call ezfio_get_hartree_fock_energy(hf_energy_ref)
  else
    hf_energy_ref = ref_bitmask_energy
  endif

  if (N_det > N_det_max) then
    psi_det = psi_det_sorted
    psi_coef = psi_coef_sorted
    N_det = N_det_max
    soft_touch N_det psi_det psi_coef
    call diagonalize_CI_dressed
    call save_wavefunction
    N_states_p = min(N_det,N_states)
  endif
  
  n_det_before = 0

  character*(8) :: pt2_string
  double precision :: threshold_selectors_save, threshold_generators_save
  threshold_selectors_save  = threshold_selectors
  threshold_generators_save = threshold_generators
  double precision :: error(N_states), energy(N_states)
  error = 0.d0

  threshold_selectors = 1.d0
  threshold_generators = 1d0 

  if (.True.) then 
    pt2_string = '(sh-Bk) '
    do while ( (N_det < N_det_max) )
      write(*,'(A)')  '--------------------------------------------------------------------------------'

      N_det_delta_ij = N_det

      do i=1,N_states
        energy(i) = psi_energy(i)+nuclear_repulsion
      enddo

      PROVIDE delta_ij_tmp
      call delta_ij_done()

      call diagonalize_ci_dressed
      do i=1,N_states
        pt2(i) = ci_energy_dressed(i) - energy(i)
      enddo

      N_states_p = min(N_det,N_states)

      print *, ''
      print '(A,I12)',  'Summary at N_det = ', N_det
      print '(A)',      '-----------------------------------'
      print *, ''
      print *, ''

      write(fmt,*) '(''# ============'',', N_states_p, '(1X,''=============================''))'
      write(*,fmt)
      write(fmt,*) '(12X,', N_states_p, '(6X,A7,1X,I6,10X))'
      write(*,fmt) ('State',k, k=1,N_states_p)
      write(fmt,*) '(''# ============'',', N_states_p, '(1X,''=============================''))'
      write(*,fmt)
      write(fmt,*) '(A12,', N_states_p, '(1X,F14.8,15X))'
      write(*,fmt) '# E          ', energy(1:N_states_p)
      if (N_states_p > 1) then
        write(*,fmt) '# Excit. (au)', energy(1:N_states_p)-energy(1)
        write(*,fmt) '# Excit. (eV)', (energy(1:N_states_p)-energy(1))*27.211396641308d0
      endif
      write(fmt,*) '(A12,', 2*N_states_p, '(1X,F14.8))'
      write(*,fmt) '# PT2'//pt2_string, (pt2(k), error(k), k=1,N_states_p)
      write(*,'(A)') '#'
      write(*,fmt) '# E+PT2      ', (energy(k)+pt2(k),error(k), k=1,N_states_p)
      if (N_states_p > 1) then
        write(*,fmt) '# Excit. (au)', ( (energy(k)+pt2(k)-energy(1)-pt2(1)), &
          dsqrt(error(k)*error(k)+error(1)*error(1)), k=1,N_states_p)
        write(*,fmt) '# Excit. (eV)', ( (energy(k)+pt2(k)-energy(1)-pt2(1))*27.211396641308d0, &
          dsqrt(error(k)*error(k)+error(1)*error(1))*27.211396641308d0, k=1,N_states_p)
      endif
      write(fmt,*) '(''# ============'',', N_states_p, '(1X,''=============================''))'
      write(*,fmt)
      print *,  ''

      print *,  'N_det             = ', N_det
      print *,  'N_states          = ', N_states

      do k=1, N_states_p
        print*,'State ',k
        print *,  'PT2             = ', pt2(k)
        print *,  'E               = ', energy(k)
        print *,  'E+PT2'//pt2_string//'   = ', energy(k)+pt2(k)
      enddo

      print *,  '-----'
      if(N_states.gt.1)then
        print *, 'Variational Energy difference (au | eV)'
        do i=2, N_states_p
          print*,'Delta E = ', (energy(i) - energy(1)), &
            (energy(i) - energy(1)) * 27.211396641308d0
        enddo
        print *,  '-----'
        print*, 'Variational + perturbative Energy difference (au | eV)'
        do i=2, N_states_p
          print*,'Delta E = ', (energy(i)+ pt2(i) - (energy(1) + pt2(1))), &
            (energy(i)+ pt2(i) - (energy(1) + pt2(1))) * 27.211396641308d0
        enddo
      endif
      call ezfio_set_shiftedbk_energy_pt2(energy(1)+pt2(1))
!      call dump_fci_iterations_value(N_det,energy,pt2) 

      n_det_before = N_det
      
      PROVIDE  psi_coef
      PROVIDE  psi_det
      PROVIDE  psi_det_sorted

      if (N_det >= N_det_max) then
        threshold_davidson = threshold_davidson_in
      end if
      call save_wavefunction
      call ezfio_set_shiftedbk_energy(energy(1))
      call ezfio_set_shiftedbk_energy_pt2(ci_energy_dressed(1))
    enddo
  endif




end

