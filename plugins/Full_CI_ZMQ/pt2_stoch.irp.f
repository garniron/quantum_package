program pt2_stoch
  implicit none
  read_wf = .True.
  SOFT_TOUCH read_wf
  PROVIDE mo_bielec_integrals_in_map
  PROVIDE psi_energy
  call run
end

subroutine run
  implicit none
  integer                        :: i,j,k
  logical, external              :: detEq
  
  double precision               :: pt2(N_states)
  integer                        :: degree
  integer                        :: n_det_before, to_select
  double precision               :: threshold_davidson_in
  
  double precision               :: E_CI_before(N_states), relative_error, error(N_states)
  
  pt2(:) = 0.d0
  
  E_CI_before(:) = psi_energy(:) + nuclear_repulsion
  threshold_selectors = 1.d0
  threshold_generators = 1.d0
  relative_error=PT2_relative_error
  
  call ZMQ_pt2(E_CI_before, pt2, relative_error, error)
  do k=1,N_states
    print *,  'State      ', k
    print *,  'N_det    = ', N_det
    print *,  'PT2      = ', pt2
    print *,  'E        = ', E_CI_before(k)
    print *,  'E+PT2    = ', E_CI_before(k)+pt2(k), ' +/- ', error(k)
    print *,  '-----'
  enddo
  call ezfio_set_full_ci_zmq_energy_pt2(E_CI_before(1)+pt2(1))
end


