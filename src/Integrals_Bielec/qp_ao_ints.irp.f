program qp_ao_ints
  use omp_lib
  implicit none
  IRP_IF MPI
    include 'mpif.h'
  IRP_ENDIF
  integer :: ierr

  BEGIN_DOC
! Increments a running calculation to compute AO integrals
  END_DOC
  integer :: i
  PROVIDE zmq_context mpi_master zmq_state zmq_context

  call switch_qp_run_to_master

  zmq_context = f77_zmq_ctx_new ()

  ! Set the state of the ZMQ
  zmq_state = 'ao_integrals'

  ! Provide everything needed
  double precision :: integral, ao_bielec_integral
  integral = ao_bielec_integral(1,1,1,1)

  do
    call wait_for_state('ao_integrals',zmq_state)
    if (zmq_state(1:7) == 'Stopped') then
      exit
    endif

    !$OMP PARALLEL DEFAULT(PRIVATE) PRIVATE(i)
    i = omp_get_thread_num()
    call ao_bielec_integrals_in_map_slave_tcp(i)
    !$OMP END PARALLEL
    IRP_IF MPI
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      if (ierr /= MPI_SUCCESS) then
        print *,  irp_here, 'error in barrier'
      endif
    IRP_ENDIF

  enddo
  IRP_IF MPI
    call MPI_finalize(i)
  IRP_ENDIF

  print *,  'Done'
end

