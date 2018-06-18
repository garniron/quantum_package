! DO NOT MODIFY BY HAND
! Created by $QP_ROOT/scripts/ezfio_interface/ei_handler.py
! from file /home/garniron/quantum_package/src/dress_zmq_det/EZFIO.cfg


BEGIN_PROVIDER [ double precision, thresh_dressed_ci  ]
  implicit none
  BEGIN_DOC
! Threshold on the convergence of the dressed CI energy
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then
    
    call ezfio_has_dress_zmq_det_thresh_dressed_ci(has)
    if (has) then
      call ezfio_get_dress_zmq_det_thresh_dressed_ci(thresh_dressed_ci)
    else
      print *, 'dress_zmq_det/thresh_dressed_ci not found in EZFIO file'
      stop 1
    endif
  endif
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST( thresh_dressed_ci, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read thresh_dressed_ci with MPI'
    endif
  IRP_ENDIF

  call write_time(6)
  if (mpi_master) then
    write(6, *) 'Read  thresh_dressed_ci'
  endif

END_PROVIDER

BEGIN_PROVIDER [ integer, n_it_max_dressed_ci  ]
  implicit none
  BEGIN_DOC
! Maximum number of dressed CI iterations
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then
    
    call ezfio_has_dress_zmq_det_n_it_max_dressed_ci(has)
    if (has) then
      call ezfio_get_dress_zmq_det_n_it_max_dressed_ci(n_it_max_dressed_ci)
    else
      print *, 'dress_zmq_det/n_it_max_dressed_ci not found in EZFIO file'
      stop 1
    endif
  endif
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST( n_it_max_dressed_ci, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read n_it_max_dressed_ci with MPI'
    endif
  IRP_ENDIF

  call write_time(6)
  if (mpi_master) then
    write(6, *) 'Read  n_it_max_dressed_ci'
  endif

END_PROVIDER
