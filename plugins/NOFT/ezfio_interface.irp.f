! DO NOT MODIFY BY HAND
! Created by $QP_ROOT/scripts/ezfio_interface/ei_handler.py
! from file /home/loos/quantum_package/src/NOFT/EZFIO.cfg


BEGIN_PROVIDER [ logical, do_jk_functionals  ]
  implicit none
  BEGIN_DOC
! Compute energies for JK-only functionals
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then
    
    call ezfio_has_noft_do_jk_functionals(has)
    if (has) then
      call ezfio_get_noft_do_jk_functionals(do_jk_functionals)
    else
      print *, 'noft/do_jk_functionals not found in EZFIO file'
      stop 1
    endif
  endif
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST( do_jk_functionals, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read do_jk_functionals with MPI'
    endif
  IRP_ENDIF

  call write_time(6)
  if (mpi_master) then
    write(6, *) 'Read  do_jk_functionals'
  endif

END_PROVIDER

BEGIN_PROVIDER [ logical, do_jkl_functionals  ]
  implicit none
  BEGIN_DOC
! Compute energies for JKL-only functionals (PNOFs)
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then
    
    call ezfio_has_noft_do_jkl_functionals(has)
    if (has) then
      call ezfio_get_noft_do_jkl_functionals(do_jkl_functionals)
    else
      print *, 'noft/do_jkl_functionals not found in EZFIO file'
      stop 1
    endif
  endif
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST( do_jkl_functionals, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read do_jkl_functionals with MPI'
    endif
  IRP_ENDIF

  call write_time(6)
  if (mpi_master) then
    write(6, *) 'Read  do_jkl_functionals'
  endif

END_PROVIDER

BEGIN_PROVIDER [ logical, do_pt2_noft  ]
  implicit none
  BEGIN_DOC
! Compute PT2 correction for NOFT
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then
    
    call ezfio_has_noft_do_pt2_noft(has)
    if (has) then
      call ezfio_get_noft_do_pt2_noft(do_pt2_noft)
    else
      print *, 'noft/do_pt2_noft not found in EZFIO file'
      stop 1
    endif
  endif
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST( do_pt2_noft, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read do_pt2_noft with MPI'
    endif
  IRP_ENDIF

  call write_time(6)
  if (mpi_master) then
    write(6, *) 'Read  do_pt2_noft'
  endif

END_PROVIDER
