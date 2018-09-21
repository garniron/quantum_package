subroutine NOFT_core(nMO,ET,EV,n)
  implicit none
  BEGIN_DOC
! Core energy for NOFT
  END_DOC

! Input variables

  integer,intent(in)            :: nMO
  double precision,intent(in)   :: n(nMO)

! Local variables

  integer                       :: i

! Output variables

  double precision,intent(out)  :: ET,EV

! Compute kinetic energy

  ET = 0d0

  do i=1,nMO
    ET = ET + n(i)*mo_kinetic_integral(i,i)
  enddo

  ET = 2d0*ET

! Compute nuclear attraction energy

  EV = 0d0

  do i=1,nMO
    EV = EV + n(i)*mo_nucl_elec_integral(i,i)
  enddo

  EV = 2d0*EV

! Dump energies

  print*, '*******************************'
  print*, '*** Core           energies ***'
  print*, '*******************************'
  print*, ''
  print*, 'Kinetic                  energy = ',ET
  print*, 'Nuclear attraction       energy = ',EV
  print*, ''

end subroutine NOFT_core

