subroutine NOFT_PT2(nMO,FL,n)
! Compute the PT2 correction from NOFT
  END_DOC

  PROVIDE mo_bielec_integrals_in_map

! Input variables

  integer,intent(in)            :: nMO,FL
  double precision,intent(in)   :: n(nMO)

! Local variables

  integer                       :: i,j,a,b
  double precision              :: EPT1,EPT2
  double precision              :: get_mo_bielec_integral

! memory allocation
  
! Useful quantities

  EPT2 = 0d0

! do i=1,FL
!   do j=1,FL
!     do a=FL+1,nMO
!       do b=FL+1,nMO

!       enddo
!     enddo
!   enddo
! enddo

! Dump energies

  print*, '*******************************'
  print*, '*** PT2 NOFT    corrections ***'
  print*, '*******************************'
  print*, ''
  print*, 'Total    PT1             energy = ',E_PT1
  print*, 'Total    PT2             energy = ',E_PT2
  print*, 'Total    PT1+PT2         energy = ',E_PT1 + E_PT2
  print*, ''

end subroutine NOFT_PT2

