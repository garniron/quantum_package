program NOFT
  implicit none
  BEGIN_DOC
! Natural orbital functional theory module
  END_DOC

  PROVIDE mo_bielec_integrals_in_map 

  integer                       :: i,j
  integer                       :: nMO,FL
  double precision              :: ET,EV
  double precision              :: integral,get_mo_bielec_integral
  double precision,allocatable  :: n(:)

  print*, ''
  print*, '*******************************'
  print*, '*** NOFT        functionals ***'
  print*, '*******************************'
  print*, ''
  print*, 'SD     = single determinant'
  print*, 'MBB    = Muller, Buijse and Baerends'
  print*, 'POWER  = Cioslowski and Pernal'
  print*, 'BCC2   = Gritsenko and coworkers'
  print*, 'CA     = Csanyi and Arias'
  print*, 'CGA    = Csanyi, Goedecker and Arias'
  print*, 'GU     = Goedecker and Umrigar'
  print*, 'ML     = Marques and Lathiotakis'
  print*, 'MLSIC  = ML with self-interaction correction'
  print*, 'PNOF2  = Piris natural orbital functional 2 (bug)'
  print*, 'PNOF3  = Piris natural orbital functional 3'
  print*, 'PNOF4  = Piris natural orbital functional 4'
  print*, 'PNOF5  = Piris natural orbital functional 5 (NYI)'
  print*, 'PNOF6x = Piris natural orbital functional 6 (x = d, u, h)'
  print*, 'PNOF7  = Piris natural orbital functional 7 (NYI)'
  print*, ''
  print*, '*******************************'
  print*, ''
  print*, '*******************************'
  print*, '*** NOFT           energies ***'
  print*, '*******************************'
  print*, ''

! Occupation numbers

  nMO = mo_tot_num
  FL = elec_num/2
  allocate(n(nMO))
  n(1:nMO) = 0.5d0*mo_occ(1:mo_tot_num)

! Compute core energies

  call NOFT_core(nMO,ET,EV,n)

! JK-only functionals

  if(do_JK_functionals) call NOFT_JKfunc(nMO,FL,ET,EV,n)

! JKL-only functionals

  if(do_JKL_functionals) call NOFT_JKLfunc(nMO,FL,ET,EV,n)

! PT2-NOFT correction

  if(do_PT2_NOFT) call NOFT_JKLfunc(nMO,FL,n)

! End

  print*, '*******************************'
  print*, ''

end

