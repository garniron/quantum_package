subroutine NOFT_JKfunc(nMO,FL,ET,EV,n)
  implicit none
  BEGIN_DOC
! JK-only functionals for NOFT
  END_DOC

  PROVIDE mo_bielec_integrals_in_map

! Input variables

  integer,intent(in)            :: nMO,FL
  double precision,intent(in)   :: ET,EV
  double precision,intent(in)   :: n(nMO)

! Local variables

  integer                       :: i,j
  double precision              :: EJ_SD
  double precision              :: EK_SD,EK_MBB,EK_POWER,EK_BBC1,EK_BBC2,EK_CA,EK_CGA,EK_GU,EK_ML,EK_MLSIC
  double precision              :: E_SD,E_MBB,E_POWER,E_BBC1,E_BBC2,E_CA,E_CGA,E_GU,E_ML,E_MLSIC
  double precision              :: alpha,a0,a1,b1
  double precision              :: f_ij,get_mo_bielec_integral
  double precision,allocatable  :: Jint(:,:),Kint(:,:)

! Memory allocation

  allocate(Jint(nMO,nMO),Kint(nMO,nMO))

! Coulomb, exchange and time-inversion integrals

  do i=1,nMO
    do j=1,nMO

      Jint(i,j) = get_mo_bielec_integral(i,j,i,j,mo_integrals_map)
      Kint(i,j) = get_mo_bielec_integral(i,j,j,i,mo_integrals_map)

    enddo
  enddo

! Compute SD Coulomb energy

  EJ_SD = 2d0*dot_product(n,matmul(Jint,n))

! Compute SD exchange energy

  EK_SD = dot_product(n,matmul(Kint,n))

! Compute MBB exchange energy

  EK_MBB = 0d0

  do i=1,nMO
    do j=1,nMO

      EK_MBB = EK_MBB + sqrt(n(i)*n(j))*Kint(i,j)

    enddo
  enddo

! Compute BBC1 exchange energy

  EK_BBC1 = 0d0

  do i=1,nMO
    do j=1,nMO

      if(i /= j .and. i > FL .and. j > FL) then
        f_ij = - sqrt(n(i)*n(j))
      else
        f_ij = sqrt(n(i)*n(j))
      endif

      EK_BBC1 = EK_BBC1 + f_ij*Kint(i,j)

    enddo
  enddo

! Compute BBC2 exchange energy

  EK_BBC2 = 0d0

  do i=1,nMO
    do j=1,i-1

      if(i > FL .and. j > FL) then
        f_ij = - sqrt(n(i)*n(j))
      elseif(i <= FL .and. j <= FL) then
        f_ij = n(i)*n(j)
      else
        f_ij = sqrt(n(i)*n(j))
      endif

      EK_BBC2 = EK_BBC2 + f_ij*Kint(i,j)

    enddo

    EK_BBC2 = EK_BBC2 + n(i)*Kint(i,i)
     
    do j=i+1,nMO

      if(i > FL .and. j > FL) then
        f_ij = - sqrt(n(i)*n(j))
      elseif(i <= FL .and. j <= FL) then
        f_ij = n(i)*n(j)
      else
        f_ij = sqrt(n(i)*n(j))
      endif

      EK_BBC2 = EK_BBC2 + f_ij*Kint(i,j)

    enddo
  enddo
  
! Compute CA exchange energy

  EK_CA = 0d0

  do i=1,nMO
    do j=1,nMO
      EK_CA = EK_CA + (sqrt(n(i)*n(j)*(1d0 - n(i))*(1d0 - n(j))) + n(i)*n(j))*Kint(i,j)
    enddo
  enddo

! Compute CGA exchange energy

  EK_CGA = 0d0

  do i=1,nMO
    do j=1,nMO
      EK_CGA = EK_CGA + 0.5d0*(sqrt(n(i)*n(j)*(2d0 - n(i))*(2d0 - n(j))) + n(i)*n(j))*Kint(i,j)
    enddo
  enddo

! Compute ML exchange energy

  EK_ML = 0d0

  a0 = 126.3101d0
  a1 = 2213.33d0
  b1 = 2338.64d0

  do i=1,nMO
    do j=1,nMO
      EK_ML = EK_ML + n(i)*n(j)*(a0 + a1*n(i)*n(j))/(1d0 + b1*n(i)*n(j))*Kint(i,j)
    enddo
  enddo

! Compute MLSIC exchange energy

  EK_MLSIC = 0d0

  a0 = 1298.78d0
  a1 = 35114.4d0
  b1 = 36412.2d0

  do i=1,nMO

    do j=1,i-1
      EK_MLSIC = EK_MLSIC + n(i)*n(j)*(a0 + a1*n(i)*n(j))/(1d0 + b1*n(i)*n(j))*Kint(i,j)
    enddo

    EK_MLSIC = EK_MLSIC + n(i)*n(i)*Kint(i,i)

    do j=i+1,nMO
      EK_MLSIC = EK_MLSIC + n(i)*n(j)*(a0 + a1*n(i)*n(j))/(1d0 + b1*n(i)*n(j))*Kint(i,j)
    enddo

  enddo

! Compute GU exchange energy

  EK_GU = 0d0

  do i=1,nMO

    do j=1,i-1
      EK_GU = EK_GU + sqrt(n(i)*n(j))*Kint(i,j)
    enddo

    EK_GU = EK_GU + n(i)*n(i)*Kint(i,j)

    do j=i+1,nMO
      EK_GU = EK_GU + sqrt(n(i)*n(j))*Kint(i,j)
    enddo

  enddo

! Compute POWER exchange energy

  EK_POWER = 0d0
  alpha = 1d0/3d0

  do i=1,nMO
    do j=1,nMO
      EK_POWER = EK_POWER + (n(i)*n(j))**alpha*Kint(i,j)
    enddo
  enddo

! Compute total energies

  E_SD    = ET + EV + EJ_SD - EK_SD
  E_MBB   = ET + EV + EJ_SD - EK_MBB
  E_BBC1  = ET + EV + EJ_SD - EK_BBC1
  E_BBC2  = ET + EV + EJ_SD - EK_BBC2
  E_CA    = ET + EV + EJ_SD - EK_CA
  E_CGA   = ET + EV + EJ_SD - EK_CGA
  E_ML    = ET + EV + EJ_SD - EK_ML
  E_MLSIC = ET + EV + EJ_SD - EK_MLSIC
  E_GU    = ET + EV + EJ_SD - EK_GU
  E_POWER = ET + EV + EJ_SD - EK_POWER

! Dump energies

  print*, '*******************************'
  print*, '*** JK  NOFT    functionals ***'
  print*, '*******************************'
  print*, ''
  print*, '*** Coulomb        energies ***'
  print*, 'Coulomb SD               energy = ',EJ_SD
  print*, ''
  print*, '*** Exchange       energies ***'
  print*, 'Exchange SD              energy = ',-EK_SD
  print*, 'Exchange MBB             energy = ',-EK_MBB
  print*, 'Exchange BBC1            energy = ',-EK_BBC1
  print*, 'Exchange BBC2            energy = ',-EK_BBC2
  print*, 'Exchange CA              energy = ',-EK_CA
  print*, 'Exchange CGA             energy = ',-EK_CGA
  print*, 'Exchange ML              energy = ',-EK_ML
  print*, 'Exchange MLSIC           energy = ',-EK_MLSIC
  print*, 'Exchange GU              energy = ',-EK_GU
  print*, 'Exchange POWER           energy = ',-EK_POWER
  print*, ''
  print*, ''
  print*, '*** Two-electron   energies ***'
  print*, 'J+K      SD              energy = ',EJ_SD - EK_SD
  print*, 'J+K      MBB             energy = ',EJ_SD - EK_MBB
  print*, 'J+K      BBC1            energy = ',EJ_SD - EK_BBC1
  print*, 'J+K      BBC2            energy = ',EJ_SD - EK_BBC2
  print*, 'J+K      CA              energy = ',EJ_SD - EK_CA
  print*, 'J+K      CGA             energy = ',EJ_SD - EK_CGA
  print*, 'J+K      ML              energy = ',EJ_SD - EK_ML
  print*, 'J+K      MLSIC           energy = ',EJ_SD - EK_MLSIC
  print*, 'J+K      GU              energy = ',EJ_SD - EK_GU
  print*, 'J+K      POWER           energy = ',EJ_SD - EK_POWER
  print*, ''
  print*, '*** Total          energies ***'
  print*, 'Total    SD              energy = ',E_SD
  print*, 'Total    MBB             energy = ',E_MBB
  print*, 'Total    BBC1            energy = ',E_BBC1
  print*, 'Total    BBC2            energy = ',E_BBC2
  print*, 'Total    CA              energy = ',E_CA
  print*, 'Total    CGA             energy = ',E_CGA
  print*, 'Total    ML              energy = ',E_ML
  print*, 'Total    MLSIC           energy = ',E_MLSIC
  print*, 'Total    GU              energy = ',E_GU
  print*, 'Total    POWER           energy = ',E_POWER
  print*, ''

end subroutine NOFT_JKfunc

