subroutine NOFT_JKLfunc(nMO,FL,ET,EV,n)
! JKL-only functionals for NOFT
  END_DOC

  PROVIDE mo_bielec_integrals_in_map

! Input variables

  integer,intent(in)            :: nMO,FL
  double precision,intent(in)   :: ET,EV
  double precision,intent(in)   :: n(nMO)

! Local variables

  integer                       :: i,j
  double precision              :: EJ_SD,EK_SD
  double precision              :: EJ_PNOF2,EJ_PNOF3,EJ_PNOF4,EJ_PNOF5,EJ_PNOF6d,EJ_PNOF6u,EJ_PNOF6h,EJ_PNOF7
  double precision              :: EK_PNOF2,EK_PNOF3,EK_PNOF4,EK_PNOF5,EK_PNOF6d,EK_PNOF6u,EK_PNOF6h,EK_PNOF7
  double precision              :: EL_PNOF2,EL_PNOF3,EL_PNOF4,EL_PNOF5,EL_PNOF6d,EL_PNOF6u,EL_PNOF6h,EL_PNOF7
  double precision              :: E_PNOF2,E_PNOF3,E_PNOF4,E_PNOF5,E_PNOF6d,E_PNOF6u,E_PNOF6h,E_PNOF7
  double precision              :: get_mo_bielec_integral

  double precision              :: SF,Sd,Su,Sh,Delta_ij,T_ij,Pi_ij

  double precision,allocatable  :: h(:),kappa(:),gam(:),Jint(:,:),Kint(:,:),Lint(:,:)

! memory allocation
  
  allocate(h(nMO),kappa(nMO),gam(nMO),Jint(nMO,nMO),Kint(nMO,nMO),Lint(nMO,nMO))

! Useful quantities

  h(1:nMO) = 1d0 - n(1:nMO)

  SF = 0d0
  do i=1,FL
    SF = SF + h(i)
  enddo

! Useful quantities for PNOF6

  do i=1,FL
    kappa(i) = h(i)*exp(-SF)
  enddo
  do i=FL+1,nMO
    kappa(i) = n(i)*exp(-SF)
  enddo

  gam(:) = 0d0
  do i=1,nMO
    do j=i,FL
      gam(i) = gam(i) + kappa(j)
    enddo
    gam(i) = n(i)*h(i) + kappa(i)*kappa(i) - kappa(i)*gam(i)
  enddo

  Sd = 0d0
  do i=1,FL
    Sd = Sd + gam(i)
  enddo

  Su = 0d0
  do i=FL+1,nMO
    Su = Su + gam(i)
  enddo

  Sh = 0.5d0*(Sd + Su)

! Coulomb, exchange and time-inversion integrals

  do i=1,nMO
    do j=1,nMO

      Jint(i,j) = get_mo_bielec_integral(i,j,i,j,mo_integrals_map)
      Kint(i,j) = get_mo_bielec_integral(i,j,j,i,mo_integrals_map)
      Lint(i,j) = get_mo_bielec_integral(i,i,j,j,mo_integrals_map)

    enddo
  enddo

!****************************************
!*** Coulomb and exchange parts of SD ***
!****************************************

  EJ_SD = +2d0*dot_product(n,matmul(Jint,n))
  EK_SD = -1d0*dot_product(n,matmul(Kint,n))

! *************
! *** PNOF2 ***
! *************

  EJ_PNOF2 = 0d0
  EK_PNOF2 = 0d0
  EL_PNOF2 = 0d0

  do i=1,nMO
    do j=1,nMO

      if(i == j) then

        Delta_ij = n(i)*n(j)
        T_ij     = 0d0
        Pi_ij    = sqrt(n(i)*n(j))

      elseif(i <= FL .and. j <= FL) then

        Delta_ij = h(i)*h(j)
        T_ij     = n(i)*n(j) - Delta_ij
        Pi_ij    = sqrt(n(i)*n(j)) + sqrt(h(i)*h(j)) + T_ij

      elseif(i <= FL .and. j > FL) then

        Delta_ij = n(j)*h(i)*(1d0 - SF)/SF
        T_ij     = n(i)*n(j) - Delta_ij
        Pi_ij    = sqrt(n(i)*n(j)) - sqrt(n(j)*h(i)) + T_ij

      elseif(i > FL .and. j <= FL) then

        Delta_ij = n(i)*h(j)*(1d0 - SF)/SF
        T_ij     = n(i)*n(j) - Delta_ij
        Pi_ij    = sqrt(n(i)*n(j)) - sqrt(n(i)*h(j)) + T_ij

      elseif(i > FL .and. j > FL) then

        Delta_ij = n(i)*n(j)
        T_ij     = n(i)*n(j) - Delta_ij
        Pi_ij    = T_ij

      else

        Delta_ij = 0d0
        T_ij     = 0d0
        Pi_ij    = 0d0

      endif

      EJ_PNOF2 = EJ_PNOF2 - 2d0*Delta_ij*Jint(i,j)
      EK_PNOF2 = EK_PNOF2 + Delta_ij*Kint(i,j)
      EL_PNOF2 = EL_PNOF2 + Pi_ij*Lint(i,j)

    enddo   
  enddo   

! *************
! *** PNOF3 ***
! *************

  EJ_PNOF3 = 0d0
  EK_PNOF3 = 0d0
  EL_PNOF3 = 0d0

  do i=1,nMO
    do j=1,nMO

      if(i == j) then

        Delta_ij = n(i)*n(j)
        Pi_ij    = sqrt(n(i)*n(j))

      elseif(i <= FL .and. j <= FL) then

        Delta_ij = h(i)*h(j)
        Pi_ij    = n(i)*n(j) - sqrt(n(i)*n(j))

      elseif(i <= FL .and. j > FL) then

        Delta_ij = n(j)*h(i)*(1d0 - SF)/SF
        Pi_ij    = n(i)*n(j) - sqrt(n(i)*n(j)) - sqrt(n(j)*h(i)) 

      elseif(i > FL .and. j <= FL) then

        Delta_ij = n(i)*h(j)*(1d0 - SF)/SF
        Pi_ij    = n(i)*n(j) - sqrt(n(i)*n(j)) - sqrt(n(i)*h(j)) 

      elseif(i > FL .and. j > FL) then

        Delta_ij = n(i)*n(j)
        Pi_ij    = n(i)*n(j) + sqrt(n(i)*n(j))

      else

        Delta_ij = 0d0
        Pi_ij    = 0d0

      endif

      EJ_PNOF3 = EJ_PNOF3 - Delta_ij*Jint(i,j)
      EK_PNOF3 = EK_PNOF3 
      EL_PNOF3 = EL_PNOF3 + Pi_ij*Lint(i,j)

    enddo   
  enddo   

! *************
! *** PNOF4 ***
! *************

  EJ_PNOF4 = 0d0
  EK_PNOF4 = 0d0
  EL_PNOF4 = 0d0

  do i=1,nMO
    do j=1,nMO

      if(i == j) then

        Delta_ij = n(i)*n(j)
        Pi_ij    = sqrt(n(i)*n(j))

      elseif(i <= FL .and. j <= FL) then

        Delta_ij = h(i)*h(j)
        Pi_ij    = - sqrt(h(i)*h(j))

      elseif(i <= FL .and. j > FL) then

        Delta_ij = n(j)*h(i)*(1d0 - SF)/SF
        Pi_ij    = - sqrt( (h(i)*n(j)/SF) * (n(i)-n(j)+h(i)*n(j)/SF))

      elseif(i > FL .and. j <= FL) then

        Delta_ij = n(i)*h(j)*(1d0 - SF)/SF
        Pi_ij    = - sqrt( (h(j)*n(i)/SF) * (n(j)-n(i)+h(j)*n(i)/SF))

      elseif(i > FL .and. j >= FL) then

        Delta_ij = n(i)*n(j)
        Pi_ij    = sqrt(n(i)*n(j))

      else

        Delta_ij = 0d0
        Pi_ij    = 0d0

      endif

      EJ_PNOF4 = EJ_PNOF4 - 2d0*Delta_ij*Jint(i,j)
      EK_PNOF4 = EK_PNOF4 + Delta_ij*Kint(i,j)
      EL_PNOF4 = EL_PNOF4 + Pi_ij*Lint(i,j)

    enddo   
  enddo   


! **************
! *** PNOF6d ***
! **************

  EJ_PNOF6d = 0d0
  EK_PNOF6d = 0d0
  EL_PNOF6d = 0d0

  do i=1,nMO
    do j=1,nMO

      if(i == j) then

        Delta_ij = n(i)*n(j)
        Pi_ij    = sqrt(n(i)*n(j))

      elseif(i <= FL .and. j <= FL) then

        Delta_ij = h(i)*h(j)*exp(-2d0*SF)
        Pi_ij    = - sqrt(h(i)*h(j))*exp(-SF)
        
      elseif(i <= FL .and. j > FL) then

        Delta_ij = gam(i)*gam(j)/Sd
        Pi_ij    = - sqrt( (n(i)*h(j) + gam(i)*gam(j)/Sd) * (n(j)*h(i) + gam(i)*gam(j)/Sd))

      elseif(i > FL .and. j <= FL) then

        Delta_ij = gam(i)*gam(j)/Sd
        Pi_ij    = - sqrt( (n(i)*h(j) + gam(i)*gam(j)/Sd) * (n(j)*h(i) + gam(i)*gam(j)/Sd))

      elseif(i > FL .and. j >= FL) then

        Delta_ij = n(i)*n(j)*exp(-2d0*SF)
        Pi_ij    = sqrt(n(i)*n(j))*exp(-SF)

      else

        Delta_ij = 0d0
        Pi_ij    = 0d0

      endif

      EJ_PNOF6d = EJ_PNOF6d - 2d0*Delta_ij*Jint(i,j)
      EK_PNOF6d = EK_PNOF6d + Delta_ij*Kint(i,j)
      EL_PNOF6d = EL_PNOF6d + Pi_ij*Lint(i,j)

    enddo   
  enddo   

! **************
! *** PNOF6u ***
! **************

  EJ_PNOF6u = 0d0
  EK_PNOF6u = 0d0
  EL_PNOF6u = 0d0

  do i=1,nMO
    do j=1,nMO

      if(i == j) then

        Delta_ij = n(i)*n(j)
        Pi_ij    = sqrt(n(i)*n(j))

      elseif(i <= FL .and. j <= FL) then

        Delta_ij = h(i)*h(j)*exp(-2d0*SF)
        Pi_ij    = - sqrt(h(i)*h(j))*exp(-SF)
        
      elseif(i <= FL .and. j > FL) then

        Delta_ij = gam(i)*gam(j)/Su
        Pi_ij    = - sqrt( (n(i)*h(j) + gam(i)*gam(j)/Su) * (n(j)*h(i) + gam(i)*gam(j)/Su))

      elseif(i > FL .and. j <= FL) then

        Delta_ij = gam(i)*gam(j)/Su
        Pi_ij    = - sqrt( (n(i)*h(j) + gam(i)*gam(j)/Su) * (n(j)*h(i) + gam(i)*gam(j)/Su))

      elseif(i > FL .and. j >= FL) then

        Delta_ij = n(i)*n(j)*exp(-2d0*SF)
        Pi_ij    = sqrt(n(i)*n(j))*exp(-SF)

      else

        Delta_ij = 0d0
        Pi_ij    = 0d0

      endif

      EJ_PNOF6u = EJ_PNOF6u - 2d0*Delta_ij*Jint(i,j)
      EK_PNOF6u = EK_PNOF6u + Delta_ij*Kint(i,j)
      EL_PNOF6u = EL_PNOF6u + Pi_ij*Lint(i,j)

    enddo   
  enddo   

! **************
! *** PNOF6h ***
! **************

  EJ_PNOF6h = 0d0
  EK_PNOF6h = 0d0
  EL_PNOF6h = 0d0

  do i=1,nMO
    do j=1,nMO

      if(i == j) then

        Delta_ij = n(i)*n(j)
        Pi_ij    = sqrt(n(i)*n(j))

      elseif(i <= FL .and. j <= FL) then

        Delta_ij = h(i)*h(j)*exp(-2d0*SF)
        Pi_ij    = - sqrt(h(i)*h(j))*exp(-SF)
        
      elseif(i <= FL .and. j > FL) then

        Delta_ij = gam(i)*gam(j)/Sh
        Pi_ij    = - sqrt( (n(i)*h(j) + gam(i)*gam(j)/Sh) * (n(j)*h(i) + gam(i)*gam(j)/Sh))

      elseif(i > FL .and. j <= FL) then

        Delta_ij = gam(i)*gam(j)/Sh
        Pi_ij    = - sqrt( (n(i)*h(j) + gam(i)*gam(j)/Sh) * (n(j)*h(i) + gam(i)*gam(j)/Sh))

      elseif(i > FL .and. j >= FL) then

        Delta_ij = n(i)*n(j)*exp(-2d0*SF)
        Pi_ij    = sqrt(n(i)*n(j))*exp(-SF)

      else

        Delta_ij = 0d0
        Pi_ij    = 0d0

      endif

      EJ_PNOF6h = EJ_PNOF6h - 2d0*Delta_ij*Jint(i,j)
      EK_PNOF6h = EK_PNOF6h + Delta_ij*Kint(i,j)
      EL_PNOF6h = EL_PNOF6h + Pi_ij*Lint(i,j)

    enddo   
  enddo   

! Add the SD part

  EJ_PNOF2  = EJ_SD + EJ_PNOF2 
  EJ_PNOF3  = EJ_SD + EJ_PNOF3 
  EJ_PNOF4  = EJ_SD + EJ_PNOF4 
! EJ_PNOF5  = EJ_SD + EJ_PNOF5 
  EJ_PNOF6d = EJ_SD + EJ_PNOF6d
  EJ_PNOF6u = EJ_SD + EJ_PNOF6u
  EJ_PNOF6h = EJ_SD + EJ_PNOF6h
! EJ_PNOF7  = EJ_SD + EJ_PNOF7 

  EK_PNOF2  = EK_SD + EK_PNOF2 
  EK_PNOF3  = EK_SD + EK_PNOF3 
  EK_PNOF4  = EK_SD + EK_PNOF4 
! EK_PNOF5  = EK_SD + EK_PNOF5 
  EK_PNOF6d = EK_SD + EK_PNOF6d
  EK_PNOF6u = EK_SD + EK_PNOF6u
  EK_PNOF6h = EK_SD + EK_PNOF6h
! EK_PNOF7  = EK_SD + EK_PNOF7 

! Compute total energies

  E_PNOF2  = ET + EV + EJ_PNOF2  + EK_PNOF2  + EL_PNOF2
  E_PNOF3  = ET + EV + EJ_PNOF3  + EK_PNOF3  + EL_PNOF3
  E_PNOF4  = ET + EV + EJ_PNOF4  + EK_PNOF4  + EL_PNOF4
! E_PNOF5  = ET + EV + EJ_PNOF5  + EK_PNOF5  + EL_PNOF5
  E_PNOF6d = ET + EV + EJ_PNOF6d + EK_PNOF6d + EL_PNOF6d
  E_PNOF6u = ET + EV + EJ_PNOF6u + EK_PNOF6u + EL_PNOF6u
  E_PNOF6h = ET + EV + EJ_PNOF6h + EK_PNOF6h + EL_PNOF6h
! E_PNOF7  = ET + EV + EJ_PNOF7  + EK_PNOF7  + EL_PNOF7

! Dump energies

  print*, '*******************************'
  print*, '*** JKL NOFT    functionals ***'
  print*, '*******************************'
  print*, ''
  print*, '*** Coulomb        energies ***'
  print*, 'Coulomb  PNOF2           energy = ',EJ_PNOF2
  print*, 'Coulomb  PNOF3           energy = ',EJ_PNOF3
  print*, 'Coulomb  PNOF4           energy = ',EJ_PNOF4
! print*, 'Coulomb  PNOF5           energy = ',EJ_PNOF5
  print*, 'Coulomb  PNOF6d          energy = ',EJ_PNOF6d
  print*, 'Coulomb  PNOF6u          energy = ',EJ_PNOF6u
  print*, 'Coulomb  PNOF6h          energy = ',EJ_PNOF6h
! print*, 'Coulomb  PNOF7           energy = ',EJ_PNOF7
  print*, ''
  print*, '*** Exchange       energies ***'
  print*, 'Exchange PNOF2           energy = ',EK_PNOF2
  print*, 'Exchange PNOF3           energy = ',EK_PNOF3
  print*, 'Exchange PNOF4           energy = ',EK_PNOF4
! print*, 'Exchange PNOF5           energy = ',EK_PNOF5
  print*, 'Exchange PNOF6d          energy = ',EK_PNOF6d
  print*, 'Exchange PNOF6u          energy = ',EK_PNOF6u
  print*, 'Exchange PNOF6h          energy = ',EK_PNOF6h
! print*, 'Exchange PNOF7           energy = ',EK_PNOF7
  print*, ''
  print*, '*** Time-inversion energies ***'
  print*, 'Time-inversion PNOF2     energy = ',EL_PNOF2
  print*, 'Time-inversion PNOF3     energy = ',EL_PNOF3
  print*, 'Time-inversion PNOF4     energy = ',EL_PNOF4
! print*, 'Time-inversion PNOF5     energy = ',EL_PNOF5
  print*, 'Time-inversion PNOF6d    energy = ',EL_PNOF6d
  print*, 'Time-inversion PNOF6u    energy = ',EL_PNOF6u
  print*, 'Time-inversion PNOF6h    energy = ',EL_PNOF6h
! print*, 'Time-inversion PNOF7     energy = ',EL_PNOF7
  print*, ''
  print*, '*** Two-electron   energies ***'
  print*, 'J+K+L          PNOF2     energy = ',EJ_PNOF2  + EK_PNOF2  + EL_PNOF2
  print*, 'J+K+L          PNOF3     energy = ',EJ_PNOF3  + EK_PNOF3  + EL_PNOF3
  print*, 'J+K+L          PNOF4     energy = ',EJ_PNOF4  + EK_PNOF4  + EL_PNOF4
! print*, 'J+K+L          PNOF5     energy = ',EJ_PNOF5  + EK_PNOF5  + EL_PNOF5
  print*, 'J+K+L          PNOF6d    energy = ',EJ_PNOF6d + EK_PNOF6d + EL_PNOF6d
  print*, 'J+K+L          PNOF6u    energy = ',EJ_PNOF6u + EK_PNOF6u + EL_PNOF6u
  print*, 'J+K+L          PNOF6h    energy = ',EJ_PNOF6h + EK_PNOF6h + EL_PNOF6h
! print*, 'J+K+L          PNOF7     energy = ',EJ_PNOF7  + EK_PNOF7  + EL_PNOF7
  print*, ''
  print*, '*** Total          energies ***'
  print*, 'Total    PNOF2           energy = ',E_PNOF2
  print*, 'Total    PNOF3           energy = ',E_PNOF3
  print*, 'Total    PNOF4           energy = ',E_PNOF4
! print*, 'Total    PNOF5           energy = ',E_PNOF5
  print*, 'Total    PNOF6d          energy = ',E_PNOF6d
  print*, 'Total    PNOF6u          energy = ',E_PNOF6u
  print*, 'Total    PNOF6h          energy = ',E_PNOF6h
! print*, 'Total    PNOF7           energy = ',E_PNOF7
  print*, ''

end subroutine NOFT_JKLfunc

