use bitmasks

BEGIN_PROVIDER [ integer, N_det_generators ]
 implicit none
 BEGIN_DOC
 ! For Single reference wave functions, the number of generators is 1 : the
 ! Hartree-Fock determinant
 END_DOC
 integer :: i
 double precision :: norm
 call write_time(output_determinants)
 norm = 0.d0
 N_det_generators = N_det
 do i=1,N_det
   norm = norm + psi_average_norm_contrib_sorted(i)
   if (norm >= threshold_generators) then
     N_det_generators = i
     exit
   endif
 enddo
 N_det_generators = max(N_det_generators,1)
 call write_int(output_determinants,N_det_generators,'Number of generators')
END_PROVIDER

 BEGIN_PROVIDER [ integer(bit_kind), psi_det_generators, (N_int,2,psi_det_size) ]
&BEGIN_PROVIDER [ double precision, psi_coef_generators, (psi_det_size,N_states) ]
 implicit none
 BEGIN_DOC
 ! For Single reference wave functions, the generator is the
 ! Hartree-Fock determinant
 END_DOC
 integer                        :: i, k
 psi_coef_generators = 0.d0
 psi_det_generators  = 0_bit_kind
 do i=1,N_det_generators
   do k=1,N_int
     psi_det_generators(k,1,i) = psi_det_sorted(k,1,i)
     psi_det_generators(k,2,i) = psi_det_sorted(k,2,i)
   enddo
   psi_coef_generators(i,:) = psi_coef_sorted(i,:)
 enddo

END_PROVIDER

BEGIN_PROVIDER [integer, degree_max_generators]
 implicit none
 BEGIN_DOC
! Max degree of excitation (respect to HF) of the generators
 END_DOC
 integer :: i,degree
  degree_max_generators = 0
  do i = 1, N_det_generators
   call get_excitation_degree(HF_bitmask,psi_det_generators(1,1,i),degree,N_int)
   if(degree .gt. degree_max_generators)then
    degree_max_generators = degree
   endif
  enddo
END_PROVIDER 

BEGIN_PROVIDER [ integer, size_select_max]
 implicit none
 BEGIN_DOC
 ! Size of the select_max array
 END_DOC
 size_select_max = 10000
END_PROVIDER

BEGIN_PROVIDER [ double precision, select_max, (size_select_max) ]
 implicit none
 BEGIN_DOC
 ! Memo to skip useless selectors
 END_DOC
 select_max = huge(1.d0)
END_PROVIDER

