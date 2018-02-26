BEGIN_PROVIDER [ logical, initialize_dress_E0_denominator ]
 implicit none
 BEGIN_DOC
 ! If true, initialize dress_E0_denominator
 END_DOC
 initialize_dress_E0_denominator = .True.
END_PROVIDER

BEGIN_PROVIDER [ double precision, dress_E0_denominator, (N_states) ]
 implicit none
 BEGIN_DOC
 ! E0 in the denominator of the dress
 END_DOC
 if (initialize_dress_E0_denominator) then
  dress_E0_denominator(1:N_states) = psi_energy(1:N_states)
! dress_E0_denominator(1:N_states) = HF_energy - nuclear_repulsion
! dress_E0_denominator(1:N_states) = barycentric_electronic_energy(1:N_states)
  call write_double(6,dress_E0_denominator(1)+nuclear_repulsion, 'dress Energy denominator')
 else
   dress_E0_denominator = -huge(1.d0)
 endif
END_PROVIDER

