

 BEGIN_PROVIDER [ double precision, dressing_column_h, (N_det,N_states) ]
&BEGIN_PROVIDER [ double precision, dressing_column_s, (N_det,N_states) ]
 implicit none
 BEGIN_DOC
 ! Null dressing vectors
 END_DOC

 integer :: i,ii,k,j,jj, l
 double precision :: f, tmp
 double precision, external :: u_dot_v

 dressing_column_h(:,:) = 0.d0
 dressing_column_s(:,:) = 0.d0
 do k=1,N_states
   l = dressed_column_idx(k)
   f = -1.d0/psi_coef(l,k)
   do jj=1,N_det_non_ref
     j = idx_non_ref(jj)
     dressing_column_h(j,k) = 2.d0*delta_ij   (k,jj) 
     dressing_column_s(j,k) = 2.d0*delta_ij_s2(k,jj)
     dressing_column_h(l,k) += psi_coef(j,k) * delta_ij(k,jj)
     dressing_column_s(l,k) += psi_coef(j,k) * delta_ij_s2(k,jj)
   enddo
   dressing_column_h(l,k) *= f
   dressing_column_s(l,k) *= f
 enddo
END_PROVIDER

