program shifted_bk
  implicit none
  BEGIN_DOC
! TODO
  END_DOC
  
  PROVIDE psi_bilinear_matrix_columns_loc psi_det_alpha_unique psi_det_beta_unique
  PROVIDE psi_bilinear_matrix_rows psi_det_sorted_gen_from_det_sorted psi_bilinear_matrix_order
  PROVIDE psi_bilinear_matrix_transp_rows_loc psi_bilinear_matrix_transp_columns
  PROVIDE psi_bilinear_matrix_transp_order


  !call diagonalize_CI()
  call dress_zmq()
end

