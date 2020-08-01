!
! List of all dynamical fields
! (so we don't need to worry about ordering of compilation)
!
module mod_fields_list
   use mod_field, only: field

   type(field) :: psi4_lin_p, psi4_lin_q, psi4_lin_f, &
                  res_lin_q, &

                  psi4_scd_p, psi4_scd_q, psi4_scd_f, &
                  res_scd_q, &

                  psi3, psi2, la, pi, muhll, hlmb, hmbmb, &

                  res_bianchi3, res_bianchi2, res_hll

end module mod_fields_list
