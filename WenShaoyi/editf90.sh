#!/bin/bash
rm .*.swp
rm tags
ctags *
vim -c "edit dfd_main.f90|tabe gems.inp|tabe dfd_fv.f90|tabe dfd_output.f90|vsp dfd_input.f90| \
  tabe dfd_linear.f90|vsp dfd_precon.f90|tabe dfd_disc.f90|vsp dfd_bound.f90|tabe dfd_data.f90| \
  tabe dfd_source.f90|tabe dfd_type.f90|vsp dfd_state.f90|tabn|tabn"
# gvim -c "edit dfd_main.f90|tabe dfd.inp|tabe dfd_fv.f90|tabe dfd_output.f90|vsp dfd_input.f90| \
#   tabe dfd_linear.f90|vsp dfd_precon.f90|tabe dfd_disc.f90|vsp dfd_bound.f90|tabe dfd_data.f90| \
#   tabe dfd_source.f90|tabe dfd_type.f90|vsp dfd_state.f90"
