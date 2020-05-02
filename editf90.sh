#!/bin/bash
rm .*.swp
rm tags
ctags *
vim -c "edit gems_main.f90|tabe gems.inp|tabe gems_fv.f90|tabe gems_output.f90|vsp gems_input.f90| \
  tabe gems_linear.f90|vsp gems_precon.f90|tabe gems_disc.f90|vsp gems_bound.f90|tabe gems_data.f90| \
  tabe gems_source.f90|tabe gems_type.f90|vsp gems_state.f90"
# gvim -c "edit gems_main.f90|tabe gems.inp|tabe gems_fv.f90|tabe gems_output.f90|vsp gems_input.f90| \
#   tabe gems_linear.f90|vsp gems_precon.f90|tabe gems_disc.f90|vsp gems_bound.f90|tabe gems_data.f90| \
#   tabe gems_source.f90|tabe gems_type.f90|vsp gems_state.f90"
