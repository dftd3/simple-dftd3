module d3mod_main
   implicit none

contains

subroutine main_run(env, mol)
   use d3def_environment
   use d3def_molecule
   class(d3_environment), intent(inout) :: env
   class(d3_molecule), intent(inout) :: mol

   call mol%print_info(env%unit)

end subroutine main_run

end module d3mod_main
