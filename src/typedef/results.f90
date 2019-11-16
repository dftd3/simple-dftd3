module d3def_results
   use iso_fortran_env, only: wp => real64
   implicit none
   public :: d3_results
   private


   type :: d3_results
      real(wp), allocatable :: energies(:)
      real(wp), allocatable :: gradient(:, :)
      real(wp), allocatable :: stress(:, :)
   end type d3_results


end module d3def_results
