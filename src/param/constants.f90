module d3par_constants
   use iso_fortran_env, only: wp => real64
   implicit none

   !> ratio of a circle's circumference to its diameter
   real(wp), parameter :: pi = 4.0_wp * atan(1.0_wp)
   !> √π
   real(wp), parameter :: sqrtpi = sqrt(pi)
   !> 2×π
   real(wp), parameter :: twopi = 2.0_wp * pi
   !> 4×π
   real(wp), parameter :: fourpi = 4.0_wp * pi
   !> π/2
   real(wp), parameter :: pihalf = 0.5_wp * pi

end module d3par_constants
