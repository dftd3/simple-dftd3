module d3mod_help
   implicit none
   public :: helpmsg

contains

subroutine helpmsg(unit)
   integer, intent(in) :: unit
   include 's-dftd3-version.fh'
   write(unit, '(a,1x,a,1x,a)') &
      "Usage:", name, "[options] <input>"
end subroutine helpmsg

end module d3mod_help
