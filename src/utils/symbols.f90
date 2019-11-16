module d3mod_utils_symbols
   use d3par_pse
   implicit none
   public :: symbol_to_number, assignment(=)

   interface assignment(=)
      module procedure :: symbol_to_number
   end interface assignment(=)

contains

pure subroutine symbol_to_number(number, symbol)
   character(len=*), intent(in) :: symbol
   integer, intent(out) :: number
   character(len=2) :: lc_symbol
   integer :: i, j, k, l
   integer, parameter :: offset = iachar('a')-iachar('A')

   number = 0
   lc_symbol = '  '

   k = 0
   do j = 1, len_trim(symbol)
      if (k > 2) exit
      l = iachar(symbol(j:j))
      if (k >= 1 .and. l == iachar(' ')) exit
      if (k >= 1 .and. l == 9) exit
      if (l >= iachar('A') .and. l <= iachar('Z')) l = l + offset
      if (l >= iachar('a') .and. l <= iachar('z')) then
         k = k+1
         lc_symbol(k:k) = achar(l)
      endif
   enddo

   do i = 1, size(pse)
      if (lc_symbol == pse(i)) then
         number = i
         exit
      endif
   enddo

end subroutine symbol_to_number

end module d3mod_utils_symbols
