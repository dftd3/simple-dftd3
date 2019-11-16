! This file is part of s-dftd3.
!
! Copyright (C) 2019 Sebastian Ehlert
!
! s-dftd3 is free software: you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! s-dftd3 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with s-dftd3.  If not, see <https://www.gnu.org/licenses/>.

submodule(d3def_molecule) molecule_reader
   use d3def_molecule
   implicit none

   interface assignment(=)
      module procedure :: symbol_to_number
   end interface assignment(=)

contains

module subroutine read_molecule_generic(self, unit, format)
   use d3mod_file_utils
   class(d3_molecule), intent(out) :: self
   integer, intent(in) :: unit
   integer, intent(in), optional :: format
   integer :: ftype
   logical :: status
   character(len=:), allocatable :: message
   if (present(format)) then
      ftype = format
   else
      ftype = self%ftype
   endif

   select case(ftype)
   case(p_ftype%xyz)
      call read_molecule_xyz(self, unit, status, iomsg=message)
   case default
      status = .false.
      message = "coordinate format unknown"
   end select

   if (.not.status) then
      self%ftype = -1
      self%name = message
   else
      self%ftype = ftype
   endif

end subroutine read_molecule_generic


subroutine read_molecule_xyz(mol, unit, status, iomsg)
   use iso_fortran_env, wp => real64
   use d3par_atomic_units
   logical, parameter :: debug = .false.
   class(d3_molecule), intent(out) :: mol
   integer, intent(in) :: unit
   logical, intent(out) :: status
   character(len=:), allocatable, intent(out) :: iomsg
   integer  :: n, iat
   real(wp) :: xyz(3)
   real(wp) :: conv

   character(len=:),allocatable :: message
   character(len=:),allocatable :: line
   character(len=10) :: chdum
   integer  :: err

   status = .false.

   conv = aatoau

   read(unit, *, iostat=err) n
   if (err /= 0) then
      iomsg = "Could not read number of atoms, check format!"
      return
   endif

   if (n < 1) then
      iomsg = "Found no atoms, cannot work without atoms!"
      return
   endif

   call mol%new(n)
   mol%npbc = 0 ! Xmol is always molecular (there are extensions to this...)
   mol%pbc = .false.

   ! drop next record
   read(unit, '(a)')

   n = 0
   do while (n < mol%n)
      read(unit, *, iostat=err) chdum, xyz(1), xyz(2), xyz(3)
      if (is_iostat_end(err)) exit
      if (err /= 0) then
         iomsg = "Could not parse coordinates from Xmol file"
         return
      endif
      if (debug) print'("->",a)',chdum
      if (debug) print'("->",3g0)',xyz

      iat = chdum
      if (debug) print'("->",g0)',iat
      if (iat > 0) then
         n = n+1
         mol%at(n) = iat
         mol%sym(n) = trim(chdum)
         mol%xyz(:,n) = xyz*conv
      else
         iomsg = "Unknown element symbol: '"//trim(chdum)//"'"
         return
      endif
   enddo

   if (n /= mol%n) then
      iomsg = "Atom number missmatch in Xmol file"
      return
   endif

   status = .true.
end subroutine read_molecule_xyz


pure elemental subroutine symbol_to_number(number, symbol)
   use d3par_pse
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

end submodule molecule_reader
