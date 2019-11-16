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

module d3mod_header
   implicit none
   public :: sdftd3_header, sdftd3_version, generic_header
   private

contains

subroutine sdftd3_header(unit)
   integer, intent(in) :: unit
   write(unit, '(32("="),">",1x,6(1x,a),2x,"<",32("="))') &
      & ['s','D','F','T','D','3']
   call sdftd3_version(unit)
end subroutine sdftd3_header

subroutine sdftd3_version(unit)
   integer, intent(in) :: unit
   include 's-dftd3-version.fh'
   write(unit, '(1x,"*",*(1x,a))') &
      & name, "version", version, "compiled by", author, "on", date
end subroutine sdftd3_version

subroutine generic_header(unit, string, width, offset)
   integer, intent(in) :: unit
   integer, intent(in) :: offset
   integer, intent(in) :: width
   character(len=*), intent(in) :: string
   character(len=width) :: dum1,dum2
   character(len=2*width) :: outstring
   character(len=width) :: formatstr
   integer :: strlen, ifront, iback
   strlen = len(string)
   ifront = (width - strlen)/2
   iback  = width - ifront - strlen
   write(dum1,*) width
   write(dum2,*) offset
   write(formatstr,'(i0,"x,a,",i0,"x")') ifront, iback
   write(outstring,'("|",'//formatstr//',"|")') string
   write(unit,'('//dum2//'x,1x,'//dum1//'("-"),1x)')
   write(unit,'('//dum2//'x,a)') trim(outstring)
   write(unit,'('//dum2//'x,1x,'//dum1//'("-"),1x)')
end subroutine generic_header

end module d3mod_header
