! This file is part of s-dftd3.
! SPDX-Identifier: LGPL-3.0-or-later
!
! s-dftd3 is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! s-dftd3 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with s-dftd3.  If not, see <https://www.gnu.org/licenses/>.

module dftd3_citation
   implicit none
   private

   public :: citation_type, author_name, new_citation
   public :: format_bibtex

   type :: author_name
      character(len=:), allocatable :: name
   end type author_name

   interface author_name
      module procedure :: new_author
   end interface author_name

   type :: citation_type
      character(len=:), allocatable :: title
      type(author_name), allocatable :: author(:)
      character(len=:), allocatable :: issue
      character(len=:), allocatable :: journal
      character(len=:), allocatable :: volume
      character(len=:), allocatable :: year
      character(len=:), allocatable :: pages
      character(len=:), allocatable :: doi
   end type citation_type

   character(len=*), parameter :: nl = new_line('a')

contains

pure function new_citation(title, author, journal, issue, volume, year, pages, doi) result(citation)
   character(len=*), intent(in) :: title
   type(author_name), intent(in) :: author(:)
   character(len=*), intent(in) :: journal
   character(len=*), intent(in), optional :: issue
   character(len=*), intent(in) :: volume
   character(len=*), intent(in) :: year
   character(len=*), intent(in) :: pages
   character(len=*), intent(in) :: doi
   type(citation_type) :: citation

   citation%title = title
   citation%author = author
   citation%journal = journal
   citation%volume = volume
   citation%year = year
   citation%pages = pages
   citation%doi = doi
   if (present(issue)) citation%issue = issue
end function new_citation

pure function new_author(name) result(author)
   character(len=*), intent(in) :: name
   type(author_name) :: author
   author%name = name
end function new_author

subroutine format_bibtex(string, citation)
   character(len=:), allocatable, intent(out) :: string
   class(citation_type), intent(in) :: citation
   integer :: idx

   if (.not.(allocated(citation%doi) &
      & .and.allocated(citation%title) &
      & .and.allocated(citation%author) &
      & .and.allocated(citation%journal) &
      & .and.allocated(citation%volume) &
      & .and.allocated(citation%year) &
      & .and.allocated(citation%pages) &
      )) then
      return
   end if

   string = &
      & "@article{" // citation%doi // "," // nl // &
      & "  title = {{" // citation%title // "}}," // nl // &
      & "  author = {" // citation%author(1)%name // nl
   do idx = 2, size(citation%author)
      string = string // &
         & "    and " // citation%author(idx)%name // nl
   end do
   string = string // &
      & "  }," // nl
   if (allocated(citation%issue)) then
      string = string // &
         & "  issue = {" // citation%issue // "}," // nl
   end if
   string = string // &
      & "  volume = {" // citation%volume // "}," // nl // &
      & "  pages = {" // citation%pages // "}," // nl // &
      & "  doi = {" // citation%doi // "}," // nl // &
      & "  url = {https://doi.org/" // citation%doi // "}" // nl // &
      & "}"
end subroutine format_bibtex

end module dftd3_citation