! This file is part of s-dftd3.
! Based on the DFTB+ implementation: https://github.com/dftbplus/dftbplus
!
! Copyright (C) 2006 - 2019  DFTB+ developers group
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


!> Implements dynamic neighbour list with iterator, like present in DFTB+.
!
!  The dynamic neighbour list does not store the entire neighbour list,
!  but creates it on the fly, allowing for a low memory footprint for large
!  neighbour lists (at the cost of speed).
!
!  To use iterate over a neighbour list, use something similar like this snippet:
!
!  ```fortran
!  type(d3_molecule) :: mol
!  type(d3_neighbourlist), target :: neighlist
!  type(d3_neighlist_iterator) :: neighiter
!  integer, parameter :: iter_chunk_size = 1000
!  integer :: iat, jat, ij, neighs
!  integer :: image(iter_chunk_size)
!  real(wp) :: rij(3), r2
!  real(wp) :: dists2(iter_chunk_size), coords(3, iter_chunk_size)
!  do iat = 1, len(mol)
!     call neighlist%get_iterator(neighiter, iat)
!     neighs = iter_chunk_size
!     do while(neighs == iter_chunk_size)
!        call neighIter%next(neighs, coords=coords, image=image, dists2=dist2)
!        do ij = 1, neighs
!           jat = image(ij)
!           rij = mol%xyz(:, iat) - coords(:, ij)
!           r2 = dists2(ij)
!           ...
!        enddo
!     enddo
!  enddo
!  ```
!
!  More complicated is the case if you want to iterate over triples using the
!  dynamic neighbour list. It is necessary to have *two* iterators like here
!
!  ```fortran
!  type(d3_molecule) :: mol
!  type(d3_neighbourlist), target :: neighlist
!  type(d3_neighlist_iterator) :: jneighiter, kneighiter
!  integer, parameter :: iter_chunk_size = 1000
!  integer :: iat, jat, kat, ij, jk, jneighs, kneighs
!  integer :: jimage(iter_chunk_size), kimage(iter_chunk_size)
!  real(wp) :: rij(3), rik(3), rjk(3), rij2, rik2, rjk2
!  real(wp) :: jdists2(iter_chunk_size), jcoords(3, iter_chunk_size)
!  real(wp) :: kdists2(iter_chunk_size), kcoords(3, iter_chunk_size)
!  do iat = 1, len(mol)
!     call neighlist%get_iterator(jneighiter, iat)
!     jneighs = iter_chunk_size
!     do while(jneighs == iter_chunk_size)
!        call jneighiter%next(jneighs, coords=jcoords, image=jimage, &
!           &                 dists2=jdist2)
!        do ij = 1, jneighs
!           jat = jimage(ij)
!           rij = mol%xyz(:, iat) - jcoords(:, ij)
!           rij2 = jdists2(ij)
!           kneighs = iter_chunk_size
!           call neighlist%get_iterator(kneighiter, jat)
!           do while(kneighs == iter_chunk_size)
!              call kneighiter%next(kneighs, coords=kcoords, image=kimage, &
!                 &                 dists2=kdist2)
!              do jk = 1, kneighs
!                 kat = kimage(jk)
!                 rjk = jcoords(:, ij) - kcoords(:, jk)
!                 rjk2 = kdists2(jk)
!                 rik = rij + rjk
!                 rik2 = rik(1)**2 + rik(2)**2 + rik(3)**2
!                 ...
!              enddo
!           enddo
!        enddo
!     enddo
!  enddo
!  ```
module d3def_neighbourlist
   use iso_fortran_env, only: wp => real64
   use d3def_latp_generator
   implicit none
   private

   public :: d3_neighbourlist
   public :: d3_neighlist_iterator
   public :: iter_chunk_size

   integer, parameter :: iter_chunk_size = 512

   !> Dynamic neighbour list
   type :: d3_neighbourlist
      private
      !> Cutoff for neighbour generation
      real(wp) :: cutoff
      !> Nr. of atoms 
      integer :: nAtom
      !> Coordinates of atoms (folded into unit cell, if periodic)
      real(wp), allocatable :: coords0(:, :)
      !> Whether system is periodic
      logical :: periodic
      !> Lattice vectors, if system is periodic
      real(wp), allocatable :: lattice(:, :)
      !> Inverse lattice vectors, if system is periodic
      real(wp), allocatable :: inv_lat(:, :)
   contains
      procedure :: new => neighlist_new
      procedure :: set_lattice => neighlist_set_lattice
      procedure :: set_coords => neighlist_set_coords
      procedure :: set_cutoff => neighlist_set_cutoff
      procedure, pass(neighlist) :: get_iterator => neighiter_new
   end type d3_neighbourlist


   !> Iterator over a dynamic neighbour list
   type :: d3_neighlist_iterator
      private
      !> Pointer to the original neighbour list
      type(d3_neighbourlist), pointer :: neighList
      !> Lattice point generator (if system is periodic)
      type(d3_latp_generator), allocatable :: latpgen
      !> Whether system is periodic
      logical :: periodic
      !> Lattice vectors, if system is periodic
      real(wp), allocatable :: lattice(:, :)
      !> Number of atoms
      integer :: nAtom
      !> Atom for which neighbours are returned by the iterator
      integer :: iAtom1
      !> Coordinates of the atom
      real(wp) :: coordsAtom1(3)
      !> Neighbour atom to be returned as next
      integer :: iAtom2
      !> Current cell shift vector (if system is periodic)
      real(wp) :: cellVec(3)
      !> Square of the cutoff radius
      real(wp) :: cutoff2
      !> Whether iterator has finished
      logical :: finished
   contains
      procedure :: new => neighiter_new
      procedure :: next => neighiter_next
   end type d3_neighlist_iterator


contains

!> Initializes a dynamic neighbour list.
subroutine neighlist_new(self, cutoff, nAtom, coords, lattice, periodic)
   !> Initialized instance on exit.
   class(d3_neighbourlist), intent(out) :: self
   !> Cutoff up to which the neighbours should be generated
   real(wp), intent(in) :: cutoff
   !> Nr. of atoms in the system.
   integer, intent(in) :: nAtom
   !> New coordinates.
   real(wp), intent(in), optional :: coords(:, :)
   !> New lattice vectors.
   real(wp), intent(in), optional :: lattice(:, :)
   !> Whether the system is periodic
   logical, intent(in) :: periodic

   self%cutoff = cutoff
   self%nAtom = nAtom
   self%periodic = periodic
   allocate(self%coords0(3, self%nAtom))
   if (present(coords)) then
      call self%set_coords(coords)
   endif

   if (self%periodic) then
      allocate(self%lattice(3, 3))
      allocate(self%inv_lat(3, 3))
      if (present(lattice)) then
         call self%set_lattice(lattice)
      end if
   end if


end subroutine neighlist_new


!> Updates the lattice vectors.
subroutine neighlist_set_lattice(self, lattice)
   use d3mod_utils_periodic
   !> Instance.
   class(d3_neighbourlist), intent(inout) :: self
   !> New lattice vectors.
   real(wp), intent(in) :: lattice(:, :)

   if (self%periodic) then
      self%lattice = lattice
      self%inv_lat = mat_inv_3x3(lattice)
   endif

end subroutine neighlist_set_lattice

!> Updates the coordiantes.
subroutine neighlist_set_coords(self, coords)
   !> Instance.
   class(d3_neighbourlist), intent(inout) :: self
   !> New coordinates.
   real(wp), intent(in) :: coords(:, :)

   self%coords0 = coords

end subroutine neighlist_set_coords

!> Updates the cutoff radius.
subroutine neighlist_set_cutoff(self, cutoff)
   use d3mod_utils_periodic
   !> Instance.
   class(d3_neighbourlist), intent(inout) :: self
   !> New lattice vectors.
   real(wp), intent(in) :: cutoff

   if (cutoff > 0.0_wp) then
      self%cutoff = cutoff
   end if

end subroutine neighlist_set_cutoff


!> Initializes an iterator for the dynamic neighbours of a given atom.
subroutine neighiter_new(self, neighList, iAtom, includeSelf)
   !> Initialized instance on exit.
   class(d3_neighlist_iterator), intent(out) :: self
   !> Dynamic neighbour list containing the basic data
   class(d3_neighbourlist), target, intent(in) :: neighList
   !> Index of the atom for which the neighbours should be generated.
   integer, intent(in) :: iAtom
   !> Whether the atom itself should be also returned as first neighbour
   !  (default: false)
   logical, intent(in), optional :: includeSelf

   logical :: includeSelf0

   if (present(includeSelf)) then
      includeSelf0 = includeSelf
   else
      includeSelf0 = .false.
   end if

   self%neighList => neighList
   self%cutoff2 = self%neighList%cutoff**2
   self%nAtom = self%neighList%nAtom
   self%periodic = self%neighList%periodic

   self%cellVec = 0.0_wp
   self%iAtom1 = iAtom
   self%coordsAtom1 = self%neighList%coords0(:,iAtom)
   if (includeSelf0) then
      self%iAtom2 = self%iAtom1
   else
      self%iAtom2 = self%iAtom1 + 1
   end if

   if (self%periodic) then
      self%lattice = self%neighList%lattice
      allocate(self%latpgen)
      call self%latpgen%new(self%lattice, self%neighList%inv_lat, &
         &                  self%neighList%cutoff, pos_extension=1, &
         &                  neg_extension=1, exclude_origin=.true.)
   end if

   self%finished = .false.

end subroutine neighiter_new


!> Returns the next group of neighbours.
pure subroutine neighiter_next(self, nNeighbour, coords, dists2, image)
   !> Instance.
   class(d3_neighlist_iterator), intent(inout) :: self
   !> Nr. of neighbours to return on entry, nr. of neighbours found on exit.
   !  When entry and exit values differ, the iterator has finished and can not
   !  return any more neighbours.
   integer, intent(inout) :: nNeighbour
   !> Coordinates of the neighbours. Shape: (3, nNeighbour)
   real(wp), intent(out), optional :: coords(:, :)
   !> Distances of the neighbours. Shape: (nNeighbour)
   real(wp), intent(out), optional :: dists2(:)
   !> Correspoinding images of the neighbours in the central cell.
   !  Shape: (nNeighbour)
   integer, intent(out), optional :: image(:)

   integer :: maxNeighs
   integer :: imageTmp(nNeighbour)
   real(wp) :: dist2
   real(wp) :: neighCoords(3)
   real(wp) :: coordsTmp(3, nNeighbour), distsTmp(nNeighbour)

   maxNeighs = nNeighbour
   nNeighbour = 0
   if (self%finished) return

   do while (nNeighbour < maxNeighs)
      if (self%iAtom2 > self%nAtom) then
         if (self%periodic) then
            call self%latpgen%get_next_point(self%cellVec, self%finished)
            self%cellVec = matmul(self%lattice, self%cellVec)
         else
            self%finished = .true.
         end if
         if (self%finished) exit
         self%iAtom2 = self%iAtom1
      end if

      neighCoords = self%neighList%coords0(:, self%iAtom2) + self%cellVec
      dist2 = sum((self%coordsAtom1 - neighCoords)**2)
      if (dist2 <= self%cutoff2) then
         nNeighbour = nNeighbour + 1
         coordsTmp(:, nNeighbour) = neighCoords
         distsTmp(nNeighbour) = dist2
         imageTmp(nNeighbour) = self%iAtom2
      end if
      self%iAtom2 = self%iAtom2 + 1
   end do

   if (present(coords)) then
      coords(:, 1:nNeighbour) = coordsTmp(:, 1:nNeighbour)
   end if
   if (present(dists2)) then
      dists2(1:nNeighbour) = distsTmp(1:nNeighbour)
   end if
   if (present(image)) then
      image(1:nNeighbour) = imageTmp(1:nNeighbour)
   end if

end subroutine neighiter_next

end module d3def_neighbourlist
