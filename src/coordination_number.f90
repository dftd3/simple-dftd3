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

!> Implementation of the coordination number describing the local chemical
!  environment.
!
!  A cutoff radius of 40 Bohr is usually totally sufficient to converge the CN in
!  real space. The exp-CN was originally used in DFT-D3 but tends to be rather
!  long-ranged, such that it can accumulate a small but significant contribution
!  from the second and third coordination shell, while the erf-CN will rapidly
!  drop to zero.
!
!  This leads to the consequence that one can usally savely divide through the
!  exp-CN (even if it is discouraged) since it almost always has some value,
!  While it certainly will result in a divide by zero exception for the erf-CN.
!
!  To use the coordination in an algorithm generate it somewhat similar to this
!  snippet:
!
!  ```fortran
!  type(d3_molecule) :: mol
!  type(d3_neighbourlist) :: neighlist
!  real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :)
!  allocate(cn(len(mol)), dcndr(3, len(mol), len(mol)), dcndL(3, 3, len(mol)))
!  ...
!  call get_coordination_number(mol, neighlist, d3_cn_type%erf, &
!     &                         cn, dcndr, dcndL)
!  call cut_coordination_number(cn, dcndr, dcndL, cn_max=8.0_wp)
!  ```
!
!  Large coordination numbers can be removed by the cutting function, note that
!  the coordination number derivatives loose all symmetry after the cutting,
!  algorithms exploiting this symmetry (for no good reason) will fail with
!  a cutted CN.
!
!  To calculate derivatives of the coordination number use the provided
!  derivatives w.r.t. to cartesian displacements (`dcndr`) and strain
!  deformations (`dcndL`). The strain derivative does not include any volume
!  normalization! The last (third) dimension can be contracted with
!  an energy derivative w.r.t. the coordination number to result in the correct
!  gradient and stress contributions:
!
!  ```fortran
!  gradient = gradient + reshape(matmul(reshape(dcndr, [3*nat, nat]), &
!     &                                 dEdcn), shape(gradient))
!  if (mol%npbc > 0) then
!     sigma = sigma + reshape(matmul(reshape(dcndL, [9, nat]), &
!        &                           dEdcn), shape(sigma))
!  endif
!  ```
!
!  or using BLAS and exploiting that arrays are stored continous in memory
!  (if they are not, than this will fail):
!
!  ```fortran
!  call dgemv('n', 3*len(mol), len(mol), 1.0_wp, dcndr, 3*len(mol), dEdcn, 1, &
!     &       1.0_wp, gradients, 1)
!  if (mol%npbc > 0) then
!     call dgemv('n', 9, len(mol), 1.0_wp, dcndL, 9, dEdcn, 1, 1.0_wp, sigma, 1)
!  endif
!  ```
module d3mod_coordination_number
   use iso_fortran_env, only: wp => real64
   implicit none
   public :: get_coordination_number
   public :: cut_coordination_number
   public :: d3_cn_type
   private

   !> Implemented kinds of counting functions for the CN.
   type :: enum_cn_type
      integer :: exp = 1
      integer :: erf = 2
   end type
   !> Enumerator for the coordination number types.
   type(enum_cn_type), parameter :: d3_cn_type = enum_cn_type()

   abstract interface
      real(wp) pure function counting_function(k, r, r0)
         import wp
         real(wp), intent(in) :: k, r, r0
      end function counting_function
   end interface

contains

!> Geometric fractional coordination number, supports both error function
!  and exponential counting functions.
subroutine get_coordination_number(mol, neighlist, cf, cn, dcndr, dcndL)
   use d3def_molecule
   use d3def_neighbourlist
   !> Molecular structure information.
   type(d3_molecule), intent(in) :: mol
   !> Neighbourlist.
   type(d3_neighbourlist), intent(in) :: neighlist
   !> Coordination number type (by counting function).
   integer, intent(in) :: cf
   !> Error function coordination number.
   real(wp), intent(out) :: cn(:)
   !> Derivative of the CN with respect to the Cartesian coordinates.
   real(wp), intent(out) :: dcndr(:, :, :)
   !> Derivative of the CN with respect to strain deformations.
   real(wp), intent(out) :: dcndL(:, :, :)

   real(wp), parameter :: kcn_exp = 16.0_wp
   real(wp), parameter :: kcn_erf = 7.5_wp

   select case(cf)
   case(d3_cn_type%exp)
      call ncoord_impl(mol, neighlist, kcn_exp, exp_count, dexp_count, &
         &             cn, dcndr, dcndL)
   case(d3_cn_type%erf)
      call ncoord_impl(mol, neighlist, kcn_erf, erf_count, derf_count, &
         &             cn, dcndr, dcndL)
   end select

end subroutine get_coordination_number

!> Actual implementation of the coordination number, takes a generic counting
!  function to return the respective CN.
subroutine ncoord_impl(mol, neighlist, kcn, cfunc, dfunc, cn, dcndr, dcndL)
   use d3def_molecule
   use d3def_neighbourlist
   use d3par_covalent_radii, only: covalent_radius => covalent_radius_d3
   !> Molecular structure information.
   type(d3_molecule), intent(in) :: mol
   !> Neighbourlist.
   type(d3_neighbourlist), target, intent(in) :: neighlist
   !> Function implementing the counting function.
   procedure(counting_function) :: cfunc
   !> Function implementing the derivative of counting function w.r.t. distance.
   procedure(counting_function) :: dfunc
   !> Steepness of counting function
   real(wp), intent(in) :: kcn
   !> Error function coordination number.
   real(wp), intent(out) :: cn(:)
   !> Derivative of the CN with respect to the Cartesian coordinates.
   real(wp), intent(out) :: dcndr(:, :, :)
   !> Derivative of the CN with respect to strain deformations.
   real(wp), intent(out) :: dcndL(:, :, :)

   type(d3_neighlist_iterator) :: neighiter

   integer :: neighs
   integer :: iat, jat, ati, atj, ij
   integer :: image(iter_chunk_size)
   real(wp) :: r2, r1, rc, rij(3), countf, countd(3), stress(3, 3)
   real(wp) :: dists2(iter_chunk_size), coords(3, iter_chunk_size)

   cn = 0.0_wp
   dcndr = 0.0_wp
   dcndL = 0.0_wp

   !$omp parallel do default(none) schedule(runtime) &
   !$omp reduction(+:cn, dcndr, dcndL) shared(mol, neighlist, kcn) &
   !$omp private(neighiter, neighs, ij, jat, ati, atj, coords, image, dists2, &
   !$omp&        r2, rij, r1, rc, countf, countd)
   do iat = 1, len(mol)
      ati = mol%at(iat)
      call neighlist%get_iterator(neighiter, iat)
      neighs = iter_chunk_size
      do while(neighs == iter_chunk_size)
         call neighiter%next(neighs, coords=coords, image=image, dists2=dists2)
         do ij = 1, neighs
            jat = image(ij)
            r2 = dists2(ij)
            rij = mol%xyz(:, iat) - coords(:, ij)
            atj = mol%at(jat)
            r1 = sqrt(r2)

            rc = covalent_radius(ati) + covalent_radius(atj)

            countf = cfunc(kcn, r1, rc)
            countd = dfunc(kcn, r1, rc) * rij/r1

            cn(iat) = cn(iat) + countf
            if (iat /= jat) then
               cn(jat) = cn(jat) + countf
            endif

            dcndr(:, iat, iat) = dcndr(:, iat, iat) + countd
            dcndr(:, jat, jat) = dcndr(:, jat, jat) - countd
            dcndr(:, iat, jat) = dcndr(:, iat, jat) + countd
            dcndr(:, jat, iat) = dcndr(:, jat, iat) - countd

            stress = spread(countd, 1, 3) * spread(rij, 2, 3)

            dcndL(:, :, iat) = dcndL(:, :, iat) + stress
            if (iat /= jat) then
               dcndL(:, :, jat) = dcndL(:, :, jat) + stress
            endif

         enddo
      enddo
   enddo
   !$omp end parallel do
end subroutine ncoord_impl

!> cutoff function for large coordination numbers
subroutine cut_coordination_number(cn, dcndr, dcndL, cn_max)
   !> on input coordination number, on output modified CN
   real(wp), intent(inout) :: cn(:)
   !> on input derivative of CN w.r.t. cartesian coordinates,
   !  on output derivative of modified CN
   real(wp), intent(inout), optional :: dcndr(:, :, :)
   !> on input derivative of CN w.r.t. strain deformation,
   !  on output derivative of modified CN
   real(wp), intent(inout), optional :: dcndL(:, :, :)
   !> maximum CN (not strictly obeyed)
   real(wp), intent(in), optional :: cn_max

   real(wp) :: cnmax
   integer :: iat

   if (present(cn_max)) then
      cnmax = cn_max
   else
      cnmax = 4.5_wp
   endif

   if (cnmax <= 0.0_wp) return

   if (present(dcndL)) then
      do iat = 1, size(dcndL, 3)
         dcndL(:, :, iat) = dcndL(:, :, iat) * dcut_cn(cn(iat), cnmax)
      enddo
   endif

   if (present(dcndr)) then
      do iat = 1, size(dcndr, 3)
         dcndr(:, :, iat) = dcndr(:, :, iat) * dcut_cn(cn(iat), cnmax)
      enddo
   endif

   do iat = 1, size(cn, 1)
      cn(iat) = cut_cn(cn(iat), cnmax)
   enddo

end subroutine cut_coordination_number

!> Exponential counting function for coordination number constributions.
pure function exp_count(k, r, r0) result(count)
   real(wp), intent(in) :: k
   real(wp), intent(in) :: r
   real(wp), intent(in) :: r0
   real(wp) :: count
   count = 1.0_wp/(1.0_wp+exp(-k*(r0/r-1.0_wp)))
end function exp_count

!> Derivative of the counting function w.r.t. the distance.
pure function dexp_count(k, r, r0) result(count)
   real(wp), intent(in) :: k
   real(wp), intent(in) :: r
   real(wp), intent(in) :: r0
   real(wp) :: count
   real(wp) :: expterm
   expterm=exp(-k*(r0/r-1._wp))
   count = -k*r0*expterm/(r**2*((expterm+1._wp)**2))
end function dexp_count

!> Error function counting function for coordination number constributions.
pure function erf_count(k, r, r0) result(count)
   real(wp), intent(in) :: k
   real(wp), intent(in) :: r
   real(wp), intent(in) :: r0
   real(wp) :: count
   count = 0.5_wp * (1.0_wp + erf(-k*(r-r0)/r0))
end function erf_count

!> Derivative of the counting function w.r.t. the distance.
pure function derf_count(k, r, r0) result(count)
   use d3par_constants, only: pi
   real(wp), parameter :: sqrtpi = sqrt(pi)
   real(wp), intent(in) :: k
   real(wp), intent(in) :: r
   real(wp), intent(in) :: r0
   real(wp) :: count
   count = -k/sqrtpi/r0*exp(-k**2*(r-r0)**2/r0**2)
end function derf_count

!> Cutting function for the coordination number.
pure elemental function cut_cn(cn, cut) result(cnp)
   real(wp), intent(in) :: cn
   real(wp), intent(in) :: cut
   real(wp) :: cnp
   cnp = log(1.0_wp + exp(cut)) - log(1.0_wp + exp(cut - cn))
end function cut_cn

!> Derivative of the cutting function w.r.t. coordination number
pure elemental function dcut_cn(cn, cut) result(dcnpdcn)
   real(wp), intent(in) :: cn
   real(wp), intent(in) :: cut
   real(wp) :: dcnpdcn
   dcnpdcn = exp(cut)/(exp(cut) + exp(cn))
end function dcut_cn

end module d3mod_coordination_number
