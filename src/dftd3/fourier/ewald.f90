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

!> Reciprocal space evaluation of the two-body dispersion energy.
!>
!> With a separable representation of the C6 coefficients the lattice sum
!>
!>    E = -1/2 sum_l lambda_l sum_T sum_ij C_li C_lj phi(|r_ij + T|)
!>
!> becomes a product of structure factors in reciprocal space. The damped pair
!> potential is bounded at the origin, so no real space complement is required and
!> the reciprocal sum alone converges exponentially. The self interaction of the
!> unrestricted double sum is removed by the value of the potential at the origin.
module dftd3_fourier_ewald
   use dftd3_cutoff, only : get_lattice_points
   use dftd3_fourier_decomposition, only : d3_lowrank_c6
   use dftd3_fourier_kernel, only : fourier_term, get_fourier_transform, &
      & get_potential_zero
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use mctc_io_constants, only : pi
   implicit none
   private

   public :: get_dispersion_ewald, get_reciprocal_lattice


contains


!> Evaluate the two-body dispersion energy by summation over the reciprocal lattice
subroutine get_dispersion_ewald(mol, lowrank, ghost, terms, nterm, kcut, gwvec, &
      & gwdcn, energies, dEdcn, gradient, sigma)

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Separable representation of the C6 coefficients
   class(d3_lowrank_c6), intent(in) :: lowrank

   !> Atoms excluded from the dispersion calculation
   logical, intent(in) :: ghost(:)

   !> Terms of the damped pair potential for each pair of species
   type(fourier_term), intent(in) :: terms(:, :, :)

   !> Number of terms for each pair of species
   integer, intent(in) :: nterm(:, :)

   !> Reciprocal space cutoff
   real(wp), intent(in) :: kcut

   !> Weighting function for the atomic reference systems
   real(wp), intent(in) :: gwvec(:, :)

   !> Derivative of the weighting function w.r.t. the coordination number
   real(wp), intent(in), optional :: gwdcn(:, :)

   !> Dispersion energy
   real(wp), intent(inout) :: energies(:)

   !> Derivative of the energy w.r.t. the coordination number
   real(wp), intent(inout), optional :: dEdcn(:)

   !> Dispersion gradient
   real(wp), intent(inout), optional :: gradient(:, :)

   !> Dispersion virial
   real(wp), intent(inout), optional :: sigma(:, :)

   logical :: grad
   integer :: nat, nid, rank, nk, ik, iat, izp, jzp, il, it, ic
   real(wp) :: vol, rec(3, 3), kvec(3), knorm, kr, erecip, phi, dphi, pval, dpval
   real(wp) :: qq, zre, zim
   complex(wp) :: zval
   real(wp), allocatable :: kpoints(:, :), c6l(:, :), dc6ldcn(:, :)
   real(wp), allocatable :: phihat(:, :), dphihat(:, :)
   real(wp), allocatable :: esum(:), dsum(:), gsum(:, :), ssum(:, :)
   real(wp), allocatable :: energies_local(:), dEdcn_local(:)
   real(wp), allocatable :: gradient_local(:, :), sigma_local(:, :)
   complex(wp), allocatable :: phase(:), sf(:, :), fvec(:, :)

   grad = present(gwdcn) .and. present(dEdcn) .and. present(gradient) &
      & .and. present(sigma)

   nat = mol%nat
   nid = mol%nid
   rank = lowrank%rank

   call get_reciprocal_lattice(mol%lattice, rec, vol)
   call get_lattice_points([.true., .true., .true.], rec, kcut, kpoints)
   nk = size(kpoints, 2)

   allocate(c6l(rank, nat))
   allocate(dc6ldcn(rank, nat), source=0.0_wp)
   if (grad) then
      call lowrank%get_weights(mol, ghost, gwvec, c6l, gwdcn, dc6ldcn)
   else
      call lowrank%get_weights(mol, ghost, gwvec, c6l)
   end if

   allocate(esum(nat), dsum(nat), source=0.0_wp)
   allocate(gsum(3, nat), source=0.0_wp)
   allocate(ssum(3, 3), source=0.0_wp)

   do iat = 1, nat
      izp = mol%id(iat)
      pval = 0.0_wp
      do it = 1, nterm(izp, izp)
         pval = pval + get_potential_zero(terms(it, izp, izp))
      end do
      do il = 1, rank
         esum(iat) = esum(iat) &
            & + 0.5_wp * lowrank%lambda(il) * c6l(il, iat)**2 * pval
         dsum(iat) = dsum(iat) &
            & + lowrank%lambda(il) * c6l(il, iat) * dc6ldcn(il, iat) * pval
      end do
   end do

   erecip = 0.0_wp

   !$omp parallel default(none) &
   !$omp shared(mol, lowrank, terms, nterm, kcut, kpoints, nk, nat, nid, rank, &
   !$omp& c6l, dc6ldcn, vol, grad, esum, dsum, gsum, ssum, erecip) &
   !$omp private(ik, iat, izp, jzp, il, it, ic, kvec, knorm, kr, phi, dphi, &
   !$omp& pval, dpval, qq, zval, zre, zim, phihat, dphihat, phase, sf, fvec, &
   !$omp& energies_local, dEdcn_local, gradient_local, sigma_local)
   allocate(phihat(nid, nid), dphihat(nid, nid))
   allocate(phase(nat), sf(rank, nid), fvec(rank, nid))
   allocate(energies_local(nat), source=0.0_wp)
   allocate(dEdcn_local(nat), source=0.0_wp)
   allocate(gradient_local(3, nat), source=0.0_wp)
   allocate(sigma_local(3, 3), source=0.0_wp)

   !$omp do schedule(runtime) reduction(+:erecip)
   do ik = 1, nk
      kvec(:) = kpoints(:, ik)
      knorm = norm2(kvec)
      if (knorm > kcut) cycle

      do izp = 1, nid
         do jzp = 1, izp
            phi = 0.0_wp
            dphi = 0.0_wp
            do it = 1, nterm(jzp, izp)
               call get_fourier_transform(terms(it, jzp, izp), knorm, pval, dpval)
               phi = phi + pval
               dphi = dphi + dpval
            end do
            phihat(jzp, izp) = phi
            phihat(izp, jzp) = phi
            dphihat(jzp, izp) = dphi
            dphihat(izp, jzp) = dphi
         end do
      end do

      sf(:, :) = cmplx(0.0_wp, 0.0_wp, wp)
      do iat = 1, nat
         kr = kvec(1)*mol%xyz(1, iat) + kvec(2)*mol%xyz(2, iat) &
            & + kvec(3)*mol%xyz(3, iat)
         phase(iat) = cmplx(cos(kr), sin(kr), wp)
         izp = mol%id(iat)
         sf(:, izp) = sf(:, izp) + c6l(:, iat) * phase(iat)
      end do

      fvec(:, :) = cmplx(0.0_wp, 0.0_wp, wp)
      do izp = 1, nid
         do jzp = 1, nid
            fvec(:, izp) = fvec(:, izp) + phihat(jzp, izp) * sf(:, jzp)
         end do
      end do

      do iat = 1, nat
         izp = mol%id(iat)
         do il = 1, rank
            zval = phase(iat) * conjg(fvec(il, izp))
            zre = real(zval, wp)
            energies_local(iat) = energies_local(iat) &
               & - 0.5_wp * lowrank%lambda(il) * c6l(il, iat) * zre / vol
            erecip = erecip - 0.5_wp * lowrank%lambda(il) * c6l(il, iat) * zre / vol
            if (grad) then
               zim = aimag(zval)
               dEdcn_local(iat) = dEdcn_local(iat) &
                  & - lowrank%lambda(il) * dc6ldcn(il, iat) * zre / vol
               gradient_local(:, iat) = gradient_local(:, iat) &
                  & + lowrank%lambda(il) * c6l(il, iat) * zim * kvec / vol
            end if
         end do
      end do

      if (grad .and. knorm > 0.0_wp) then
         qq = 0.0_wp
         do izp = 1, nid
            do jzp = 1, nid
               do il = 1, rank
                  qq = qq + lowrank%lambda(il) * dphihat(jzp, izp) &
                     & * real(sf(il, izp) * conjg(sf(il, jzp)), wp)
               end do
            end do
         end do
         do ic = 1, 3
            sigma_local(ic, :) = sigma_local(ic, :) &
               & + 0.5_wp * qq * kvec(ic) * kvec(:) / (knorm * vol)
         end do
      end if
   end do
   !$omp end do

   !$omp critical (get_dispersion_ewald_)
   esum(:) = esum(:) + energies_local(:)
   dsum(:) = dsum(:) + dEdcn_local(:)
   gsum(:, :) = gsum(:, :) + gradient_local(:, :)
   ssum(:, :) = ssum(:, :) + sigma_local(:, :)
   !$omp end critical (get_dispersion_ewald_)
   !$omp end parallel

   energies(:) = energies(:) + esum(:)
   if (grad) then
      do ic = 1, 3
         ssum(ic, ic) = ssum(ic, ic) - erecip
      end do
      dEdcn(:) = dEdcn(:) + dsum(:)
      gradient(:, :) = gradient(:, :) + gsum(:, :)
      sigma(:, :) = sigma(:, :) + ssum(:, :)
   end if

end subroutine get_dispersion_ewald


!> Reciprocal lattice vectors and cell volume of a periodic structure
pure subroutine get_reciprocal_lattice(lattice, rec, vol)

   !> Direct lattice vectors
   real(wp), intent(in) :: lattice(:, :)

   !> Reciprocal lattice vectors
   real(wp), intent(out) :: rec(3, 3)

   !> Volume of the unit cell
   real(wp), intent(out) :: vol

   real(wp) :: cross(3, 3), det

   call crossproduct(lattice(:, 2), lattice(:, 3), cross(:, 1))
   call crossproduct(lattice(:, 3), lattice(:, 1), cross(:, 2))
   call crossproduct(lattice(:, 1), lattice(:, 2), cross(:, 3))

   det = lattice(1, 1)*cross(1, 1) + lattice(2, 1)*cross(2, 1) &
      & + lattice(3, 1)*cross(3, 1)

   vol = abs(det)
   rec(:, :) = 2.0_wp * pi * cross / det

end subroutine get_reciprocal_lattice


pure subroutine crossproduct(a, b, c)
   real(wp), intent(in) :: a(3), b(3)
   real(wp), intent(out) :: c(3)

   c(1) = a(2)*b(3) - b(2)*a(3)
   c(2) = a(3)*b(1) - b(3)*a(1)
   c(3) = a(1)*b(2) - b(1)*a(2)

end subroutine crossproduct


end module dftd3_fourier_ewald
