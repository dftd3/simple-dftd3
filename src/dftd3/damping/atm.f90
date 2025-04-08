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

module dftd3_damping_atm
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   implicit none

   public :: get_atm_dispersion, get_atm_pairwise_dispersion


contains


!> Evaluation of the dispersion energy expression
subroutine get_atm_dispersion(mol, trans, cutoff, s9, rs9, alp, rvdw, c6, dc6dcn, &
      & energy, dEdcn, gradient, sigma)

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Lattice points
   real(wp), intent(in) :: trans(:, :)

   !> Real space cutoff
   real(wp), intent(in) :: cutoff

   !> Scaling for dispersion coefficients
   real(wp), intent(in) :: s9

   !> Scaling for van-der-Waals radii in damping function
   real(wp), intent(in) :: rs9

   !> Exponent of zero damping function
   real(wp), intent(in) :: alp

   !> Van-der-Waals radii for all element pairs
   real(wp), intent(in) :: rvdw(:, :)

   !> C6 coefficients for all atom pairs.
   real(wp), intent(in) :: c6(:, :)

   !> Derivative of the C6 w.r.t. the coordination number
   real(wp), intent(in), optional :: dc6dcn(:, :)

   !> Dispersion energy
   real(wp), intent(inout) :: energy(:)

   !> Derivative of the energy w.r.t. the coordination number
   real(wp), intent(inout), optional :: dEdcn(:)

   !> Dispersion gradient
   real(wp), intent(inout), optional :: gradient(:, :)

   !> Dispersion virial
   real(wp), intent(inout), optional :: sigma(:, :)

   logical :: grad

   if (abs(s9) < epsilon(1.0_wp)) return
   grad = present(dc6dcn) .and. present(dEdcn) .and. present(gradient) &
      & .and. present(sigma)

   if (grad) then
      call get_atm_dispersion_derivs(mol, trans, cutoff, s9, rs9, alp, rvdw, c6, dc6dcn, &
      & energy, dEdcn, gradient, sigma)
   else
      call get_atm_dispersion_energy(mol, trans, cutoff, s9, rs9, alp, rvdw, c6, energy)
   end if

end subroutine get_atm_dispersion


!> Evaluation of the dispersion energy expression
subroutine get_atm_dispersion_energy(mol, trans, cutoff, s9, rs9, alp, rvdw, c6, energy)

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Lattice points
   real(wp), intent(in) :: trans(:, :)

   !> Real space cutoff
   real(wp), intent(in) :: cutoff

   !> Scaling for dispersion coefficients
   real(wp), intent(in) :: s9

   !> Scaling for van-der-Waals radii in damping function
   real(wp), intent(in) :: rs9

   !> Exponent of zero damping function
   real(wp), intent(in) :: alp

   !> Van-der-Waals radii for all element pairs
   real(wp), intent(in) :: rvdw(:, :)

   !> C6 coefficients for all atom pairs.
   real(wp), intent(in) :: c6(:, :)

   !> Dispersion energy
   real(wp), intent(inout) :: energy(:)

   integer :: iat, jat, kat, izp, jzp, kzp, jtr, ktr
   real(wp) :: vij(3), vjk(3), vik(3), r2ij, r2jk, r2ik, c6ij, c6jk, c6ik, triple
   real(wp) :: r0ij, r0jk, r0ik, r0, r1, r2, r3, r5, rr, fdmp, ang
   real(wp) :: cutoff2, c9, dE

   ! Thread-private arrays for reduction
   ! Set to 0 explicitly as the shared variants are potentially non-zero (inout)
   real(wp), allocatable :: energy_local(:)

   cutoff2 = cutoff*cutoff

   !$omp parallel default(none) &
   !$omp shared(mol, trans, c6, s9, rs9, alp, rvdw, cutoff2) &
   !$omp private(iat, jat, kat, izp, jzp, kzp, jtr, ktr, vij, vjk, vik, &
   !$omp& r2ij, r2jk, r2ik, c6ij, c6jk, c6ik, triple, r0ij, r0jk, r0ik, r0, &
   !$omp& r1, r2, r3, r5, rr, fdmp, ang, c9, dE) &
   !$omp shared(energy) &
   !$omp private(energy_local)
   allocate(energy_local(size(energy, 1)), source=0.0_wp)
   !$omp do schedule(runtime)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat
         jzp = mol%id(jat)
         c6ij = c6(jat, iat)
         r0ij = rs9 * rvdw(jzp, izp)
         do jtr = 1, size(trans, 2)
            vij(:) = mol%xyz(:, jat) + trans(:, jtr) - mol%xyz(:, iat)
            r2ij = vij(1)*vij(1) + vij(2)*vij(2) + vij(3)*vij(3)
            if (r2ij > cutoff2 .or. r2ij < epsilon(1.0_wp)) cycle
            do kat = 1, jat
               kzp = mol%id(kat)
               c6ik = c6(kat, iat)
               c6jk = c6(kat, jat)
               c9 = -s9 * sqrt(abs(c6ij*c6ik*c6jk))
               r0ik = rs9 * rvdw(kzp, izp)
               r0jk = rs9 * rvdw(kzp, jzp)
               r0 = r0ij * r0ik * r0jk
               triple = triple_scale(iat, jat, kat)
               do ktr = 1, size(trans, 2)
                  vik(:) = mol%xyz(:, kat) + trans(:, ktr) - mol%xyz(:, iat)
                  r2ik = vik(1)*vik(1) + vik(2)*vik(2) + vik(3)*vik(3)
                  if (r2ik > cutoff2 .or. r2ik < epsilon(1.0_wp)) cycle
                  vjk(:) = mol%xyz(:, kat) + trans(:, ktr) - mol%xyz(:, jat) &
                     & - trans(:, jtr)
                  r2jk = vjk(1)*vjk(1) + vjk(2)*vjk(2) + vjk(3)*vjk(3)
                  if (r2jk > cutoff2 .or. r2jk < epsilon(1.0_wp)) cycle
                  r2 = r2ij*r2ik*r2jk
                  r1 = sqrt(r2)
                  r3 = r2 * r1
                  r5 = r3 * r2

                  fdmp = 1.0_wp / (1.0_wp + 6.0_wp * (r0 / r1)**(alp / 3.0_wp))
                  ang = 0.375_wp*(r2ij + r2jk - r2ik)*(r2ij - r2jk + r2ik)&
                     & *(-r2ij + r2jk + r2ik) / r5 + 1.0_wp / r3

                  rr = ang*fdmp

                  dE = rr * c9 * triple / 3.0_wp
                  energy_local(iat) = energy_local(iat) - dE
                  energy_local(jat) = energy_local(jat) - dE
                  energy_local(kat) = energy_local(kat) - dE
               end do
            end do
         end do
      end do
   end do
   !$omp end do
   !$omp critical (get_atm_dispersion_energy_)
   energy(:) = energy(:) + energy_local(:)
   !$omp end critical (get_atm_dispersion_energy_)
   deallocate(energy_local)
   !$omp end parallel

end subroutine get_atm_dispersion_energy


!> Evaluation of the dispersion energy expression
subroutine get_atm_dispersion_derivs(mol, trans, cutoff, s9, rs9, alp, rvdw, c6, dc6dcn, &
      & energy, dEdcn, gradient, sigma)

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Lattice points
   real(wp), intent(in) :: trans(:, :)

   !> Real space cutoff
   real(wp), intent(in) :: cutoff

   !> Scaling for dispersion coefficients
   real(wp), intent(in) :: s9

   !> Scaling for van-der-Waals radii in damping function
   real(wp), intent(in) :: rs9

   !> Exponent of zero damping function
   real(wp), intent(in) :: alp

   !> Van-der-Waals radii for all element pairs
   real(wp), intent(in) :: rvdw(:, :)

   !> C6 coefficients for all atom pairs.
   real(wp), intent(in) :: c6(:, :)

   !> Derivative of the C6 w.r.t. the coordination number
   real(wp), intent(in) :: dc6dcn(:, :)

   !> Dispersion energy
   real(wp), intent(inout) :: energy(:)

   !> Derivative of the energy w.r.t. the coordination number
   real(wp), intent(inout) :: dEdcn(:)

   !> Dispersion gradient
   real(wp), intent(inout) :: gradient(:, :)

   !> Dispersion virial
   real(wp), intent(inout) :: sigma(:, :)

   integer :: iat, jat, kat, izp, jzp, kzp, jtr, ktr
   real(wp) :: vij(3), vjk(3), vik(3), r2ij, r2jk, r2ik, c6ij, c6jk, c6ik, triple
   real(wp) :: r0ij, r0jk, r0ik, r0, r1, r2, r3, r5, rr, fdmp, dfdmp, ang, dang
   real(wp) :: cutoff2, c9, dE, dGij(3), dGjk(3), dGik(3), dS(3, 3)

   ! Thread-private arrays for reduction
   ! Set to 0 explicitly as the shared variants are potentially non-zero (inout)
   real(wp), allocatable :: energy_local(:)
   real(wp), allocatable :: dEdcn_local(:)
   real(wp), allocatable :: gradient_local(:, :)
   real(wp), allocatable :: sigma_local(:, :)

   cutoff2 = cutoff*cutoff

   !$omp parallel default(none) &
   !$omp shared(mol, trans, c6, s9, rs9, alp, rvdw, cutoff2, dc6dcn) &
   !$omp private(iat, jat, kat, izp, jzp, kzp, jtr, ktr, vij, vjk, vik, &
   !$omp& r2ij, r2jk, r2ik, c6ij, c6jk, c6ik, triple, r0ij, r0jk, r0ik, r0, &
   !$omp& r1, r2, r3, r5, rr, fdmp, dfdmp, ang, dang, c9, dE, dGij, dGjk, &
   !$omp& dGik, dS) &
   !$omp shared(energy, gradient, sigma, dEdcn) &
   !$omp private(energy_local, gradient_local, sigma_local, dEdcn_local)
   allocate(energy_local(size(energy, 1)), source=0.0_wp)
   allocate(dEdcn_local(size(dEdcn, 1)), source=0.0_wp)
   allocate(gradient_local(size(gradient, 1), size(gradient, 2)), source=0.0_wp)
   allocate(sigma_local(size(sigma, 1), size(sigma, 2)), source=0.0_wp)
   !$omp do schedule(runtime)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat
         jzp = mol%id(jat)
         c6ij = c6(jat, iat)
         r0ij = rs9 * rvdw(jzp, izp)
         do jtr = 1, size(trans, 2)
            vij(:) = mol%xyz(:, jat) + trans(:, jtr) - mol%xyz(:, iat)
            r2ij = vij(1)*vij(1) + vij(2)*vij(2) + vij(3)*vij(3)
            if (r2ij > cutoff2 .or. r2ij < epsilon(1.0_wp)) cycle
            do kat = 1, jat
               kzp = mol%id(kat)
               c6ik = c6(kat, iat)
               c6jk = c6(kat, jat)
               c9 = -s9 * sqrt(abs(c6ij*c6ik*c6jk))
               r0ik = rs9 * rvdw(kzp, izp)
               r0jk = rs9 * rvdw(kzp, jzp)
               r0 = r0ij * r0ik * r0jk
               triple = triple_scale(iat, jat, kat)
               do ktr = 1, size(trans, 2)
                  vik(:) = mol%xyz(:, kat) + trans(:, ktr) - mol%xyz(:, iat)
                  r2ik = vik(1)*vik(1) + vik(2)*vik(2) + vik(3)*vik(3)
                  if (r2ik > cutoff2 .or. r2ik < epsilon(1.0_wp)) cycle
                  vjk(:) = mol%xyz(:, kat) + trans(:, ktr) - mol%xyz(:, jat) &
                     & - trans(:, jtr)
                  r2jk = vjk(1)*vjk(1) + vjk(2)*vjk(2) + vjk(3)*vjk(3)
                  if (r2jk > cutoff2 .or. r2jk < epsilon(1.0_wp)) cycle
                  r2 = r2ij*r2ik*r2jk
                  r1 = sqrt(r2)
                  r3 = r2 * r1
                  r5 = r3 * r2

                  fdmp = 1.0_wp / (1.0_wp + 6.0_wp * (r0 / r1)**(alp / 3.0_wp))
                  ang = 0.375_wp*(r2ij + r2jk - r2ik)*(r2ij - r2jk + r2ik)&
                     & *(-r2ij + r2jk + r2ik) / r5 + 1.0_wp / r3

                  rr = ang*fdmp

                  dfdmp = -2.0_wp * alp * (r0 / r1)**(alp / 3.0_wp) * fdmp**2

                  ! d/drij
                  dang = -0.375_wp * (r2ij**3 + r2ij**2 * (r2jk + r2ik)&
                     & + r2ij * (3.0_wp * r2jk**2 + 2.0_wp * r2jk*r2ik&
                     & + 3.0_wp * r2ik**2)&
                     & - 5.0_wp * (r2jk - r2ik)**2 * (r2jk + r2ik)) / r5
                  dGij(:) = c9 * (-dang*fdmp + ang*dfdmp) / r2ij * vij

                  ! d/drik
                  dang = -0.375_wp * (r2ik**3 + r2ik**2 * (r2jk + r2ij)&
                     & + r2ik * (3.0_wp * r2jk**2 + 2.0_wp * r2jk * r2ij&
                     & + 3.0_wp * r2ij**2)&
                     & - 5.0_wp * (r2jk - r2ij)**2 * (r2jk + r2ij)) / r5
                  dGik(:) = c9 * (-dang * fdmp + ang * dfdmp) / r2ik * vik

                  ! d/drjk
                  dang = -0.375_wp * (r2jk**3 + r2jk**2*(r2ik + r2ij)&
                     & + r2jk * (3.0_wp * r2ik**2 + 2.0_wp * r2ik * r2ij&
                     & + 3.0_wp * r2ij**2)&
                     & - 5.0_wp * (r2ik - r2ij)**2 * (r2ik + r2ij)) / r5
                  dGjk(:) = c9 * (-dang * fdmp + ang * dfdmp) / r2jk * vjk

                  dE = rr * c9 * triple
                  energy_local(iat) = energy_local(iat) - dE/3.0_wp
                  energy_local(jat) = energy_local(jat) - dE/3.0_wp
                  energy_local(kat) = energy_local(kat) - dE/3.0_wp

                  gradient_local(:, iat) = gradient_local(:, iat) - dGij - dGik
                  gradient_local(:, jat) = gradient_local(:, jat) + dGij - dGjk
                  gradient_local(:, kat) = gradient_local(:, kat) + dGik + dGjk

                  dS(:, :) = spread(dGij, 1, 3) * spread(vij, 2, 3)&
                     & + spread(dGik, 1, 3) * spread(vik, 2, 3)&
                     & + spread(dGjk, 1, 3) * spread(vjk, 2, 3)

                  sigma_local(:, :) = sigma_local + dS * triple

                  dEdcn_local(iat) = dEdcn_local(iat) - dE * 0.5_wp &
                     & * (dc6dcn(iat, jat) / c6ij + dc6dcn(iat, kat) / c6ik)
                  dEdcn_local(jat) = dEdcn_local(jat) - dE * 0.5_wp &
                     & * (dc6dcn(jat, iat) / c6ij + dc6dcn(jat, kat) / c6jk)
                  dEdcn_local(kat) = dEdcn_local(kat) - dE * 0.5_wp &
                     & * (dc6dcn(kat, iat) / c6ik + dc6dcn(kat, jat) / c6jk)
               end do
            end do
         end do
      end do
   end do
   !$omp end do
   !$omp critical (get_atm_dispersion_derivs_)
   energy(:) = energy(:) + energy_local(:)
   dEdcn(:) = dEdcn(:) + dEdcn_local(:)
   gradient(:, :) = gradient(:, :) + gradient_local(:, :)
   sigma(:, :) = sigma(:, :) + sigma_local(:, :)
   !$omp end critical (get_atm_dispersion_derivs_)
   deallocate(energy_local)
   deallocate(dEdcn_local)
   deallocate(gradient_local)
   deallocate(sigma_local)
   !$omp end parallel

end subroutine get_atm_dispersion_derivs


!> Evaluation of the dispersion energy expression
subroutine get_atm_pairwise_dispersion(mol, trans, cutoff, s9, rs9, alp, rvdw, c6, &
      & energy)

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Lattice points
   real(wp), intent(in) :: trans(:, :)

   !> Real space cutoff
   real(wp), intent(in) :: cutoff

   !> Scaling for dispersion coefficients
   real(wp), intent(in) :: s9

   !> Scaling for van-der-Waals radii in damping function
   real(wp), intent(in) :: rs9

   !> Exponent of zero damping function
   real(wp), intent(in) :: alp

   !> Van-der-Waals radii for all element pairs
   real(wp), intent(in) :: rvdw(:, :)

   !> C6 coefficients for all atom pairs.
   real(wp), intent(in) :: c6(:, :)

   !> Dispersion energy
   real(wp), intent(inout) :: energy(:, :)

   integer :: iat, jat, kat, izp, jzp, kzp, jtr, ktr
   real(wp) :: vij(3), vjk(3), vik(3), r2ij, r2jk, r2ik, c6ij, c6jk, c6ik, triple
   real(wp) :: r0ij, r0jk, r0ik, r0, r1, r2, r3, r5, rr, fdmp, ang
   real(wp) :: cutoff2, c9, dE

   ! Thread-private arrays for reduction
   ! Set to 0 explicitly as the shared variants are potentially non-zero (inout)
   real(wp), allocatable :: energy_local(:, :)

   if (abs(s9) < epsilon(1.0_wp)) return
   cutoff2 = cutoff*cutoff

   !$omp parallel default(none) &
   !$omp shared(mol, trans, c6, cutoff2, s9, rs9, alp, rvdw) &
   !$omp private(iat, jat, kat, izp, jzp, kzp, jtr, ktr, vij, vjk, vik, &
   !$omp& r2ij, r2jk, r2ik, c6ij, c6jk, c6ik, triple, r0ij, r0jk, r0ik, r0, &
   !$omp& r1, r2, r3, r5, rr, fdmp, ang, c9, dE) &
   !$omp shared(energy) &
   !$omp private(energy_local)
   allocate(energy_local(size(energy, 1), size(energy, 2)), source=0.0_wp)
   !$omp do schedule(runtime)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat
         jzp = mol%id(jat)
         c6ij = c6(jat, iat)
         r0ij = rs9 * rvdw(jzp, izp)
         do jtr = 1, size(trans, 2)
            vij(:) = mol%xyz(:, jat) + trans(:, jtr) - mol%xyz(:, iat)
            r2ij = vij(1)*vij(1) + vij(2)*vij(2) + vij(3)*vij(3)
            if (r2ij > cutoff2 .or. r2ij < epsilon(1.0_wp)) cycle
            do kat = 1, jat
               kzp = mol%id(kat)
               c6ik = c6(kat, iat)
               c6jk = c6(kat, jat)
               c9 = -s9 * sqrt(abs(c6ij*c6ik*c6jk))
               r0ik = rs9 * rvdw(kzp, izp)
               r0jk = rs9 * rvdw(kzp, jzp)
               r0 = r0ij * r0ik * r0jk
               triple = triple_scale(iat, jat, kat)
               do ktr = 1, size(trans, 2)
                  vik(:) = mol%xyz(:, kat) + trans(:, ktr) - mol%xyz(:, iat)
                  r2ik = vik(1)*vik(1) + vik(2)*vik(2) + vik(3)*vik(3)
                  if (r2ik > cutoff2 .or. r2ik < epsilon(1.0_wp)) cycle
                  vjk(:) = mol%xyz(:, kat) + trans(:, ktr) - mol%xyz(:, jat) &
                     & - trans(:, jtr)
                  r2jk = vjk(1)*vjk(1) + vjk(2)*vjk(2) + vjk(3)*vjk(3)
                  if (r2jk > cutoff2 .or. r2jk < epsilon(1.0_wp)) cycle
                  r2 = r2ij*r2ik*r2jk
                  r1 = sqrt(r2)
                  r3 = r2 * r1
                  r5 = r3 * r2

                  fdmp = 1.0_wp / (1.0_wp + 6.0_wp * (r0 / r1)**(alp / 3.0_wp))
                  ang = 0.375_wp*(r2ij + r2jk - r2ik)*(r2ij - r2jk + r2ik)&
                     & *(-r2ij + r2jk + r2ik) / r5 + 1.0_wp / r3

                  rr = ang*fdmp

                  dE = rr * c9 * triple / 6.0_wp
                  energy_local(jat, iat) = energy_local(jat, iat) - dE
                  energy_local(kat, iat) = energy_local(kat, iat) - dE
                  energy_local(iat, jat) = energy_local(iat, jat) - dE
                  energy_local(kat, jat) = energy_local(kat, jat) - dE
                  energy_local(iat, kat) = energy_local(iat, kat) - dE
                  energy_local(jat, kat) = energy_local(jat, kat) - dE
               end do
            end do
         end do
      end do
   end do
   !$omp end do
   !$omp critical (get_atm_pairwise_dispersion_)
   energy(:, :) = energy(:, :) + energy_local(:, :)
   !$omp end critical (get_atm_pairwise_dispersion_)
   deallocate(energy_local)
   !$omp end parallel

end subroutine get_atm_pairwise_dispersion


!> Logic exercise to distribute a triple energy to atomwise energies.
elemental function triple_scale(ii, jj, kk) result(triple)

   !> Atom indices
   integer, intent(in) :: ii, jj, kk

   !> Fraction of energy
   real(wp) :: triple

   if (ii == jj) then
      if (ii == kk) then
         ! ii'i" -> 1/6
         triple = 1.0_wp/6.0_wp
      else
         ! ii'j -> 1/2
         triple = 0.5_wp
      end if
   else
      if (ii /= kk .and. jj /= kk) then
         ! ijk -> 1 (full)
         triple = 1.0_wp
      else
         ! ijj' and iji' -> 1/2
         triple = 0.5_wp
      end if
   end if

end function triple_scale


end module dftd3_damping_atm
