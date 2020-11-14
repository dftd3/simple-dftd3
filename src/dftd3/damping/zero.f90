! This file is part of s-dftd3.
! SPDX-Identifier: LGLP-3.0-or-later
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

module dftd3_damping_zero
   use dftd3_damping, only : damping_param
   use dftd3_data, only : get_r4r2_val, get_vdw_rad
   use dftd3_param, only : d3_param
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   implicit none

   public :: zero_damping_param, new_zero_damping


   !> Zero (Chai-Head-Gordon) damping model
   type, extends(damping_param) :: zero_damping_param
      real(wp) :: s6
      real(wp) :: s8
      real(wp) :: s9
      real(wp) :: rs6
      real(wp) :: rs8
      real(wp) :: alp
      real(wp), allocatable :: r4r2(:)
      real(wp), allocatable :: rvdw(:, :)
   contains

      !> Evaluate pairwise dispersion energy expression
      procedure :: get_dispersion2

      !> Evaluate ATM three-body dispersion energy expression
      procedure :: get_dispersion3

   end type zero_damping_param


contains


!> Create new zero damping model
subroutine new_zero_damping(self, param, num)

   !> Zero damping parameters
   type(zero_damping_param), intent(out) :: self

   !> Parameters
   type(d3_param), intent(in) :: param

   !> Atomic numbers
   integer, intent(in) :: num(:)

   integer :: isp, jsp, izp, jzp

   self%s6 = param%s6   
   self%s8 = param%s8   
   self%s9 = param%s9   
   self%rs6 = param%rs6  
   self%rs8 = param%rs8  
   self%alp = param%alp  

   allocate(self%r4r2(size(num)))
   do isp = 1, size(num)
      izp = num(isp)
      self%r4r2(isp) = get_r4r2_val(izp)
   end do

   allocate(self%rvdw(size(num), size(num)))
   do isp = 1, size(num)
      izp = num(isp)
      do jsp = 1, isp
         jzp = num(jsp)
         self%rvdw(jsp, isp) = get_vdw_rad(izp, jzp)
         self%rvdw(isp, jsp) = self%rvdw(jsp, isp)
      end do
   end do

end subroutine new_zero_damping


!> Evaluation of the dispersion energy expression
subroutine get_dispersion2(self, mol, trans, cutoff, c6, dc6dcn, &
      & energy, dEdcn, gradient, sigma)

   !> Damping parameters
   class(zero_damping_param), intent(in) :: self

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Lattice points
   real(wp), intent(in) :: trans(:, :)

   !> Real space cutoff
   real(wp), intent(in) :: cutoff

   !> C6 coefficients for all atom pairs.
   real(wp), intent(in) :: c6(:, :)

   !> Derivative of the C6 w.r.t. the coordination number
   real(wp), intent(in), optional :: dc6dcn(:, :)

   !> Dispersion energy
   real(wp), intent(inout) :: energy

   !> Derivative of the energy w.r.t. the coordination number
   real(wp), intent(inout), optional :: dEdcn(:)

   !> Dispersion gradient
   real(wp), intent(inout), optional :: gradient(:, :)

   !> Dispersion virial
   real(wp), intent(inout), optional :: sigma(:, :)

   logical :: grad
   integer :: iat, jat, izp, jzp, jtr
   real(wp) :: vec(3), r2, r1, r6, r8, t6, t8, d6, d8, f6, f8, alp6, alp8
   real(wp) :: edisp, gdisp, cutoff2, r0ij, rrij, c6ij, dE, dG(3), dS(3, 3)

   if (abs(self%s6) < epsilon(1.0_wp) .and. abs(self%s8) < epsilon(1.0_wp)) return
   grad = present(dc6dcn) .and. present(dEdcn) .and. present(gradient) &
      & .and. present(sigma)
   cutoff2 = cutoff*cutoff
   alp6 = self%alp
   alp8 = self%alp + 2.0_wp

   if (grad) then
      do iat = 1, mol%nat
         izp = mol%id(iat)
         do jat = 1, iat
            jzp = mol%id(jat)
            rrij = 3*self%r4r2(izp)*self%r4r2(jzp)
            r0ij = self%rvdw(jzp, izp)
            c6ij = c6(jat, iat)
            do jtr = 1, size(trans, 2)
               vec(:) = mol%xyz(:, iat) - (mol%xyz(:, jat) + trans(:, jtr))
               r2 = vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3)
               if (r2 > cutoff2 .or. r2 < epsilon(1.0_wp)) cycle
               r1 = sqrt(r2)

               r6 = r2*r2*r2
               r8 = r6*r2

               t6 = (self%rs6*r0ij/r1)**alp6
               t8 = (self%rs8*r0ij/r1)**alp8

               f6 = 1.0_wp / (1.0_wp + 6.0_wp*t6)
               f8 = 1.0_wp / (1.0_wp + 6.0_wp*t8)

               d6 = -6.0_wp * f6 / r2 &
                  & + 6.0_wp * alp6 * t6 * f6**2 / r2
               d8 = -8.0_wp * f8 / r2 &
                  & + 6.0_wp * alp8 * t8 * f8**2 / r2

               edisp = self%s6 * f6 / r6 + self%s8 * rrij * f8 / r8
               gdisp = self%s6 * d6 / r6 + self%s8 * rrij * d8 / r8

               dE = -c6ij*edisp * 0.5_wp
               dG(:) = -c6ij*gdisp*vec
               dS(:, :) = spread(dG, 1, 3) * spread(vec, 2, 3) * 0.5_wp

               energy = energy + dE
               dEdcn(iat) = dEdcn(iat) - dc6dcn(iat, jat) * edisp
               sigma(:, :) = sigma + dS
               if (iat /= jat) then
                  energy = energy + dE
                  dEdcn(jat) = dEdcn(jat) - dc6dcn(jat, iat) * edisp
                  gradient(:, iat) = gradient(:, iat) + dG
                  gradient(:, jat) = gradient(:, jat) - dG
                  sigma(:, :) = sigma + dS
               end if
            end do
         end do
      end do
   else
      do iat = 1, mol%nat
         izp = mol%id(iat)
         do jat = 1, iat
            jzp = mol%id(jat)
            rrij = 3*self%r4r2(izp)*self%r4r2(jzp)
            r0ij = self%rvdw(jzp, izp)
            c6ij = c6(jat, iat)
            do jtr = 1, size(trans, 2)
               vec(:) = mol%xyz(:, iat) - (mol%xyz(:, jat) + trans(:, jtr))
               r2 = vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3)
               if (r2 > cutoff2 .or. r2 < epsilon(1.0_wp)) cycle
               r1 = sqrt(r2)

               r6 = r2*r2*r2
               r8 = r6*r2

               t6 = (self%rs6*r0ij/r1)**alp6
               t8 = (self%rs8*r0ij/r1)**alp8

               f6 = 1.0_wp / (1.0_wp + 6.0_wp*t6)
               f8 = 1.0_wp / (1.0_wp + 6.0_wp*t8)

               edisp = self%s6 * f6 / r6 + self%s8 * rrij * f8 / r8

               dE = -c6ij*edisp * 0.5_wp

               energy = energy + dE
               if (iat /= jat) then
                  energy = energy + dE
               end if
            end do
         end do
      end do
   end if
end subroutine get_dispersion2


!> Evaluation of the dispersion energy expression
subroutine get_dispersion3(self, mol, trans, cutoff, c6, dc6dcn, &
      & energy, dEdcn, gradient, sigma)

   !> Damping parameters
   class(zero_damping_param), intent(in) :: self

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Lattice points
   real(wp), intent(in) :: trans(:, :)

   !> Real space cutoff
   real(wp), intent(in) :: cutoff

   !> C6 coefficients for all atom pairs.
   real(wp), intent(in) :: c6(:, :)

   !> Derivative of the C6 w.r.t. the coordination number
   real(wp), intent(in), optional :: dc6dcn(:, :)

   !> Dispersion energy
   real(wp), intent(inout) :: energy

   !> Derivative of the energy w.r.t. the coordination number
   real(wp), intent(inout), optional :: dEdcn(:)

   !> Dispersion gradient
   real(wp), intent(inout), optional :: gradient(:, :)

   !> Dispersion virial
   real(wp), intent(inout), optional :: sigma(:, :)

   logical :: grad
   integer :: iat, jat, kat, izp, jzp, kzp, jtr, ktr
   real(wp) :: vij(3), vjk(3), vik(3), r2ij, r2jk, r2ik, c6ij, c6jk, c6ik, triple
   real(wp) :: r0ij, r0jk, r0ik, r0, r1, r2, r3, r5, rr, fdmp, dfdmp, ang, dang
   real(wp) :: cutoff2, alp, c9, dE, dGij(3), dGjk(3), dGik(3), dS(3, 3)

   if (abs(self%s9) < epsilon(1.0_wp)) return
   grad = present(dc6dcn) .and. present(dEdcn) .and. present(gradient) &
      & .and. present(sigma)
   cutoff2 = cutoff*cutoff
   alp = self%alp + 2.0_wp

   if (grad) then
      do iat = 1, mol%nat
         izp = mol%id(iat)
         do jat = 1, iat
            jzp = mol%id(jat)
            c6ij = c6(jat, iat)
            r0ij = 4.0_wp / 3.0_wp * self%rvdw(jzp, izp)
            do jtr = 1, size(trans, 2)
               vij(:) = mol%xyz(:, jat) + trans(:, jtr) - mol%xyz(:, iat)
               r2ij = vij(1)*vij(1) + vij(2)*vij(2) + vij(3)*vij(3)
               if (r2ij > cutoff2 .or. r2ij < epsilon(1.0_wp)) cycle
               do kat = 1, jat
                  kzp = mol%id(kat)
                  c6ik = c6(kat, iat)
                  c6jk = c6(kat, jat)
                  c9 = -self%s9 * sqrt(abs(c6ij*c6ik*c6jk))
                  r0ik = 4.0_wp / 3.0_wp * self%rvdw(kzp, izp)
                  r0jk = 4.0_wp / 3.0_wp * self%rvdw(kzp, jzp)
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
                     energy = energy - dE

                     gradient(:, iat) = gradient(:, iat) - dGij - dGik
                     gradient(:, jat) = gradient(:, jat) + dGij - dGjk
                     gradient(:, kat) = gradient(:, kat) + dGik + dGjk

                     dS(:,:) = spread(dGij, 1, 3) * spread(vij, 2, 3)&
                        & + spread(dGik, 1, 3) * spread(vik, 2, 3)&
                        & + spread(dGjk, 1, 3) * spread(vjk, 2, 3)

                     sigma(:,:) = sigma + dS * triple

                     dEdcn(iat) = dEdcn(iat) - dE * 0.5_wp &
                        & * (dc6dcn(iat, jat) / c6ij + dc6dcn(iat, kat) / c6ik)
                     dEdcn(jat) = dEdcn(jat) - dE * 0.5_wp &
                        & * (dc6dcn(jat, iat) / c6ij + dc6dcn(jat, kat) / c6jk)
                     dEdcn(kat) = dEdcn(kat) - dE * 0.5_wp &
                        & * (dc6dcn(kat, iat) / c6ik + dc6dcn(kat, jat) / c6jk)
                  end do
               end do
            end do
         end do
      end do
   else
      do iat = 1, mol%nat
         izp = mol%id(iat)
         do jat = 1, iat
            jzp = mol%id(jat)
            c6ij = c6(jat, iat)
            r0ij = 4.0_wp / 3.0_wp * self%rvdw(jzp, izp)
            do jtr = 1, size(trans, 2)
               vij(:) = mol%xyz(:, jat) + trans(:, jtr) - mol%xyz(:, iat)
               r2ij = vij(1)*vij(1) + vij(2)*vij(2) + vij(3)*vij(3)
               if (r2ij > cutoff2 .or. r2ij < epsilon(1.0_wp)) cycle
               do kat = 1, jat
                  kzp = mol%id(kat)
                  c6ik = c6(kat, iat)
                  c6jk = c6(kat, jat)
                  c9 = -self%s9 * sqrt(abs(c6ij*c6ik*c6jk))
                  r0ik = 4.0_wp / 3.0_wp * self%rvdw(kzp, izp)
                  r0jk = 4.0_wp / 3.0_wp * self%rvdw(kzp, jzp)
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

                     dE = rr * c9 * triple
                     energy = energy - dE
                  end do
               end do
            end do
         end do
      end do
   end if

end subroutine get_dispersion3


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


end module dftd3_damping_zero
