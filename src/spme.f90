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

!> Smooth particle mesh evaluation of the two-body dispersion energy.
!>
!> The structure factors of the reciprocal space sum are interpolated from an
!> equispaced mesh using cardinal B-splines, which replaces the explicit sum over
!> atoms at every wave vector by a fast Fourier transform. The remaining sum runs
!> over the mesh instead of a sphere of wave vectors, so the cost grows as
!> O(N log N) with the number of atoms rather than O(N^2).
!>
!> Because the damping radius depends on the pair of species, the kernel cannot be
!> folded into a single mesh: the structure factors are resolved by species and the
!> transforms are repeated for every term of the low-rank expansion.
module dftd3_fourier_spme
   use dftd3_fourier_decomposition, only : d3_lowrank_c6
   use dftd3_fourier_fft, only : fft_mesh, new_fft_mesh, fft_3d, valid_mesh_size
   use dftd3_fourier_kernel, only : fourier_term, get_fourier_transform, &
      & get_potential_zero
   use dftd3_partition, only : work_partition, owns_index
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use mctc_io_constants, only : pi
   use mctc_io_math, only : matdet_3x3, matinv_3x3
   implicit none
   private

   public :: get_dispersion_spme, get_spme_mesh, spline_order


   !> Order of the cardinal B-spline interpolation.
   !>
   !> Six gives an interpolation error well below the mesh discretisation error for
   !> the meshes of interest, while keeping the spreading cost at 216 points per atom.
   integer, parameter :: spline_order = 6


contains


!> Number of mesh points resolving the requested reciprocal space cutoff
pure function get_spme_mesh(lattice, kcut) result(mesh)

   !> Direct lattice vectors
   real(wp), intent(in) :: lattice(:, :)

   !> Reciprocal space cutoff
   real(wp), intent(in) :: kcut

   !> Number of mesh points along each direction
   integer :: mesh(3)

   integer :: idir

   ! the mesh has to resolve oscillations of wavelength 2*pi/kcut along each
   ! lattice vector, the transform rounds the request up to a supported size
   do idir = 1, 3
      mesh(idir) = valid_mesh_size(ceiling(kcut * norm2(lattice(:, idir)) / pi))
   end do

end function get_spme_mesh


!> Evaluate the two-body dispersion energy on a particle mesh
subroutine get_dispersion_spme(mol, lowrank, ghost, terms, nterm, mesh, gwvec, &
      & gwdcn, energies, dEdcn, gradient, sigma, partition)

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

   !> Number of mesh points along each direction
   integer, intent(in) :: mesh(3)

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

   !> Work partition over the terms of the low-rank expansion
   type(work_partition), intent(in), optional :: partition

   logical :: grad
   type(fft_mesh) :: fmesh
   integer :: nat, nid, rank, iat, izp, jzp, il, it, ic, jc
   integer :: i1, i2, i3, j1, j2, j3, g1, g2, g3
   real(wp) :: vol, rec(3, 3), invlat(3, 3), kvec(3), knorm, erecip
   real(wp) :: pval, dpval, phi, dphi, wgt, dwgt(3), sfac, qq, escale
   real(wp) :: dudr(3, 3), tvec(3), dtvec(3)
   integer, allocatable :: base(:, :)
   real(wp), allocatable :: theta(:, :, :), dtheta(:, :, :)
   real(wp), allocatable :: bfac(:, :), c6l(:, :), dc6ldcn(:, :)
   real(wp), allocatable :: phik(:, :, :, :), dphik(:, :, :, :)
   real(wp), allocatable :: pot(:, :), esum(:), dsum(:), gsum(:, :), ssum(:, :)
   complex(wp), allocatable :: qgrid(:, :, :, :), tgrid(:, :, :)

   grad = present(gwdcn) .and. present(dEdcn) .and. present(gradient) &
      & .and. present(sigma)

   nat = mol%nat
   nid = mol%nid
   rank = lowrank%rank

   invlat(:, :) = matinv_3x3(mol%lattice)
   rec(:, :) = 2.0_wp * pi * transpose(invlat)
   vol = abs(matdet_3x3(mol%lattice))

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

   ! the unrestricted double sum contains the interaction of every atom with
   ! itself, which the value of the potential at the origin removes again
   do iat = 1, nat
      izp = mol%id(iat)
      pval = 0.0_wp
      do it = 1, nterm(izp, izp)
         pval = pval + get_potential_zero(terms(it, izp, izp))
      end do
      do il = 1, rank
         if (.not.owns_index(partition, il)) cycle
         esum(iat) = esum(iat) &
            & + 0.5_wp * lowrank%lambda(il) * c6l(il, iat)**2 * pval
         dsum(iat) = dsum(iat) &
            & + lowrank%lambda(il) * c6l(il, iat) * dc6ldcn(il, iat) * pval
      end do
   end do

   call new_fft_mesh(fmesh, mesh)
   call get_spline_weights(mol, mesh, base, theta, dtheta)
   call get_euler_factors(mesh, bfac)

   ! the transform of the pair potential depends on the species pair but not on
   ! the term of the expansion, so it is tabulated once for the whole mesh
   call get_mesh_kernel(mesh, rec, terms, nterm, nid, grad, phik, dphik)

   do ic = 1, 3
      do jc = 1, 3
         dudr(jc, ic) = mesh(ic) * rec(jc, ic) / (2.0_wp * pi)
      end do
   end do

   allocate(qgrid(mesh(1), mesh(2), mesh(3), nid))
   allocate(tgrid(mesh(1), mesh(2), mesh(3)))
   allocate(pot(nat, nid))

   erecip = 0.0_wp
   escale = -0.5_wp / vol

   ! the terms of the expansion decouple, so only one set of meshes is live and
   ! a partition can hand out whole terms
   do il = 1, rank
      if (.not.owns_index(partition, il)) cycle

      qgrid(:, :, :, :) = cmplx(0.0_wp, 0.0_wp, wp)
      do iat = 1, nat
         if (ghost(iat)) cycle
         izp = mol%id(iat)
         do j3 = 1, spline_order
            g3 = modulo(base(3, iat) + j3, mesh(3)) + 1
            do j2 = 1, spline_order
               g2 = modulo(base(2, iat) + j2, mesh(2)) + 1
               do j1 = 1, spline_order
                  g1 = modulo(base(1, iat) + j1, mesh(1)) + 1
                  wgt = theta(j1, 1, iat) * theta(j2, 2, iat) * theta(j3, 3, iat)
                  qgrid(g1, g2, g3, izp) = qgrid(g1, g2, g3, izp) &
                     & + cmplx(c6l(il, iat) * wgt, 0.0_wp, wp)
               end do
            end do
         end do
      end do

      do izp = 1, nid
         call fft_3d(fmesh, qgrid(:, :, :, izp), +1)
      end do

      ! contract the structure factors with the kernel to obtain the potential
      ! each species experiences, and accumulate energy and virial on the way
      do izp = 1, nid
         tgrid(:, :, :) = cmplx(0.0_wp, 0.0_wp, wp)
         !$omp parallel do schedule(runtime) default(none) collapse(2) &
         !$omp shared(mesh, nid, izp, phik, bfac, qgrid, tgrid) &
         !$omp private(i1, i2, i3, jzp, sfac)
         do i3 = 1, mesh(3)
            do i2 = 1, mesh(2)
               do i1 = 1, mesh(1)
                  sfac = bfac(i1, 1) * bfac(i2, 2) * bfac(i3, 3)
                  do jzp = 1, nid
                     tgrid(i1, i2, i3) = tgrid(i1, i2, i3) + sfac &
                        & * phik(pair_index(izp, jzp), i1, i2, i3) &
                        & * qgrid(i1, i2, i3, jzp)
                  end do
               end do
            end do
         end do

         call fft_3d(fmesh, tgrid, -1)

         ! gather the potential back onto the atoms of this species
         do iat = 1, nat
            pot(iat, izp) = 0.0_wp
         end do
         do iat = 1, nat
            if (ghost(iat) .or. mol%id(iat) /= izp) cycle
            wgt = 0.0_wp
            dwgt(:) = 0.0_wp
            do j3 = 1, spline_order
               g3 = modulo(base(3, iat) + j3, mesh(3)) + 1
               do j2 = 1, spline_order
                  g2 = modulo(base(2, iat) + j2, mesh(2)) + 1
                  do j1 = 1, spline_order
                     g1 = modulo(base(1, iat) + j1, mesh(1)) + 1
                     sfac = real(tgrid(g1, g2, g3), wp)
                     tvec(1) = theta(j1, 1, iat)
                     tvec(2) = theta(j2, 2, iat)
                     tvec(3) = theta(j3, 3, iat)
                     wgt = wgt + tvec(1) * tvec(2) * tvec(3) * sfac
                     if (grad) then
                        dtvec(1) = dtheta(j1, 1, iat) * tvec(2) * tvec(3)
                        dtvec(2) = tvec(1) * dtheta(j2, 2, iat) * tvec(3)
                        dtvec(3) = tvec(1) * tvec(2) * dtheta(j3, 3, iat)
                        dwgt(:) = dwgt(:) + matmul(dudr, dtvec) * sfac
                     end if
                  end do
               end do
            end do
            pot(iat, izp) = wgt
            esum(iat) = esum(iat) + escale * lowrank%lambda(il) * c6l(il, iat) * wgt
            erecip = erecip + escale * lowrank%lambda(il) * c6l(il, iat) * wgt
            if (grad) then
               dsum(iat) = dsum(iat) &
                  & + 2.0_wp * escale * lowrank%lambda(il) * dc6ldcn(il, iat) * wgt
               gsum(:, iat) = gsum(:, iat) &
                  & + 2.0_wp * escale * lowrank%lambda(il) * c6l(il, iat) * dwgt(:)
            end if
         end do
      end do

      if (grad) then
         !$omp parallel do schedule(runtime) default(none) collapse(2) &
         !$omp shared(mesh, nid, rec, bfac, dphik, qgrid, lowrank, il, vol) &
         !$omp private(i1, i2, i3, j1, j2, j3, izp, jzp, kvec, knorm, qq, sfac, ic) &
         !$omp reduction(+:ssum)
         do i3 = 1, mesh(3)
            do i2 = 1, mesh(2)
               do i1 = 1, mesh(1)
                  j1 = i1 - 1 - mesh(1) * ((2*(i1 - 1)) / mesh(1))
                  j2 = i2 - 1 - mesh(2) * ((2*(i2 - 1)) / mesh(2))
                  j3 = i3 - 1 - mesh(3) * ((2*(i3 - 1)) / mesh(3))
                  kvec(:) = rec(:, 1)*j1 + rec(:, 2)*j2 + rec(:, 3)*j3
                  knorm = norm2(kvec)
                  if (knorm <= 0.0_wp) cycle
                  sfac = bfac(i1, 1) * bfac(i2, 2) * bfac(i3, 3)
                  qq = 0.0_wp
                  do izp = 1, nid
                     do jzp = 1, nid
                        qq = qq + lowrank%lambda(il) &
                           & * dphik(pair_index(izp, jzp), i1, i2, i3) * sfac &
                           & * real(qgrid(i1, i2, i3, izp) &
                           &        * conjg(qgrid(i1, i2, i3, jzp)), wp)
                     end do
                  end do
                  do ic = 1, 3
                     ssum(ic, :) = ssum(ic, :) &
                        & + 0.5_wp * qq * kvec(ic) * kvec(:) / (knorm * vol)
                  end do
               end do
            end do
         end do
      end if
   end do

   energies(:) = energies(:) + esum(:)
   if (grad) then
      do ic = 1, 3
         ssum(ic, ic) = ssum(ic, ic) - erecip
      end do
      dEdcn(:) = dEdcn(:) + dsum(:)
      gradient(:, :) = gradient(:, :) + gsum(:, :)
      sigma(:, :) = sigma(:, :) + ssum(:, :)
   end if

end subroutine get_dispersion_spme


!> Packed index of a species pair
pure function pair_index(izp, jzp) result(ij)

   !> Species indices
   integer, intent(in) :: izp, jzp

   !> Index into the packed pair dimension
   integer :: ij

   ij = max(izp, jzp) * (max(izp, jzp) - 1) / 2 + min(izp, jzp)

end function pair_index


!> Tabulate the transform of the pair potential on the mesh
subroutine get_mesh_kernel(mesh, rec, terms, nterm, nid, grad, phik, dphik)

   !> Number of mesh points along each direction
   integer, intent(in) :: mesh(3)

   !> Reciprocal lattice vectors
   real(wp), intent(in) :: rec(3, 3)

   !> Terms of the damped pair potential for each pair of species
   type(fourier_term), intent(in) :: terms(:, :, :)

   !> Number of terms for each pair of species
   integer, intent(in) :: nterm(:, :)

   !> Number of species
   integer, intent(in) :: nid

   !> Whether the derivative is required
   logical, intent(in) :: grad

   !> Transform of the pair potential
   real(wp), allocatable, intent(out) :: phik(:, :, :, :)

   !> Derivative of the transform w.r.t. the wave number
   real(wp), allocatable, intent(out) :: dphik(:, :, :, :)

   integer :: i1, i2, i3, j1, j2, j3, izp, jzp, it, npair
   real(wp) :: kvec(3), knorm, phi, dphi, pval, dpval

   npair = nid * (nid + 1) / 2
   allocate(phik(npair, mesh(1), mesh(2), mesh(3)))
   if (grad) then
      allocate(dphik(npair, mesh(1), mesh(2), mesh(3)))
   else
      allocate(dphik(0, 0, 0, 0))
   end if

   !$omp parallel do schedule(runtime) default(none) collapse(2) &
   !$omp shared(mesh, rec, terms, nterm, nid, grad, phik, dphik) &
   !$omp private(i1, i2, i3, j1, j2, j3, izp, jzp, it, kvec, knorm, phi, dphi, &
   !$omp& pval, dpval)
   do i3 = 1, mesh(3)
      do i2 = 1, mesh(2)
         do i1 = 1, mesh(1)
            ! the mesh index wraps to the negative half of the Brillouin zone
            j1 = i1 - 1 - mesh(1) * ((2*(i1 - 1)) / mesh(1))
            j2 = i2 - 1 - mesh(2) * ((2*(i2 - 1)) / mesh(2))
            j3 = i3 - 1 - mesh(3) * ((2*(i3 - 1)) / mesh(3))
            kvec(:) = rec(:, 1)*j1 + rec(:, 2)*j2 + rec(:, 3)*j3
            knorm = norm2(kvec)
            do izp = 1, nid
               do jzp = 1, izp
                  phi = 0.0_wp
                  dphi = 0.0_wp
                  do it = 1, nterm(jzp, izp)
                     call get_fourier_transform(terms(it, jzp, izp), knorm, pval, dpval)
                     phi = phi + pval
                     dphi = dphi + dpval
                  end do
                  phik(pair_index(izp, jzp), i1, i2, i3) = phi
                  if (grad) dphik(pair_index(izp, jzp), i1, i2, i3) = dphi
               end do
            end do
         end do
      end do
   end do

end subroutine get_mesh_kernel


!> Cardinal B-spline weights and their derivatives for every atom
subroutine get_spline_weights(mol, mesh, base, theta, dtheta)

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Number of mesh points along each direction
   integer, intent(in) :: mesh(3)

   !> Index of the first mesh point contributing to an atom
   integer, allocatable, intent(out) :: base(:, :)

   !> Spline weights
   real(wp), allocatable, intent(out) :: theta(:, :, :)

   !> Derivative of the spline weights w.r.t. the mesh coordinate
   real(wp), allocatable, intent(out) :: dtheta(:, :, :)

   integer :: iat, idir
   real(wp) :: invlat(3, 3), frac(3), uval

   allocate(base(3, mol%nat))
   allocate(theta(spline_order, 3, mol%nat))
   allocate(dtheta(spline_order, 3, mol%nat))

   invlat(:, :) = matinv_3x3(mol%lattice)

   do iat = 1, mol%nat
      frac(:) = matmul(invlat, mol%xyz(:, iat))
      do idir = 1, 3
         uval = mesh(idir) * (frac(idir) - floor(frac(idir)))
         base(idir, iat) = floor(uval) - spline_order
         call fill_bspline(uval - floor(uval), theta(:, idir, iat), &
            & dtheta(:, idir, iat))
      end do
   end do

end subroutine get_spline_weights


!> Cardinal B-spline of order spline_order and its derivative at a mesh offset
pure subroutine fill_bspline(wval, theta, dtheta)

   !> Position within the mesh cell, between zero and one
   real(wp), intent(in) :: wval

   !> Spline weights
   real(wp), intent(out) :: theta(:)

   !> Derivative of the spline weights
   real(wp), intent(out) :: dtheta(:)

   integer :: korder, jorder
   real(wp) :: div

   theta(:) = 0.0_wp
   theta(1) = 1.0_wp - wval
   theta(2) = wval

   do korder = 3, spline_order - 1
      div = 1.0_wp / (korder - 1)
      theta(korder) = div * wval * theta(korder-1)
      do jorder = 1, korder - 2
         theta(korder-jorder) = div * ((wval + jorder) * theta(korder-jorder-1) &
            & + (korder - jorder - wval) * theta(korder-jorder))
      end do
      theta(1) = div * (1.0_wp - wval) * theta(1)
   end do

   ! the derivative follows from the splines of one order lower
   dtheta(1) = -theta(1)
   do jorder = 2, spline_order
      dtheta(jorder) = theta(jorder-1) - theta(jorder)
   end do

   div = 1.0_wp / (spline_order - 1)
   theta(spline_order) = div * wval * theta(spline_order-1)
   do jorder = 1, spline_order - 2
      theta(spline_order-jorder) = div * ((wval + jorder) * theta(spline_order-jorder-1) &
         & + (spline_order - jorder - wval) * theta(spline_order-jorder))
   end do
   theta(1) = div * (1.0_wp - wval) * theta(1)

end subroutine fill_bspline


!> Squared modulus of the Euler exponential spline factors along each direction
subroutine get_euler_factors(mesh, bfac)

   !> Number of mesh points along each direction
   integer, intent(in) :: mesh(3)

   !> Correction for the interpolation of the structure factors
   real(wp), allocatable, intent(out) :: bfac(:, :)

   integer :: idir, im, jm
   real(wp) :: arg, spl(spline_order), dspl(spline_order)
   complex(wp) :: denom

   call fill_bspline(0.0_wp, spl, dspl)

   allocate(bfac(maxval(mesh), 3), source=0.0_wp)
   do idir = 1, 3
      do im = 1, mesh(idir)
         denom = cmplx(0.0_wp, 0.0_wp, wp)
         do jm = 1, spline_order - 1
            arg = 2.0_wp * pi * (im - 1) * (jm - 1) / mesh(idir)
            denom = denom + spl(jm) * cmplx(cos(arg), sin(arg), wp)
         end do
         bfac(im, idir) = 1.0_wp / real(denom * conjg(denom), wp)
      end do
   end do

end subroutine get_euler_factors


end module dftd3_fourier_spme
