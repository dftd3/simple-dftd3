module d3mod_dftd3
   use iso_fortran_env, only: wp => real64
   implicit none
   public :: d3bj_eg
   private

contains 

!> Calculate the weights of the reference system and the derivatives w.r.t.
!  coordination number for later use.
subroutine weight_references(nat, numbers, weighting_factor, cn, gwvec, gwdcn)
   use d3par_dftd3
   !> Nr. of atoms (without periodic images)
   integer, intent(in) :: nat
   !> Atomic numbers of every atom.
   integer, intent(in) :: numbers(:)
   !> Exponent for the Gaussian weighting.
   real(wp), intent(in) :: weighting_factor
   !> Coordination number of every atom.
   real(wp), intent(in) :: cn(:)
   !> weighting for the atomic reference systems
   real(wp), intent(out) :: gwvec(:, :)
   !> derivative of the weighting function w.r.t. the coordination number
   real(wp), intent(out) :: gwdcn(:, :)

   integer :: iat, ati, iref, icount
   real(wp) :: norm, dnorm, wf, gw, expw, expd, gwk, dgwk

   gwVec = 0.0_wp
   gwVec = 0.0_wp

   do iat = 1, nat
      ati = numbers(iat)
      norm = 0.0_wp
      dnorm = 0.0_wp
      do iref = 1, number_of_references(ati)
         do icount = 1, reference_count(iref, ati)
            wf = icount * weighting_factor
            gw = weight_CN(wf, cn(iat), reference_cn(iref, ati))
            norm = norm + gw
            dnorm = dnorm + 2*wf*(reference_cn(iref, ati) - cn(iat)) * gw
         end do
      end do
      norm = 1.0_wp / norm
      do iref = 1, number_of_references(ati)
         expw = 0.0_wp
         expd = 0.0_wp
         do icount = 1, reference_count(iref, ati)
            wf = icount * weighting_factor
            gw = weight_CN(wf, cn(iat), reference_cn(iref, ati))
            expw = expw + gw
            expd = expd + 2*wf*(reference_cn(iref, ati) - cn(iat)) * gw
         enddo

         gwk = expw * norm
         if (gwk /= gwk) then
            if (maxval(reference_cn(:number_of_references(ati), ati)) &
               & == reference_cn(iref, ati)) then
               gwk = 1.0_wp
            else
               gwk = 0.0_wp
            endif
         endif
         gwVec(iref, iat) = gwk

         dgwk = expd*norm-expw*dnorm*norm**2
         if (dgwk /= dgwk) then
            dgwk = 0.0_wp
         endif
         gwdcn(iref, iat) = dgwk

      end do
   end do

end subroutine weight_references

!> calculate atomic dispersion coefficients and their derivatives w.r.t.
!  the coordination number.
subroutine get_atomic_c6(nat, numbers, gwvec, gwdcn, c6, dc6dcn)
   use d3par_dftd3
   !> Nr. of atoms (without periodic images)
   integer, intent(in) :: nat
   !> numbers of every atom.
   integer, intent(in) :: numbers(:)
   !> weighting function for the atomic reference systems
   real(wp), intent(in) :: gwvec(:, :)
   !> derivative of the weighting function w.r.t. the coordination number
   real(wp), intent(in) :: gwdcn(:, :)
   !> C6 coefficients for all atom pairs.
   real(wp), intent(out) :: c6(:, :)
   !> derivative of the C6 w.r.t. the coordination number
   real(wp), intent(out) :: dc6dcn(:, :)

   integer :: iat, jat, ati, atj, iref, jref
   real(wp) :: refc6, dc6, dc6dcn1, dc6dcn2

   c6 = 0.0_wp
   dc6dcn = 0.0_wp

   do iat = 1, nat
      ati = numbers(iat)
      do jat = 1, iat
         atj = numbers(jat)
         dc6 = 0.0_wp
         dc6dcn1 = 0.0_wp
         dc6dcn2 = 0.0_wp
         do iref = 1, number_of_references(ati)
            do jref = 1, number_of_references(atj)
               refc6 = get_c6(iref, jref, ati, atj)
               dc6 = dc6 + gwvec(iref, iat) * gwvec(jref, jat) * refc6
               dc6dcn1 = dc6dcn1 + gwdcn(iref, iat) * gwvec(jref, jat) * refc6
               dc6dcn2 = dc6dcn2 + gwvec(iref, iat) * gwdcn(jref, jat) * refc6
            end do
         end do
         c6(iat, jat) = dc6
         c6(jat, iat) = dc6
         dc6dcn(iat, jat) = dc6dcn1
         dc6dcn(jat, iat) = dc6dcn2
      end do
   end do
end subroutine get_atomic_c6

subroutine d3bj_eg(mol, neighlist, par, weighting_factor, &
      &            cn, dcndr, dcndL, energies, gradient, sigma)
   use d3def_molecule
   use d3def_neighbourlist
   use d3def_damping_parameters
   use d3par_dftd3
   use d3par_r4r2, only: r4r2 => sqrt_z_r4_over_r2

   class(d3_molecule), intent(in) :: mol
   class(d3_neighbourlist), intent(in) :: neighlist
   class(d3_damping_parameters), intent(in) :: par

   real(wp), intent(in) :: weighting_factor
   real(wp), intent(in) :: cn(:)
   real(wp), intent(in) :: dcndr(:, :, :)
   real(wp), intent(in) :: dcndL(:, :, :)

   real(wp), intent(inout) :: energies(:)
   real(wp), intent(inout) :: gradient(:, :)
   real(wp), intent(inout) :: sigma(:, :)

   type(d3_neighlist_iterator) :: neighiter

   integer :: nat, neighs, max_ref
   integer :: iat, jat, ati, atj, ij
   integer :: image(iter_chunk_size)

   real(wp) :: r4r2ij, r0, rij(3), r2, t6, t8, t10, d6, d8, d10
   real(wp) :: dE, dG(3), disp, ddisp
   real(wp) :: dists2(iter_chunk_size), coords(3, iter_chunk_size)

   real(wp), allocatable :: gw(:, :), dgwdcn(:, :)
   real(wp), allocatable :: c6(:, :), dc6dcn(:, :)
   real(wp), allocatable :: dEdcn(:)

   if (.not.allocated(reference_c6)) call copy_c6(reference_c6)

   nat = len(mol)
   max_ref = maxval(number_of_references(mol%at))
   allocate(gw(max_ref, nat), dgwdcn(max_ref, nat), c6(nat, nat), &
      &     dc6dcn(nat, nat), dEdcn(nat), source=0.0_wp)

   call weight_references(nat, mol%at, weighting_factor, cn, gw, dgwdcn)

   call get_atomic_c6(nat, mol%at, gw, dgwdcn, c6, dc6dcn)

   !$omp parallel do default(none) schedule(runtime) &
   !$omp reduction(+:energies, gradient, sigma, dEdcn) &
   !$omp shared(mol, neighlist, par, c6) &
   !$omp private(neighiter, neighs, ij, jat, ati, atj, coords, image, dists2, &
   !$omp&        r2, rij, r4r2ij, r0, t6, t8, t10, d6, d8, d10, disp, ddisp, &
   !$omp&        dE, dG)
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

            r4r2ij = 3*r4r2(ati)*r4r2(atj)
            r0 = par%a1*sqrt(r4r2ij) + par%a2

            t6 = 1._wp/(r2**3+r0**6)
            t8 = 1._wp/(r2**4+r0**8)
            t10 = 1._wp/(r2**5+r0**10)

            d6 = -6*r2**2*t6**2
            d8 = -8*r2**3*t8**2
            d10 = -10*r2**4*t10**2

            disp = (par%s6*t6 + par%s8*r4r2ij*t8 &
               &  + par%s10*49.0_wp/40.0_wp*r4r2ij**2*t10)*pair_scale(iat, jat)
            ddisp= par%s6*d6 + par%s8*r4r2ij*d8 &
               & + par%s10*49.0_wp/40.0_wp*r4r2ij**2*d10

            dE = -c6(iat, jat)*disp
            dG = -c6(iat, jat)*ddisp*rij

            energies(iat) = energies(iat) + 0.5_wp*dE
            energies(jat) = energies(jat) + 0.5_wp*dE

            dEdcn(iat) = dEdcn(iat) - dc6dcn(iat, jat) * disp
            dEdcn(jat) = dEdcn(jat) - dc6dcn(jat, iat) * disp

            gradient(:, iat) = gradient(:, iat) + dG
            gradient(:, jat) = gradient(:, jat) - dG

            sigma = sigma + spread(dG, 1, 3)*spread(rij, 2, 3)*pair_scale(iat, jat)

         enddo
      enddo
   enddo
   !$omp end parallel do

   gradient = gradient + reshape(matmul(reshape(dcndr, [3*nat, nat]), &
      &                                 dEdcn), shape(gradient))
   if (mol%npbc > 0) then
      sigma = sigma + reshape(matmul(reshape(dcndL, [9, nat]), &
         &                           dEdcn), shape(sigma))
   endif

end subroutine d3bj_eg

real(wp) pure elemental function weight_cn(wf,cn,cnref) result(cngw)
   real(wp),intent(in) :: wf, cn, cnref
   intrinsic :: exp
   cngw = exp ( -wf * ( cn - cnref )**2 )
end function weight_cn

real(wp) pure elemental function pair_scale(iat, jat) result(scale)
   integer, intent(in) :: iat, jat
   if (iat == jat) then
      scale = 0.5_wp
   else
      scale = 1.0_wp
   endif
end function pair_scale

end module d3mod_dftd3
