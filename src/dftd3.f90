module d3mod_dftd3
   use iso_fortran_env, only: wp => real64
   implicit none
   public :: d3bj_eg
   private

contains 

subroutine d3bj_eg(mol, par, weighting_factor, cn, dcndr, dcndL, r_thr, &
      &            energy, gradient, sigma)
   use d3def_molecule
   use d3def_damping_parameters
   use d3par_dftd3
   use d3par_r4r2, only: r4r2 => sqrt_z_r4_over_r2
   use d3mod_utils_periodic, only: get_realspace_cutoff

   class(d3_molecule), intent(in) :: mol
   class(d3_damping_parameters), intent(in) :: par

   real(wp), intent(in) :: weighting_factor
   real(wp), intent(in) :: cn(:)
   real(wp), intent(in) :: dcndr(:, :, :)
   real(wp), intent(in) :: dcndL(:, :, :)
   real(wp), intent(in) :: r_thr

   real(wp), intent(inout) :: energy
   real(wp), intent(inout) :: gradient(:, :)
   real(wp), intent(inout) :: sigma(:, :)

   integer :: iat, jat, ati, atj, iref, jref, ic, max_ref
   integer :: rep(3), tx, ty, tz

   real(wp) :: norm, dnorm, tgw, dgw
   real(wp) :: r4r2ij, r0, t(3), rij(3), r2, r, oor6, oor8, oor10
   real(wp) :: door6, door8, door10, disp, ddisp, dtmp(3)
   real(wp) :: c6ij, dic6ij, djc6ij, c6ref, gtmp(3)

   real(wp), allocatable :: gw(:, :), dgwdcn(:, :), dEdcn(:)

   if (mol%npbc > 0) &
   call get_realspace_cutoff(mol%lattice,r_thr,rep)
   where(.not.mol%pbc) rep = 0

   if (.not.allocated(reference_c6)) call copy_c6(reference_c6)
   
   max_ref = maxval(number_of_references(mol%at))
   allocate(gw(max_ref, len(mol)), dgwdcn(max_ref, len(mol)), dEdcn(len(mol)), &
      &     source=0.0_wp)

   ! omp parallel default(none) &
   ! omp shared(mol, par, cn, rep, r_thr, gw, dgwdcn, weighting_factor, &
   ! omp&       reference_c6) &
   ! omp private(iat, jat, ati, atj, iref, jref, ic, tx, ty, tz, &
   ! omp&        norm, dnorm, tgw, dgw, r4r2ij, r0, t, rij, r2, r, &
   ! omp&        oor6, oor8, oor10, door6, door8, door10, gtmp, &
   ! omp&        c6ij, dic6ij, djc6ij, c6ref, disp, ddisp, dtmp) &
   ! omp reduction(+:energy,gradient,sigma,dEdcn)
   ! omp do schedule(runtime)
   do iat = 1, len(mol)
      ati = mol%at(iat)
      norm = 0.0_wp
      dnorm = 0.0_wp
      do iref = 1, number_of_references(ati)
         tgw = cngw(weighting_factor, cn(iat), reference_cn(iref, ati))
         norm = norm + tgw
         dnorm = dnorm + 2*weighting_factor*(reference_cn(iref, ati)-cn(iat))*tgw
      enddo
      norm = 1.0_wp/(norm + 1.0e-14_wp)
      do iref = 1, number_of_references(ati)
         tgw = cngw(weighting_factor, cn(iat), reference_cn(iref, ati))
         dgw =2*weighting_factor*(reference_cn(iref, ati)-cn(iat))*tgw

         gw(iref, iat) = tgw * norm
         dgwdcn(iref, iat) = dgw*norm - tgw*dnorm*norm**2
      enddo
   enddo
   ! omp enddo

   ! omp do schedule(runtime)
   do iat = 1, len(mol)
      ati = mol%at(iat)
      c6ij = 0.0_wp
      dic6ij = 0.0_wp
      djc6ij = 0.0_wp
      do iref = 1, number_of_references(ati)
         do jref = 1, number_of_references(ati)
            ic = ati * (1 + ati)/2
            c6ref = reference_c6(iref, jref, ic)
            c6ij = c6ij + gw(iref, iat) * gw(jref, iat) * c6ref
            dic6ij = dic6ij + dgwdcn(iref, iat) * gw(jref, iat) * c6ref
            djc6ij = djc6ij + gw(iref, iat) * dgwdcn(jref, iat) * c6ref
         enddo
      enddo
      r4r2ij = 3*r4r2(ati)*r4r2(ati)
      r0 = par%a1*sqrt(r4r2ij) + par%a2
      do concurrent(tx = -rep(1):rep(1), &
            &       ty = -rep(2):rep(2), &
            &       tz = -rep(3):rep(3), &
            &       tx /= 0 .or. ty /= 0 .or. tz /= 0)
         ! cycle iat with iat interaction in same cell
         t = [tx,ty,tz]
         rij = matmul(mol%lattice,t)
         r2  = sum( rij**2 )
         if (r2.gt.r_thr) cycle
         r   = sqrt(r2)
         oor6 = 1._wp/(r2**3+r0**6)
         oor8 = 1._wp/(r2**4+r0**8)
         oor10 = 1._wp/(r2**5+r0**10)
         door6 = -6*r2**2*r*oor6**2
         door8 = -8*r2**3*r*oor8**2
         door10 = -10*r2**4*r*oor10**2
         disp = par%s6*oor6 + par%s8*r4r2ij*oor8   &
            & + par%s10*49.0_wp/40.0_wp*r4r2ij**2*oor10
         ddisp= par%s6*door6 + par%s8*r4r2ij*door8 &
            & + par%s10*49.0_wp/40.0_wp*r4r2ij**2*door10
         energy = energy - c6ij*disp
         dtmp = c6ij*ddisp*rij/r
         dEdcn(iat) = dEdcn(iat) + (dic6ij + djc6ij)*disp
         sigma = sigma - spread(dtmp, 1, 3)*spread(rij, 2, 3)
      enddo
      do jat = 1, iat-1
         atj = mol%at(jat)
         ! temps
         gtmp = 0.0_wp
         c6ij = 0.0_wp
         dic6ij = 0.0_wp
         djc6ij = 0.0_wp
         ! all refs
         do iref = 1, number_of_references(ati)
            do jref = 1, number_of_references(atj)
               if (ati > atj) then
                  ic = atj + ati*(ati-1)/2
                  c6ref = reference_c6(iref, jref, ic)
               else
                  ic = ati + atj*(atj-1)/2
                  c6ref = reference_c6(jref, iref, ic)
               endif
               c6ij   = c6ij   + gw(iref, iat) * gw(jref, jat) * c6ref
               dic6ij = dic6ij + dgwdcn(iref, iat) * gw(jref, jat) * c6ref
               djc6ij = djc6ij + gw(iref, iat) * dgwdcn(jref, jat) * c6ref
            enddo
         enddo
         r4r2ij = 3*r4r2(ati)*r4r2(atj)
         r0 = par%a1*sqrt(r4r2ij) + par%a2
         do concurrent(tx = -rep(1):rep(1), &
               &       ty = -rep(2):rep(2), &
               &       tz = -rep(3):rep(3))
            t = [tx,ty,tz]
            rij = mol%xyz(:,iat) - mol%xyz(:,jat) + matmul(mol%lattice,t)
            r2 = sum(rij**2)
            if (r2.gt.r_thr) cycle
            r = sqrt(r2)
            oor6 = 1._wp/(r2**3+r0**6)
            oor8 = 1._wp/(r2**4+r0**8)
            oor10 = 1._wp/(r2**5+r0**10)
            door6 = -6*r2**2*r*oor6**2
            door8 = -8*r2**3*r*oor8**2
            door10 = -10*r2**4*r*oor10**2
            disp = par%s6*oor6 + par%s8*r4r2ij*oor8 &
               & + par%s10*49.0_wp/40.0_wp*r4r2ij**2*oor10
            ddisp= par%s6*door6 + par%s8*r4r2ij*door8 &
               & + par%s10*49.0_wp/40.0_wp*r4r2ij**2*door10
            energy = energy - c6ij*disp
            ! save this
            dEdcn(iat) = dEdcn(iat) + dic6ij *disp
            dEdcn(jat) = dEdcn(jat) + djc6ij *disp
            dtmp = c6ij*ddisp*rij/r
            gtmp = gtmp - dtmp
            sigma = sigma - spread(dtmp, 1, 3)*spread(rij, 2, 3)
         enddo
         gradient(:, iat) = gradient(:, iat) + gtmp
         gradient(:, jat) = gradient(:, jat) - gtmp
      enddo
   enddo
   ! omp enddo
   ! omp endparallel

   gradient = gradient - reshape(matmul(reshape(dcndr, [3*len(mol), len(mol)]), &
      &                                 dEdcn), shape(gradient))
   if (mol%npbc > 0) then
      sigma = sigma - reshape(matmul(reshape(dcndL, [9, len(mol)]), &
         &                           dEdcn), shape(sigma))
   endif
   
end subroutine d3bj_eg

integer pure elemental function lin(i, j)
   integer, intent(in) :: i, j
   integer :: ii, jj
   intrinsic :: max, min
   ii = max(i, j)
   jj = min(i, j)
   lin = jj + ii*(ii-1)/2
end function lin

real(wp) pure elemental function cngw(wf,cn,cnref)
   real(wp),intent(in) :: wf,cn,cnref
   intrinsic :: exp
   cngw = exp ( -wf * ( cn - cnref )**2 )
end function cngw

end module d3mod_dftd3
