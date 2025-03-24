
! epilogue: 12-6 LF terms

do im_new = 1,icount
   j = cache_bckptr(im_new)

   dfee = cache_df(im_new)
   delx = cache_x(im_new)
   dely = cache_y(im_new)
   delz = cache_z(im_new)
   delr2inv = cache_r2(im_new)
   !write(6,*) i
   !write(6,*) iaci
   !write(6,*) j
   !write(6,*) iac(j)
   ic = ico(iaci+iac(j))
   !write(6,*) ic 
   !write(6,*) " "
   r6 = delr2inv*delr2inv*delr2inv
   delr12inv = r6 * r6
#ifdef LES 
   lfac=lesfac(lestmp+lestyp(j))
   f6 = cn2(ic)*r6*lfac
   f12 = cn1(ic)*delr12inv*lfac

   if(ipimd>0) then
      if(cnum(i).eq.0.and.cnum(j).eq.0) then
         nrg_all(1:nbead)=nrg_all(1:nbead) + (f12-f6)/nbead
      else
         if(cnum(i).ne.0) then
            nrg_all(cnum(i)) = nrg_all(cnum(i)) + f12-f6
         else
            nrg_all(cnum(j)) = nrg_all(cnum(j)) + f12-f6
         endif
      endif
   endif
#else
   if ( vdwmodel == 0 ) then
      if (fswitch < 0) then
         f6 = cn2(ic)*r6
         f12 = cn1(ic)*delr12inv
      else
         delr3inv = delr2inv * sqrt(delr2inv)
         if (delr2inv < fswitch2inv) then
            p12 = cut6 / (cut6 - fswitch6)
            p6 = cut3 / (cut3 - fswitch3)
            f12 = cn1(ic) * p12 * (r6 - cut6inv) * (r6 - cut6inv)
            f6 = cn2(ic) * p6 * (delr3inv - cut3inv) * (delr3inv - cut3inv)
            df12 = 12.d0 * p12 * delr2inv * r6 * (r6 - cut6inv)
            df6 = 6.d0 * p6 * (delr3inv - cut3inv) * delr3inv * delr2inv
         else
            f12 = cn1(ic) * delr12inv - cn1(ic) * invfswitch6cut6
            f6 = cn2(ic) * r6 - cn2(ic) * invfswitch3cut3
            df12 = 12.d0 * delr2inv * delr12inv
            df6 = 6.d0 * delr2inv * r6
         end if
         df = dfee + df12 * cn1(ic) - df6 * cn2(ic)
      end if
   else
      f6 = cn5(ic)*r6
      f12 = cn4(ic)/exp(cn3(ic)/sqrt(delr2inv))
   endif
#endif

   if  (lj1264 == 1 .and. plj1264 == 1) then !New2021
      mr4 = delr2inv*delr2inv
      f4 = cn6(ic)*mr4
      do im_c4 = 1, C4Pairwise
         if(MIN(i,j) == (MIN(cn7(im_c4*3-2), cn7(im_c4*3-1))+3)/3 .and. MAX(i,j) == (MAX(cn7(im_c4*3-2), cn7(im_c4*3-1))+3)/3) then
            f4 = f4 + cn8(im_c4)*mr4
         end if
      end do
      if(decpr .and. idecomp > 0) &
         call decpair(3,i,j,(f12 - f6 - f4)/(nstlim/ntpr))
      evdw = evdw + f12 - f6 - f4
   else if (lj1264 == 1 .and. plj1264 /= 1) then 
      mr4 = delr2inv*delr2inv
      f4 = cn6(ic)*mr4
      if(decpr .and. idecomp > 0) &
         call decpair(3,i,j,(f12 - f6 - f4)/(nstlim/ntpr))
      evdw = evdw + f12 - f6 - f4
   else if (lj1264 /= 1 .and. plj1264 == 1) then
      mr4 = delr2inv*delr2inv
      do im_c4 = 1, C4Pairwise
         if(MIN(i,j) == (MIN(cn7(im_c4*3-2), cn7(im_c4*3-1))+3)/3 .and. MAX(i,j) == (MAX(cn7(im_c4*3-2), cn7(im_c4*3-1))+3)/3) then
            f4 = cn8(im_c4)*mr4
         end if
      end do
      if(decpr .and. idecomp > 0) &
         call decpair(3,i,j,(f12 - f6 - f4)/(nstlim/ntpr))
      evdw = evdw + f12 - f6 - f4
   else
   ! -- ti decomp
      if(decpr .and. idecomp > 0) call decpair(3,i,j,(f12 - f6)/(nstlim/ntpr))
      evdw = evdw + f12 - f6
   endif

   if (lj1264 == 1 .or. plj1264 == 1) then !New2021
      df = dfee + (12.d0*f12 - 6.d0*f6 - 4.d0*f4 )*delr2inv
   else if (fswitch < 0) then
      df = dfee + (12.d0*f12 - 6.d0*f6)*delr2inv
   endif

   dfx = delx*df
   dfy = dely*df
   dfz = delz*df
#ifndef noVIRIAL
   vxx = vxx - dfx*delx
   vxy = vxy - dfx*dely
   vxz = vxz - dfx*delz
   vyy = vyy - dfy*dely
   vyz = vyz - dfy*delz
   vzz = vzz - dfz*delz
#endif
   dumx = dumx + dfx
   dumy = dumy + dfy
   dumz = dumz + dfz
   force(1,j) = force(1,j) + dfx
   force(2,j) = force(2,j) + dfy
   force(3,j) = force(3,j) + dfz
   f4 = 0.0 ! clear f4 to avoid repetitive calculation
end do  !  im_new = 1,icount

