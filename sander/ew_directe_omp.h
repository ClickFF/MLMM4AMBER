
! epilogue: 12-6 LF terms

!$omp parallel do  &
!$omp&  private(im_c4,im_new,j,dfee,delx,dely,delz,delr2inv,ic,r6,f6,f12,df,dfx,dfy,dfz)  &
!$omp&  reduction(+:evdw,dumx,dumy,dumz)
!New2021
do im_new = 1,icount
   j = cache_bckptr(im_new)

   dfee = cache_df(im_new)
   delx = cache_x(im_new)
   dely = cache_y(im_new)
   delz = cache_z(im_new)
   delr2inv = cache_r2(im_new)

   ic = ico(iaci+iac(j))
   r6 = delr2inv*delr2inv*delr2inv
   delr12inv = r6 * r6
   f6 = cn2(ic)*r6
   f12 = cn1(ic)*delr12inv

   evdw = evdw + f12 - f6

   df = dfee + (12.d0*f12 - 6.d0*f6)*delr2inv

   dfx = delx*df
   dfy = dely*df
   dfz = delz*df
#if 0
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
end do  !  im_new = 1,icount
