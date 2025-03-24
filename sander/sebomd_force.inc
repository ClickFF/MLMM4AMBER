!  this file is included in force.F90

    ! step 2/2 to save forces from nfe_on_force
    call sebomd_save_forces(2, natom, f, xx(gradsebomd), xx(grad1tmp), &
                            xx(grad2tmp))
    sebomd_obj%diverror = 0
    if (ntp > 0 .or. ntb == 2) then
      do i = 1, 3
        sebomd_obj%dbox(i) = box(i)
      enddo
    endif

    ! First call to SEBOMD to initialize the parameters
400 if (sebomd_obj%idcflag == 0) then
      if (ntb == 0) then
        sebomd_obj%dbox(1) = 0.d0
        sebomd_obj%dbox(2) = 0.d0
        sebomd_obj%dbox(3) = 0.d0
      else
        sebomd_obj%dbox(1) = box(1)
        sebomd_obj%dbox(2) = box(2)
        sebomd_obj%dbox(3) = box(3)
!       ! force wrapping (needed by sebomd)
!       call wrap_molecules(nspm,ix(i70),x)
      endif
      
      if (sebomd_obj%chtype > 0) then
!       call sebomd_readatchg(xx(divchg))
      else
        do i = 1, natom
          xx(divchg + i - 1) = xx(l15 + i - 1)/18.2223d0
        enddo
      endif
      call sebomd_energy(sebomd_obj%idcflag, sebomd_obj%iflagch, &
                         sebomd_obj%iflagbo, ntb, natom, vir, &
                         sebomd_obj%esebomd, sebomd_obj%dbox, x, &
                         xx(gradsebomd), nres, ix(i02), xx(divchg), ih(m06), &
                         nspm, ix(i70), sebomd_obj%pdmx)
    endif
    ! End of first call to sebomd, as indicated by idcflag == 0

    ! Checking that everything is ok, exiting if not
    if (sebomd_obj%idcflag == -1) then
      write (6,*) 'Error during SEBOMD calculation'
      call mexit(6, 1)
    endif
    
    sebomd_obj%idcflag = 1
    
!    if (ntb.ne.0) then
!       ! force wrapping (needed by sebomd)
!       call wrap_molecules(nspm,ix(i70),x)
!    endif

    ! second call to SEBOMD to calculate the potential energy
    call sebomd_energy(sebomd_obj%idcflag, sebomd_obj%iflagch, &
                       sebomd_obj%iflagbo, ntb, natom, vir, &
                       sebomd_obj%esebomd, sebomd_obj%dbox, &
                       x, xx(gradsebomd), nres, ix(i02), xx(divchg), ih(m06), &
                       nspm, ix(i70), sebomd_obj%pdmx)
!    pot%scf = escf
    
    ! Checking that everything is ok, initializing again if not
    if (sebomd_obj%idcflag == -1) then
      if (sebomd_obj%diverror == 0) then
        sebomd_obj%idcflag = 0
        sebomd_obj%diverror = 1
        write (6,*) 'Error during SEBOMD calculation, restarting with a ', &
               'new initial density matrix'
        goto 400
      else
        write (6,*) 'Error during SEBOMD calculation'
        call mexit(6, 1)
      endif
    endif

    ! SEBOMD is currently computing gradient (not forces)
    do i = 1, 3*natom
      xx(gradsebomd+i-1) = -xx(gradsebomd+i-1)
    enddo

    ! SEBOMD gradient computation by numerical differentiation (debugging)
    if (sebomd_obj%debugforces /= 0) then
      dx = 1.0d-5
      do i = 1, 3*natom
        ! central energy E_SEBOMD(x)
        xdx = x(i)
        call sebomd_save_forces(0, natom, f, xx(gradsebomd), xx(grad1tmp), &
                                xx(grad2tmp))
        call sebomd_energy(sebomd_obj%idcflag, sebomd_obj%iflagch, &
                           sebomd_obj%iflagbo, ntb, natom, vir, &
                           sebomd_obj%esebomd, sebomd_obj%dbox, x, &
                           xx(gradsebomd), nres, ix(i02) ,xx(divchg), &
                           ih(m06), nspm,ix(i70), sebomd_obj%pdmx)
        e0 = sebomd_obj%esebomd
        fxx = -xx(gradsebomd+i-1)

        ! E_SEBOMD(x-\delta x)
        x(i) = xdx-dx
        call sebomd_save_forces(0, natom, f, xx(gradsebomd), xx(grad1tmp), &
                                xx(grad2tmp))
        call sebomd_energy(sebomd_obj%idcflag, sebomd_obj%iflagch, &
                           sebomd_obj%iflagbo, ntb, natom, vir, &
                           sebomd_obj%esebomd, sebomd_obj%dbox, x, &
                           xx(gradsebomd), nres, ix(i02), xx(divchg), &
                           ih(m06), nspm,ix(i70), sebomd_obj%pdmx)
        em1 = sebomd_obj%esebomd

        ! E_SEBOMD(x-2\delta x)
        x(i) = xdx-dx-dx
        call sebomd_save_forces(0, natom, f, xx(gradsebomd), xx(grad1tmp), &
                                xx(grad2tmp))
        call sebomd_energy(sebomd_obj%idcflag, sebomd_obj%iflagch, &
                           sebomd_obj%iflagbo, ntb, natom, vir, &
                           sebomd_obj%esebomd, sebomd_obj%dbox, x, &
                           xx(gradsebomd), nres, ix(i02), xx(divchg), &
                           ih(m06), nspm,ix(i70), sebomd_obj%pdmx)
        em2 = sebomd_obj%esebomd

        ! E_SEBOMD(x+\delta x)
        x(i) = xdx+dx
        call sebomd_save_forces(0, natom, f, xx(gradsebomd), xx(grad1tmp), &
                                xx(grad2tmp))
        call sebomd_energy(sebomd_obj%idcflag, sebomd_obj%iflagch, &
                           sebomd_obj%iflagbo, ntb, natom, vir, &
                           sebomd_obj%esebomd, sebomd_obj%dbox, x, &
                           xx(gradsebomd), nres, ix(i02), xx(divchg), &
                           ih(m06), nspm,ix(i70), sebomd_obj%pdmx)
        ep1 = sebomd_obj%esebomd

        ! E_SEBOMD(x+2\delta x)
        x(i) = xdx+dx+dx
        call sebomd_save_forces(0, natom, f, xx(gradsebomd), xx(grad1tmp), &
                                xx(grad2tmp))
        call sebomd_energy(sebomd_obj%idcflag, sebomd_obj%iflagch, &
                           sebomd_obj%iflagbo, ntb, natom, vir, &
                           sebomd_obj%esebomd, sebomd_obj%dbox, x, &
                           xx(gradsebomd), nres, ix(i02), xx(divchg), &
                           ih(m06), nspm,ix(i70), sebomd_obj%pdmx)
        ep2 = sebomd_obj%esebomd

        fx1 = -(ep1-em1)/(2.0d0*dx)
        fx2 = (1.0d0/ 12.0d0*(ep2-em2) + 2.0d0/  3.0d0*(em1-ep1))/dx

        write(0,'("gradients ",2i4,3f25.16)') (i-1)/3+1,i,fxx,fx2,fx1
        write(0,'("delta grad",2i4,2f25.16)') (i-1)/3+1,i,fxx-fx2,fxx-fx1

        ! Finally, restore the current coordinates
        x(i) = xdx
      enddo
    endif
