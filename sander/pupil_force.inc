!  this file is included in force.F90
#ifdef PUPIL_SUPPORT

  !*****************************************************
  !     Getting the Quantum forces with PUPIL package
  !*****************************************************

  ! Reconstruct the simulation cell if there is any change
  ! call inipupcell(natms,qcell,cell,xxx,yyy,zzz)
  do iPup = 1, 3    !vector loop
    do jPup= 1, 3  !Component loop
      qcell((iPup-1)*3+jPup) = ucell(jPup,iPup)
    enddo
  enddo
  ! Minimum point of the box ..... we assume (0,0,0) ????
  qcell(10) = 0.0d0
  qcell(11) = 0.0d0
  qcell(12) = 0.0d0

  ! Temporary vector to wrap the real coordinates to pass through
  ! PUPIL interface.  This stack will be deallocated after
  ! recovering all PUPIL forces.
  call get_stack(l_puptmp, 3*natom, routine)
  if (.not. rstack_ok) then
    deallocate(r_stack)
    allocate(r_stack(1:lastrst), stat=alloc_ier)
    call reassign_rstack(routine)
  endif
  REQUIRE(rstack_ok)
  do iPup=1,3*natom
    r_stack(l_puptmp + iPup - 1) = x(iPup)
  end do

  if(ntb > 0) then
    call wrap_molecules(nspm, ix(i70), r_stack(l_puptmp))
    if (ifbox == 2) then
      call wrap_to(nspm, ix(i70), r_stack(l_puptmp), box)
    end if
  end if

  ! Preparing the coordinates, velocity and classic forces
  ! to get quantum force
  do iPup=1,natom
    bs1 = (iPup-1)*9
    bs2 = (iPup-1)*3
    do jPup=1,3
      qcdata(bs1+jPup) = r_stack(l_puptmp + bs2 + jPup - 1)
      qcdata(bs1+3+jPup) = realStack(lvel+bs2+jPup-1)
      qcdata(bs1+6+jPup) = f(bs2+jPup)
    enddo
  enddo

  ! We are going to use the qmmm_nml and qmmm_struct variables
  ! to skip quantum atoms in the force calculation
  qmmm_nml%ifqnt = .true.
  if (pupStep .eq. 0) then

    ! To keep initial values from the MD step 1
    ierr = 0
    allocate ( pupnb14(numnb14*3),stat=ierr)
    REQUIRE(ierr == 0)
    allocate ( pupbonh(nbonh*3),stat=ierr)
    REQUIRE(ierr == 0)
    allocate ( pupbona(nbona*3),stat=ierr)
    REQUIRE(ierr == 0)
    allocate ( puptheth(ntheth*4),stat=ierr)
    REQUIRE(ierr == 0)
    allocate ( puptheta(ntheta*4),stat=ierr)
    REQUIRE(ierr == 0)
    allocate ( pupphih(nphih*5),stat=ierr)
    REQUIRE(ierr == 0)
    allocate ( pupphia(nphia*5),stat=ierr)
    REQUIRE(ierr == 0)
    pupnbonh   = nbonh
    pupnbona   = nbona
    pupntheth  = ntheth
    pupntheta  = ntheta
    pupnphih   = nphih
    pupnphia   = nphia
    pupnumnb14 = numnb14
    call copy_14nb(ix(inb_14),pupnb14,numnb14)
    do iPup = 1,nbonh
      bs1 = iPup-1
      bs2 = bs1*3
      pupbonh(bs2+1) = ix(iibh+bs1)
      pupbonh(bs2+2) = ix(ijbh+bs1)
      pupbonh(bs2+3) = ix(iicbh+bs1)
    enddo
    do iPup = 1,nbona
      bs1 = iPup-1
      bs2 = bs1*3
      pupbona(bs2+1) = ix(iiba+bs1)
      pupbona(bs2+2) = ix(ijba+bs1)
      pupbona(bs2+3) = ix(iicba+bs1)
    enddo
    do iPup = 1,ntheth
      bs1 = iPup-1
      bs2 = bs1*4
      puptheth(bs2+1) = ix(i24+bs1)
      puptheth(bs2+2) = ix(i26+bs1)
      puptheth(bs2+3) = ix(i28+bs1)
      puptheth(bs2+4) = ix(i30+bs1)
    enddo
    do iPup = 1,ntheta
      bs1 = iPup-1
      bs2 = bs1*4
      puptheta(bs2+1) = ix(i32+bs1)
      puptheta(bs2+2) = ix(i34+bs1)
      puptheta(bs2+3) = ix(i36+bs1)
      puptheta(bs2+4) = ix(i38+bs1)
    enddo
    do iPup = 1,nphih
      bs1 = iPup-1
      bs2 = bs1*5
      pupphih(bs2+1) = ix(i40+bs1)
      pupphih(bs2+2) = ix(i42+bs1)
      pupphih(bs2+3) = ix(i44+bs1)
      pupphih(bs2+4) = ix(i46+bs1)
      pupphih(bs2+5) = ix(i48+bs1)
    enddo
    do iPup = 1,nphia
      bs1 = iPup-1
      bs2 = bs1*5
      pupphia(bs2+1) = ix(i50+bs1)
      pupphia(bs2+2) = ix(i52+bs1)
      pupphia(bs2+3) = ix(i54+bs1)
      pupphia(bs2+4) = ix(i56+bs1)
      pupphia(bs2+5) = ix(i58+bs1)
    enddo
  endif

  ! Getting the quantum forces for a specific quantum domain
  pupStep  = pupStep + 1
  puperror = 0
  pupLevelData = 3
  call getquantumforces(natom,pupLevelData,pupStep,puperror,qcdata,qcell)
  if (puperror .ne. 0) then
    write (6,*) 'Fatal error: Could not obtain quantum forces!'
    call mexit(6,1)
  endif

  ! Quantum energy treatment....
  pot%scf = qmEnergy

  ! Delete interactions between the classical and quantum
  ! zones if a new list of quantum atoms is given.
  if (pupQZchange .ne. 0) then

    ! Rebuild the nonbonded 1:4 list.  Delete connectivity between
    ! the QM atoms and reinitialize internal NB 1:4 list structures
    ! from the beginning.
    numnb14 = pupnumnb14
    nbonh  = pupnbonh
    nbona  = pupnbona
    ntheth = pupntheth
    ntheta = pupntheta
    nphih  = pupnphih
    nphia  = pupnphia
    call copy_14nb(pupnb14, ix(inb_14), pupnumnb14)
    do iPup = 1,nbonh
      bs1 = iPup-1
      bs2 = bs1*3
      ix(iibh +bs1) = pupbonh(bs2+1)
      ix(ijbh +bs1) = pupbonh(bs2+2)
      ix(iicbh+bs1) = pupbonh(bs2+3)
    enddo
    do iPup = 1,nbona
      bs1 = iPup-1
      bs2 = bs1*3
      ix(iiba+bs1)  = pupbona(bs2+1)
      ix(ijba+bs1)  = pupbona(bs2+2)
      ix(iicba+bs1) = pupbona(bs2+3)
    enddo
    do iPup = 1,ntheth
      bs1 = iPup-1
      bs2 = bs1*4
      ix(i24+bs1) = puptheth(bs2+1)
      ix(i26+bs1) = puptheth(bs2+2)
      ix(i28+bs1) = puptheth(bs2+3)
      ix(i30+bs1) = puptheth(bs2+4)
    enddo
    do iPup = 1,ntheta
      bs1 = iPup-1
      bs2 = bs1*4
      ix(i32+bs1) = puptheta(bs2+1)
      ix(i34+bs1) = puptheta(bs2+2)
      ix(i36+bs1) = puptheta(bs2+3)
      ix(i38+bs1) = puptheta(bs2+4)
    enddo
    do iPup = 1,nphih
      bs1 = iPup-1
      bs2 = bs1*5
      ix(i40+bs1) = pupphih(bs2+1)
      ix(i42+bs1) = pupphih(bs2+2)
      ix(i44+bs1) = pupphih(bs2+3)
      ix(i46+bs1) = pupphih(bs2+4)
      ix(i48+bs1) = pupphih(bs2+5)
    enddo
    do iPup = 1,nphia
      bs1 = iPup-1
      bs2 = bs1*5
      ix(i50+bs1) = pupphia(bs2+1)
      ix(i52+bs1) = pupphia(bs2+2)
      ix(i54+bs1) = pupphia(bs2+3)
      ix(i56+bs1) = pupphia(bs2+4)
      ix(i58+bs1) = pupphia(bs2+5)
    enddo

    call deleting_qm_atoms()
    qsetup = .true.

    ! Setting as current quantum zone
    pupQZchange = 0

  endif

  ! For PUPIL, rebuild the neighbour list and
  ! zero the charges on QM atoms at every step
  if (igb == 0 .and. ipb == 0 .and. iyammp == 0) then
         
    ! (for GB: do all nonbondeds together below)
    call timer_start(TIME_NONBON)
    call timer_start(TIME_LIST)

    ! do_list_update = .true.         
    call nonbond_list(x, ix(i04), ix(i06), ix(i08), ix(i10), ntypes, &
                      natom/am_nbead, xx, ix, ipairs, ntnb, ix(ibellygp), &
                      belly,newbalance, qsetup, do_list_update)

    ! call qm_zero_charges(x(L15))
    call timer_stop(TIME_LIST)
    call timer_stop(TIME_NONBON)
  end if

  ! charge reassign here !
  if (ifcr /= 0) then
    call cr_reassign_charge( x, f, pot%ct, xx(l15), natom )
  end if
#endif /*PUPIL_SUPPORT*/
