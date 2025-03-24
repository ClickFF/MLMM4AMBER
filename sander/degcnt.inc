
  ! Call DEGCNT to get the actual number of degrees of freedom for the
  ! solute and solvent. This call returns the correct numbers for belly
  ! simulations and simulations with separate solute/solvent scaling -- dap
  ! "IDUMAR" is dummy array. Used since this routine was also used w/ GIBBS.
  
#ifdef LES
  ! Return LES and non-LES degrees of freedom, since separate solvent coupling
  ! is no longer used.  Large changes to degcnt were made: cnum is now passed
  ! (LES copy number of each atom)
  call degcnt(ibelly, nr, ix(ibellygp), nsolut, nbonh, nbona, 0, ix(iibh), &
              ix(ijbh), ix(iiba), ix(ijba), idumar, idumar, ntc, idumar, 0, &
              0, 0, idumar, rndfp, rndfles, cnum, temp0les)
   
  ! RNDFP = # degrees of freedom for solute
  ! RNDFS = # degrees of freedom for solvent
  ! RNDF = total number of degrees of freedom.
  ! RNDFLES = # degrees of freedom for LES groups
  
  ! temp0les was init to negative number to signify not to use a LES bath
  ! just do standard code (meaning use solute/solvent baths)
  ! any positive (or zero) means to use LES bath with that target
  
  ! degcnt returns rndfs or rndfles in the rndfles variable
  ! depending on whether a LES bath was specified
  ! do this instead of duplicating call with rndfs or rndfles
   
  if (temp0les < 0.d0) then
    rndfs=rndfles
    rndfles=0.d0
  else
    rndfs=0.d0
  end if

  if (master) then
    write (6,'(a,f8.0)') "# degrees of freedom in non-LES region: ", rndfp
    write (6,'(a,f8.0)') "# degrees of freedom in     LES region: ", rndfles
  end if
   
  ! Modify RNDFP to reflect NDFMIN (set in mdread)
  rndfp = rndfp - ndfmin
  if (temp0les < 0.d0) then
    rndf = rndfp+rndfs
  else
    rndf = rndfp+rndfles
  end if
   
#else 
  call degcnt(ibelly, nr, ix(ibellygp), nsolut, nbonh, nbona, 0, ix(iibh), &
              ix(ijbh), ix(iiba), ix(ijba), idumar, idumar, ntc, idumar, 0, &
              0, 0, idumar, rndfp, rndfs)

  ! RNDFP = # degrees of freedom for solute
  ! RNDFS = # degrees of freedom for solvent
  ! RNDF = total number of degrees of freedom.
#  ifdef MPI
  if (mdloop .eq. 0 .and. master) then
#  else
  if (master) then
#  endif /* MPI */
    if (abfqmmm_param%abfqmmm /= 1 .or. &
        (abfqmmm_param%qmstep == 1 .and. abfqmmm_param%system == 1)) then
      write (6,'(a,f8.0)') "|  # of SOLUTE  degrees of freedom (RNDFP): ", &
        rndfp
      write (6,'(a,f8.0)') "|  # of SOLVENT degrees of freedom (RNDFS): ", &
        rndfs
    end if
  end if

  ! qtw - substract the number of overlapping noshake QM atoms in noshakemask
  rndfp = rndfp - qmmm_struct%noshake_overlap

  ! Modify RNDFP to reflect NDFMIN (set in mdread) and num_noshake
  rndfp = rndfp - ndfmin + num_noshake
  rndf = rndfp + rndfs
#  ifdef MPI
  if (mdloop .eq. 0 .and. master) then
#  else
  if (master) then
#  endif /* MPI */
    if (abfqmmm_param%abfqmmm /= 1 .or. &
        (abfqmmm_param%qmstep == 1 .and. abfqmmm_param%system == 1)) then
      if (qmmm_nml%ifqnt) then
        write (6,'(a,i6)') "|  QMSHAKE_NOSHAKEMASK_OVERLAP = ", &
          qmmm_struct%noshake_overlap
      end if
      write (6,'(a,f8.0,a,i6,a,f8.0)') "|  NDFMIN = ", rndfp, &
        "     NUM_NOSHAKE = ", num_noshake, "     CORRECTED RNDFP = ", rndfp
      write (6,'(a,f8.0)') "|  TOTAL # of degrees of freedom (RNDF) = ", rndf
    end if
  end if
#endif /* LES */
  ! Correct the degree of freedom count for extra points.  Note
  ! that rndfp, rndfs are uncorrected in an extra points context!
  call fix_degree_count(rndf)  
#ifndef LES
  if (isgld > 0) then

    ! Number of degrees of freedom in the SGLD part
    if (isgsta == 1) then
      sgsta_rndfp = 0
    else
      call degcnt(ibelly, nr, ix(ibellygp), isgsta-1, nbonh, nbona, 0, &
                  ix(iibh), ix(ijbh), ix(iiba), ix(ijba), idumar, idumar, &
                  ntc, idumar, 0, 0, 0, idumar, sgsta_rndfp, ignore_solvent)
    end if
    if (isgend == nr) then
      sgend_rndfp = rndf
    else
      call degcnt(ibelly, nr, ix(ibellygp), isgend, nbonh, nbona, 0, &
                  ix(iibh), ix(ijbh), ix(iiba), ix(ijba), idumar, &
                  idumar, ntc, idumar, 0, 0, 0, idumar, sgend_rndfp, &
                  ignore_solvent)
    end if

    ! Note that the solute ndf outputs above from degcnt are uncorrected
    ! for qmmm_struct%noshake_overlap, num_noshake, and extra points;
    ! also, ndfmin is not always being handled.
    call sg_fix_degree_count(sgsta_rndfp, sgend_rndfp, ndfmin, rndf)
  end if
#endif /* LES */

#ifdef MPI
  ! Soft-core TI: count additional degrees of freedom.
  if (ifsc /= 0) then
    call sc_degrees_o_freedom(ndfmin)
  end if
#endif /* MPI */
  ! End of degrees of freedom setup.
