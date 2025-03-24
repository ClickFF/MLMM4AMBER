
    plumed_stopflag=0
    call plumed_f_gcmd("setStep"//char(0), nstep)
    call plumed_f_gcmd("setPositions"//char(0), x)
    call plumed_f_gcmd("setMasses"//char(0), amass)
    call plumed_f_gcmd("setCharges"//char(0), xx(l15))
    call plumed_f_gcmd("setEnergy"//char(0), ener%pot)
    call plumed_f_gcmd("setForces"//char(0), f)
    call plumed_f_gcmd("setStopFlag"//char(0), plumed_stopflag)
    plumed_box = 0.0
    if (ifbox == 0) then
      continue
    else if (ifbox == 1) then
      plumed_box(1,1) = box(1)
      plumed_box(2,2) = box(2)
      plumed_box(3,3) = box(3)
    else if (ifbox == 2) then
      plumed_box(1,1) = box(1)
      plumed_box(2,1) = 0.0
      plumed_box(3,1) = 0.0
      plumed_box(1,2) = -box(1)/3.0
      plumed_box(2,2) = sqrt(8.0/9.0)*box(1)
      plumed_box(3,2) = 0.0
      plumed_box(1,3) = -box(1)/3.0
      plumed_box(2,3) = -sqrt(2.0/9.0)*box(1)
      plumed_box(3,3) = sqrt(2.0/3.0)*box(1)
    else
      write (6,*) "!!!!! PLUMED ERROR: Only orthorhombic and truncted &
                   &octahedron cells are supported in this release."
      write (6,*) "!!!!! ABORTING RUN"
      call mexit(6, 1)
    endif
    plumed_virial=0.0

    ! It's not completely clear where the factor 2.0 comes from.
    ! Anyway, I was able to match a change in press of 1000 bar with
    ! a corresponding SLOPE=66.02 added to VOLUME CV in PLUMED GB.
    plumed_virial(1,1)=2.0*ener%vir(1)
    plumed_virial(2,2)=2.0*ener%vir(2)
    plumed_virial(3,3)=2.0*ener%vir(3)
    call plumed_f_gcmd("setVirial"//char(0),plumed_virial)
    call plumed_f_gcmd("setBox"//char(0),plumed_box)
    call plumed_f_gcmd("calc"//char(0),0);
#ifdef MPI
    ! This is required since PLUMED only updates virial on master processor
    call mpi_bcast(plumed_virial, 9, MPI_DOUBLE_PRECISION, 0, commsander, ierr)
#endif
    ener%vir(1)=0.5*plumed_virial(1,1)
    ener%vir(2)=0.5*plumed_virial(2,2)
    ener%vir(3)=0.5*plumed_virial(3,3)
