!  to be included in runmd.F90

    call plumed_f_gcreate()
    call plumed_f_gcmd("setRealPrecision"//char(0), 8)
    call plumed_f_gcmd("getApiVersion"//char(0), plumed_version)
    if (plumed_version>1) then
      plumed_kbt = 2.0 * temp0 * boltz2
      call plumed_f_gcmd("setKbT"//char(0), plumed_kbt)
    endif
    plumed_energyUnits = 4.184
    plumed_lengthUnits = 0.1
    plumed_timeUnits = 1.0
    plumed_chargeUnits = 1.0/18.2223
    call plumed_f_gcmd("setMDEnergyUnits"//char(0), plumed_energyUnits)
    call plumed_f_gcmd("setMDLengthUnits"//char(0), plumed_lengthUnits)
    call plumed_f_gcmd("setMDTimeUnits"//char(0), plumed_timeUnits)
    if (plumed_version > 3) then
      call plumed_f_gcmd("setMDChargeUnits"//char(0), plumed_chargeUnits)
    endif
    call plumed_f_gcmd("setPlumedDat"//char(0), &
                       trim(adjustl(plumedfile))//char(0))
    call plumed_f_gcmd("setNatoms"//char(0), nr)
    call plumed_f_gcmd("setMDEngine"//char(0), "amber")
    call plumed_f_gcmd("setTimestep"//char(0), dt)
#  ifdef MPI
    call plumed_f_gcmd("setMPIFComm"//char(0), commsander)
    if (numgroup>1) then
      call plumed_f_gcmd("GREX setMPIFIntracomm"//char(0), commsander)
      if (master) then
        call plumed_f_gcmd("GREX setMPIFIntercomm"//char(0), commmaster)
      endif
      call plumed_f_gcmd("GREX init"//char(0), 0)
    endif
#  endif
    call plumed_f_gcmd("init"//char(0), 0)

