!------------------------------------------------------------------------------
! This file is not meant to be compiled.  It is to be included in runmd.
!------------------------------------------------------------------------------

  ! The generalized canonical-isokinetic algorithm

    ! Propagation of the physical and virtual velocities as well as chain
    ! variables related to the generalized canonical-isokinetic ensemble
    ! using a specified set of atoms for each processor 
    iii = 3*(istart-1)
    do j = istart, iend

      ! The contribution of thermostat forces (first half step):
      do kija = 1, nkija
        glcl(j) = glcl(j) + (s2(j) + s3(j) - (s1(j)**2)*s2(j)*erlxt2)*sinsh4
        s3(j) = s3(j) + (s2(j)**2 - erlixt2)*sinsh8
        if (s3(j) == 0.d0) then
          s2(j) = s2(j) + (s1(j)**2 - erlixt2)*sinsh4
        else
          s2(j) = s2(j) * dexp(-s3(j)*sinsh4) + &
                  (s1(j)**2 - erlixt2)*(1.d0 - dexp(-s3(j)*sinsh4))/s3(j)
        end if
        s3(j) = s3(j) + (s2(j)**2 - erlixt2)*sinsh8
        hkin = amass(j)*(v(iii+1)**2 + v(iii+2)**2 + v(iii+3)**2) + &
               3.d0 * tktk * erlxt2 * (s1(j)**2) / 4.d0
        erlst2 = erlxt2 * 3.d0 * tktk / (hkin * 2.d0)
        esi = dexp(-s2(j) * sinsh2)
        esim = dsqrt((s1(j)**2) * erlst2 * (esi**2 - 1.d0)/2.d0 + 1.d0)
        v(iii+1) = v(iii+1) / esim
        v(iii+2) = v(iii+2) / esim
        v(iii+3) = v(iii+3) / esim
        s1(j) = s1(j) * esi / esim
        s3(j) = s3(j) + (s2(j)**2 - erlixt2)*sinsh8
        if (s3(j) == 0.d0) then
          s2(j) = s2(j) + (s1(j)**2 - erlixt2)*sinsh4
        else
          s2(j) = s2(j) * dexp(-s3(j)*sinsh4) + &
                  (s1(j)**2 - erlixt2)*(1.d0 - dexp(-s3(j)*sinsh4))/s3(j)
        end if
        s3(j) = s3(j) + (s2(j)**2 - erlixt2)*sinsh8
        glcl(j) = glcl(j) + (s2(j) + s3(j) - (s1(j)**2)*s2(j)*erlxt2)*sinsh4
      end do

      ! The contribution of physical forces (full step)
      hkin = amass(j) * (v(iii+1)**2 + v(iii+2)**2 + v(iii+3)**2) + &
             3.d0 * tktk * erlxt2 * (s1(j)**2) / 4.d0
      aaa = f(iii+1)*v(iii+1) + f(iii+2)*v(iii+2) + f(iii+3)*v(iii+3)
      bbb = (f(iii+1)**2 + f(iii+2)**2 + f(iii+3)**2)/amass(j)
      aaa = aaa / hkin
      bbb = bbb / hkin
      tt = dtx * dsqrt(bbb)
      st = aaa * (dcosh(tt) - 1.d0)/bbb + dsinh(tt)/dsqrt(bbb)
      st1 = dcosh(tt) + aaa*dsinh(tt)/dsqrt(bbb)
      s1(j) = s1(j) / st1
      v(iii+1) = (v(iii+1) + f(iii+1)/amass(j)*st) / st1
      v(iii+2) = (v(iii+2) + f(iii+2)/amass(j)*st) / st1
      v(iii+3) = (v(iii+3) + f(iii+3)/amass(j)*st) / st1

      ! The contribution of thermostat forces (second half step)
      do kija = 1, nkija
        glcl(j) = glcl(j) + (s2(j) + s3(j) - (s1(j)**2)*s2(j)*erlxt2)*sinsh4
        s3(j) = s3(j) + (s2(j)**2 - erlixt2)*sinsh8
        if (s3(j) == 0.d0) then
          s2(j) = s2(j) + (s1(j)**2 - erlixt2)*sinsh4
        else
          s2(j) = s2(j) * dexp(-s3(j) * sinsh4) + &
                  (s1(j)**2 - erlixt2)*(1.d0 - dexp(-s3(j)*sinsh4))/s3(j)
        end if
        s3(j) = s3(j) + (s2(j)**2 - erlixt2)*sinsh8
        hkin = amass(j) * (v(iii+1)**2 + v(iii+2)**2 + v(iii+3)**2) + &
               3.d0 * tktk * erlxt2 * (s1(j)**2) / 4.d0
        erlst2 = erlxt2 * 3.d0 * tktk / (hkin * 2.d0)
        esi = dexp(-s2(j) * sinsh2)
        esim = dsqrt((s1(j)**2) * erlst2 * (esi**2 - 1.d0)/2.d0 + 1.d0)
        v(iii+1) = v(iii+1) / esim
        v(iii+2) = v(iii+2) / esim
        v(iii+3) = v(iii+3) / esim
        s1(j) = s1(j) * esi / esim
        s3(j) = s3(j) + (s2(j)**2 - erlixt2)*sinsh8
        if (s3(j) == 0.d0) then
          s2(j) = s2(j) + (s1(j)**2 - erlixt2)*sinsh4
        else
          s2(j) = s2(j) * dexp(-s3(j) * sinsh4) + &
                  (s1(j)**2 - erlixt2)*(1.d0 - dexp(-s3(j)*sinsh4))/s3(j)
        end if
        s3(j) = s3(j) + (s2(j)**2 - erlixt2)*sinsh8
        glcl(j) = glcl(j) + (s2(j) + s3(j) - (s1(j)**2)*s2(j)*erlxt2)*sinsh4
      end do
      iii = iii + 3
    end do

    ! Calculation of some quantities using parallel processing
    iii = 3*(istart - 1)
    do j = istart, iend

      ! Measuring the maximal kinetic constraints deviations
      tkkk = (amass(j) * (v(iii+1)**2 + v(iii+2)**2 + v(iii+3)**2) + &
              3.d0 * tktk * erlxt2 * (s1(j)**2) / 4.d0) / (3.d0 * boltz)
      if (dabs(tkkk - temp0) > davalev) then
        davalev = dabs(tkkk - temp0)
      end if

      ! Measuring the conservation laws fluctuations:
      etlci = (s2(j)**2 + s3(j)**2)*erlxt2/2.d0 + glcl(j)
      clfs = clfs + (etlci - etlc(j))**2
      iii = iii+3
    end do

#ifdef MPI
    ! Extracting the quantities from each processor
    ! and collecting them in the root
    call MPI_Gather(davalev, 1, MPI_DOUBLE_PRECISION, sdavalev, 1, &
                    MPI_DOUBLE_PRECISION, 0, commsander, ierr);

    ! Summing up over all processors and sending the result to the master root
    call mpi_reduce(clfs, sclfs, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
                    commsander, ierr)
    if (master) then
      davalevs = 0.d0
      do inum = 1, numtasks
        if (sdavalev(inum) > davalevs) then
          davalevs = sdavalev(inum)
        end if
      end do
    end if
#endif /* MPI */
    iconfs = iconfs + 1

    ! Writing the results only once using the master processor:
#ifdef MPI
    if (master) then
#endif /* MPI */
      if ((nstep+1) == ntpr*((nstep+1)/ntpr)) then
#ifdef MPI
        write(6,256) '| Max deviations in the kinetic constraints:', &
          davalevs / temp0 * 100.d0, ' %'
        write(6,266) '| Mean fluctuations in the conservation laws:', &
          dsqrt(sclfs / (iconfs*natom))   
#else
        write(6,256) '| Max deviations in the kinetic constraints:', &
          davalev / temp0 * 100.d0, ' %'
        write(6,266) '| Mean fluctuations in the conservation laws:', &
          dsqrt(clfs / (iconfs*natom))   
#endif /* MPI */
  256   format(a44,e12.3,a2)
  266   format(a45,e12.3)
      end if
#ifdef MPI
    end if
#endif /* MPI */

    ! Calculation of velocity and thermostat distribution functions
    if (idistr .ne. 0) then
      if ((nstep+1) == idistr * ((nstep+1)/idistr)) then
        txtr = 0.d0
        iii = 3*(istart - 1)
        do j = istart, iend
          txtr = txtr + amass(j)*(v(iii+1)**2 + v(iii+2)**2 + v(iii+3)**2)
          iii = iii+3
        end do
#ifdef MPI

        ! Summing up the partial terms over all processors and 
        ! sending the result to all
        call mpi_allreduce(txtr, stxtr, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                           commsander, ierr)
        txtr = stxtr / (natom * 3.d0)
#else
        txtr = txtr / (natom * 3.d0)
#endif
        iii = 3*(istart - 1)
        do j = istart, iend

          ! Distribution of atomic velocities:
          kves = idint(dabs(dsqrt(amass(j))*v(iii+1))/(dsqrt(txtr)*dargvs)) + 1
          if (v(iii+1) < 0.d0) then
            kves = -kves
          end if
          if (iabs(kves) <= lep) then
            distrs(kves) = distrs(kves) + 1.d0
          end if
          kves = idint(dabs(dsqrt(amass(j))*v(iii+2))/(dsqrt(txtr)*dargvs)) + 1
          if (v(iii+2) < 0.d0) then
            kves = -kves
          end if
          if (iabs(kves) <= lep) then
            distrs(kves) = distrs(kves) + 1.d0
          end if
          kves = idint(dabs(dsqrt(amass(j))*v(iii+3))/(dsqrt(txtr)*dargvs)) + 1
          if (v(iii+3) < 0.d0) then
            kves = -kves
          end if
          if (iabs(kves) <= lep) then
            distrs(kves) = distrs(kves) + 1.d0
          end if
 
          ! And their thermostat counterparts
          kves = idint(s1(j) * erlxt / dargvs) + 1
          if (kves <= lep) then
            distrh1(kves) = distrh1(kves) + 1.d0
          end if

          ! Include chain variables:
          kves = idint(dabs(s2(j)*erlxt)/dargvs) + 1
          if (s2(j) < 0.d0) then
            kves = -kves
          end if
          if (iabs(kves) <= lep) then
            distrh2(kves) = distrh2(kves) + 1.d0
          end if
          kves = idint(dabs(s3(j)*erlxt)/dargvs) + 1
          if (s3(j) < 0.d0) then
            kves = -kves
          end if
          if (iabs(kves) <= lep) then
            distrh3(kves) = distrh3(kves) + 1.d0
          end if
          iii = iii+3
        end do
        iconfd = iconfd + 1
      end if
    end if

#ifdef MPI
    ! Summing up over all processors and sending the results to the master root
    call mpi_reduce(distrs, sdistrs, 2*lep+1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                    0, commsander, ierr)
    call mpi_reduce(distrh1, sdistrh1, lep, MPI_DOUBLE_PRECISION, MPI_SUM, &
                    0, commsander, ierr)
    call mpi_reduce(distrh2, sdistrh2, 2*lep+1, MPI_DOUBLE_PRECISION, &
                    MPI_SUM, 0, commsander, ierr)
    call mpi_reduce(distrh3, sdistrh3, 2*lep+1, MPI_DOUBLE_PRECISION, &
                    MPI_SUM, 0, commsander, ierr)

    ! Writing the results only once using the master processor
    if (master) then
#endif /* MPI */

    ! Writing the velocity distribution functions:
    if (iconfd .ne. 0 .and. idistr .ne. 0) then
      if ((nstep+1) == ntpr*((nstep+1)/ntpr)) then
        open(26, file='velocds.dat')
        do jl = -lep, lep
          if (jl .ne. 0) then
#ifdef MPI
            distr = sdistrs(jl) / (iconfd * natom * dargvs * 3.d0)
#else
            distr = distrs(jl) / (iconfd * natom * dargvs * 3.d0)
#endif /* MPI */
            distrmaxwel = ((1.d0 / (2.d0*pif1))**(1.d0/2.d0)) * &
                           dexp(-0.5d0 * argvel(jl)**2)
            write(26, 604) argvel(jl), distr, distrmaxwel
  604       format(1x, f10.5, 2(2x, f10.5))
          end if
        end do
        close(26)
        ! Writing the accompanied distribution function:
        open(27, file='velocds1.dat')
        do jl = 1, lep
#ifdef MPI
          distr = sdistrh1(jl) / (iconfd * natom * dargvs)
#else
          distr = distrh1(jl) / (iconfd * natom * dargvs)
#endif /* MPI */
          write(27, 605) argvel(jl), distr
  605     format(1x, f10.5, 2x, f10.5)
        end do
        close(27)
        ! Writing the chain distributions:
        open(27, file='velocd23.dat')
        do jl = -lep, lep
          if (jl .ne. 0) then
#ifdef MPI
            distr = sdistrh2(jl) / (iconfd * natom * dargvs)
            distra = sdistrh3(jl) / (iconfd * natom * dargvs)
#else
            distr = distrh2(jl) / (iconfd * natom * dargvs)
            distra = distrh3(jl) / (iconfd * natom * dargvs)
#endif /* MPI */
            write(27, 604) argvel(jl), distr, distra
          end if
        end do
        close(27)
      end if
    end if
#ifdef MPI
    end if
    ! End of contingency for writing out velocity distributions,
    ! which occurs only on the master process.  The presence of
    ! this contingency keeps us at two levels of indentation.
#endif /* MPI */

    ! Writing the atomic velocities and chain-thermostat variables
    ! every ntwr steps
    if ((nstep+1) == ntwr*((nstep+1)/ntwr) .or. nstep+1 == nstlim) then
#ifdef MPI
      sv = 0.d0
      s1a = 0.d0
      s2a = 0.d0
      s3a = 0.d0
      iii = 3*(istart-1)
      do j = istart, iend
        sv(iii+1) = v(iii+1)
        sv(iii+2) = v(iii+2)
        sv(iii+3) = v(iii+3)
        iii = iii + 3
        s1a(j) = s1(j)
        s2a(j) = s2(j)
        s3a(j) = s3(j)
      end do
      svs = 0.d0
      call mpi_reduce(sv, svs, 3*natom, MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
                      commsander, ierr)
      s1s = 0.d0
      s2s = 0.d0
      s3s = 0.d0
      call mpi_reduce(s1a, s1s, natom, MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
                      commsander, ierr)
      call mpi_reduce(s2a, s2s, natom, MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
                      commsander, ierr)
      call mpi_reduce(s3a, s3s, natom, MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
                      commsander, ierr)
      if (master) then
        open(578, file='vfreez.rst')
        iii = 0
        do j = 1, natom
          write(578, *) svs(iii+1), svs(iii+2), svs(iii+3)
          iii = iii+3
        end do
        close(578)
        open(579, file='tfreez.rst')
        do j = 1, natom
          write(579, *) s1s(j), s2s(j), s3s(j)
        end do
        close(579)
      end if
#else
      open(578, file='vfreez.rst')
      iii = 0
      do j = 1, natom
        write(578, *) v(iii+1), v(iii+2), v(iii+3)
        iii = iii+3
      end do
      close(578)
      open(579, file='tfreez.rst')
      do j = 1, natom
        write(579,*) s1(j),s2(j),s3(j)
      end do
      close(579)
#endif /* MPI */
    end if

