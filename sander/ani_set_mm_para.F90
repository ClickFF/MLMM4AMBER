module ani_set_mm_para

    implicit none

    contains
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !+ Get rid of bond force constant for QM/MM or non-belly atom pairs.
    subroutine ani_setbon(nb,ib,jb,icb,igrp)
    use ani_module, only : ani_nml, ani_struct
    use parms,       only : numbnd, orig_numbnd, add_qmmm_bonds, req
    use constants,   only : zero
    implicit none
    ! Passed arguments
    integer, intent(inout) :: nb
    integer, intent(inout) :: ib(nb), jb(nb), icb(nb)
    integer, intent(in)    :: igrp(*)

    ! Internal variables
    integer nba, i, iat, jat
    logical iq, jq
    integer j,iatnum,jatnum
    real(8), dimension(nb) :: new_rk
    real(8), dimension(nb) :: new_req
    integer mlp_shake
    
    ! numbnd = number of bond parameters (changed)
    ! orig_numbnd = number of bond parameters (unchanged)
    ! numbnd - orig_numbnd = number of bond parameters added during this subroutine
    
    orig_numbnd = numbnd
    
    mlp_shake = ani_nml%mlp_shake
    nba = 0
    do i = 1, nb
        iat = ib(i)/3+1
        jat = jb(i)/3+1
        iq = .FALSE.
        jq = .FALSE.
        if (ani_nml%ifmlp) then
            ! remove quantum - quantum bond pairs from the bond list. This is done by rebuilding
            ! the list without these pairs.
            ! Note if lnk_method == 2 then we will treat the MMLink pair atom as being a QM atom so
            ! we also delete any MML-QM bonds.
            do j=1, ani_nml%ani_natom
            if (iat==ani_struct%ani_mask_index(j)) then
                iq = .true.
                iatnum = ani_struct%ani_species(j)
            elseif (jat==ani_struct%ani_mask_index(j)) then
                jq = .true.
                jatnum = ani_struct%ani_species(j)
            end if
            end do
        endif
        iq = iq .and. jq
        if((igrp(iat) > 0 .or. igrp(jat) > 0)) then
            if (iq) then !Both are QM atoms (or also QM-MMLink if lnk_method==2)
            !In order to shake QM atoms we need to add the bonds to the bond list
            !but we need to make the force constant zero while preserving the eqm bond length.
            !We will do this by creating a new type for each of the QM-QM bonds.
            !For the moment only do QM-H bonds for shake NTC=2
            !Only do this if qmshake is set.
                if ((iatnum == 1 .OR. jatnum == 1) .AND. mlp_shake == 1) then
                    !at least one of them is a hydrogen
                    !We need to make a new bond type here.
                    numbnd = numbnd+1
                    new_rk(numbnd - orig_numbnd) = zero  ! Set force constant to zero
                    new_req(numbnd - orig_numbnd) = req(icb(i)) !Preserve the eqm distance
                    nba = nba+1
                    ib(nba) = ib(i)
                    jb(nba) = jb(i)
                    icb(nba) = numbnd
                    !else we do nothing and the bond doesn't get added to the bond list.
                    write(6,*) 'Found ANI bond pair, set ',iat,jat,' to zero force constant'
                else
                    write(6,*) 'Found ANI bond pair, neglect ',iat,jat
                end if 
            else
            !We add the current pair to the list if both atoms are not bellied AND/OR
            !both in the QM region
                nba = nba+1
                ib(nba) = ib(i)
                jb(nba) = jb(i)
                icb(nba) = icb(i)
            end if
        end if
    end do
    nb = nba

    ! Add these additional bond parameters to the global rk/req arrays
    call add_qmmm_bonds(new_rk, new_req)

    return

    end subroutine ani_setbon 

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !+ Remove QM/MM or non-belly angle parameters
    subroutine ani_setang(nt,it,jt,kt,ict,igrp)
    use ani_module, only : ani_nml,ani_struct
    implicit NONE
    integer nt
    integer it(nt),jt(nt),kt(nt),ict(nt),igrp(*)
    integer nta, i, iat, jat, kat
    logical iq, jq, kq
    integer j
    
    nta = 0
    do i = 1,nt
        iat = it(i)/3+1
        jat = jt(i)/3+1
        kat = kt(i)/3+1
        iq = .FALSE.
        jq = .FALSE.
        kq = .FALSE.
        if (ani_nml%ifmlp) then
            ! remove quantum - quantum - quantum angle triplets from the angle list. This is done by rebuilding
            ! the list without these triplets.
            ! if lnk_method=2 we treat the MMlink pair atoms as being QM so we also remove ML-QM-QM angle triplets.
            do j=1, ani_nml%ani_natom
            iq = iq.or.iat==ani_struct%ani_mask_index(j)
            jq = jq.or.jat==ani_struct%ani_mask_index(j)
            kq = kq.or.kat==ani_struct%ani_mask_index(j)
            end do
            iq = iq .and. jq .and. kq
            !iq will be true if iat,jat and kat are all in the QM region

            if (iq) then
                !write(6,*) 'Found ANI angle triplet ',iat,jat,kat
            end if
        end if
        if((igrp(iat) > 0 .or. igrp(jat) > 0 .or. igrp(kat) > 0) .and. (.not. iq) ) then
            !We add the current triplet to the list if all atoms are not bellied AND/OR
            !in the QM region
            nta = nta+1
            it(nta) = it(i)
            jt(nta) = jt(i)
            kt(nta) = kt(i)
            ict(nta) = ict(i)
        else
            write(6,*) 'Found ANI angle triplet ',iat,jat,kat
        end if
    end do
    nt = nta
    return
    end subroutine ani_setang 

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !+ [Enter a one-line description of subroutine setdih here]
    subroutine ani_setdih(np,ip,jp,kp,lp,icp,igrp)
    use ani_module, only : ani_nml,ani_struct
    implicit NONE
    integer np
    integer ip(np),jp(np),kp(np),lp(np),icp(np),igrp(*)
    integer npa, i, iat, jat, kat, lat
    logical iq, jq, kq, lq
    integer j
    
    npa = 0
    do i = 1,np
        iat = ip(i)/3+1
        jat = jp(i)/3+1
        kat = iabs(kp(i))/3+1
        lat = iabs(lp(i))/3+1
        iq = .FALSE.
        jq = .FALSE.
        kq = .FALSE.
        lq = .FALSE.
        if (ani_nml%ifmlp) then
            ! remove quantum - quantum - quantum - quantum dihedrals from the dihedra; list. This is done by rebuilding
            ! the list without these dihedrals.
            do j=1, ani_nml%ani_natom
            iq = iq.or.iat==ani_struct%ani_mask_index(j)
            jq = jq.or.jat==ani_struct%ani_mask_index(j)
            kq = kq.or.kat==ani_struct%ani_mask_index(j)
            lq = lq.or.lat==ani_struct%ani_mask_index(j)
            end do
            iq = iq .and. jq .and. kq .and. lq
            !iq will be true if iat,jat, kat and lat are all in the QM region
            if (iq) then
                !write(6,*) 'Found ANI dihedral quadruplet ',iat,jat,kat,lat
            endif
        end if
        if((igrp(iat) > 0 .or. igrp(jat) > 0 .or. igrp(kat) > 0 .or. igrp(lat) > 0) .and. (.not. iq)) then
            !add current dihedral set to the list since all atoms are not bellied or in the QM region
            npa = npa+1
            ip(npa) = ip(i)
            jp(npa) = jp(i)
            kp(npa) = kp(i)
            lp(npa) = lp(i)
            icp(npa) = icp(i)
        else
            write(6,*) 'Found ANI dihedral quadruplet ',iat,jat,kat,lat
        end if
    end do
    np = npa
    return
    end subroutine ani_setdih 

end module ani_set_mm_para
