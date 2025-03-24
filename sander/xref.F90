! <compile=optimized>
#include "copyright.h"
#include "../include/assert.fh"
#include "../include/dprec.fh"

!------------------------------------------------------------------------------
! cns_ref:  Interface to the Axel Brunger's CNS (Crystallography and NMR
!           System) program.
!------------------------------------------------------------------------------
module cns_xref
   
  ! Description:
  ! Module for interfacing to the Crystallography & NMR System
  ! (CNS) software suite:
  ! http://cns.csb.yale.edu
  
  ! History:
  ! $Id: xref.f,v 10.4 2010/03/19 18:13:05 rcw Exp $
  use state
  implicit none

  private

  logical, public :: is_xref_on
  public :: cns_xref_run
  public :: write_cns_xref_md_energies
  public :: write_cns_xref_min_energies
  character(*), parameter :: output_token = "CNS_XREF: "

  contains
  ! public routines are given in alphabetical order

  !---------------------------------------------------------------------------
  ! cns_xref_run: run cns_solve and process the output.
  !
  ! Arguments:
  !   natom:    the number of atoms in the system
  !   ih:       global Hollerith array for atom names, residue names,
  !             atom types, and other character information
  !   x:        system coordinates for all atoms
  !   fg:       
  !---------------------------------------------------------------------------
  subroutine cns_xref_run(natom, ih, x, fg, ener)

      implicit none
      character(len=4), intent(in)    :: ih(*)
      _REAL_,           intent(inout) :: x(*)
      _REAL_,           intent(inout) :: fg(*)
      type(state_rec),  intent(inout) :: ener

      integer MAXATREF
      parameter (MAXATREF=20000)
      character(len=20) fcnsgen
      integer natom, itimes, natcns, medcns(MAXATREF)

      !Local scratch
      _REAL_ :: ener_scratch(4)


      data itimes /0/, medcns /MAXATREF*0/
      save fcnsgen, itimes, natcns, medcns

      if ( is_xref_active() ) then

        if (itimes .eq. 0) call rdxrefin(natom,ih,x,fcnsgen,natcns,medcns)

        itimes = itimes + 1
        call wrtpdbcns(fcnsgen, medcns, x)
#ifdef NO_SYSTEM_INTERFACE
        call sander_bomb('cns_xref_run', 'system calls are not supported on &
                         &this architecture', 'CNS XREF requires system call &
                         &support to run.')
#else
        call system("cns_solve < minimize.inp > minimize.out")
        call rdxraycns(natcns, medcns, ener, ener_scratch, fg)
        call system("mv -f minimize.pdb cnsout.pdb")
        call system("mv -f mmen2.dat cnsout.ene")
        call system("mv -f force2.dat cnsout.grd")
#endif

        ! Disable xref calculations due to the complete absence of test caes
        ! and strange, undocumented use of the ener arrays which prevented this
        ! being successfully converted over to the new energy module structure.
        call sander_bomb('cns_xref_run', 'XREF calculations are currently &
                         &disabled due to their incompatibility with the new &
                         &state module structure.', 'Code updates and test &
                         &cases are required to for this to work.')

        !ener(16) = ener_scratch(1)
        !ener(18) = ener_scratch(3)
        !ener(19) = ener_scratch(4)
      end if
     return
  end subroutine cns_xref_run

 
  !----------------------------------------------------------------------------
  ! write_cns_xref_md_energies: emit energies during molecular dynamics.
  !
  ! Arguments:
  !   energy:   a special state_rec struct to store components of the energy
  !             (see state.F90, written by Mark James Williamson and inspired
  !             by similar features in pmemd)
  !----------------------------------------------------------------------------
  subroutine write_cns_xref_md_energies( energy )

    implicit none
    type(state_rec),  intent(in) :: energy

    character(*), parameter :: frmt = "(1x,'EXRAY = ', f14.4,2x, 'R     = &
                                       &',f14.4,2x, 'RFREE     = ',f14.4)"

    if (is_xref_active()) then
      write(6, frmt) energy%exray, energy%eptot, energy%rfree
    end if

    return
  end subroutine write_cns_xref_md_energies

  !----------------------------------------------------------------------------
  ! write_cns_xref_min_energies: emit energies during minimization, similar to
  !                              write_cns_xref_md_energies above.
  !
  ! Arguments:
  !   energy:   a special state_rec struct to store components of the energy
  !             (see state.F90)
  !----------------------------------------------------------------------------
  subroutine write_cns_xref_min_energies( energy )

    implicit none
    type(state_rec), intent(in) :: energy
    character(*), parameter :: frmt = "(1x,'EXRAY   = ',f13.4,2x,'R       = &
                                       &',f13.4,2x, 'RFREE      = ',f13.4)"

    if (is_xref_active()) then
      write(6, frmt) energy%exray, energy%eptot, energy%rfree
      write(6, "(1x,'EAMBER  = ',f13.4)" ) energy%eptot !i mjw maybe wrong
    end if

    return
  end subroutine write_cns_xref_min_energies

  ! Private routines in alphabetical order

  !----------------------------------------------------------------------------
  ! is_xref_active: returns true if xref is active, false if xref is unactive
  !----------------------------------------------------------------------------
  function is_xref_active()

    implicit none
    logical :: is_xref_active

    is_xref_active = is_xref_on
  end function is_xref_active

  !----------------------------------------------------------------------------
  ! rdxrefin: read input for the xref module
  !
  ! Arguments:
  !   natom:    the number of atoms in the system
  !   ih:       global Hollerith array for atom names, residue names,
  !             atom types, and other character information
  !   x:        system coordinates for all atoms
  !   fcnsgen:
  !   natcns:
  !   medcns:
  !----------------------------------------------------------------------------
  subroutine rdxrefin(natom, ih, x, fcnsgen, natcns, medcns)

    implicit none
    character(len=4) ih(*)
    _REAL_ x(*)
    character(len=20) fcnsgen
    integer natom, ixref, natcns, iread, medcns(*)
    character(len=80) line
    logical there
    parameter ( ixref = 1 )

    fcnsgen='generate_easy.pdb'
    inquire(file='xref.in',exist=there)
    if (.not.there) then
      call sander_bomb('rdxrefin<xref.f>', ' file xref.in does not exist', &
                       ' construct an xref.in file.')
    endif
    open(ixref,file='xref.in',status='old')
    read(ixref,'(A80)',end=1000) line
    do while (line(1:7) .ne. 'MAPPING')
      fcnsgen = line(1:20)
      inquire(file=fcnsgen,exist=there)
      if (.not.there) then
        call sander_bomb('rdxrefin<xref.f>', ' file cnsgen does not exist', &
                         ' check file xref.in.')
      endif
      read(ixref,'(A80)',end=1000) line
    end do

    ! If xref.in contains a mapping table, read it in.
    natcns = 0
    do while (.true.)
      read(ixref, '(A80)', end=2000) line
      if (line(1:11) == 'END_MAPPING') then
        close(ixref)
        return
      end if
      backspace(ixref)
      read(ixref, '(I5)') iread
      natcns = natcns + 1
      medcns(natcns) = iread
    end do

    ! If xref.in does not contain a mapping table, we will construct it.
 1000 close(ixref)
    call map_table(natom,ih,x,fcnsgen,natcns,medcns)
    return

 2000 close(ixref)
    call sander_bomb('rdxrefin<xref.f>', ' file xref.in incomplete', &
                     ' make sure xref.in contains MAPPING and END_MAPPING &
                     &tags.')
  end subroutine rdxrefin

  !----------------------------------------------------------------------------
  ! map_table: compare the atoms in the CNS pdb file against those in Amber,
  !            then build the mapping table.  If an atom in the CNS PDB file
  !            is not matched, it could be because 1.) it is in an alternate
  !            conformation of a disordered residue, or, 2.) it is related to
  !            another atom already matched by an NCS symmetry.  A warning
  !            will be generated but the program will proceed.
  !
  ! Arguments:
  !   natom:    the number of atoms in the system
  !   ih:       global Hollerith array for atom names, residue names,
  !             atom types, and other character information
  !   x:        system coordinates for all atoms
  !   fcnsgen:
  !   natcns:
  !   medcns:
  !----------------------------------------------------------------------------
  subroutine map_table(natom, ih, x, fcnsgen, natcns, medcns)

    implicit none
    character(len=4) ih(*),jnam
    _REAL_ x(*),xi,yi,zi,occi,bfaci,xtol
    character(len=20) fcnsgen
    integer natom,ipdb,natcns,medcns(*),jatm,jj
    parameter ( ipdb=1, xtol=1.1D-3 )
    character line*80,str1*13,inam*1,str2*16,str3*14

    open(ipdb, file=fcnsgen, status='old')
    natcns = 0
    do while (.true.)
      read(ipdb, '(A80)', end=1000) line
      if (line(1:4) == 'ATOM' .or. line(1:6) == 'HETATM') then
        natcns=natcns+1
        backspace(ipdb)
        read(ipdb, '(a13,a1,a16,f8.3,f8.3,f8.3,f6.2,f6.2,a14)', &
             end=1000) str1, inam, str2, xi, yi, zi, occi, bfaci, str3
        if (str3(7:9) == 'AC2') then
          medcns(natcns) = 0
        else
          do jatm = 1, natom
            jj = 3*(jatm - 1)
            jnam = ih(jatm)
            if (abs(x(jj+1)-xi) <= xtol .and. abs(x(jj+2)-yi) <= xtol .and. &
                abs(x(jj+3)-zi) <= xtol) then
              medcns(natcns)=jatm
              goto 2000
            end if
          end do
          write(6,'("WARNING: NO MATCH FOUND FOR - ",a13)') str1
        end if
 2000   continue
      end if
    end do
 1000 close(ipdb)
    call amflsh(6)

  end subroutine map_table

  !----------------------------------------------------------------------------
  ! wrtpdbcns: write a PDB file from CNS coordinates
  !
  ! Arguments:
  !   fcnsgen:
  !   medcns:
  !   x:        Amber coordinates array
  !----------------------------------------------------------------------------
  subroutine wrtpdbcns(fcnsgen, medcns, x)

    implicit none
    character(len=20) fcnsgen
    integer ipdb,iout,iatm,jatm,medcns(*),j,jbase
    _REAL_ x(*)
    character(len=80) line
    parameter ( ipdb=1, iout=2 )

    open(ipdb,file=fcnsgen,status='old')
    open(iout,file='cnsin.pdb')
    rewind(iout)
    iatm = 0
    do while (.true.)
      read(ipdb,'(A80)',end=1000) line
      if (line(1:4) == 'ATOM' .or. line(1:6) == 'HETATM') then
        if (line(73:74) /= 'AC' .or. line(73:75) == 'AC1') then
          iatm = iatm + 1
          jatm = medcns(iatm)
          jbase = 3*(jatm - 1)
          write(iout, '(A30,3(F8.3),A26)') line(1:30), &
                                           (x(jbase+j), j = 1, 3), line(55:80)
        else
          write(iout,'(A80)') line
        end if
      else
        write(iout,'(A80)') line
      end if
    end do
 1000 close(ipdb)
    close(iout)

  end subroutine wrtpdbcns

  !----------------------------------------------------------------------------
  ! rdxraycns: read a CNS Xray data file
  !
  ! Arguments:
  !   natcns:
  !   medcns:
  !   ener:
  !   ener_scratch:
  !   fg:
  !----------------------------------------------------------------------------
  subroutine rdxraycns(natcns, medcns, ener, ener_scratch, fg)

    use state
    implicit none
    type(state_rec) ::  ener
    _REAL_ fg(*), xgrd1, xgrd2, xgrd3
    _REAL_, intent(out) :: ener_scratch(4)
    integer natcns, medcns(*), ipdb, iread, iatm, jatm, jbase
    character(len=80) line
    parameter (ipdb=1)

    open(ipdb, file='mmen2.dat', status='old')
    read(ipdb, '(A80)') line
    read(ipdb, '(E20.14)') ener_scratch(1)
    read(ipdb, '(A80)') line
    read(ipdb, '(E20.14)') ener_scratch(2)
    read(ipdb, '(A80)') line
    read(ipdb, '(F8.6,1X,F8.6)') ener_scratch(4), ener_scratch(3)
    close(ipdb)
    ener%tot = ener%tot + ener_scratch(1)

    open(ipdb, file='force2.dat', status='old')
    read(ipdb, '(A80)') line
    read(ipdb, '(I6)') iread
    if (iread .ne. natcns) then
      call sander_bomb('rdxraycns<xref.f>', ' number of atoms in forces2.dat &
                       &incorrect', '.')
    end if
    do iatm = 1, natcns
      read(ipdb,'(3(E20.14,1X))') xgrd1, xgrd2, xgrd3
      jatm = medcns(iatm)
      if (jatm > 0) then
        jbase = 3*(jatm - 1)
        fg(jbase+1) = fg(jbase+1) - xgrd1
        fg(jbase+2) = fg(jbase+2) - xgrd2
        fg(jbase+3) = fg(jbase+3) - xgrd3
      end if
    end do
    close(ipdb)
  end subroutine rdxraycns

end module cns_xref
