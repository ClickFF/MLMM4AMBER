#include "../include/dprec.fh"

! This is a .F03 file.
!
! module for interfacing with PuReMD code (ReaxFF+EEM in QM/MM mode)
module qm2_extern_reaxff_puremd_module
        use, intrinsic :: iso_c_binding
        implicit none

        private

        interface
                type(c_ptr) function setup_qmmm &
                        (num_qm_atoms, qm_symbols, qm_pos, &
                        num_mm_atoms, mm_symbols, mm_pos_q, sim_box_info, &
                        ffield_filename, control_filename) &
                        bind(C, name='setup_qmmm') 
                        use, intrinsic :: iso_c_binding
                        implicit none
                        integer (c_int), value :: num_qm_atoms
                        type(c_ptr), value :: qm_symbols
                        type(c_ptr), value :: qm_pos
                        integer (c_int), value :: num_mm_atoms
                        type(c_ptr), value :: mm_symbols
                        type(c_ptr), value :: mm_pos_q
                        type(c_ptr), value :: sim_box_info
                        type(c_ptr), value :: ffield_filename
                        type(c_ptr), value :: control_filename
                end function setup_qmmm

                integer (c_int) function reset_qmmm &
                        (handle, num_qm_atoms, qm_symbols, qm_pos, &
                        num_mm_atoms, mm_symbols, mm_pos_q, sim_box_info, &
                        ffield_filename, control_filename) &
                        bind(C, name='reset_qmmm') 
                        use, intrinsic :: iso_c_binding
                        implicit none
                        type(c_ptr), value :: handle
                        integer (c_int), value :: num_qm_atoms
                        type(c_ptr), value :: qm_symbols
                        type(c_ptr), value :: qm_pos
                        integer (c_int), value :: num_mm_atoms
                        type(c_ptr), value :: mm_symbols
                        type(c_ptr), value :: mm_pos_q
                        type(c_ptr), value :: sim_box_info
                        type(c_ptr), value :: ffield_filename
                        type(c_ptr), value :: control_filename
                end function reset_qmmm

                integer (c_int) function simulate &
                        (handle) &
                        bind(C, name='simulate') 
                        use, intrinsic :: iso_c_binding
                        implicit none
                        type(c_ptr), value :: handle
                end function simulate

                integer (c_int) function cleanup &
                        (handle) &
                        bind(C, name='cleanup') 
                        use, intrinsic :: iso_c_binding
                        implicit none
                        type(c_ptr), value :: handle
                end function cleanup

                integer (c_int) function set_output_enabled &
                        (handle, is_enabled) &
                        bind(C, name='set_output_enabled') 
                        use, intrinsic :: iso_c_binding
                        implicit none
                        type(c_ptr), value :: handle
                        integer (c_int), value :: is_enabled
                end function set_output_enabled

                integer (c_int) function set_control_parameter &
                        (handle, keyword, values) &
                        bind(C, name='set_control_parameter') 
                        use, intrinsic :: iso_c_binding
                        implicit none
                        type(c_ptr), value :: handle
                        type(c_ptr), value :: keyword
                        type(c_ptr) :: values
                end function set_control_parameter

                integer (c_int) function set_contiguous_charge_constraints &
                        (handle, num_char_const_contig, char_const_contig_start, &
                        char_const_contig_end, char_const_contig_value) &
                        bind(C, name='set_contiguous_charge_constraints') 
                        use, intrinsic :: iso_c_binding
                        implicit none
                        type(c_ptr), value :: handle
                        integer (c_int), value :: num_char_const_contig
                        type(c_ptr), value :: char_const_contig_start
                        type(c_ptr), value :: char_const_contig_end
                        type(c_ptr), value :: char_const_contig_value
                end function set_contiguous_charge_constraints

                integer (c_int) function set_custom_charge_constraints &
                        (handle, num_char_const_custom, char_const_custom_count, &
                        char_const_custom_atom_index, &
                        char_const_custom_coeff, char_const_custom_rhs) &
                        bind(C, name='set_custom_charge_constraints') 
                        use, intrinsic :: iso_c_binding
                        implicit none
                        type(c_ptr), value :: handle
                        integer (c_int), value :: num_char_const_custom
                        type(c_ptr), value :: char_const_custom_count
                        type(c_ptr), value :: char_const_custom_atom_index
                        type(c_ptr), value :: char_const_custom_coeff
                        type(c_ptr), value :: char_const_custom_rhs
                end function set_custom_charge_constraints

                integer (c_int) function get_atom_forces_qmmm &
                        (handle, qm_f, mm_f) &
                        bind(C, name='get_atom_forces_qmmm') 
                        use, intrinsic :: iso_c_binding
                        implicit none
                        type(c_ptr), value :: handle
                        type(c_ptr), value :: qm_f
                        type(c_ptr), value :: mm_f
                end function get_atom_forces_qmmm

                integer (c_int) function get_atom_charges_qmmm &
                        (handle, qm_q, mm_q) &
                        bind(C, name='get_atom_charges_qmmm') 
                        use, intrinsic :: iso_c_binding
                        implicit none
                        type(c_ptr), value :: handle
                        type(c_ptr), value :: qm_q
                        type(c_ptr), value :: mm_q
                end function get_atom_charges_qmmm

                integer (c_int) function get_system_info &
                        (handle, e_potential, e_kinetic, e_total, temperature, &
                        volume, pressure) &
                        bind(C, name='get_system_info') 
                        use, intrinsic :: iso_c_binding
                        implicit none
                        type(c_ptr), value :: handle
                        type(c_ptr), value :: e_potential
                        type(c_ptr), value :: e_kinetic
                        type(c_ptr), value :: e_total
                        type(c_ptr), value :: temperature
                        type(c_ptr), value :: volume
                        type(c_ptr), value :: pressure
                end function get_system_info
        end interface

        public :: get_reaxff_puremd_forces, reaxff_puremd_finalize

        character(len=*), parameter, public :: module_name = "qm2_extern_reaxff_puremd_module"
        type(c_ptr), save, private :: handle = C_NULL_PTR
        
        type reaxff_nml_type
             ! control file name
             character(len=256) :: control
             ! force field file name
             character(len=256) :: ffield
             ! atom start number for charge constraint range
             integer, dimension(256) :: char_const_contig_start
             ! atom end number for charge constraint range
             integer, dimension(256) :: char_const_contig_end
             ! charge constraint for specified atom range
             real, dimension(256) :: char_const_contig_value
             ! counts for custom charge constraint
             integer, dimension(256) :: char_const_custom_count
             ! atom indices for custom charge constraints
             integer, dimension(1024) :: char_const_custom_atom_index
             ! coefficients for custom charge constraints
             real, dimension(1024) :: char_const_custom_coeff
             ! right-hand size (RHS) constants for custom charge constraints
             real, dimension(256) :: char_const_custom_rhs
             ! number of threads for reaxff-puremd (OpenMP only)
             integer :: numthreads
             ! charge method
             integer :: charge_method
             ! charge solver tolerance
             real :: solvtol
             ! max number of iterations for the charge solver
             integer :: solvmaxit
             ! preconditioner type for the charge solver
             integer :: solvprecond
             ! bonded interaction cutoff
             real :: nbrcut
             ! hydrogen bonding interaction cutoff
             real :: hbondcut
             ! valence cutoff (based on bond-order)
             real :: thbcut
             ! 1 = include polarization energy in ReaxFF, 0 = otherwise
             integer :: include_polar_energy
        end type reaxff_nml_type


contains

        subroutine get_reaxff_puremd_forces( num_qm_atoms, qm_pos, qm_types, &
                        qm_q, num_mm_atoms, mm_pos_q, mm_types, e_total, &
                        qm_f, mm_f, qmcharge )
                use, intrinsic :: iso_c_binding
                use ElementOrbitalIndex, only : elementSymbol
                implicit none

                integer (c_int), intent(in) :: num_qm_atoms                     ! number of QM atoms
                _REAL_, intent(in) :: qm_pos(3,num_qm_atoms)                    ! QM atom coordinates
                integer (c_int), target, intent(in) :: qm_types(num_qm_atoms)   ! QM atom types
                _REAL_, intent(inout) :: qm_q(num_qm_atoms)                     ! QM atom charges (nuclear charge in au)
                integer (c_int), intent(in) :: num_mm_atoms                     ! number of MM atoms
                _REAL_, intent(in) :: mm_pos_q(4,num_mm_atoms)                  ! MM atom coordinates and charges (nuclear charge in au)
                integer (c_int), target, intent(in) :: mm_types(num_mm_atoms)   ! MM atom types
                _REAL_, intent(out) :: e_total                                  ! SCF energy (kcal/mol)
                _REAL_, intent(out) :: qm_f(3,num_qm_atoms)                     ! SCF QM force (AMU * Angstroms / ps^2)
                _REAL_, intent(out) :: mm_f(3,num_mm_atoms)                     ! SCF MM force (AMU * Angstroms / ps^2)
                integer (c_int), intent(in) :: qmcharge                         ! total charge of the QM region

                ! double precision types required for PuReMD interface
                real (c_double), target :: qm_pos_(3,num_qm_atoms)
                real (c_double), target :: qm_q_(num_qm_atoms)
                real (c_double), target :: mm_pos_q_(4,num_mm_atoms)
                real (c_double), target :: e_total_
                real (c_double), target :: qm_f_(3,num_qm_atoms)
                real (c_double), target :: mm_f_(3,num_mm_atoms)

                logical, save :: first_call = .true.
                type(reaxff_nml_type), save :: reaxff_nml
                logical :: control_exists, ffield_exists
                integer :: ix
                integer (c_int) :: ret, num_char_const_contig, num_char_const_custom
                character(kind=c_char, len=1024), target :: control_filename, ffield_filename
                character(kind=c_char, len=1024), target :: keyword, values
                ! triplets for lengths and angles of QM region simulation box (Angstroms and degrees)
                real (c_double), target :: sim_box_info(6)
                character(len=2), dimension(num_qm_atoms), target :: qm_symbols
                character(len=2), dimension(num_mm_atoms), target :: mm_symbols
                integer (c_int), dimension(256), target :: char_const_contig_start, char_const_contig_end
                real (c_double), dimension(256), target :: char_const_contig_value
                integer (c_int), dimension(256), target :: char_const_custom_count
                integer (c_int), dimension(1024), target :: char_const_custom_atom_index
                real (c_double), dimension(1024), target :: char_const_custom_coeff
                real (c_double), dimension(256), target :: char_const_custom_rhs

#if defined(HAVE_REAXFF_PUREMD)
                ! NOTE: PuReMD must run with periodic boundary conditions (PBCs) ON,
                !       so to compensate the simulation box will have void space added around it
                !       (20 angstroms, as the long-range cut-off is 10 angstroms) in order
                !       negate the effect of PBCs for QMMM
                do ix = 1, 3
                        sim_box_info(ix) = MAX(MAXVAL(mm_pos_q(ix,:)), MAXVAL(qm_pos(ix,:))) &
                                - MIN(MINVAL(mm_pos_q(ix,:)), MINVAL(qm_pos(ix,:))) + 20.0
                end do
                ! orthogonal simulation box
                sim_box_info(4:6) = 90.0

                ! get character representations of atoms
                do ix = 1, num_qm_atoms
                        if ( qm_types(ix) >= 1 ) then
                                qm_symbols(ix) = elementSymbol(qm_types(ix))
                        else if ( qm_types(ix) == -1 ) then
                                ! map Amber dummy atom representation to reaxff-puremd representation
                                qm_symbols(ix) = 'X' // C_NULL_CHAR
                        else
                                call sander_bomb("get_reaxff_puremd_forces", &
                                        "ERROR: Invalid QM atom type!", &
                                        "Will exit now")
                        endif
                end do
                do ix = 1, num_mm_atoms
                        if ( mm_types(ix) >= 1 ) then
                                mm_symbols(ix) = elementSymbol(mm_types(ix))
                        else if ( mm_types(ix) == -1 ) then
                                ! map Amber dummy atom representation to reaxff-puremd representation
                                mm_symbols(ix) = 'X' // C_NULL_CHAR
                        else
                                call sander_bomb("get_reaxff_puremd_forces", &
                                        "ERROR: Invalid MM atom type!", &
                                        "Will exit now")
                        endif
                end do

                ! copy Amber mixed precision arrays into double precision arrays
                qm_pos_ = qm_pos
                qm_q_ = qm_q
                mm_pos_q_ = mm_pos_q
                num_char_const_contig = 0
                num_char_const_custom = 0

                if ( first_call ) then
                        first_call = .false.

                        write (6, '(/,a,/)') '  >>> Running calculations with ReaxFF <<<'

                        call get_namelist_reaxff(reaxff_nml)
                        call print_namelist_reaxff(reaxff_nml)

                        control_filename = trim(adjustl(reaxff_nml%control)) // C_NULL_CHAR
                        ffield_filename = trim(adjustl(reaxff_nml%ffield)) // C_NULL_CHAR
                        char_const_contig_start = reaxff_nml%char_const_contig_start
                        char_const_contig_end = reaxff_nml%char_const_contig_end
                        char_const_contig_value = reaxff_nml%char_const_contig_value
                        ! count atom ranges for contiguous constraints (1-based numbering of atoms)
                        num_char_const_contig = count(char_const_contig_start >= 1 .and. char_const_contig_end >= 1)
                        ! adjust constraint ends based on total num. atoms
                        do ix = 1, num_char_const_contig
                                char_const_contig_end(ix) = MIN( char_const_contig_end(ix), num_qm_atoms + num_mm_atoms )
                        enddo
                        char_const_custom_count = reaxff_nml%char_const_custom_count
                        char_const_custom_atom_index = reaxff_nml%char_const_custom_atom_index
                        char_const_custom_coeff = reaxff_nml%char_const_custom_coeff
                        char_const_custom_rhs = reaxff_nml%char_const_custom_rhs
                        ! count custom constraints
                        num_char_const_custom = count(char_const_custom_count >= 1)
 
                        inquire( file=control_filename, exist=control_exists )
                        inquire( file=ffield_filename, exist=ffield_exists )

                        ! check if the force field file (default or specified by user) exists
                        if ( ffield_exists .eqv. .false. ) then
                                call sander_bomb("get_reaxff_puremd_forces", &
                                        "ERROR: Specified force field file does not exist!", &
                                        "Will exit now")
                        endif

                        ! check if the control file specified by user exists
                        if ( (len_trim(reaxff_nml%control) > 0) .and. (control_exists .eqv. .false.) ) then
                                call sander_bomb("get_reaxff_puremd_forces", &
                                        "ERROR: Specified control file does not exist!", &
                                        "Will exit now")
                        endif

                        if ( (len_trim(reaxff_nml%control) > 0) .and. control_exists ) then
                                handle = setup_qmmm( num_qm_atoms, c_loc(qm_symbols), &
                                        c_loc(qm_pos_), num_mm_atoms, c_loc(mm_symbols), &
                                        c_loc(mm_pos_q_), c_loc(sim_box_info), &
                                        c_loc(ffield_filename), c_loc(control_filename) )
                        else
                                handle = setup_qmmm( num_qm_atoms, c_loc(qm_symbols), &
                                        c_loc(qm_pos_), num_mm_atoms, c_loc(mm_symbols), &
                                        c_loc(mm_pos_q_), c_loc(sim_box_info), &
                                        c_loc(ffield_filename), C_NULL_PTR )
                        endif

                        ! NVE ensemble
                        keyword = "ensemble_type" // C_NULL_CHAR
                        values = "0" // C_NULL_CHAR
                        ret = set_control_parameter( handle, c_loc(keyword), c_loc(values) )
                        if ( ret /= 0_c_int ) call sander_bomb("get_reaxff_puremd_forces", &
                                "ERROR: get_reaxff_puremd_forces::set_control_parameter", &
                                "Will exit now")

                        ! MD steps
                        keyword = "nsteps" // C_NULL_CHAR
                        values = "0" // C_NULL_CHAR
                        ret = set_control_parameter( handle, c_loc(keyword), c_loc(values) )
                        if ( ret /= 0_c_int ) call sander_bomb("get_reaxff_puremd_forces", &
                                "ERROR: get_reaxff_puremd_forces::set_control_parameter", &
                                "Will exit now")

                        ! time step length (in fs)
                        keyword = "dt" // C_NULL_CHAR
                        values = "0.25" // C_NULL_CHAR
                        ret = set_control_parameter( handle, c_loc(keyword), c_loc(values) )
                        if ( ret /= 0_c_int ) call sander_bomb("get_reaxff_puremd_forces", &
                                "ERROR: get_reaxff_puremd_forces::set_control_parameter", &
                                "Will exit now")

                        ! OpenMP number of threads
                        keyword = "num_threads" // C_NULL_CHAR
                        write (values, *) reaxff_nml%numthreads
                        values = trim(adjustl(values)) // C_NULL_CHAR
                        ret = set_control_parameter( handle, c_loc(keyword), c_loc(values) )
                        if ( ret /= 0_c_int ) call sander_bomb("get_reaxff_puremd_forces", &
                                "ERROR: get_reaxff_puremd_forces::set_control_parameter", &
                                "Will exit now")

                        ! enable periodic boundary conditions
                        keyword = "periodic_boundaries" // C_NULL_CHAR
                        values = "1" // C_NULL_CHAR
                        ret = set_control_parameter( handle, c_loc(keyword), c_loc(values) )
                        if ( ret /= 0_c_int ) call sander_bomb("get_reaxff_puremd_forces", &
                                "ERROR: get_reaxff_puremd_forces::set_control_parameter", &
                                "Will exit now")

                        ! do not remap atom coordinates within simulation box boundaries
                        keyword = "reposition_atoms" // C_NULL_CHAR
                        values = "0" // C_NULL_CHAR
                        ret = set_control_parameter( handle, c_loc(keyword), c_loc(values) )
                        if ( ret /= 0_c_int ) call sander_bomb("get_reaxff_puremd_forces", &
                                "ERROR: get_reaxff_puremd_forces::set_control_parameter", &
                                "Will exit now")

                        ! recompute Verlet neighbor lists at every (1) MD step
                        keyword = "reneighbor" // C_NULL_CHAR
                        values = "1" // C_NULL_CHAR
                        ret = set_control_parameter( handle, c_loc(keyword), c_loc(values) )
                        if ( ret /= 0_c_int ) call sander_bomb("get_reaxff_puremd_forces", &
                                "ERROR: get_reaxff_puremd_forces::set_control_parameter", &
                                "Will exit now")

                        ! disable force and energy tabulation for Coulomb interactions
                        keyword = "tabulate_long_range" // C_NULL_CHAR
                        values = "0" // C_NULL_CHAR
                        ret = set_control_parameter( handle, c_loc(keyword), c_loc(values) )
                        if ( ret /= 0_c_int ) call sander_bomb("get_reaxff_puremd_forces", &
                                "ERROR: get_reaxff_puremd_forces::set_control_parameter", &
                                "Will exit now")

                        ! calculate energies at every (1) MD step
                        keyword = "energy_update_freq" // C_NULL_CHAR
                        values = "1" // C_NULL_CHAR
                        ret = set_control_parameter( handle, c_loc(keyword), c_loc(values) )
                        if ( ret /= 0_c_int ) call sander_bomb("get_reaxff_puremd_forces", &
                                "ERROR: get_reaxff_puremd_forces::set_control_parameter", &
                                "Will exit now")

                        ! add a 2.5 Angstrom buffer to Verlet neighbor list cut-off
                        keyword = "vlist_buffer" // C_NULL_CHAR
                        values = "2.5" // C_NULL_CHAR
                        ret = set_control_parameter( handle, c_loc(keyword), c_loc(values) )
                        if ( ret /= 0_c_int ) call sander_bomb("get_reaxff_puremd_forces", &
                                "ERROR: get_reaxff_puremd_forces::set_control_parameter", &
                                "Will exit now")

                        ! 5.0 Angstrom bond interaction cut-off
                        keyword = "nbrhood_cutoff" // C_NULL_CHAR
                        write (values, *) reaxff_nml%nbrcut
                        values = trim(adjustl(values)) // C_NULL_CHAR
                        ret = set_control_parameter( handle, c_loc(keyword), c_loc(values) )
                        if ( ret /= 0_c_int ) call sander_bomb("get_reaxff_puremd_forces", &
                                "ERROR: get_reaxff_puremd_forces::set_control_parameter", &
                                "Will exit now")

                        ! 0.001 threshold for valence angle interactions
                        keyword = "thb_cutoff" // C_NULL_CHAR
                        write (values, *) reaxff_nml%thbcut
                        values = trim(adjustl(values)) // C_NULL_CHAR
                        ret = set_control_parameter( handle, c_loc(keyword), c_loc(values) )
                        if ( ret /= 0_c_int ) call sander_bomb("get_reaxff_puremd_forces", &
                                "ERROR: get_reaxff_puremd_forces::set_control_parameter", &
                                "Will exit now")

                        ! 7.5 Angstrom hydrogen bond interaction cut-off
                        keyword = "hbond_cutoff" // C_NULL_CHAR
                        write (values, *) reaxff_nml%hbondcut
                        values = trim(adjustl(values)) // C_NULL_CHAR
                        ret = set_control_parameter( handle, c_loc(keyword), c_loc(values) )
                        if ( ret /= 0_c_int ) call sander_bomb("get_reaxff_puremd_forces", &
                                "ERROR: get_reaxff_puremd_forces::set_control_parameter", &
                                "Will exit now")

                        ! 0.3 Angstrom bond graph calculation cut-off
                        keyword = "bond_graph_cutoff" // C_NULL_CHAR
                        values = "0.3" // C_NULL_CHAR
                        ret = set_control_parameter( handle, c_loc(keyword), c_loc(values) )
                        if ( ret /= 0_c_int ) call sander_bomb("get_reaxff_puremd_forces", &
                                "ERROR: get_reaxff_puremd_forces::set_control_parameter", &
                                "Will exit now")

                        ! EEM model (full system) for charge calculations
                        keyword = "charge_method" // C_NULL_CHAR
                        write (values, *) reaxff_nml%charge_method
                        values = trim(adjustl(values)) // C_NULL_CHAR
                        ret = set_control_parameter( handle, c_loc(keyword), c_loc(values) )
                        if ( ret /= 0_c_int ) call sander_bomb("get_reaxff_puremd_forces", &
                                "ERROR: get_reaxff_puremd_forces::set_control_parameter", &
                                "Will exit now")

                        ! net charge for system (in Coulombs)
                        keyword = "cm_q_net" // C_NULL_CHAR
                        write (values, *) qmcharge
                        values = trim(adjustl(values)) // C_NULL_CHAR
                        ret = set_control_parameter( handle, c_loc(keyword), c_loc(values) )
                        if ( ret /= 0_c_int ) call sander_bomb("get_reaxff_puremd_forces", &
                                "ERROR: get_reaxff_puremd_forces::set_control_parameter", &
                                "Will exit now")

                        ! GMRES algorithm in charge solver
                        keyword = "cm_solver_type" // C_NULL_CHAR
                        values = "0" // C_NULL_CHAR
                        ret = set_control_parameter( handle, c_loc(keyword), c_loc(values) )
                        if ( ret /= 0_c_int ) call sander_bomb("get_reaxff_puremd_forces", &
                                "ERROR: get_reaxff_puremd_forces::set_control_parameter", &
                                "Will exit now")

                        ! max. iterations before restarting (GMRES) in charge solver
                        keyword = "cm_solver_restart" // C_NULL_CHAR
                        values = "200" // C_NULL_CHAR
                        ret = set_control_parameter( handle, c_loc(keyword), c_loc(values) )
                        if ( ret /= 0_c_int ) call sander_bomb("get_reaxff_puremd_forces", &
                                "ERROR: get_reaxff_puremd_forces::set_control_parameter", &
                                "Will exit now")

                        ! max. iterations in charge solver
                        keyword = "cm_solver_max_iters" // C_NULL_CHAR
                        write (values, *) reaxff_nml%solvmaxit
                        values = trim(adjustl(values)) // C_NULL_CHAR
                        ret = set_control_parameter( handle, c_loc(keyword), c_loc(values) )
                        if ( ret /= 0_c_int ) call sander_bomb("get_reaxff_puremd_forces", &
                                "ERROR: get_reaxff_puremd_forces::set_control_parameter", &
                                "Will exit now")

                        ! tolerance in charge solver
                        keyword = "cm_solver_q_err" // C_NULL_CHAR
                        write (values, *) reaxff_nml%solvtol
                        values = trim(adjustl(values)) // C_NULL_CHAR
                        ret = set_control_parameter( handle, c_loc(keyword), c_loc(values) )
                        if ( ret /= 0_c_int ) call sander_bomb("get_reaxff_puremd_forces", &
                                "ERROR: get_reaxff_puremd_forces::set_control_parameter", &
                                "Will exit now")

                        ! Jacobi preconditioner in charge solver
                        keyword = "cm_solver_pre_comp_type" // C_NULL_CHAR
                        write (values, *) reaxff_nml%solvprecond
                        values = trim(adjustl(values)) // C_NULL_CHAR
                        ret = set_control_parameter( handle, c_loc(keyword), c_loc(values) )
                        if ( ret /= 0_c_int ) call sander_bomb("get_reaxff_puremd_forces", &
                                "ERROR: get_reaxff_puremd_forces::set_control_parameter", &
                                "Will exit now")

                        ! polarization energy calculation within ReaxFF
                        keyword = "include_polarization_energy" // C_NULL_CHAR
                        write (values, *) reaxff_nml%include_polar_energy
                        values = trim(adjustl(values)) // C_NULL_CHAR
                        ret = set_control_parameter( handle, c_loc(keyword), c_loc(values) )
                        if ( ret /= 0_c_int ) call sander_bomb("get_reaxff_puremd_forces", &
                                "ERROR: get_reaxff_puremd_forces::set_control_parameter", &
                                "Will exit now")

                        ! disable file I/O
                        ret = set_output_enabled( handle, 0 )
                        if ( ret /= 0_c_int ) call sander_bomb("get_reaxff_puremd_forces", &
                                "ERROR: get_reaxff_puremd_forces::set_output_enabled", &
                                "Will exit now")

                        ! set contiguous charge constraints
                        ret = set_contiguous_charge_constraints( handle, &
                                num_char_const_contig, c_loc(char_const_contig_start), &
                                c_loc(char_const_contig_end), c_loc(char_const_contig_value) )
                        if ( ret /= 0_c_int ) call sander_bomb("get_reaxff_puremd_forces", &
                                "ERROR: get_reaxff_puremd_forces::set_contiguous_charge_constraints", &
                                "Will exit now")

                        ! set custom charge constraints
                        ret = set_custom_charge_constraints( handle, &
                                num_char_const_custom, c_loc(char_const_custom_count), &
                                c_loc(char_const_custom_atom_index), &
                                c_loc(char_const_custom_coeff), &
                                c_loc(char_const_custom_rhs) )
                        if ( ret /= 0_c_int ) call sander_bomb("get_reaxff_puremd_forces", &
                                "ERROR: get_reaxff_puremd_forces::set_custom_charge_constraints", &
                                "Will exit now")

                        ret = simulate( handle )
                        if ( ret /= 0_c_int ) call sander_bomb("get_reaxff_puremd_forces", &
                                "ERROR: get_reaxff_puremd_forces::simulate", &
                                "Will exit now")

                        ret = get_atom_forces_qmmm( handle, c_loc(qm_f_), c_loc(mm_f_) )
                        if ( ret /= 0_c_int ) call sander_bomb("get_reaxff_puremd_forces", &
                                "ERROR: get_reaxff_puremd_forces::get_atom_forces_qmmm", &
                                "Will exit now")

                        ! disregard MM atom charges, as static (input-only)
                        ret = get_atom_charges_qmmm( handle, c_loc(qm_q_), C_NULL_PTR )
                        if ( ret /= 0_c_int ) call sander_bomb("get_reaxff_puremd_forces", &
                                "ERROR: get_reaxff_puremd_forces::get_atom_charges_qmmm", &
                                "Will exit now")

                        ! disregard all values except total energy
                        ret = get_system_info( handle, C_NULL_PTR, C_NULL_PTR, &
                                c_loc(e_total_), C_NULL_PTR, C_NULL_PTR, C_NULL_PTR )
                        if ( ret /= 0_c_int ) call sander_bomb("get_reaxff_puremd_forces", &
                                "ERROR: get_reaxff_puremd_forces::get_system_info", &
                                "Will exit now")
                else
                        char_const_contig_start = reaxff_nml%char_const_contig_start
                        char_const_contig_end = reaxff_nml%char_const_contig_end
                        char_const_contig_value = reaxff_nml%char_const_contig_value
                        ! count range ranges (1-based numbering of atoms)
                        num_char_const_contig = count(char_const_contig_start >= 1 .and. char_const_contig_end >= 1)
                        ! adjust constraint ends based on total num. atoms
                        do ix = 1, num_char_const_contig
                                char_const_contig_end(ix) = MIN( char_const_contig_end(ix), num_qm_atoms + num_mm_atoms )
                        enddo
                        char_const_custom_count = reaxff_nml%char_const_custom_count
                        char_const_custom_atom_index = reaxff_nml%char_const_custom_atom_index
                        char_const_custom_coeff = reaxff_nml%char_const_custom_coeff
                        char_const_custom_rhs = reaxff_nml%char_const_custom_rhs
                        ! count custom constraints
                        num_char_const_custom = count(char_const_custom_count >= 1)

                        ret = reset_qmmm( handle, num_qm_atoms, c_loc(qm_symbols), &
                                c_loc(qm_pos_), num_mm_atoms, c_loc(mm_symbols), &
                                c_loc(mm_pos_q_), c_loc(sim_box_info), &
                                C_NULL_PTR, C_NULL_PTR )
                        if ( ret /= 0_c_int ) call sander_bomb("get_reaxff_puremd_forces", &
                                "ERROR: get_reaxff_puremd_forces::reset_qmmm", &
                                "Will exit now")

                        ! net charge for system (in Coulombs)
                        keyword = "cm_q_net" // C_NULL_CHAR
                        write (values, *) qmcharge
                        values = trim(adjustl(values)) // C_NULL_CHAR
                        ret = set_control_parameter( handle, c_loc(keyword), c_loc(values) )
                        if ( ret /= 0_c_int ) call sander_bomb("get_reaxff_puremd_forces", &
                                "ERROR: get_reaxff_puremd_forces::set_control_parameter", &
                                "Will exit now")

                        ! set contiguous charge constraints
                        ret = set_contiguous_charge_constraints( handle, &
                                num_char_const_contig, c_loc(char_const_contig_start), &
                                c_loc(char_const_contig_end), c_loc(char_const_contig_value) )
                        if ( ret /= 0_c_int ) call sander_bomb("get_reaxff_puremd_forces", &
                                "ERROR: get_reaxff_puremd_forces::set_contiguous_charge_constraints", &
                                "Will exit now")

                        ! set custom charge constraints
                        ret = set_custom_charge_constraints( handle, &
                                num_char_const_custom, c_loc(char_const_custom_count), &
                                c_loc(char_const_custom_atom_index), &
                                c_loc(char_const_custom_coeff), &
                                c_loc(char_const_custom_rhs) )
                        if ( ret /= 0_c_int ) call sander_bomb("get_reaxff_puremd_forces", &
                                "ERROR: get_reaxff_puremd_forces::set_custom_charge_constraints", &
                                "Will exit now")

                        ret = simulate( handle )
                        if ( ret /= 0_c_int ) call sander_bomb("get_reaxff_puremd_forces", &
                                "ERROR: get_reaxff_puremd_forces::simulate", &
                                "Will exit now")

                        ret = get_atom_forces_qmmm( handle, c_loc(qm_f_), c_loc(mm_f_) )
                        if ( ret /= 0_c_int ) call sander_bomb("get_reaxff_puremd_forces", &
                                "ERROR: get_reaxff_puremd_forces::get_atom_forces_qmmm", &
                                "Will exit now")

                        ! disregard MM atom charges, as static (input-only)
                        ret = get_atom_charges_qmmm( handle, c_loc(qm_q_), C_NULL_PTR )
                        if ( ret /= 0_c_int ) call sander_bomb("get_reaxff_puremd_forces", &
                                "ERROR: get_reaxff_puremd_forces::get_atom_charges_qmmm", &
                                "Will exit now")

                        ! disregard all values except total energy
                        ret = get_system_info( handle, C_NULL_PTR, C_NULL_PTR, &
                                c_loc(e_total_), C_NULL_PTR, C_NULL_PTR, C_NULL_PTR )
                        if ( ret /= 0_c_int ) call sander_bomb("get_reaxff_puremd_forces", &
                                "ERROR: get_reaxff_puremd_forces::get_system_info", &
                                "Will exit now")
                end if

                ! copy double precision arrays into Amber mixed precision arrays 
                qm_q = qm_q_
                e_total = e_total_
                qm_f = qm_f_
                mm_f = mm_f_
#else
                call sander_bomb('get_reaxff_puremd_forces','reaxff-puremd is not enabled', &
                        'Check your installation or reconfigure with the -reaxff-puremd option.')
#endif
        end subroutine get_reaxff_puremd_forces


        subroutine reaxff_puremd_finalize( )
                use, intrinsic :: iso_c_binding
                implicit none

                integer (c_int) :: ret

#if defined(HAVE_REAXFF_PUREMD)
!                if ( c_associated(handle) ) then
!                        ret = cleanup( handle )
!                        if ( ret /= 0_c_int ) call sander_bomb("reaxff_puremd_finalize", &
!                                "ERROR: reaxff_puremd_finalize::cleanup", &
!                                "Will exit now")
!                endif
#else
                call sander_bomb('reaxff_puremd_finalize','reaxff-puremd is not enabled', &
                        'Check your installation or reconfigure with the -reaxff-puremd option.')
#endif
        end subroutine reaxff_puremd_finalize


        subroutine get_namelist_reaxff( reaxff_nml )
                implicit none

                type(reaxff_nml_type), intent(out) :: reaxff_nml

                ! control file name
                character(len=256) :: control
                ! force field file name
                character(len=256) :: ffield
                ! atom start number (1-based) of range for charge constraint
                integer, dimension(256) :: char_const_contig_start
                ! atom end number (1-based) of range for charge constraint
                integer, dimension(256) :: char_const_contig_end
                ! charge constraint for specified consecutive atom range
                real, dimension(256) :: char_const_contig_value
                ! counts for custom charge constraints
                integer, dimension(256) :: char_const_custom_count
                ! atom indices (1-based) for custom charge constraints
                integer, dimension(1024) :: char_const_custom_atom_index
                ! coefficients for custom charge constraints
                real, dimension(1024) :: char_const_custom_coeff
                ! right-hand size (RHS) constants for custom charge constraints
                real, dimension(256) :: char_const_custom_rhs
                ! number of threads for reaxff-puremd (OpenMP only)
                integer :: numthreads
                ! charge method
                integer :: charge_method
                ! charge solver tolerance
                real :: solvtol
                ! max number of iterations for the charge solver
                integer :: solvmaxit
                ! preconditioner type for the charge solver
                integer :: solvprecond
                ! bonded interaction cutoff
                real :: nbrcut
                ! hydrogen bonding interaction cutoff
                real :: hbondcut
                ! valence cutoff (based on bond-order)
                real :: thbcut
                ! 1 = include polarization energy in ReaxFF, 0 = otherwise
                integer :: include_polar_energy
                integer :: ierr

                namelist /reaxff/ control, ffield, char_const_contig_start, &
                        char_const_contig_end, char_const_contig_value, &
                        char_const_custom_count, char_const_custom_atom_index, &
                        char_const_custom_coeff, char_const_custom_rhs, &
                        numthreads, charge_method, solvtol, solvmaxit, solvprecond, &
                        nbrcut, hbondcut, thbcut, include_polar_energy

                ! default values
                control = '' 
                ffield = 'ffield.reaxff'
                char_const_contig_start = 0
                char_const_contig_end = 0
                char_const_contig_value = 0.0
                char_const_custom_count = 0
                char_const_custom_atom_index = 0
                char_const_custom_coeff = 0.0
                char_const_custom_rhs = 0.0
                numthreads = 1
                charge_method = 1
                solvtol = 1.0e-8
                solvmaxit = 200
                solvprecond = 1 ! Jacobi preconditioner for charge solver
                nbrcut = 5.0 ! Angstroms
                hbondcut = 7.5 ! Angstroms
                thbcut = 0.001 
                include_polar_energy = 1

                ! Read namelist
                rewind 5
                read(5, nml=reaxff, iostat=ierr)
                if ( ierr > 0 ) then
                        call sander_bomb('get_namelist_reaxff (qm2_extern_reaxff_puremd_module)', &
                        '&reaxff namelist read error', &
                        'Please check your input.')
                else if ( ierr < 0 ) then
                        write(6,'(a/a)') '&reaxff namelist read encountered end of file', &
                                'Please check your input if the calculation encounters a problem'
                end if

                reaxff_nml%control = control
                reaxff_nml%ffield = ffield
                reaxff_nml%char_const_contig_start = char_const_contig_start
                reaxff_nml%char_const_contig_end = char_const_contig_end
                reaxff_nml%char_const_contig_value = char_const_contig_value
                reaxff_nml%char_const_custom_count = char_const_custom_count
                reaxff_nml%char_const_custom_atom_index = char_const_custom_atom_index
                reaxff_nml%char_const_custom_coeff = char_const_custom_coeff
                reaxff_nml%char_const_custom_rhs = char_const_custom_rhs
                reaxff_nml%numthreads = numthreads
                reaxff_nml%charge_method = charge_method
                reaxff_nml%solvtol = solvtol
                reaxff_nml%solvmaxit = solvmaxit
                reaxff_nml%solvprecond = solvprecond
                reaxff_nml%nbrcut = nbrcut
                reaxff_nml%hbondcut = hbondcut
                reaxff_nml%thbcut = thbcut
                reaxff_nml%include_polar_energy = include_polar_energy
        end subroutine get_namelist_reaxff


        subroutine print_namelist_reaxff( reaxff_nml )
                implicit none

                type(reaxff_nml_type), intent(in) :: reaxff_nml        

                write(6,'(a)') '---------ReaxFF options-------'
                ! Print the control file name if it is set (advanced option)
                if ( len_trim(reaxff_nml%control) > 0 ) then
                        write(6,'(2a)') ' control            ', reaxff_nml%control
                endif
                write(6,'(2a)')         ' ffield               ', reaxff_nml%ffield
                write(6,'(A,I7)')       ' charge_method        ', reaxff_nml%charge_method
                write(6,"(A,E7.1)")     ' solvtol              ', reaxff_nml%solvtol
                write(6,"(A,I7)")       ' solvmaxit            ', reaxff_nml%solvmaxit
                write(6,"(A,I7)")       ' solvprecond          ', reaxff_nml%solvprecond
                write(6,"(A,F7.2)")     ' nbrcut               ', reaxff_nml%nbrcut
                write(6,"(A,F7.2)")     ' hbondcut             ', reaxff_nml%hbondcut
                write(6,"(A,F7.5)")     ' thbcut               ', reaxff_nml%thbcut
                write(6,"(A,I7)")       ' include_polar_energy ', reaxff_nml%include_polar_energy
                write(6,"(A,I7,A)")     ' numthreads           ', reaxff_nml%numthreads, '(OpenMP only)'
                write(6,'(a)')          '-----end ReaxFF options-------'
        end subroutine print_namelist_reaxff
end module qm2_extern_reaxff_puremd_module
