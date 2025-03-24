module ani_module
    use findmask
    use lmod_driver
    use iso_c_binding

    implicit none
    private
    public :: ani_nml, ani_struct, init_isani, init_ani_coord_species

    integer, dimension(:), pointer :: isani
    integer :: ani_error

    type ani_nml_type
        logical :: ifmlp
        integer :: anitype
        character(len=256) :: animask
        integer :: ani_natom
        integer :: mlp_model
        integer :: mlp_shake
        integer :: mlp_confine
        real :: confine1, confine2, confinek
    end type ani_nml_type

    type ani_struct_type
        real(c_float), allocatable, dimension(:) :: ani_coordinates  
        integer(c_int), allocatable, dimension(:) :: ani_species  
        double precision, allocatable, dimension(:) :: ani_gradient
        real, allocatable, dimension(:) :: ani_mass
        real, allocatable :: ani_link_atom(:)
        real :: ani_reference_energy
        real :: ani_energy
        integer, dimension(:), pointer :: ani_isani
        integer, allocatable, dimension(:) :: ani_mask_index
        real :: ani_virial(4)
        real :: ani_mom_p(3)
        real :: ani_ekrot, ani_ekcm
        real :: ani_mass_center(3)
        double precision :: confine_energy

        
    end type ani_struct_type

    type(ani_nml_type), save :: ani_nml
    type(ani_struct_type), save :: ani_struct

contains

    subroutine init_isani(natom)
        integer, intent(in) :: natom
        
        allocate(isani( natom ), stat=ani_error)
        if (ani_error /= 0) then
            print *, "Error in allocation"
            stop
        endif
        ani_struct%ani_isani => isani
    end subroutine init_isani

    subroutine init_ani_coord_species(ani_atom)
        integer, intent(in) :: ani_atom
        allocate(ani_struct%ani_coordinates(ani_atom*3), stat=ani_error)
        allocate(ani_struct%ani_species(ani_atom), stat=ani_error)
        allocate(ani_struct%ani_gradient(ani_atom*3), stat=ani_error)  
        allocate(ani_struct%ani_mask_index(ani_atom), stat=ani_error)
        allocate(ani_struct%ani_mass(ani_atom), stat=ani_error)
        if (ani_error /= 0) then
            print *, "Error in allocation"
            stop
        endif
    end subroutine init_ani_coord_species

    
end module ani_module

module ml_shared
  double precision, dimension(:), allocatable :: ml_force
  integer, dimension(:), allocatable :: ani_mask_index
end module ml_shared