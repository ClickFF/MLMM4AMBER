type :: const_ph_info
   sequence
   integer :: num_states, first_atom, num_atoms, first_state, first_charge
end type const_ph_info

type (const_ph_info), parameter :: NULL_CPH_INFO = const_ph_info(0,0,0,0,0)
