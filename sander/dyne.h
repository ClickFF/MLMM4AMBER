type :: const_e_info
   sequence
   integer :: num_states, first_atom, num_atoms, first_state, first_charge
end type const_e_info

type (const_e_info), parameter :: NULL_CE_INFO = const_e_info(0,0,0,0,0)
