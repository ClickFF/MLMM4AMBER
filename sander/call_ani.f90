module call_ani

    use iso_c_binding
    implicit none

    interface
        subroutine mlp_potential_forward(input_species, input_coordinates, atom_num) bind(C, name="mlp_potential_forward")
            use iso_c_binding
            integer(c_int), intent(in) :: input_species(*)
            real(c_float), intent(in) :: input_coordinates(*)
            integer(c_int), intent(in) :: atom_num
        end subroutine mlp_potential_forward

        subroutine init_ani(atom_num, mlp_model) bind(C, name="init_mlp_potential")
            use iso_c_binding
            integer(c_int), intent(in) :: atom_num
            integer(c_int), intent(in) :: mlp_model
        end subroutine init_ani

        subroutine ani_synchronize() bind(C, name="mlp_potential_synchronize")
            use iso_c_binding
        end subroutine ani_synchronize

        subroutine mlp_get_energy(energy,force) bind(C, name="mlp_potential_process")
            use iso_c_binding
            real(c_float), intent(out) :: energy
            double precision, intent(out) :: force(*)
        end subroutine mlp_get_energy

    end interface

    contains

    subroutine get_ani_virial(virial, force, input_coordinates, atom_num)
        use iso_c_binding
        real, intent(inout) :: virial(4)
        double precision, intent(in) :: force(*)
        real(c_float), intent(in) :: input_coordinates(*)
        integer, intent(in) :: atom_num
        integer :: i

        virial = 0.0

        do i = 1, atom_num
            virial(1) = virial(1) + force(3*i-2) * input_coordinates(3*i-2)
            virial(2) = virial(2) + force(3*i-1) * input_coordinates(3*i-1)
            virial(3) = virial(3) + force(3*i) * input_coordinates(3*i)
        end do

    end subroutine get_ani_virial

    subroutine ani_zero_grad(grad, natom)
        use iso_c_binding
        double precision, intent(inout) :: grad(*)
        integer, intent(in) :: natom
        double precision :: sum_x, sum_y, sum_z
        integer :: i

        ! 初始化总和变量
        sum_x = 0.0d0
        sum_y = 0.0d0
        sum_z = 0.0d0

        ! 计算x, y, z三个方向上的总和
        do i = 1, natom
            sum_x = sum_x + grad(3*i - 2)
            sum_y = sum_y + grad(3*i - 1)
            sum_z = sum_z + grad(3*i)
        end do

        grad(1) = grad(1) - sum_x
        grad(2) = grad(2) - sum_y
        grad(3) = grad(3) - sum_z


    end subroutine ani_zero_grad

end module call_ani