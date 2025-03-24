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

        subroutine mlp_get_energy(energy,force) bind(C, name="mlp_potential_process")
            use iso_c_binding
            real(c_float), intent(out) :: energy
            double precision, intent(out) :: force(*)
        end subroutine mlp_get_energy

    end interface

    contains

    subroutine get_ani_virial(virial, gradient, input_coordinates, atom_num)
        use iso_c_binding
        real, intent(inout) :: virial(4)
        double precision, intent(in) :: gradient(*)
        real(c_float), intent(in) :: input_coordinates(*)
        integer, intent(in) :: atom_num
        integer :: i

        virial = 0.0

        do i = 1, atom_num
            virial(1) = virial(1) + gradient(3*i-2) * input_coordinates(3*i-2)
            virial(2) = virial(2) + gradient(3*i-1) * input_coordinates(3*i-1)
            virial(3) = virial(3) + gradient(3*i) * input_coordinates(3*i)
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

    subroutine ani_get_mass_center(mass_center, input_coordinates, mass, atom_num)
        use iso_c_binding
        real, intent(out) :: mass_center(3)
        real(c_float), intent(in) :: input_coordinates(*), mass(*)
        integer, intent(in) :: atom_num
        integer :: i

        mass_center = 0.0

        do i = 1, atom_num
            mass_center(1) = mass_center(1) + input_coordinates(3*i-2)
            mass_center(2) = mass_center(2) + input_coordinates(3*i-1)
            mass_center(3) = mass_center(3) + input_coordinates(3*i)
        end do

        mass_center = mass_center / atom_num

    end subroutine ani_get_mass_center

    ! 计算限制势能和力的子程序
    subroutine compute_confining_potential(force, coordinates, mass_center, atom_num, r1, r2, k, energy)
        implicit none
        ! 输入参数
        real, intent(in) :: mass_center(3)
        real(c_float), intent(in) :: coordinates(3*atom_num)
        real, intent(in) :: r1, r2, k
        integer, intent(in) :: atom_num
        double precision, intent(inout) :: force(3*atom_num)
        
        ! 输出参数
        double precision, intent(out) :: energy

        ! 局部变量
        integer :: i
        real :: dx, dy, dz, r
        real :: ux, uy, uz
        real :: fx, fy, fz

        real :: constant_r
        double precision :: total_energy

        constant_r = r2 - r1
        total_energy = 0.0d0

        do i = 1, atom_num
            ! 计算与质心的位移向量
            dx = dble(coordinates(3*i-2)) - mass_center(1)
            dy = dble(coordinates(3*i-1)) - mass_center(2)
            dz = dble(coordinates(3*i))   - mass_center(3)

            ! 计算距离 r
            r = sqrt(dx*dx + dy*dy + dz*dz)
            if (r < r1) then
                ! 区域一：不施加额外的力和势能
                cycle
            else if (r <= r2) then

                ! 区域二：线性增长力 F = -k * (r - r1) * unit_vector
                ! 势能 V = 0.5 * k * (r - r1)^2

                ! 计算单位向量
                if (r /= 0.0) then
                    ux = dx / r
                    uy = dy / r
                    uz = dz / r
                else
                    ux = 0.0
                    uy = 0.0
                    uz = 0.0
                end if

                ! 计算力
                fx = -k * (r - r1) * ux
                fy = -k * (r - r1) * uy
                fz = -k * (r - r1) * uz

                ! 累加到总力
                force(3*i-2) = force(3*i-2) + fx
                force(3*i-1) = force(3*i-1) + fy
                force(3*i)   = force(3*i)   + fz

                ! 计算并累加势能
                total_energy = total_energy + 0.5d0 * k * (r - r1)**2

            else
                ! 区域三：恒定力 F = -k * (r2 - r1) * unit_vector
                ! 势能 V = 0.5 * k * (constant_r)^2 + k * constant_r * (r - r2)

                ! 计算单位向量
                if (r /= 0.0) then
                    ux = dx / r
                    uy = dy / r
                    uz = dz / r
                else
                    ux = 0.0
                    uy = 0.0
                    uz = 0.0
                end if

                ! 计算力
                fx = -k * constant_r * ux
                fy = -k * constant_r * uy
                fz = -k * constant_r * uz

                ! 累加到总力
                force(3*i-2) = force(3*i-2) + fx
                force(3*i-1) = force(3*i-1) + fy
                force(3*i)   = force(3*i)   + fz

                ! 计算并累加势能
                total_energy = total_energy + 0.5d0 * k * (constant_r)**2 + k * constant_r * (r - r2)
            end if
        end do

        ! 将总能量赋值给输出参数
        energy = total_energy

    end subroutine compute_confining_potential

end module call_ani