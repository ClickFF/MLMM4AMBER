#include "mlp_potential.hpp"
#include "ani.hpp"
#include "mace.hpp"

torch::jit::script::Module model;
auto device = torch::Device(torch::kCUDA, 0);

InitFunc initFunc = nullptr;
ForwardFunc forwardFunc = nullptr;
ProcessFunc processFunc = nullptr;

extern "C" {
    void init_mlp_potential(int* atom_num, int* mlp_model) {
        try {
            // 选择模型
            if (*mlp_model == 0 or *mlp_model==2) {
                initFunc = init_ani;
                forwardFunc = ani_forward;
                processFunc = ani_process; // 初始化 processFunc
            } else if (*mlp_model == 1 or *mlp_model==3 or *mlp_model==4) {
                initFunc = init_mace;
                forwardFunc = mace_forward;
                processFunc = mace_process; // 初始化 processFunc
            } else {
                throw std::invalid_argument("Unknown model type");
            }
        } catch (const std::invalid_argument& e) {
            std::cerr << "Error: " << e.what() << std::endl;
            return;
        }

        try {
            // 初始化模型
            if (initFunc) {
                initFunc(atom_num, mlp_model);
            }
        } catch (const std::exception& e) {
            std::cerr << "Initialization error: " << e.what() << std::endl;
            return;
        }
    }

    void mlp_potential_forward(int* input_species, float* input_coordinates, int* atom_num) {
        if (forwardFunc) {
            forwardFunc(input_species, input_coordinates, atom_num);
        } else {
            std::cerr << "Error: forwardFunc is not initialized." << std::endl;
        }
    }

    void mlp_potential_process(float* output_energy, double* output_force) {
        if (processFunc) {
            processFunc(output_energy, output_force);
        } else {
            std::cerr << "Error: processFunc is not initialized." << std::endl;
        }
    }
}
