#include <torch/script.h>
#include <torch/torch.h>
#include <iostream>
#include <algorithm>
#include <cstdlib> // for std::getenv
#include <string>
#include "ani.hpp"
#include <future> 

extern torch::jit::script::Module model;
extern torch::Device device;

// ANI output in unit of Hartree
// Transform to kcal/mol
const double ani_unit_factor = 627.5094740631;
torch::Tensor energy_t;
torch::Tensor coordinates;

torch::Tensor grad;
torch::Tensor grad_1d;

c10::intrusive_ptr<c10::ivalue::Tuple> outputs;

// 异步执行模型推理的函数
c10::intrusive_ptr<c10::ivalue::Tuple> async_forward(const torch::Tensor& species, const torch::Tensor& coordinates, torch::jit::script::Module& model) {
    // 在异步任务中执行前向传播
    outputs = model.forward({species, coordinates}).toTuple();
}

void init_ani(int* atom_num, int* model_type)
{
    const char* model_env_path = std::getenv("MODEL_PATH");
    if (model_env_path == nullptr) {
        std::cerr << "Environment variable MODEL_PATH is not set." << std::endl;
        return;
    }

    std::string model_name;
    std::string model_path;

    // 根据 model_type 参数选择不同的模型路径
    switch (*model_type) {
        case 0:
            model_path = std::string(model_env_path) + "/ani2x_model.pt";
            model_name = "ANI-2x";
            break;
        case 2:
            model_path = std::string(model_env_path) + "/ani1_xnr.pt";
            model_name = "ANI-1xnr";
            break;
        default:
            std::cerr << "Unknown model type: " << *model_type << std::endl;
            return;
    }

    // 加载模型并赋值给全局或类成员变量 model
    try {
        model = torch::jit::load(model_path);  // 直接赋值给全局 model 变量
        model.to(torch::kCUDA);  // 如果有 GPU 支持，将模型转到 GPU
        std::cout << "Model " << model_name << " loaded successfully!\n";
    }
    catch (const c10::Error& e) {
        std::cerr << "Error loading the model: " << e.msg() << std::endl;
    }
}

void ani_forward(int* input_species, float* input_coordinates, int* atom_num) {

    /* Used for debugging
    
    // 将input_species转换为std::vector
    std::vector<int> cpp_species(input_species, input_species + *atom_num);

    // 将input_coordinates转换为std::vector
    std::vector<float> cpp_coordinates(input_coordinates, input_coordinates + *atom_num * 3);

    // 打印species
    std::cout<<"species: "<<std::endl;
    for (int i = 0; i < cpp_species.size(); ++i) {
        std::cout << cpp_species[i] << " ";
    }
    std::cout << std::endl;

    // 打印coordinates
    std::cout<<"coordinates: "<<std::endl;
    for (int i = 0; i < cpp_coordinates.size(); ++i) {
        std::cout << cpp_coordinates[i] << " ";
    }
    std::cout << std::endl;
    */

    // 创建 torch::Tensor
    torch::Tensor species = torch::from_blob(input_species, {1, *atom_num}, torch::kInt32).to(device);

    //std::cout << "species: " << species << std::endl;

    // 将坐标数据转换为 torch::Tensor，直接使用输入的坐标数据
    coordinates = torch::from_blob(input_coordinates, {1, *atom_num, 3}, torch::kFloat32).to(device).set_requires_grad(true);
    //coordinates.set_requires_grad(true);
    //std::cout << "coordinates: " << coordinates << std::endl;
    /* Type test    
    std::cout << "Total Test: " << std::endl;        
    std::cout << "species sizes: " << species.sizes() << std::endl;
    std::cout << "species dtype: " << species.dtype() << std::endl;
    std::cout << "coordinates sizes: " << coordinates.sizes() << std::endl;
    std::cout << "coordinates dtype: " << coordinates.dtype() << std::endl;
    */
    
    outputs = model.forward({species, coordinates}).toTuple();

    /*
    if (!energy_t.defined()) throw std::runtime_error("Energy tensor is undefined.");
    
    std::cout << "Pre-backward Energy (CPP): " << *output_energy << std::endl;
    */
    // if (!grad.defined()) throw std::runtime_error("Gradient tensor is undefined after backward.");
}

void ani_process(float* output_energy, double* output_gradient) {
    torch::cuda::synchronize();
    energy_t = outputs->elements()[0].toTensor();
    energy_t.backward();
    *output_energy = energy_t.item<float>()*ani_unit_factor;
    torch::Tensor grad = coordinates.grad();

    grad = grad.cpu().neg() * ani_unit_factor;

    // Efficiently copy gradient data
    std::copy(grad.data_ptr<float>(), grad.data_ptr<float>() + grad.numel(), output_gradient);
    /*
    grad = grad.to(torch::kCPU);

    // turn gradient to force
    grad = grad * -1.0;

    grad_1d = grad.view({-1});
    grad_1d *= ani_unit_factor;
    std::vector<float> grad_vector(grad_1d.data_ptr<float>(), grad_1d.data_ptr<float>() + grad_1d.numel());
    std::copy(grad_vector.begin(), grad_vector.end(), output_gradient);
    
    
    std::cout << "Energy (CPP): " << *output_energy << std::endl;
    
    
    std::cout << std::endl;
    std::cout << "Gradient(CPP): " << grad << std::endl;
    */



    // Copy the energy and gradient back to Fortran
    
    /*
    std::cout << "Energy finished!" << std::endl;
    */
    }
        
