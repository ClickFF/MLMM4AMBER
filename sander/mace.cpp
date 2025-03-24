#include <torch/script.h>
#include <torch/torch.h>
#include <iostream>
#include <vector>
#include <map>
#include <stdexcept>
#include <cstdlib> // for std::getenv
#include <string>
#include "mace.hpp"

extern torch::jit::script::Module model;
extern torch::Device device;

// MACE output in units of eV
// Transform to kcal/mol
const double mace_unit_factor = 23.060541945329;

std::vector<int> AtomicNumberTable = {1, 6, 7, 8, 9, 15, 16, 17, 35, 53};
c10::IValue output;

void init_mace(int* atom_num, int* model_type) {
    const char* model_env_path = std::getenv("MODEL_PATH");
    
    if (model_env_path == nullptr) {
        std::cerr << "Environment variable MODEL_PATH is not set." << std::endl;
        return;
    }
    
    std::string model_name;
    std::string model_path;

    // 根据 model_type 参数选择不同的模型路径
    switch (*model_type) {
        case 1:
            model_path = std::string(model_env_path) + "/MACE-OFF23_small.pt";
            model_name = "MACE-OFF23(S)";
            break;
        case 3:
            model_path = std::string(model_env_path) + "/MACE-OFF23_medium.pt";
            model_name = "MACE-OFF23(M)";
            break;
        case 4:
            model_path = std::string(model_env_path) + "/MACE-OFF23_large.pt";
            model_name = "MACE-OFF23(L)";
            break;
        default:
            std::cerr << "Unknown model type: " << *model_type << std::endl;
            return;
    }

    // 尝试加载模型并打印信息
    try {
        model = torch::jit::load(model_path);  // 直接赋值给全局 model 变量
        model.to(device);  // 将模型转移到设备上（如 GPU）
        std::cout << "Model " << model_name << " loaded successfully!\n";
    }
    catch (const c10::Error& e) {
        std::cerr << "Error loading the model: " << e.msg() << std::endl;
    }
}



void calculate_edge_index_and_shifts(int* atom_num, torch::Tensor atom_coords, torch::Tensor& edge_index, torch::Tensor& shifts, double r_max_squared) {
    int n_nodes = *atom_num;
    int n_edges = 0;
    
    // 第一步：计算边的数量
    std::vector<int> n_edges_vec(n_nodes, 0);
    #pragma omp parallel for reduction(+:n_edges)
    for (int ii = 0; ii < n_nodes; ++ii) {
        for (int jj = 0; jj < n_nodes; ++jj) {
            if (ii != jj) {
                double delx = atom_coords[ii][0].item<double>() - atom_coords[jj][0].item<double>();
                double dely = atom_coords[ii][1].item<double>() - atom_coords[jj][1].item<double>();
                double delz = atom_coords[ii][2].item<double>() - atom_coords[jj][2].item<double>();
                double rsq = delx * delx + dely * dely + delz * delz;

                if (rsq < r_max_squared) {
                    n_edges_vec[ii]++;
                    n_edges++;
                }
            }
        }
    }

    // 第二步：构造 first_edge 数组
    std::vector<int> first_edge(n_nodes);
    first_edge[0] = 0;
    for (int ii = 0; ii < n_nodes - 1; ++ii) {
        first_edge[ii + 1] = first_edge[ii] + n_edges_vec[ii];
    }

    // 第三步：初始化 edge_index 和 shifts（全为0）
    edge_index = torch::empty({2, n_edges}, torch::dtype(torch::kInt64));
    shifts = torch::zeros({n_edges, 3}, torch::dtype(torch::kFloat64));

    // 填充 edge_index
    #pragma omp parallel for
    for (int ii = 0; ii < n_nodes; ++ii) {
        int k = first_edge[ii];
        for (int jj = 0; jj < n_nodes; ++jj) {
            if (ii != jj) {
                double delx = atom_coords[ii][0].item<double>() - atom_coords[jj][0].item<double>();
                double dely = atom_coords[ii][1].item<double>() - atom_coords[jj][1].item<double>();
                double delz = atom_coords[ii][2].item<double>() - atom_coords[jj][2].item<double>();
                double rsq = delx * delx + dely * dely + delz * delz;

                if (rsq < r_max_squared) {
                    edge_index[0][k] = ii;
                    edge_index[1][k] = jj;
                    k++;
                }
            }
        }
    }
}

// 主函数：调用两次计算
void mace_forward(int* input_species, float* input_coordinates, int* atom_num) {
    torch::Tensor species = torch::from_blob(input_species, {*atom_num}, torch::kInt64);
    // 创建一个空的 coordinates tensor
    torch::Tensor coordinates = torch::empty({*atom_num, 3}, torch::kFloat32);
    
    // 使用 OpenMP 并行将 input_coordinates 复制到 coordinates tensor 中
    #pragma omp parallel for
    for (int ii = 0; ii < *atom_num; ++ii) {
        coordinates[ii][0] = input_coordinates[ii * 3 + 0];  // 复制 x 坐标
        coordinates[ii][1] = input_coordinates[ii * 3 + 1];  // 复制 y 坐标
        coordinates[ii][2] = input_coordinates[ii * 3 + 2];  // 复制 z 坐标
    }

    // 定义 r_max_squared，用于判断原子之间的距离
    double r_max_squared = 25.0; // 举例，3Å的最大半径
    
    // 定义 edge_index 和 shifts tensors
    torch::Tensor edge_index, shifts;

    // 调用 calculate_edge_index_and_shifts 函数
    calculate_edge_index_and_shifts(atom_num, coordinates, edge_index, shifts, r_max_squared);

    // 将计算得到的 edge_index 和 shifts（全为0）传递到设备上
    edge_index = edge_index.to(device);
    shifts = shifts.to(device);
    coordinates = coordinates.to(device);

    // 如果 node_attrs 需要更新
    auto node_attrs = torch::zeros({*atom_num, static_cast<int32_t>(AtomicNumberTable.size())}, torch::dtype(torch::kFloat64));
    torch::Tensor updated_node_attrs = node_attrs.clone(); 
    for (int i = 0; i < *atom_num; ++i) {
        int atomic_number = input_species[i];
        auto it = std::find(AtomicNumberTable.begin(), AtomicNumberTable.end(), atomic_number);
        if (it != AtomicNumberTable.end()) {
            int idx = it - AtomicNumberTable.begin();
            updated_node_attrs[i][idx] = 1.0;
        } else {
            throw std::runtime_error("Unknown atomic number");
        }
    }
    node_attrs = updated_node_attrs;
    node_attrs.to(device);

    // 初始化Tensor
    auto batch = torch::zeros({*atom_num}, torch::dtype(torch::kInt64));
    auto cell = torch::zeros({3, 3}, torch::dtype(torch::kFloat64));
    auto energy = torch::tensor({0.0}, torch::dtype(torch::kFloat64));
    auto force = torch::zeros({*atom_num, 3}, torch::dtype(torch::kFloat64));
    auto ptr = torch::empty({2}, torch::dtype(torch::kInt64));
    auto unit_shifts = torch::zeros({*atom_num * (*atom_num - 1), 3}, torch::dtype(torch::kFloat64)); // 容量假设
    auto weight = torch::tensor({1.0}, torch::dtype(torch::kFloat64));
    auto local_or_ghost = torch::ones({1}, torch::kFloat64);
    ptr[0] = 0;
    ptr[1] = *atom_num;

    batch = batch.to(device);
    cell = cell.to(device);
    energy = energy.to(device);
    force = force.to(device);
    node_attrs = node_attrs.to(device);
    ptr = ptr.to(device);
    unit_shifts = unit_shifts.to(device);
    weight = weight.to(device);
    local_or_ghost = local_or_ghost.to(device);

    // 重新准备 data_dict
    c10::Dict<std::string, torch::Tensor> data_dict;

    data_dict.insert("node_attrs", node_attrs);
    data_dict.insert("batch", batch);
    data_dict.insert("cell", cell);
    data_dict.insert("energy", energy);
    data_dict.insert("forces", force);
    data_dict.insert("ptr", ptr);
    data_dict.insert("unit_shifts", unit_shifts);
    data_dict.insert("weight", weight);

    data_dict.insert("positions", coordinates);
    data_dict.insert("edge_index", edge_index.to(device));
    data_dict.insert("shifts", shifts.to(device));

    //std::cout<<data_dict.at("positions").sizes()<<std::endl;
    //std::cout<<data_dict.at("edge_index").sizes()<<std::endl;
    //std::cout<<data_dict.at("shifts").sizes()<<std::endl;


    
    //std::cout<<data_dict.at("positions")<<std::endl;

    //std::cout<<"Before forward"<<std::endl;
    // 执行 forward 操作
    output = model.forward({data_dict, local_or_ghost});
    //std::cout<<"After forward"<<std::endl;
}

void mace_process(float* output_energy, double* output_force) {
    //torch::cuda::synchronize();
    auto output_dict = output.toGenericDict();
    //std::cout<<output_dict.at("energy").toTensor().cpu()<<std::endl;
    auto out_forces = output_dict.at("forces").toTensor().cpu().to(torch::kDouble);
    std::copy(out_forces.data_ptr<double>(), out_forces.data_ptr<double>() + out_forces.numel(), output_force);
    for (int i = 0; i < out_forces.numel(); ++i) {
        output_force[i] *= mace_unit_factor;  // forces转换为kcal/mol
    }
    auto out_energy = output_dict.at("energy").toTensor().cpu();
    //std::cout << out_energy << std::endl;
    *output_energy = out_energy.item<float>() * mace_unit_factor;  // energy转换为kcal/mol
}
