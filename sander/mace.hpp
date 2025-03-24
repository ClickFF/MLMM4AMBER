#ifndef MACE_MODULE_HPP
#define MACE_MODULE_HPP

#include <torch/script.h>
#include <torch/torch.h>
#include <iostream>
#include <vector>
#include <map>
#include <stdexcept>
#include <cstdlib> // for std::getenv
#include <string>

void init_mace(int* atom_num, int* model_type);
void mace_forward(int* input_species, float* input_coordinates, int* atom_num);
void mace_process(float* output_energy, double* output_gradient);
#endif // MACE_MODULE_HPP
