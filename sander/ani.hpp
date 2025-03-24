#ifndef ANI_MODULE_HPP
#define ANI_MODULE_HPP

#include <torch/script.h>
#include <torch/torch.h>
#include <iostream>
#include <algorithm>
#include <cstdlib> // for std::getenv
#include <string>

void init_ani(int* atom_num, int* model_type);
void ani_forward(int* input_species, float* input_coordinates, int* atom_num);
void ani_process(float* output_energy, double* output_gradient);
#endif // ANI_MODULE_HPP
