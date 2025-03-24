#ifndef MLP_POTENTIAL_HPP
#define MLP_POTENTIAL_HPP

#include <vector>
#include <string>

extern "C" {
    void init_mlp_potential(int* atom_num, int* mlp_model);
    void mlp_potential_forward(int* input_species, float* input_coordinates, int* atom_num);
    void mlp_potential_process(float* output_energy, double* output_force);
}

typedef void (*InitFunc)(int*, int*);
typedef void (*ForwardFunc)(int*, float*, int*);
typedef void (*ProcessFunc)(float*, double*);

extern InitFunc initFunc;
extern ForwardFunc forwardFunc;
extern ProcessFunc processFunc;

#endif // MLP_POTENTIAL_HPP
