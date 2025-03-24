# Hybrid ML/MM Molecular Dynamics Engine in AMBER

This repository contains the implementation of a hybrid Machine Learning / Molecular Mechanics (ML/MM) molecular dynamics simulation engine, based on the AMBER 2023 software suite. This tool enables multiscale simulations by integrating machine-learned interatomic potentials (MLIPs) into the AMBER MD engine, specifically SANDER, to perform accurate and efficient free energy calculations.

## 🚀 Features

- ⚙️ Built into AMBER's SANDER engine with Fortran/C++ interfacing
- 🤖 Integrated support for popular MLIPs:
  - ANI series (ANI-1x, ANI-1ccx, ANI-2x)
  - MACE series (MACE-OFF23(S), MACE-OFF23(M), MACE-OFF23(L))
- 🔄 Supports inference + auto-gradient in C++ for force calculation
- 🧠 ML on GPU, MD on CPU for parallel asynchronous computing
- 📈 Validated with energy & momentum conservation tests
- 💧 Supports thermodynamic integration (TI) for hydration free energy calculations
- 🔬 Compatible with MM-PBSA for protein-ligand binding free energy estimations

## 📚 Reference

This work is based on the methods described in the paper:

**"Accurate Free Energy Calculation via Multiscale Simulations Driven by Hybrid Machine Learning and Molecular Mechanics Potentials"**

> Xujian Wang, Xiongwu Wu, Bernard R. Brooks, and Junmei Wang  
> ChemRxiv, 2024, DOI: [10.26434/chemrxiv-2024-zq975-v2](https://doi.org/10.26434/chemrxiv-2024-zq975-v2)

## 📬 Contact

For questions or collaboration, please contact:

- hsuchein0126@outlook.com  
- juw79@pitt.edu

---

© 2025. Released under [CC BY 4.0 License](https://creativecommons.org/licenses/by/4.0/)
