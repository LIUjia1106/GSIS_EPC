# GSIS for e-p BTEs

## Introduction
The project utilize GSIS solve coupled electron-phonon Boltzmann equation, implemented with DG.

## Usage
### Build Requirement
- Dependency:
    - [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) (You don't have to install it since it has been involved in _'./includes'_)
    - [OpenMPI](https://www.open-mpi.org/)
- Compiler: GCC 11.4.0
### Build & Run
Please fork or download the project from the **main** branch
- To configure: _(in the root directory)_

    `cmake -S . -B build`

- To build: _(in the root directory)_

    `cmake --build build`

- To run: _(in the root directory)_

    `./AGEFEM_Poros`

### Personalization
- The change of the parameters of the phonon transport model (such as Kndusen numerbs, etc.), boundary conditions, and the parameters of the numerical computation can be conducted in _'./src/AGEFEM_Poros.cpp'_;
- The change of the computation domain (as well as the mesh) requires the modification to the mesh file (in _'./Mesh'_), _SpatialMesh_ and _'./src/AGEFEM_Poros.cpp'_.
