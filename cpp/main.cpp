#include <chrono>
#include <fstream>
#include <iostream>
#include "lj_direct_summation.h"
#include <verlet.h>

const double DT{0.001};
const double SIGMA{1.44};
const size_t NUMBER_OF_RUNS{1};
const size_t NUMBER_OF_TIMESTEPS{100};
const size_t NB_ATOMS_MAX{5000};
const size_t NB_ATOMS_STEP{250};

/// @brief Measure the execution time for the regular  Lennard Jones lattice
/// with direct summation for a given number of atoms and simulation steps
/// @param nb_atoms number of atoms in the lattice
/// @param direct whether to use direct summation or the truncated and shifted
/// potential
/// @return number of microseconds of execution time
int64_t run_timed(size_t nb_atoms)
{
    Atoms atoms{Atoms(nb_atoms, SIGMA * pow(2, 1. / 6.))};
    // time the execution from here
    auto start = std::chrono::high_resolution_clock::now();
    (void)lj_direct_summation(atoms, 1.0, SIGMA);
    for (size_t i{0}; i < NUMBER_OF_TIMESTEPS; i++)
    {
        verlet_step1(atoms.positions, atoms.velocities, atoms.forces, DT, 1.0);
        atoms.forces.setZero();
        (void)(lj_direct_summation(atoms, 1.0, SIGMA));
        verlet_step2(atoms.velocities, atoms.forces, DT, 1.0);
    }
    auto stop = std::chrono::high_resolution_clock::now();
    // compute the duration of the loop
    return std::chrono::duration_cast<std::chrono::microseconds>(stop - start)
        .count();
}

int main(int argc, char *argv[])
{
    // open a csv file and write the header describing the stored data
    std::ofstream file("runtimes.csv");
    file << "nb_atoms,runtime_micros" << std::endl;

    // run timed simulations
    for (size_t nb_atoms{std::max((size_t)2, NB_ATOMS_STEP)};
         nb_atoms <= NB_ATOMS_MAX; nb_atoms += NB_ATOMS_STEP)
    {
        // the core loop where simulation runs are timed:
        for (size_t run{0}; run < NUMBER_OF_RUNS; run++)
        {
            int64_t micros{run_timed(nb_atoms) / NUMBER_OF_TIMESTEPS};
            file << nb_atoms << "," << micros << std::endl;
            // output to stddout to keep track of the process
            std::cout << nb_atoms << " atoms took " << micros << " μs/iter" << std::endl;
        }
    }
    return 0;
}