#include "atoms.h"

/// @brief Query the number of atoms in the system.
/// @return the number of atoms
size_t Atoms::nb_atoms() const
{
    return positions.cols();
};

/// @brief Initialize positions on a regular lattice, given a spacing between
/// grid points. The bounding volume of the lattice should be close to
/// cube shaped and centred around the origin.
/// @param n number of atoms
/// @param spacing spacing between atoms in the lattice
void initialize_lattice(Positions_t &positions, size_t n, double spacing)
{
    size_t cube_length{(size_t)std::ceil(std::cbrt(n))};
    size_t i{0};
    for (size_t x{0}; x < cube_length + 1; x++)
    {
        for (size_t y{0}; y < cube_length + 1; y++)
        {
            for (size_t z{0}; z < cube_length + 1; z++)
            {
                if (i >= n)
                {
                    return;
                }
                // place the particles
                positions(0, i) = (x - cube_length * 0.5) * spacing;
                positions(1, i) = (y - cube_length * 0.5) * spacing;
                positions(2, i) = (z - cube_length * 0.5) * spacing;
                i++;
            }
        }
    }
}

/// @brief Initialize a set of `n` atoms with masses of `1` and zero velocities
/// on a regular lattice, given a spacing between grid points.
/// @param n number of atoms
/// @param spacing spacing between atoms in the lattice
Atoms::Atoms(const size_t n, double spacing)
    : positions{3, n},
      velocities{3, n},
      forces{3, n}
{
    positions.setZero();
    velocities.setZero();
    forces.setZero();
    initialize_lattice(positions, n, spacing);
};
