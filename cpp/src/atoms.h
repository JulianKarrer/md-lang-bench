#ifndef __ATOMS_H
#define __ATOMS_H

#include "types.h"

/// @brief  A container class for atoms, holding relevant information such as
/// positions, velocities, currently acting forces, masses etc. Provides
/// convenience functions for querying properties of the system such as kinetic
/// energy or atom count.
class Atoms {
  public:
    Positions_t positions;
    Velocities_t velocities;
    Forces_t forces;

    Atoms(const size_t n, double spacing);

    size_t nb_atoms() const;
};
#endif // __ATOMS_H