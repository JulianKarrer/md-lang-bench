#ifndef __VERLET_H
#define __VERLET_H

#include <Eigen/Dense>
#include <types.h>

void verlet_step1(Positions_t &positions, Velocities_t &velocities,
                  const Forces_t &forces, double dt, double m);
void verlet_step2(Velocities_t &velocities, const Forces_t &forces, double dt,
                  double m);

#endif // __VERLET_H
