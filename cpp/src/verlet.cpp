#include "verlet.h"

/// The predictor step of a Velocity-Verlet time integration scheme
void verlet_step1(Positions_t &positions, Velocities_t &velocities,
                  const Forces_t &forces, double dt, double m)
{
    velocities += forces / m * dt * .5;
    positions += velocities * dt;
}

/// The corrector step of a velocity-Verlet time integration scheme.
/// Use this after verlet_step1 AND a subsequent update to fx,fy,fz using
/// the positions obtained through verlet_step1
void verlet_step2(Velocities_t &velocities, const Forces_t &forces, double dt,
                  double m)
{
    velocities += forces / m * dt * .5;
}
