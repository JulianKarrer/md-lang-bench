use nalgebra::Vector3;
use std::{
    fs::{File, OpenOptions},
    io::{BufWriter, Write},
    time::Instant,
};

// Define types
/// Type alias for typechecking positions
type PositionT = Vec<Vector3<f64>>;
/// Type alias for typechecking velocities
type VelocityT = Vec<Vector3<f64>>;
/// Type alias for typechecking forces
type ForceT = Vec<Vector3<f64>>;

/// Structure that holds buffers for the positions, velocities, forces, masses etc. of all atoms.
/// We use a structure of arrays instead of an array of structs for performance reasons.
struct Atoms {
    pos: PositionT,
    vel: VelocityT,
    frc: ForceT,
}

impl Atoms {
    /// Initialize a set of n grid points with given spacing on a regular lattice
    /// around the origin that should be as cubic as possible
    fn initialize_lattice(nb_atoms: usize, spacing: f64) -> PositionT {
        let cube_length: usize = (nb_atoms as f64).cbrt().ceil() as usize;
        let mut res: PositionT = vec![Vector3::new(0.0, 0.0, 0.0); nb_atoms];
        let mut i = 0;
        for x in 0..cube_length + 1 {
            for y in 0..cube_length + 1 {
                for z in 0..cube_length + 1 {
                    if i >= nb_atoms {
                        break;
                    }
                    res[i] = Vector3::new(
                        (x as f64 - cube_length as f64 * 0.5) * spacing,
                        (y as f64 - cube_length as f64 * 0.5) * spacing,
                        (z as f64 - cube_length as f64 * 0.5) * spacing,
                    );
                    i += 1;
                }
            }
        }
        res
    }

    /// Initialize a set of `nb_atoms` on an oblique lattice with unit mass and no velocities or forces
    pub fn new(nb_atoms: usize, sigma: f64) -> Self {
        let pos = Atoms::initialize_lattice(nb_atoms, sigma * 2.0f64.powf(1. / 6.));
        return Self {
            pos,
            vel: vec![Vector3::zeros(); nb_atoms],
            frc: vec![Vector3::zeros(); nb_atoms],
        };
    }

    /// Compute the full lennard-Jones potential by iterating over all pairs of atoms in
    /// O(n^2)-time
    pub fn lj_direct_sum(&mut self, epsilon: f64, sigma: f64) -> f64 {
        self.frc.fill(Vector3::zeros());
        let n = self.pos.len();
        let mut potential = 0.0;
        for i in 0..n {
            for j in i + 1..n {
                let dist_vec = self.pos[i] - self.pos[j];
                let distance = dist_vec.norm();
                let (energy, force_mag) = lj_pair(distance, epsilon, sigma);
                potential += energy;
                let f = force_mag * dist_vec / distance;
                self.frc[i] += f;
                self.frc[j] -= f;
            }
        }
        return potential;
    }
}

/// The Lennard-Jones potential: returns the energy and force as a tuple
fn lj_pair(distance: f64, epsilon: f64, sigma: f64) -> (f64, f64) {
    let sd = sigma / distance;
    let sd6 = sd.powi(6);
    let sd12 = sd6 * sd6;
    (
        4.0 * epsilon * (sd12 - sd6),
        24.0 * epsilon * (2.0 * sd12 - sd6) / distance,
    )
}

/// Predictor step of the Velocity Verlet time integration scheme
fn velocity_verlet_step1(
    pos: &mut PositionT,
    vel: &mut VelocityT,
    frc: &ForceT,
    mas: f64,
    dt: f64,
) {
    vel.iter_mut().zip(pos).zip(frc).for_each(|((v, x), f)| {
        *v += f / mas * dt * 0.5;
        *x += *v * dt;
    });
}

/// Corrector step of the Velocity Verlet time integration scheme
fn velocity_verlet_step2(vel: &mut VelocityT, frc: &ForceT, mas: f64, dt: f64) {
    vel.iter_mut().zip(frc).for_each(|(v, f)| {
        *v += *f / mas * dt * 0.5;
    });
}

// define parameters
const DT: f64 = 0.001;
const SIGMA: f64 = 1.44;
const EPSILON: f64 = 1.;

const NB_RUNS: usize = 1;
const NB_TIMESTEPS: usize = 100;
const NB_ATOMS_MAX: usize = 5000;
const NB_ATOMS_STEP: usize = 250;

/// Run the simulation for `nb_atoms` many atoms, using LJDS or LJTS, and return
/// the runtime in microseconds.
fn run_timed(nb_atoms: usize) -> u128 {
    let mut atoms = Atoms::new(nb_atoms, SIGMA);

    // time the execution of a simulation as seen in milestone 5
    let start = Instant::now();
    atoms.lj_direct_sum(EPSILON, SIGMA);
    // main simulation loop
    for _ in 0..NB_TIMESTEPS {
        // compute lj forces and integrate with velocity verlet
        velocity_verlet_step1(&mut atoms.pos, &mut atoms.vel, &atoms.frc, 1.0, DT);
        atoms.lj_direct_sum(EPSILON, SIGMA);
        velocity_verlet_step2(&mut atoms.vel, &atoms.frc, 1.0, DT);
    }
    start.elapsed().as_micros()
}

fn main() {
    // open csv file and write header
    let mut file: BufWriter<File> = BufWriter::new(
        OpenOptions::new()
            .write(true)
            .create(true)
            .truncate(true)
            .open("runtimes.csv")
            .unwrap(),
    );
    write!(file, "nb_atoms,runtime_micros\n").unwrap();

    for n in (2.max(NB_ATOMS_STEP)..=NB_ATOMS_MAX).step_by(NB_ATOMS_STEP) {
        for _ in 0..NB_RUNS {
            let micros = run_timed(n) / NB_TIMESTEPS as u128;
            write!(file, "{},{}\n", n, micros).unwrap();
            println!("{} atoms took {} μs", n, micros);
        }
    }
}
