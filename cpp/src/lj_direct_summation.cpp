/*
 * Copyright 2021 Lars Pastewka
 *
 * ### MIT license
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include <iostream>

#include <Eigen/Dense>

#include "lj_direct_summation.h"

inline std::tuple<double, double> lj_pair(double distance, double epsilon,
                                          double sigma)
{
    double sd = sigma / distance;
    double sd2 = sd * sd;
    double sd6 = sd2 * sd2 * sd2;
    double sd12 = sd6 * sd6;
    return {4 * epsilon * (sd12 - sd6),
            24 * epsilon * (2 * sd12 - sd6) / distance};
}

double lj_direct_summation(Atoms &atoms, double epsilon, double sigma)
{
    double potential_energy{0};

    for (int i{0}; i < atoms.nb_atoms(); ++i)
    {
        for (int j{i + 1}; j < atoms.nb_atoms(); ++j)
        {
            Eigen::Array3d distance_vector{atoms.positions.col(i) -
                                           atoms.positions.col(j)};
            double distance{
                std::sqrt((distance_vector * distance_vector).sum())};
            auto [pair_energy, pair_force]{lj_pair(distance, epsilon, sigma)};

            potential_energy += pair_energy;

            Eigen::Array3d force_vector{pair_force * distance_vector /
                                        distance};
            atoms.forces.col(i) += force_vector;
            atoms.forces.col(j) -= force_vector;
        }
    }

    return potential_energy;
}
