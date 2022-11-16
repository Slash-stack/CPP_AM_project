//
// Created by theo on 3/24/22.
//

#ifndef TP2_UNIVERS_H
#define TP2_UNIVERS_H

#include <deque>
#include <vector>
#include <cmath>
#include "Particule.h"
#include <iostream>
#include <fstream>

/**
 * Class for a 1, 2 or 3-dimensional universe of \a Particles.
 * Supports simulation according to Strömer-Verlet's algorithm for electrostatic forces and output to .vtk file for visualisation.
 *
 * The particles are stored inside deques into a grid (1, 2 or 3 D) this grid is strored in a 1D vector.
 * To factorise code as much as possible we had to store the x axis into the last dimension of the grid so that for example grid[y * grid_size[1] + x]
 * gets the deque (x, y) in 2D and the deque x in 1D when grid_size[1] is 0
 */
class Univers {
public:
    Univers(int dimension, double r_cut, std::vector<double> l_d);

    void print_state();

    void simulation(double dt, double t_end);

    void update_positions(double dt);

    void add_particule(const Particule &part);

    void write_vtk(int);

    void write_vtk1(int);

private:
    void update_forces();

    std::vector<int> size_grille;
    std::vector<std::deque<Particule>> grille;
    int nbr_part;
    double t;
    int dim;
    double r_cut;
    std::vector<double> l_d;

};

#endif //TP2_UNIVERS_H
