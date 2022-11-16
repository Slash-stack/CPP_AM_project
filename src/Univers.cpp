//
// Created by theo on 3/24/22.
//

#include "Univers.h"
#include <iostream>

const double EPSILON = 1.0;
const double SIGMA = 1.0;

/**
 * Computes the electrostatic force of \a p2 over \a p1 and adds it to \a p1.force.
 * @param p1 Reference particle
 * @param p2 Particle from which to compute the force
 */
void compute_force(Particule &p1, Particule &p2) {
    double d = p1.distance(p2);
    auto v = Vector(p2.get_pos() - p1.get_pos());
    //p1.add_force(p1.get_m() * p2.get_m() / d / d / d * v);

    double constante = 1 / d / d * pow(SIGMA / d, 6) * (1 - 2 * pow(SIGMA / d, 6));
    p1.add_force(constante * v);
}

/**
 * Computes and updates forces of all the particles for the current state of the universe.
 */
void Univers::update_forces() {
    //pour chaque case (i, j, k)
    //  pour chaque particule de cette case
    //      pour chaque casse voisine (l, m, n)
    //          pour chaque particule de cette case voisine
    //              calculer et set_force de p1

    if (dim == 1) {
        for (int i = 0; i < size_grille[0]; ++i) {
            for (auto &p1: grille[i]) {
                p1.set_force({});
                for (int l = std::max(0, i - 1); l < std::min(i + 2, size_grille[0]); ++l) {
                    for (auto &p2: grille[l]) {
                        if (p1.get_id() != p2.get_id()) {
                            compute_force(p1, p2);
                        }
                    }
                    p1.set_force(24 * EPSILON * p1.get_force());
                }
            }
        }
    } else if (dim == 2) {
        for (int i = 0; i < size_grille[0]; ++i) {
            for (int j = 0; j < size_grille[1]; ++j) {

                for (auto &p1: grille[i + size_grille[0] * j]) {
                    p1.set_force({});

                    for (int l = std::max(0, i - 1); l < std::min(size_grille[0], i + 2); ++l) {
                        for (int m = std::max(0, j - 1); m < std::min(size_grille[1], j + 2); ++m) {

                            for (auto &p2: grille[l + size_grille[0] * m]) {
                                if (p1.get_id() != p2.get_id()) {
                                    compute_force(p1, p2);
                                }
                            }
                            p1.set_force(24 * EPSILON * p1.get_force());
                        }
                    }
                }
            }
        }
    } else if (dim == 3) {
        for (int i = 0; i < size_grille[0]; ++i) {
            for (int j = 0; j < size_grille[1]; ++j) {
                for (int k = 0; k < size_grille[2]; ++k) {

                    for (auto &p1: grille[i + size_grille[0] * j + size_grille[0] * size_grille[1] * k]) {
                        p1.set_force({});

                        for (int l = std::max(0, i - 1); l < std::min(size_grille[0], i + 2); ++l) {
                            for (int m = std::max(0, j - 1); m < std::min(size_grille[1], j + 2); ++m) {
                                for (int n = std::max(0, k - 1); n < std::min(size_grille[2], k + 2); ++n) {

                                    for (auto &p2: grille[l + size_grille[0] * m +
                                                          size_grille[0] * size_grille[1] * n]) {
                                        if (p1.get_id() != p2.get_id()) {
                                            compute_force(p1, p2);
                                        }
                                    }
                                    p1.set_force(24 * EPSILON * p1.get_force());
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

/**
 * @brief \n Simulates the evolution of the universe.
 * @details At every step \a dt of time, updates positions of particles according to Störmer-Verlet's algorithm.
 *          \n i.e. calls update_forces then update_positions and adds dt to t until t reaches t_end.
 *          \n Also writes universe's state to a .vtk file.
 * @param dt time step.
 * @param t_end end time.
 */
void Univers::simulation(double dt, double t_end) {
    t = 0;
    update_forces();
    int step = 0;
    while (t < t_end) {
        update_positions(dt);
        if (step % 10 == 0) {
            this->write_vtk(step);
        }
        step++;
        t += dt;
    }
}


/**
 * Constructor for a Universe.
 * @param dimension the dimention of the universe, can be 1, 2 or 3.
 * @param r_cut the distance from which interaction forces between particles can be neglected.
 * @param l_d vector of size 1, 2 or 3 (according to \a dimension) that describes the size of the universe along the axis.
 */
Univers::Univers(int dimension, double r_cut, std::vector<double> l_d) : nbr_part(0), t(0), dim(dimension),
                                                                         r_cut(r_cut), l_d(l_d) {

    size_grille = {0, 0, 0};
    grille = {};
    switch (dimension) {
        case 3:
            size_grille[2] = ceil(l_d[2] / r_cut);
        case 2:
            size_grille[1] = ceil(l_d[1] / r_cut);
        case 1:
            size_grille[0] = ceil(l_d[0] / r_cut);
            break;
        default:
            std::cout << "Mauvaise dimension dans le constructeur" << std::endl;
    }
    grille.resize(size_grille[0] * std::max(1, size_grille[1]) * std::max(1, size_grille[2]));
}

/**
 * Adds a particle to the universe.
 * @param part particle to add.
 */
void Univers::add_particule(const Particule &part) {
    nbr_part++;
    int x = 0, y = 0, z = 0;
    switch (dim) {
        case 3:
            z = (part.get_pos()[2] + l_d[2] / 2) / r_cut;
        case 2:
            y = (part.get_pos()[1] + l_d[1] / 2) / r_cut;
        case 1:
            x = (part.get_pos()[0] + l_d[0] / 2) / r_cut;
            break;
        default:
            std::cout << "erreur dimension ajout particule" << std::endl;
    }
    grille[z * size_grille[1] * size_grille[2] + y * size_grille[1] + x].push_back(part);
}

/**
 * Prints the position of every particles.
 */
void Univers::print_state() {
    for (auto &i: grille) {
        for (const auto &part: i) {
            std::cout << "Particule " << part.get_id() << ":" << std::endl
                      << "\tPosition " << part.get_pos()[0] << " " << part.get_pos()[1] << " " << part.get_pos()[2]
                      << std::endl;
        }
    }
}

/**
 * Saves the state of the universe as a .vtk file
 * @param step step of the simulation that we want to save
 */
void Univers::write_vtk(int step) {
    // On créée le fichier .vtk qui aura pour nom time_{step}.vtk
    std::ofstream myfile;
    std::string myfile_name = "time_" + std::to_string(step) + ".vtk";
    myfile.open(myfile_name);
    // On écrit l'entête du fichier
    myfile << "# vtk DataFile Version 2.0\n";
    myfile << myfile_name << std::endl;
    myfile << "ASCII\n";
    myfile << "DATASET UNSTRUCTURED_GRID\n";
    myfile << "POINTS " << nbr_part << " float\n";
    // On écrit les  positions des différentes particules
    for (auto &i: grille) {
        for (const auto &part: i) {
            myfile << part.get_pos()[0] << ' '
                   << part.get_pos()[1] << ' '
                   << part.get_pos()[2] << std::endl;
        }
    }
    myfile.close();
}

/**
 * Saves the state of the universe as a .vtk file (xml format)
 * @param step step of the simulation that we want to save
 * @bug
 */
void Univers::write_vtk1(int step) {
    std::ofstream myfile;
    std::string myfile_name = "time_" + std::to_string(step) + ".vtk";
    myfile.open(myfile_name);

    myfile << "<VTKFile type=\"UnstructuredGrid\" version=\"2.0\" byte_order=\"LittleEndian\">\n";
    myfile << "<UnstructuredGrid>\n";
    myfile << "<Piece NumberOfPoints=\"" << nbr_part << "\" NumberOfCells=\"0\">\n";
    myfile << "<Points>\n";
    myfile << R"(<DataArray name="Position" type="Float32" NumberOfComponents=")" << dim << "\" format=\"ascii\">\n";

    for (auto &i: grille) {
        for (const auto &part: i) {
            myfile << part.get_pos()[0];
            if (dim >= 2) {
                myfile << ' ' << part.get_pos()[1];
            }
            if (dim == 3) {
                myfile << ' ' << part.get_pos()[2];
            }
            myfile << std::endl;
        }
    }

    myfile << "</DataArray>\n";
    myfile << "</Points>\n";
    myfile << "<PointData Vectors=\"vector\">\n";
    myfile << R"(<DataArray name="Velocity" type="Float32" NumberOfComponents=")" << dim << "\" format=\"ascii\">\n";

    for (auto &i: grille) {
        for (const auto &part: i) {
            myfile << part.get_speed()[0];
            if (dim >= 2) {
                myfile << ' ' << part.get_speed()[1];
            }
            if (dim == 3) {
                myfile << ' ' << part.get_speed()[2];
            }
            myfile << std::endl;
        }
    }

    myfile << "</DataArray>\n";
    myfile << R"(<DataArray name="Masse" type="Float32" format="ascii">)" << std::endl;

    for (auto &i: grille) {
        for (const auto &part: i) {
            myfile << part.get_m() << std::endl;
        }
    }

    myfile << "</DataArray>\n";
    myfile << "</PointData>\n";

    myfile << "<Cells>\n"
              "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n"
              "</DataArray>\n"
              "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n"
              "</DataArray>\n"
              "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n"
              "</DataArray>\n"
              "</Cells>\n";

    myfile << "</Piece>\n";
    myfile << "</UnstructuredGrid>\n";
    myfile << "</VTKFile>\n";
    myfile.close();
}

/** Updates the positions of the particles in our universe
 * @brief we compute the new positions of the particles then we identify the particles that need to move (we erase these
 * particles and re-add them in the right cell.
 * @param dt
 */
void Univers::update_positions(double dt) {
    for (auto &cellule: grille) {
        for (auto &part: cellule) {
            part.set_position(part.get_pos() + dt * (part.get_force() * dt * (0.5 / part.get_m()) + part.get_speed()));
            part.set_force_old(part.get_force());
        }
    }

    std::deque<Particule> doit_bouger = {};
    int x, y, z;
    for (int i = 0; i < size_grille[0]; ++i) {
        for (int j = 0; j < size_grille[1]; ++j) {
            for (int k = 0; k < size_grille[2]; ++k) {
                auto iter = grille[i + size_grille[0] * j + size_grille[0] * size_grille[1] * k].begin();
                while (iter != grille[i + size_grille[0] * j + size_grille[0] * size_grille[1] * k].end()) {
                    auto part = *iter;
                    x = 0, y = 0, z = 0;
                    switch (dim) {
                        case 3:
                            z = (part.get_pos()[2] + l_d[2] / 2) / r_cut;
                        case 2:
                            y = (part.get_pos()[1] + l_d[1] / 2) / r_cut;
                        case 1:
                            x = (part.get_pos()[0] + l_d[0] / 2) / r_cut;
                            break;
                        default:
                            std::cout << "erreur dimension ajout particule" << std::endl;
                    }
                    if (x != i || y != j || z != k) {
                        iter = grille[i + size_grille[0] * j + size_grille[0] * size_grille[1] * k].erase(iter);
                        doit_bouger.push_back(part);
                    } else { iter++; }
                }
            }
        }
    }

    for (auto part: doit_bouger) {
        this->add_particule(part);
        nbr_part--;
    }

    update_forces();

    for (auto &cellule: grille) {
        for (auto &part: cellule) {
            part.set_speed(part.get_speed() + dt * (0.5 / part.get_m()) * (part.get_force() + part.get_force_old()));
        }
    }
}