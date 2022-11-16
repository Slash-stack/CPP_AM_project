#include <iostream>
#include "Particule.h"
#include <random>
#include <vector>
#include <chrono>
#include "Vector.h"
#include "Univers.h"


int main() {
    std::vector<double> l_carac = {250., 40., 0};
    Univers u = Univers(2, 2.5, l_carac);

    double distance_entre_particules = pow(2, 1 / 6);

    for (int i = 0; i < 40; i++) {
        for (int j = 0; j < 40; j++) {
            auto p = Particule(
                    {60 * distance_entre_particules + i * distance_entre_particules,
                     22 * distance_entre_particules + j * distance_entre_particules, 0},
                    {0, -10, 0},
                    {}, {}, 1., 40 * i + j, "rouge");
            u.add_particule(p);
        }
    }
    for (int i = 0; i < 160; i++) {
        for (int j = 0; j < 40; j++) {
            auto p = Particule(
                    {i * distance_entre_particules, -22 * distance_entre_particules + j * distance_entre_particules, 0},
                    {0, 0, 0},
                    {}, {}, 1., 1600 + 40 * i + j, "bleu");
            u.add_particule(p);
        }
    }
    u.simulation(0.00005, 2.5);
    return 0;
}


int main5() {
    Vector v;
    Vector v2(v);
    v += Vector(1, 1, 1);
    v2[0] = 3;
    v2[1] = 4;
    v2 = 3 * v2;
    std::cout << "vecteur2: " << v2[0] << "  vecteur2 length: " << v2.distance(v) << "\n";
    return 0;
}


int main4() {
    std::vector<Particule> particules;

    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    auto start = std::chrono::steady_clock::now();
    int iter = 1 << 20;
    for (int i = 0; i < iter; ++i) {
        Particule part{};
        part.set_position({dist(mt), dist(mt), dist(mt)});
        part.set_speed({dist(mt), dist(mt), dist(mt)});
        part.set_force({dist(mt), dist(mt), dist(mt)});
        particules.push_back(part);
    }
    for (auto p: particules) {
        p.set_position({});
    }
    std::chrono::duration<double> time = std::chrono::steady_clock::now() - start;
    std::cout << "temps: " << time.count() << "\n";
    return 0;
}



/* Main correcpondant au premier algorithme
int main4() {
    auto u = Univers();
    auto soleil = Particule({0, 0, 0}, {0, 0, 0}, {}, {}, 1., 0, "soleil");
    u.add_particule(soleil);
    auto terre = Particule({0, 1., 0}, {-1., 0, 0}, {}, {}, 3.0e-6, 1, "soleil");
    u.add_particule(terre);
    auto jupiter = Particule({0, 5.36, 0}, {-0.425, 0, 0}, {}, {}, 9.55e-4, 2, "soleil");
    u.add_particule(jupiter);
    auto haley = Particule({34.75, 0, 0}, {0, 0.0296, 0}, {}, {}, 1.e-14, 3, "soleil");
    u.add_particule(haley);
    u.simulation(0.015, 468.5, true);
}
 */