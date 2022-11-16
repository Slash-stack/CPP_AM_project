//
// Created by theo on 2/10/22.
//

#ifndef TP2_PARTICULE_H
#define TP2_PARTICULE_H

#include <array>
#include "Vector.h"
#include <string>

/**
 * Class for a 3-dimentional particle.
 */
class Particule {
private:
    Vector position;
    Vector speed;
    Vector force;
    Vector force_old;
    std::string categorie;
    double m;
    int id;

    static int next_id;
public:
    Particule();

    Particule(const Vector &pos, const Vector &speed, const Vector &force, const Vector &force_old, double m, int id,
              std::string s);

    [[nodiscard]] double distance(const Particule &) const;

    [[nodiscard]] double get_m() const;

    [[nodiscard]] Vector get_pos() const;

    [[nodiscard]] Vector get_force() const;

    [[nodiscard]] Vector get_force_old() const;

    [[nodiscard]] Vector get_speed() const;

    [[nodiscard]] int get_id() const;

    void set_position(const Vector &);

    void set_speed(const Vector &);

    void set_force(const Vector &);

    void add_force(const Vector &);

    void set_force_old(const Vector &);

};


#endif //TP2_PARTICULE_H
