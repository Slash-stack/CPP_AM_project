//
// Created by theo on 2/10/22.
//

#include "Particule.h"
#include <random>

#include <utility>

int Particule::next_id = 0;

Particule::Particule(const Vector &pos, const Vector &speed, const Vector &force, const Vector &force_old, double m,
                     int id, std::string s) {
    this->position = pos;
    this->speed = speed;
    this->force = force;
    this->force_old = force_old;
    this->categorie = std::move(s);
    this->m = m;
    this->id = id;
}

double Particule::distance(const Particule &p) const {
    return p.position.distance(this->position);
}

/**
 * Contructor for a randomly created particle between 0 and 1 in 3 directions with mass 1
 */
Particule::Particule() : Particule({}, {}, {}, {}, 1, 0, "") {
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    position[0] = dist(mt);
    position[1] = dist(mt);
    position[2] = dist(mt);
    this->id = next_id++;
}

double Particule::get_m() const {
    return m;
}

Vector Particule::get_pos() const {
    return position;
}

Vector Particule::get_force() const {
    return force;
}

Vector Particule::get_force_old() const {
    return force_old;
}

Vector Particule::get_speed() const {
    return speed;
}

int Particule::get_id() const {
    return id;
}

void Particule::set_position(const Vector &v) {
    this->position.set_vec(v);
}

void Particule::set_speed(const Vector &v) {
    this->speed.set_vec(v);
}

void Particule::set_force(const Vector &v) {
    this->force.set_vec(v);
}

void Particule::set_force_old(const Vector &v) {
    this->force_old.set_vec(v);
}

void Particule::add_force(const Vector &v) {
    this->set_force(this->force + v);
}




