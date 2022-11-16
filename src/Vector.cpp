//
// Created by theo on 2/17/22.
//

#include "Vector.h"
#include "cassert"
#include "cmath"


Vector::Vector(double x, double y, double z) {
    vec = {};
    vec[0] = x;
    vec[1] = y;
    vec[2] = z;
}

/**
 * Base vector constructor, initial position is (0, 0 ,0)
 */
Vector::Vector() : Vector(0, 0, 0) {}

Vector::Vector(const Vector &v) : Vector() {
    vec[0] = v.vec[0];
    vec[1] = v.vec[1];
    vec[2] = v.vec[2];
}

/**
 * Computes the length of the vector (origin to the point)
 * @return length of the vector.
 */
double Vector::length() const {
    return sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
}

double &Vector::operator[](int i) {
    assert(i >= 0 && i <= 2);
    return vec[i];
}

Vector Vector::operator*(double s) {
    return {s * vec[0], s * vec[1], s * vec[2]};
}

Vector &Vector::operator+=(const Vector &v) {
    vec[0] += v.vec[0];
    vec[1] += v.vec[1];
    vec[2] += v.vec[2];
    return *this;
}

Vector &Vector::operator=(const Vector &v) {
    if (this != &v) {
        vec[0] = v.vec[0];
        vec[1] = v.vec[1];
        vec[2] = v.vec[2];
    }
    return *this;
}

/**
 * Computes the distance between points \a this and \a v
 * @param v other vector
 * @return length of the vector \a this to \a v
 */
double Vector::distance(const Vector &v) const {
    Vector v2 = *this - v;
    return v2.length();
}

Vector Vector::operator-(const Vector &v) const {
    return Vector(this->vec[0] - v.vec[0],
                  this->vec[1] - v.vec[1],
                  this->vec[2] - v.vec[2]);
}

double Vector::operator[](int i) const {
    assert(i >= 0 && i <= 2);
    return vec[i];
}

Vector Vector::operator+(const Vector &v) const {
    return Vector(this->vec[0] + v.vec[0],
                  this->vec[1] + v.vec[1],
                  this->vec[2] + v.vec[2]);
}

Vector Vector::get_vec() const {
    return reinterpret_cast<const Vector &>(vec);
}

void Vector::set_vec(double x, double y, double z) {
    this->vec[0] = x;
    this->vec[1] = y;
    this->vec[2] = z;
}

void Vector::set_vec(const Vector &v) {
    this->set_vec(v[0], v[1], v[2]);
}


Vector operator*(double s, Vector v) {
    return {s * v[0], s * v[1], s * v[2]};
}

Vector operator-(Vector v) {
    return {-v[0], -v[1], -v[2]};
}
