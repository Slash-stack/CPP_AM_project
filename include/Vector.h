//
// Created by theo on 2/17/22.
//

#include <array>

#ifndef TP2_VECTOR_H
#define TP2_VECTOR_H

/**
 * Base class for a 3 dimensional vector, overriding several operators.
 */
class Vector {
public:
    Vector();

    Vector(double x, double y, double z);

    Vector(const Vector &v);

    double &operator[](int i);

    double operator[](int i) const;

    [[nodiscard]] double length() const;

    [[nodiscard]] double distance(const Vector &) const;

    Vector operator*(double s);

    Vector &operator+=(const Vector &);

    Vector &operator=(const Vector &);

    Vector operator-(const Vector &) const;

    Vector operator+(const Vector &) const;

    [[nodiscard]] Vector get_vec() const;

    void set_vec(double, double, double);

    void set_vec(const Vector &);

private:
    std::array<double, 3> vec;
};

Vector operator*(double s, Vector v);

Vector operator-(Vector v);

#endif //TP2_VECTOR_H
