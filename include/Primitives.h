/*
 * Primitives.h
 *
 *  Created on: 29.04.2016
 *  Author: veske & scratchapixel
 *  http://www.scratchapixel.com/code.php?id=9&origin=/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle&src=1
 */

#ifndef PRIMITIVES_H_
#define PRIMITIVES_H_

#include <math.h>
#include <deal.II/numerics/vector_tools.h>

/** Template class to define the 3D vector with its operations */
template<typename T>
class Vec3 {
public:
    /** Vec3 constructors */
    Vec3() : x(T(0)), y(T(0)), z(T(0)) {}
    Vec3(T xx) : x(xx), y(xx), z(xx) {}
    Vec3(T xx, T yy, T zz) : x(xx), y(yy), z(zz) {}

    /** Addition of two vectors */
    Vec3 operator +(const Vec3 &v) const {
        return Vec3(x + v.x, y + v.y, z + v.z);
    }

    /** Subtraction of two vectors */
    Vec3 operator -(const Vec3 &v) const {
        return Vec3(x - v.x, y - v.y, z - v.z);
    }
    Vec3 operator -() const {
        return Vec3(-x, -y, -z);
    }

    /** Scalar multiplication of vector with a scalar and with another vector */
    Vec3 operator *(const T &r) const {
        return Vec3(x * r, y * r, z * r);
    }
    Vec3 operator *(const Vec3 &v) const {
        return Vec3(x * v.x, y * v.y, z * v.z);
    }
    Vec3& operator *=(const T &r) {
        x *= r, y *= r, z *= r;
        return *this;
    }
    friend Vec3 operator *(const T &r, const Vec3 &v) {
        return Vec3<T>(v.x * r, v.y * r, v.z * r);
    }

    /** Equals operator */
    bool operator ==(const Vec3 &v) const {
        return x == v.x && y == v.y && z == v.z;
    }

    /** Scalar division of vector with a scalar or with another vector */
    Vec3 operator /(const T &r) const {
        return Vec3(x / r, y / r, z / r);
    }
    Vec3& operator /=(const T &r) {
        x /= r, y /= r, z /= r;
        return *this;
    }
    friend Vec3 operator /(const T &r, const Vec3 &v) {
        return Vec3<T>(r / v.x, r / v.y, r / v.z);
    }

    /** Dot product, cross product, norm and length */
    T dotProduct(const Vec3<T> &v) const {
        return x * v.x + y * v.y + z * v.z;
    }
    Vec3 crossProduct(const Vec3<T> &v) const {
        return Vec3<T>(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x);
    }
    T norm() const {
        return x * x + y * y + z * z;
    }
    T length() const {
        return sqrt(norm());
    }

    /** Function to normalize the vector */
    Vec3& normalize() {
        T n = norm();
        if (n > 0) {
            T factor = 1 / sqrt(n);
            x *= factor, y *= factor, z *= factor;
        }
        return *this;
    }

    /**
     Define access operators or accessors.
     The Vec3 coordinates can be accessed that way v[0], v[1], v[2], rather than v.x, v.y, v.z.
     This is useful in loops: the coordinates can be accessed with the loop index (e.g. v[i]).
     */
    const T& operator [](uint8_t i) const {
        return (&x)[i];
    }
    T& operator [](uint8_t i) {
        return (&x)[i];
    }

    /** Defining the behaviour of cout */
    friend std::ostream& operator <<(std::ostream &s, const Vec3<T> &v) {
        return s << v.x << ' ' << v.y << ' ' << v.z;
    }

    T x, y, z;
};

typedef Vec3<double> Vec3d; //!> 3D vector class with double values
typedef Vec3<float> Vec3f;  //!> 3D vector class with float values
typedef Vec3<int> Vec3i;    //!> 3D vector class with integer values


/** Class to define elementary operations between points */
template<typename T>
class Point3{
public:

    /** Constructors of Point3 class */
    Point3() : x(0), y(0), z(0) {}
    Point3(T xx) : x(xx), y(xx), z(xx) {}
    Point3(T xx, T yy, T zz) : x(xx), y(yy), z(zz) {}

    /** Defining the distance between two Point3-s */
    double distance(const Point3 &p) {
        Vec3d dif(x - p.x, y - p.y, z - p.z);
        return dif.length();
    }

    /** Defining the distance between a Point3 and dealii::Point<3> */
    double distance(const dealii::Point<3> &p) {
        Vec3d dif(x - p[0], y - p[1], z - p[2]);
        return dif.length();
    }

    /** Defining the comparison operator between two Point3-s */
    bool operator ==(const Point3 &p) {
        return x == p.x && y == p.y && z == p.z;
    }

    /** Defining the comparison operator between a Point3 and dealii::Point<3> */
    bool operator ==(const dealii::Point<3> &p) {
        return x == p[0] && y == p[1] && z == p[2];
    }

    /** Define accessors */
    const T& operator [](uint8_t i) const {
        return (&x)[i];
    }
    T& operator [](uint8_t i) {
        return (&x)[i];
    }

    /** Defining the behaviour of cout */
    friend std::ostream& operator <<(std::ostream &s, const Point3 &p) {
        return s << p.x << ' ' << p.y << ' ' << p.z;
    }

    T x, y, z;
};

typedef Point3<double> Point3d; //!> 3D point class with double values
typedef Point3<float> Point3f;  //!> 3D point class with float values
typedef Point3<int> Point3i;    //!> 3D point class with integer values

#endif /* PRIMITIVES_H_ */
