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

/** Template class to define the 3-dimensional vector with its operations */
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

/** Class to define elementary operations between 3-dimensional points */
template<typename T>
class Point3{
public:

    /** Constructors of Point3 class */
    Point3() : x(0), y(0), z(0) {}
    Point3(T xx) : x(xx), y(xx), z(xx) {}
    Point3(T xx, T yy, T zz) : x(xx), y(yy), z(zz) {}

    /** Squared distance between two Point3-s */
    T distance2(const Point3<T> &p) const {
        T xx = x - p.x;
        T yy = y - p.y;
        T zz = z - p.z;
        return (T) (xx * xx + yy * yy + zz * zz);
    }

    /** Squared distance between a Point3 and dealii::Point<3> */
    T distance2(const dealii::Point<3> &p) const {
        T xx = x - p[0];
        T yy = y - p[1];
        T zz = z - p[2];
        return (T) (xx * xx + yy * yy + zz * zz);
    }

    /** Distance between two Point3-s */
    T distance(const Point3<T> &p) const {
        T xx = x - p.x;
        T yy = y - p.y;
        T zz = z - p.z;
        return (T) sqrt(xx * xx + yy * yy + zz * zz);
    }

    /** Distance between a Point3 and dealii::Point<3> */
    T distance(const dealii::Point<3> &p) const {
        T xx = x - p[0];
        T yy = y - p[1];
        T zz = z - p[2];
        return (T) sqrt(xx * xx + yy * yy + zz * zz);
    }

    /** Function to figure out whether point is near some coordinate.
     * Non-negative return value indicates the index of coordinate near the Point,
     * -1 indicates that Point is not near any given coordinate.
     */
    const int near(const T &p0, const T eps) const {
        if(fabs(x - p0) <= eps)
            return 0;
        else
            return -1;
    }
    const int near(const T &p0, const T &p1, const T eps) const {
        if(fabs(x - p0) <= eps)
            return 0;
        else if (fabs(y - p1) <= eps)
            return 1;
        else
            return -1;
    }
    const int near(const T &p0, const T &p1, const T &p2, const T eps) const {
        if(fabs(x - p0) <= eps)
            return 0;
        else if (fabs(y - p1) <= eps)
            return 1;
        else if (fabs(z - p2) <= eps)
            return 2;
        else
            return -1;
    }

    /** Comparison operator between two Point3-s */
    bool operator ==(const Point3<T> &p) const {
        return x == p.x && y == p.y && z == p.z;
    }

    /** Comparison operator between a Point3 and dealii::Point<3> */
    bool operator ==(const dealii::Point<3> &p) {
        return x == p[0] && y == p[1] && z == p[2];
    }

    /** Point3 accessors */
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

/** Class to define elementary operations between 3-dimensional points */
template<typename T>
class Point2{
public:

    /** Constructors of Point2 class */
    Point2() : x(0), y(0) {}
    Point2(T xx) : x(xx), y(xx) {}
    Point2(T xx, T yy) : x(xx), y(yy) {}

    /** Function to figure out whether point is near some coordinate.
     * Non-negative return value indicates the index of coordinate near the Point,
     * -1 indicates that Point is not near any given coordinate.
     */
    const int near(const T &p0, const T eps) const {
        if(fabs(x - p0) <= eps)
            return 0;
        else
            return -1;
    }
    const int near(const T &p0, const T &p1, const T eps) const {
        if(fabs(x - p0) <= eps)
            return 0;
        else if (fabs(y - p1) <= eps)
            return 1;
        else
            return -1;
    }

    /** Distance between two Point2-s */
    double distance(const Point2<T> &p) const {
        T xx = x - p.x;
        T yy = y - p.y;
        return sqrt(xx * xx + yy * yy);
    }

    /** Squared distance between two Point2-s; it's a bit faster than distance */
    T distance2(const Point2<T> &p) const {
        T xx = x - p.x;
        T yy = y - p.y;
        return xx * xx + yy * yy;
    }

    /** Comparison operator between two Point2-s */
    bool operator ==(const Point2<T> &p) const {
        return x == p.x && y == p.y;
    }

    /** Point2 accessors */
    const T& operator [](uint8_t i) const {
        return (&x)[i];
    }
    T& operator [](uint8_t i) {
        return (&x)[i];
    }

    /** Defining the behaviour of cout */
    friend std::ostream& operator <<(std::ostream &s, const Point2 &p) {
        return s << p.x << ' ' << p.y;
    }

    T x, y;
};

typedef Vec3<double> Vec3d; //!> 3D vector class with double values
typedef Vec3<float> Vec3f;  //!> 3D vector class with float values
typedef Vec3<int> Vec3i;    //!> 3D vector class with integer values

typedef Point3<double> Point3d; //!> 3-dimensional point class with double values
typedef Point3<float> Point3f;  //!> 3-dimensional point class with float values
typedef Point3<int> Point3i;    //!> 3-dimensional point class with integer values

typedef Point2<double> Point2d; //!> 2-dimensional point class with double values
typedef Point2<float> Point2f;  //!> 2-dimensional point class with float values
typedef Point2<int> Point2i;    //!> 2-dimensional point class with integer values

#endif /* PRIMITIVES_H_ */
