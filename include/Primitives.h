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

namespace femocs {

/** Template class to define the operations between the indices of face nodes */
template<typename T>
class SimpleFace_T {
public:
    /** SimpleFace constructors */
    SimpleFace_T() : n1(0), n2(0), n3(0){}
    SimpleFace_T(T n) : n1(n), n2(n), n3(n) {}
    SimpleFace_T(T nn1, T nn2, T nn3) : n1(nn1), n2(nn2), n3(nn3) {}

    /** Number of vertices in SimpleFace */
    const T n_verts() const {
        return 3;
    }

    /** Equals operator */
    bool operator ==(const SimpleFace_T &t) const {
        return (n1 == t.n1 || n1 == t.n2 || n1 == t.n3)
                && (n2 == t.n1 || n2 == t.n2 || n2 == t.n3)
                && (n3 == t.n1 || n3 == t.n2 || n3 == t.n3);
    }

    /** Less than, less than or equal, bigger than, bigger than or equal operators */
    vector<bool> operator <=(const T &t) const {
        return {n1 <= t, n2 <= t, n3 <= t};
    }
    vector<bool> operator <(const T &t) const {
        return {n1 < t, n2 < t, n3 < t};
    }
    vector<bool> operator >=(const T &t) const {
        return {n1 >= t, n2 >= t, n3 >= t};
    }
    vector<bool> operator >(const T &t) const {
        return {n1 > t, n2 > t, n3 > t};
    }

    /** Define access operators or accessors. */
    const T& operator [](uint8_t i) const {
        require(i < 3, "Invalid index!");
        return (&n1)[i];
    }
    T& operator [](uint8_t i) {
        require(i < 3, "Invalid index!");
        return (&n1)[i];
    }

    /** Defining the behaviour of cout */
    friend std::ostream& operator <<(std::ostream &s, const SimpleFace_T<T> &t) {
        return s << t.n1 << ' ' << t.n2 << ' ' << t.n3;
    }

    /** Transform SimpleFace to vector */
    vector<int> to_vector() const {
        return vector<int> {n1, n2, n3};
    }

    T n1, n2, n3;
};

/** Template class to define the operations between the indices of element nodes */
template<typename T>
class SimpleElement_T {
public:
    /** SimpleFace constructors */
    SimpleElement_T() : n1(0), n2(0), n3(0), n4(0){}
    SimpleElement_T(T n) : n1(n), n2(n), n3(n), n4(n) {}
    SimpleElement_T(T nn1, T nn2, T nn3, T nn4) : n1(nn1), n2(nn2), n3(nn3), n4(nn4) {}

    /** Number of vertices in tetrahedron */
    const T n_verts() const {
        return 4;
    }

    /** Transform SimpleElement to vector */
    vector<int> to_vector() const {
        return vector<int> {n1, n2, n3, n4};
    }

    /** Equals operator */
    bool operator ==(const SimpleElement_T &t) const {
        return (n1 == t.n1 || n1 == t.n2 || n1 == t.n3 || n1 = t.n4)
                && (n2 == t.n1 || n2 == t.n2 || n2 == t.n3 || n2 = t.n4)
                && (n3 == t.n1 || n3 == t.n2 || n3 == t.n3 || n3 = t.n4)
                && (n4 == t.n1 || n4 == t.n2 || n4 == t.n3 || n4 = t.n4);
        }

    /** Less than, less than or equal, bigger than, bigger than or equal operators */
    vector<bool> operator <=(const T &t) const {
        return {n1 <= t, n2 <= t, n3 <= t, n4 <= t};
    }
    vector<bool> operator <(const T &t) const {
        return {n1 < t, n2 < t, n3 < t, n4 < t};
    }
    vector<bool> operator >=(const T &t) const {
        return {n1 >= t, n2 >= t, n3 >= t, n4 >= t};
    }
    vector<bool> operator >(const T &t) const {
        return {n1 > t, n2 > t, n3 > t, n4 > t};
    }

    /** Define access operators or accessors */
    const T& operator [](uint8_t i) const {
        require(i < 4, "Invalid index!");
        return (&n1)[i];
    }
    T& operator [](uint8_t i) {
        require(i < 4, "Invalid index!");
        return (&n1)[i];
    }

    /** Defining the behaviour of cout */
    friend std::ostream& operator <<(std::ostream &s, const SimpleElement_T &t) {
        return s << t.n1 << ' ' << t.n2 << ' ' << t.n3 << ' ' << t.n4;
    }

    T n1, n2, n3, n4;
};

/** Template class to define the 3-dimensional vector with its operations */
template<typename T>
class Vec3_T {
public:
    /** Vec3 constructors */
    Vec3_T() : x(T(0)), y(T(0)), z(T(0)) {}
    Vec3_T(T xx) : x(xx), y(xx), z(xx) {}
    Vec3_T(T xx, T yy, T zz) : x(xx), y(yy), z(zz) {}

    /** Dimensionality of vector */
    const int size() const {
        return 3;
    }

    /** Addition of two vectors */
    Vec3_T operator +(const Vec3_T &v) const {
        return Vec3_T(x + v.x, y + v.y, z + v.z);
    }

    /** Subtraction of two vectors */
    Vec3_T operator -(const Vec3_T &v) const {
        return Vec3_T(x - v.x, y - v.y, z - v.z);
    }
    Vec3_T operator -() const {
        return Vec3_T(-x, -y, -z);
    }

    /** Scalar multiplication of vector with a scalar and with another vector */
    Vec3_T operator *(const T &r) const {
        return Vec3_T(x * r, y * r, z * r);
    }
    Vec3_T operator *(const Vec3_T &v) const {
        return Vec3_T(x * v.x, y * v.y, z * v.z);
    }
    Vec3_T& operator *=(const T &r) {
        x *= r, y *= r, z *= r;
        return *this;
    }
    friend Vec3_T operator *(const T &r, const Vec3_T &v) {
        return Vec3_T<T>(v.x * r, v.y * r, v.z * r);
    }

    /** Equals operator */
    bool operator ==(const Vec3_T &v) const {
        return x == v.x && y == v.y && z == v.z;
    }

    /** Scalar division of vector with a scalar or with another vector */
    Vec3_T operator /(const T &r) const {
        return Vec3_T(x / r, y / r, z / r);
    }
    Vec3_T& operator /=(const T &r) {
        x /= r, y /= r, z /= r;
        return *this;
    }
    friend Vec3_T operator /(const T &r, const Vec3_T &v) {
        return Vec3_T<T>(r / v.x, r / v.y, r / v.z);
    }

    /** Dot product, cross product, norm and length */
    T dotProduct(const Vec3_T<T> &v) const {
        return x * v.x + y * v.y + z * v.z;
    }
    Vec3_T crossProduct(const Vec3_T<T> &v) const {
        return Vec3_T<T>(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x);
    }
    T norm() const {
        return x * x + y * y + z * z;
    }
    T length() const {
        return sqrt(norm());
    }

    /** Function to normalize the vector */
    Vec3_T& normalize() {
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
    friend std::ostream& operator <<(std::ostream &s, const Vec3_T<T> &v) {
        return s << v.x << ' ' << v.y << ' ' << v.z;
    }

    T x, y, z;
};

/** Class to define elementary operations between 3-dimensional points */
template<typename T>
class Point3_T{
public:

    /** Constructors of Point3 class */
    Point3_T() : x(0), y(0), z(0) {}
    Point3_T(T xx) : x(xx), y(xx), z(xx) {}
    Point3_T(T xx, T yy, T zz) : x(xx), y(yy), z(zz) {}

    /** Dimensionality of Point */
    const int size() const {
        return 3;
    }

    /** Squared distance between two Point3-s */
    const double distance2(const Point3_T<T> &p) const {
        T xx = x - p.x;
        T yy = y - p.y;
        T zz = z - p.z;
        return (T) (xx * xx + yy * yy + zz * zz);
    }

    /** Squared distance between a Point3 and dealii::Point<3> */
    const double distance2(const dealii::Point<3> &p) const {
        T xx = x - p[0];
        T yy = y - p[1];
        T zz = z - p[2];
        return (T) (xx * xx + yy * yy + zz * zz);
    }

    /** Distance between two Point3-s */
    const double distance(const Point3_T<T> &p) const {
        T xx = x - p.x;
        T yy = y - p.y;
        T zz = z - p.z;
        return (T) sqrt(xx * xx + yy * yy + zz * zz);
    }

    /** Distance between a Point3 and dealii::Point<3> */
    const double distance(const dealii::Point<3> &p) const {
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

    /** Subtraction of two Point3-d */
    Point3_T<T> operator -(const Point3_T<T> &p) const {
        return Point3_T(x - p.x, y - p.y, z - p.z);
    }

    /** Adding a point to existing one */
    Point3_T& operator +=(const Point3_T<T> &p) {
        x += p.x, y += p.y, z += p.z;
        return *this;
    }
    /** Subtracting a point from existing one */
    Point3_T& operator -=(const Point3_T<T> &p) {
        x -= p.x, y -= p.y, z -= p.z;
        return *this;
    }
    /** Multiplying a Point with constant */
    Point3_T& operator *=(const T &r) {
        x *= r, y *= r, z *= r;
        return *this;
    }
    /** Dividing a Point with constant */
    Point3_T& operator /=(const T &r) {
        x /= r, y /= r, z /= r;
        return *this;
    }



    /** Comparison operator between two Point3-s */
    const bool operator ==(const Point3_T<T> &p) const {
        return x == p.x && y == p.y && z == p.z;
    }
    bool operator ==(const Point3_T<T> &p) {
        return x == p.x && y == p.y && z == p.z;
    }

    /** Comparison operator between a Point3 and dealii::Point<3> */
    const bool operator ==(const dealii::Point<3> &p) const {
        return x == p[0] && y == p[1] && z == p[2];
    }
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
    friend std::ostream& operator <<(std::ostream &s, const Point3_T &p) {
        return s << p.x << ' ' << p.y << ' ' << p.z;
    }

    T x, y, z;
};

/** Class to define elementary operations between 2-dimensional points */
template<typename T>
class Point2_T{
public:

    /** Constructors of Point2 class */
    Point2_T() : x(0), y(0) {}
    Point2_T(T xx) : x(xx), y(xx) {}
    Point2_T(T xx, T yy) : x(xx), y(yy) {}

    /** Dimensionality of Point */
    const int size() const {
        return 2;
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

    /** Distance between two Point2-s */
    const double distance(const Point2_T<T> &p) const {
        T xx = x - p.x;
        T yy = y - p.y;
        return sqrt(xx * xx + yy * yy);
    }

    /** Squared distance between two Point2-s; it's a bit faster than distance */
    const double distance2(const Point2_T<T> &p) const {
        T xx = x - p.x;
        T yy = y - p.y;
        return xx * xx + yy * yy;
    }

    /** Subtraction of two Point2-d */
    Point2_T<T> operator -(const Point2_T<T> &p) const {
        return Point2_T(x - p.x, y - p.y);
    }

    /** Comparison operator between two Point2-s */
    bool operator ==(const Point2_T<T> &p) const {
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
    friend std::ostream& operator <<(std::ostream &s, const Point2_T &p) {
        return s << p.x << ' ' << p.y;
    }

    T x, y;
};

typedef SimpleFace_T<unsigned int> SimpleFace;          //!> face class without Point data
typedef SimpleElement_T<unsigned int> SimpleElement;    //!> element class without Point data

typedef Vec3_T<double> Vec3;     //!> 3-dimensional vector class with double values
typedef Point3_T<double> Point3; //!> 3-dimensional point class with double values
typedef Point2_T<double> Point2; //!> 2-dimensional point class with double values

} // namaspace femocs

#endif /* PRIMITIVES_H_ */
