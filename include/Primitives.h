/*
 * Primitives.h
 *
 *  Created on: 29.04.2016
 *  Author: veske & scratchapixel
 *  http://www.scratchapixel.com/code.php?id=9&origin=/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle&src=1
 */

#ifndef PRIMITIVES_H_
#define PRIMITIVES_H_

#include "Macros.h"
#include "deal.II/numerics/vector_tools.h"
#include <math.h>

namespace femocs {

/** Template class to enable iterators and range-based for-loops in other Primitives */
template<typename T, typename R>
class Iterator {
public:
    Iterator (const T* dat, uint8_t i) : data(dat), iter_indx(i) {}
    /** Function to compare two iterators */
    bool operator!= (const Iterator& rhs) const { return iter_indx != rhs.iter_indx; }
    /** Function to access member of iterator */
    R operator* () const { return (*data)[iter_indx]; }
    /** Function to advance iterator */
    const Iterator& operator++ () {
        ++iter_indx;
        return *this;
    }

private:
    uint8_t iter_indx;
    const T* data;
};

/** Template class for finite element cell */
template <size_t dim>
class SimpleCell {
public:

    /** SimpleCell constructors */
    SimpleCell() { std::fill_n(node, dim, int()); }
    SimpleCell(const unsigned int &nn) { std::fill_n(node, dim, nn); }

    /** Less than, bigger than, less than or equal, bigger than or equal operators */
    vector<bool> operator <(const unsigned int &t) const {
        vector<bool> v; v.reserve(dim);
        for (unsigned int n : node) v.push_back(n < t);
        return v;
    }
    vector<bool> operator >(const unsigned int &t) const {
        vector<bool> v; v.reserve(dim);
        for (unsigned int n : node) v.push_back(n > t);
        return v;
    }
    vector<bool> operator <=(const unsigned int &t) const {
        vector<bool> v; v.reserve(dim);
        for (unsigned int n : node) v.push_back(n <= t);
        return v;
    }
    vector<bool> operator >=(const unsigned int &t) const {
        vector<bool> v; v.reserve(dim);
        for (unsigned int n : node) v.push_back(n >= t);
        return v;
    }

    /** Define access operators or accessors */
    const unsigned int operator [](size_t i) const {
        require(i >= 0 && i < dim, "Invalid index: " + to_string(i));
        return node[i];
    }
    unsigned int operator [](size_t i) {
        require((i >= 0) && (i < dim), "Invalid index: " + to_string(i));
        return node[i];
    }

    /** Define the behaviour of string stream */
    friend std::ostream& operator <<(std::ostream &s, const SimpleCell<dim> &t) {
        for (unsigned int nn : t.node) s << nn << ' ';
        return s;
    }

    /** Return data as string */
    string to_str() const { stringstream ss; ss << this; return ss.str(); }

    /** Transform SimpleCell to vector */
    vector<int> to_vector() const { return vector<int>(std::begin(node), std::end(node)); }

    /** Attach iterator */
    typedef Iterator<SimpleCell, unsigned int> iterator;
    iterator begin() const { return iterator(this, 0); }
    iterator end() const { return iterator(this, dim); }

    unsigned int node[dim]; //!< vertices of SimpleCell
};

/** Node class without Point data */
class SimpleNode: public SimpleCell<1> {
public:
    /** SimpleNode constructors */
    SimpleNode() : SimpleCell<1>() {}
    SimpleNode(const unsigned int &n1) : SimpleCell<1>(n1) {}
};

/** Edge class without Point data */
class SimpleEdge: public SimpleCell<2> {
public:
    /** SimpleEdge constructors */
    SimpleEdge() : SimpleCell<2>() {}
    SimpleEdge(const unsigned int &n1) : SimpleCell<2>(n1) {}
    SimpleEdge(const unsigned int &n1, const unsigned int &n2) {
        node[0] = n1; node[1] = n2;
    }
};

/** Face class without Point data */
class SimpleFace : public SimpleCell<3> {
public:
    /** SimpleFace constructors */
    SimpleFace() : SimpleCell<3>() {}
    SimpleFace(const unsigned int &n1) : SimpleCell<3>(n1) {}
    SimpleFace(const unsigned int &n1, const unsigned int &n2, const unsigned int &n3) {
        node[0] = n1; node[1] = n2; node[2] = n3;
    }

    /** Get i-th edge of the face */
    const SimpleEdge edge(const unsigned int i) const {
        if (i <= 0) return SimpleEdge(node[0], node[1]);
        if (i == 1) return SimpleEdge(node[0], node[2]);
        else        return SimpleEdge(node[1], node[2]);
    }
};

/** Element class without Point data */
class SimpleElement : public SimpleCell<4> {
public:
    /** SimpleEdge constructors */
    SimpleElement() : SimpleCell<4>() {}
    SimpleElement(const unsigned int &n1) : SimpleCell<4>(n1) {}
    SimpleElement(const unsigned int &n1, const unsigned int &n2, const unsigned int &n3, const unsigned int &n4) {
        node[0] = n1; node[1] = n2; node[2] = n3; node[3] = n4;
    }

    /** Get i-th edge of the element */
    const SimpleEdge edge(const unsigned int i) const {
        if (i <= 0) return SimpleEdge(node[0], node[1]);
        if (i == 1) return SimpleEdge(node[0], node[2]);
        if (i == 2) return SimpleEdge(node[0], node[3]);
        if (i == 3) return SimpleEdge(node[1], node[2]);
        if (i == 4) return SimpleEdge(node[1], node[3]);
        else        return SimpleEdge(node[2], node[3]);
    }

    /** Get i-th face of the element */
    const SimpleFace face(const unsigned int i) const {
        if (i <= 0) return SimpleFace(node[0], node[1], node[2]);
        if (i == 1) return SimpleFace(node[0], node[1], node[3]);
        if (i == 2) return SimpleFace(node[1], node[2], node[3]);
        else        return SimpleFace(node[2], node[0], node[1]);
    }
};

/** Template class to define elementary operations between 3-dimensional points */
template<typename T>
class Point3_T{
public:

    /** Constructors of Point3 class */
    Point3_T() : x(0), y(0), z(0), r(0) {}
    Point3_T(T xx) : x(xx), y(xx), z(xx), r(0) {}
    Point3_T(T xx, T yy, T zz) : x(xx), y(yy), z(zz), r(0) {}

    /** Dimensionality of SimpleCell_T */
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

    /** Squared distance between a Point3 and dealii::SimpleCell_T<3> */
    const double distance2(const dealii::Point<3> &p) const {
        T xx = x - p[0];
        T yy = y - p[1];
        T zz = z - p[2];
        return (T) (xx * xx + yy * yy + zz * zz);
    }

    /** Squared distance between two Point3-s taking into account the simulation cell periodicity.
     * Period == 0 in some direction gives the result without periodicity in that direction.
     */
    const double periodic_distance2(const Point3_T<T> &p, double period_x = 0, double period_y = 0,
            double period_z = 0) const {
        T dx = fabs(x - p.x);
        T dy = fabs(y - p.y);
        T dz = fabs(z - p.z);

        dx = min(dx, fabs(dx - period_x)); // apply periodic boundary condition in x-direction
        dy = min(dy, fabs(dy - period_y)); // apply periodic boundary condition in y-direction
        dz = min(dz, fabs(dz - period_z)); // apply periodic boundary condition in z-direction

        return dx * dx + dy * dy + dz * dz;
    }

    /** Distance between two Point3-s */
    const double distance(const Point3_T<T> &p) const {
        return sqrt(distance2(p));
    }

    /** Distance between a Point3 and dealii::SimpleCell_T<3> */
    const double distance(const dealii::Point<3> &p) const {
        return sqrt(distance2(p));
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
    /** Multiplying a SimpleCell_T with constant */
    /** Scalar multiplication of vector with a scalar and with another vector */
    Point3_T operator *(const T &r) const {
        return Point3_T(x * r, y * r, z * r);
    }
    Point3_T& operator *=(const T &r) {
        x *= r, y *= r, z *= r;
        return *this;
    }
    /** Dividing a SimpleCell_T with constant */
    Point3_T operator /(const T &r) const {
        return Point3_T(x / r, y / r, z / r);
    }
    Point3_T& operator /=(const T &r) {
        x /= r, y /= r, z /= r;
        return *this;
    }

    /** Comparison operator between two Point3-s */
    const bool operator ==(const Point3_T<T> &p) const
        { return x == p.x && y == p.y && z == p.z; }
    bool operator ==(const Point3_T<T> &p)
        { return x == p.x && y == p.y && z == p.z; }

    const bool operator >(const Point3_T<T> &p) const { return r > p.r; }
    bool operator >(const Point3_T<T> &p) { return r > p.r; }
    const bool operator <(const Point3_T<T> &p) const { return r < p.r; }
    bool operator <(const Point3_T<T> &p) { return r < p.r; }
    const bool operator >=(const Point3_T<T> &p) const { return r >= p.r; }
    bool operator >=(const Point3_T<T> &p) { return r >= p.r; }
    const bool operator <=(const Point3_T<T> &p) const { return r <= p.r; }
    bool operator <=(const Point3_T<T> &p) { return r <= p.r; }

    /** Comparison operator between a Point3 and dealii::SimpleCell_T<3> */
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

    /** Attach iterator */
    typedef Iterator<Point3_T, T> iterator;
    iterator begin() const { return iterator(this, 0); }
    iterator end() const { return iterator(this, size()); }

    /** Define the behaviour of string stream */
    friend std::ostream& operator <<(std::ostream &s, const Point3_T &p) {
        return s << p.x << ' ' << p.y << ' ' << p.z;
    }

    /** Return data as string */
    string to_str() const { stringstream ss; ss << (*this); return ss.str(); }

    T x, y, z, r;
};

/** Template class to define elementary operations between 2-dimensional points */
template<typename T>
class Point2_T{
public:
    /** Constructors of Point2 class */
    Point2_T() : x(0), y(0), r(0) {}
    Point2_T(T xx) : x(xx), y(xx), r(0) {}
    Point2_T(T xx, T yy) : x(xx), y(yy), r(0) {}

    /** Dimensionality of SimpleCell_T */
    const int size() const {
        return 2;
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
    bool operator ==(const Point2_T<T> &p) const
        { return x == p.x && y == p.y; }
    bool operator ==(const Point2_T<T> &p)
        { return x == p.x && y == p.y; }

    const bool operator >(const Point2_T<T> &p) const { return r > p.r; }
    bool operator >(const Point2_T<T> &p) { return r > p.r; }
    const bool operator <(const Point2_T<T> &p) const { return r < p.r; }
    bool operator <(const Point2_T<T> &p) { return r < p.r; }
    const bool operator >=(const Point2_T<T> &p) const { return r >= p.r; }
    bool operator >=(const Point2_T<T> &p) { return r >= p.r; }
    const bool operator <=(const Point2_T<T> &p) const { return r <= p.r; }
    bool operator <=(const Point2_T<T> &p) { return r <= p.r; }

    /** Point2 accessors */
    const T& operator [](uint8_t i) const {
        return (&x)[i];
    }
    T& operator [](uint8_t i) {
        return (&x)[i];
    }

    /** Attach iterator */
    typedef Iterator<Point2_T, T> iterator;
    iterator begin() const { return iterator(this, 0); }
    iterator end() const { return iterator(this, size()); }

    /** Define the behaviour of string stream */
    friend std::ostream& operator <<(std::ostream &s, const Point2_T &p) {
        return s << p.x << ' ' << p.y;
    }

    /** Return data as string */
    string to_str() const { stringstream ss; ss << (*this); return ss.str(); }

    T x, y, r;
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
    Vec3_T& operator +=(const Vec3_T &v) {
        x += v.x, y += v.y, z += v.z;
        return *this;
    }

    /** Subtraction of two vectors */
    Vec3_T operator -(const Vec3_T &v) const {
        return Vec3_T(x - v.x, y - v.y, z - v.z);
    }
    Vec3_T& operator -=(const Vec3_T &v) {
        x -= v.x, y -= v.y, z -= v.z;
        return *this;
    }

    /** Scalar multiplication of vector with a scalar and with another vector */
    Vec3_T operator *(const T &r) const {
        return Vec3_T(x * r, y * r, z * r);
    }
    Vec3_T operator *(const Vec3_T &v) const {
        return Vec3_T(x * v.x, y * v.y, z * v.z);
    }
    friend Vec3_T operator *(const T &r, const Vec3_T &v) {
        return Vec3_T<T>(v.x * r, v.y * r, v.z * r);
    }
    Vec3_T& operator *=(const T &r) {
        x *= r, y *= r, z *= r;
        return *this;
    }

    /** Scalar division of vector with a scalar or with another vector */
    Vec3_T operator /(const T &r) const {
        return Vec3_T(x / r, y / r, z / r);
    }
    Vec3_T operator /(const Vec3_T &v) const {
        return Vec3_T(x / v.x, y / v.y, z / v.z);
    }
    friend Vec3_T operator /(const T &r, const Vec3_T &v) {
        return Vec3_T<T>(r / v.x, r / v.y, r / v.z);
    }
    Vec3_T& operator /=(const T &r) {
        x /= r, y /= r, z /= r;
        return *this;
    }

    /** Equals operator */
    bool operator ==(const Vec3_T &v) const {
        return x == v.x && y == v.y && z == v.z;
    }

    /** Dot product, cross product, norm and length */
    T dotProduct(const Vec3_T<T> &v) const {
        return x * v.x + y * v.y + z * v.z;
    }
    Vec3_T crossProduct(const Vec3_T<T> &v) const {
        return Vec3_T<T>(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x);
    }
    T length2() const {
        return x * x + y * y + z * z;
    }
    T length() const {
        return sqrt(length2());
    }

    /** Function to normalize the vector */
    Vec3_T& normalize() {
        T n = length();
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
    const T& operator [](uint8_t i) const { return (&x)[i]; }
    T& operator [](uint8_t i) { return (&x)[i]; }

    /** Attach iterator */
    typedef Iterator<Vec3_T, T> iterator;
    iterator begin() const { return iterator(this, 0); }
    iterator end() const { return iterator(this, size()); }

    /** Define the behaviour of string stream */
    friend std::ostream& operator <<(std::ostream &s, const Vec3_T<T> &v) {
        return s << v.x << ' ' << v.y << ' ' << v.z;
    }

    /** Return data as string */
    string to_str() const { stringstream ss; ss << (*this); return ss.str(); }

    T x, y, z;
};

/** Template class to define the 4-dimensional vector with its operations */
template<typename T>
class Vec4_T {
public:
    /** Vec4 constructors */
    Vec4_T() : x(T(0)), y(T(0)), z(T(0)), w(T(0)) {}
    Vec4_T(T xx) : x(xx), y(xx), z(xx), w(xx) {}
    Vec4_T(T xx, T yy, T zz, T ww) : x(xx), y(yy), z(zz), w(ww) {}
    Vec4_T(const Vec3_T<T> &v, T ww) : x(v.x), y(v.y), z(v.z), w(ww) {}
    Vec4_T(const Point3_T<T> &p, T ww) : x(p.x), y(p.y), z(p.z), w(ww) {}

    /** Dimensionality of vector */
    const int size() const {
        return 4;
    }

    /** Addition of two vectors */
    Vec4_T operator +(const Vec4_T &v) const {
        return Vec4_T(x + v.x, y + v.y, z + v.z, w + v.w);
    }
    Vec4_T& operator +=(const Vec4_T &v) {
        x += v.x, y += v.y, z += v.z, w += v.w;
        return *this;
    }

    /** Subtraction of two vectors */
    Vec4_T operator -(const Vec4_T &v) const {
        return Vec4_T(x - v.x, y - v.y, z - v.z, w - v.w);
    }
    Vec4_T& operator -=(const Vec4_T &v) {
        x -= v.x, y -= v.y, z -= v.z, w -= v.w;
        return *this;
    }

    /** Scalar multiplication of vector with a scalar and with another vector */
    Vec4_T operator *(const T &r) const {
        return Vec4_T(x * r, y * r, z * r, w * r);
    }
    Vec4_T operator *(const Vec4_T &v) const {
        return Vec4_T(x * v.x, y * v.y, z * v.z, w * v.w);
    }
    Vec4_T& operator *=(const T &r) {
        x *= r, y *= r, z *= r, w *= r;
        return *this;
    }
    Vec4_T& operator *=(const Vec4_T &v) const {
        x *= v.x, y *= v.y, z *= v.z, w *= v.w;
        return *this;
    }

    /** Scalar division of vector with a scalar or with another vector */
    Vec4_T operator /(const T &r) const {
        return Vec4_T(x / r, y / r, z / r, w / r);
    }
    Vec4_T operator /(const Vec4_T &v) const {
        return Vec4_T(x / v.x, y / v.y, z / v.z, w / v.w);
    }
    Vec4_T& operator /=(const T &r) {
        x /= r, y /= r, z /= r, w /= r;
        return *this;
    }
    Vec4_T& operator /=(const Vec4_T &v) const {
        x /= v.x, y /= v.y, z /= v.z, w /= v.w;
        return *this;
    }

    /** Equals operator */
    bool operator ==(const Vec4_T &v) const {
        return x == v.x && y == v.y && z == v.z && w == v.w;
    }

    /** Dot product and length */
    T dotProduct(const Vec4_T<T> &v) const {
        return x * v.x + y * v.y + z * v.z + w * v.w;
    }
    T length2() const {
        return x * x + y * y + z * z + w * w;
    }
    T length() const {
        return sqrt(length2());
    }

    /** Function to normalize the vector */
    Vec4_T& normalize() {
        T n = length();
        if (n > 0) {
            T factor = 1 / sqrt(n);
            x *= factor, y *= factor, z *= factor;
        }
        return *this;
    }

    /** Define access operators */
    const T& operator [](uint8_t i) const { return (&x)[i]; }
    T& operator [](uint8_t i) { return (&x)[i]; }

    /** Attach iterator */
    typedef Iterator<Vec4_T, T> iterator;
    iterator begin() const { return iterator(this, 0); }
    iterator end() const { return iterator(this, size()); }

    /** Define the behaviour of string stream */
    friend std::ostream& operator <<(std::ostream &s, const Vec4_T<T> &v) {
        return s << v.x << ' ' << v.y << ' ' << v.z << ' ' << v.w;
    }

    /** Return data as string */
    string to_str() const { stringstream ss; ss << (*this); return ss.str(); }

    T x, y, z, w;
};

typedef Vec4_T<double> Vec4;     //!> 4-dimensional vector class with double values
typedef Vec3_T<double> Vec3;     //!> 3-dimensional vector class with double values
typedef Point3_T<double> Point3; //!> 3-dimensional point class with double values
typedef Point2_T<double> Point2; //!> 2-dimensional point class with double values


/** Class to define elementary operations between atoms */
class Atom {
public:

    /** Constructors of Atom class */
    Atom() : id(0), point(0), coord(0), sort_indx(0) {}
    Atom(int i, Point3 &p, int c) : id(i), point(p.x, p.y, p.z), coord(c), sort_indx(0) {}
    Atom(const int i, const Point3 &p, const int c) : id(i), point(p.x, p.y, p.z), coord(c), sort_indx(0) {}

    /** Comparison operator between two Atom-s */
    const bool operator ==(const Atom &a) const
        { return id == a.id && point == a.point && coord == a.coord; }
    bool operator ==(const Atom &a)
        { return id == a.id && point == a.point && coord == a.coord; }

    /** Define the behaviour of string stream */
    friend std::ostream& operator <<(std::ostream &s, const Atom &a) {
        return s << a.id << ' ' << a.point << ' ' << a.coord;
    }

    /** Functor for sorting atoms in ascending order by their ID */
    struct sort_id {
        inline bool operator() (const Atom& a1, const Atom& a2) { return a1.id < a2.id; }
    };

    /** Functor for sorting atoms in ascending order by their x, y, z or radial coordinate */
    struct sort_up {
        int I;           //!< coordinate along which the atoms are sorted; 0=x, 1=y, 2=z, 3=radial
        sort_up(int d = 3) : I(d) {}
        inline bool operator() (const Atom& a1, const Atom& a2) { return a1.point[I] < a2.point[I]; }
    };

    /** Functor for sorting atoms in descending order by their x, y, z or radial coordinate */
    struct sort_down {
        int I;           //!< coordinate along which the atoms are sorted; 0=x, 1=y, 2=z, 3=radial
        sort_down(int d = 3) : I(d) {}
        inline bool operator() (const Atom& a1, const Atom& a2) { return a1.point[I] > a2.point[I]; }
    };

    /** Functor for sorting atoms in ascending order first by their x, y, z or radial coordinate
     * and then in ascending order by another x, y, z or radial coordinate */
    struct sort_up2 {
        int I1, I2;      //!< coordinates along which the atoms are sorted; 0=x, 1=y, 2=z, 3=radial
        sort_up2(int coord1, int coord2) : I1(coord1), I2(coord2) {}
        inline bool operator() (const Atom& a1, const Atom& a2) {
            return (a1.point[I1] < a2.point[I1]) || ((a1.point[I1] == a2.point[I1]) && (a1.point[I2] < a2.point[I2]));
        }
    };

    /** Functor for sorting atoms in descending order first by their x, y, z or radial coordinate
     * and then in descending order by another x, y, z or radial coordinate */
    struct sort_down2 {
        int I1, I2;      //!< coordinates along which the atoms are sorted; 0=x, 1=y, 2=z, 3=radial
        sort_down2(int coord1, int coord2) : I1(coord1), I2(coord2) {}
        inline bool operator() (const Atom& a1, const Atom& a2) {
            return (a1.point[I1] > a2.point[I1]) || ((a1.point[I1] == a2.point[I1]) && (a1.point[I2] > a2.point[I2]));
        }
    };

    /** Functor for sorting atoms in ascending order first by their x, y, z or radial coordinate
     * and then in descending order by another x, y, z or radial coordinate */
    struct sort_up_down {
        int I1, I2;      //!< coordinates along which the atoms are sorted; 0=x, 1=y, 2=z, 3=radial
        sort_up_down(int coord1, int coord2) : I1(coord1), I2(coord2) {}
        inline bool operator() (const Atom& a1, const Atom& a2) {
            return (a1.point[I1] < a2.point[I1]) || ((a1.point[I1] == a2.point[I1]) && (a1.point[I2] > a2.point[I2]));
        }
    };

    /** Functor for sorting atoms in descending order first by their x, y, z or radial coordinate
     * and then in ascending order by another x, y, z or radial coordinate */
    struct sort_down_up {
        int I1, I2;      //!< coordinates along which the atoms are sorted; 0=x, 1=y, 2=z, 3=radial
        sort_down_up(int coord1, int coord2) : I1(coord1), I2(coord2) {}
        inline bool operator() (const Atom& a1, const Atom& a2) {
            return (a1.point[I1] > a2.point[I1]) || ((a1.point[I1] == a2.point[I1]) && (a1.point[I2] < a2.point[I2]));
        }
    };

    Point3 point;
    int id;
    int coord;
    int sort_indx;
};

/** Class to hold solution data and its operations */
class Solution {
public:
    Solution() : elfield(Vec3(0)), el_norm(0), potential(0), sort_indx(0) {}
    Solution(double d) : elfield(Vec3(d)), el_norm(d), potential(d), sort_indx(0) {}
    Solution(Vec3& v, double en, double pot) : elfield(v), el_norm(en), potential(pot), sort_indx(0) {}
    Solution(const Vec3& v, const double en, const double pot) : elfield(v), el_norm(en), potential(pot), sort_indx(0) {}

    /** Define the behaviour of string stream */
    friend std::ostream& operator <<(std::ostream &ss, const Solution &sol) {
        return ss << sol.elfield << ' ' << sol.el_norm << ' ' << sol.potential;
    }

    /** Functors used to sort vector of Solution into same order as vector of Atom */
    struct sort_up {
        inline bool operator() (const Solution& lh, const Solution& rh) { return lh.sort_indx < rh.sort_indx; }
    };
    struct sort_down {
        inline bool operator() (const Solution& lh, const Solution& rh) { return lh.sort_indx > rh.sort_indx; }
    };

    Vec3 elfield;
    double el_norm;
    double potential;
    int sort_indx;
};

} // namaspace femocs

#endif /* PRIMITIVES_H_ */
