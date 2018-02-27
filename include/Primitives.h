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
#include "Globals.h"

namespace femocs {

/** Template class to enable iterators and range-based for-loops in other Primitives */
template<typename T, typename R>
class Iterator {
public:
    Iterator (const T* dat, unsigned int i) : iter_indx(i), data(dat) {}
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
    unsigned int iter_indx;
    const T* data;
};

template <size_t dim>
class Vec {
public:

    /** Vec constructors */
    Vec() { std::fill_n(node, dim, double()); }
    Vec(const double &nn) { std::fill_n(node, dim, nn); }

    /** Dimensionality of vector */
    int size() const { return dim; }

    /** Check if one of the nodes equals to the one of interest */
    bool operator ==(const double &t) const {
        for (double n : node) if (n == t) return true;
        return false;
    }

    bool operator ==(const Vec<dim> &t) const {
        for (size_t i = 0; i < dim; ++i) if ((*this)[i] != t[i]) return false;
        return true;
    }

    /** Addition of two vectors */
    Vec<dim>& operator +=(const Vec<dim> &p) { for(size_t i = 0; i < dim; ++i) node[i] += p[i]; return *this; }

    /** Subtraction of two vectors */
    Vec<dim>& operator -=(const Vec<dim> &p) { for(size_t i = 0; i < dim; ++i) node[i] -= p[i]; return *this; }

    /** Multiplication of two vectors */
    Vec<dim>& operator *=(const Vec<dim> &p) { for(size_t i = 0; i < dim; ++i) node[i] *= p[i]; return *this; }

    /** Division of two vectors */
    Vec<dim>& operator /=(const Vec<dim> &p) { for(size_t i = 0; i < dim; ++i) node[i] /= p[i]; return *this; }

    /** Adding a scalar to the vector */
    Vec<dim>& operator +=(const double &t) { for(double &n : node) n += t; return *this; }

    /** Subtracting a scalar from the vector */
    Vec<dim>& operator -=(const double &t) { for(double &n : node) n -= t; return *this; }

    /** Multiplying vector with a scalar */
    Vec<dim>& operator *=(const double &t) { for(double &n : node) n *= t; return *this; }

    /** Dividing vector with a scalar */
    Vec<dim>& operator /=(const double &t) { for(double &n : node) n /= t; return *this; }

    /** Dot product of two vectors */
    double dotProduct(const Vec<dim> &v) const {
        double product = 0;
        for (size_t i = 0; i < dim; ++i)
            product += node[i] *= v[i];
        return product;
    }

    /** Define the behaviour of string stream */
    friend std::ostream& operator <<(std::ostream &s, const Vec<dim> &t) {
        for (double nn : t.node) s << nn << ' ';
        return s;
    }

    /** Return data as string */
    string to_str() const { stringstream ss; ss << this; return ss.str(); }

    /** Accessor for accessing the i-th node */
    const double& operator [](const size_t i) const {
        require(i >= 0 && i < dim, "Invalid index: " + to_string(i));
        return node[i];
    }

    /** Iterator to access the cell nodes */
    typedef Iterator<Vec<dim>, double> iterator;
    iterator begin() const { return iterator(this, 0); }
    iterator end() const { return iterator(this, dim); }

    double node[dim]; ///< Vec data
};

/** Class to define basic operations with 2-dimensional points */
class Point2 {
public:
    /** Constructors of Point2 class */
    Point2() : x(0), y(0) {}
    Point2(const double xx) : x(xx), y(xx) {}
    Point2(const double xx, const double yy) : x(xx), y(yy) {}

    /** Dimensionality of point */
    int size() const { return 2; }

    /** Squared distance between two points; it's a bit faster than distance */
    double distance2(const Point2 &p) const {
        const double xx = x - p.x;
        const double yy = y - p.y;
        return xx * xx + yy * yy;
    }

    /** Distance between two points */
    double distance(const Point2 &p) const { return sqrt(distance2(p)); }

    /** Accessor for accessing the i-th coordinate */
    const double& operator [](size_t i) const { return (&x)[i]; }

    /** Iterator for accessing the coordinates */
    typedef Iterator<Point2, double> iterator;
    iterator begin() const { return iterator(this, 0); }
    iterator end() const { return iterator(this, size()); }

    /** Define the behaviour of string stream */
    friend std::ostream& operator <<(std::ostream &s, const Point2 &p) {
        return s << p.x << ' ' << p.y;
    }

    /** Return data as string */
    string to_str() const { stringstream ss; ss << (*this); return ss.str(); }

    void periodic(Point2 &pmax, Point2 &pmin){
        x = periodic_image(x, pmax.x, pmin.x);
        y = periodic_image(y, pmax.y, pmin.y);
    }

    double x, y;    ///< Cartesian coordinates
};

/** Common operations for 3-dimensional point and vector */
class Vector3Data {
public:
    Vector3Data() : x(0), y(0), z(0) {}
    Vector3Data(const double xx) : x(xx), y(xx), z(xx) {}
    Vector3Data(const double xx, const double yy, const double zz) : x(xx), y(yy), z(zz) {}
    Vector3Data(const Vector3Data& v) : x(v.x), y(v.y), z(v.z) {}

    /** Dimensionality of vector */
    int size() const { return 3; }

    Vector3Data& operator +=(const Point2 &p) { x += p.x, y += p.y; return *this; }

    /** Addition of two vectors */
    Vector3Data operator +(const Vector3Data &p) const { return Vector3Data(x+p.x, y+p.y, z+p.z); }
    Vector3Data& operator +=(const Vector3Data &p) { x += p.x, y += p.y, z += p.z; return *this; }

    /** Subtraction of two vectors */
    Vector3Data operator -(const Vector3Data &p) const { return Vector3Data(x-p.x, y-p.y, z-p.z); }
    Vector3Data& operator -=(const Vector3Data &p) { x -= p.x, y -= p.y, z -= p.z; return *this; }

    /** Multiplication of two vectors */
    Vector3Data operator *(const Vector3Data &p) const { return Vector3Data(x*p.x, y*p.y, z*p.z); }
    Vector3Data& operator *=(const Vector3Data &p) { x *= p.x, y *= p.y, z *= p.z; return *this; }

    /** Division of two vectors */
    Vector3Data operator /(const Vector3Data &p) const { return Vector3Data(x/p.x, y/p.y, z/p.z); }
    Vector3Data& operator /=(const Vector3Data &p) { x /= p.x, y /= p.y, z /= p.z; return *this; }

    /** Adding a scalar to the vector */
    Vector3Data operator +(const double &r) const { return Vector3Data(x + r, y + r, z + r); }
    Vector3Data& operator +=(const double &r) { x += r, y += r, z += r; return *this; }

    /** Subtracting a scalar from the vector */
    Vector3Data operator -(const double &r) const { return Vector3Data(x - r, y - r, z - r); }
    Vector3Data& operator -=(const double &r) { x -= r, y -= r, z -= r; return *this; }

    /** Multiplying vector with a scalar */
    Vector3Data operator *(const double &r) const { return Vector3Data(x * r, y * r, z * r); }
    Vector3Data& operator *=(const double &r) { x *= r, y *= r, z *= r; return *this; }

    /** Dividing vector with a scalar */
    Vector3Data operator /(const double &r) const { return Vector3Data(x / r, y / r, z / r); }
    Vector3Data& operator /=(const double &r) { x /= r, y /= r, z /= r; return *this; }

    /** Comparison operators */
    bool operator ==(const Vector3Data &p) const { return x == p.x && y == p.y && z == p.z; }
    bool operator ==(const double d) const { return x == d && y == d && z == d; }

    /** Accessor for accessing the i-th coordinate */
    const double& operator [](const size_t i) const { return (&x)[i]; }
    double& operator [](const size_t i) { return (&x)[i]; }

    /** Iterator for accessing the coordinates */
    typedef Iterator<Vector3Data, double> iterator;
    iterator begin() const { return iterator(this, 0); }
    iterator end() const { return iterator(this, size()); }

    /** Define the behaviour of string stream */
    friend std::ostream& operator <<(std::ostream &s, const Vector3Data &p) {
        return s << p.x << ' ' << p.y << ' ' << p.z;
    }

    /** Return data as string */
    string to_str() const { stringstream ss; ss << (*this); return ss.str(); }

    double x, y, z; ///< VectorData data
};

/** Basic operations with 3-dimensional point */
// TODO: For some reason not working properly
class Point3 : public Vector3Data {
public:
    /** Constructors of Point3 class */
    Point3() : Vector3Data() {}
    Point3(const double xx) : Vector3Data(xx) {}
    Point3(const double xx, const double yy, const double zz) : Vector3Data(xx, yy, zz) {}
    Point3(const double xx, const double yy) : Vector3Data(xx, yy, 0) {}
    Point3(const Vector3Data &v) : Vector3Data(v) {}
    Point3(const dealii::Point<3> &p) : Vector3Data(p[0], p[1], p[2]) {}

    /** Squared distance between two points */

    double distance2(const Point2 &p) const {
        const double xx = x - p.x;
        const double yy = y - p.y;
        return xx * xx + yy * yy;
    }

    double distance2(const Point3 &p) const {
        double xx = x - p.x;
        double yy = y - p.y;
        double zz = z - p.z;
        return xx * xx + yy * yy + zz * zz;
    }

    /** Distance between two points */
    double distance(const Point3 &p) const { return sqrt(distance2(p)); }

    /** Squared distance between two Point3-s taking into account the simulation cell periodicity.
     * Period == 0 in some direction gives the result without periodicity in that direction. */
    double periodic_distance2(const Point3 &p, const double period_x = 0,
            const double period_y = 0, const double period_z = 0) const {

        double dx = fabs(x - p.x);
        double dy = fabs(y - p.y);
        double dz = fabs(z - p.z);

        if (MODES.PERIODIC) {
            if (dx > period_x * 0.5) dx = period_x - dx; // apply periodic boundary condition in x-direction
            if (dy > period_y * 0.5) dy = period_y - dy; // apply periodic boundary condition in y-direction
            if (dz > period_z * 0.5) dz = period_z - dz; // apply periodic boundary condition in z-direction
        }

        return dx * dx + dy * dy + dz * dz;
    }

    void periodic(Point3 &pmax, Point3 &pmin){
        x = periodic_image(x, pmax.x, pmin.x);
        y = periodic_image(y, pmax.y, pmin.y);
        z = periodic_image(z, pmax.z, pmin.z);
    }
};

/** Basic operations with 3-dimensional vector */
class Vec3 : public Vector3Data {
public:
    /** Vec3 constructors */
    Vec3() : Vector3Data() {}
    Vec3(const double xx) : Vector3Data(xx) {}
    Vec3(const double xx, const double yy, const double zz) : Vector3Data(xx, yy, zz) {}
    Vec3(const Vector3Data &v) : Vector3Data(v) {}
    Vec3(const dealii::Tensor<1,3>& t) : Vector3Data(t[0], t[1], t[2])  {}

    /** Vector norm */
    double norm2() const { return x * x + y * y + z * z; }
    double norm() const { return sqrt(norm2()); }

    /** Normalize the vector */
    Vec3& normalize() {
        double n = norm();
        if (n > 0) {
            double factor = 1 / n;
            x *= factor, y *= factor, z *= factor;
        }
        return *this;
    }

    /** Dot product of two vectors */
    double dotProduct(const Vec3 &v) const { return x * v.x + y * v.y + z * v.z; }

    /** Cross product of two vectors */
    Vec3 crossProduct(const Vec3 &v) const {
        return Vec3(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x);
    }
};

/** Class to define basic operations with 4-dimensional vector */
class Vec4 {
public:
    /** Vec4 constructors */
    Vec4() : x(0), y(0), z(0), w(0) {}
    Vec4(const double xx) : x(xx), y(xx), z(xx), w(xx) {}
    Vec4(const double xx, double yy, double zz, double ww) : x(xx), y(yy), z(zz), w(ww) {}
    Vec4(const Vec3 &v, const double ww) : x(v.x), y(v.y), z(v.z), w(ww) {}

    /** Dimensionality of vector */
    int size() const { return 4; }

    /** Define equality of two vectors */
    bool operator ==(const Vec4 &v) const { return x == v.x && y == v.y && z == v.z && w == v.w; }

    /** Addition of two vectors */
    Vec4 operator +(const Vec4 &v) const { return Vec4(x + v.x, y + v.y, z + v.z, w + v.w); }
    Vec4& operator +=(const Vec4 &v) { x += v.x, y += v.y, z += v.z, w += v.w; return *this; }

    /** Subtraction of two vectors */
    Vec4 operator -(const Vec4 &v) const { return Vec4(x - v.x, y - v.y, z - v.z, w - v.w); }
    Vec4& operator -=(const Vec4 &v) { x -= v.x, y -= v.y, z -= v.z, w -= v.w; return *this; }

    /** Scalar multiplication of vector with a scalar */
    Vec4 operator *(const double &r) const { return Vec4(x * r, y * r, z * r, w * r); }
    Vec4& operator *=(const double &r) { x *= r, y *= r, z *= r, w *= r; return *this; }

    /** Scalar division of vector with a scalar */
    Vec4 operator /(const double &r) const { return Vec4(x / r, y / r, z / r, w / r); }
    Vec4& operator /=(const double &r) {  x /= r, y /= r, z /= r, w /= r; return *this; }

    /** Dot product of two vectors */
    double dotProduct(const Vec4 &v) const { return x * v.x + y * v.y + z * v.z + w * v.w; }

    /** Vector norm */
    double norm2() const { return x * x + y * y + z * z + w * w; }
    double norm() const { return sqrt(norm2()); }

    /** Function to normalize the vector */
    Vec4& normalize() {
        double n = norm();
        if (n > 0) {
            double factor = 1 / sqrt(n);
            x *= factor, y *= factor, z *= factor, w *= factor;
        }
        return *this;
    }

    /** Define access operator */
    const double& operator [](const size_t i) const { return (&x)[i]; }

    /** Attach iterator for accessing the vector components */
    typedef Iterator<Vec4, double> iterator;
    iterator begin() const { return iterator(this, 0); }
    iterator end() const { return iterator(this, size()); }

    /** Define the behaviour of string stream */
    friend std::ostream& operator <<(std::ostream &s, const Vec4 &v) {
        return s << v.x << ' ' << v.y << ' ' << v.z << ' ' << v.w;
    }

    /** Return data as string */
    string to_str() const { stringstream ss; ss << (*this); return ss.str(); }

    double x, y, z, w;    ///< Cartesian coordinates
};

/** Class to define elementary operations between atoms */
class Atom {
public:
    /** Constructors of Atom class */
    Atom() : id(0), point(0), marker(0) {}
    Atom(const int i, const Point3 &p, const int c) : id(i), point(p.x, p.y, p.z), marker(c) {}

    /** Comparison operator between two Atom-s */
    bool operator ==(const Atom &a) const {
        return id == a.id && point == a.point && marker == a.marker;
    }

    /** Define the behaviour of string stream */
    friend std::ostream& operator <<(std::ostream &s, const Atom &a) {
        return s << a.id << ' ' << a.point << ' ' << a.marker;
    }

    /** Functor for sorting atoms spatially in ascending order; used in CGAL Hilbert sorting */
    struct sort_spatial {
        typedef Atom Point_3;

        struct LessX {
            bool operator() (const Atom& a1, const Atom& a2) { return a1.point.x < a2.point.x; }
        };
        struct LessY {
            bool operator() (const Atom& a1, const Atom& a2) { return a1.point.y < a2.point.y; }
        };
        struct LessZ {
            bool operator() (const Atom& a1, const Atom& a2) { return a1.point.z < a2.point.z; }
        };
        struct ComputeX {
            double operator() (const Atom& a) { return a.point.x; }
        };
        struct ComputeY {
            double operator() (const Atom& a) { return a.point.y; }
        };
        struct ComputeZ {
            double operator() (const Atom& a) { return a.point.z; }
        };

        LessX less_x_3_object() const { return LessX(); }
        LessY less_y_3_object() const { return LessY(); }
        LessZ less_z_3_object() const { return LessZ(); }

        ComputeX compute_x_3_object() const { return ComputeX(); }
        ComputeY compute_y_3_object() const { return ComputeY(); }
        ComputeZ compute_z_3_object() const { return ComputeZ(); }
    };

    /** Functor for sorting atoms in ascending order by their ID */
    struct sort_id {
        inline bool operator() (const Atom& a1, const Atom& a2) { return a1.id < a2.id; }
    };

    /** Functor for sorting atoms in ascending order by their marker */
    struct sort_marker_up {
        inline bool operator() (const Atom& a1, const Atom& a2) { return a1.marker < a2.marker; }
    };

    /** Functor for sorting atoms in descending order by their marker */
    struct sort_marker_down {
        inline bool operator() (const Atom& a1, const Atom& a2) { return a1.marker > a2.marker; }
    };

    /** Functor for sorting atoms in ascending order by their x, y, z or radial coordinate */
    struct sort_up {
        int I;           //!< coordinate along which the atoms are sorted; 0=x, 1=y, 2=z
        sort_up(int d) : I(d) {}
        inline bool operator() (const Atom& a1, const Atom& a2) { return a1.point[I] < a2.point[I]; }
    };

    /** Functor for sorting atoms in descending order by their x, y, z or radial coordinate */
    struct sort_down {
        int I;           //!< coordinate along which the atoms are sorted; 0=x, 1=y, 2=z
        sort_down(int d) : I(d) {}
        inline bool operator() (const Atom& a1, const Atom& a2) { return a1.point[I] > a2.point[I]; }
    };

    /** Functor for sorting atoms in ascending order first by their x, y, z or radial coordinate
     * and then in ascending order by another x, y, z or radial coordinate */
    struct sort_up2 {
        int I1, I2;      //!< coordinates along which the atoms are sorted; 0=x, 1=y, 2=z
        sort_up2(int coord1, int coord2) : I1(coord1), I2(coord2) {}
        inline bool operator() (const Atom& a1, const Atom& a2) {
            return (a1.point[I1] < a2.point[I1]) || ((a1.point[I1] == a2.point[I1]) && (a1.point[I2] < a2.point[I2]));
        }
    };

    /** Functor for sorting atoms in descending order first by their x, y, z or radial coordinate
     * and then in descending order by another x, y, z or radial coordinate */
    struct sort_down2 {
        int I1, I2;      //!< coordinates along which the atoms are sorted; 0=x, 1=y, 2=z
        sort_down2(int coord1, int coord2) : I1(coord1), I2(coord2) {}
        inline bool operator() (const Atom& a1, const Atom& a2) {
            return (a1.point[I1] > a2.point[I1]) || ((a1.point[I1] == a2.point[I1]) && (a1.point[I2] > a2.point[I2]));
        }
    };

    /** Functor for sorting atoms in ascending order first by their x, y, z or radial coordinate
     * and then in descending order by another x, y, z or radial coordinate */
    struct sort_up_down {
        int I1, I2;      //!< coordinates along which the atoms are sorted; 0=x, 1=y, 2=z
        sort_up_down(int coord1, int coord2) : I1(coord1), I2(coord2) {}
        inline bool operator() (const Atom& a1, const Atom& a2) {
            return (a1.point[I1] < a2.point[I1]) || ((a1.point[I1] == a2.point[I1]) && (a1.point[I2] > a2.point[I2]));
        }
    };

    /** Functor for sorting atoms in descending order first by their x, y, z or radial coordinate
     * and then in ascending order by another x, y, z or radial coordinate */
    struct sort_down_up {
        int I1, I2;      //!< coordinates along which the atoms are sorted; 0=x, 1=y, 2=z
        sort_down_up(int coord1, int coord2) : I1(coord1), I2(coord2) {}
        inline bool operator() (const Atom& a1, const Atom& a2) {
            return (a1.point[I1] > a2.point[I1]) || ((a1.point[I1] == a2.point[I1]) && (a1.point[I2] < a2.point[I2]));
        }
    };

    int id;
    Point3 point;
    int marker;
};

/** Class to hold solution data and its operations */
class Solution {
public:

    /** Constructors of Solution */
    Solution() : id(0), vector(Vec3(0)), norm(0), scalar(0) {}
    Solution(const double d) : id(0), vector(Vec3(d)), norm(d), scalar(d) {}
    Solution(const Vec3& v) : id(0), vector(v), norm(v.norm()), scalar(0) {}
    Solution(const Vec3& v, const double s) : id(0), vector(v), norm(v.norm()), scalar(s) {}
    Solution(const Vec3& v, const double n, const double s) : id(0), vector(v), norm(n), scalar(s) {}

    /** Define the behaviour of string stream */
    friend std::ostream& operator <<(std::ostream &ss, const Solution &sol) {
        return ss << sol.vector << ' ' << sol.norm << ' ' << sol.scalar;
    }

    /** Functors used to sort vector of Solution into same order as vector of Atom */
    struct sort_up {
        inline bool operator() (const Solution& lh, const Solution& rh) { return lh.id < rh.id; }
    };
    struct sort_down {
        inline bool operator() (const Solution& lh, const Solution& rh) { return lh.id > rh.id; }
    };

    int id;
    Vec3 vector;
    double norm;
    double scalar;
};

/** Super particles for PIC calculations */
class SuperParticle {
public:
    SuperParticle() : pos(Point3(0)), vel(Vec3(0)), cell(0) {}
    SuperParticle(const Point3 &p, const Vec3 &v, const int c) : pos(p), vel(v), cell(c) {}

    Point3 pos;    ///< particle position [Å]
    Vec3 vel;      ///< particle velocity [Å/fs]
    int cell;      ///< mesh hexahedron ID, -1 if outside the domain

    /** Define sorting order */
    bool operator < (const SuperParticle &p) { return (cell < p.cell); }

    /** Define the behaviour of string stream */
    friend std::ostream& operator <<(std::ostream &s, const SuperParticle &p) {
        return s << p.pos << ' ' << p.vel << ' ' << p.cell;
    }
};

/** Template class for finite element cell */
template <size_t dim>
class SimpleCell {
public:

    /** SimpleCell constructors */
    SimpleCell() { std::fill_n(node, dim, int()); }
    SimpleCell(const unsigned int &nn) { std::fill_n(node, dim, nn); }

    /** Dimensionality of cell */
    int size() const { return dim; }

    /** Check if SimpleCell contains the node of interest */
    bool operator ==(const unsigned int &t) const {
        for (unsigned int n : node) if (n == t) return true;
        return false;
    }

    /** Check if SimpleCell does not contain the node of interest */
    bool operator !=(const unsigned int &t) const {
        for (unsigned int n : node) if (n == t) return false;
        return true;
    }

    /** Check if two SimpleCells are equal to one-another */
    bool operator ==(const SimpleCell<dim> &sc) const {
        for (unsigned int n : node) if (sc != n) return false;
        return true;
    }

    /** Check if two SimpleCells are not equal to one-another */
    bool operator !=(const SimpleCell<dim> &sc) const {
        for (unsigned int n : node) if (sc == n) return false;
        return true;
    }

    /** Less than, bigger than, less than or equal, bigger than or equal operators */
    bool operator <(const unsigned int &t) const {
        for (unsigned int n : node) if (n >= t) return false;
        return true;
    }
    bool operator >(const unsigned int &t) const {
        for (unsigned int n : node) if (n <= t) return false;
        return true;
    }
    bool operator <=(const unsigned int &t) const {
        for (unsigned int n : node) if (n > t) return false;
        return true;
    }
    bool operator >=(const unsigned int &t) const {
        for (unsigned int n : node) if (n < t) return false;
        return true;
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

    /** Accessor for accessing the i-th node */
    const unsigned int& operator [](size_t i) const {
        require(i >= 0 && i < dim, "Invalid index: " + to_string(i));
        return node[i];
    }

    /** Iterator to access the cell nodes */
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
    SimpleNode(const SimpleCell<1> &s) { node[0] = s.node[0]; }
};

/** Edge class without Point data */
class SimpleEdge: public SimpleCell<2> {
public:
    /** SimpleEdge constructors */
    SimpleEdge() : SimpleCell<2>() {}
    SimpleEdge(const unsigned int &n1) : SimpleCell<2>(n1) {}
    SimpleEdge(const unsigned int &n1, const unsigned int &n2) { node[0] = n1; node[1] = n2; }
    SimpleEdge(const SimpleCell<2> &s) { node[0] = s.node[0]; node[1] = s.node[1]; }

    /** Get i-th vertex of the edge */
    SimpleNode vert(const int i) const {
        require(i >= 0 && i <= 1, "Invalid index: " + to_string(i));
        return SimpleNode(node[i]);
    }

    /** Check whether two edges share a node */
    bool neighbor(const SimpleEdge& s) const {
        const bool b1 = s == node[0];
        const bool b2 = s == node[1];
        return b1 + b2 == 1;
    }

    /** Check whether edge contains a node */
    bool contains(const SimpleNode& s) const {
        return s == node[0] || s == node[1];
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
    SimpleFace(const SimpleCell<3> &s) {
        std::copy( std::begin(s.node), std::end(s.node), std::begin(node) );
    }

    /** Get i-th edge of the face */
    SimpleEdge edge(const int i) const {
        require(i >= 0 && i <= 2, "Invalid index: " + to_string(i));
        if (i == 0) return SimpleEdge(node[0], node[1]);
        if (i == 1) return SimpleEdge(node[0], node[2]);
        else        return SimpleEdge(node[1], node[2]);
    }

    /** Check whether two triangles share an edge */
    bool edge_neighbor(const SimpleFace& s) const {
        int sum = s == node[0];
        sum += s == node[1];
        sum += s == node[2];
        return sum == 2;
    }

    /** Check whether two triangles share a node or an edge */
    bool node_neighbor(const SimpleFace& s) const {
        int sum = s == node[0];
        sum += s == node[1];
        sum += s == node[2];
        return sum == 1 || sum == 2 ;
    }

    /** Check whether triangle contains an edge */
    bool contains(const SimpleEdge& s) const {
        const bool b1 = s == node[0];
        const bool b2 = s == node[1];
        const bool b3 = s == node[2];
        return b1 + b2 + b3 == 2;
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
    SimpleElement(const SimpleCell<4> &s) {
        std::copy( std::begin(s.node), std::end(s.node), std::begin(node) );
    }

    /** Get i-th edge of the element */
    SimpleEdge edge(const int i) const {
        require(i >= 0 && i <= 5, "Invalid index: " + to_string(i));
        if (i == 0) return SimpleEdge(node[0], node[1]);
        if (i == 1) return SimpleEdge(node[0], node[2]);
        if (i == 2) return SimpleEdge(node[0], node[3]);
        if (i == 3) return SimpleEdge(node[1], node[2]);
        if (i == 4) return SimpleEdge(node[1], node[3]);
        else        return SimpleEdge(node[2], node[3]);
    }

    /** Get i-th face of the element */
    SimpleFace face(const int i) const {
        require(i >= 0 && i <= 3, "Invalid index: " + to_string(i));
        if (i == 0) return SimpleFace(node[0], node[1], node[2]);
        if (i == 1) return SimpleFace(node[0], node[1], node[3]);
        if (i == 2) return SimpleFace(node[0], node[2], node[3]);
        else        return SimpleFace(node[1], node[2], node[3]);
    }

    /** Check whether two tetrahedra share a triangle */
    bool neighbor(const SimpleElement& s) const {
        int sum = s == node[0];
        sum += s == node[1];
        sum += s == node[2];
        sum += s == node[3];
        return sum == 3;
    }

    /** Check whether tetrahedron contains a triangle */
    bool contains(const SimpleFace& s) const {
        const bool b1 = s == node[0];
        const bool b2 = s == node[1];
        const bool b3 = s == node[2];
        const bool b4 = s == node[3];
        return b1 + b2 + b3 + b4 == 3;
    }
};

/** Quadrangle class without Point data */
class SimpleQuad: public SimpleCell<4> {
public:
    /** SimpleQuad constructors */
    SimpleQuad() : SimpleCell<4>() {}
    SimpleQuad(const unsigned int &n1) : SimpleCell<4>(n1) {}
    SimpleQuad(const unsigned int &n1, const unsigned int &n2, const unsigned int &n3, const unsigned int &n4) {
        node[0] = n1; node[1] = n2; node[2] = n3; node[3] = n4;
    }
    SimpleQuad(const SimpleCell<4> &s) {
        std::copy( std::begin(s.node), std::end(s.node), std::begin(node) );
    }
};

/** Hexahedron class without Point data */
class SimpleHex: public SimpleCell<8> {
public:
    /** SimpleHex constructors */
    SimpleHex() : SimpleCell<8>() {}
    SimpleHex(const unsigned int &n1) : SimpleCell<8>(n1) {}
    SimpleHex(const unsigned int &n1, const unsigned int &n2, const unsigned int &n3, const unsigned int &n4,
            const unsigned int &n5, const unsigned int &n6, const unsigned int &n7, const unsigned int &n8) {
        node[0] = n1; node[1] = n2; node[2] = n3; node[3] = n4;
        node[4] = n5; node[5] = n6; node[6] = n7; node[7] = n8;
    }
    SimpleHex(const SimpleCell<8> &s) {
        std::copy( std::begin(s.node), std::end(s.node), std::begin(node) );
    }
};

/** Quadratic triangle class without Point data */
class QuadraticTri: public SimpleCell<6> {
public:
    QuadraticTri() : SimpleCell<6>() {}
    QuadraticTri(const unsigned &n1) : SimpleCell<6>(n1) {}

    QuadraticTri(const unsigned &n1, const unsigned &n2, const unsigned &n3,
            const unsigned &n4, const unsigned &n5, const unsigned &n6) {
        node[0] = n1; node[1] = n2; node[2] = n3; node[3] = n4; node[4] = n5; node[5] = n6;
    }

    QuadraticTri(const SimpleFace &s, const unsigned &n4, const unsigned &n5, const unsigned &n6) {
        std::copy( std::begin(s.node), std::end(s.node), std::begin(node) );
        node[3] = n4, node[4] = n5; node[5] = n6;
    }

    QuadraticTri(const SimpleCell<6> &s) {
        std::copy( std::begin(s.node), std::end(s.node), std::begin(node) );
    }
};

/** Quadratic tetrahedron class without Point data */
class QuadraticTet: public SimpleCell<10> {
public:
    QuadraticTet() : SimpleCell<10>() {}
    QuadraticTet(const unsigned &n1) : SimpleCell<10>(n1) {}

    QuadraticTet(const unsigned &n1, const unsigned &n2, const unsigned &n3, const unsigned &n4,
            const unsigned &n5, const unsigned &n6, const unsigned &n7,
            const unsigned &n8, const unsigned &n9, const unsigned &n10) {
        node[0] = n1; node[1] = n2; node[2] = n3; node[3] = n4; node[4] = n5;
        node[5] = n6; node[6] = n7; node[7] = n8; node[8] = n9; node[9] = n10;
    }

    QuadraticTet(const SimpleElement &s, const unsigned &n5, const unsigned &n6, const unsigned &n7,
            const unsigned &n8, const unsigned &n9, const unsigned &n10) {
        std::copy( std::begin(s.node), std::end(s.node), std::begin(node) );
        node[4] = n5; node[5] = n6; node[6] = n7; node[7] = n8; node[8] = n9; node[9] = n10;
    }

    QuadraticTet(const SimpleCell<10> &s) {
        std::copy( std::begin(s.node), std::end(s.node), std::begin(node) );
    }
};

} // namaspace femocs

#endif /* PRIMITIVES_H_ */
