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

template <typename T, size_t dim>
class VectorData {
public:

    /** SimpleCell constructors */
    VectorData() { std::fill_n(node, dim, T()); }
    VectorData(const T &nn) { std::fill_n(node, dim, nn); }

    /** Dimensionality of vector */
    int size() const { return dim; }

    /** Check if one of the nodes equals to the one of interest */
    bool operator ==(const T &t) const {
        for (T n : node) if (n == t) return true;
        return false;
    }

    bool operator ==(const VectorData<T, dim> &t) const {
        for (size_t i = 0; i < dim; ++i) if ((*this)[i] != t[i]) return false;
        return true;
    }

    /** Define the behaviour of string stream */
    friend std::ostream& operator <<(std::ostream &s, const VectorData<T, dim> &t) {
        for (T nn : t.node) s << nn << ' ';
        return s;
    }

    /** Return data as string */
    string to_str() const { stringstream ss; ss << this; return ss.str(); }

    /** Accessor for accessing the i-th node */
    const T& operator [](const size_t i) const {
        require(i >= 0 && i < dim, "Invalid index: " + to_string(i));
        return node[i];
    }

    /** Iterator to access the cell nodes */
    typedef Iterator<VectorData<T, dim>, T> iterator;
    iterator begin() const { return iterator(this, 0); }
    iterator end() const { return iterator(this, dim); }

    T node[dim]; ///< VectorData data
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

    /** Subtraction of two points */
    Point2 operator -(const Point2 &p) const {
        return Point2(x - p.x, y - p.y);
    }

    /** Comparison operators between two points */
    bool operator ==(const Point2 &p) const { return x == p.x && y == p.y; }

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

    double x, y;    ///< Cartesian coordinates
};

/** Class to define basic operations with 3-dimensional points */
class Point3 {
public:
    /** Constructors of Point3 class */
    Point3() : x(0), y(0), z(0) {}
    Point3(const double xx) : x(xx), y(xx), z(xx) {}
    Point3(const double xx, const double yy, const double zz) : x(xx), y(yy), z(zz) {}
    Point3(const Point3& p) : x(p.x), y(p.y), z(p.z) {}

    /** Dimensionality of point */
    int size() const { return 3; }

    /** Squared distance between two Point3-s */
    double distance2(const Point3 &p) const {
        double xx = x - p.x;
        double yy = y - p.y;
        double zz = z - p.z;
        return xx * xx + yy * yy + zz * zz;
    }

    /** Distance between two Point3-s */
    double distance(const Point3 &p) const { return sqrt(distance2(p)); }

    /** Squared distance between a Point3 and dealii::Point<3> */
    double distance2(dealii::Point<3> &p) const {
        double xx = x - p[0];
        double yy = y - p[1];
        double zz = z - p[2];
        return xx * xx + yy * yy + zz * zz;
    }

    /** Distance between a Point3 and dealii::SimpleCell_T<3> */
    double distance(dealii::Point<3> &p) const { return sqrt(distance2(p)); }

    /** Squared distance between two Point3-s taking into account the simulation cell periodicity.
     * Period == 0 in some direction gives the result without periodicity in that direction. */
    double periodic_distance2(const Point3 &p, const double period_x = 0,
            const double period_y = 0, const double period_z = 0) const {

        double dx = fabs(x - p.x);
        double dy = fabs(y - p.y);
        double dz = fabs(z - p.z);

//        dx = min(dx, fabs(dx - period_x)); // apply periodic boundary condition in x-direction
//        dy = min(dy, fabs(dy - period_y)); // apply periodic boundary condition in y-direction
//        dz = min(dz, fabs(dz - period_z)); // apply periodic boundary condition in z-direction

        if (MODES.PERIODIC) {
            if (dx > period_x * 0.5) dx = period_x - dx; // apply periodic boundary condition in x-direction
            if (dy > period_y * 0.5) dy = period_y - dy; // apply periodic boundary condition in y-direction
            if (dz > period_z * 0.5) dz = period_z - dz; // apply periodic boundary condition in z-direction
        }

        return dx * dx + dy * dy + dz * dz;
    }

    /** Addition of two points */
    Point3 operator +(const Point3 &p) const { return Point3(x+p.x, y+p.y, z+p.z); }
    Point3& operator +=(const Point2 &p) { x += p.x, y += p.y; return *this; }
    Point3& operator +=(const Point3 &p) { x += p.x, y += p.y, z += p.z; return *this; }

    /** Subtraction of two points */
    Point3 operator -(const Point3 &p) const { return Point3(x-p.x, y-p.y, z-p.z); }
    Point3& operator -=(const Point3 &p) { x -= p.x, y -= p.y, z -= p.z; return *this; }

    /** Multiplying point with a scalar */
    Point3 operator *(const double &r) const { return Point3(x * r, y * r, z * r); }
    Point3& operator *=(const double &r) { x *= r, y *= r, z *= r; return *this; }

    /** Dividing point with a scalar */
    Point3 operator /(const double &r) const { return Point3(x / r, y / r, z / r); }
    Point3& operator /=(const double &r) { x /= r, y /= r, z /= r; return *this; }

    /** Comparison operators between two 3D points */
    bool operator ==(const Point3 &p) const { return x == p.x && y == p.y && z == p.z; }
    bool operator ==(const dealii::Point<3> &p) const {
        return x == p[0] && y == p[1] && z == p[2];
    }

    /** Accessor for accessing the i-th coordinate */
    const double& operator [](const size_t i) const { return (&x)[i]; }

    /** Iterator for accessing the coordinates */
    typedef Iterator<Point3, double> iterator;
    iterator begin() const { return iterator(this, 0); }
    iterator end() const { return iterator(this, size()); }

    /** Define the behaviour of string stream */
    friend std::ostream& operator <<(std::ostream &s, const Point3 &p) {
        return s << p.x << ' ' << p.y << ' ' << p.z;
    }

    /** Return data as string */
    string to_str() const { stringstream ss; ss << (*this); return ss.str(); }

    double x, y, z; ///< Cartesian coordinates
};

/** Class to define basic operations with 3-dimensional vector */
class Vec3 {
public:
    /** Vec3 constructors */
    Vec3() : x(0), y(0), z(0) {}
    Vec3(const double xx) : x(xx), y(xx), z(xx) {}
    Vec3(const double xx, const double yy, const double zz) : x(xx), y(yy), z(zz) {}
    Vec3(const Point3 &p) : x(p.x), y(p.y), z(p.z) {}

    /** Dimensionality of vector */
    int size() const { return 3; }

    /** Compare vector with a scalar */
    bool operator ==(const double d) const { return x == d && y == d && z == d; }
    
    /** Compare one vector with another */
    bool operator ==(const Vec3 &v) const { return x == v.x && y == v.y && z == v.z; }

    /** Addition of two vectors */
    Vec3 operator +(const Vec3 &v) const { return Vec3(x + v.x, y + v.y, z + v.z); }
    Vec3& operator +=(const Vec3 &v) { x += v.x, y += v.y, z += v.z; return *this; }

    /** Subtraction of two vectors */
    Vec3 operator -(const Vec3 &v) const { return Vec3(x - v.x, y - v.y, z - v.z); }
    Vec3& operator -=(const Vec3 &v) { x -= v.x, y -= v.y, z -= v.z; return *this; }

    /** Multiplication of two vectors */
    Vec3 operator *(const Vec3 &v) const { return Vec3(x * v.x, y * v.y, z * v.z); }
    Vec3& operator *=(const Vec3 &v) { x *= v.x, y *= v.y, z *= v.z; return *this; }

    /** Scalar multiplication of vector with a scalar */
    Vec3 operator *(const double &r) const { return Vec3(x * r, y * r, z * r); }
    Vec3& operator *=(const double &r) { x *= r, y *= r, z *= r; return *this; }

    /** Division of two vectors */
    Vec3 operator /(const Vec3 &v) const { return Vec3(x / v.x, y / v.y, z / v.z); }
    Vec3& operator /=(const Vec3 &v) { x /= v.x, y /= v.y, z /= v.z; return *this; }

    /** Scalar division of vector with a scalar */
    Vec3 operator /(const double &r) const { return Vec3(x / r, y / r, z / r); }
    Vec3& operator /=(const double &r) { x /= r, y /= r, z /= r; return *this; }

    /** Dot product and cross product of two vectors */
    double dotProduct(const Vec3 &v) const { return x * v.x + y * v.y + z * v.z; }
    Vec3 crossProduct(const Vec3 &v) const {
        return Vec3(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x);
    }

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

    /**
     Define access operators or accessors.
     The Vec3 coordinates can be accessed that way v[0], v[1], v[2], rather than v.x, v.y, v.z.
     This is useful in loops: the coordinates can be accessed with the loop index (e.g. v[i]).
     */
    const double& operator [](const size_t i) const { return (&x)[i]; }

    /** Attach iterator for accessing the vector components */
    typedef Iterator<Vec3, double> iterator;
    iterator begin() const { return iterator(this, 0); }
    iterator end() const { return iterator(this, size()); }

    /** Define the behaviour of string stream */
    friend std::ostream& operator <<(std::ostream &s, const Vec3 &v) {
        return s << v.x << ' ' << v.y << ' ' << v.z;
    }

    /** Return data as string */
    string to_str() const { stringstream ss; ss << (*this); return ss.str(); }

    double x, y, z; ///< Cartesian coordinates
};

/** Class to define basic operations with 4-dimensional vector */
class Vec4 {
public:
    /** Vec4 constructors */
    Vec4() : x(0), y(0), z(0), w(0) {}
    Vec4(const double xx) : x(xx), y(xx), z(xx), w(xx) {}
    Vec4(const double xx, double yy, double zz, double ww) : x(xx), y(yy), z(zz), w(ww) {}
    Vec4(const Vec3 &v, const double ww) : x(v.x), y(v.y), z(v.z), w(ww) {}
    Vec4(const Point3 &p, const double ww) : x(p.x), y(p.y), z(p.z), w(ww) {}

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

/** Template class for finite element cell */
template <size_t dim>
class SimpleCell {
public:

    /** SimpleCell constructors */
    SimpleCell() { std::fill_n(node, dim, int()); }
    SimpleCell(const unsigned int &nn) { std::fill_n(node, dim, nn); }

    /** Dimensionality of cell */
    int size() const { return dim; }

    /** Check if one of the nodes equals to the one of interest */
    bool operator ==(const unsigned int &t) const {
        for (unsigned int n : node) if (n == t) return true;
        return false;
    }

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

    /** Compare node with scalar */
    bool operator ==(const unsigned int n) const {
        return node[0] == n;
    }
};

/** Edge class without Point data */
class SimpleEdge: public SimpleCell<2> {
public:
    /** SimpleEdge constructors */
    SimpleEdge() : SimpleCell<2>() {}
    SimpleEdge(const unsigned int &n1) : SimpleCell<2>(n1) {}
    SimpleEdge(const unsigned int &n1, const unsigned int &n2) { node[0] = n1; node[1] = n2; }
    SimpleEdge(const SimpleCell<2> &s) { node[0] = s.node[0]; node[1] = s.node[1]; }

    /** Check whether edge contains node */
    bool operator ==(const unsigned int n) const {
        return node[0] == n || node[1] == n;
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

    /** Check whether face contains node */
    bool operator ==(const unsigned int n) const {
        return node[0] == n || node[1] == n || node[2] == n;
    }

    /** Get i-th edge of the face */
    SimpleEdge edge(const unsigned int i) const {
        require(i >= 0 && i <= 2, "Invalid index: " + to_string(i));
        if (i == 0) return SimpleEdge(node[0], node[1]);
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
    SimpleElement(const SimpleCell<4> &s) {
        std::copy( std::begin(s.node), std::end(s.node), std::begin(node) );
    }
    /** Check whether element contains node */
    bool operator ==(const unsigned int n) const {
        return node[0] == n || node[1] == n || node[2] == n || node[3] == n;
    }

    /** Get i-th edge of the element */
    SimpleEdge edge(const unsigned int i) const {
        require(i >= 0 && i <= 5, "Invalid index: " + to_string(i));
        if (i == 0) return SimpleEdge(node[0], node[1]);
        if (i == 1) return SimpleEdge(node[0], node[2]);
        if (i == 2) return SimpleEdge(node[0], node[3]);
        if (i == 3) return SimpleEdge(node[1], node[2]);
        if (i == 4) return SimpleEdge(node[1], node[3]);
        else        return SimpleEdge(node[2], node[3]);
    }

    /** Get i-th face of the element */
    SimpleFace face(const unsigned int i) const {
        require(i >= 0 && i <= 3, "Invalid index: " + to_string(i));
        if (i == 0) return SimpleFace(node[0], node[1], node[2]);
        if (i == 1) return SimpleFace(node[0], node[1], node[3]);
        if (i == 2) return SimpleFace(node[1], node[2], node[3]);
        else        return SimpleFace(node[2], node[0], node[1]);
    }
};

/** Hexahedron class without Point data */
class SimpleHex: public SimpleCell<8> {
public:
    /** SimpleEdge constructors */
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

} // namaspace femocs

#endif /* PRIMITIVES_H_ */
