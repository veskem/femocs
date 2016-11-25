/*
 * Coarsener.h
 *
 *  Created on: 7.10.2016
 *      Author: veske
 */

#ifndef COARSENER_H_
#define COARSENER_H_

#include "Primitives.h"
#include "Medium.h"
#include "Config.h"

using namespace std;
namespace femocs {

/** Virtual class to get data needed to coarsen surface atoms */
class Coarsener {
public:
    Coarsener();
    Coarsener(const Point3 &origin, double radius, double A, double r0_min=0, double r0_max=1e20);

    virtual ~Coarsener() {};

    /** Points INSIDE the region will be coarsened */
    virtual void pick_cutoff(const Point3 &point) {
        if (in_region(point))
            cutoff2 = get_const_cutoff();
        else
            cutoff2 = get_inf_cutoff();
    }

    bool nearby(const Point3 &p1, const Point3 &p2) const {
        return p1.distance2(p2) <= cutoff2;
    }

    /** Number of points written to the vtk file */
    virtual int get_n_points() const { return 0; }

    /** Number of polygons written to the vtk file */
    virtual int get_n_polygons() const { return 0; }

    /** Get the points for writing to vtk file */
    virtual vector<Point3> get_points() { return vector<Point3>{}; }

    /** Get nodes of polygons for writing to vtk file */
    virtual vector<vector<int>> get_polygons() { return vector<vector<int>>{vector<int>{}}; }

protected:
    Point3 origin3d;  //!< centre of the coarsener
    double cutoff2;   //!< squared cut off radius
    double radius;    //!< radius of the system
    double radius2;   //!< squared radius of the system
    double r0_min;    //!< cut off radius for the points on the edge of the system
    double r0_max;    //!< cut off radius for the points far away from of the system
    double A;         //!< coarsening factor

    const int n_nodes_per_circle = 50; //!< number of nodes in vtk circle-polygon
    const int n_nodes_per_line = 2;    //!< number of nodes in vtk line-polygon
    const int n_circles = 0;           //!< number of circles written to vtk file
    const int n_lines = 0;             //!< number of lines written to vtk file

    /** Get constant squared cut off radius */
    double get_const_cutoff() const { return r0_min * r0_min; }

    /** Get cut off radius that increases with distance from origin */
    double get_increasing_cutoff(const Point3 &point) const {
        double cutoff = max(0.0, origin3d.distance(point) - radius);
        cutoff = min( r0_max, A * sqrt(cutoff) + r0_min );
        return cutoff * cutoff;
    }

    /** Get cut off radius that is smaller than any possible distance between atoms */
    double get_inf_cutoff() const { return -1e20; }

    /** Point in region? */
    virtual inline bool in_region(const Point3 &point) const { return false; }
};


/** Class to coarsen whole area uniformly */
class ConstCoarsener: public Coarsener {
public:
    ConstCoarsener();
    ConstCoarsener(double r0_min);

private:
    /** Point in anywhere? */
    inline bool in_region(const Point3 &point) const { return true; }
};

/** Class to coarsen surface outside one infinite vertical cylinder */
class FlatlandCoarsener: public Coarsener {
public:
    FlatlandCoarsener();
    FlatlandCoarsener(const Point3 &origin, double radius, double A, double r0_min=0, double r0_max=1e20);

    /** Points OUTSIDE the region will be coarsened */
    void pick_cutoff(const Point3 &point) {
        if (in_region(point))
            cutoff2 = get_increasing_cutoff(point);
        else
            cutoff2 = get_inf_cutoff();
    }

private:
    Point2 origin2d;

    /** Point outside infinitely high vertical cylinder? */
    inline bool in_region(const Point3 &point) const {
        return origin2d.distance2(Point2(point.x, point.y)) > radius2;
    }
};

/** Class to coarsen surface inside one infinite vertical cylinder */
class CylinderCoarsener: public Coarsener {
public:
    CylinderCoarsener();
    CylinderCoarsener(const Point2 &base, double radius, double r0_cylinder=0);

    /** Points INSIDE the nanotip will be coarsened */
    void pick_cutoff(const Point3 &point) {
        if (in_region(point))
            cutoff2 = get_const_cutoff();
        else
            cutoff2 = get_inf_cutoff();
    }

    vector<Point3> get_points();
    vector<vector<int>> get_polygons();
    int get_n_points() const { return n_lines * n_nodes_per_line + n_circles * n_nodes_per_circle; }
    int get_n_polygons() const { return n_lines + n_circles; };

protected:
    /** Point in infinitely high vertical cylinder? */
    inline bool in_region(const Point3 &point) const {
        return origin2d.distance2(Point2(point.x, point.y)) <= radius2;
    }

private:
    const int n_circles = 1;
    const int n_lines = 4;
    Point2 origin2d;
};

/** Class to coarsen surface inside one infinite vertical nanotip */
class NanotipCoarsener: public Coarsener {
public:
    NanotipCoarsener();
    NanotipCoarsener(const Point3 &apex, double radius, double A, double r0_apex, double r0_cylinder);

    /** Points INSIDE the nanotip will be coarsened */
    void pick_cutoff(const Point3 &point) {
        if (in_region(point))
            cutoff2 = get_increasing_cutoff(point);
        else
            cutoff2 = get_inf_cutoff();
    }

    vector<Point3> get_points();
    vector<vector<int>> get_polygons();
    int get_n_points() const { return n_lines * n_nodes_per_line + n_circles * n_nodes_per_circle; }
    int get_n_polygons() const { return n_lines + n_circles; };

protected:
    /** Point in infinitely high vertical cylinder? */
    inline bool in_region(const Point3 &point) const {
        return origin2d.distance2(Point2(point.x, point.y)) <= radius2;
    }

private:
    const int n_circles = 4;
    const int n_lines = 6;
    Point2 origin2d;
};

/** Class to coarsen surface inside one infinite tilted nanotip */
class TiltedNanotipCoarsener: public NanotipCoarsener {
public:
    TiltedNanotipCoarsener();
    TiltedNanotipCoarsener(const Point3 &apex, const Point3 &base, double radius, double A, double r0_apex, double r0_cylinder);

private:
    const int n_circles = 4;
    const int n_lines = 6;
    Vec3 bottom, axis;
    double height2;

    /** Point in infinite tilted cylinder? */
    inline bool in_region(const Point3 &point) const {
        Vec3 testpoint(point.x, point.y, point.z);  // convert point to vector
        Vec3 bot2point = testpoint - bottom;        // vector from centre of bottom face to test point

        double dot = axis.dotProduct(bot2point);

        // point between the ends of cylinder?
        //        if( dot < 0.0 || dot > height2 )
        //            return false;

        // check the squared distance to the cylinder axis
        return (bot2point.norm2() - dot*dot / height2) <= radius2;
    }
};


/** Class for binding together coarseners */
class Coarseners {
public:
    /** Coarseners constructor */
    Coarseners() : r0_inf(0), zmean(0) {}

    /** Append coarsener to the array of all the coarseners */
    void attach_coarsener(shared_ptr<Coarsener> c) {
        coarseners.push_back(c);
    }

    /** Pick cut off radiuses for coarseners */
    void pick_cutoff(const Point3 &point) {
        for (auto &c : coarseners)
            c->pick_cutoff(point);
    }

    /** Specify the collective action of coarseners */
    bool nearby(const Point3 &p1, const Point3 &p2) {
        bool near = false;
        for (auto &c : coarseners)
            near |= c->nearby(p1, p2);

        return near;
    }

    /** Write the contours of coarseners to file in .vtk format */
    void write(const string &file_name);

    /** Generate coarseners for one nanotip system */
    void generate(Medium &medium, const double radius, const Config::CoarseFactor &cf, const double latconst);

    double r0_inf;
    double zmean;

private:
    vector<shared_ptr<Coarsener>> coarseners;
};

} // namespace femocs

#endif /* COARSENER_H_ */
