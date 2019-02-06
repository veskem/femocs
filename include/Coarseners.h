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
#include "FileWriter.h"
#include <memory>

using namespace std;
namespace femocs {

/** Virtual class to get data needed to coarsen surface atoms */
class Coarsener {
public:
    Coarsener();
    Coarsener(const Point3 &origin, const double exp, const double radius,
            const double A, const double r0_min=0, const double r0_max=1e20);

    virtual ~Coarsener() {};

    /** Points INSIDE the region will be coarsened */
    virtual void pick_cutoff(const Point3 &point) {
        if (in_region(point))
            cutoff2 = get_const_cutoff();
        else
            cutoff2 = get_inf_cutoff();
    }

    /** Determine if two points are within their cut-off radius */
    bool nearby(const Point3 &p1, const Point3 &p2) const {
        return p1.distance2(p2) <= cutoff2;
    }

    /** Number of points written to the vtk file */
    virtual int get_n_points() const { return 0; }

    /** Number of polygons written to the vtk file */
    virtual int get_n_polygons() const { return 0; }

    /** Get the points for writing to vtk file
     * @param zmin minimum z coordinate*/
    virtual vector<Point3> get_points(const double) { return vector<Point3>{}; }

    /** Get nodes of polygons for writing to vtk file */
    virtual vector<vector<int>> get_polygons() { return vector<vector<int>>{vector<int>{}}; }

    /** Return active cut-off radius */
    double get_cutoff2(const Point3& p) {
        pick_cutoff(p);
        return cutoff2;
    }

protected:
    Point3 origin3d;  ///< centre of the coarsener
    double exponential; ///< exponential factor that determines coarsening rate outside the region of interest
    double cutoff2;   ///< squared cut off radius
    double radius;    ///< radius of the system
    double radius2;   ///< squared radius of the system
    double r0_min;    ///< cut off radius for the points on the edge of the system
    double r0_max;    ///< cut off radius for the points far away from of the system
    double A;         ///< coarsening factor

    const int n_nodes_per_circle = 50; ///< number of nodes in vtk circle-polygon
    const int n_nodes_per_line = 2;    ///< number of nodes in vtk line-polygon
    const int n_circles = 0;           ///< number of circles written to vtk file
    const int n_lines = 0;             ///< number of lines written to vtk file

    /** Get constant squared cut off radius */
    double get_const_cutoff() const { return r0_min * r0_min; }

    /** Get cut off radius that increases with distance from origin */
    double get_increasing_cutoff(const Point3 &point) const {
        double cutoff = max(0.0, origin3d.distance(point) - radius);
        cutoff = min( r0_max, A * pow(cutoff, exponential) + r0_min );
        return cutoff * cutoff;
    }

    /** Get cut off radius that is smaller than any possible distance between atoms */
    double get_inf_cutoff() const { return -1e100; }

    /** Point in region? */
    virtual inline bool in_region(const Point3 &) const { return false; }
};


/** Class to coarsen whole area uniformly */
class ConstCoarsener: public Coarsener {
public:
    ConstCoarsener();
    ConstCoarsener(const double r0_min);

private:
    /** Point in anywhere? */
    inline bool in_region(const Point3 &) const { return true; }
};

/** Class to coarsen surface outside one infinite vertical cylinder */
class FlatlandCoarsener: public Coarsener {
public:
    FlatlandCoarsener();
    FlatlandCoarsener(const Point3 &origin, const double exponential, const double radius,
            const double A, const double r0_min=0, const double r0_max=1e20);

    /** Points OUTSIDE the region will be coarsened */
    void pick_cutoff(const Point3 &point) {
        if (in_region(point))
            cutoff2 = get_increasing_cutoff(point);
        else
            cutoff2 = get_inf_cutoff();
    }

private:
    Point2 origin2d;

    /** Point outside cone with big apex angle? */
    inline bool in_region(const Point3 &point) const {
        constexpr double min_angle = M_PI / 6.0;
        Vec3 diff(point - origin3d);
        return asin(diff.z / diff.norm()) < min_angle;
    }
};

/** Class to coarsen surface inside one infinite vertical cylinder */
class CylinderCoarsener: public Coarsener {
public:
    CylinderCoarsener();
    CylinderCoarsener(const Point2 &base, const double exponential, const double radius,
            const double r0_cylinder=0);

    /** Points INSIDE the nanotip will be coarsened */
    void pick_cutoff(const Point3 &point) {
        if (in_region(point))
            cutoff2 = get_const_cutoff();
        else
            cutoff2 = get_inf_cutoff();
    }

    vector<Point3> get_points(const double zmin);
    vector<vector<int>> get_polygons();
    int get_n_points() const { return n_lines * n_nodes_per_line + n_circles * n_nodes_per_circle; }
    int get_n_polygons() const { return n_lines + n_circles; };

protected:
    /** Point in infinitely high vertical cylinder? */
    inline bool in_region(const Point3 &point) const {
        return point.distance2(origin2d) <= radius2;
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
    NanotipCoarsener(const Point3 &apex, const double exponential, const double radius,
            const double A, const double r0_apex, const double r0_cylinder);

    /** Points INSIDE the nanotip will be coarsened */
    void pick_cutoff(const Point3 &point) {
        if (in_region(point))
            cutoff2 = get_increasing_cutoff(point);
        else
            cutoff2 = get_inf_cutoff();
    }

    virtual vector<Point3> get_points(const double zmin);
    vector<vector<int>> get_polygons();
    int get_n_points() const { return n_lines * n_nodes_per_line + n_circles * n_nodes_per_circle; }
    int get_n_polygons() const { return n_lines + n_circles; };

protected:
    /** Point in infinitely high vertical cylinder? */
    virtual inline bool in_region(const Point3 &point) const {
        return point.distance2(bottom) <= radius2;
    }

    const int n_circles = 4;
    const int n_lines = 6;
    Point2 bottom;
};

/** Class to coarsen surface inside one vertical cone with a sphere on top */
class ConeCoarsener: public NanotipCoarsener {
public:
    ConeCoarsener();
    ConeCoarsener(const Point3 &apex, const Point3 &base, double theta,
            double exponential, double radius, double A, double r0_apex, double r0_cone);

    vector<Point3> get_points(const double zmin);

protected:
    /** Point in infinitely high vertical cone? */
    inline bool in_region(const Point3 &point) const {
        double r = max(0.0, r_bottom + (point.z - z_bottom) * tan_theta);
        return point.distance2(bottom) <= r*r;
    }

private:
    double tan_theta;  ///< tangent of cone apex angle
    double z_bottom;   ///< z-coordinate of nanotip-substrate junction
    double r_bottom;
};

/** Class to coarsen surface inside one infinite tilted nanotip */
class TiltedNanotipCoarsener: public NanotipCoarsener {
public:
    TiltedNanotipCoarsener();
    TiltedNanotipCoarsener(const Point3 &apex, const Point3 &base, const double exponential,
            const double radius, const double A, const double r0_apex, const double r0_cylinder);

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
class Coarseners: public FileWriter {
public:
    /** Coarseners constructor */
    Coarseners() : tan_theta(0), radius(0), amplitude(0), r0_cylinder(0) {}

    /** Append coarsener to the array of all the coarseners */
    void attach_coarsener(shared_ptr<Coarsener> c) {
        coarseners.push_back(c);
    }

    /** Pick cut off radiuses for coarseners */
    void pick_cutoff(const Point3 &point) {
        for (auto &c : coarseners)
            c->pick_cutoff(point);
    }

    /** Calculate the cut off radius for given point */
    double get_cutoff(const Point3 &point) {
        for (auto &c : coarseners) {
            double cutoff = c->get_cutoff2(point);
            if (cutoff >= 0) return sqrt(cutoff);
        }
        return -1;
    }

    /** Specify the collective action of coarseners */
    bool nearby(const Point3 &p1, const Point3 &p2) {
        bool near = false;
        for (auto &c : coarseners)
            near |= c->nearby(p1, p2);

        return near;
    }

    /** Generate coarseners for one nanotip system */
    void generate(const Medium &medium, const Config::Geometry &conf, const Config::CoarseFactor &cf);

    /** Get the distance between atoms on the edge of simulation cell */
    double get_r0_inf(const Medium::Sizes &s);

    /** Radius of coarsening cylinder */
    double get_radius() const { return radius; }

    /** Return whether point is inside region-of-interest (nanotip or similar) */
    bool inside_roi(const Point3& p) const;

    Point3 centre;
private:
    double tan_theta;
    double radius;
    double amplitude;
    double r0_cylinder;
    vector<shared_ptr<Coarsener>> coarseners;

    /** Get histogram for atom z-coordinates */
    void get_histogram(vector<int> &bins, vector<double> &bounds, const Medium& medium);

    /** Get the average z-coordinate of substrate atoms */
    double get_z_mean(const Medium& medium);

    /** Write the contours of coarseners to file in .vtk format */
    void write_vtk(ofstream &out) const;

    /** Specify file types that can be written */
    bool valid_extension(const string &ext) const {
        return ext == "vtk" ||  ext == "vtks";
    }
};

} // namespace femocs

#endif /* COARSENER_H_ */
