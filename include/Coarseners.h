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
    Coarsener(const Point3 &origin, double exp, double radius,
            double A, double r0_min, double r0_max);

    virtual ~Coarsener() {};

    /** Choose clearance distance for atoms inside the region */
    virtual void pick_cutoff(const Point3 &point);

    /** Determine if two points are within their cut-off radius */
    bool nearby(const Point3 &p1, const Point3 &p2) const;

    /** Return active cut-off radius */
    double get_cutoff2(const Point3& p) const { return cutoff2; }

    /** Number of points written to the vtk file */
    virtual int get_n_points() const { return 0; }

    /** Number of polygons written to the vtk file */
    virtual int get_n_polygons() const { return 0; }

    /** Get the points for writing to vtk file
     * @param zmin minimum z coordinate*/
    virtual vector<Point3> get_points() const {
        return vector<Point3>{};
    }

    /** Get nodes of polygons for writing to vtk file */
    virtual vector<vector<int>> get_polygons() const {
        return vector<vector<int>>{vector<int>{}};
    }

protected:
    static constexpr int n_nodes_per_circle = 50; ///< number of nodes in vtk circle-polygon
    static constexpr int n_nodes_per_line = 2;    ///< number of nodes in vtk line-polygon

    Point3 origin3d;  ///< centre of the coarsener
    double exponential; ///< exponential factor that determines coarsening rate outside the region of interest
    double cutoff2;   ///< squared cut off radius
    double radius;    ///< radius of the system
    double radius2;   ///< squared radius of the system
    double r0_min;    ///< cut off radius for the points on the edge of the system
    double r0_max;    ///< cut off radius for the points far away from of the system
    double A;         ///< coarsening factor

    /** Point in region? */
    virtual bool in_region(const Point3 &) const { return false; }

    /** Get constant squared cut off radius */
    double get_const_cutoff() const;

    /** Get cut off radius that increases with distance from origin */
    double get_increasing_cutoff(const Point3 &point) const;

    /** Get cut off radius that is smaller than any possible distance between atoms */
    double get_inf_cutoff() const { return -1e100; }
};


/** Class to coarsen whole area uniformly */
class ConstCoarsener: public Coarsener {
public:
    ConstCoarsener();
    ConstCoarsener(double r0_min);

private:
    /** Point in anywhere? */
    bool in_region(const Point3 &) const { return true; }
};

/** Class to coarsen the surface atoms outside the nanotip */
class FlatlandCoarsener: public Coarsener {
public:
    FlatlandCoarsener();
    FlatlandCoarsener(const Point3 &origin, double exponential, double radius,
            double A, double r0_min=0, double r0_max=1e20);

    /** Choose clearance distance for atoms on the flat surface */
    void pick_cutoff(const Point3 &point);

private:
    /** Point outside a cone with big apex angle? */
    bool in_region(const Point3 &point) const;
};

/** Class to coarsen the surface atoms inside one vertical nanotip
 * with conical walls and sphere on top */
class NanotipCoarsener: public Coarsener {
public:
    NanotipCoarsener();
    NanotipCoarsener(const Point3 &apex, const Point3 &base, double theta,
            double exponential, double radius, double A, double r0_apex, double r0_cone);

    /** Choose clearance distance for atoms inside the nanotip */
    void pick_cutoff(const Point3 &point);

    /** Generate points for vtk file */
    vector<Point3> get_points() const;

    /** Generate polygons for vtk file */
    vector<vector<int>> get_polygons() const;

    /** Nr of points added to vtk file */
    int get_n_points() const {
        return n_lines * n_nodes_per_line + n_circles * n_nodes_per_circle;
    }

    /** Nr of polygons added to vtk file */
    int get_n_polygons() const {
        return n_lines + n_circles;
    }

protected:
    static constexpr int n_circles = 4; ///< nr on circles in vtk file
    static constexpr int n_lines = 6;   ///< nr of lines in vtk file
    Point2 bottom;     ///< centre of a cone
    double tan_theta;  ///< tangent of cone apex angle
    double z_bottom;   ///< z-coordinate of nanotip-substrate junction
    double r_bottom;   ///< radius of cone on nanotip-substrate junction

    /** Point inside the nanotip? */
    inline bool in_region(const Point3 &point) const;
};

/** Class for binding together coarseners */
class Coarseners: public FileWriter {
public:
    /** Coarseners constructor */
    Coarseners() : tan_theta(0), radius(0), amplitude(0), r0_wall(0) {}

    /** Generate coarseners for one nanotip system */
    void generate(const Medium &medium, const Config::Geometry &conf, const Config::CoarseFactor &cf);

    /** Pick cut-off radiuses for coarseners */
    void pick_cutoff(const Point3 &point);

    /** Specify the collective action of coarseners */
    bool nearby(const Point3 &p1, const Point3 &p2) const;

    /** Return whether point is inside region-of-interest (nanotip or similar) */
    bool inside_roi(const Point3& p) const;

    /** Calculate the cut-off radius for given point */
    double get_cutoff(const Point3 &point);

    /** Get the distance between atoms on the edge of simulation cell */
    double get_r0_inf(const Medium::Sizes &s) const;

    /** Radius of coarsening cylinder */
    double get_radius() const { return radius; }

    /** Average z-coordinate of flat region */
    double get_z_mean() const { return centre.z; }

private:
    Point3 centre;   ///< coordinates of nanotip-substrate junction
    double tan_theta; ///< tangent of apex angle of coarsening cone
    double radius;    ///< radius of nanotip-substrate junction
    double amplitude; ///< coarsening amplitude aka factor
    double r0_wall;   ///< min distance between atoms on nanotip wall
    vector<shared_ptr<Coarsener>> coarseners;  ///< coarseners of various regions

    /** Calculate histogram for atom z-coordinates */
    void calc_histogram(vector<int> &bins, vector<double> &bounds, const Medium& medium) const;

    /** Get the average z-coordinate of substrate atoms */
    double calc_z_mean(const Medium& medium) const;

    /** Write the contours of coarseners to file in .vtk format */
    void write_vtk(ofstream &out) const;

    /** Specify file types that can be written */
    bool valid_extension(const string &ext) const {
        return ext == "vtk" ||  ext == "vtks";
    }
};

} // namespace femocs

#endif /* COARSENER_H_ */
