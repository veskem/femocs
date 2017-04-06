/*
 * Coarsener.cpp
 *
 *  Created on: 7.10.2016
 *      Author: veske
 */

#include "Coarseners.h"
#include <fstream>
#include <numeric>
#include <algorithm>
#include <float.h>

using namespace std;
namespace femocs {

// Coarsener constructors
Coarsener::Coarsener() :
        origin3d(Point3(0)), cutoff2(0), radius(0), radius2(0), r0_min(0), r0_max(0), A(0) {}

Coarsener::Coarsener(const Point3 &origin, const double radius, const double A, const double r0_min, const double r0_max) :
        origin3d(origin), cutoff2(0), radius(radius), radius2(radius*radius), r0_min(r0_min), r0_max(r0_max), A(A) {}

// ConstCoarsener for coarsening area uniformly
ConstCoarsener::ConstCoarsener() : Coarsener() {}

ConstCoarsener::ConstCoarsener(const double r0_min) : Coarsener(Point3(), 0, 0, r0_min) {}

// FlatlandCoarsener for coarsening area outside one infinite cylinder
FlatlandCoarsener::FlatlandCoarsener() : Coarsener(), origin2d(Point2(0.0)) {}

FlatlandCoarsener::FlatlandCoarsener(const Point3 &origin, const double radius, const double A,
        const double r0_min, const double r0_max) :
        Coarsener(origin, radius, A, r0_min, r0_max),
        origin2d(Point2(origin.x, origin.y)) {}

// CylinderCoarsener for coarsening area inside one infinite cylinder
CylinderCoarsener::CylinderCoarsener() : Coarsener() {}

CylinderCoarsener::CylinderCoarsener(const Point2 &base, const double radius, const double r0_cylinder) :
        Coarsener(Point3(), radius, 0, r0_cylinder), origin2d(base) {}

// NanotipCoarsener for coarsening area inside one infinite vertical nanotip
NanotipCoarsener::NanotipCoarsener() : Coarsener(), origin2d(Point2()) {}

NanotipCoarsener::NanotipCoarsener(const Point3 &apex, const double radius, const double A,
        const double r0_apex, const double r0_cylinder) :
        Coarsener(apex, radius, A, r0_apex, r0_cylinder), origin2d(Point2(apex.x, apex.y)) {}

// TiltedNanotipCoarsener for coarsening area inside one infinite tilted nanotip
TiltedNanotipCoarsener::TiltedNanotipCoarsener() : NanotipCoarsener(), bottom(Vec3()), axis(Vec3()), height2(0) {}

TiltedNanotipCoarsener::TiltedNanotipCoarsener(const Point3 &apex, const Point3 &base, const double radius,
        const double A, const double r0_apex, const double r0_cylinder) :
        NanotipCoarsener(apex, radius, A, r0_apex, r0_cylinder), bottom(Vec3(base.x, base.y, base.z)) {

    Vec3 top(apex.x, apex.y, apex.z);
    axis = top - bottom;
    height2 = axis.norm2();
}

vector<vector<int>> CylinderCoarsener::get_polygons() {
    // Reserve memory for nodes
    vector<vector<int>> polys;
    polys.resize(get_n_polygons());

    // Make nodes for circle
    polys[0].resize(n_nodes_per_circle);
    std::iota (begin(polys[0]), end(polys[0]), 0);

    // Make nodes for lines
    for (int i = 0; i < n_lines; ++i) {
        int node = 2 * i + n_circles * n_nodes_per_circle;
        polys[n_circles+i] = vector<int>{node, node+1};
    }

    return polys;
}

vector<Point3> CylinderCoarsener::get_points() {
    const double point_res = 2 * M_PI / n_nodes_per_circle;
    const double line_res = 2 * M_PI / n_lines;

    // Reserve memory for points
    vector<Point3> points;
    points.reserve(get_n_points());

    // Make points for circle
    for (double a = 0; a < 2*M_PI; a += point_res)
        points.push_back( Point3(radius*cos(a), radius*sin(a), 0) );

    // Make points for lines
    for (double a = 0; a < 2*M_PI; a += line_res) {
        points.push_back( Point3(radius*cos(a), radius*sin(a), 1) );
        points.push_back( Point3(radius*cos(a), radius*sin(a),-1) );
    }

    // Shift points to the origin in x-y plane
    for (unsigned i = 0; i < points.size(); ++i)
        points[i] += origin2d;

    return points;
}

vector<vector<int>> NanotipCoarsener::get_polygons() {
    // Reserve memory for nodes
    vector<vector<int>> polys;
    polys.resize(get_n_polygons());

    // Make nodes for apex circle
    polys[0].resize(n_nodes_per_circle);
    std::iota (begin(polys[0]), end(polys[0]), 0);

    // Make nodes for bottom circle
    polys[1].resize(n_nodes_per_circle);
    std::iota (begin(polys[1]), end(polys[1]), 1*n_nodes_per_circle);

    // Make points for circle in y-z plane
    polys[2].resize(n_nodes_per_circle);
    std::iota (begin(polys[2]), end(polys[2]), 2*n_nodes_per_circle);

    // Make points for circle in x-z plane
    polys[3].resize(n_nodes_per_circle);
    std::iota (begin(polys[3]), end(polys[3]), 3*n_nodes_per_circle);

    // Make nodes for lines
    for (int i = 0; i < n_lines; ++i) {
        int node = 2 * i + n_circles*n_nodes_per_circle;
        polys[n_circles+i] = vector<int>{node, node+1};
    }

    return polys;
}

vector<Point3> NanotipCoarsener::get_points() {
    const double circle_res = 2 * M_PI / n_nodes_per_circle;
    const double line_res = 2 * M_PI / 4;
    const double zbottom = -2.0 * radius;

    // Reserve memory for points
    vector<Point3> points;
    points.reserve(get_n_points());

    // Make points for apex circle
    for (double a = 0; a < 2*M_PI; a += circle_res)
        points.push_back(Point3(radius*cos(a), radius*sin(a), 0));

    // Make points for bottom circle
    for (double a = 0; a < 2*M_PI; a += circle_res)
        points.push_back(Point3(radius*cos(a), radius*sin(a), zbottom));

    // Make points for circle in y-z plane
    for (double a = 0; a < 2*M_PI; a += circle_res)
        points.push_back(Point3(0, radius*cos(a), radius*sin(a)));

    // Make points for circle in x-z plane
    for (double a = 0; a < 2*M_PI; a += circle_res)
        points.push_back(Point3(radius*cos(a), 0, radius*sin(a)));

    // Make points for vertical lines
    for (double a = 0; a < 2*M_PI; a += line_res) {
        points.push_back( Point3(radius*cos(a), radius*sin(a), 0) );
        points.push_back( Point3(radius*cos(a), radius*sin(a), zbottom-1) );
    }

    // Make points for horizontal lines
    points.push_back( Point3(radius, 0, 0) );
    points.push_back( Point3(-radius,0, 0) );
    points.push_back( Point3(0, radius, 0) );
    points.push_back( Point3(0,-radius, 0) );

    // Shift points according to the origin
    for (unsigned i = 0; i < points.size(); ++i)
        points[i] += origin3d;

    return points;
}

// Generate coarseners for one nanotip system
void Coarseners::generate(const Medium &medium, const double radius, const Config::CoarseFactor &cf, const double latconst) {
    const int n_atoms = medium.size();
    require(n_atoms > 0, "Not enough points to generate coarseners: " + to_string(n_atoms));

    const double z_bot = get_z_mean(medium);  // medium.sizes.zmin;
    const double z_top = max(z_bot, medium.sizes.zmax - 0.5*radius);

    centre = Point3(medium.sizes.xmid, medium.sizes.ymid, z_bot);
    Point3 apex(medium.sizes.xmid, medium.sizes.ymid, z_top);

    this->radius = radius;                 // store the coarsener radius
    amplitude = cf.amplitude * latconst;
    r0_cylinder = max(0.0, cf.r0_cylinder) * latconst;
    const double r0_sphere = max(0.0, cf.r0_sphere) * latconst;
    const double r0_flat = min(amplitude*1e20, r0_cylinder);

    coarseners.clear();
    attach_coarsener( make_shared<NanotipCoarsener>(apex, radius, amplitude, r0_sphere, r0_cylinder) );
    attach_coarsener( make_shared<FlatlandCoarsener>(centre, radius, amplitude, r0_flat) );
}

// Return the average z-coordinate of substrate atoms
double Coarseners::get_z_mean(const Medium& medium) {
    const int n_atoms = medium.size();
    const int n_bins = 100;

    // make the histogram for all the surface atoms
    vector<int> bins(n_bins, 0);
    vector<double> bounds(n_bins+1);
    get_histogram(bins, bounds, medium);

    // locate the bin with maximum amount of entries
    int max_bin = distance( bins.begin(), max_element(bins.begin(), bins.end()) );
    double zmin = bounds[max_bin];
    double zmax = bounds[max_bin+1];

    // For systems with small amount of atoms in the substrate, use the minimum coordinate
    if ( (zmin + zmax) / 2 > medium.sizes.zmean )
        return medium.sizes.zmin;

    // average z-coordinate == average for the entries of most popular bin
    double zmean = 0;
    int n_average = 0;
    for (int i = 0; i< n_atoms; ++i) {
        double z = medium.get_point(i).z;
        if (z >= zmin && z < zmax) {
            zmean += z;
            n_average++;
        }
    }

    return zmean / n_average;
}

// Get histogram for atom z-coordinates
void Coarseners::get_histogram(vector<int> &bins, vector<double> &bounds, const Medium& medium) {
    const int n_atoms = medium.size();
    const int n_bins = bins.size();
    const int n_bounds = bounds.size();
    const double eps = 1e-5;

    // Find minimum and maximum values from all z-coordinates
    double value_min = DBL_MAX;
    double value_max =-DBL_MAX;
    double value;
    for (int i = 0; i < n_atoms; ++i) {
        value = medium.get_point(i).z;
        value_min = min(value_min, value);
        value_max = max(value_max, value);
    }

    // Fill the bounds with values value_min:value_step:(value_max + epsilon)
    // Epsilon is added to value_max to include the maximum value in the up-most bin
    double value_step = (value_max - value_min) / n_bins;
    for (int i = 0; i < n_bounds; ++i)
        bounds[i] = value_min + value_step * i;
    bounds[n_bounds-1] += eps;

    for (int i = 0; i < n_atoms; ++i) {
        value = medium.get_point(i).z;
        for (int j = 0; j < n_bins; ++j)
            if (value >= bounds[j] && value < bounds[j+1]) {
                bins[j]++;
                break;
            }
    }
}

// Get the distance between atoms on the edge of simulation cell
double Coarseners::get_r0_inf(const Medium::Sizes &s) {
    const double max_distance = centre.distance(Point3(s.xmin, s.ymin, s.zmin));
    if ((max_distance - radius) > 0)
        return 1.1 * amplitude * sqrt(max_distance - radius) + r0_cylinder;
    else
        return r0_cylinder;
}

// Write the contours of coarseners to file in .vtk format
void Coarseners::write(const string &file_name) {
    if (!MODES.WRITEFILE) return;

    string file_type = get_file_type(file_name);
    require(file_type == "vtk", "Unimplemented file type: " + file_type);

    std::ofstream out(file_name.c_str());
    require(out, "File " + file_name + " cannot be opened for writing!");

    out.precision(8);

    out << "# vtk DataFile Version 3.0\n";
    out << "# Coarsener data\n";
    out << "ASCII\n";
    out << "DATASET POLYDATA\n";

    // Calculate total number of points
    size_t n_points = 0;
    for (auto &c : coarseners)
        n_points += c->get_n_points();

    out << "\nPOINTS " << n_points << " float\n";

    // Write points to file
    for (auto &c : coarseners)
        for (Point3 p : c->get_points())
            out << p << "\n";

    // Calculate total number of polygons and nodes
    size_t n_polygons = 0;
    size_t n_total = 0;
    for (auto &c : coarseners) {
        n_polygons += c->get_n_polygons();
        n_total += c->get_n_polygons() + c->get_n_points();
    }

    out << "\nPOLYGONS " << n_polygons << " " << n_total << "\n";

    // Write polygons to file
    int offset = 0;
    for (auto &c : coarseners) {
        for (vector<int> polygon : c->get_polygons()) {
            out << polygon.size() << " ";  // number of nodes
            for (int node : polygon)
                out << (node + offset) << " ";        // index of node
            out << "\n";
        }
        offset += c->get_n_points();
    }
}

} // namespace femocs
