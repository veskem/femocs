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
#include <math.h>

using namespace std;
namespace femocs {

Coarsener::Coarsener() :
        origin3d(Point3(0)), exponential(0.5), cutoff2(0), radius(0), radius2(0), r0_min(0), r0_max(0), A(0) {}

Coarsener::Coarsener(const Point3 &origin, const double exp, const double radius,
        const double A, const double r0_min, const double r0_max) :
        origin3d(origin), exponential(exp), cutoff2(0), radius(radius), radius2(radius*radius),
        r0_min(r0_min), r0_max(r0_max), A(A) {}

ConstCoarsener::ConstCoarsener() : Coarsener() {}

ConstCoarsener::ConstCoarsener(const double r0_min) : Coarsener(Point3(), 0, 0, 0, r0_min) {}

FlatlandCoarsener::FlatlandCoarsener() : Coarsener(), origin2d(0.0) {}

FlatlandCoarsener::FlatlandCoarsener(const Point3 &origin, const double exp, const double radius,
        const double A, const double r0_min, const double r0_max) :
        Coarsener(origin, exp, radius, A, r0_min, r0_max),
        origin2d(origin.x, origin.y) {}

CylinderCoarsener::CylinderCoarsener() : Coarsener(), origin2d(0.0) {}

CylinderCoarsener::CylinderCoarsener(const Point2 &base, const double exp, const double radius,
        const double r0_cylinder) :
        Coarsener(Point3(), exp, radius, 0, r0_cylinder), origin2d(base.x, base.y) {}

NanotipCoarsener::NanotipCoarsener() : Coarsener(), origin2d(0.0) {}

NanotipCoarsener::NanotipCoarsener(const Point3 &apex, const double exp, const double radius,
        const double A, const double r0_apex, const double r0_cylinder) :
        Coarsener(apex, exp, radius, A, r0_apex, r0_cylinder), origin2d(apex.x, apex.y) {}

TiltedNanotipCoarsener::TiltedNanotipCoarsener() : NanotipCoarsener(), bottom(Vec3()), axis(Vec3()), height2(0) {}

TiltedNanotipCoarsener::TiltedNanotipCoarsener(const Point3 &apex, const Point3 &base, const double exp,
        const double radius, const double A, const double r0_apex, const double r0_cylinder) :
        NanotipCoarsener(apex, exp, radius, A, r0_apex, r0_cylinder), bottom(Vec3(base.x, base.y, base.z)) {

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

vector<Point3> CylinderCoarsener::get_points(const double zmin) {
    const double point_res = 2 * M_PI / n_nodes_per_circle;
    const double line_res = 2 * M_PI / n_lines;

    // Reserve memory for points
    vector<Point3> points;
    points.reserve(get_n_points());

    // Make points for circle
    for (double a = 0; a < 2*M_PI; a += point_res)
        points.push_back( Point3(radius*cos(a), radius*sin(a), zmin) );

    // Make points for lines
    for (double a = 0; a < 2*M_PI; a += line_res) {
        points.push_back( Point3(radius*cos(a), radius*sin(a), 1.1*zmin) );
        points.push_back( Point3(radius*cos(a), radius*sin(a), 0.9*zmin) );
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

vector<Point3> NanotipCoarsener::get_points(const double zmin) {
    const double circle_res = 2 * M_PI / n_nodes_per_circle;
    const double line_res = 2 * M_PI / 4;
    const double zbottom = zmin - origin3d.z;

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
        points.push_back( Point3(radius*cos(a), radius*sin(a), zbottom) );
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

void Coarseners::generate(const Medium &medium, const double radius,
    const Config::CoarseFactor &cf, const double latconst)
{
    const int n_atoms = medium.size();
    require(n_atoms > 0, "Not enough points to generate coarseners.");
    require(cf.r0_cylinder >= 0 && cf.r0_sphere >= 0, "Coarsening factors must be non-negative!");
    require(cf.r0_cylinder >= cf.r0_sphere, "Coarsening factor in cylinder wall must be >= coarsening factor in apex!");
    require(cf.exponential > 0, "Coarsening rate must be positive!");

    const double z_bot = get_z_mean(medium);
    const double z_top = max(z_bot, medium.sizes.zmax - 0.5*radius);

    centre = Point3(medium.sizes.xmid, medium.sizes.ymid, z_bot);
    Point3 apex(medium.sizes.xmid, medium.sizes.ymid, z_top);

    this->radius = radius;                 // store the coarsener radius
    amplitude = cf.amplitude * latconst;
    r0_cylinder = cf.r0_cylinder * 0.25 * latconst;
    const double r0_sphere = cf.r0_sphere * 0.25 * latconst;
    const double r0_flat = min(amplitude*1e20, r0_cylinder);

    coarseners.clear();
    attach_coarsener( make_shared<NanotipCoarsener>(apex, cf.exponential, radius, amplitude, r0_sphere, r0_cylinder) );
    attach_coarsener( make_shared<FlatlandCoarsener>(centre, cf.exponential, radius, amplitude, r0_flat) );
}

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
        if (z >= zmin && z <= zmax) {
            zmean += z;
            n_average++;
        }
    }

    require(n_average > 0, "Error calculating the average z of a surface atoms!");
    return zmean / n_average;
}

void Coarseners::get_histogram(vector<int> &bins, vector<double> &bounds, const Medium& medium) {
    const int n_atoms = medium.size();
    const int n_bins = bins.size();
    const int n_bounds = bounds.size();

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
    bounds[n_bounds-1] += 1e-5 * value_step;

    for (int i = 0; i < n_atoms; ++i) {
        value = medium.get_point(i).z;
        for (int j = 0; j < n_bins; ++j)
            if (value >= bounds[j] && value < bounds[j+1]) {
                bins[j]++;
                break;
            }
    }
}

double Coarseners::get_r0_inf(const Medium::Sizes &s) {
    const double max_distance = centre.distance(Point3(s.xmin, s.ymin, s.zmin));
    if ((max_distance - radius) > 0)
        return 1.1 * amplitude * sqrt(max_distance - radius) + r0_cylinder;
    else
        return r0_cylinder;
}

bool Coarseners::inside_interesting_region(const Point3& p) const {
    const double min_angle = M_PI / 6.0;
    Vec3 diff(p - centre);
    bool inside_tip = diff.x * diff.x + diff.y * diff.y <= radius * radius;
    bool inside_cone = asin(diff.z / diff.norm()) >= min_angle;
    return inside_tip && inside_cone;
}

bool Coarseners::valid_extension(const string &ext) const {
    return ext == "vtk";
}

void Coarseners::write_vtk(ofstream &out) const {
    out << "# vtk DataFile Version 3.0\n";
    out << "# " + class_name() + "\n";
    out << "ASCII\n";
    out << "DATASET POLYDATA\n";

    // Calculate total number of points
    size_t n_points = 0;
    for (auto &c : coarseners)
        n_points += c->get_n_points();

    out << "POINTS " << n_points << " float\n";

    // Write points to file
    for (auto &c : coarseners)
        for (Point3 p : c->get_points(centre.z))
            out << p << "\n";

    // Calculate total number of polygons and nodes
    size_t n_polygons = 0;
    size_t n_total = 0;
    for (auto &c : coarseners) {
        n_polygons += c->get_n_polygons();
        n_total += c->get_n_polygons() + c->get_n_points();
    }

    out << "POLYGONS " << n_polygons << " " << n_total << "\n";

    // Write polygons to file
    int offset = 0;
    for (auto &c : coarseners) {
        for (vector<int> polygon : c->get_polygons()) {
            int n_nodes = polygon.size(); // number of nodes in polygon
            if (n_nodes > 0) {
                out << n_nodes;
                for (int node : polygon)
                    out << " " << (node + offset);   // index of node
                out << "\n";
            }
        }
        offset += c->get_n_points();
    }
}

} // namespace femocs
