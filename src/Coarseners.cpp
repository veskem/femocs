/*
 * Coarsener.cpp
 *
 *  Created on: 7.10.2016
 *      Author: veske
 */

#include "Coarseners.h"
#include <fstream>

using namespace std;
namespace femocs {

// Coarsener constructors
Coarsener::Coarsener() :
        origin3d(Point3(0)), cutoff2(0), radius(0), radius2(0), A(0), r0_min(0), r0_max(0) {}

Coarsener::Coarsener(const Point3 &origin, double radius, double A, double r0_min, double r0_max) :
        origin3d(origin), cutoff2(0), radius(radius), radius2(radius*radius), A(A), r0_min(r0_min), r0_max(r0_max) {}

// ConstCoarsener for coarsening area uniformly
ConstCoarsener::ConstCoarsener() : Coarsener() {}

ConstCoarsener::ConstCoarsener(double r0_min) : Coarsener(Point3(), 0, 0, r0_min) {}

// FlatlandCoarsener for coarsening area outside one infinite cylinder
FlatlandCoarsener::FlatlandCoarsener() : Coarsener(), origin2d(Point2(0.0)) {}

FlatlandCoarsener::FlatlandCoarsener(const Point3 &origin, double radius, double A, double r0_min, double r0_max) :
        Coarsener(origin, radius, A, r0_min, r0_max),
        origin2d(Point2(origin.x, origin.y)) {}

// CylinderCoarsener for coarsening area inside one infinite cylinder
CylinderCoarsener::CylinderCoarsener() : Coarsener() {}

CylinderCoarsener::CylinderCoarsener(const Point2 &base, double radius, double r0_cylinder) :
        Coarsener(Point3(), radius, 0, r0_cylinder), origin2d(base) {}

// NanotipCoarsener for coarsening area inside one infinite vertical nanotip
NanotipCoarsener::NanotipCoarsener() : Coarsener(), origin2d(Point2()) {}

NanotipCoarsener::NanotipCoarsener(const Point3 &apex, double radius, double A, double r0_apex, double r0_cylinder) :
        Coarsener(apex, radius, A, r0_apex, r0_cylinder), origin2d(Point2(apex.x, apex.y)) {}

// TiltedNanotipCoarsener for coarsening area inside one infinite tilted nanotip
TiltedNanotipCoarsener::TiltedNanotipCoarsener() : NanotipCoarsener(), bottom(Vec3()), axis(Vec3()), height2(0) {}

TiltedNanotipCoarsener::TiltedNanotipCoarsener(const Point3 &apex, const Point3 &base, double radius,double A, double r0_apex, double r0_cylinder) :
        NanotipCoarsener(apex, radius, A, r0_apex, r0_cylinder), bottom(Vec3(base.x, base.y, base.z)) {

    Vec3 top(apex.x, apex.y, apex.z);
    axis = top - bottom;
    height2 = axis.length2();
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
    for (int i = 0; i < points.size(); ++i)
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
    const double zbottom = origin3d.z - 2*radius;

    // Reserve memory for points
    vector<Point3> points;
    points.reserve(get_n_points());

    // Make points for apex circle
    for (double a = 0; a < 2*M_PI; a += circle_res)
        points.push_back(Point3(radius*cos(a), radius*sin(a), origin3d.z));

    // Make points for bottom circle
    for (double a = 0; a < 2*M_PI; a += circle_res)
        points.push_back(Point3(radius*cos(a), radius*sin(a), zbottom));

    // Make points for circle in y-z plane
    for (double a = 0; a < 2*M_PI; a += circle_res)
        points.push_back(Point3(0, radius*cos(a), origin3d.z+radius*sin(a)));

    // Make points for circle in x-z plane
    for (double a = 0; a < 2*M_PI; a += circle_res)
        points.push_back(Point3(radius*cos(a), 0, origin3d.z+radius*sin(a)));

    // Make points for vertical lines
    for (double a = 0; a < 2*M_PI; a += line_res) {
        points.push_back( Point3(radius*cos(a), radius*sin(a), origin3d.z) );
        points.push_back( Point3(radius*cos(a), radius*sin(a), zbottom-1) );
    }

    // Make points for horizontal lines
    points.push_back( Point3(radius, 0, origin3d.z) );
    points.push_back( Point3(-radius,0, origin3d.z) );
    points.push_back( Point3(0, radius, origin3d.z) );
    points.push_back( Point3(0,-radius, origin3d.z) );

    // Shift points to the origin
    for (int i = 0; i < points.size(); ++i)
        points[i] += origin2d;

    return points;
}

void Coarseners::write(const string &file_name) {
#if not DEBUGMODE
    return;
#endif

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
