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

void CylinderCoarsener::get_write_data(vector<Point3> &points, vector<vector<int>> &nodes) const {
    const int n_nodes_per_circle = 50;
    const int n_nodes_per_line = 2;
    const int n_circles = 1;
    const int n_lines = 4;
    const double point_res = 2 * M_PI / n_nodes_per_circle;
    const double line_res = 2 * M_PI / n_lines;

    // Reserve memory for points
    points.reserve(n_circles * n_nodes_per_circle + n_lines * n_nodes_per_line);

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

    // Reserve memory for nodes
    nodes.resize(n_circles + n_lines);

    // Make nodes for circle
    nodes[0].resize(n_nodes_per_circle);
    std::iota (begin(nodes[0]), end(nodes[0]), 0);

    // Make nodes for lines
    for (int i = 0; i < n_lines; ++i) {
        int node = 2 * i + n_circles * n_nodes_per_circle;
        nodes[n_circles+i] = vector<int>{node, node+1};
    }
}

void NanotipCoarsener::get_write_data(vector<Point3> &points, vector<vector<int>> &nodes) const {
    const int n_nodes_per_circle = 50;
    const int n_nodes_per_line = 2;
    const int n_nodes_per_arc = 50;
    const int n_circles = 2;
    const int n_lines = 4;
    const int n_arcs = 2;
    const double point_res = 2 * M_PI / n_nodes_per_circle;
    const double line_res = 2 * M_PI / n_lines;
    const double arc_res = M_PI / (n_nodes_per_arc-1);

    // Reserve memory for points
    points.reserve(n_circles*n_nodes_per_circle + n_lines*n_nodes_per_line + n_arcs*n_nodes_per_arc);

    // Make points for apex circle
    for (double a = 0; a < 2*M_PI; a += point_res)
        points.push_back(Point3(radius*cos(a), radius*sin(a), origin3d.z));

    // Make points for bottom circle
    for (double a = 0; a < 2*M_PI; a += point_res)
        points.push_back(Point3(radius*cos(a), radius*sin(a), 0));

    // Make points for arc in y-z plane
    for (double a = 0; a < (M_PI+1e-5); a += arc_res)
        points.push_back(Point3(0, radius*cos(a), origin3d.z+radius*sin(a)));

    // Make points for arc in x-z plane
    for (double a = 0; a < (M_PI+1e-5); a += arc_res)
        points.push_back(Point3(radius*cos(a), 0, origin3d.z+radius*sin(a)));

    // Make points for lines
    for (double a = 0; a < 2*M_PI; a += line_res) {
        points.push_back( Point3(radius*cos(a), radius*sin(a), origin3d.z) );
        points.push_back( Point3(radius*cos(a), radius*sin(a), -1) );
    }

    // Shift points to the origin
    for (int i = 0; i < points.size(); ++i)
        points[i] += origin2d;

    // Reserve memory for nodes
    nodes.resize(n_circles + n_lines + n_arcs);

    // Make nodes for apex circle
    nodes[0].resize(n_nodes_per_circle);
    std::iota (begin(nodes[0]), end(nodes[0]), 0);

    // Make nodes for bottom circle
    nodes[1].resize(n_nodes_per_circle);
    std::iota (begin(nodes[1]), end(nodes[1]), n_nodes_per_circle);

    // Make points for arc in y-z plane
    nodes[2].resize(n_nodes_per_arc);
    std::iota (begin(nodes[2]), end(nodes[2]), n_circles * n_nodes_per_circle);

    // Make points for arc in x-z plane
    nodes[3].resize(n_nodes_per_arc);
    std::iota (begin(nodes[3]), end(nodes[3]), n_circles * n_nodes_per_circle + n_nodes_per_arc);

    // Make nodes for lines
    for (int i = 0; i < n_lines; ++i) {
        int node = 2 * i + n_circles*n_nodes_per_circle + n_arcs*n_nodes_per_arc;
        nodes[n_circles+n_arcs+i] = vector<int>{node, node+1};
    }
}

void Coarsener::write(const string &file_name) const {
#if not DEBUGMODE
    return;
#endif

    string file_type = get_file_type(file_name);
    require(file_type == "vtk", "Unimplemented file type: " + file_type);

    // Generate data to be written to file
    vector<Point3> points;
    vector<vector<int>> nodes;
    get_write_data(points, nodes);

    size_t n_points = points.size();
    size_t n_polygons = nodes.size();
    size_t n_total = n_polygons;
    for (int i = 0; i < n_polygons; ++i)
        n_total += nodes[i].size();

    std::ofstream out(file_name.c_str());
    require(out, "File " + file_name + " cannot be opened for writing!");

//    out.setf(std::ios::scientific);
    out.precision(8);

    out << "# vtk DataFile Version 3.0\n";
    out << "# Coarsener data\n";
    out << "ASCII\n";
    out << "DATASET POLYDATA\n";

    out << "\nPOINTS " << n_points << " float\n";
    for (Point3 p : points)
        out << p << "\n";

    out << "\nPOLYGONS " << n_polygons << " " << n_total << "\n";
    for (size_t i = 0; i < n_polygons; ++i) {
        size_t n_nodes = nodes[i].size();
        out << n_nodes << " "; // number of nodes
        for (int node : nodes[i])
            out << node << " ";
        out << "\n";
    }
}

} // namespace femocs
