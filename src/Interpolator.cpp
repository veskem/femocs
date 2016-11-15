/*
 * Interpolator.cpp
 *
 *  Created on: 19.10.2016
 *      Author: veske
 */

#include "Interpolator.h"
#include "Macros.h"

#include <fstream>
#include <float.h>

using namespace std;
namespace femocs {

Interpolator::Interpolator() : solution(NULL), mesh(NULL) {
    reserve(0);
    reserve_precompute(0);
};

Interpolator::Interpolator(SolutionReader* sr) : solution(sr), mesh(sr->get_mesh()) {
    reserve(0);
    reserve_precompute(0);
};

// Get average electric field around I-th solution point
const Vec3 Interpolator::get_average_solution(const int I, const double smooth_factor, const double r_cut) {
    const double r_cut2 = r_cut * r_cut;
    const double a = 1.0;

    Vec3 average(0.0);
    Point3 point1 = get_point(I);
    double w_sum = 0.0;

    for (int i = 0; i < get_n_atoms(); ++i)
        if (i != I) {
            double dist2 = point1.distance2(get_point(i));
            if (dist2 > r_cut2 || interpolation[i].el_norm >= error_field) continue;

            double w = a * exp(-1.0 * sqrt(dist2) / smooth_factor);
            w_sum += w;
            average += interpolation[i].elfield * w;
        }

    average *= (1.0 / w_sum);
    return average;
}

// Get histogram for electric field x,y,z component or for its norm
const void Interpolator::get_histogram(vector<int> &bins, vector<double> &bounds, const int coordinate) {
    require(coordinate >= 0 && coordinate <= 3, "Invalid component: " + to_string(coordinate));

    const int n_atoms = get_n_atoms();
    const int n_bins = bins.size();
    const int n_bounds = bounds.size();
    const double eps = 1e-5;

    // Find minimum and maximum values from all non-error values
    double value_min = DBL_MAX;
    double value_max =-DBL_MAX;
    double value;
    for (int i = 0; i < n_atoms; ++i) {
        if (coordinate == 3) value = interpolation[i].el_norm;
        else                 value = interpolation[i].elfield[coordinate];

        if (abs(value) < error_field) {
            value_min = min(value_min, value);
            value_max = max(value_max, value);
        }
    }

    // Fill the bounds with values value_min:value_step:(value_max + epsilon)
    // Epsilon is added to value_max to include the maximum value in the up-most bin
    double value_step = (value_max - value_min) / n_bins;
    for (int i = 0; i < n_bounds; ++i)
        bounds[i] = value_min + value_step * i;
    bounds[n_bounds-1] += eps;

    for (int i = 0; i < n_atoms; ++i)
        for (int j = 0; j < n_bins; ++j) {
            if (coordinate == 3) value = interpolation[i].el_norm;
            else                 value = interpolation[i].elfield[coordinate];

            if (value >= bounds[j] && value < bounds[j+1]) {
                bins[j]++;
                continue;
            }
        }
}

// Function to clean the result from peaks
const void Interpolator::clean(const int coordinate, const int n_bins, const double smooth_factor, const double r_cut) {
    require(coordinate >= 0 && coordinate <= 3, "Invalid coordinate: " + to_string(coordinate));
    if (n_bins <= 1) return;

    const int n_atoms = get_n_atoms();

    vector<int> bins(n_bins, 0);
    vector<double> bounds(n_bins+1);
    get_histogram(bins, bounds, coordinate);

    // Find the first bin with zero entries from positive edge of bounds;
    // this will determine the maximum allowed elfield value
    double value_max = bounds[n_bins];
    for (int i = n_bins-1; i >= 0; --i) {
        if (bounds[i] < 0) break;
        if (bins[i] == 0) value_max = bounds[i];
    }

    // Find the last bin with zero entries from negative edge of bounds;
    // this will determine the minimum allowed elfield value
    double value_min = bounds[0];
    for (int i = 0; i < n_bins; ++i) {
        if (bounds[i+1] > 0) break;
        if (bins[i] == 0) value_min = bounds[i+1];
    }

    require(value_min <= value_max, "Error in histogram cleaner!");

//    cout.precision(3);
//    cout << endl << coordinate << " " << value_min << " " << value_max << endl;
//    for (int i = 0; i < bins.size(); ++i) cout << bins[i] << " ";
//    cout << endl;
//    for (int i = 0; i < bounds.size(); ++i) cout << bounds[i] << " ";
//    cout << endl;

    // If all the bins are filled, no blocking will be applied
    if (value_min == bounds[0] && value_max == bounds[n_bins])
        return;

    double value;
    for (int i = 0; i < n_atoms; ++i) {
        if (coordinate == 3) value = abs(interpolation[i].el_norm);
        else                 value = abs(interpolation[i].elfield[coordinate]);

        if (value < value_min || value > value_max) {
            interpolation[i].elfield = get_average_solution(i, smooth_factor, r_cut);
            interpolation[i].el_norm = interpolation[i].elfield.norm();
        }
    }
}

// Compile data string from the data vectors for file output
const string Interpolator::get_data_string(const int i) {
    if (i < 0) return "Interpolator data: id x y z coordination Ex Ey Ez Enorm potential";

    ostringstream strs;
    strs << atoms[i] << " " << interpolation[i];
    return strs.str();
}

// Reserve memory for interpolation data
const void Interpolator::reserve(const int N) {
    atoms.clear();
    interpolation.clear();

    Medium::reserve(N);
    interpolation.reserve(N);
}

// Reserve memory for pre-compute data
const void Interpolator::reserve_precompute(const int N) {
    centroid.clear();
    det0.clear();
    det1.clear();
    det2.clear();
    det3.clear();
    det4.clear();
    tet_not_valid.clear();

    centroid.reserve(N);
    det0.reserve(N);
    det1.reserve(N);
    det2.reserve(N);
    det3.reserve(N);
    det4.reserve(N);
    tet_not_valid.reserve(N);
}

// Precompute the data to tetrahedra to make later bcc calculations faster
const void Interpolator::precompute_tetrahedra() {
    const int n_elems = mesh->elems.size();
    double d0, d1, d2, d3, d4;

    reserve_precompute(n_elems);

    // Calculate centroids of elements
    for (int i = 0; i < n_elems; ++i)
        centroid.push_back(mesh->elems.get_centroid(i));

    /* Calculate main and minor determinants for 1st, 2nd, 3rd and 4th
     * barycentric coordinate of tetrahedra using the relations below */
    for (SimpleElement se : mesh->elems) {
        Vec3 v1 = mesh->nodes.get_vec(se[0]);
        Vec3 v2 = mesh->nodes.get_vec(se[1]);
        Vec3 v3 = mesh->nodes.get_vec(se[2]);
        Vec3 v4 = mesh->nodes.get_vec(se[3]);

        /* =====================================================================================
         * det0 = |x1 y1 z1 1|
                  |x2 y2 z2 1|
                  |x3 y3 z3 1|
                  |x4 y4 z4 1|  */
        d0 = determinant(v1, v2, v3, v4);
        tet_not_valid.push_back(fabs(d0) < epsilon);
        det0.push_back(1.0 / d0);

        /* =====================================================================================
         * det1 = |x  y  z  1| = + x * |y2 z2 1| - y * |x2 z2 1| + z * |x2 y2 1| - |x2 y2 z2|
                  |x2 y2 z2 1|         |y3 z3 1|       |x3 z3 1|       |x3 y3 1|   |x3 y3 z3|
                  |x3 y3 z3 1|         |y4 z4 1|       |x4 z4 1|       |x4 y4 1|   |x4 y4 z4|
                  |x4 y4 z4 1|  */
        d1 = determinant(Vec3(v2.y, v3.y, v4.y), Vec3(v2.z, v3.z, v4.z));
        d2 = determinant(Vec3(v2.x, v3.x, v4.x), Vec3(v2.z, v3.z, v4.z));
        d3 = determinant(Vec3(v2.x, v3.x, v4.x), Vec3(v2.y, v3.y, v4.y));
        d4 = determinant(Vec3(v2.x, v3.x, v4.x), Vec3(v2.y, v3.y, v4.y), Vec3(v2.z, v3.z, v4.z));

        det1.push_back(Vec4(d1, -d2, d3, -d4));

        /* =====================================================================================
         * det2 = |x1 y1 z1 1| = - x * |y1 z1 1| + y * |x1 z1 1| - z * |x1 y1 1| + |x1 y1 z1|
                  |x  y  z  1|         |y3 z3 1|       |x3 z3 1|       |x3 y3 1|   |x3 y3 z3|
                  |x3 y3 z3 1|         |y4 z4 1|       |x4 z4 1|       |x4 y4 1|   |x4 y4 z4|
                  |x4 y4 z4 1|  */
        d1 = determinant(Vec3(v1.y, v3.y, v4.y), Vec3(v1.z, v3.z, v4.z));
        d2 = determinant(Vec3(v1.x, v3.x, v4.x), Vec3(v1.z, v3.z, v4.z));
        d3 = determinant(Vec3(v1.x, v3.x, v4.x), Vec3(v1.y, v3.y, v4.y));
        d4 = determinant(Vec3(v1.x, v3.x, v4.x), Vec3(v1.y, v3.y, v4.y), Vec3(v1.z, v3.z, v4.z));

        det2.push_back(Vec4(-d1, d2, -d3, d4));

        /* =====================================================================================
         * det3 = |x1 y1 z1 1| = + x * |y1 z1 1| - y * |x1 z1 1| + z * |x1 y1 1| - |x1 y1 z1|
                  |x2 y2 z2 1|         |y2 z2 1|       |x2 z2 1|       |x2 y2 1|   |x2 y2 z2|
                  |x  y  z  1|         |y4 z4 1|       |x4 z4 1|       |x4 y4 1|   |x4 y4 z4|
                  |x4 y4 z4 1|  */
        d1 = determinant(Vec3(v1.y, v2.y, v4.y), Vec3(v1.z, v2.z, v4.z));
        d2 = determinant(Vec3(v1.x, v2.x, v4.x), Vec3(v1.z, v2.z, v4.z));
        d3 = determinant(Vec3(v1.x, v2.x, v4.x), Vec3(v1.y, v2.y, v4.y));
        d4 = determinant(Vec3(v1.x, v2.x, v4.x), Vec3(v1.y, v2.y, v4.y), Vec3(v1.z, v2.z, v4.z));

        det3.push_back(Vec4(d1, -d2, d3, -d4));

        /* =====================================================================================
         * det4 = |x1 y1 z1 1| = - x * |y1 z1 1| + y * |x1 z1 1| - z * |x1 y1 1| + |x1 y1 z1|
                  |x2 y2 z2 1|         |y2 z2 1|       |x2 z2 1|       |x2 y2 1|   |x2 y2 z2|
                  |x3 y3 z3 1|         |y3 z3 1|       |x3 z3 1|       |x3 y3 1|   |x3 y3 z3|
                  |x  y  z  1|  */

        d1 = determinant(Vec3(v1.y, v2.y, v3.y), Vec3(v1.z, v2.z, v3.z));
        d2 = determinant(Vec3(v1.x, v2.x, v3.x), Vec3(v1.z, v2.z, v3.z));
        d3 = determinant(Vec3(v1.x, v2.x, v3.x), Vec3(v1.y, v2.y, v3.y));
        d4 = determinant(v1, v2, v3);

        det4.push_back(Vec4(-d1, d2, -d3, d4));
    }
}

// Check with barycentric coordinates whether the point is inside the i-th tetrahedron
const bool Interpolator::point_in_tetrahedron(const Point3 &point, const int i) {
    expect(i >= 0 && i < det0.size(), "Index out of bounds: " + to_string(i));

    // Ignore co-planar tetrahedra
//    if (tet_not_valid[i]) return false;

    Vec4 pt(point, 1);

    // If one of the barycentric coordinates is < zero, the point is outside the tetrahedron
    // Source: http://steve.hollasch.net/cgindex/geometry/ptintet.html
    if (det0[i] * pt.dotProduct(det1[i]) < zero) return false;
    if (det0[i] * pt.dotProduct(det2[i]) < zero) return false;
    if (det0[i] * pt.dotProduct(det3[i]) < zero) return false;
    if (det0[i] * pt.dotProduct(det4[i]) < zero) return false;

    // All bcc-s are >= 0, so point is inside the tetrahedron
    return true;
}

// Calculate barycentric coordinates for point
const Vec4 Interpolator::get_bcc(const Point3 &point, const int elem) {
    expect(elem >= 0 && elem < det0.size(), "Index out of bounds: " + to_string(elem));

    Vec4 pt(point, 1);

    double bcc1 = det0[elem] * pt.dotProduct(det1[elem]);
    double bcc2 = det0[elem] * pt.dotProduct(det2[elem]);
    double bcc3 = det0[elem] * pt.dotProduct(det3[elem]);
    double bcc4 = det0[elem] * pt.dotProduct(det4[elem]);

    return Vec4(bcc1, bcc2, bcc3, bcc4);
}

// Find the element which contains the point or is the closest to it
const int Interpolator::locate_element(const Point3 &point, const int elem_guess) {
    // Check first the guessed element
    if (point_in_tetrahedron(point, elem_guess)) return elem_guess;

    // Check then the neighbours of guessed element
    for (int nbr : mesh->elems.get_neighbours(elem_guess))
        if (nbr >= 0 && point_in_tetrahedron(point, nbr))
            return nbr;

    const int n_elems = det0.size();
    double min_distance2 = DBL_MAX;
    int min_index = 0;

    // If no success, loop through all the elements
    for (int elem = 0; elem < n_elems; ++elem) {
        // If correct element is found, we're done
        if (point_in_tetrahedron(point, elem)) return elem;

        // Otherwise look for the element whose centroid distance to the point is the smallest
        else {
            double distance2 = point.distance2(centroid[elem]);
            if (distance2 < min_distance2) {
                min_distance2 = distance2;
                min_index = elem;
            }
        }
    }

    // If no perfect element found, return the best; indicate the imperfectness with the minus sign
    return -min_index;
}

// Calculate interpolation for point inside or near the elem-th tetrahedron
const Solution Interpolator::get_interpolation(const Point3 &point, const int elem) {
    expect(elem >= 0 && elem < mesh->elems.size(), "Index out of bounds: " + to_string(elem));

    // Get barycentric coordinates of point in tetrahedron
    Vec4 bcc = get_bcc(point, elem);

    SimpleElement selem = mesh->elems[elem];

    // Interpolate electric field
    Vec3 elfield_i(0.0);
    for (int i = 0; i < mesh->elems.DIM; ++i)
        elfield_i += solution->get_solution(selem[i]).elfield * bcc[i];

    // Interpolate potential
    double potential_i(0.0);
    for (int i = 0; i < mesh->elems.DIM; ++i)
        potential_i += solution->get_solution(selem[i]).potential * bcc[i];

    return Solution(elfield_i, elfield_i.norm(), potential_i);
}

// Linearly interpolate solution on Medium atoms
const void Interpolator::extract_interpolation(const Medium &medium) {
    const int n_atoms = medium.get_n_atoms();
    reserve(n_atoms);

    int elem = 0;

    for (int i = 0; i < n_atoms; ++i) {
        Atom atom = medium.get_atom(i);
        add_atom(atom);

        // Find the element that contains (elem >= 0) or is closest (elem < 0) to the point
        elem = locate_element(atom.point, abs(elem));

        // Store the data whether the point is in- of outside the mesh
        if (elem < 0) set_coordination(i, 1);
        else          set_coordination(i, 0);

        // Calculate the interpolation
        interpolation.push_back( get_interpolation(atom.point, abs(elem)) );
    }
}

// Linearly interpolate electric field on the set of points
const void Interpolator::extract_elfield(int n_points, double* x, double* y, double* z,
        double* Ex, double* Ey, double* Ez, double* Enorm, int* flag) {

    Medium medium(n_points);
    for (int i = 0; i < n_points; ++i)
        medium.add_atom(Point3(x[i], y[i], z[i]));

    // Interpolate solution and store the index of first point outside the mesh
    extract_interpolation(medium);

    for (int i = 0; i < n_points; ++i) {
        Ex[i] = interpolation[i].elfield.x;
        Ey[i] = interpolation[i].elfield.y;
        Ez[i] = interpolation[i].elfield.z;
        Enorm[i] = interpolation[i].el_norm;
        flag[i] = atoms[i].coord;
    }
}

// Linearly interpolate electric potential on the set of points
const void Interpolator::extract_potential(int n_points, double* x, double* y, double* z, double* phi, int* flag) {
    Medium medium(n_points);
    for (int i = 0; i < n_points; ++i)
        medium.add_atom(Point3(x[i], y[i], z[i]));

    // Interpolate solution and store the index of first point outside the mesh
    extract_interpolation(medium);

    for (int i = 0; i < n_points; ++i) {
        phi[i] = interpolation[i].potential;
        flag[i] = atoms[i].coord;
    }
}

// Determinant of 3x3 matrix which's last column consists of ones
const double Interpolator::determinant(const Vec3 &v1, const Vec3 &v2) {
    return v1.x * (v2.y - v2.z)
         - v1.y * (v2.x - v2.z)
         + v1.z * (v2.x - v2.y);
}

// Determinant of 3x3 matrix which's columns consist of Vec3-s
const double Interpolator::determinant(const Vec3 &v1, const Vec3 &v2, const Vec3 &v3) {
    return v1.x * (v2.y * v3.z - v3.y * v2.z)
         - v2.x * (v1.y * v3.z - v3.y * v1.z)
         + v3.x * (v1.y * v2.z - v2.y * v1.z);
}

// Determinant of 4x4 matrix which's last column consists of ones
const double Interpolator::determinant(const Vec3 &v1, const Vec3 &v2, const Vec3 &v3, const Vec3 &v4) {
    const double det1 = determinant(v2, v3, v4);
    const double det2 = determinant(v1, v3, v4);
    const double det3 = determinant(v1, v2, v4);
    const double det4 = determinant(v1, v2, v3);

    return det4 - det3 + det2 - det1;
}

// Determinant of 4x4 matrix which's columns consist of Vec4-s
const double Interpolator::determinant(const Vec4 &v1, const Vec4 &v2, const Vec4 &v3, const Vec4 &v4) {
    double det1 = determinant(Vec3(v1.y,v1.z,v1.w), Vec3(v2.y,v2.z,v2.w), Vec3(v3.y,v3.z,v3.w));
    double det2 = determinant(Vec3(v1.x,v1.z,v1.w), Vec3(v2.x,v2.z,v2.w), Vec3(v3.x,v3.z,v3.w));
    double det3 = determinant(Vec3(v1.x,v1.y,v1.w), Vec3(v2.x,v2.y,v2.w), Vec3(v3.x,v3.y,v3.w));
    double det4 = determinant(Vec3(v1.x,v1.y,v1.z), Vec3(v2.x,v2.y,v2.z), Vec3(v3.x,v3.y,v3.z));

    return v4.w * det4 - v4.z * det3 + v4.y * det2 - v4.x * det1;
}

// Export interpolated electric field to helmod
const void Interpolator::export_helmod(int n_atoms, double* Ex, double* Ey, double* Ez, double* Enorm) {
    //require(n_atoms <= get_n_atoms(), "Solution vector is shorter than requested: " + to_string(get_n_atoms()));

    // Initially pass the zero electric field for all the atoms
    for (int i = 0; i < n_atoms; ++i) {
        Ex[i] = 0;
        Ey[i] = 0;
        Ez[i] = 0;
        Enorm[i] = 0;
    }

    // Pass the the real electric field for surface atoms that were in the mesh
    for (int i = 0; i < get_n_atoms(); ++i) {
        int identifier = get_id(i);
        if (identifier < 0) continue;


        require(identifier < n_atoms, "Mismatch between import and export data sizes detected!");

        Ex[identifier] = interpolation[i].elfield.x;
        Ey[identifier] = interpolation[i].elfield.y;
        Ez[identifier] = interpolation[i].elfield.z;
        Enorm[identifier] = interpolation[i].el_norm;
    }
}

} /* namespace femocs */
