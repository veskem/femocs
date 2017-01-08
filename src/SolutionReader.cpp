/*
 * SolutionReader.cpp
 *
 *  Created on: 6.6.2016
 *      Author: veske
 */

#include "SolutionReader.h"
#include "Macros.h"

#include <float.h>

using namespace std;
namespace femocs {

// Initialize data vectors
SolutionReader::SolutionReader() : interpolator(NULL) {
    reserve(0);
}

SolutionReader::SolutionReader(Interpolator* ip) : interpolator(ip) {
    reserve(0);
}

// Linearly interpolate solution on Medium atoms
const void SolutionReader::interpolate(const Medium &medium) {
    const int n_atoms = medium.get_n_atoms();
    reserve(n_atoms);   
    
    // Copy the atoms
    for (int i = 0; i < n_atoms; ++i)
        add_atom(medium.get_atom(i));
    
    // Sort atoms into sequential order to speed up interpolation
    sort_spatial();
    
    int elem = 0;
    for (int i = 0; i < n_atoms; ++i) {
        Point3 point = get_point(i);
        // Find the element that contains (elem >= 0) or is closest (elem < 0) to the point
        elem = interpolator->locate_element(point, abs(elem));

        // Store whether the point is in- or outside the mesh
        if (elem < 0) set_coordination(i, 1);
        else          set_coordination(i, 0);

        // Calculate the interpolation
        interpolation.push_back( interpolator->get_interpolation(point, abs(elem)) );
    }
}

// Linearly interpolate electric field on a set of points
const void SolutionReader::export_elfield(int n_points, double* x, double* y, double* z,
        double* Ex, double* Ey, double* Ez, double* Enorm, int* flag) {

    Medium medium(n_points);
    for (int i = 0; i < n_points; ++i)
        medium.add_atom(Point3(x[i], y[i], z[i]));

    // Interpolate solution and store the index of first point outside the mesh
    interpolate(medium);

    for (int i = 0; i < n_points; ++i) {
        Ex[i] = interpolation[i].elfield.x;
        Ey[i] = interpolation[i].elfield.y;
        Ez[i] = interpolation[i].elfield.z;
        Enorm[i] = interpolation[i].el_norm;
        flag[i] = atoms[i].coord;
    }
}

// Linearly interpolate electric potential on a set of points
const void SolutionReader::export_potential(int n_points, double* x, double* y, double* z, double* phi, int* flag) {
    Medium medium(n_points);
    for (int i = 0; i < n_points; ++i)
        medium.add_atom(Point3(x[i], y[i], z[i]));

    // Interpolate solution and store the index of first point outside the mesh
    interpolate(medium);

    for (int i = 0; i < n_points; ++i) {
        phi[i] = interpolation[i].potential;
        flag[i] = atoms[i].coord;
    }
}

// Export interpolated electric field
const void SolutionReader::export_solution(int n_atoms, double* Ex, double* Ey, double* Ez, double* Enorm) {
    // Initially pass the zero electric field for all the atoms
    for (int i = 0; i < n_atoms; ++i) {
        Ex[i] = 0;
        Ey[i] = 0;
        Ez[i] = 0;
        Enorm[i] = 0;
    }

    // Pass the the calculated electric field for stored atoms
    for (int i = 0; i < get_n_atoms(); ++i) {
        int identifier = get_id(i);
        if (identifier < 0 || identifier >= n_atoms) continue;

        Ex[identifier] = interpolation[i].elfield.x;
        Ey[identifier] = interpolation[i].elfield.y;
        Ez[identifier] = interpolation[i].elfield.z;
        Enorm[identifier] = interpolation[i].el_norm;
    }
}

// Get average electric field around I-th solution point
const Solution SolutionReader::get_average_solution(const int I, const double smooth_factor, const double r_cut) {
    const double r_cut2 = r_cut * r_cut;
    const double a = 1.0;

    Vec3 elfield(0.0);
    double potential = 0.0;

    Point3 point1 = get_point(I);
    double w_sum = 0.0;

    // E_I = sum i!=I [E_i * a*exp( -b*distance(i,I) )] / sum i!=I [a*exp( -b*distance(i,I) )]
    for (int i = 0; i < get_n_atoms(); ++i)
        if (i != I) {
            double dist2 = point1.distance2(get_point(i));
            if (dist2 > r_cut2 || interpolation[i].el_norm >= interpolator->error_field) continue;

            double w = a * exp(-1.0 * sqrt(dist2) / smooth_factor);
            w_sum += w;
            elfield += interpolation[i].elfield * w;
            potential += interpolation[i].potential * w;
        }

    elfield *= (1.0 / w_sum); potential *= (1.0 / w_sum);
    return Solution(elfield, elfield.norm(), potential);
}

// Get histogram for electric field x,y,z component or for its norm
const void SolutionReader::get_histogram(vector<int> &bins, vector<double> &bounds, const int coordinate) {
    require(coordinate >= 0 && coordinate <= 4, "Invalid component: " + to_string(coordinate));

    const int n_atoms = get_n_atoms();
    const int n_bins = bins.size();
    const int n_bounds = bounds.size();
    const double eps = 1e-5;

    // Find minimum and maximum values from all non-error values
    double value_min = DBL_MAX;
    double value_max =-DBL_MAX;
    double value;
    for (int i = 0; i < n_atoms; ++i) {
        if (coordinate == 4) value = interpolation[i].potential;
        else if (coordinate == 3) value = interpolation[i].el_norm;
        else                 value = interpolation[i].elfield[coordinate];

        if (abs(value) < interpolator->error_field) {
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
            if (coordinate == 4) value = interpolation[i].potential;
            else if (coordinate == 3) value = interpolation[i].el_norm;
            else                 value = interpolation[i].elfield[coordinate];

            if (value >= bounds[j] && value < bounds[j+1]) {
                bins[j]++;
                continue;
            }
        }
}

// Clean the interpolation from peaks
const void SolutionReader::clean(const int coordinate, const int n_bins, const double smooth_factor, const double r_cut) {
    require(coordinate >= 0 && coordinate <= 4, "Invalid coordinate: " + to_string(coordinate));
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
        if (bounds[i+1] >= 0) break;
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
        if (coordinate == 4) value = abs(interpolation[i].potential);
        else if (coordinate == 3) value = abs(interpolation[i].el_norm);
        else                 value = abs(interpolation[i].elfield[coordinate]);

        if (value < value_min || value > value_max)
            interpolation[i] = get_average_solution(i, smooth_factor, r_cut);
    }
}

// Reserve memory for solution vectors
const void SolutionReader::reserve(const int n_nodes) {
    require(n_nodes >= 0, "Invalid number of nodes: " + to_string(n_nodes));

    atoms.clear();
    interpolation.clear();

    atoms.reserve(n_nodes);
    interpolation.reserve(n_nodes);
}

// Compile data string from the data vectors for file output
const string SolutionReader::get_data_string(const int i) const{
    if (i < 0) return "SolutionReader data: id x y z out_of_mesh Ex Ey Ez Enorm potential";

    ostringstream strs;
    strs << atoms[i] << " " << interpolation[i];
    return strs.str();
}

const void SolutionReader::get_point_data(ofstream& out) const {
    const int n_atoms = get_n_atoms();

    // output vector magnitude explicitly to make it possible to apply filters in ParaView
    out << "SCALARS Elfield.norm double\nLOOKUP_TABLE default\n";
    for (size_t i = 0; i < n_atoms; ++i)
        out << interpolation[i].el_norm << "\n";

    // output scalar (electric potential, temperature etc)
    out << "SCALARS Potential double\nLOOKUP_TABLE default\n";
    for (size_t i = 0; i < n_atoms; ++i)
        out << interpolation[i].potential << "\n";

    // output vector data (electric field, current density etc)
    out << "VECTORS Electric_field double\n";
    for (size_t i = 0; i < n_atoms; ++i)
        out << interpolation[i].elfield << "\n";
}

} // namespace femocs
