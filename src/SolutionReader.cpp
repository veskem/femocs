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

/* ==========================================
 * ============= SOLUTION READER ============
 * ========================================== */

// Initialize SolutionReader
SolutionReader::SolutionReader() : interpolator(NULL), vec_label("vec"), vec_norm_label("vec_norm"), scalar_label("scalar") {
    reserve(0);
}

SolutionReader::SolutionReader(LinearInterpolator* ip, const string& vec_lab, const string& vec_norm_lab, const string& scal_lab) :
        interpolator(ip), vec_label(vec_lab), vec_norm_label(vec_norm_lab), scalar_label(scal_lab) {
    reserve(0);
}

// Reserve memory for solution vectors
void SolutionReader::reserve(const int n_nodes) {
    require(n_nodes >= 0, "Invalid number of nodes: " + to_string(n_nodes));

    atoms.clear();
    interpolation.clear();

    atoms.reserve(n_nodes);
    interpolation.reserve(n_nodes);
}

// Compile data string from the data vectors for file output
string SolutionReader::get_data_string(const int i) const{
    if (i < 0) return "SolutionReader properties=id:I:1:pos:R:3:marker:I:1:force:R:3:" + vec_norm_label + ":R:1:" + scalar_label + ":R:1";

    ostringstream strs; strs << fixed;
    strs << atoms[i] << " " << interpolation[i];
    return strs.str();
}

void SolutionReader::get_point_data(ofstream& out) const {
    const int n_atoms = size();

    // write IDs of atoms
    out << "SCALARS id int\nLOOKUP_TABLE default\n";
    for (size_t i = 0; i < n_atoms; ++i)
        out << atoms[i].id << "\n";

    // write atom markers
    out << "SCALARS marker int\nLOOKUP_TABLE default\n";
    for (size_t i = 0; i < n_atoms; ++i)
        out << atoms[i].marker << "\n";

    // output scalar (electric potential, temperature etc)
    out << "SCALARS " << scalar_label << " double\nLOOKUP_TABLE default\n";
    for (size_t i = 0; i < n_atoms; ++i)
        out << interpolation[i].scalar << "\n";

    // output vector magnitude explicitly to make it possible to apply filters in ParaView
    out << "SCALARS " << vec_norm_label << " double\nLOOKUP_TABLE default\n";
    for (size_t i = 0; i < n_atoms; ++i)
        out << interpolation[i].norm << "\n";

    // output vector data (electric field, current density etc)
    out << "VECTORS " << vec_label << " double\n";
    for (size_t i = 0; i < n_atoms; ++i)
        out << interpolation[i].vector << "\n";
}

// Get average electric field around I-th solution point
Solution SolutionReader::get_average_solution(const int I, const double r_cut) {
    // Cut off weights after 5 sigma
    const double r_cut2 = pow(5 * r_cut, 2);
    const double smooth_factor = r_cut / 5.0;

    Vec3 elfield(0.0);
    double potential = 0.0;

    Point3 point1 = get_point(I);
    double w_sum = 0.0;

    // E_I = sum i!=I [E_i * a*exp( -b*distance(i,I) )] / sum i!=I [a*exp( -b*distance(i,I) )]
    for (int i = 0; i < size(); ++i)
        if (i != I) {
            double dist2 = point1.distance2(get_point(i));
            if (dist2 > r_cut2 || interpolation[i].norm >= interpolator->error_field) continue;

            double w = exp(-1.0 * sqrt(dist2) / smooth_factor);
            w_sum += w;
            elfield += interpolation[i].vector * w;
            potential += interpolation[i].scalar * w;
        }

    elfield *= (1.0 / w_sum); potential *= (1.0 / w_sum);
    return Solution(elfield, potential);
}

// Get histogram for electric field x,y,z component or for its norm
void SolutionReader::get_histogram(vector<int> &bins, vector<double> &bounds, const int coordinate) {
    require(coordinate >= 0 && coordinate <= 4, "Invalid component: " + to_string(coordinate));

    const int n_atoms = size();
    const int n_bins = bins.size();
    const int n_bounds = bounds.size();
    const double eps = 1e-5;

    // Find minimum and maximum values from all non-error values
    double value_min = DBL_MAX;
    double value_max =-DBL_MAX;
    double value;
    for (int i = 0; i < n_atoms; ++i) {
        if (coordinate == 4) value = interpolation[i].scalar;
        else if (coordinate == 3) value = interpolation[i].norm;
        else                 value = interpolation[i].vector[coordinate];

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
            if (coordinate == 4) value = interpolation[i].scalar;
            else if (coordinate == 3) value = interpolation[i].norm;
            else                 value = interpolation[i].vector[coordinate];

            if (value >= bounds[j] && value < bounds[j+1]) {
                bins[j]++;
                continue;
            }
        }
}

// Clean the interpolation from peaks using histogram cleaner
void SolutionReader::clean(const int coordinate, const double r_cut) {
    require(coordinate >= 0 && coordinate <= 4, "Invalid coordinate: " + to_string(coordinate));
    const int n_bins = (int) size() / 250;

    if (n_bins <= 1 || r_cut < 0.1) return;

    const int n_atoms = size();

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
        if (coordinate == 4) value = abs(interpolation[i].scalar);
        else if (coordinate == 3) value = abs(interpolation[i].norm);
        else                 value = abs(interpolation[i].vector[coordinate]);

        if (value < value_min || value > value_max)
            interpolation[i] = get_average_solution(i, r_cut);
    }
}

/* ==========================================
 * ============== FIELD READER ==============
 * ========================================== */

FieldReader::FieldReader() : SolutionReader() {}
FieldReader::FieldReader(LinearInterpolator* ip) : SolutionReader(ip, "elfield", "elfield_norm", "potential") {}

// Print statistics about interpolated solution
void FieldReader::print_statistics() {
    if (!MODES.VERBOSE) return;

    const int n_atoms = size();
    Vec3 elfield(0);
    Vec3 rms_elfield(0);
    double potential = 0;
    double rms_potential = 0;
    int n_points = 0;

    for (int i = 0; i < n_atoms; ++i) {
        double pot = interpolation[i].scalar;
        if (pot >= interpolator->error_field) continue;
        Vec3 ef = interpolation[i].vector;

        elfield += ef;
        rms_elfield += ef * ef;
        potential += pot;
        rms_potential += pot * pot;
        n_points++;
    }

    elfield *= (1.0 / n_points);
    rms_elfield = Vec3(sqrt(rms_elfield.x), sqrt(rms_elfield.y), sqrt(rms_elfield.z)) * (1.0 / n_points);
    potential = potential / n_points;
    rms_potential = sqrt(rms_potential) / n_points;

    cout << "  mean elfield: \t" << elfield << endl;
    cout << "   rms elfield: \t" << rms_elfield << endl;
    cout << "  mean potential: \t" << potential << endl;
    cout << "   rms potential: \t" << rms_potential << endl;
}

// Linearly interpolate solution on Medium atoms
void FieldReader::interpolate(const Medium &medium, const double r_cut, const int component, const bool srt) {
    require(component >= 0 && component <= 2, "Invalid interpolation component: " + to_string(component));
    require(interpolator, "NULL interpolator cannot be used!");

    const int n_atoms = medium.size();
    reserve(n_atoms);   
    
    // Copy the atoms
    for (int i = 0; i < n_atoms; ++i)
        append(Atom(i, medium.get_point(i), 0));

    // Sort atoms into sequential order to speed up interpolation
    if (srt) sort_spatial();

    int elem = 0;
    for (int i = 0; i < n_atoms; ++i) {
        Point3 point = get_point(i);
        // Find the element that contains (elem >= 0) or is closest (elem < 0) to the point
        elem = interpolator->locate_element(point, abs(elem));

        // Store whether the point is in- or outside the mesh
        if (elem < 0) set_marker(i, 1);
        else          set_marker(i, 0);

        // Calculate the interpolation
        if (component == 0) interpolation.push_back( interpolator->get_solution(point, abs(elem)) );
        if (component == 1) interpolation.push_back( interpolator->get_vector(point, abs(elem)) );
        if (component == 2) interpolation.push_back( interpolator->get_scalar(point, abs(elem)) );
    }

    clean(0, r_cut);  // clean by vector x-component
    clean(1, r_cut);  // clean by vector y-component
    clean(2, r_cut);  // clean by vector z-component
    clean(3, r_cut);  // clean by vector norm
    clean(4, r_cut);  // clean by scalar

    if (srt) {
        // sort atoms back to their initial order
        for (int i = 0; i < n_atoms; ++i)
            interpolation[i].id = atoms[i].id;
        sort(interpolation.begin(), interpolation.end(), Solution::sort_up());
        sort(atoms.begin(), atoms.end(), Atom::sort_id());
    }
}

// Linearly interpolate electric field on a set of points
void FieldReader::interpolate(const int n_points, const double* x, const double* y, const double* z,
        const double r_cut, const int component, const bool sort) {
    Medium medium(n_points);
    for (int i = 0; i < n_points; ++i)
        medium.append(Point3(x[i], y[i], z[i]));

    // Interpolate solution
    interpolate(medium, r_cut, component, sort);
}

// Linearly interpolate electric field on a set of points
void FieldReader::export_elfield(const int n_points, double* Ex, double* Ey, double* Ez, double* Enorm, int* flag) {
    require(n_points == size(), "Invalid array size: " + to_string(n_points));

    for (int i = 0; i < n_points; ++i) {
        Ex[i] = interpolation[i].vector.x;
        Ey[i] = interpolation[i].vector.y;
        Ez[i] = interpolation[i].vector.z;
        Enorm[i] = interpolation[i].norm;
        flag[i] = atoms[i].marker;
    }
}

// Linearly interpolate electric potential on a set of points
void FieldReader::export_potential(const int n_points, double* phi, int* flag) {
    require(n_points == size(), "Invalid array size: " + to_string(n_points));

    for (int i = 0; i < n_points; ++i) {
        phi[i] = interpolation[i].scalar;
        flag[i] = atoms[i].marker;
    }
}

// Export interpolated electric field
void FieldReader::export_solution(const int n_atoms, double* Ex, double* Ey, double* Ez, double* Enorm) {
    // Initially pass the zero electric field for all the atoms
    for (int i = 0; i < n_atoms; ++i) {
        Ex[i] = 0;
        Ey[i] = 0;
        Ez[i] = 0;
        Enorm[i] = 0;
    }

    // Pass the the calculated electric field for stored atoms
    for (int i = 0; i < size(); ++i) {
        int identifier = get_id(i);
        if (identifier < 0 || identifier >= n_atoms) continue;

        Ex[identifier] = interpolation[i].vector.x;
        Ey[identifier] = interpolation[i].vector.y;
        Ez[identifier] = interpolation[i].vector.z;
        Enorm[identifier] = interpolation[i].norm;
    }
}

Vec3 FieldReader::get_elfield(const int i) const {
    require(i >= 0 && i < size(), "Invalid index: " + to_string(i));
    return interpolation[i].vector;
}

double FieldReader::get_potential(const int i) const {
    require(i >= 0 && i < size(), "Invalid index: " + to_string(i));
    return interpolation[i].scalar;
}

/* ==========================================
 * =============== HEAT READER ==============
 * ========================================== */

HeatReader::HeatReader() : SolutionReader() {}
HeatReader::HeatReader(LinearInterpolator* ip) : SolutionReader(ip, "rho", "rho_norm", "temperature") {}

// Linearly interpolate solution on Medium atoms
void HeatReader::interpolate(const Medium &medium, const double r_cut) {
    require(interpolator, "NULL interpolator cannot be used!");

    const int n_atoms = medium.size();
    reserve(n_atoms);

    // Copy the atoms
    for (int i = 0; i < n_atoms; ++i)
        append(Atom(i, medium.get_point(i), 0));

    // Sort atoms into sequential order to speed up interpolation
    sort_spatial();

    int elem = 0;
    for (int i = 0; i < n_atoms; ++i) {
        Point3 point = get_point(i);
        // Find the element that contains (elem >= 0) or is closest (elem < 0) to the point
        elem = interpolator->locate_element(point, abs(elem));

        // Store whether the point is in- or outside the mesh
        if (elem < 0) set_marker(i, 1);
        else          set_marker(i, 0);

        // Calculate the interpolation
        interpolation.push_back( interpolator->get_solution(point, abs(elem)) );
    }

    // sort atoms back to their initial order
    for (int i = 0; i < n_atoms; ++i)
        interpolation[i].id = atoms[i].id;
    sort(interpolation.begin(), interpolation.end(), Solution::sort_up());
    sort(atoms.begin(), atoms.end(), Atom::sort_id());
}

Vec3 HeatReader::get_rho(const int i) const {
    require(i >= 0 && i < size(), "Invalid index: " + to_string(i));
    return interpolation[i].vector;
}

double HeatReader::get_temperature(const int i) const {
    require(i >= 0 && i < size(), "Invalid index: " + to_string(i));
    return interpolation[i].scalar;
}

/* ==========================================
 * ============== CHARGE READER =============
 * ========================================== */

ChargeReader::ChargeReader() : SolutionReader() {}
ChargeReader::ChargeReader(LinearInterpolator* ip) : SolutionReader(ip, "elfield", "area", "charge") {}

// Calculate charges on surface faces
void ChargeReader::calc_interpolated_charges(const TetgenMesh& mesh) {
    require(interpolator, "NULL interpolator cannot be used!");

    const int n_faces = mesh.faces.size();
    reserve(n_faces);

    // Store the centroids of the triangles
    for (int i = 0; i < n_faces; ++i)
        append( Atom(i, mesh.faces.get_centroid(i), 0) );

    // Sort centroids into sequential order to speed up interpolation
    sort_spatial();

    int elem = 0;
    for (int i = 0; i < n_faces; ++i) {
        Point3 point = get_point(i);
        int face = get_id(i);

        // Interpolate electric field in the centroid
        elem = interpolator->locate_element(point, abs(elem));
        Vec3 elfield = interpolator->get_vector(point, abs(elem));

        double area = mesh.faces.get_area(face);
        double charge = area * elfield.norm();  // * eps0;
        interpolation.push_back(Solution(elfield, area, charge));
    }

    // sort atoms back to their initial order
    for (int i = 0; i < n_faces; ++i)
        interpolation[i].id = atoms[i].id;

    sort(interpolation.begin(), interpolation.end(), Solution::sort_up());
    sort(atoms.begin(), atoms.end(), Atom::sort_id());
}

void ChargeReader::calc_charges(const TetgenMesh& mesh) {
    const int n_faces = mesh.faces.size();
    reserve(n_faces);

    // Store the centroids of the triangles
    for (int i = 0; i < n_faces; ++i)
        append( Atom(i, mesh.faces.get_centroid(i), 0) );

    // Find the electric fields in the centroids of the triangles
    vector<Vec3> elfields;
    get_elfields(mesh, elfields);
    require(elfields.size() == n_faces, "Electric fields were not extracted for every face!");

    // Calculate the charges for the triangles
    for (int face = 0; face < n_faces; ++face) {
        double area = mesh.faces.get_area(face);
        double charge = area * elfields[face].norm();  // * eps0;
        interpolation.push_back(Solution(elfields[face], area, charge));
    }
}

Vec3 ChargeReader::get_elfield(const int i) const {
    require(i >= 0 && i < size(), "Invalid index: " + to_string(i));
    return interpolation[i].vector;
}

double ChargeReader::get_area(const int i) const {
    require(i >= 0 && i < size(), "Invalid index: " + to_string(i));
    return interpolation[i].norm;
}

double ChargeReader::get_charge(const int i) const {
    require(i >= 0 && i < size(), "Invalid index: " + to_string(i));
    return interpolation[i].scalar;
}

void ChargeReader::print_statistics(const Medium::Sizes& sizes, const double E0) {
    if (!MODES.VERBOSE) return;

    double Q = E0 * sizes.xbox * sizes.ybox;
    double q = 0;
    for (int i = 0; i < size(); ++i)
        q += get_charge(i);

    cout << "  Q / sum(q) = " << Q << " / " << q << " = " << Q/q << endl;
}

void ChargeReader::get_elfields(const TetgenMesh& mesh, vector<Vec3> &elfields) {
    const int n_hexs = mesh.hexahedra.size();
    const int n_faces = mesh.faces.size();
    const int n_nodes = mesh.nodes.size();

    // Mark the nodes that are connected to the surface triangles
    vector<bool> on_face(n_nodes);
    for (int face = 0; face < n_faces; ++face)
        for (int node : mesh.faces[face])
            on_face[node] = true;

    // Store the indices of hexahedra connected to surface nodes
    vector<vector<int>> node2hex(n_nodes);
    for (int hex = 0; hex < n_hexs; ++hex)
        for (int node : mesh.hexahedra[hex]) {
            if (on_face[node])
                node2hex[node].push_back(hex);
        }

    // Extract the electric fields in the centroids of the triangles
    elfields.clear(); elfields.reserve(n_faces);
    for (int face = 0; face < n_faces; ++face) {
        Point3 centroid = get_point(face);
        int vert = mesh.faces[face][0];

        int centroid_indx = -1;
        double min_dist = 1e100;

        // Loop through all the hexahedra connected to the first vertex of triangle
        for (int hex : node2hex[vert])
            // Loop through all the vertices of hexahedron
            for (int node : mesh.hexahedra[hex]) {
                // Determine the node that is closest to the centroid of a face
                double dist = centroid.distance2(mesh.nodes[node]);
                if (dist < min_dist) {
                    min_dist = dist;
                    centroid_indx = node;
                }
            }

        elfields.push_back(interpolator->get_vector(centroid_indx));
    }
}

/* ==========================================
 * ============== FORCE READER ==============
 * ========================================== */

ForceReader::ForceReader() : SolutionReader() {}
ForceReader::ForceReader(LinearInterpolator* ip) : SolutionReader(ip, "force", "force_norm", "charge") {}

void ForceReader::calc_forces(const FieldReader &fields, const ChargeReader& faces, const Medium::Sizes& sizes, const double r_cut) {
    const int n_atoms = fields.size();
    const int n_faces = faces.size();

    // Copy the atom data
    reserve(n_atoms);
    for (int i = 0; i < n_atoms; ++i)
        append(fields.get_atom(i));

    /* Distribute the charges on surface faces between surface atoms.
     * If q_i and Q_j is the charge on i-th atom and j-th face, respectively, then
     *     q_i = sum_j(w_ij * Q_j),  sum_i(w_ij) = 1 for every j
     * where w_ij is the weight of charge on j-th face for the i-th atom. */

    vector<double> charges(n_atoms);
    vector<double> weights;
    for (int face = 0; face < n_faces; ++face) {
        Point3 point1 = faces.get_point(face);
        double q_face = faces.get_charge(face);

        double r_cut2 = faces.get_area(face) * 50.0;
        double smooth_factor = sqrt(r_cut2) / 10.0;

        // Find weights and normalization factor for all the atoms for given face
        weights = vector<double>(n_atoms);
        double w_sum = 0.0;
        for (int atom = 0; atom < n_atoms; ++atom) {
            double dist2 = point1.periodic_distance2(get_point(atom), sizes.xbox, sizes.ybox);
            if (dist2 > r_cut2) continue;

            double w = exp(-1.0 * sqrt(dist2) / smooth_factor);
            weights[atom] = w;
            w_sum += w;
        }

        // Store the partial charges on atoms
        w_sum = 1.0 / w_sum;
        for (int atom = 0; atom < n_atoms; ++atom)
            if (weights[atom] > 0)
                charges[atom] += weights[atom] * w_sum * q_face;
    }

    for (int atom = 0; atom < n_atoms; ++atom) {
        Vec3 force = fields.get_elfield(atom) * charges[atom];
        interpolation.push_back(Solution(force, charges[atom]));
    }
}

Vec3 ForceReader::get_force(const int i) const {
    require(i >= 0 && i < size(), "Invalid index: " + to_string(i));
    return interpolation[i].vector;
}

double ForceReader::get_charge(const int i) const {
    require(i >= 0 && i < size(), "Invalid index: " + to_string(i));
    return interpolation[i].scalar;
}

void ForceReader::print_statistics(const Medium::Sizes& sizes, const double E0) {
    if (!MODES.VERBOSE) return;

    double Q = E0 * sizes.xbox * sizes.ybox;
    double q = 0;
    for (int i = 0; i < size(); ++i)
        q += get_charge(i);

    cout << "  Q / sum(q) = " << Q << " / " << q << " = " << Q/q << endl;
}

} // namespace femocs
