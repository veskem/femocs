/*
 * SolutionReader.cpp
 *
 *  Created on: 6.6.2016
 *      Author: veske
 */

#include "SolutionReader.h"
#include "Macros.h"
#include "getelec.h"
#include "Config.h"

#include <float.h>
#include <stdio.h>

using namespace std;
namespace femocs {

/* ==========================================
 * ============= SOLUTION READER ============
 * ========================================== */

// Initialize SolutionReader
SolutionReader::SolutionReader() :
        vec_label("vec"), vec_norm_label("vec_norm"), scalar_label("scalar"),
        limit_min(0), limit_max(0), sort_atoms(false), dim(0), rank(0), interpolator(NULL)
{
    reserve(0);
}

SolutionReader::SolutionReader(Interpolator* i, const string& vec_lab, const string& vec_norm_lab, const string& scal_lab) :
        vec_label(vec_lab), vec_norm_label(vec_norm_lab), scalar_label(scal_lab),
        limit_min(0), limit_max(0), sort_atoms(false), dim(0), rank(0), interpolator(i)
{
    reserve(0);
}

void SolutionReader::calc_interpolation() {
    require(interpolator, "NULL interpolator cannot be used!");

    const int n_atoms = size();
    if (interpolation.size() != 0) {
        interpolation.clear();
        interpolation.reserve(n_atoms);
    }

    // Sort atoms into sequential order to speed up interpolation
    if (sort_atoms) sort_spatial();
    array<double,8> sf;

    int cell = 0;
    for (int i = 0; i < n_atoms; ++i) {
        Point3 point = get_point(i);

        // Depending on interpolation dimension and rank, pick corresponding functions
        if (dim == 2) {
            if (rank == 1) {
                cell = interpolator->lintri.locate_cell(point, abs(cell));
                append_interpolation(interpolator->lintri.interp_solution(point, cell));
            } else if (rank == 2) {
                cell = interpolator->quadtri.locate_cell(point, abs(cell));
                append_interpolation(interpolator->quadtri.interp_solution(point, cell));
            } else if (rank == 3) {
                cell = interpolator->linquad.locate_cell(point, abs(cell));
                append_interpolation(interpolator->linquad.interp_solution(point, cell));
            }
        } else {
            if (rank == 1) {
                cell = interpolator->lintet.locate_cell(point, abs(cell));
                append_interpolation(interpolator->lintet.interp_solution(point, cell));
            } else if (rank == 2) {
                cell = interpolator->quadtet.locate_cell(point, abs(cell));
                append_interpolation(interpolator->quadtet.interp_solution(point, cell));
            } else if (rank == 3) {
                cell = interpolator->linhex.locate_cell(point, abs(cell));
                append_interpolation(interpolator->linhex.interp_solution(point, cell));
            }
        }

        set_marker(i, cell);
    }

    // Sort atoms back to their initial order
    if (sort_atoms) {
        for (int i = 0; i < n_atoms; ++i)
            interpolation[i].id = atoms[i].id;
        sort( interpolation.begin(), interpolation.end(), Solution::sort_up() );
        sort( atoms.begin(), atoms.end(), Atom::sort_id() );
    }
}

void SolutionReader::calc_interpolation(vector<int>& atom2cell) {
    require(interpolator, "NULL interpolator cannot be used!");

    const int n_atoms = size();
    const bool cells_known = atom2cell.size() == n_atoms;

    // are the atoms already mapped against the triangles?
    if (cells_known) {
        // ...yes, no need to calculate them again, just interpolate
        for (int i = 0; i < n_atoms; ++i) {
            // locate the face
            int cell = atom2cell[i];
            set_marker(i, cell);

            // calculate the interpolation
            if (dim == 2) {
                if (rank == 1)
                    append_interpolation(interpolator->lintri.interp_solution(get_point(i), cell));
                else if (rank == 2)
                    append_interpolation(interpolator->quadtri.interp_solution(get_point(i), cell));
                else if (rank == 3)
                    append_interpolation(interpolator->linquad.interp_solution(get_point(i), cell));
            } else {
                if (rank == 1)
                    append_interpolation(interpolator->lintet.interp_solution(get_point(i), cell));
                else if (rank == 2)
                    append_interpolation(interpolator->quadtet.interp_solution(get_point(i), cell));
                else if (rank == 3)
                    append_interpolation(interpolator->linhex.interp_solution(get_point(i), cell));
            }
        }

    // ...nop, do it and interpolate
    } else {
        calc_interpolation();
        atom2cell = vector<int>(n_atoms);
        for (int i = 0; i < n_atoms; ++i)
            atom2cell[i] = abs(get_marker(i));
    }
}

// Reserve memory for solution vectors
void SolutionReader::reserve(const int n_nodes) {
    require(n_nodes >= 0, "Invalid number of nodes: " + to_string(n_nodes));

    atoms.clear();
    atoms.reserve(n_nodes);
    interpolation.clear();
    interpolation.reserve(n_nodes);
}

// Append solution
void SolutionReader::append_interpolation(const Solution& s) {
    expect(interpolation.size() < interpolation.capacity(), "Allocated vector size exceeded!");
    interpolation.push_back(s);
}

// Get pointer to interpolation vector
vector<Solution>* SolutionReader::get_interpolations() {
    return &interpolation;
}

// Get i-th Solution
Solution SolutionReader::get_interpolation(const int i) const {
    require(i >= 0 && i < static_cast<int>(interpolation.size()), "Index out of bounds: " + to_string(i));
    return interpolation[i];
}

// Set i-th Solution
void SolutionReader::set_interpolation(const int i, const Solution& s) {
    require(i >= 0 && i < static_cast<int>(interpolation.size()), "Index out of bounds: " + to_string(i));
    interpolation[i] = s;
}

// Compile data string from the data vectors for file output
string SolutionReader::get_data_string(const int i) const{
    if (i < 0) return "SolutionReader properties=id:I:1:pos:R:3:marker:I:1:force:R:3:" + vec_norm_label + ":R:1:" + scalar_label + ":R:1";

    ostringstream strs; strs << fixed;
    strs << atoms[i] << " " << interpolation[i];
    return strs.str();
}

// Get information about data vectors for .vtk file
void SolutionReader::get_point_data(ofstream& out) const {
    const int n_atoms = size();

    // write IDs of atoms
    out << "SCALARS id int\nLOOKUP_TABLE default\n";
    for (int i = 0; i < n_atoms; ++i)
        out << atoms[i].id << "\n";

    // write atom markers
    out << "SCALARS marker int\nLOOKUP_TABLE default\n";
    for (int i = 0; i < n_atoms; ++i)
        out << atoms[i].marker << "\n";

    // output scalar (electric potential, temperature etc)
    out << "SCALARS " << scalar_label << " double\nLOOKUP_TABLE default\n";
    for (int i = 0; i < n_atoms; ++i)
        out << interpolation[i].scalar << "\n";

    // output vector magnitude explicitly to make it possible to apply filters in ParaView
    out << "SCALARS " << vec_norm_label << " double\nLOOKUP_TABLE default\n";
    for (int i = 0; i < n_atoms; ++i)
        out << interpolation[i].norm << "\n";

    // output vector data (electric field, current density etc)
    out << "VECTORS " << vec_label << " double\n";
    for (int i = 0; i < n_atoms; ++i)
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
            if (dist2 > r_cut2) continue;

            double w = exp(-1.0 * sqrt(dist2) / smooth_factor);
            w_sum += w;
            elfield += interpolation[i].vector * w;
            potential += interpolation[i].scalar * w;
        }

    if (w_sum > 0) {
        elfield *= (1.0 / w_sum); potential /= w_sum;
        return Solution(elfield, potential);
    }

    expect(false, "Node " + to_string(I) + " can't be averaged!");
    return(interpolation[I]);
}

// Get histogram for electric field x,y,z component or for its norm
void SolutionReader::get_histogram(vector<int> &bins, vector<double> &bounds, const int coordinate) {
    require(coordinate >= 0 && coordinate <= 4, "Invalid component: " + to_string(coordinate));

    const int n_atoms = size();
    const int n_bins = bins.size();
    const int n_bounds = bounds.size();

    // Find minimum and maximum values from all non-error values
    double value_min = DBL_MAX;
    double value_max =-DBL_MAX;
    double value;
    for (int i = 0; i < n_atoms; ++i) {
        if (coordinate == 4) value = interpolation[i].scalar;
        else if (coordinate == 3) value = interpolation[i].norm;
        else                 value = interpolation[i].vector[coordinate];

        value_min = min(value_min, value);
        value_max = max(value_max, value);
    }

    // Fill the bounds with values value_min:value_step:(value_max + epsilon)
    // Epsilon is added to value_max to include the maximum value in the up-most bin
    double value_step = (value_max - value_min) / n_bins;
    for (int i = 0; i < n_bounds; ++i)
        bounds[i] = value_min + value_step * i;
    bounds[n_bounds-1] += 1e-5 * value_step;

    for (int i = 0; i < n_atoms; ++i)
        for (int j = 0; j < n_bins; ++j) {
            if (coordinate == 4) value = interpolation[i].scalar;
            else if (coordinate == 3) value = interpolation[i].norm;
            else                 value = interpolation[i].vector[coordinate];

            if (value >= bounds[j] && value < bounds[j+1]) {
                bins[j]++;
                break;
            }
        }
}

// Clean the interpolation from peaks using histogram cleaner
void SolutionReader::histogram_clean(const int coordinate, const double r_cut) {
    require(coordinate >= 0 && coordinate <= 4, "Invalid coordinate: " + to_string(coordinate));
    const int n_atoms = size();
    const int n_bins = static_cast<int>(n_atoms / 250);

    if (n_bins <= 1 || r_cut < 0.1) return;

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
        if (coordinate < 3)  value = fabs(interpolation[i].vector[coordinate]);
        else if (coordinate == 3) value = fabs(interpolation[i].norm);
        else if (coordinate == 4) value = fabs(interpolation[i].scalar);

        if (value < value_min || value > value_max)
            interpolation[i] = get_average_solution(i, r_cut);
    }
}

bool SolutionReader::clean(const double r_cut, const bool use_hist_clean) {
    const int n_atoms = size();

    // Apply histogram cleaner for the solution
    if (use_hist_clean) {
        histogram_clean(0, r_cut);  // clean by vector x-component
        histogram_clean(1, r_cut);  // clean by vector y-component
        histogram_clean(2, r_cut);  // clean by vector z-component
        histogram_clean(3, r_cut);  // clean by vector norm
        histogram_clean(4, r_cut);  // clean by scalar
    }

    // replace the NaN-s with average solution
    bool fail = false;
    for (int i = 0; i < n_atoms; ++i) {
        double s = interpolation[i].scalar;
        if (s != s) {
            expect(false, "Replacing NaN at interpolation point " + to_string(i));
            interpolation[i] = get_average_solution(i, r_cut);
            fail = true;
        }
    }
    return fail;
}

// Initialise statistics about the solution
void SolutionReader::init_statistics() {
    Medium::init_statistics();
    stat.vec_norm_min = stat.scal_min = DBL_MAX;
    stat.vec_norm_max = stat.scal_max = -DBL_MAX;
}

// Calculate statistics about the solution
void SolutionReader::calc_statistics() {
    init_statistics();
    Medium::calc_statistics();

    for (int i = 0; i < size(); ++i) {
        double norm = interpolation[i].norm;
        double scalar = interpolation[i].scalar;
        stat.vec_norm_max = max(stat.vec_norm_max, norm);
        stat.vec_norm_min = min(stat.vec_norm_min, norm);
        stat.scal_max = max(stat.scal_max, scalar);
        stat.scal_min = min(stat.scal_min, scalar);
    }
}

// Print statistics about interpolated solution
void SolutionReader::print_statistics() {
    if (!MODES.VERBOSE) return;

    const int n_atoms = size();
    Vec3 vec(0), rms_vec(0);
    double scalar = 0, rms_scalar = 0;

    for (int i = 0; i < n_atoms; ++i) {
        Vec3 v = interpolation[i].vector;
        double s = interpolation[i].scalar;

        vec += v; rms_vec += v * v;
        scalar += s; rms_scalar += s * s;
    }

    vec *= (1.0 / n_atoms);
    rms_vec = Vec3(sqrt(rms_vec.x), sqrt(rms_vec.y), sqrt(rms_vec.z)) * (1.0 / n_atoms);
    scalar = scalar / n_atoms;
    rms_scalar = sqrt(rms_scalar) / n_atoms;

    stringstream stream;
    stream << "mean " << vec_label << ": \t" << vec;
    stream << "\n   rms " << vec_label << ": \t" << rms_vec;
    stream << "\n  mean & rms " << scalar_label << ": " << scalar << "\t" << rms_scalar;
    write_verbose_msg(stream.str());
}

int SolutionReader::export_results(const int n_points, const string &data_type, const bool append, double* data) {
    check_return(size() == 0, "No " + data_type + " to export!");

    // Check whether the already excising data must survive or not
    if (!append) {
        // No, clear the previous data
        int export_type_len = 1;
        if (data_type == vec_label)
            export_type_len = 3;
        for (int i = 0; i < export_type_len * n_points; ++i)
            data[i] = 0;
    }

    // Pass the desired solution
    for (int i = 0; i < size(); ++i) {
        int id = get_id(i);
        if (id < 0 || id >= n_points) continue;

        if (data_type == vec_label) {
            id *= 3;
            for (double v : interpolation[i].vector)
                data[id++] = v;
        }
        else if (data_type == scalar_label)
            data[id] = interpolation[i].scalar;
        else
            data[id] = interpolation[i].norm;
    }

    return 0;
}

int SolutionReader::interpolate_results(const int n_points, const string &data_type, const double* x,
        const double* y, const double* z, double* data) {
    check_return(size() == 0, "No " + data_type + " to interpolate!");

    // transfer coordinates
    SolutionReader sr(interpolator, vec_label, vec_norm_label, scalar_label);
    sr.set_preferences(sort_atoms, dim, rank);
    sr.reserve(n_points);
    for (int i = 0; i < n_points; ++i)
        sr.append(Atom(i, Point3(x[i], y[i], z[i]), 0));

    // interpolate solution
    sr.calc_interpolation();

    // export solution
    if (data_type == vec_label) {
        for (int i = 0; i < n_points; ++i) {
            int j = 3*i;
            for (double v : sr.interpolation[i].vector)
                data[j++] = v;
        }
    } else if (data_type == scalar_label) {
        for (int i = 0; i < n_points; ++i)
            data[i] = sr.interpolation[i].scalar;
    } else {
        for (int i = 0; i < n_points; ++i)
            data[i] = sr.interpolation[i].norm;
    }

    return 0;
}

/* ==========================================
 * ============== FIELD READER ==============
 * ========================================== */

FieldReader::FieldReader(Interpolator* i) :
        SolutionReader(i, "elfield", "elfield_norm", "potential"),
        E0(0), radius1(0), radius2(0) {}

void FieldReader::compare_interpolators(fch::PoissonSolver<3> &poisson, const Medium::Sizes &sizes) {
    const double x = sizes.xmid;
    const double y = sizes.ymid;
    const double zmin = 1.0 + sizes.zmax;
    const double step = 0.001;
    const int n_points = 10000;

    
    reserve(n_points);
    for (int i = 0; i < n_points; ++i)
        append(Point3(x + 0.1*i*step, y + 0.1*i*step, zmin + i * step));

    int cell_index;
    array<double,8> shape_functions;

    double t0;

    start_msg(t0, "poisson");
    cell_index = 0;
    for (int i = 0; i < n_points; ++i) {
        Point3 p = get_point(i);
        dealii::Point<3> deal_point(p.x, p.y, p.z);

        cell_index = interpolator->linhex.locate_cell(p, cell_index);
        int cell_index_in_deal = interpolator->linhex.femocs2deal(cell_index);
        if (cell_index_in_deal < 0)
            append_interpolation(Solution(0));
        else {
            double val1 = poisson.probe_efield_norm(deal_point,cell_index_in_deal);
            double val2 = poisson.probe_potential(deal_point, cell_index_in_deal);
            append_interpolation(Solution(Vec3(0), val1, val2));
        }
    }
    end_msg(t0);

    write("out/potential1.xyz");

    interpolation.clear();
    interpolation.reserve(n_points);

    start_msg(t0, "linhexs");
    cell_index = 0;
    for (int i = 0; i < n_points; ++i) {
        Point3 p = get_point(i);
        cell_index = interpolator->linhex.locate_cell(p, cell_index);
        append_interpolation(interpolator->linhex.interp_solution_v2(p, cell_index));
    }
    end_msg(t0);
    write("out/potential2.xyz");

    interpolation.clear();
    interpolation.reserve(n_points);

    start_msg(t0, "lintets");
    cell_index = 0;
    for (int i = 0; i < n_points; ++i) {
        Point3 p = get_point(i);
        cell_index = interpolator->lintet.locate_cell(p, cell_index);
        append_interpolation(interpolator->lintet.interp_solution_v2(p, cell_index));
    }
    end_msg(t0);
    write("out/potential3.xyz");

    interpolation.clear();
    interpolation.reserve(n_points);

    start_msg(t0, "quadtets");
    cell_index = 0;
    for (int i = 0; i < n_points; ++i) {
        Point3 p = get_point(i);
        cell_index = interpolator->quadtet.locate_cell(p, cell_index);
        append_interpolation(interpolator->quadtet.interp_solution_v2(p, cell_index));
    }
    end_msg(t0);
    write("out/potential4.xyz");
}

void FieldReader::test_corners(const TetgenMesh& mesh) const {
    int cell_index;
    array<double,8> shape_functions;
    array<double,4> bcc;

    int cntr = 0;
    for (int i = 0; i < mesh.tets.size(); ++i)
        if (mesh.tets.get_marker(i) == TYPES.VACUUM && cntr++ >= 0) {
            cell_index = i;
            break;
        }

    SimpleElement stet = mesh.tets[cell_index];

    Point3 n1 = mesh.nodes[stet[0]];
    Point3 n2 = mesh.nodes[stet[1]];
    Point3 n3 = mesh.nodes[stet[2]];
    Point3 n4 = mesh.nodes[stet[3]];

    vector<Point3> points;

    points.push_back(n1);
    points.push_back(n2);
    points.push_back(n3);
    points.push_back(n4);
    points.push_back((n1 + n2) / 2.0);
    points.push_back((n1 + n3) / 2.0);
    points.push_back((n1 + n4) / 2.0);
    points.push_back((n2 + n3) / 2.0);
    points.push_back((n2 + n4) / 2.0);
    points.push_back((n3 + n4) / 2.0);
    points.push_back((n1 + n2 + n3) / 3.0);
    points.push_back((n1 + n2 + n4) / 3.0);
    points.push_back((n1 + n3 + n4) / 3.0);
    points.push_back((n2 + n3 + n4) / 3.0);
    points.push_back((n1 + n2 + n3 + n4) / 4.0);

    vector<string> labels = {"n1","n2","n3","n4","c12","c13","c14","c23","c24","c34","c123","c124","c134","c234","c1234"};
    require(labels.size() == points.size(), "Incompatible vectors!");
    
    cout << "results for tet=" << cell_index << ", 4tet=" << 4*cell_index << endl;

    for (int i = 0; i < 4; ++i) {
        cell_index = interpolator->linhex.locate_cell(points[i], 0);
        interpolator->lintet.get_shape_functions(bcc, points[i], 0);
        cout << endl << labels[i] << ":\t" << cell_index << "\t";
        for (double b : bcc)
            cout << ", " << b;
    }

    for (int i = 4; i < labels.size(); ++i) {
        cell_index = interpolator->linhex.locate_cell(points[i], 0);
        interpolator->lintet.get_shape_functions(bcc, points[i], abs(int(cell_index/4)));
        cout << endl << labels[i] << ":\t" << cell_index << "\t";
        for (double b : bcc)
            cout << ", " << b;
    }
}

// Interpolate electric field and potential on a set of points
void FieldReader::interpolate(const int n_points, const double* x, const double* y, const double* z) {
    // store the point coordinates
    reserve(n_points);
    for (int i = 0; i < n_points; ++i)
        append(Atom(i, Point3(x[i], y[i], z[i]), 0));

    // interpolate solution
    calc_interpolation();
}

// Interpolate electric field and potential on a Medium atoms
void FieldReader::interpolate(const Medium &medium) {
    const int n_atoms = medium.size();

    // store the atom coordinates
    reserve(n_atoms);
    for (int i = 0; i < n_atoms; ++i)
        append( Atom(i, medium.get_point(i), 0) );

    // interpolate solution
    calc_interpolation();

    // store the original atom id-s
    for (int i = 0; i < n_atoms; ++i)
        atoms[i].id = medium.get_id(i);
}

void FieldReader::interpolate(const fch::DealSolver<3>& solver) {
    // import the surface nodes the solver needs
    vector<dealii::Point<3>> nodes;
    solver.get_surface_nodes(nodes);

    const int n_atoms = nodes.size();

    // store the node coordinates
    reserve(n_atoms);
    int i = 0;
    for (dealii::Point<3>& node : nodes)
        append( Atom(i++, Point3(node[0], node[1], node[2]), 0) );

    // interpolate or assign solution on the atoms
    calc_interpolation();
}

// Return electric field in i-th interpolation point
Vec3 FieldReader::get_elfield(const int i) const {
    require(i >= 0 && i < size(), "Invalid index: " + to_string(i));
    return interpolation[i].vector;
}

// Return electric field norm in i-th interpolation point
double FieldReader::get_elfield_norm(const int i) const {
    require(i >= 0 && i < size(), "Invalid index: " + to_string(i));
    return interpolation[i].norm;
}

// Return electric potential in i-th interpolation point
double FieldReader::get_potential(const int i) const {
    require(i >= 0 && i < size(), "Invalid index: " + to_string(i));
    return interpolation[i].scalar;
}

// Analytical potential for i-th point near the hemisphere
double FieldReader::get_analyt_potential(const int i, const Point3& origin) const {
    require(i >= 0 && i < size(), "Invalid index: " + to_string(i));

    Point3 point = get_point(i);
    point -= origin;
    double r = point.distance(Point3(0));
    return -E0 * point.z * (1 - pow(radius1 / r, 3.0));
}

// Analytical electric field for i-th point near the hemisphere
Vec3 FieldReader::get_analyt_field(const int i, const Point3& origin) const {
    require(i >= 0 && i < size(), "Invalid index: " + to_string(i));

    Point3 point = get_point(i);
    point -= origin;
    double r5 = pow(point.x * point.x + point.y * point.y + point.z * point.z, 2.5);
    double r3 = pow(radius1, 3.0);
    double f = point.x * point.x + point.y * point.y - 2.0 * point.z * point.z;

    double Ex = 3 * E0 * r3 * point.x * point.z / r5;
    double Ey = 3 * E0 * r3 * point.y * point.z / r5;
    double Ez = E0 * (1.0 - r3 * f / r5);

    return Vec3(Ex, Ey, Ez);
}

// Analytical field enhancement for ellipsoidal nanotip
double FieldReader::get_analyt_enhancement() const {
    expect(radius1 > 0, "Invalid nanotip minor semi-axis: " + to_string(radius1));

    if ( radius2 <= radius1 )
        return 3.0;
    else {
        double nu = radius2 / radius1;
        double zeta = sqrt(nu*nu - 1);
        return pow(zeta, 3.0) / (nu * log(zeta + nu) - zeta);
    }
}

// Compare the analytical and calculated field enhancement
bool FieldReader::check_limits(const vector<Solution>* solutions) const {
    if (limit_min == limit_max)
        return false;

    double Emax = -1e100;
    if (solutions) {
        for (Solution s : *solutions)
            Emax = max(Emax, s.norm);
    } else {
        for (Solution s : interpolation)
            Emax = max(Emax, s.norm);
    }

    const double gamma1 = fabs(Emax / E0);
    const double gamma2 = get_analyt_enhancement();
    const double beta = fabs(gamma1 / gamma2);

    stringstream stream;
    stream << fixed << setprecision(3);
    stream << "field enhancements:  (F)emocs:" << gamma1
            << "  (A)nalyt:" << gamma2
            << "  F-A:" << gamma1 - gamma2
            << "  F/A:" << gamma1 / gamma2;

    write_verbose_msg(stream.str());
    return beta < limit_min || beta > limit_max;
}

// Set parameters for calculating analytical solution
void FieldReader::set_check_params(const double E0, const double limit_min, const double limit_max,
        const double radius1, const double radius2) {
    this->E0 = E0;
    this->limit_min = limit_min;
    this->limit_max = limit_max;
    this->radius1 = radius1;
    if (radius2 > radius1)
        this->radius2 = radius2;
    else
        this->radius2 = radius1;
}

// Export interpolated electric field
int FieldReader::export_elfield(const int n_atoms, double* Ex, double* Ey, double* Ez, double* Enorm) {
    if (n_atoms <= 0) return 0;
    check_return(size() == 0, "No field to export!");

    // Initially pass the zero electric field for all the atoms
    for (int i = 0; i < n_atoms; ++i) {
        Ex[i] = 0;
        Ey[i] = 0;
        Ez[i] = 0;
        Enorm[i] = 0;
    }

    // Pass the the calculated electric field for stored atoms
    for (int i = 0; i < size(); ++i) {
        int id = get_id(i);
        if (id < 0 || id >= n_atoms) continue;

        Ex[id] = interpolation[i].vector.x;
        Ey[id] = interpolation[i].vector.y;
        Ez[id] = interpolation[i].vector.z;
        Enorm[id] = interpolation[i].norm;
    }
    return 0;
}

int FieldReader::interpolate_surface_elfield(const int n_points, const double* x, const double* y, const double* z,
        double* Ex, double* Ey, double* Ez, double* Enorm, int* flag) {
    if (n_points <= 0) return 0;
    check_return(interpolator->nodes.size() == 0, "No surface field to export!");

    FieldReader fr(interpolator);
    fr.set_preferences(false, 2, rank);
    fr.interpolate(n_points, x, y, z);

    for (int i = 0; i < n_points; ++i) {
        Ex[i] = fr.interpolation[i].vector.x;
        Ey[i] = fr.interpolation[i].vector.y;
        Ez[i] = fr.interpolation[i].vector.z;
        Enorm[i] = fr.interpolation[i].norm;
        flag[i] = fr.atoms[i].marker < 0;
    }

    fr.write("out/interpolation_E_surf.movie");
    return check_limits(fr.get_interpolations());
}

int FieldReader::interpolate_elfield(const int n_points, const double* x, const double* y, const double* z,
        double* Ex, double* Ey, double* Ez, double* Enorm, int* flag) {
    if (n_points <= 0) return 0;
    check_return(interpolator->nodes.size() == 0, "No electric field to interpolate!");

    FieldReader fr(interpolator);
    fr.set_preferences(false, 3, rank);
    fr.interpolate(n_points, x, y, z);

    for (int i = 0; i < n_points; ++i) {
        Ex[i] = fr.interpolation[i].vector.x;
        Ey[i] = fr.interpolation[i].vector.y;
        Ez[i] = fr.interpolation[i].vector.z;
        Enorm[i] = fr.interpolation[i].norm;
        flag[i] = fr.atoms[i].marker < 0;
    }

    fr.write("out/interpolation_E.movie");
    return check_limits(fr.get_interpolations());
}

int FieldReader::interpolate_phi(const int n_points, const double* x, const double* y, const double* z,
        double* phi, int* flag) {
    if (n_points <= 0) return 0;
    check_return(interpolator->nodes.size() == 0, "No electric potential to interpolate!");

    FieldReader fr(interpolator);
    fr.set_preferences(false, 3, 2);
    fr.interpolate(n_points, x, y, z);

    for (int i = 0; i < n_points; ++i) {
        phi[i] = fr.interpolation[i].scalar;
        flag[i] = fr.atoms[i].marker < 0;
    }

    return 0;
}

/* ==========================================
 * =============== HEAT READER ==============
 * ========================================== */

HeatReader::HeatReader(Interpolator* i) :
        SolutionReader(i, "rho", "rho_norm", "temperature") {}

// Interpolate solution on Medium atoms
void HeatReader::interpolate(const Medium &medium) {
    const int n_atoms = medium.size();

    // store the atom coordinates
    reserve(n_atoms);
    for (int i = 0; i < n_atoms; ++i)
        append( Atom(i, medium.get_point(i), 0) );

    // interpolate or assign solution
    calc_interpolation();

    // restore the original atom id-s
    for (int i = 0; i < n_atoms; ++i)
        atoms[i].id = medium.get_id(i);
}

void HeatReader::interpolate(fch::DealSolver<3>& solver) {
    // import the surface nodes the solver needs
    vector<dealii::Point<3>> nodes;
    solver.get_surface_nodes(nodes);

    const int n_atoms = nodes.size();

    // store the node coordinates
    reserve(n_atoms);
    int i = 0;
    for (dealii::Point<3>& node : nodes)
        append( Atom(i++, Point3(node[0], node[1], node[2]), 0) );

    // interpolate or assign solution on the atoms
    calc_interpolation();
}

// Export interpolated temperature
int HeatReader::export_temperature(const int n_atoms, double* T) {
    if (n_atoms <= 0) return 0;
    check_return(size() == 0, "No temperature to export!");

    // Pass the the calculated temperature for stored atoms
    for (int i = 0; i < size(); ++i) {
        int identifier = get_id(i);
        if (identifier < 0 || identifier >= n_atoms) continue;
        T[identifier] = get_temperature(i);
    }

    return 0;
}

Vec3 HeatReader::get_rho(const int i) const {
    require(i >= 0 && i < size(), "Invalid index: " + to_string(i));
    return interpolation[i].vector;
}

double HeatReader::get_rho_norm(const int i) const {
    require(i >= 0 && i < size(), "Invalid index: " + to_string(i));
    return interpolation[i].norm;
}

double HeatReader::get_temperature(const int i) const {
    require(i >= 0 && i < size(), "Invalid index: " + to_string(i));
    return interpolation[i].scalar;
}

/* ==========================================
 * ============= EMISSION READER ============
 * ========================================== */

EmissionReader::EmissionReader(const FieldReader *fr, const HeatReader *hr,
        const fch::PoissonSolver<3> *p, Interpolator* i) :
        SolutionReader(i, "none", "rho_norm", "temperature"),
        fields(fr), heat(hr), mesh(NULL), poisson(p)
{}

void EmissionReader::initialize(const TetgenMesh* m) {
    mesh = m;

    int n_nodes = fields->size();
    require(n_nodes > 0, "EmissionReader can't use empty dependencies!");

    atoms = fields->atoms;

    // deallocate and allocate currents data
    current_densities.resize(n_nodes);
    nottingham.resize(n_nodes);
    currents.resize(n_nodes);
    field_loc.resize(n_nodes);

    //deallocate and allocate lines
    rline.resize(n_lines);
    Vline.resize(n_lines);

    // find Fmax
    get_field_loc();

    //Initialise data
    global_data.Jmax = 0.;
    global_data.Frep = global_data.Fmax;
    global_data.Jrep = 0.;
    global_data.multiplier = 1.;
}

void EmissionReader::emission_line(const Point3& point, const Vec3& direction, const double rmax) {
    const int interpolation_rank = 3;
    const double nm_per_angstrom = 0.1;
    const double rmin = 1.e-5 * rmax;
    Point3 pfield(direction.x, direction.y, direction.z);

    FieldReader fr(interpolator);
    fr.set_preferences(false, 3, interpolation_rank);
    fr.reserve(n_lines);

    for (int i = 0; i < n_lines; i++){
        rline[i] = rmin + ((rmax - rmin) * i) / (n_lines - 1);
        fr.append(point - pfield * rline[i]);
    }
    fr.calc_interpolation();


    for (int i = 0; i < n_lines; i++){
        int hex = fr.get_marker(i);
        int hex_deal = interpolator->linhex.femocs2deal(hex);

        Point3 p = fr.get_point(i);
        dealii::Point<3> p_deal(p.x, p.y, p.z);

        Vline[i] = global_data.multiplier * poisson->probe_potential(p_deal, hex_deal);
        rline[i] *= nm_per_angstrom;

    }

    for (int i = 0; i < n_lines; i++){
        Vline[i] -= Vline[0];
        rline[i] -= rline[0];
    }

    // Check data condition (should be monotonous)
    for (int i = 1; i < n_lines; ++i) { // go through points
        if (Vline[i] < Vline[i-1]) { // if decreasing at a point
            double dVdx = 0.0;
            int j;
            for (j = i + 1; j < n_lines; ++j) {
                if (Vline[j] > Vline[i-1]) {
                    dVdx = (Vline[j] - Vline[i-1]) / (rline[j] - rline[i-1]);
                    break;
                }
            }

            if (dVdx == 0.0) {
                if (i > 1)
                    dVdx = (Vline[i-1] - Vline[i-2]) / (rline[i-1] - rline[i-2]);
                else
                    write_verbose_msg("Non-monotonous Vline could not be recovered at i = "
                            + to_string(i));
            }
            for (int k = 0; k <= j; ++k)
                Vline[k] =  Vline[i-1] + (rline[k] - rline[i-1]) * dVdx;
        }
    }
}

void EmissionReader::calc_representative() {
    double area = 0.; // total emitting (FWHM) area
    double FJ = 0.; // int_FWHMarea (F*J)dS
    global_data.I_tot = 0;
    global_data.I_fwhm = 0;

    for (int i = 0; i < currents.size(); ++i){ // go through face centroids
        int tri = mesh->quads.to_tri(abs(fields->get_marker(i)));
        // quadrangle area is 1/3 of corresponding triangle area
        double face_area = mesh->tris.get_area(tri) / 3.;
        currents[i] = face_area * current_densities[i];
        global_data.I_tot += currents[i];

        if (current_densities[i] > global_data.Jmax * 0.5){ //if point eligible
            area += face_area; // increase total area
            global_data.I_fwhm += currents[i]; // increase total current
            FJ += currents[i] * field_loc[i].norm();
        }
    }

    global_data.Jrep = global_data.I_fwhm / area;
    global_data.Frep = global_data.multiplier * FJ / global_data.I_fwhm;
}

void EmissionReader::inject_electrons(double delta_t, double Wsp, vector<Point3> &pos, vector<int> &cells) {

    for (int i = 0; i < fields->size(); ++i){ // go through face centroids
        double current = currents[i] * electrons_per_fs;
        double charge = current * delta_t; //in e
        double n_sps = charge / Wsp;

        int intpart = (int) floor(n_sps);
        double frpart = n_sps - intpart;
        int n_electrons_sp = intpart;

        if ((double) rand() / RAND_MAX < frpart)
            n_electrons_sp++;

        if (n_electrons_sp == 0) continue;

        //TODO : fix rng

        int quad = abs(fields->get_marker(i));
        int tri = mesh->quads.to_tri(quad);
        int hex = mesh->quad2hex(quad, TYPES.VACUUM);
        hex = interpolator->linhex.femocs2deal(hex);
        SimpleQuad squad = mesh->quads[quad];

        // generate desired amount of electrons
        // that are uniformly distributed on a given quadrangle
        for (int j = 0; j < n_electrons_sp; j++){
            Point3 position = interpolator->linquad.get_rnd_point(quad);
            // push point little bit inside the vacuum mesh
            position += mesh->tris.get_norm(tri) * (mesh->tris.stat.edgemin * 1.e-5);

            pos.push_back(position);
			cells.push_back(hex);
        }
    }
}

void EmissionReader::get_field_loc(){

    global_data.Fmax = 0;

    for (int i = 0; i < fields->size(); ++i) { // go through all face centroids
        int quad = fields->get_marker(i);
        int hex_femocs = mesh->quad2hex(quad, TYPES.VACUUM);
        int hex_deal = interpolator->linhex.femocs2deal(hex_femocs);

        Point3 centroid = fields->get_point(i);
        dealii::Point<3> point(centroid.x, centroid.y, centroid.z);
        dealii::Tensor<1, 3, double> Fdealii = poisson->probe_efield(point, hex_deal);
        field_loc[i] = Vec3(Fdealii);
        global_data.Fmax = max(global_data.Fmax, field_loc[i].norm());
    }

}

void EmissionReader::emission_cycle(double workfunction, bool blunt, bool cold) {

    struct emission gt;
    gt.W = workfunction;    // set workfuntion, must be set in conf. script
    gt.R = 1000.0;   // radius of curvature (overrided by femocs potential distribution)
    gt.gamma = 10;  // enhancement factor (overrided by femocs potential distribution)
    double F, J;    // Local field and current density in femocs units (Angstrom)

    get_field_loc();

    for (int i = 0; i < fields->size(); ++i) { // go through all face centroids

        Vec3 field = field_loc[i];

        F = global_data.multiplier * field.norm();
        gt.mode = 0;
        gt.F = angstrom_per_nm * F;
        gt.Temp = heat->get_temperature(i);
        set_marker(i, 0); // set marker for output emission xyz file. Means No full calculation

        if (F > 0.6 * global_data.Fmax && !blunt){ // Full calculation with line only for high field points
            field.normalize(); // get line direction

            emission_line(get_point(i), field, 1.6 * workfunction / F); //get emission line data

            gt.Nr = n_lines;
            gt.xr = &rline[0];
            gt.Vr = &Vline[0];
            gt.mode = -21; // set mode to potential input data
            set_marker(i, 1); //marker = 1, emission calculated with line
        }
        gt.approx = 0; // simple GTF approximation
        cur_dens_c(&gt); // calculate emission
        if (gt.ierr != 0 )
            write_verbose_msg("GETELEC 1st call returned with error, ierr = " + to_string(gt.ierr));
        J = gt.Jem * nm2_per_angstrom2; // current density in femocs units

        if (J > 0.1 * global_data.Jmax && !cold){ // If J is worth it, calculate with full energy integration
            gt.approx = 1;
            cur_dens_c(&gt);
            if (gt.ierr != 0 )
                write_verbose_msg("GETELEC 2nd call returned with error, ierr = " + to_string(gt.ierr));
            J = gt.Jem * nm2_per_angstrom2;
            set_marker(i, 2);
        }

        global_data.Jmax = max(global_data.Jmax, J); // output data
        current_densities[i] = J;
        nottingham[i] = nm2_per_angstrom2 * gt.heat;
    }
}

double EmissionReader::calc_emission(const double multiplier, const Config::Emission &conf,
        double Vappl) {

    global_data.multiplier = multiplier;
    calc_emission(conf, Vappl);
    return global_data.multiplier;
}

void EmissionReader::calc_emission(const Config::Emission &conf, double Vappl) {

    double theta_old = global_data.multiplier;
    double err_fact = 0.5, error;
    double Fmax_0 = global_data.Fmax;

    for (int i = 0; i < 20; ++i){ // SC calculation loop
        global_data.Jmax = 0.;
        global_data.Fmax = global_data.multiplier * Fmax_0;

        emission_cycle(conf.work_function, conf.blunt, conf.cold);
        calc_representative();

        if (conf.omega_SC <= 0.) break; // if Vappl<=0, SC is ignored
        if (i > 5) err_fact *= 0.5; // if not converged in first 6 steps, reduce factor

        // calculate SC multiplier (function coming from getelec)
        global_data.multiplier = theta_SC(global_data.Jrep / nm2_per_angstrom2,
                conf.omega_SC * Vappl, angstrom_per_nm * global_data.Frep);
        error = global_data.multiplier - theta_old;
        global_data.multiplier = theta_old + error * err_fact;
        theta_old = global_data.multiplier;

        printf("SC cycle #%d, theta=%f, Jrep=%e, Frep=%e, Itot=%e\n", i, global_data.multiplier,
                global_data.Jrep, global_data.Frep, global_data.I_tot);

        // if converged break
        if (abs(error) < conf.SC_error) break;
    }
}

string EmissionReader::get_data_string(const int i) const {
    if (i < 0) return "EmissionReader properties=id:I:1:pos:R:3:marker:I:1:force:R:3:" + vec_norm_label + ":R:1:" + scalar_label + ":R:1";

    ostringstream strs; strs << fixed;
    strs << atoms[i] << ' ' << field_loc[i]
             << ' ' << log(current_densities[i]) << ' ' << log(fabs(nottingham[i]));

    return strs.str();
}

string EmissionReader::get_global_data(const bool first_line) const {
    ostringstream strs;

    // uncomment to add header to the file
//    if (first_line) strs << "time Itot Jrep Frep Jmax Fmax multiplier";

    strs << fixed << setprecision(2)
            << GLOBALS.TIME;
    strs << scientific << setprecision(6)
            << " " << global_data.I_tot << " " << global_data.Jrep << " " << global_data.Frep
            << " " << global_data.Jmax << " " << global_data.Fmax << " " << global_data.multiplier;

    return strs.str();
}

void EmissionReader::export_emission(fch::CurrentHeatSolver<3>& ch_solver) {
    ch_solver.current.set_bc(current_densities);
    ch_solver.heat.set_bc(nottingham);
}

/* ==========================================
 * ============== CHARGE READER =============
 * ========================================== */

ChargeReader::ChargeReader(Interpolator* i) :
        SolutionReader(i, "elfield", "area", "charge"), Q_tot(0) {}

void ChargeReader::calc_charges(const TetgenMesh& mesh, const double E0) {
    const double sign = fabs(E0) / E0;
    const int n_faces = mesh.tris.size();
    const int n_quads_per_triangle = 3;

    // Store the centroids of the triangles
    reserve(n_faces);
    for (int i = 0; i < n_faces; ++i)
        append( Atom(i, mesh.tris.get_centroid(i), 0) );

    // create triangle index to its centroid index mapping
    vector<int> tri2centroid(n_faces);
    for (int face = 0; face < n_faces; ++face) {
        for (int node : mesh.quads[n_quads_per_triangle * face])
            if (mesh.nodes.get_marker(node) == TYPES.FACECENTROID) {
                tri2centroid[face] = node;
                break;
            }
    }

    // Calculate the charges for the triangles
    for (int face = 0; face < n_faces; ++face) {
        double area = mesh.tris.get_area(face);
        Vec3 elfield = interpolator->nodes.get_vector(tri2centroid[face]);
        double charge = eps0 * area * elfield.norm() * sign;
        append_interpolation(Solution(elfield, area, charge));
    }
}

// Remove the atoms and their solutions outside the box
void ChargeReader::clean(const Medium::Sizes& sizes, const double latconst) {
    const int n_atoms = size();
    const int eps = latconst / 2.0;;
    vector<bool> in_box; in_box.reserve(n_atoms);

    // Check the locations of points
    for (int i = 0; i < n_atoms; ++i) {
        Point3 point = get_point(i);
        const bool box_x = point.x >= (sizes.xmin - eps) && point.x <= (sizes.xmax + eps);
        const bool box_y = point.y >= (sizes.ymin - eps) && point.y <= (sizes.ymax + eps);
        const bool box_z = point.z >= (sizes.zmin - eps) && point.z <= (sizes.zmax + eps);
        in_box.push_back(box_x && box_y && box_z);
    }

    const int n_box = vector_sum(in_box);
    vector<Atom> _atoms; _atoms.reserve(n_box);
    vector<Solution> _interpolation; _interpolation.reserve(n_box);

    // Copy the solutions and atoms that remain into box
    for (int i = 0; i < n_atoms; ++i)
        if (in_box[i]) {
            _atoms.push_back(atoms[i]);
            _interpolation.push_back(interpolation[i]);
        }
    atoms = _atoms;
    interpolation = _interpolation;

    // Calculate the new "correct" sum of charges in the remaining region
    Q_tot = 0;
    for (Solution s : interpolation)
        Q_tot += s.scalar;
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

// Check whether charge is conserved within specified limits
bool ChargeReader::check_limits(const vector<Solution>* solutions) const {
    if (limit_min == limit_max)
        return false;

    double q = 0;
    if (solutions) {
        for (Solution s : *solutions)
            q += s.scalar;
    } else {
        for (Solution s : interpolation)
            q += s.scalar;
    }

    stringstream stream;
    stream << fixed << setprecision(3);
    stream << "Q_tot / sum(charge) = " << Q_tot << " / " << q << " = " << Q_tot / q;
    write_verbose_msg(stream.str());

    q = Q_tot / q;
    return q < limit_min || q > limit_max;
}

// Set parameters for calculating analytical solution
void ChargeReader::set_check_params(const double Q_tot, const double limit_min, const double limit_max) {
    this->Q_tot = Q_tot * eps0;
    this->limit_min = limit_min;
    this->limit_max = limit_max;
}

/* ==========================================
 * ============== FORCE READER ==============
 * ========================================== */

ForceReader::ForceReader(Interpolator* i) :
        SolutionReader(i, "force", "pair_potential", "charge") {}

int ForceReader::get_nanotip(Medium& nanotip, const double radius) {
    const int n_atoms = size();
    const double radius2 = radius * radius;
    Medium::calc_statistics();

    // Make map for atoms in nanotip
    Point2 centre(sizes.xmid, sizes.ymid);
    vector<bool> atom_in_nanotip(n_atoms);
    for (int i = 0; i < n_atoms; ++i)
        atom_in_nanotip[i] = centre.distance2(get_point2(i)) <= radius2;

    const int n_nanotip_atoms = vector_sum(atom_in_nanotip);

    // Separate nanotip from substrate
    nanotip.reserve(n_nanotip_atoms);
    for (int i = 0; i < n_atoms; ++i)
        if (atom_in_nanotip[i]) {
            nanotip.append(Atom(i, get_point(i), TYPES.SURFACE));
            set_marker(i, 1);
        } else
            set_marker(i, 0);

    nanotip.calc_statistics();
    return n_nanotip_atoms;
}

void ForceReader::clean_voro_faces(VoronoiMesh& mesh) {
    const int nanotip_end = mesh.nodes.indxs.surf_end;

    // calculate mean area of surface faces
    vector<double> areas(mesh.vfaces.size());
    int n_surface_faces = 0;
    double mean_area = 0;
    for (int cell = 0; cell <= nanotip_end; ++cell)
        for (VoronoiFace face : mesh.voros[cell])
            if (mesh.vfaces.get_marker(face.id) == TYPES.SURFACE) {
                double area = face.area();
                areas[face.id] = area;
                mean_area += area;
                n_surface_faces++;
            }

    mean_area /= n_surface_faces;

    // calculate standard deviation of the areas
    double std = 0;
    for (double area : areas)
        if (area > 0) {
            area -= mean_area;
            std += area * area;
        }
    std = sqrt(std / (n_surface_faces - 1));

    // remove cells whose face area is bigger than threshold
    const double max_area = 5.0 * std + mean_area;

    for (int cell = 0; cell <= nanotip_end; ++cell)
        for (VoronoiFace face1 : mesh.voros[cell])
            if (areas[face1.id] > max_area) {
                for (VoronoiFace face2 : mesh.voros[cell])
                    mesh.vfaces.set_marker(face2.id, TYPES.NONE);
                break;
            }

}

int ForceReader::calc_voronois(VoronoiMesh& mesh, const vector<int>& atom2face,
        const double radius, const double latconst, const string& mesh_quality)
{
    require(interpolator, "NULL interpolator cannot be used!");
    const int n_atoms = size();
    const bool faces_known = n_atoms == atom2face.size();
    // TODO put those values to Config, because they affect heavily how the voronoi charges will look like
    const double max_distance_from_surface = 0.5 * latconst;
    const double shift_distance = 1.0 * latconst;

    Medium nanotip;
    const int n_nanotip_atoms = get_nanotip(nanotip, radius);

    // calculate support points for the nanotip
    // by moving the nanotip points in direction of its corresponding triangle norm by r_cut
    Medium support(n_nanotip_atoms);
    int face = 0;
    for (int i = 0; i < n_atoms; ++i)
        if (get_marker(i)) {
            Point3 point = get_point(i);
            if (faces_known) face = atom2face[i];
            else face = abs( interpolator->lintri.locate_cell(point, face) );

            if (interpolator->lintri.fast_distance(point, face) < max_distance_from_surface) {
                point += interpolator->lintri.get_norm(face) * shift_distance;
                support.append(Atom(i, point, TYPES.VACANCY));
            }
        }

    nanotip += support;
    nanotip.calc_statistics();

    // Generate Voronoi cells around the nanotip
    // r - reconstruct, v - output Voronoi cells, Q - quiet, q - mesh quality
    int err_code = mesh.generate(nanotip, latconst, "rQq" + mesh_quality, "vQ");
    if (err_code) return err_code;

    // Clean the mesh from faces and cells that have node in the infinity
    mesh.clean();

    require(mesh.nodes.size() > 0, "Empty Voronoi mesh cannot be handled!");
    require(mesh.voros.size() > 0, "Empty Voronoi mesh cannot be handled!");

    const int nanotip_end = n_nanotip_atoms - 1;
    const int support_start = n_nanotip_atoms;
    const int support_end = n_nanotip_atoms + support.size() - 1;

    // specify the location of Voronoi faces
    for (int cell = 0; cell <= nanotip_end; ++cell)
        for (VoronoiFace face : mesh.voros[cell]) {
            const int nborcell = face.nborcell(cell);
            if (nborcell < support_start)
                mesh.vfaces.set_marker(face.id, TYPES.PERIMETER);
            else if (nborcell >= support_start && nborcell <= support_end)
                mesh.vfaces.set_marker(face.id, TYPES.SURFACE);
            else
                mesh.vfaces.set_marker(face.id, TYPES.NONE);
        }

    return 0;
}

void ForceReader::calc_charge_and_lorentz(const VoronoiMesh& mesh, const FieldReader& fields) {
    const int n_atoms = size();
    int cell = -1;
    for (int i = 0; i < n_atoms; ++i)
        if (get_marker(i)) {
            Vec3 field = fields.get_elfield(i);
            double charge = 0;
            for (VoronoiFace face : mesh.voros[++cell]) {
                if (mesh.vfaces.get_marker(face.id) == TYPES.SURFACE) {
                    const double dp = field.dotProduct(face.norm(cell));
                    if (dp < 0)
                        charge += dp * face.area();
                }
            }
            charge *= eps0;
            interpolation[i] = Solution(field * charge, 0, charge);
        }
}

// Calculate forces from atomic electric fields and face charges
void ForceReader::calc_lorentz(const FieldReader &fields) {
    const int n_atoms = fields.size();

    // Copy the atom data
    reserve(n_atoms);
    atoms = fields.atoms;

    // Calculate the charges by ensuring the that the total sum of it remains conserved
    vector<double> charges;
    interpolator->lintri.interp_conserved(charges, atoms);

    // calculate forces and store them
    for (int i = 0; i < n_atoms; ++i) {
        Vec3 force = fields.get_elfield(i) * (charges[i] * force_factor);   // [e*V/A]
        interpolation.push_back(Solution(force, 0, charges[i]));
    }
}

void ForceReader::calc_coulomb(const double r_cut) {
    const double r_cut2 = r_cut * r_cut;
    const int n_atoms = size();

    // it is more efficinet to calculate forces from linked list than from neighbour list,
    // as in that ways it is possible to avoid double calculation of the distances
    calc_linked_list(r_cut);
    require(list.size() == n_atoms, "Invalid linked list size: " + to_string(list.size()));
    require(head.size() == nborbox_size[0]*nborbox_size[1]*nborbox_size[2],
            "Invalid linked list header size: " + to_string(head.size()));

    // loop through the atoms
    for (int i = 0; i < n_atoms; ++i) {
        array<int,3> i_atom = nborbox_indices[i];
        Point3 point = atoms[i].point;

        // loop through the boxes where the neighbours are located; there are up to 3^3=27 boxes
        // no periodicity needed, as the charge on simubox boundary is very small
        for (int iz = i_atom[2]-1; iz <= i_atom[2]+1; ++iz) {
            // some of the iterations are be skipped if the box is on simubox boundary
            if (iz < 0 || iz >= nborbox_size[2]) continue;
            for (int iy = i_atom[1]-1; iy <= i_atom[1]+1; ++iy) {
                if (iy < 0 || iy >= nborbox_size[1]) continue;
                for (int ix = i_atom[0]-1; ix <= i_atom[0]+1; ++ix) {
                    if (ix < 0 || ix >= nborbox_size[0]) continue;

                    // transform volumetric neighbour box index to linear one
                    int i_cell = (iz * nborbox_size[1] + iy) * nborbox_size[0] + ix;
                    require(i_cell >= 0 && i_cell < head.size(), "Invalid neighbouring cell index: " + to_string(i_cell));

                    // get the index of first atom in given neighbouring cell and loop through neighbours
                    int j = head[i_cell];
                    while(j >= 0) {
                        // avoid double looping over atom-pairs
                        if (i < j) {
                            Vec3 displacement = point - get_point(j);
                            const double r_squared = displacement.norm2();
                            if (r_squared <= r_cut2) {
                                double r = sqrt(r_squared);
                                double V = exp(-q_screen * r) * couloumb_constant *
                                        get_charge(i) * get_charge(j) / r;

                                Vec3 force = displacement * (V / r_squared);
                                interpolation[i].vector += force;
                                interpolation[j].vector -= force;
                                interpolation[i].norm += 0.5 * V;
                                interpolation[j].norm += 0.5 * V;
                            }
                        }
                        j = list[j];
                    }
                }
            }
        }
    }
}

void ForceReader::distribute_charges(const FieldReader &fields, const ChargeReader& faces,
        const double r_cut, const double smooth_factor) {

    const int n_atoms = fields.size();
    const int n_faces = faces.size();

    // Copy the atom data
    reserve(n_atoms);
    for (int i = 0; i < n_atoms; ++i)
        append(fields.get_atom(i));

    calc_statistics();

    /* Distribute the charges on surface faces between surface atoms.
     * If q_i and Q_j is the charge on i-th atom and j-th face, respectively, then
     *     q_i = sum_j(w_ij * Q_j),  sum_i(w_ij) = 1 for every j
     * where w_ij is the weight of charge on j-th face for the i-th atom. */

    vector<double> charges(n_atoms);
    vector<double> weights;
    for (int face = 0; face < n_faces; ++face) {
        Point3 point1 = faces.get_point(face);
        double q_face = faces.get_charge(face);

        double r_cut2 = faces.get_area(face) * 100.0;
        double sf = smooth_factor * sqrt(r_cut2) / 10.0;

        // Find weights and normalization factor for all the atoms for given face
        // Get the charge for real surface atoms
        weights = vector<double>(n_atoms);
        double w_sum = 0.0;
        for (int atom = 0; atom < n_atoms; ++atom) {
            double dist2 = point1.periodic_distance2(get_point(atom), sizes.xbox, sizes.ybox);
            if (dist2 > r_cut2) continue;

            double w = exp(-1.0 * sqrt(dist2) / sf);
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
        Vec3 force = fields.get_elfield(atom) * (charges[atom] * force_factor);   // [e*V/A]
        interpolation.push_back(Solution(force, charges[atom]));
    }

    histogram_clean(0, r_cut);  // clean by vector x-component
    histogram_clean(1, r_cut);  // clean by vector y-component
    histogram_clean(2, r_cut);  // clean by vector z-component
    histogram_clean(3, r_cut);  // clean by vector norm
    histogram_clean(4, r_cut);  // clean by scalar
}

int ForceReader::export_charge_and_force(const int n_atoms, double* xq) const {
    if (n_atoms < 0) return 0;
    check_return(size() == 0, "No charge & force to export!");

    // Initially pass the zero force and charge for all the atoms
    for (int i = 0; i < 4*n_atoms; ++i)
        xq[i] = 0;

    // Pass the the calculated electric field for stored atoms
    for (int i = 0; i < size(); ++i) {
        int identifier = get_id(i);
        if (identifier < 0 || identifier >= n_atoms) continue;

        identifier *= 4;
        xq[identifier++] = interpolation[i].scalar;
        for (double x : interpolation[i].vector)
            xq[identifier++] = x;
    }

    return 0;
}

int ForceReader::export_force_and_pairpot(const int n_atoms, double* xnp, double* Epair, double* Vpair) const {
    if (n_atoms < 0) return 0;
    check_return(size() == 0, "No force & pair potential to export!");

    Vec3 box(sizes.xbox, sizes.ybox, sizes.zbox);

    for (int i = 0; i < size(); ++i) {
        int id = get_id(i);
        if (id < 0 || id >= n_atoms) continue;

        double V = interpolation[i].norm;
        Epair[id] += V;
        Vpair[0] += V;

        id *= 3;
        Vec3 force = interpolation[i].vector;
        for (int j = 0; j < 3; ++j)
            xnp[id+j] += force[j] / box[j];
    }

    return 0;
}

} // namespace femocs
