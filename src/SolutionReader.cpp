/*
 * SolutionReader.cpp
 *
 *  Created on: 6.6.2016
 *      Author: veske
 */

#include "SolutionReader.h"
#include "Macros.h"
#include "getelec.h"

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
        limit_min(0), limit_max(0), sort_atoms(false), dim(0), rank(0), empty_val(0), interpolator(NULL)
{
    reserve(0);
}

SolutionReader::SolutionReader(Interpolator* i, const string& vec_lab, const string& vec_norm_lab, const string& scal_lab) :
        vec_label(vec_lab), vec_norm_label(vec_norm_lab), scalar_label(scal_lab),
        limit_min(0), limit_max(0), sort_atoms(false), dim(0), rank(0), empty_val(0), interpolator(i)
{
    reserve(0);
}

void SolutionReader::calc_interpolation() {
    require(interpolator, "NULL interpolator cannot be used!");

    const int n_atoms = size();

    const bool b1 = dim == 2 && rank == 1 && interpolator->lintris.size() == 0;
    const bool b2 = dim == 3 && rank == 1 && interpolator->lintets.size() == 0;
    const bool b3 = dim == 2 && rank == 2 && interpolator->quadtris.size() == 0;
    const bool b4 = dim == 3 && rank == 2 && interpolator->quadtets.size() == 0;
    if (b1 || b2 || b3 || b4) {
        interpolation = vector<Solution>(n_atoms, Solution(empty_val));
        return;
    }

    // Sort atoms into sequential order to speed up interpolation
    if (sort_atoms) sort_spatial();
    array<double,8> sf;

    int cell = 0;
    for (int i = 0; i < n_atoms; ++i) {
        Point3 point = get_point(i);

        // Depending on interpolation dimension and rank, pick corresponding functions
        if (dim == 2 && rank == 1) {
            cell = interpolator->lintris.locate_cell(point, abs(cell));
            append_interpolation(interpolator->lintris.interp_solution(point, cell));
        } else if (dim == 2 && rank == 2) {
            cell = interpolator->quadtris.locate_cell(point, abs(cell));
            append_interpolation(interpolator->quadtris.interp_solution(point, cell));
        } else if (dim == 3 && rank == 1) {
            cell = interpolator->lintets.locate_cell(point, abs(cell));
            append_interpolation(interpolator->lintets.interp_solution(point, cell));
        } else if (dim == 3 && rank == 2) {
            cell = interpolator->quadtets.locate_cell(point, abs(cell));
            append_interpolation(interpolator->quadtets.interp_solution(point, cell));
        } else if (dim == 3 && rank == 3) {
            cell = interpolator->linhexs.locate_cell(point, abs(cell));
            append_interpolation(Solution(0));
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
    if (cells_known)
        // ...yes, no need to calculate them again, just interpolate
        for (int i = 0; i < n_atoms; ++i) {
            // locate the face
            int cell = atom2cell[i];
            set_marker(i, cell);

            // calculate the interpolation
            if (dim == 2 && rank == 1)
                append_interpolation(interpolator->lintris.interp_solution(get_point(i), cell));
            else if (dim == 2 && rank == 2)
                append_interpolation(interpolator->quadtris.interp_solution(get_point(i), cell));
            else if (dim == 3 && rank == 1)
                append_interpolation(interpolator->lintets.interp_solution(get_point(i), cell));
            else if (dim == 3 && rank == 2)
                append_interpolation(interpolator->quadtets.interp_solution(get_point(i), cell));
            else if (dim == 3 && rank == 3)
                append_interpolation(Solution(0));
        }

    // ...nop, do it and interpolate
    else {
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

// Separate cylindrical region from substrate region
int SolutionReader::get_nanotip(Medium& nanotip, vector<bool>& atom_in_nanotip, const double radius) {
    const int n_atoms = size();
    const double radius2 = radius * radius;
    Medium::calc_statistics();

    // Make map for atoms in nanotip
    Point2 centre(sizes.xmid, sizes.ymid);
    atom_in_nanotip = vector<bool>(n_atoms);
    for (int i = 0; i < n_atoms; ++i)
        atom_in_nanotip[i] = centre.distance2(get_point2(i)) <= radius2;

    const int n_nanotip_atoms = vector_sum(atom_in_nanotip);

    // Separate nanotip from substrate
    nanotip.reserve(n_nanotip_atoms);
    for (int i = 0; i < n_atoms; ++i)
        if (atom_in_nanotip[i])
            nanotip.append(Atom(i, get_point(i), TYPES.SURFACE));

    nanotip.calc_statistics();
    return n_nanotip_atoms;
}

// Initialise statistics about the solution
void SolutionReader::init_statistics() {
    stat.vec_norm_min = stat.scal_min = DBL_MAX;
    stat.vec_norm_max = stat.scal_max = -DBL_MAX;
}

// Calculate statistics about the solution
void SolutionReader::calc_statistics() {
    init_statistics();

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

/* ==========================================
 * ============== FIELD READER ==============
 * ========================================== */

FieldReader::FieldReader(Interpolator* i) :
        SolutionReader(i, "elfield", "elfield_norm", "potential"),
        E0(0), radius1(0), radius2(0) {}

void FieldReader::test_pic(fch::Laplace<3>* laplace, const Medium& medium) {
    const double x = 0;
    const double y = 0;
    const double zmin = medium.sizes.zmax;
    const double step = 0.5;

    const int n_points = 30;

    reserve(n_points);
    for (int i = 0; i < n_points; ++i)
        append(Point3(x, y, zmin + i * step));

    cout << "\nprobing potential:\n";

    array<double,8> shape_functions;

    int hex_index = 0;
    for (int i = 0; i < n_points; ++i) {
        Point3 p = get_point(i);
        hex_index = interpolator->linhexs.locate_cell(p, hex_index);
        dealii::Point<3> deal_point(p.x, p.y, p.z);
        double val1 = laplace->probe_potential(deal_point,interpolator->linhexs.femocs2deal(hex_index));
        double val2 = laplace->probe_potential(deal_point);
        printf("%.2f, %e, %e, %e\n", p.z, val1, val2, val1 - val2);
    }

    cout << "\nprobing elfield:\n";

    hex_index = 0;
    for (int i = 0; i < n_points; ++i) {
        Point3 p = get_point(i);
        hex_index = interpolator->linhexs.locate_cell(p, hex_index);
        dealii::Point<3> deal_point(p.x, p.y, p.z);

        double val1 = laplace->probe_efield_norm(deal_point,interpolator->linhexs.femocs2deal(hex_index));
        double val2 = laplace->probe_efield_norm(deal_point);
        printf("%.2f, %e, %e, %e\n", p.z, val1, val2, val1 - val2);
    }

    cout << "\nstoring interpolation:\n";

    hex_index = 0;
    for (int i = 0; i < n_points; ++i) {
        Point3 p = get_point(i);
        hex_index = interpolator->linhexs.locate_cell(p, hex_index);
        dealii::Point<3> deal_point(p.x, p.y, p.z);

        double val1 = laplace->probe_efield_norm(deal_point,interpolator->linhexs.femocs2deal(hex_index));
        double val2 = laplace->probe_potential(deal_point,interpolator->linhexs.femocs2deal(hex_index));
        append_interpolation(Solution(Vec3(0), val1, val2));
    }
}

void FieldReader::test_pic_vol2(fch::Laplace<3>* laplace, const Medium& medium, const TetgenMesh& mesh) {
    const double x = 0;
    const double y = 0;
    const double zmin = 0.5 + medium.sizes.zmax;
    const double step = 0.0005;

    const int n_points = 10000;

    reserve(n_points);
    for (int i = 0; i < n_points; ++i)
        append(Point3(x, y, zmin + i * step));

    int cell_index;
    array<double,8> shape_functions;

//    cout << "\nprobing potential:\n";
//
//    cell_index = 0;
//    for (int i = 0; i < n_points; ++i) {
//        Point3 p = get_point(i);
//        cell_index = interpolator->linhexs.locate_cell(p, cell_index); // necessary to per
//        double val1 = interpolator->linhexs.interp_solution(p, cell_index).scalar;
//        double val2 = laplace->probe_potential(dealii::Point<3>(p.x, p.y, p.z));
//        printf("%.2f, %e, %e, %.6f\n", p.z, val1, val2, val1/val2);
//    }
//
//    cout << "\nprobing elfield:\n";
//
//    cell_index = 0;
//    for (int i = 0; i < n_points; ++i) {
//        Point3 p = get_point(i);
//        cell_index = interpolator->linhexs.locate_cell(p, cell_index);
//        double val1 = interpolator->linhexs.interp_solution(p, cell_index).norm;
//        double val2 = laplace->probe_efield_norm(dealii::Point<3>(p.x, p.y, p.z));
//        printf("%.2f, %e, %e, %.6f\n", p.z, val1, val2, val1/val2);
//    }

//    cout << "\ntesting shape functions:\n";
//    cout << setprecision(3) << fixed;
//
//    for (int i = 0; i < mesh.elems.size(); ++i)
//        if (mesh.elems.get_marker(i) == TYPES.VACUUM) {
//            cell_index = 4*i;
//            break;
//        }
//
//    SimpleHex shex = mesh.hexahedra[cell_index];
//    for (int i : shex) {
//        Point3 p = interpolator->nodes.get_vertex(i);
//        interpolator->linhexs.get_shape_functions(shape_functions, p, cell_index);
//
//        double shape_sum = 0;
//        for (double sf : shape_functions) {
//            cout << fabs(sf) << ", ";
//            shape_sum += sf;
//        }
//        cout << "sum=" << shape_sum << endl;
//    }
//
//    cout << endl;

    cout << "\ncomparing interpolations:\n";
    double t0;

    start_msg(t0, "laplace");
    cell_index = 0;
    for (int i = 0; i < n_points; ++i) {
        Point3 p = get_point(i);
        dealii::Point<3> deal_point(p.x, p.y, p.z);

        cell_index = interpolator->linhexs.locate_cell(p, cell_index);
        int cell_index_in_deal = interpolator->linhexs.femocs2deal(cell_index);
        if (cell_index_in_deal < 0)
            append_interpolation(Solution(0));
        else {
            double val1 = laplace->probe_efield_norm(deal_point,cell_index_in_deal);
            double val2 = laplace->probe_potential(deal_point, cell_index_in_deal);
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
        cell_index = interpolator->linhexs.locate_cell(p, cell_index);
        append_interpolation(interpolator->linhexs.interp_solution(p, cell_index));
    }
    end_msg(t0);
    write("out/potential2.xyz");

    interpolation.clear();
    interpolation.reserve(n_points);

    start_msg(t0, "lintets");
    cell_index = 0;
    for (int i = 0; i < n_points; ++i) {
        Point3 p = get_point(i);
        cell_index = interpolator->lintets.locate_cell(p, cell_index);
        append_interpolation(interpolator->lintets.interp_solution(p, cell_index));
    }
    end_msg(t0);
    write("out/potential3.xyz");

    interpolation.clear();
    interpolation.reserve(n_points);

    start_msg(t0, "quadtets");
    cell_index = 0;
    for (int i = 0; i < n_points; ++i) {
        Point3 p = get_point(i);
        cell_index = interpolator->quadtets.locate_cell(p, cell_index);
        append_interpolation(interpolator->quadtets.interp_solution(p, cell_index));
    }
    end_msg(t0);
    write("out/potential4.xyz");
}

void FieldReader::test_pic_vol3(const TetgenMesh& mesh) const {
    int cell_index;
    array<double,8> shape_functions;
    array<double,4> bcc;

    int cntr = 0;
    for (int i = 0; i < mesh.elems.size(); ++i)
        if (mesh.elems.get_marker(i) == TYPES.VACUUM && cntr++ >= 0) {
            cell_index = i;
            break;
        }

    SimpleElement stet = mesh.elems[cell_index];


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

    cout << "results for tet=" << cell_index << ", 4tet=" << 4*cell_index << endl;
    cout << "with neighbours ";
    for (int nbor : mesh.elems.get_neighbours(cell_index))
        cout << nbor << " (" << 4*nbor << ")   ";
    cout << endl;// << fixed << setprecision(3);

    require(labels.size() == points.size(), "Incompatible vectors!");

    for (int i = 0; i < 4; ++i) {
        cell_index = interpolator->linhexs.locate_cell(points[i], 0);
        interpolator->lintets.get_shape_functions(bcc, points[i], 0);
        cout << endl << labels[i] << ":\t" << cell_index+1 << "\t";
        for (double b : bcc)
            cout << ", " << b;
    }

    for (int i = 4; i < labels.size(); ++i) {
        cell_index = interpolator->linhexs.locate_cell(points[i], 0);
        interpolator->lintets.get_shape_functions(bcc, points[i], abs(int(cell_index/4)));
        cout << endl << labels[i] << ":\t" << cell_index+1 << "\t";
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

void FieldReader::interpolate(vector<double>& elfields, const vector<dealii::Point<3>>& nodes) {
    // import the surface nodes the solver needs
    const int n_nodes = nodes.size();

    // store the node coordinates
    reserve(n_nodes);
    int i = 0;
    for (const dealii::Point<3>& node : nodes)
        append( Atom(i++, Point3(node[0], node[1], node[2]), 0) );

    // interpolate solution on the nodes and clean peaks
    calc_interpolation();

    // export electric field norms to the solver
    elfields = vector<double>(n_nodes);
    for (int i = 0; i < n_nodes; ++i)
        elfields[i] = 10.0 * get_elfield_norm(i);
}

// Linearly interpolate electric field for the currents and temperature solver
void FieldReader::transfer_elfield(fch::CurrentsAndHeatingStationary<3>* ch_solver) {
    // import the surface nodes the solver needs
    vector<dealii::Point<3>> nodes;
    ch_solver->get_surface_nodes(nodes);

    // interpolate electric fields and potentials
    vector<double> elfields;
    interpolate(elfields, nodes);

    // transfer electric fields to solver
    ch_solver->set_electric_field_bc(elfields);
}

// Linearly interpolate electric field for the currents and temperature solver
void FieldReader::transfer_elfield(fch::CurrentsAndHeating<3>& ch_solver) {
    // import the surface nodes the solver needs
    vector<dealii::Point<3>> nodes;
    ch_solver.get_surface_nodes(nodes);

    // interpolate electric fields and potentials
    vector<double> elfields;
    interpolate(elfields, nodes);

    // transfer electric fields to solver
    ch_solver.set_electric_field_bc(elfields);
}

// Linearly interpolate electric field on a set of points
void FieldReader::export_elfield(const int n_points, double* Ex, double* Ey, double* Ez, double* Enorm, int* flag) {
    require(n_points == size(), "Invalid query size: " + to_string(n_points));
    for (int i = 0; i < n_points; ++i) {
        Ex[i] = interpolation[i].vector.x;
        Ey[i] = interpolation[i].vector.y;
        Ez[i] = interpolation[i].vector.z;
        Enorm[i] = interpolation[i].norm;
        flag[i] = atoms[i].marker < 0;
    }
}

// Linearly interpolate electric potential on a set of points
void FieldReader::export_potential(const int n_points, double* phi, int* flag) {
    require(n_points == size(), "Invalid query size: " + to_string(n_points));
    for (int i = 0; i < n_points; ++i) {
        phi[i] = interpolation[i].scalar;
        flag[i] = atoms[i].marker < 0;
    }
}

// Export interpolated electric field
void FieldReader::export_solution(const int n_atoms, double* Ex, double* Ey, double* Ez, double* Enorm) {
    if (n_atoms <= 0) return;

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

// Interpolate temperature for the currents and temperature solver
void HeatReader::interpolate(fch::CurrentsAndHeating<3>& ch_solver) {
    // import the surface nodes the solver needs
    vector<dealii::Point<3>> nodes;
    ch_solver.get_surface_nodes(nodes);

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
void HeatReader::export_temperature(const int n_atoms, double* T) {
    if (n_atoms <= 0) return;

    // Pass the the calculated temperature for stored atoms
    for (int i = 0; i < size(); ++i) {
        int identifier = get_id(i);
        if (identifier < 0 || identifier >= n_atoms) continue;
        T[identifier] = get_temperature(i);
    }
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

EmissionReader::EmissionReader(const FieldReader& fr, const HeatReader& hr, const TetgenFaces& f,
        Interpolator* i) :
        SolutionReader(i, "none", "rho_norm", "temperature"),
        fields(fr), heat(hr), faces(f) {
    initialize();
}

void EmissionReader::initialize() {
    int n_nodes = fields.size();

    // deallocate and allocate currents data
    current_densities.clear();
    nottingham.clear();
    currents.clear();
    current_densities.reserve(n_nodes);
    nottingham.reserve(n_nodes);
    currents.reserve(n_nodes);

    //deallocate and allocate lines
    rline.clear();
    Vline.clear();
    rline.reserve(n_lines);
    Vline.reserve(n_lines);

    // find Fmax
    for (int i = 0; i < n_nodes; ++i)
        Fmax = max(Fmax, fields.get_elfield_norm(i));

    //Initialise data
    Jmax = 0.;
    Frep = Fmax;
    Jrep = 0.;
    multiplier = 1.;
}

void EmissionReader::emission_line(const Point3& point, const Vec3& direction, const double rmax) {
    const int interpolation_rank = 1;
    const double nm_per_angstrom = 0.1;
    const double rmin = 0.0;
    Point3 pfield(direction.x, direction.y, direction.z);

    FieldReader fr(interpolator);
    fr.set_preferences(false, 3, interpolation_rank, 0);
    fr.reserve(n_lines);

    for (int i = 0; i < n_lines; i++){
        rline[i] = rmin + ((rmax - rmin) * i) / (n_lines - 1);
        fr.append(point - pfield * rline[i]);
    }
    fr.calc_interpolation();
    for (int i = 0; i < n_lines; i++){
        Vline[i] = multiplier * fr.get_potential(i);
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
    double I_fwhm = 0.; // total emitted current within FWHM emitting area
    double FJ = 0.; // int_FWHMarea (F*J)dS
    double I_tot = 0;

    for (int i = 0; i < fields.size(); ++i){ // go through face centroids
        double face_area = faces.get_area(abs(fields.get_marker(i))) / 3.;
        currents[i] = face_area * current_densities[i];
        I_tot += currents[i];

        if (current_densities[i] > Jmax * 0.5){ //if point eligible
            //quadrangle face area is 1/3 of corresponding triangle face area
            double face_area = faces.get_area(abs(fields.get_marker(i))) / 3.;
            area += face_area; // increase total area
            I_fwhm += currents[i]; // increase total current
            FJ += currents[i] * fields.get_elfield_norm(i);
        }
    }

    if (MODES.VERBOSE)
        printf("I_fwhm = %e, I_tot = %e [Amps]\n", I_fwhm, I_tot);
    Jrep = I_fwhm / area;
    Frep = multiplier * FJ / I_fwhm;
}

void EmissionReader::calc_emission(double workfunction, bool blunt){

    struct emission gt;
    gt.W = workfunction;    // set workfuntion, must be set in conf. script
    gt.R = 1000.0;   // radius of curvature (overrided by femocs potential distribution)
    gt.gamma = 10;  // enhancement factor (overrided by femocs potential distribution)
    double F, J;    // Local field and current density in femocs units (Angstrom)

    for (int i = 0; i < fields.size(); ++i) { // go through all face centroids

        Vec3 field = fields.get_elfield(i); //local field
        F = multiplier * field.norm();
        gt.mode = 0;
        gt.F = angstrom_per_nm * F;
        gt.Temp = heat.get_temperature(i);
        set_marker(i, 0); // set marker for output emission xyz file. Means No full calculation

        if (F > 0.6 * Fmax && !blunt){ // Full calculation with line only for high field points
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

        if (J > 0.1 * Jmax){ // If J is worth it, calculate with full energy integration
            gt.approx = 1;
            cur_dens_c(&gt);
            if (gt.ierr != 0 )
                write_verbose_msg("GETELEC 2nd call returned with error, ierr = " + to_string(gt.ierr));
            J = gt.Jem * nm2_per_angstrom2;
            set_marker(i, 2);
        }

        Jmax = max(Jmax, J); // output data
        current_densities[i] = J;
        nottingham[i] = nm2_per_angstrom2 * gt.heat;
    }

}

void EmissionReader::transfer_emission(fch::CurrentsAndHeating<3>& ch_solver,
        const double workfunction, const double Vappl, bool blunt) {

    const int n_nodes = fields.size();

    reserve(n_nodes);
    atoms = fields.atoms;

    double theta_old = multiplier;
    double err_fact = 0.5, error;

    double Fmax_0 = Fmax;

    for (int i = 0; i < 20; ++i){ // SC calculation loop
        Jmax = 0.;
        Fmax = multiplier * Fmax_0;

        calc_emission(workfunction, blunt);
        calc_representative();

        if (MODES.VERBOSE)
            printf("\nSC j= %d th= %f Jmax= %e Jrep= %e Fmax= %f Frep= %f\n", i, multiplier, Jmax,
                    Jrep , Fmax, Frep);

        if (Vappl <= 0) break; // if Vappl<=0, SC is ignored
        if (i > 5) err_fact *= 0.5; // if not converged in first 6 steps, reduce factor

        // calculate SC multiplier (function coming from getelec)
        multiplier = theta_SC(Jrep / nm2_per_angstrom2, Vappl, angstrom_per_nm * Frep);
        error = multiplier - theta_old;
        multiplier = theta_old + error * err_fact;
        theta_old = multiplier;
        if (abs(error) < 1.e-3) break; //if converged break

    }

    for (int i = 0; i < n_nodes; i++) // append data for surface emission xyz file
        append_interpolation( Solution(Vec3(0), log(current_densities[i]), log(fabs(nottingham[i]))));

    ch_solver.set_emission_bc(current_densities, nottingham); // output data for heat BCs
}

/* ==========================================
 * ============== CHARGE READER =============
 * ========================================== */

ChargeReader::ChargeReader(Interpolator* i) :
        SolutionReader(i, "elfield", "area", "charge"), Q_tot(0) {}

void ChargeReader::calc_charges(const TetgenMesh& mesh, const double E0) {
    const double sign = fabs(E0) / E0;
    const int n_faces = mesh.faces.size();
    const int n_quads_per_triangle = 3;

    // Store the centroids of the triangles
    reserve(n_faces);
    for (int i = 0; i < n_faces; ++i)
        append( Atom(i, mesh.faces.get_centroid(i), 0) );

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
        double area = mesh.faces.get_area(face);
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
        SolutionReader(i, "force", "force_norm", "charge") {}

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

int ForceReader::calc_voronois(VoronoiMesh& mesh, vector<bool>& atom_in_nanotip,
        const vector<int>& atom2face, const double radius, const double latconst, const string& mesh_quality)
{
    require(interpolator, "NULL interpolator cannot be used!");
    const int n_atoms = size();
    const bool faces_known = n_atoms == atom2face.size();
    // TODO put those values to Config, because they affect heavily how the voronoi charges will look like
    const double max_distance_from_surface = 0.5 * latconst;
    const double shift_distance = 1.0 * latconst;

    Medium nanotip;
    const int n_nanotip_atoms = get_nanotip(nanotip, atom_in_nanotip, radius);

    // calculate support points for the nanotip
    // by moving the nanotip points in direction of its corresponding triangle norm by r_cut
    Medium support(n_nanotip_atoms);
    int face = 0;
    for (int i = 0; i < n_atoms; ++i)
        if (atom_in_nanotip[i]) {
            Point3 point = get_point(i);
            if (faces_known)
                face = atom2face[i];
            else if (rank == 1)
                face = abs( interpolator->lintris.locate_cell(point, face) );
            else if (rank == 2)
                face = abs( interpolator->quadtris.locate_cell(point, face) );

            if (interpolator->lintris.fast_distance(point, face) < max_distance_from_surface) {
                Vec3 shift = interpolator->lintris.get_norm(face) * shift_distance;
                point += Point3(shift.x, shift.y, shift.z);
                support.append(Atom(i, point, TYPES.VACANCY));
            }
        }

    nanotip += support;
    nanotip.calc_statistics();
    nanotip.write("out/nanotip.xyz");

    // Generate Voronoi cells around the nanotip
    // r - reconstruct, v - output Voronoi cells, Q - quiet, q - mesh quality
    int err_code = mesh.generate(nanotip, latconst, "rQq" + mesh_quality, "vQ");
    if (err_code) return err_code;

    // Clean the mesh from faces and cells that have node in the infinity
    mesh.clean();

    // Save the location of nanotip and support atoms to speed up later calculations
    mesh.nodes.save_indices(n_nanotip_atoms, 0, support.size());

    return 0;
}

int ForceReader::calc_voronoi_charges(VoronoiMesh& mesh, const vector<int>& atom2face, const FieldReader& fields,
         const double radius, const double latconst, const string& mesh_quality)
{
    // Extract nanotip and generate Voronoi cells around it
    vector<bool> atom_in_nanotip;
    int err_code = calc_voronois(mesh, atom_in_nanotip, atom2face, radius, latconst, mesh_quality);
    if (err_code) return err_code;

    const int nanotip_end = mesh.nodes.indxs.surf_end;
    const int support_start = mesh.nodes.indxs.vacuum_start;
    const int support_end = mesh.nodes.indxs.vacuum_end;

    require(mesh.nodes.size() > 0, "Empty Voronoi mesh cannot be handled!");
    require(mesh.voros.size() > 0, "Empty Voronoi mesh cannot be handled!");

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

    // remove too big surface faces
//    clean_voro_faces(mesh);

    // calculate charge on Voronoi cell
    int cell = -1;
    for (int i = 0; i < size(); ++i)
        if (atom_in_nanotip[i]) {
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
            interpolation[i] = Solution(field * charge, charge);
        }

    return 0;
}

// Calculate forces from atomic electric fields and face charges
void ForceReader::calc_forces(const FieldReader &fields) {
    const int n_atoms = fields.size();

    // Copy the atom data
    reserve(n_atoms);
    atoms = fields.atoms;

    // Calculate the charges by ensuring the that the total sum of it remains conserved
    vector<double> charges;
    interpolator->lintris.interp_conserved(charges, atoms);

    // calculate forces and store them
    for (int i = 0; i < n_atoms; ++i) {
        Vec3 force = fields.get_elfield(i) * (charges[i] * force_factor);   // [e*V/A]
        interpolation.push_back(Solution(force, charges[i]));
    }
}

// Calculate forces from atomic electric fields and face charges
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

// Export the induced charge and force on imported atoms
void ForceReader::export_force(const int n_atoms, double* xq) {
    if (n_atoms <= 0) return;

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
}

Vec3 ForceReader::get_force(const int i) const {
    require(i >= 0 && i < size(), "Invalid index: " + to_string(i));
    return interpolation[i].vector;
}

double ForceReader::get_force_norm(const int i) const {
    require(i >= 0 && i < size(), "Invalid index: " + to_string(i));
    return interpolation[i].norm;
}

double ForceReader::get_charge(const int i) const {
    require(i >= 0 && i < size(), "Invalid index: " + to_string(i));
    return interpolation[i].scalar;
}

/* ==========================================
 * ============= CoulombReader =============
 * ========================================== */

/**
 * @brief Module containing subrotuines for calculating the Coloumb force between atoms
 * @brief Add forces due to Coloumb interaction and Lorentz force to atoms
 * The subroutine adds the Coulomb force between atoms based on the charge in the "first column" in @p xq.
 * The Lorentz force is calculated elsewhere, but the force itself is added to the total force here.
 * The Lorentz force components are also stored in @p xq.
 */

/*

void CoulombReader::init(vector<int>& charged, double* xnp,
        const double* xq, const double* _box, const double* pbc, const double qrcut, const int natoms) {

    static bool firsttime = true;
    charged.reserve(natoms / 3);  // List of atoms with charge

    // Figure out size of image cells and neighbour cells
    if (firsttime) {
        box = Vec3(_box[0], _box[1], _box[2]);

        imxmax = ceil(qrcut / box[0]);
        if (pbc[0] == 0) imxmax = 0;
        imymax = ceil(qrcut / box[1]);
        if (pbc[1] == 0) imymax = 0;
        imzmax = ceil(qrcut / box[2]);
        if (pbc[2] == 0) imzmax = 0;

        neigh_cell_size = {2.01*qrcut, 2.01*qrcut, 2.01*qrcut};
        // Assume max 50% higher atomic density than in crystal
        max_in_cell = ceil(1.5*0.85 * neigh_cell_size[0] * neigh_cell_size[1] * neigh_cell_size[2]);

        for (int j = 0; j < 3; ++j)
          ncell[j] = ceil(box[j] / neigh_cell_size[j]);

        if(MODES.VERBOSE) printf("  qforces initial range of image cells %i %i %i", imxmax,imymax,imzmax);
        if(MODES.VERBOSE) printf("neigh_cells %i, %i, %i, %i", ncell[0], ncell[1], ncell[2], max_in_cell);

        neigh_cells[0].reserve(ncell[0]);
        neigh_cells[1].reserve(ncell[1]);
        neigh_cells[2].reserve(ncell[2]);
        neigh_cells[3].reserve(max_in_cell);
    }

    for (int i = 0; i < natoms; ++i) {
        int i3 = i * 3;
        int i4 = i * 4;

        // Transform any outer forces from Parcas units into eV/A for force routines
        for (int j = 0; j < 3; ++j)
            xnp[i3+j] *= box[j];

        //  Make a list of charged atoms. Later we only need to loop over these
        if(abs(xq[i4]) > 0.0)
            charged.push_back(i);
    }

    firsttime = false;
}

void CoulombReader::qforces(
        const double* x0,    ///< Atom positions (parcas units)
        double* xnp,         ///<  Forces on atoms (parcas units)
        const double* _box,  ///<  Simulation box size (Å)
        const double* pbc,   ///<  Periodic boundaries
        double* Epair,       ///<  Potential energy per atom
        const double* xq,    ///<  Charges on atoms (unit charges) and Lorentz force components
        double Vpair,        ///<  Total potential energy of atoms. Pot. due to Coloumb forces are added here. NOTE: Lorentz is missing!
        double Vqq,          ///<  Potnetial energy due to coloumb interaction
        const double qrcut,  ///<  Cut-off for Coloumb force
        const double qscreen, ///<  Screening factor for Coulomb force
        const int natoms     ///<  Number of atoms
        ) {

    const int update_neighbours_every = 10;
    const double couloumb_constant = 14.399758;
    const double r_cut_square = qrcut * qrcut;

    static int times_called = 0;

    double t0;
    array<int,3> cell, nei;
    array<double,3> cellc; // Center-point of cell
    vector<int> charged;

    start_msg(t0, "Initializing Coulomb forces...");
    init(charged, xnp, xq, _box, pbc, qrcut, natoms);
    Vqq = 0;
    end_msg(t0);

    // Divide system into cells, update which atoms belong where
    if (times_called++ % update_neighbours_every == 0) {
        start_msg(t0, "Building neighbour list for qforces cells...");
        calc_nborlist(charged, x0);
        end_msg(t0);
    }

    start_msg(t0, "Looping over image cells...");

    for (int imx = -imxmax; imx <= imxmax; ++imx) {
    for (int imy = -imymax; imy <= imymax; ++imy) {
    for (int imz = -imzmax; imz <= imzmax; ++imz) {
        // Determine whether we should include atoms in this image cell
        // (idea similar to fig. 5.7 in Allen-Tildesley)
        // by checking whether cell _center_ is farther than qrcut
        // from the central cell

        double rsq = imx*imx*box[0]*box[0] + imy*imy*box[1]*box[1] + imz*imz*box[2]*box[2];
        if (rsq > r_cut_square) continue;

        // Once we have chosen to handle one image cell, take
        // care of all atoms in there

//      $OMP PARALLEL DO DEFAULT(SHARED), &
//      $OMP PRIVATE(i,i3,i4,xi,yi,zi,t,t3,t4,&
//      $OMP xtim,ytim,ztim,dx,dy,dz,rsq,r,kCoulomb_q1q2,qscreenfact,V,dVdr), &
//      $OMP REDUCTION(+:Vqq), REDUCTION(+:Vpair)

        for (int i : charged) {
            int i3 = i*3;
            int i4 = i*4;

            Vec3 xx(x0[i3+0], x0[i3+1], x0[i3+2]);
            xx *= box;

            for (int j = 0; j < 3; ++j) {
                cell[j] = ceil( (x0[i3+j]*box[j] + box[j]/2.0) / neigh_cell_size[j] );
                cellc[j] = cell[j] * neigh_cell_size[j] - neigh_cell_size[j] / 2.0;
                if (cell[j] > ncell[j]) cell[j] = ncell[j];
            }

            for (int ineiz = -1; ineiz <= 1; ++ineiz) {
                if (check_limits(xx, cellc, nei, ncell,cell, ineiz, 2)) continue;

                for (int ineiy = -1; ineiy <= 1; ++ineiy) {
                    if (check_limits(xx, cellc, nei, ncell,cell, ineiy, 1)) continue;

                    for (int ineix = -1; ineix <= 1; ++ineix) {
                        if (check_limits(xx, cellc, nei, ncell,cell, ineix, 0)) continue;

                        // loop over atoms s which are neighbours of i
                        for (int cht = 1; cht < neigh_cells(nei[0],nei[1],nei[2],0); ++cht) {
                            int t = neigh_cells(nei[0],nei[1],nei[2],cht);
                            int t4 = t*4;
                            int t3 = t*3;
                            if (xq[t4] == 0) continue;
                            if (i==t && imx==0 && imy==0 && imz==0) continue;

                            // Find distance. Because image cells are used, periodics are not needed.
                            Vec3 point(x0[t3+0]+imx, x0[t3+1]+imy, x0[t3+2]+imz);
                            point *= box;
                            point -= xx;

                            rsq = point.norm2();
                            double r=sqrt(rsq);

                            // With cutoff 20.0 Vqq differs by 0.0002% from w.o. cutoff
                            if (r > qrcut) continue;

                            double V = exp(-qscreen * r) * couloumb_constant * xq[i4] * xq[t4] / r;
                            Epair[i] += 0.5 * V;
                            Vqq += 0.5 * V;
                            Vpair += 0.5 * V;

                            double dVdr=-V/rsq;
                            // rsq because one 1/r from derivative, another from
                            // vector projection

                            for (int j = 0; j < 3; ++j)
                                xnp[i3+j] += dVdr * point[j];
                        }
                    }
                }
            }
        }
//      $OMP END PARALLEL DO
    }
    }
    }

    end_msg(t0);

//    $OMP PARALLEL DO PRIVATE(i3,i4)
    // TODO check whether unit transformation from eV-A system to parcas units is needed
    for (int i = 0; i < natoms; ++i) {
        int i3 = i * 3;
        int i4 = i * 4;

        for (int j = 0; j < 3; ++j)
            xnp[i3+j] = (xnp[i3+j] + xq[i4+j+1]) / box[j];
    }
//    $OMP END PARALLEL DO
}

bool CoulombReader::check_limits(Vec3& xx, array<double,3>& cellc,
        array<int,3>& nei, array<int,3>& ncell, array<int,3>& cell, int inei, int i)
{
    // check if past midpoint for not looking the other way
    if (xx[i] > cellc[i] && inei == -1) return true;
    if (xx[i] < cellc[i] && inei == 1) return true;

    nei[i] = cell[i] + inei;

    // check periodic boundaries
    if (i == 2)
        return nei[i] > ncell[i] || nei[i] < 1;

    if (nei[i] > ncell[i]) nei[i] = 1;
    if (nei[i] < 1) nei[i] = ncell[i];

    return false;
}

void CoulombReader::calc_nborlist(const bool lateral_periodic) {
    const int n_atoms = size();
    const double r_cut = 6.0;
    Point3 simubox_size(sizes.xbox, sizes.ybox, sizes.zbox);

    vector<unsigned> list(n_atoms);
    vector<unsigned> head(n_atoms);

    array<int,3> M;
    for (int j = 0; j < 3; ++j)
        M[j] = ceil(simubox_size[j] / r_cut);

    // calculate linked list for the atoms
    for (int i = 0; i < n_atoms; ++i) {
        Point3 dx = atoms[i].point + simubox_size / 2;
        // Check that we are inside lateral boundaries
        if (lateral_periodic)
            for (int j = 0; j < 2; ++j) {
                if (dx[j] < 0) dx[j] += simubox_size[j];
                else if (dx[j] > simubox_size[j]) dx[j] -= simubox_size[j];
            }

        array<int,3> point_index;
        for (int j = 0; j < 3; ++j)
            point_index[j] = int( (dx[j] / simubox_size[j]) * M[j] );

        // If not periodic, let border cells continue to infinity
        if (!lateral_periodic)
            for (int j = 0; j < 2; ++j) {
                point_index[j] = max(0, point_index[j]);
                point_index[j] = min(M[j]-1, point_index[j]);
            }

        int i_cell = (point_index[2] * M[1] + point_index[1]) * M[0] + point_index[0];
        set_marker(i, i_cell);
        list[i] = head[i_cell];
        head[i_cell] = i;
    }

    for (int i = 0; i < n_atoms; ++i) {
        Point3 &p1 = atoms[i].point;
        for (int ix = 0; ix < M[0]; ++ix) {
            for (int iy = 0; iy < M[1]; ++iy) {
                for (int iz = 0; iz < M[2]; ++iz) {
                    int i_cell = (iz * M[1] + iy) * M[0] + ix;
                    int j = head[i_cell];
                    while(true) {
                        if (j == 0) break;
                        double r2 = p1.distance2(get_point(j));
                        if (r2 < r_cut2) {
                            printf("%i and %i are neighbours!\n", i, j);
                        }
                        j = list[j];
                    }
                }
            }
        }
    }
}

void CoulombReader::check_neighbours(const int i)

void CoulombReader::calc_nborlist(const vector<int>& charged, const double* x0) {
    neigh_cells = -1;

    array<int,3> cell;
    for (int i : charged) {
        int i3 = i * 3;
        int i4 = i * 4;

        for (int j = 0; j < 3; ++j) {
            cell[j] = ceil( (x0[i3+j] * box[j] + box[j] / 2.0) / neigh_cell_size[j] );
            if (cell[j] > ncell[j]) cell[j] = ncell[j];
            if(cell[j] < 1) cell[j] = 1; // Can happen if atom is exactly at border?
        }

        int atoms_in_cell = neigh_cells(cell[0], cell[1], cell[2], 0);
        if (atoms_in_cell >= max_in_cell) {
            printf("Qforces: too many atoms in cell %i, %i", atoms_in_cell,max_in_cell);
            return;
        }

        atoms_in_cell = max(atoms_in_cell, 0);

        atoms_in_cell++;
        neigh_cells(cell[0], cell[1], cell[2], 0) = atoms_in_cell;
        neigh_cells(cell[0], cell[1], cell[2], atoms_in_cell) = i;
    }
}

//*/

} // namespace femocs
