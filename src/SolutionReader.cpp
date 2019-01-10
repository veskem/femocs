/*
 * SolutionReader.cpp
 *
 *  Created on: 6.6.2016
 *      Author: veske
 */

#include "SolutionReader.h"
#include "Macros.h"
#include "Config.h"
#include "FileWriter.h"

#include <float.h>
#include <stdio.h>

using namespace std;
namespace femocs {

/* ==========================================
 * ============= SOLUTION READER ============
 * ========================================== */

SolutionReader::SolutionReader() :
        vec_label("vec"), norm_label("vec_norm"), scalar_label("scalar"),
        limit_min(0), limit_max(0), sort_atoms(false), dim(0), rank(0), interpolator(NULL)
{
    reserve(0);
}

SolutionReader::SolutionReader(Interpolator* i, const string& vec_lab, const string& vec_norm_lab, const string& scal_lab) :
        vec_label(vec_lab), norm_label(vec_norm_lab), scalar_label(scal_lab),
        limit_min(0), limit_max(0), sort_atoms(false), dim(0), rank(0), interpolator(i)
{
    reserve(0);
}

// that function can't be directly moved into calc_interpolation, as OpenMP can't handle reference to cell
int SolutionReader::update_interpolation(const int i, int cell) {
    Atom &atom = atoms[i];

    // Depending on interpolation dimension and rank, pick corresponding functions
    if (dim == 2) {
        if (rank == 1)
            interpolation[i] = interpolator->lintri.locate_interpolate(atom.point, cell);
        else if (rank == 2)
            interpolation[i] = interpolator->quadtri.locate_interpolate(atom.point, cell);
        else if (rank == 3)
            interpolation[i] = interpolator->linquad.locate_interpolate(atom.point, cell);
    } else {
        if (rank == 1)
            interpolation[i] = interpolator->lintet.locate_interpolate(atom.point, cell);
        else if (rank == 2)
            interpolation[i] = interpolator->quadtet.locate_interpolate(atom.point, cell);
        else if (rank == 3)
            interpolation[i] = interpolator->linhex.locate_interpolate(atom.point, cell);
    }

    if (!sort_atoms) atom.marker = cell;
    return cell;
}

void SolutionReader::calc_full_interpolation() {
    require(interpolator, "NULL interpolator cannot be used!");
    const int n_atoms = size();

    // Sort atoms into sequential order to speed up interpolation
    if (sort_atoms) sort_spatial();

    int cell = -1;

#pragma omp parallel for private(cell)
    for (int i = 0; i < n_atoms; ++i)
        cell = update_interpolation(i, cell);

    // Sort atoms back to their initial order
    if (sort_atoms) {
        for (int i = 0; i < n_atoms; ++i)
            interpolation[i].id = atoms[i].id;
        std::sort( interpolation.begin(), interpolation.end(), Solution::sort_up() );
        std::sort( atoms.begin(), atoms.end(), Atom::sort_id() );
    }

    atoms_mapped_to_cells = true;
}

void SolutionReader::calc_interpolation() {
    require(interpolator, "NULL interpolator cannot be used!");

    const int n_atoms = size();

    // are the atoms already mapped against the triangles?
    if (!atoms_mapped_to_cells) {
        // ...nop, do the mapping, interpolate and output mapping
        calc_full_interpolation();
        return;
    }

    // ...yes, no need to calculate them again, just interpolate
    int cell;

#pragma omp parallel for private(cell)
    for (int i = 0; i < n_atoms; ++i) {
        // locate the face
        Atom &atom = atoms[i];
        cell = abs(atom.marker);

        // calculate the interpolation
        if (dim == 2) {
            if (rank == 1)
                interpolation[i] = interpolator->lintri.interp_solution(atom.point, cell);
            else if (rank == 2)
                interpolation[i] = interpolator->quadtri.interp_solution(atom.point, cell);
            else if (rank == 3)
                interpolation[i] = interpolator->linquad.interp_solution(atom.point, cell);
        } else {
            if (rank == 1)
                interpolation[i] = interpolator->lintet.interp_solution(atom.point, cell);
            else if (rank == 2)
                interpolation[i] = interpolator->quadtet.interp_solution(atom.point, cell);
            else if (rank == 3)
                interpolation[i] = interpolator->linhex.interp_solution(atom.point, cell);
        }
    }
}

void SolutionReader::reserve(const int n_nodes) {
    require(n_nodes >= 0, "Invalid number of nodes: " + to_string(n_nodes));

    atoms.clear();
    atoms.reserve(n_nodes);
    interpolation.resize(n_nodes);
    atoms_mapped_to_cells = false;
}

void SolutionReader::write_xyz(ofstream &out) const {
    // write the start of xyz header
    FileWriter::write_xyz(out);

    // write Ovito header
    out << "properties=id:I:1:pos:R:3:marker:I:1:force:R:3:"
            + vec_label + ":R:1:" + norm_label + ":R:1:" + scalar_label + ":R:1\n";

    // write data
    const int n_atoms = size();
    for (int i = 0; i < n_atoms; ++i)
        out << atoms[i] << " " << interpolation[i] << "\n";
}

void SolutionReader::write_vtk_point_data(ofstream& out) const {
    // write point data header, ID-s and markers
    Medium::write_vtk_point_data(out);

    const int n_atoms = size();

    // output scalar (electric potential, temperature etc)
    out << "SCALARS " << scalar_label << " double\nLOOKUP_TABLE default\n";
    for (int i = 0; i < n_atoms; ++i)
        out << interpolation[i].scalar << "\n";

    // output vector magnitude explicitly to make it possible to apply filters in ParaView
    out << "SCALARS " << norm_label << " double\nLOOKUP_TABLE default\n";
    for (int i = 0; i < n_atoms; ++i)
        out << interpolation[i].norm << "\n";

    // output vector data (electric field, current density etc)
    out << "VECTORS " << vec_label << " double\n";
    for (int i = 0; i < n_atoms; ++i)
        out << interpolation[i].vector << "\n";
}

void SolutionReader::init_statistics() {
    Medium::init_statistics();
    stat.vec_norm_min = stat.scal_min = DBL_MAX;
    stat.vec_norm_max = stat.scal_max = -DBL_MAX;
}

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

int SolutionReader::contains(const string& data_label) const {
    if (data_label == vec_label) return 1;
    if (data_label == norm_label) return 2;
    if (data_label == scalar_label) return 3;

    return 0;
}

int SolutionReader::export_results(const int n_points, const string &data_type, double* data) const {
    check_return(size() == 0, "No " + data_type + " to export!");

    int dtype = contains(data_type);
    require(dtype > 0 && dtype <= 3, "SolutionReader does not contain " + data_type);

    // Pass the desired solution
    for (int i = 0; i < size(); ++i) {
        int id = get_id(i);
        if (id < 0 || id >= n_points) continue;

        switch(dtype) {
            case 1:
                id *= 3;
                for (double v : interpolation[i].vector)
                    data[id++] = v;
                continue;
            case 2:
                data[id] = interpolation[i].norm;
                continue;
            default:
                data[id] = interpolation[i].scalar;
        }
    }

    return 0;
}

int SolutionReader::interpolate_results(const int n_points, const string &data_type, const double* x,
        const double* y, const double* z, double* data) {
    check_return(size() == 0, "No " + data_type + " to interpolate!");

    // transfer coordinates
    SolutionReader sr(interpolator, vec_label, norm_label, scalar_label);
    sr.set_preferences(sort_atoms, dim, rank);
    sr.reserve(n_points);
    for (int i = 0; i < n_points; ++i)
        sr.append(Atom(i, Point3(x[i], y[i], z[i]), 0));

    // interpolate solution
    sr.calc_interpolation();

    // export interpolation
    return sr.export_results(n_points, data_type, data);
}

void SolutionReader::interpolate(const DealSolver<3>& solver) {
    solver.export_surface_centroids(*this);
    calc_interpolation();
}

void SolutionReader::interpolate(const int n_points, const double* x, const double* y, const double* z) {
    // store the point coordinates
    reserve(n_points);
    for (int i = 0; i < n_points; ++i)
        append(Atom(i, Point3(x[i], y[i], z[i]), 0));

    // interpolate solution
    calc_interpolation();
}

void SolutionReader::interpolate(const Medium &medium) {
    const int n_atoms = medium.size();
    reserve(n_atoms);

    // store atom coordinates
    for (int i = 0; i < n_atoms; ++i)
        append( Atom(i, medium.get_point(i), 0) );

    // interpolate solution
    calc_interpolation();

    // restore original atom id-s
    for (int i = 0; i < n_atoms; ++i)
        atoms[i].id = medium.get_id(i);
}

void SolutionReader::interpolate(const AtomReader &reader) {
    require(!sort_atoms, "Atom sorting not allowed here!");
    const int n_atoms = reader.size();
    reserve(n_atoms);

    // store atom coordinates
    for (int i = 0; i < n_atoms; ++i)
        if (reader.get_marker(i) != TYPES.FIXED)
            append(reader.get_atom(i));

    // interpolate solution
    calc_interpolation();
}

/* ==========================================
 * ============== FIELD READER ==============
 * ========================================== */

FieldReader::FieldReader(Interpolator* i) :
        SolutionReader(i, LABELS.elfield, LABELS.elfield_norm, LABELS.potential),
        E0(0), radius1(0), radius2(0), beta(0) {}

double FieldReader::get_analyt_potential(const int i, const Point3& origin) const {
    require(i >= 0 && i < size(), "Invalid index: " + to_string(i));

    Point3 point = get_point(i);
    point -= origin;
    double r = point.distance(Point3(0));
    return -E0 * point.z * (1 - pow(radius1 / r, 3.0));
}

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

bool FieldReader::check_limits(const vector<Solution>* solutions, bool verbose) {
    if (limit_min == limit_max)
        return false;

    double Emax = calc_max_field(solutions);
    const double gamma1 = fabs(Emax / E0);
    const double gamma2 = get_analyt_enhancement();
    beta = fabs(gamma1 / gamma2);

    bool out_of_limits = beta < limit_min || beta > limit_max;
    if (verbose || out_of_limits) {
        stringstream stream;
        stream << fixed << setprecision(3);
        stream << "field enhancements:  (M)easured:" << gamma1
                << "  (A)nalyt:" << gamma2
                << "  M/A:" << gamma1 / gamma2;
        write_verbose_msg(stream.str());
    }
    return out_of_limits;
}

double FieldReader::calc_max_field(const vector<Solution>* solutions) const {
    double E_max = -1e100;
    if (solutions) {
        for (Solution s : *solutions)
            E_max = max(E_max, s.norm);
    } else {
        for (Solution s : interpolation)
            E_max = max(E_max, s.norm);
    }
    return E_max;
}

void FieldReader::set_check_params(const Config& conf, double radius,
        double tip_height, double box_height)
{
    require(conf.field.anode_BC == "neumann" || conf.field.anode_BC == "dirichlet",
            "Invalid anode boundary condition: " + conf.field.anode_BC);

    if (conf.field.anode_BC == "neumann")
        E0 = conf.field.E0;
    else
        E0 = conf.field.V0 / box_height;
    limit_min = conf.tolerance.field_min;
    limit_max = conf.tolerance.field_max;

    if (tip_height > radius) {
        radius2 = tip_height;
        radius1 = sqrt(tip_height * radius);
    } else {
        radius2 = radius;
        radius1 = radius;
    }
}

/* ==========================================
 * =============== HEAT READER ==============
 * ========================================== */

HeatReader::HeatReader(Interpolator* i) :
        SolutionReader(i, LABELS.rho, LABELS.potential, LABELS.temperature)
{}

void HeatReader::sort_spatial() {
    std::sort( atoms.begin(), atoms.end(), Atom::sort_marker_up() );
}

void HeatReader::sort_spatial(const TetgenMesh* mesh) {
    const int n_nodes = mesh->nodes.size();
    const int n_tets = mesh->tets.size();
    const int n_dealii_nodes = size();

    // determine tetrahedra that are connected to node
    vector<vector<int>> node2tets(n_nodes);
    for (int tet = 0; tet < n_tets; ++tet) {
        for (int node : mesh->tets[tet])
            node2tets[node].push_back(tet);
    }

    // obtain sort indices for Femocs mesh nodes that are also present in Deal.II
    vector<int> indices(n_nodes, -1);
    int cntr = 0;

    for (int i = 0; i < n_nodes; ++i)
        for (int tet : node2tets[i]) {
            // skip tetrahedra that are located in  vacuum
            if (mesh->tets.get_marker(tet) == TYPES.VACUUM) continue;

            // enumerate nodes of tetrahedron
            for (int node : mesh->tets[tet]) {
                if (indices[node] < 0)
                    indices[node] = cntr++;
            }

            // enumerate remaining nodes of hexahedra
            // that are located inside the tetrahedron
            for (int hex : mesh->tets.to_hexs(tet)) {
                for (int node : mesh->hexs[hex])
                    if (indices[node] < 0)
                        indices[node] = cntr++;
            }
        }

    // must be below 1 to ensure nodes are shrinked into the mesh
    // if > 1, nodes are expanded and therefore some nodes will move outside the mesh making interpolating slower
    const double shrink_factor = 0.99;

    cntr = 0;
    for (int i : indices)
        if (i >= 0) {
            require(cntr < n_dealii_nodes, "Index of Deal.II vertex exceeds # of vertices: "
                    + d2s(cntr) + " >= " + d2s(n_dealii_nodes));
            // store sort index
            atoms[cntr].marker = i;
            // move point little bit to ensure it doesn't overlap with previous mesh node
            atoms[cntr].point *= shrink_factor;
            cntr++;
        }
}

void HeatReader::interpolate_dofs(CurrentHeatSolver<3>& solver, const TetgenMesh* mesh) {
    solver.heat.export_vertices(*this);
    if (sort_atoms) sort_spatial(mesh);
    calc_interpolation();

    // transfer temperatures and potentials into separate vectors
    const int n_points = size();
    temperatures.resize(n_points);
    for (int i = 0; i < n_points; ++i)
        temperatures[i] = interpolation[i].scalar;

    solver.heat.import_solution(&temperatures);
}

void HeatReader::precalc_berendsen() {
    const int n_tets = interpolator->lintet.size();

    // calculate mapping between atoms and tets that surround them
    tet2atoms = vector<vector<int>>(n_tets);
    for (int i = 0; i < size(); ++i)
        tet2atoms[abs(get_marker(i))].push_back(i);

    // calculate average temperature inside tetrahedron
    fem_temp = vector<double>(n_tets);
    for (int tet = 0; tet < n_tets; ++tet) {
        int n_atoms_in_tet = tet2atoms[tet].size();
        if (n_atoms_in_tet > 0) {
            for (int node : interpolator->lintet.get_cell(tet))
                fem_temp[tet] += interpolator->nodes.get_scalar(node);
            fem_temp[tet] /= n_nodes_per_tet;
        }
    }
}

int HeatReader::scale_berendsen(double* x1, const int n_atoms,
        const Vec3& parcas2si, const Config& conf)
{
    check_return(size() == 0, "No " + LABELS.parcas_velocity + " to export!");

    // PARCAS internal time step in fs
    // needed to transfer velocity from PARCAS units to Angstrom / fs
    const double delta_t = conf.behaviour.timestep_fs / (10.1805*sqrt(conf.behaviour.mass));
    // MD time step over Berendsen tau
    timestep_over_tau = conf.behaviour.timestep_fs / conf.heating.tau;

    const unsigned int n_tets = tet2atoms.size();
    require(n_tets > 0, "Data is missing for long Berendsen thermostat!");
    require(fem_temp.size() == n_tets, "Mismatch between vector sizes: "
            + d2s(n_tets) + ", " + d2s(fem_temp.size()));

    // convert velocities from Parcas to SI units
    vector<Vec3> velocities;
    calc_SI_velocities(velocities, n_atoms, parcas2si / delta_t, x1);

    for (unsigned int tet = 0; tet < n_tets; ++tet) {
        int n_atoms_in_tet = tet2atoms[tet].size();
        if (n_atoms_in_tet) {

            // calculate average temperature of atoms inside a tetrahedron
            double md_temp = 0;
            for (int atom : tet2atoms[tet]) {
                int id = get_id(atom);
                if (id < 0 || id >= n_atoms) continue;
                md_temp += velocities[id].norm2();
            }
            md_temp *= heat_factor / n_atoms_in_tet;

            // calculate scaling factor
            double lambda = calc_lambda(md_temp, fem_temp[tet]);

            // scale the velocities towards tetrahedron temperature
            for (int atom : tet2atoms[tet]) {
                int id = get_id(atom);
                if (id < 0 || id >= n_atoms) continue;

                int I = 3 * id;
                for (int j = 0; j < 3; ++j)
                    x1[I+j] *= lambda;

                // store scaled temperature
                interpolation[atom].scalar = velocities[id].norm2() * (heat_factor * lambda);
            }
        }
    }

    write("out/berendsen.movie");

    return 0;
}

inline double HeatReader::calc_lambda(const double T_start, const double T_end) const {
    return sqrt( 1.0 + timestep_over_tau * (T_end / T_start - 1.0) );
}

void HeatReader::calc_SI_velocities(vector<Vec3>& velocities,
        const int n_atoms, const Vec3& parcas2si, double* x1) {

    // perform conversion
    velocities.clear();
    velocities.reserve(n_atoms);
    for (int i = 0; i < n_atoms; ++i) {
        int I = 3*i;
        velocities.push_back(Vec3(x1[I], x1[I+1], x1[I+2]) * parcas2si);
    }
}

/* ==========================================
 * ============== CHARGE READER =============
 * ========================================== */

ChargeReader::ChargeReader(Interpolator* i) :
        SolutionReader(i, LABELS.elfield, LABELS.area, LABELS.charge), Q_tot(0) {}

void ChargeReader::calc_charges(const TetgenMesh& mesh, const double E0) {
    const double sign = fabs(E0) / E0;
    const int n_faces = mesh.tris.size();

    // Store the centroids of the triangles
    reserve(n_faces);
    for (int i = 0; i < n_faces; ++i)
        append( Atom(i, mesh.tris.get_centroid(i), 0) );

    // create triangle index to its centroid index mapping
    vector<int> tri2centroid(n_faces);
    for (int face = 0; face < n_faces; ++face) {
        for (int node : mesh.quads[n_quads_per_tri * face])
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
        interpolation[face] = Solution(elfield, area, charge);
    }
}

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

void ChargeReader::set_check_params(const double Q_tot, const double limit_min, const double limit_max) {
    this->Q_tot = Q_tot * eps0;
    this->limit_min = limit_min;
    this->limit_max = limit_max;
}

/* ==========================================
 * ============== FORCE READER ==============
 * ========================================== */

ForceReader::ForceReader(Interpolator* i) :
        SolutionReader(i, LABELS.force, LABELS.pair_potential, LABELS.charge) {}

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
            nanotip.append(Atom(TYPES.SURFACE, get_point(i), get_marker(i)));
            set_marker(i, 1);
        }
        else set_marker(i, 0);

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

int ForceReader::calc_voronois(VoronoiMesh& mesh,
        const Config::Geometry& conf, const string& mesh_quality)
{
    require(interpolator, "NULL interpolator cannot be used!");
    const int n_atoms = size();
    // TODO put those values to Config, because they affect heavily how the voronoi charges will look like
    const double max_distance_from_surface = 0.5 * conf.latconst;
    const double shift_distance = 1.0 * conf.latconst;

    Medium nanotip;
    const int n_nanotip_atoms = get_nanotip(nanotip, conf.radius);

    // calculate support points for the nanotip by moving the nanotip points
    // in direction of its corresponding triangle norm by shift_distance
    Medium support(n_nanotip_atoms);
    for (int i = 0; i < n_nanotip_atoms; ++i) {
        Point3 point = nanotip.get_point(i);
        int face = abs( nanotip.get_marker(i) );

        if (interpolator->lintri.fast_distance(point, face) < max_distance_from_surface) {
            point += interpolator->lintri.get_norm(face) * shift_distance;
            support.append(Atom(TYPES.VACANCY, point, face));
        }
    }

    nanotip += support;

    // Generate Voronoi cells around the nanotip
    // r - reconstruct, v - output Voronoi cells, Q - quiet, q - mesh quality
    int err_code = mesh.generate(nanotip, conf.latconst, "rQq" + mesh_quality, "vQ");
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
            interpolation[i] = Solution(field * (charge * force_factor), 0, charge);
        }
}

void ForceReader::recalc_lorentz(const FieldReader &fields) {
    // the assumption is that the amount and IDs of field points
    // from previous full timestep have not changed
    const int n_atoms = fields.size();
    require(n_atoms == size(), "Invalid input fields!");

    for (int i = 0; i < n_atoms; ++i) {
        double charge = get_charge(i);
        Vec3 force = fields.get_elfield(i) * (charge * force_factor);
        interpolation[i] = Solution(force, 0, charge);
    }
}

void ForceReader::calc_coulomb(const double r_cut) {
    const double r_cut2 = r_cut * r_cut;
    const int n_atoms = size();

    // it is more efficinet to calculate forces from linked list than from neighbour list,
    // as in that ways it is possible to avoid double calculation of the distances
    calc_linked_list(r_cut);
    require((int)list.size() == n_atoms, "Invalid linked list size: " + d2s(list.size()));
    require((int)head.size() == nborbox_size[0]*nborbox_size[1]*nborbox_size[2],
            "Invalid linked list header size: " + d2s(head.size()));

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
                    require(i_cell >= 0 && i_cell < (int)head.size(), "Invalid neighbouring cell index: " + d2s(i_cell));

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

void ForceReader::distribute_charges(const FieldReader &fields, const ChargeReader& faces, const double smooth_factor) {
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
        interpolation[atom] = Solution(force, charges[atom]);
    }
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

int ForceReader::export_parcas(const int n_points, const string &data_type, const Vec3& si2parcas,
        double* data) const
{
    check_return(size() == 0, "No " + data_type + " to export!");

    // export total pair-potential ?
    if (data_type == LABELS.pair_potential_sum) {
        data[0] = 0;
        for (int i = 0; i < size(); ++i)
            data[0] += get_pairpot(i);
    }

    // export force perturbation in reduced units (used in Parcas)
    else if (data_type == LABELS.parcas_force) {
        for (int i = 0; i < size(); ++i) {
            int id = get_id(i);
            if (id < 0 || id >= n_points) continue;

            Vec3 reduced_force = get_force(i) * si2parcas;
            id *= 3;
            for (double f : reduced_force)
                data[id++] += f;
        }
    }

    // export charge and force in SI units ?
    else if (data_type == LABELS.charge_force) {
        for (int i = 0; i < size(); ++i) {
            int id = get_id(i);
            if (id < 0 || id >= n_points) continue;

            id *= 4;
            data[id++] = get_charge(i);
            for (double f : get_force(i))
                data[id++] = f;
        }
    }

    else
        require(false, "Unimplemented type of export data: " + data_type);

    return 0;
}

} // namespace femocs
