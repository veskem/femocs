/*
 * AtomReader.cpp
 *
 *  Created on: 14.12.2015
 *      Author: veske
 */

#include "AtomReader.h"

#include <cfloat>
#include <fstream>
#include <sstream>

#include <omp.h>

using namespace std;
namespace femocs {

// AtomReader constructor
AtomReader::AtomReader() : Medium(), simu_type("") {}

// Reserve memory for data vectors
const void AtomReader::reserve(const int n_atoms) {
    require(n_atoms > 0, "Invalid # atoms: " + to_string(n_atoms));
    atoms.clear();
    type.clear();

    Medium::reserve(n_atoms);
    type.reserve(n_atoms);
}

// Add atom with its coordinates and type
const void AtomReader::add_atom(const int id, const Point3 &point, const int type) {
    Medium::add_atom( Atom(id, point, 0) );
    this->type.push_back(type);
}

const bool AtomReader::equals_previous_run(const double eps) {
    const int n_atoms = get_n_atoms();
    const double eps2 = eps * eps;

    if (n_atoms != previous_point.size())
        return false;

    for (int i = 0; i < n_atoms; ++i)
        if ( get_point(i).distance2(previous_point[i]) > eps )
            return false;

    return true;
}

const void AtomReader::save_current_run_points() {
    const int n_atoms = get_n_atoms();

    if (n_atoms != previous_point.size())
        previous_point.resize(n_atoms);

    for (int i = 0; i < n_atoms; ++i)
        previous_point[i] = get_point(i);
}

const void AtomReader::calc_coordination(const int nnn, const double cutoff, const int* nborlist) {
    require(cutoff > 0, "Invalid cutoff: " + to_string(cutoff));
    require(nnn >= 0, "Invalid # nearest neighbors: " + to_string(nnn));

    const int n_atoms = get_n_atoms();
    const double cutoff2 = cutoff * cutoff;
    int nbor_indx;

    nbor_indx = 0;
    for (int i = 0; i < n_atoms; ++i) {
        Point3 point1 = get_point(i);

        // Loop through neighbours in neighbour list
        int n_nbors = nborlist[nbor_indx++];
        for (int j = 0; j < n_nbors; ++j) {
            int nbor = nborlist[nbor_indx + j] - 1;
            require (nbor >= 0 && nbor < n_atoms, "Invalid neighbour: " + to_string(nbor));

            double r2 = point1.periodic_distance2(get_point(nbor), sizes.xbox, sizes.ybox);

            if (r2 <= cutoff2) {
                atoms[i].coord++;
                atoms[nbor].coord++;
//                if(atoms[i].coord < nnn) atoms[i].coord++;
//                if(atoms[nbor].coord < nnn) atoms[nbor].coord++;
            }
        }
        nbor_indx += n_nbors;
    }
    check_coordination();
}

// Calculate coordination (# nearest neighbours within cutoff radius) for all the atoms
const void AtomReader::calc_coordination(const int nnn, const double cutoff) {
    require(nnn > 0, "Invalid number of nearest neighbors: " + to_string(nnn));
    require(simu_type == "md" || simu_type == "kmc", "Incorrect simulation type: " + simu_type);

    if (simu_type == "md")
        calc_slow_coordination(cutoff);
    else if (simu_type == "kmc")
        calc_dummy_coordination(nnn);

    check_coordination();
}

// Check for evaporated atoms
const void AtomReader::check_coordination() {
    const int coord_max = 1;
    const int coord_min = 0;
    for (Atom a : atoms)
        if ((a.coord < coord_max) && (a.coord >= coord_min))
            expect(false, "WARNING: evaporated atom detected!");
}

// Calculate without neighbour list the coordination for atoms coming from MD simulations
const void AtomReader::calc_slow_coordination(const double cutoff) {
    require(cutoff > 0, "Invalid cutoff: " + to_string(cutoff));

    expect(false, "WARNING: Running slow coordination calculation!");

    const int n_atoms = get_n_atoms();
    const double cutoff2 = cutoff * cutoff;

    vector<int> coordinations(n_atoms, 0);

    for (int i = 0; i < n_atoms - 1; ++i) {
        Point3 point1 = get_point(i);
        for (int j = i + 1; j < n_atoms; ++j) {
            double r2 = point1.periodic_distance2(get_point(j), sizes.xbox, sizes.ybox);
            if (r2 <= cutoff2) {
                coordinations[i]++;
                coordinations[j]++;
            }
        }
    }

    for (int i = 0; i < n_atoms; ++i)
        set_coordination(i, coordinations[i]);
}

// Calculate coordination for atoms coming from KMC simulations
const void AtomReader::calc_dummy_coordination(const int nnn) {
    require(nnn > 0, "Invalid number of nearest neighbors!");
    const int n_atoms = get_n_atoms();

    for (int i = 0; i < n_atoms; ++i) {
        if (type[i] == TYPES.BULK)
            set_coordination(i, nnn);
        else if (type[i] == TYPES.SURFACE)
            set_coordination(i, (int) nnn / 2);
        else
            set_coordination(i, 0);
    }
}

const void AtomReader::extract_types(const int nnn, const double latconst) {
    const int n_atoms = get_n_atoms();
    calc_statistics();
    crys_struct.latconst = latconst;

    for (int i = 0; i < n_atoms; ++i) {
        if( get_coordination(i) >= (nnn - 1) )
            type[i] = TYPES.BULK;
        else if(get_point(i).z < (sizes.zmin + 2.0 * crys_struct.latconst))
            type[i] = TYPES.FIXED;
        else
            type[i] = TYPES.SURFACE;
    }
}

// Redefine simubox size for example to insert electric field height
const void AtomReader::resize_box(const double zmin, const double zmax) {
    require(zmin <= zmax, "Invalid size for simulation box: " + to_string(zmin) + ", " + to_string(zmax));
    sizes.zminbox = zmin;
    sizes.zmaxbox = zmax;
    sizes.zbox = zmax - zmin;
}

// =================================
// *** GETTERS: ***************

const int AtomReader::get_type(const int i) {
    require(i >= 0 && i < get_n_atoms(), "Invalid index!");
    return type[i];
}

// Compile data string from the data vectors
const string AtomReader::get_data_string(const int i) {
    if (i < 0) return "AtomReader data: id x y z coordination type";

    ostringstream strs;
    strs << atoms[i] << " " << get_type(i);
    return strs.str();
}

// =================================
// *** IMPORTERS: ***************

const void AtomReader::import_kimocs() {
    require(false, "AtomReader::import_kimocs() not implemented!");
    simu_type = "kmc";
}

const void AtomReader::import_helmod(int n_atoms, double* x, double* y, double* z, int* types) {
    require(n_atoms > 0, "Zero input atoms detected!");
    simu_type = "md";
    reserve(n_atoms);
    for (int i = 0; i < n_atoms; ++i)
        add_atom(i, Point3(x[i], y[i], z[i]), types[i]);

    calc_statistics();
}

const void AtomReader::import_parcas(int n_atoms, const double* xyz, const double* box) {
    require(n_atoms > 0, "Zero input atoms detected!");
    simu_type = "md";
    reserve(n_atoms);
    for (int i = 0; i < 3*n_atoms; i+=3)
        add_atom(i/3, Point3(xyz[i+0]*box[0], xyz[i+1]*box[1], xyz[i+2]*box[2]), TYPES.BULK);

    calc_statistics();
}

const void AtomReader::import_file(const string file_name) {
    string file_type = get_file_type(file_name);

    if (file_type == "xyz") {
        simu_type = "md";
        import_xyz(file_name);
    }

    else if (file_type == "ckx") {
        simu_type = "kmc";
        import_ckx(file_name);
    }

    else if (file_type == "dump") {
        simu_type = "md";
        import_dump(file_name);
    }

    else
        require(false, "Unsupported file type: " + file_type);

    calc_statistics();
}

const void AtomReader::import_xyz(const string file_name) {
    ifstream in_file(file_name, ios::in);
    require(in_file.is_open(), "Did not find a file " + file_name);

    double x, y, z;
    int type, id;
    string elem, line, dummy;
    size_t n_atoms;
    istringstream iss;

    getline(in_file, line);     // Read number of atoms
    iss.clear();
    iss.str(line);
    iss >> n_atoms;
    reserve(n_atoms);

    getline(in_file, line);     // Skip comments line

    id = 0;
    // keep storing values from the text file as long as data exists:
    while ((--n_atoms > 0) && getline(in_file, line)) {
        iss.clear();
        iss.str(line);
        iss >> elem >> x >> y >> z >> type;
        add_atom(id++, Point3(x, y, z), type);
    }
}

const void AtomReader::import_ckx(const string file_name) {
    ifstream in_file(file_name, ios::in);
    require(in_file.is_open(), "Did not find a file " + file_name);

    double x, y, z;
    int type, id;
    string line, dummy;
    size_t n_atoms;
    istringstream iss;

    getline(in_file, line); 	// Read number of atoms
    iss.clear();
    iss.str(line);
    iss >> n_atoms;

    reserve(n_atoms);

    getline(in_file, line);    // Skip comments line

    id = 0;
    // keep storing values from the text file as long as data exists:
    while ((--n_atoms > 0) && getline(in_file, line)) {
        iss.clear();
        iss.str(line);
        iss >> type >> x >> y >> z;
        add_atom(id++, Point3(x, y, z), type);
    }
}

const void AtomReader::import_dump(const string file_name) {
    require(false, "AtomReader::import_dump not implemented!");
}

} /* namespace femocs */
