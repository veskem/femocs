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
AtomReader::AtomReader() : Medium() {}

// Reserve memory for data vectors
const void AtomReader::reserve(const int n_atoms) {
    require(n_atoms > 0, "Invalid # atoms: " + to_string(n_atoms));
    atoms.clear();
    types.clear();

    atoms.reserve(n_atoms);
    types.reserve(n_atoms);
}

// Add atom with its coordinates and type
const void AtomReader::add_atom(const int id, const Point3 &point, const int type) {
    Medium::add_atom( Atom(id, point, 0) );
    this->types.push_back(type);
}

const bool AtomReader::equals_previous_run(const double eps) {
    if (eps < 1e-5)
        return false;

    const int n_atoms = get_n_atoms();
    const double eps2 = eps * eps;

    if (n_atoms != previous_point.size())
        return false;

     for (int i = 0; i < n_atoms; ++i)
         if ( get_point(i).distance2(previous_point[i]) > eps2 )
             return false;
     return true;
}

// Calculate root mean square of the distances atoms have moved after previous run
const double AtomReader::get_rms_distance(const double eps) {
    if (eps < 1e-5) return DBL_MAX;

    const int n_atoms = get_n_atoms();
    if (n_atoms != previous_point.size())
        return DBL_MAX;
    
    double sum = 0;
    for (int i = 0; i < n_atoms; ++i)
        sum += get_point(i).distance2(previous_point[i]);
    
    return sqrt(sum / n_atoms);
}

const void AtomReader::save_current_run_points(const double eps) {
    if (eps < 1e-5) return;
    const int n_atoms = get_n_atoms();

    if (n_atoms != previous_point.size())
        previous_point.resize(n_atoms);

    for (int i = 0; i < n_atoms; ++i)
        previous_point[i] = get_point(i);
}

// Check for evaporated atoms
const void AtomReader::check_coordination() {
    const int coord_max = 1;
    const int coord_min = 0;
    for (Atom a : atoms)
        if (a.coord >= coord_min && a.coord <= coord_max)
            expect(false, "Evaporated atom detected!");
}

// Calculate coordination for all the atoms
const void AtomReader::calc_coordination(const int nnn, const double cutoff, const int* nborlist) {
    require(cutoff > 0, "Invalid cutoff: " + to_string(cutoff));
    require(nnn >= 0, "Invalid # nearest neighbors: " + to_string(nnn));

    const int n_atoms = get_n_atoms();
    const double cutoff2 = cutoff * cutoff;

    int nbor_indx = 0;
    // Loop through all the atoms
    for (int i = 0; i < n_atoms; ++i) {
        Point3 point1 = get_point(i);

        // Loop through atom neighbours
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

// Calculate coordination for all the atoms using brute force technique
const void AtomReader::calc_coordination(const double cutoff) {
    require(cutoff > 0, "Invalid cutoff: " + to_string(cutoff));
    expect(false, "Running slow coordination calculation!");
    
    const int n_atoms = get_n_atoms();
    const double cutoff2 = cutoff * cutoff;
    vector<int> coordinations(n_atoms, 0);

    int percentage = 0;
    // Loop through all the atoms
    for (int i = 0; i < n_atoms - 1; ++i) {
        // show progress (for big systems the calculation may take a lot of time)
        if (MODES.VERBOSE && i % (n_atoms / 44) == 0) {
            cout << "\r" << string(percentage, '*') << string(44-percentage, ' ') << " " << flush;
            percentage++;
        }
        
        Point3 point1 = get_point(i);
        // Loop through all the possible neighbours of the atom
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
    check_coordination();
}

// Calculate coordination for all the atoms using the atom types
const void AtomReader::calc_coordination(const int nnn) {
    require(nnn > 0, "Invalid number of nearest neighbors!");
    const int n_atoms = get_n_atoms();

    for (int i = 0; i < n_atoms; ++i) {
        if (types[i] == TYPES.BULK)
            set_coordination(i, nnn);
        else if (types[i] == TYPES.SURFACE)
            set_coordination(i, (int) nnn / 2);
        else if (types[i] == TYPES.VACANCY)
            set_coordination(i, -1);
        else
            set_coordination(i, 0);
    }

    check_coordination();
}

// Extract atom types from calculated coordinations
const void AtomReader::extract_types(const int nnn, const double latconst) {
    const int n_atoms = get_n_atoms();
    const int nnn_eps = 1;
    calc_statistics();

    for (int i = 0; i < n_atoms; ++i) {
        if ( (nnn - get_coordination(i)) <= nnn_eps )
            types[i] = TYPES.BULK;
        else if (get_point(i).z < (sizes.zmin + latconst))
            types[i] = TYPES.FIXED;
        else
            types[i] = TYPES.SURFACE;
    }
}

// Redefine simubox size for example to insert electric field height
const void AtomReader::resize_box(const double zmin, const double zmax) {
    require(zmin <= zmax, "Invalid size for simulation box: " + to_string(zmin) + ", " + to_string(zmax));
    sizes.zminbox = zmin;
    sizes.zmaxbox = zmax;
    sizes.zbox = zmax - zmin;
}

// Redefine the min and max values for x, y and z - coordinates
const void AtomReader::resize_box(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax) {
    require(xmin <= xmax, "Invalid x-size for simulation box: " + to_string(xmin) + ", " + to_string(xmax));
    require(ymin <= ymax, "Invalid y-size for simulation box: " + to_string(ymin) + ", " + to_string(ymax));
    require(zmin <= zmax, "Invalid z-size for simulation box: " + to_string(zmin) + ", " + to_string(zmax));
    
    sizes.xmin = xmin; sizes.xmax = xmax;
    sizes.ymin = ymin; sizes.ymax = ymax;
    sizes.zmin = zmin; sizes.zmax = zmax;
    
    // Define size of simubox
    sizes.xbox = xmax - xmin;
    sizes.ybox = ymax - ymin;
    sizes.zbox = zmax - zmin;
    sizes.zminbox = zmin;
    sizes.zmaxbox = zmax;

    // Define the centre of simubox
    sizes.xmid = (xmax + xmin) / 2;
    sizes.ymid = (ymax + ymin) / 2;
    sizes.zmid = (zmax + zmin) / 2;
}

// =================================
// *** GETTERS: ***************

const int AtomReader::get_type(const int i) const {
    require(i >= 0 && i < get_n_atoms(), "Invalid index!");
    return types[i];
}

// Compile data string from the data vectors
const string AtomReader::get_data_string(const int i) const {
    if (i < 0) return "AtomReader properties=id:R:1:pos:R:3:coordination:R:1:type:R:1";

    ostringstream strs;
    strs << atoms[i] << " " << get_type(i);
    return strs.str();
}

// =================================
// *** IMPORTERS: ***************

const void AtomReader::import_kimocs() {
    require(false, "AtomReader::import_kimocs() not implemented!");
}

const void AtomReader::import_helmod(int n_atoms, double* x, double* y, double* z, int* types) {
    require(n_atoms > 0, "Zero input atoms detected!");
    reserve(n_atoms);
    for (int i = 0; i < n_atoms; ++i)
        add_atom(i, Point3(x[i], y[i], z[i]), types[i]);

    calc_statistics();
}

const void AtomReader::import_parcas(int n_atoms, const double* xyz, const double* box) {
    require(n_atoms > 0, "Zero input atoms detected!");
    reserve(n_atoms);
    for (int i = 0; i < 3*n_atoms; i+=3)
        add_atom(i/3, Point3(xyz[i+0]*box[0], xyz[i+1]*box[1], xyz[i+2]*box[2]), TYPES.BULK);

    calc_statistics();
}

const void AtomReader::import_file(const string file_name) {
    string file_type = get_file_type(file_name);

    if (file_type == "xyz")
        import_xyz(file_name);
    else if (file_type == "ckx")
        import_ckx(file_name);
    else if (file_type == "dump")
        import_dump(file_name);
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
