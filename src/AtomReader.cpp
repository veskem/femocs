/*
 * AtomReader.cpp
 *
 *  Created on: 14.12.2015
 *      Author: veske
 */

#include "AtomReader.h"

#include <cfloat>
#include <fstream>

using namespace std;
namespace femocs {

// AtomReader constructor
AtomReader::AtomReader() : Medium(), rms_distance(0) {}

// Reserve memory for data vectors
void AtomReader::reserve(const int n_atoms) {
    require(n_atoms >= 0, "Invalid # atoms: " + to_string(n_atoms));
    atoms.clear();

    atoms.reserve(n_atoms);
    cluster = vector<int>(n_atoms, 0);
    coordination = vector<int>(n_atoms, 0);
}

bool AtomReader::equals_previous_run(const double eps) {
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
double AtomReader::calc_rms_distance(const double eps) {
    if (eps < 1e-5) return DBL_MAX;

    const int n_atoms = get_n_atoms();
    if (n_atoms != previous_point.size())
        return DBL_MAX;
    
    double sum = 0;
    for (int i = 0; i < n_atoms; ++i)
        sum += get_point(i).distance2(previous_point[i]);
    
    rms_distance = sqrt(sum / n_atoms);
    return rms_distance;
}

void AtomReader::save_current_run_points(const double eps) {
    if (eps < 1e-5) return;
    const int n_atoms = get_n_atoms();

    if (n_atoms != previous_point.size())
        previous_point.resize(n_atoms);

    for (int i = 0; i < n_atoms; ++i)
        previous_point[i] = get_point(i);
}

// Calculate list of close neighbours using Parcas diagonal neighbour list
void AtomReader::calc_nborlist(const int nnn, const double r_cut, const int* parcas_nborlist) {
    require(r_cut > 0, "Invalid cut-off radius: " + to_string(r_cut));

    const int n_atoms = get_n_atoms();
    const double r_cut2 = r_cut * r_cut;

    // Initialise list of closest neighbours
    nborlist = vector<vector<int>>(n_atoms);
    for (int i = 0; i < n_atoms; ++i)
        nborlist[i].reserve(nnn);

    // Loop through all the atoms
    int nbor_indx = 0;
    for (int i = 0; i < n_atoms; ++i) {
        Point3 point1 = get_point(i);

        int n_nbors = parcas_nborlist[nbor_indx++];

        // Loop through atom neighbours
        for (int j = 0; j < n_nbors; ++j) {
            int nbr = parcas_nborlist[nbor_indx++] - 1;
            if ( r_cut2 >= point1.periodic_distance2(get_point(nbr), sizes.xbox, sizes.ybox) ) {
                nborlist[i].push_back(nbr);
                nborlist[nbr].push_back(i);
            }
        }
    }
}

// Calculate list of close neighbours using brute force technique
void AtomReader::calc_nborlist(const int nnn, const double r_cut) {
    require(r_cut > 0, "Invalid cut-off radius: " + to_string(r_cut));

    const int n_atoms = get_n_atoms();
    const double r_cut2 = r_cut * r_cut;

    // Initialise list of closest neighbours
    nborlist = vector<vector<int>>(n_atoms);
    for (int i = 0; i < n_atoms; ++i)
        nborlist[i].reserve(nnn);

    int percentage = 0;
    if (MODES.VERBOSE) cout << endl;

    // Loop through all the atoms
    for (int i = 0; i < n_atoms - 1; ++i) {
        // show the progress (for big systems the calculation may take a lot of time)
        if (MODES.VERBOSE && i % (n_atoms / 44) == 0) {
            cout << "\r" << string(percentage, '*') << string(44-percentage, ' ') << " " << flush;
            percentage++;
        }

        Point3 point1 = get_point(i);
        // Loop through all the possible neighbours of the atom
        for (int j = i + 1; j < n_atoms; ++j) {
            if ( r_cut2 >= point1.periodic_distance2(get_point(j), sizes.xbox, sizes.ybox) ) {
                nborlist[i].push_back(j);
                nborlist[j].push_back(i);
            }
        }
    }
}

// Check for clustered atoms
void AtomReader::check_clusters() {
    const int n_clusters = vector_sum(vector_not(&cluster, 0));
    if (n_clusters > 0)
        cout << "FEMOCS: # evaporated or clustered atoms: " + to_string(n_clusters) << endl;
}

// Group atoms into clusters
void AtomReader::calc_clusters() {
    const int n_atoms = get_n_atoms();
    vector<int> neighbours;
    vector<bool> not_marked(n_atoms, true);
    int c = -1;

    // for each unvisited point P in all the points
    for (int i = 0; i < n_atoms; ++i) {
        if (not_marked[i]) {
            // Mark P as visited & expand cluster
            not_marked[i] = false;
            cluster[i] = ++c;

            neighbours = nborlist[i];

            for (int j = 0; j < neighbours.size(); ++j) {
                int nbor = neighbours[j];
                // if P' is not visited
                if (not_marked[nbor]) {
                    // Mark P' as visited & expand cluster
                    not_marked[nbor] = false;
                    cluster[nbor] = c;
                    neighbours.insert(neighbours.end(),nborlist[nbor].begin(),nborlist[nbor].end());
                }
            }
        }
    }
}

// Check for evaporated atoms
void AtomReader::check_coordinations() {
    const int coord_max = 1;
    const int coord_min = 0;
    for (int a : coordination)
        if (a >= coord_min && a <= coord_max)
            cout << "FEMOCS: Evaporated atom detected!" << endl;
}

// Calculate coordination for all the atoms
void AtomReader::calc_coordinations() {
    for (int i = 0; i < get_n_atoms(); ++i)
        coordination[i] = nborlist[i].size();
}

// Calculate pseudo-coordination for all the atoms using the atom types
void AtomReader::calc_coordinations(const int nnn) {
    require(nnn > 0, "Invalid number of nearest neighbors!");
    const int n_atoms = get_n_atoms();

    for (int i = 0; i < n_atoms; ++i) {
        if (atoms[i].marker == TYPES.BULK)
            coordination[i] = nnn;
        else if (atoms[i].marker == TYPES.SURFACE)
            coordination[i] = nnn/2;
        else if (atoms[i].marker == TYPES.VACANCY)
            coordination[i] = -1;
        else
            coordination[i] = 0;
    }
}

// Extract atom types from calculated coordinations
void AtomReader::extract_types(const int nnn, const double latconst) {
    const int n_atoms = get_n_atoms();
    const int nnn_eps = 3;
    calc_statistics();

    for (int i = 0; i < n_atoms; ++i) {
        if (cluster[i] != 0)
            atoms[i].marker = TYPES.CLUSTER;
        else if (get_point(i).z < (sizes.zmin + 0.49*latconst))
            atoms[i].marker = TYPES.FIXED;
        else if ( (nnn - coordination[i]) >= nnn_eps )
            atoms[i].marker = TYPES.SURFACE;
        else
            atoms[i].marker = TYPES.BULK;
    }
}

// Redefine simubox size for example to insert electric field height
void AtomReader::resize_box(const double zmin, const double zmax) {
    require(zmin <= zmax, "Invalid size for simulation box: " + to_string(zmin) + ", " + to_string(zmax));
    sizes.zminbox = zmin;
    sizes.zmaxbox = zmax;
    sizes.zbox = zmax - zmin;
}

// Redefine the min and max values for x, y and z - coordinates
void AtomReader::resize_box(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax) {
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

// Compile data string from the data vectors
string AtomReader::get_data_string(const int i) const {
    if (i < 0) return "AtomReader properties=id:I:1:pos:R:3:type:I:1:coordination:I:1:cluster:I:1";

    ostringstream strs; strs << fixed;
    strs << atoms[i] << " " << coordination[i] << " " << cluster[i];
    return strs.str();
}

// Get the closest neighbours of i-th atom
const vector<int>& AtomReader::get_neighbours(const int i) const {
    require(i >= 0 && i < get_n_atoms(), "Invalid index: " + to_string(i));
    require(nborlist.size() == get_n_atoms(), "Query from invalid neighbour list!");
    return nborlist[i];
}

// =================================
// *** IMPORTERS: ***************

void AtomReader::import_kimocs() {
    require(false, "AtomReader::import_kimocs() not implemented!");
}

void AtomReader::import_helmod(const int n_atoms, const double* x, const double* y, const double* z, const int* types) {
    require(n_atoms > 0, "Zero input atoms detected!");
    reserve(n_atoms);
    for (int i = 0; i < n_atoms; ++i)
        add_atom( Atom(i, Point3(x[i], y[i], z[i]), types[i]) );

    calc_statistics();
}

void AtomReader::import_parcas(const int n_atoms, const double* xyz, const double* box) {
    require(n_atoms > 0, "Zero input atoms detected!");
    reserve(n_atoms);
    for (int i = 0; i < 3*n_atoms; i+=3)
        add_atom( Atom(i/3, Point3(xyz[i+0]*box[0], xyz[i+1]*box[1], xyz[i+2]*box[2]), TYPES.BULK) );

    calc_statistics();
}

void AtomReader::import_file(const string &file_name) {
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

void AtomReader::import_xyz(const string &file_name) {
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
        add_atom( Atom(id++, Point3(x, y, z), type) );
    }
}

void AtomReader::import_ckx(const string &file_name) {
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
        add_atom( Atom(id++, Point3(x, y, z), type) );
    }
}

void AtomReader::import_dump(const string &file_name) {
    require(false, "AtomReader::import_dump not implemented!");
}

} /* namespace femocs */
