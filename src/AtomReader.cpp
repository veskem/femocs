/*
 * AtomReader.cpp
 *
 *  Created on: 14.12.2015
 *      Author: veske
 */

#include "AtomReader.h"

#include <cfloat>
#include <fstream>
#include <math.h>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

using namespace std;
namespace femocs {

AtomReader::AtomReader() : Medium(), rms_distance(0) {}

void AtomReader::reserve(const int n_atoms) {
    require(n_atoms >= 0, "Invalid # atoms: " + to_string(n_atoms));
    atoms.clear();

    atoms.reserve(n_atoms);
    cluster = vector<int>(n_atoms, 0);
    coordination = vector<int>(n_atoms, 0);
}

void AtomReader::extract(Surface& surface, const int type, const bool invert) {
    const int coord_min = 2;
    const int n_atoms = size();
    vector<bool> is_type(n_atoms);

    // Get number and locations of atoms of desired type
    if (!invert) {
        for (int i = 0; i < n_atoms; ++i)
            is_type[i] = get_marker(i) == type;
    } else {
        for (int i = 0; i < n_atoms; ++i)
            is_type[i] = get_marker(i) != type;
    }

    // Clean lonely atoms; atom is considered lonely if its coordination is lower than coord_min
    if (get_nborlist_size() == n_atoms)
        for (int i = 0; i < n_atoms; ++i)
            if (is_type[i]) {
                int n_nbors = 0;
                for (int nbor : get_neighbours(i)) {
                    require(nbor >= 0 && nbor < n_atoms, "Invalid index: " + d2s(nbor));
                    if (is_type[nbor]) n_nbors++;
                }

                is_type[i] = n_nbors >= coord_min;
            }

    // Store the atoms
    surface.reserve(vector_sum(is_type));
    for (int i = 0; i < n_atoms; ++i)
        if (is_type[i])
            surface.append(get_atom(i));

    surface.calc_statistics();
}

double AtomReader::calc_rms_distance(const double eps) {
    if (eps <= 0) return DBL_MAX;

    const size_t n_atoms = size();
    if (n_atoms != previous_points.size())
        return DBL_MAX;
    
    double sum = 0;
    for (size_t i = 0; i < n_atoms; ++i) {
        if (previous_types[i] != TYPES.CLUSTER && previous_types[i] != TYPES.EVAPORATED)
            sum += get_point(i).distance2(previous_points[i]);
    }

    rms_distance = sqrt(sum / n_atoms);
    return rms_distance;
}

void AtomReader::save_current_run_points(const double eps) {
    if (eps <= 0) return;
    const int n_atoms = size();

    if (n_atoms != (int)previous_points.size()) {
        previous_points.resize(n_atoms);
        previous_types.resize(n_atoms);
    }

    for (int i = 0; i < n_atoms; ++i) {
        previous_points[i] = get_point(i);
        previous_types[i] = atoms[i].marker;
    }
}

void AtomReader::calc_nborlist(const int nnn, const double r_cut, const int* parcas_nborlist) {
    require(r_cut > 0, "Invalid cut-off radius: " + to_string(r_cut));
    require(nnn > 0, "Invalid # nearest neighbours: " + to_string(nnn));

    const int n_atoms = size();
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

void AtomReader::recalc_nborlist(const double r_cut) {
    require(r_cut > 0, "Invalid cut-off radius: " + to_string(r_cut));

    const int n_atoms = size();
    const double r_cut2 = r_cut * r_cut;

    // Initialise list of new neighbourlist
    vector<vector<int>> new_nborlist(n_atoms);
    for (int i = 0; i < n_atoms; ++i)
        new_nborlist[i].reserve(nborlist[i].size());

    // Loop through all the atoms
    for (int i = 0; i < n_atoms; ++i) {
        Point3 point1 = get_point(i);
        // Loop through all the previously found neighbours of the atom
        for (int nbor : nborlist[i]) {
            if ( r_cut2 >= point1.periodic_distance2(get_point(nbor), sizes.xbox, sizes.ybox) )
                new_nborlist[i].push_back(nbor);
        }
    }

    nborlist = new_nborlist;
}

void AtomReader::calc_coordinations() {
    for (int i = 0; i < size(); ++i)
        coordination[i] = nborlist[i].size();
}

void AtomReader::calc_coordinations(int& nnn, double& latconst, double& coord_cutoff, const int* parcas_nborlist) {
    const double rdf_cutoff = 2.0 * latconst;

    calc_nborlist(nnn, rdf_cutoff, parcas_nborlist);
    calc_rdf(nnn, latconst, coord_cutoff, 200, rdf_cutoff);

    require(coord_cutoff <= rdf_cutoff, "Invalid cut-off: " + to_string(coord_cutoff));

    recalc_nborlist(coord_cutoff);
    calc_coordinations();
}

void AtomReader::calc_coordinations(const int nnn, const double coord_cutoff, const int* parcas_nborlist) {
    calc_nborlist(nnn, coord_cutoff, parcas_nborlist);
    calc_coordinations();
}

// Calculate coordination for all the atoms using neighbour list
void AtomReader::calc_coordinations(int& nnn, double& latconst, double& coord_cutoff) {
    const double rdf_cutoff = 2.0 * latconst;

    calc_verlet_nborlist(nborlist, rdf_cutoff, true);
    calc_rdf(nnn, latconst, coord_cutoff, 200, rdf_cutoff);

    require(coord_cutoff <= rdf_cutoff, "Invalid cut-off: " + to_string(coord_cutoff));

    recalc_nborlist(coord_cutoff);
    calc_coordinations();
}

// Calculate coordination for all the atoms using neighbour list
void AtomReader::calc_coordinations(const int nnn, const double r_cut) {
    calc_verlet_nborlist(nborlist, r_cut, true);
    for (int i = 0; i < size(); ++i)
        coordination[i] = nborlist[i].size();
}

// Calculate pseudo-coordination for all the atoms using the atom types
void AtomReader::calc_coordinations(const int nnn) {
    require(nnn > 0, "Invalid number of nearest neighbors!");
    const int n_atoms = size();

    for (int i = 0; i < n_atoms; ++i) {
        if (atoms[i].marker == TYPES.BULK)
            coordination[i] = nnn;
        else if (atoms[i].marker == TYPES.SURFACE)
            coordination[i] = nnn / 2;
        else if (atoms[i].marker == TYPES.VACANCY)
            coordination[i] = -1;
        else
            coordination[i] = 0;
    }
}

// Check for clustered and evaporated atoms
int AtomReader::check_clusters(const bool print) {
    const int n_detached = vector_sum(vector_not(&cluster, 0));
    if (print && n_detached > 0) {
        const int n_evaporated = vector_sum(vector_less(&cluster, 0));
        write_silent_msg("# evaporated|clustered atoms: " + d2s(n_evaporated) +
                "|" + d2s(n_detached - n_evaporated));
    }

    return n_detached;
}

// Group atoms into clusters
void AtomReader::calc_clusters() {
    const int n_atoms = size();
    require(get_nborlist_size() == n_atoms, "Clusters cannot be calculated if neighborlist is missing!");

    cluster = vector<int>(n_atoms, -1);
    vector<int> n_cluster_types;
    int c = -1;

    // for each unvisited point P in all the points
    for (int i = 0; i < n_atoms; ++i)
        if (cluster[i] < 0) {
            // mark P as visited & expand cluster
            cluster[i] = ++c;

            vector<int> neighbours = nborlist[i];

            int c_counter = 1;
            for (unsigned j = 0; j < neighbours.size(); ++j) {
                int nbor = neighbours[j];
                // if P' is not visited, connect it into the cluster
                if (cluster[nbor] < 0) {
                    c_counter++;
                    cluster[nbor] = c;
                    neighbours.insert(neighbours.end(),nborlist[nbor].begin(),nborlist[nbor].end());
                }
            }
            n_cluster_types.push_back(c_counter);
        }

    // mark clusters with one element (i.e. evaporated atoms) with minus sign
    for (int& cl : cluster) {
        if (n_cluster_types[cl] <= 1)
            cl *= -1;
    }

}

void AtomReader::calc_clusters(const double cluster_cutoff, const double coord_cutoff) {
    if (cluster_cutoff > 0 && cluster_cutoff != coord_cutoff) {
        if (cluster_cutoff < coord_cutoff)
            recalc_nborlist(cluster_cutoff);
        else
            calc_verlet_nborlist(nborlist, cluster_cutoff, true);
    }
    calc_clusters();
}

// Recalculate list of closest neighbours using Parcas neighbourlist and group atoms into clusters
void AtomReader::calc_clusters(const int nnn, const double cluster_cutoff,
        const double coord_cutoff, const int* parcas_nborlist) {

    if (cluster_cutoff > 0 && cluster_cutoff != coord_cutoff) {
        if (cluster_cutoff < coord_cutoff)
            recalc_nborlist(cluster_cutoff);
        else
            calc_nborlist(nnn, cluster_cutoff, parcas_nborlist);
    }
    calc_clusters();
}

//Calculate the radial distribution function (rdf) in a periodic isotropic system.
void AtomReader::calc_rdf(int & nnn, double& latconst, double& coord_cutoff, const int n_bins, const double r_cut) {
    require(r_cut > 0, "Invalid cut-off radius: " + to_string(r_cut));
    require(n_bins > 1, "Invalid # histogram bins: " + to_string(n_bins));

    const int n_atoms = size();
    const double bin_width = r_cut / n_bins;
    // Only valid for single-atom isotropic system
    const double norm_factor = 4.0/3.0 * M_PI * n_atoms * n_atoms / (sizes.xbox * sizes.ybox * sizes.zbox);

    // calculate the rdf histogram
    vector<double> rdf(n_bins);
    for (int i = 0; i < n_atoms; ++i) {
        Point3 point = get_point(i);
        for (int nbor : nborlist[i]) {
            const double distance2 = point.periodic_distance2(get_point(nbor), sizes.xbox, sizes.ybox);
            rdf[size_t(sqrt(distance2) / bin_width)]++;
        }
    }

    // Normalise rdf histogram by the volume and # atoms
    // Normalisation is done in a ways OVITO does
    double rdf_max = -1.0;
    for (int i = 0; i < n_bins; ++i) {
        double r = bin_width * i;
        double r2 = r + bin_width;
        rdf[i] /= norm_factor * (r2*r2*r2 - r*r*r);
        rdf_max = max(rdf_max, rdf[i]);
    }

    // Normalise rdf histogram by the maximum peak and remove noise
    for (double& r : rdf) {
        r /= rdf_max;
        if (r < 0.05) r = 0;
    }

    vector<double> peaks;
    calc_rdf_peaks(peaks, rdf, bin_width);
    require(peaks.size() >= 5, "Not enough peaks in RDF: " + to_string(peaks.size()));

    latconst = peaks[1];
    coord_cutoff = peaks[4];
    nnn = 48;
}

void AtomReader::calc_rdf_peaks(vector<double>& peaks, const vector<double>& rdf, const double bin_width) {
    const int grad_length = rdf.size() - 1;
    vector<double> gradients(grad_length);
    for (int i = 0; i < grad_length; ++i)
        gradients[i] = rdf[i+1] - rdf[i];

    peaks.clear();
    for (int i = 0; i < grad_length - 1; ++i)
        if (gradients[i] * gradients[i+1] < 0 && gradients[i] > gradients[i+1])
            peaks.push_back((i+1.5) * bin_width);
}

// Extract atom types from calculated coordinations
void AtomReader::extract_types(const int nnn, const double latconst) {
    const int n_atoms = size();
    calc_statistics();

    for (int i = 0; i < n_atoms; ++i) {
        if (cluster[i] > 0)
            atoms[i].marker = TYPES.CLUSTER;
        else if (cluster[i] < 0)
            atoms[i].marker = TYPES.EVAPORATED;
        else if (get_point(i).z < (sizes.zmin + 0.49*latconst))
            atoms[i].marker = TYPES.FIXED;
        else if (coordination[i] < nnn)
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

string AtomReader::get_data_string(const int i) const {
    if (i < 0) return "AtomReader properties=id:I:1:pos:R:3:type:I:1:coordination:I:1";

    ostringstream strs; strs << fixed;
    strs << atoms[i] << " " << coordination[i];
    return strs.str();
}

const vector<int>& AtomReader::get_neighbours(const int i) const {
    require(i >= 0 && i < size(), "Invalid index: " + to_string(i));
    require(size() == (int)nborlist.size(), "Query from invalid neighbour list!");
    return nborlist[i];
}

// Return the size of neighbour list
int AtomReader::get_nborlist_size() const {
    return nborlist.size();
}

Vec3 AtomReader::get_si2parcas_box() const {
    require(simubox.x > 0 && simubox.y > 0 && simubox.z > 0, "Invalid simubox dimensions: " + d2s(simubox));
    return Vec3(1/simubox.x, 1/simubox.y, 1/simubox.z);
}

Vec3 AtomReader::get_parcas2si_box() const {
    require(simubox.x > 0 && simubox.y > 0 && simubox.z > 0, "Invalid simubox dimensions: " + d2s(simubox));
    return simubox;
}

// =================================
// *** IMPORTERS: ***************

void AtomReader::generate_nanotip(const double h, const double radius, const double latconst) {
    const double latconst2 = latconst * latconst;
    const double tau = 2 * M_PI;
    const double box_width = 1.5*radius;
    const double height = fabs(h) * radius;

    // Over estimate the number of generated points and reserve memory for them
    const int n_nanotip = tau * (radius + latconst) * (height + radius) / latconst2;
    const int n_substrate = (4 * box_width * box_width - M_PI * pow(radius-latconst, 2)) / latconst2;
    reserve(n_nanotip + n_substrate);

    double x, y, z, r, phi, theta, dPhi, dTheta, dZ;
    unsigned int id = 0;

    dTheta = 0.5 * M_PI / round(0.5 * M_PI * radius / latconst);

    // Add the topmost atom
    append( Atom(id++, Point3(0, 0, height + radius), TYPES.SURFACE) );

    // Add the apex
    for (theta = dTheta; theta < 0.5 * M_PI; theta += dTheta) {
        z = height + radius * cos(theta);
        dPhi = tau / round(tau * radius * sin(theta) / latconst);

        for (phi = 0; phi < tau; phi += dPhi) {
            x = radius * sin(theta) * cos(phi);
            y = radius * sin(theta) * sin(phi);
            append( Atom(id++, Point3(x, y, z), TYPES.SURFACE) );
        }
    }

    dPhi = tau / round(tau * radius / latconst);
    dZ = height / round(height / latconst);
    // Add the sides of the cylinder
    for (z = height; z >= 0; z -= dZ) {
        for (phi = 0; phi < tau; phi += dPhi) {
            x = radius * cos(phi);
            y = radius * sin(phi);
            append( Atom(id++, Point3(x, y, z), TYPES.SURFACE) );
        }
    }

    // Add cylindrical substrate
    for (r = radius + latconst; r < box_width*sqrt(2); r += latconst) {
        dPhi = tau / round(tau * r / latconst);
        for (phi = 0; phi < tau; phi += dPhi) {
            x = r * cos(phi);
            y = r * sin(phi);
            if ( fabs(x) <= box_width && fabs(y) <= box_width)
                append( Atom(id++, Point3(x, y, 0), TYPES.SURFACE) );
        }
    }

    // in case of negative height, turn nanotip into opened nanovoid
    if (h < 0) {
        for (int i = 0; i < size(); ++i)
            set_z(i, -1.0*get_point(i).z);
    }

    calc_statistics();
}

void AtomReader::import_kimocs() {
    require(false, "AtomReader::import_kimocs() not implemented!");
}

void AtomReader::import_helmod(const int n_atoms, const double* x, const double* y, const double* z, const int* types) {
    require(n_atoms > 0, "Zero input atoms detected!");
    reserve(n_atoms);
    for (int i = 0; i < n_atoms; ++i)
        append( Atom(i, Point3(x[i], y[i], z[i]), types[i]) );

    calc_statistics();
}

void AtomReader::import_parcas(const int n_atoms, const double* xyz, const double* box) {
    require(n_atoms > 0, "Zero input atoms detected!");

    simubox = Vec3(box[0], box[1], box[2]);
    require(simubox.x > 0 && simubox.y > 0 && simubox.z > 0, "Invalid simubox dimensions: " + d2s(simubox));

    reserve(n_atoms);
    for (int i = 0; i < 3*n_atoms; i+=3)
        append( Atom(i/3, Point3(xyz[i+0]*box[0], xyz[i+1]*box[1], xyz[i+2]*box[2]), TYPES.BULK) );

    calc_statistics();
}

void AtomReader::import_file(const string &file_name, const bool add_noise) {
    string file_type = get_file_type(file_name);

    if (file_type == "xyz")
        import_xyz(file_name);
    else if (file_type == "ckx")
        import_ckx(file_name);
    else if (file_type == "dump")
        import_dump(file_name);
    else
        require(false, "Unsupported file type: " + file_type);

    if (add_noise) {
        // initialize random seed
        srand (time(NULL));
        // add random number that is close to the distance between nn atoms to all the coordinates
        const double eps = 0.1 * get_point(0).distance(get_point(1));
        for (int i = 0; i < size(); ++i)
            atoms[i].point += eps * static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
    }

    calc_statistics();
}

void AtomReader::import_xyz(const string &file_name) {
    ifstream in_file(file_name, ios::in);
    require(in_file.is_open(), "Did not find a file " + file_name);

    double x, y, z;
    int type, n_atoms;
    string elem, line, dummy;
    istringstream iss;

    getline(in_file, line);     // Read number of atoms
    iss.clear();
    iss.str(line);
    iss >> n_atoms;
    reserve(n_atoms);

    getline(in_file, line);     // Skip comments line

    int id = 0;
    // keep storing values from the text file as long as data exists:
    while (id < n_atoms && getline(in_file, line)) {
        iss.clear();
        iss.str(line);
        iss >> elem >> x >> y >> z >> type;
        append( Atom(id++, Point3(x, y, z), type) );
    }
}

void AtomReader::import_ckx(const string &file_name) {
    ifstream in_file(file_name, ios::in);
    require(in_file.is_open(), "Did not find a file " + file_name);

    double x, y, z;
    int type, n_atoms;
    string line, dummy;
    istringstream iss;

    getline(in_file, line); 	// Read number of atoms
    iss.clear();
    iss.str(line);
    iss >> n_atoms;

    reserve(n_atoms);

    getline(in_file, line);    // Skip comments line

    int id = 0;
    // keep storing values from the text file as long as data exists:
    while (id < n_atoms && getline(in_file, line)) {
        iss.clear();
        iss.str(line);
        iss >> type >> x >> y >> z;
        append( Atom(id++, Point3(x, y, z), type) );
    }
}

void AtomReader::import_dump(const string &) {
    require(false, "AtomReader::import_dump not implemented!");
}

} /* namespace femocs */
