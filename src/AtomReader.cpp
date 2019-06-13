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

AtomReader::AtomReader() : Medium(), conf(NULL) {}

AtomReader::AtomReader(const Config::Geometry *c) : Medium(), conf(c)
{}

void AtomReader::reserve(const int n_atoms) {
    require(n_atoms >= 0, "Invalid # atoms: " + to_string(n_atoms));
    require(conf, "NULL conf can't be used!");
    atoms.clear();

    atoms.reserve(n_atoms);
    cluster = vector<int>(n_atoms, 0);
    coordination = vector<int>(n_atoms, 0);
}

void AtomReader::extract(Surface& surface, const int type, const bool invert) {
    const unsigned int coord_min = 2;
    const unsigned int n_atoms = size();
    vector<bool> is_type(n_atoms);

    // Get number and locations of atoms of desired type
    if (!invert) {
        for (unsigned int i = 0; i < n_atoms; ++i)
            is_type[i] = get_marker(i) == type;
    } else {
        for (unsigned int i = 0; i < n_atoms; ++i)
            is_type[i] = get_marker(i) != type;
    }

    // Clean lonely atoms; atom is considered lonely if its coordination is lower than coord_min
    if (nborlist.size() == n_atoms)
        for (unsigned i = 0; i < n_atoms; ++i)
            if (is_type[i]) {
                unsigned int n_nbors = 0;
                for (int nbor : nborlist[i]) {
                    require(nbor >= 0 && nbor < (int)n_atoms, "Invalid index: " + d2s(nbor));
                    if (is_type[nbor]) n_nbors++;
                }

                is_type[i] = n_nbors >= coord_min;
            }

    // Store the atoms
    surface.reserve(vector_sum(is_type));
    for (unsigned i = 0; i < n_atoms; ++i)
        if (is_type[i])
            surface.append(get_atom(i));

    surface.calc_statistics();
}

bool AtomReader::calc_rms_distance() {
    data.rms_distance = DBL_MAX;

    const size_t n_atoms = size();
    if (n_atoms != previous_points.size())
        return true;
    
    double sum = 0;
    for (size_t i = 0; i < n_atoms; ++i) {
        if (previous_types[i] != TYPES.CLUSTER &&
                previous_types[i] != TYPES.EVAPORATED &&
                previous_types[i] != TYPES.FIXED)
            sum += get_point(i).distance2(previous_points[i]);
    }

    data.rms_distance = sqrt(sum / n_atoms);
    return data.rms_distance >= conf->distance_tol;
}

void AtomReader::save_current_run_points() {
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

void AtomReader::calc_nborlist(const double r_cut, const int* parcas_nborlist) {
    require(r_cut > 0, "Invalid cut-off radius: " + to_string(r_cut));
    require(conf->nnn > 0, "Invalid # nearest neighbours: " + to_string(conf->nnn));

    const int n_atoms = size();
    const double r_cut2 = r_cut * r_cut;

    // Initialise list of closest neighbours
    nborlist = vector<vector<int>>(n_atoms);
    for (int i = 0; i < n_atoms; ++i)
        nborlist[i].reserve(conf->nnn);

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

void AtomReader::calc_rdf_coordinations(const int* parcas_nborlist) {
    if (data.latconst <= 0) {
        data.coord_cutoff = conf->coordination_cutoff;
        data.latconst = conf->latconst;
    }

    const double rdf_cutoff = 2.0 * data.latconst;

    if (parcas_nborlist)
        calc_nborlist(rdf_cutoff, parcas_nborlist);
    else
        calc_verlet_nborlist(nborlist, rdf_cutoff, true);

    calc_rdf(200, rdf_cutoff);
    require(data.coord_cutoff <= rdf_cutoff, "Invalid cut-off: " + to_string(data.coord_cutoff));

    recalc_nborlist(data.coord_cutoff);
    for (int i = 0; i < size(); ++i)
        coordination[i] = nborlist[i].size();
}

void AtomReader::calc_coordinations(const int* parcas_nborlist) {
    data.coord_cutoff = conf->coordination_cutoff;
    data.latconst = conf->latconst;

    if (parcas_nborlist)
        calc_nborlist(conf->coordination_cutoff, parcas_nborlist);
    else
        calc_verlet_nborlist(nborlist, conf->coordination_cutoff, true);

    for (int i = 0; i < size(); ++i)
        coordination[i] = nborlist[i].size();
}

void AtomReader::calc_pseudo_coordinations() {
    require(conf->nnn > 0, "Invalid # nearest neighbours: " + to_string(conf->nnn));
    const int n_atoms = size();

    for (int i = 0; i < n_atoms; ++i) {
        if (atoms[i].marker == TYPES.BULK)
            coordination[i] = conf->nnn;
        else if (atoms[i].marker == TYPES.SURFACE)
            coordination[i] = conf->nnn / 2;
        else if (atoms[i].marker == TYPES.VACANCY)
            coordination[i] = -1;
        else
            coordination[i] = 0;
    }
}

void AtomReader::calc_clusters(const int* parcas_nborlist) {
    // if needed, update neighbor list
    if (conf->cluster_cutoff > 0 && conf->cluster_cutoff != data.coord_cutoff) {
        if (conf->cluster_cutoff < data.coord_cutoff)
            recalc_nborlist(conf->cluster_cutoff);
        else if (parcas_nborlist)
            calc_nborlist(conf->cluster_cutoff, parcas_nborlist);
        else
            calc_verlet_nborlist(nborlist, conf->cluster_cutoff, true);
    }

    const unsigned int n_atoms = size();
    require(nborlist.size() == n_atoms, "Clusters cannot be calculated if neighborlist is missing!");

    // group atoms into clusters, i.e perform cluster analysis

    cluster = vector<int>(n_atoms, -1);
    vector<int> n_cluster_types;
    int c = -1;

    // for each unvisited point P in all the points
    for (unsigned int i = 0; i < n_atoms; ++i)
        if (cluster[i] < 0) {
            // mark P as visited & expand cluster
            cluster[i] = ++c;

            vector<int> neighbours = nborlist[i];

            int c_counter = 1;
            for (unsigned int j = 0; j < neighbours.size(); ++j) {
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

    // calculate statistics about clustered atoms
    data.n_detached = vector_sum(vector_not(&cluster, 0));
    data.n_evaporated = data.n_detached - vector_sum(vector_less(&cluster, 0));
}

void AtomReader::calc_rdf(const int n_bins, const double r_cut) {
    require(r_cut > 0, "Invalid cut-off radius: " + to_string(r_cut));
    require(n_bins > 1, "Invalid # histogram bins: " + to_string(n_bins));

    const int n_atoms = size();
    const double bin_width = r_cut / n_bins;
    // factor to normalize RDF with respect to ideal gas
    const double norm_factor = 4.0/3.0 * M_PI * n_atoms * n_atoms / (sizes.xbox * sizes.ybox * sizes.zbox);

    // calculate the rdf histogram
    vector<double> rdf(n_bins);
    for (int i = 0; i < n_atoms; ++i)
        if (get_marker(i) != TYPES.FIXED) {
            Point3 point = get_point(i);
            for (int nbor : nborlist[i]) {
                const double distance2 = point.periodic_distance2(get_point(nbor), sizes.xbox, sizes.ybox);
                rdf[size_t(sqrt(distance2) / bin_width)]++;
            }
        }

    // Normalise rdf histogram by with respect to ideal gas
    // Also find the location of first neighbouring cell
    // that is located on the first peak of rdf
    double rdf_max = -1.0;
    int sigma_indx = 0;
    for (int i = 0; i < n_bins; ++i) {
        double r1 = bin_width * i;
        double r2 = r1 + bin_width;
        rdf[i] /= norm_factor * (r2*r2*r2 - r1*r1*r1);
        if (rdf_max < rdf[i]) {
            rdf_max = rdf[i];
            sigma_indx = i;
        }
    }

    double sigma = sigma_indx * bin_width;
    data.latconst = sqrt(2) * sigma;
    data.coord_cutoff = sqrt(5) * sigma;
}

void AtomReader::extract_types() {
    const int n_atoms = size();
    calc_statistics();

    for (int i = 0; i < n_atoms; ++i) {
        if (cluster[i] > 0)
            atoms[i].marker = TYPES.CLUSTER;
        else if (cluster[i] < 0)
            atoms[i].marker = TYPES.EVAPORATED;
        else if (get_point(i).z < (sizes.zmin + 0.49*conf->latconst))
            atoms[i].marker = TYPES.FIXED;
        else if (coordination[i] < conf->nnn)
            atoms[i].marker = TYPES.SURFACE;
        else
            atoms[i].marker = TYPES.BULK;
    }
}

string AtomReader::get_data_string(const int i) const {
    if (i < 0) return "AtomReader properties=id:I:1:pos:R:3:type:I:1:coordination:I:1";

    ostringstream strs; strs << fixed;
    strs << atoms[i] << " " << coordination[i];
    return strs.str();
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

void AtomReader::generate_nanotip(double h, double radius, double latconst) {
    radius -= 0.05 * latconst; // make actual radius a bit smaller to prevent problems during coarsening
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

bool AtomReader::import_atoms(const int n_atoms, const double* x, const double* y, const double* z, const int* types) {
    require(n_atoms > 0, "Zero input atoms detected!");

    atoms.clear();
    atoms.reserve(n_atoms);
    for (int i = 0; i < n_atoms; ++i)
        append( Atom(i, Point3(x[i], y[i], z[i]), types[i]) );

    calc_statistics();

    if (calc_rms_distance()) {
        cluster = vector<int>(n_atoms, 0);
        coordination = vector<int>(n_atoms, 0);
        return true;
    }
    return false;
}

bool AtomReader::import_parcas(const int n_atoms, const double* xyz,
        const double* vel, const double* box, const Config& conf)
{
    require(n_atoms > 0, "Zero input atoms detected!");
    simubox = Vec3(box[0], box[1], box[2]);
    require(simubox.x > 0 && simubox.y > 0 && simubox.z > 0, "Invalid simubox dimensions: " + d2s(simubox));
    require(conf.behaviour.timestep_fs > 0, "Invalid MD time step: " + d2s(conf.behaviour.timestep_fs));
    Vec3 parcas2si = simubox / conf.behaviour.timestep_fs;

    atoms.resize(n_atoms);
    velocities.resize(n_atoms);
    for (int i = 0; i < n_atoms; ++i) {
        int I = 3*i;
        atoms[i] = Atom(i, Point3(xyz[I], xyz[I+1], xyz[I+2]) * simubox, TYPES.BULK);
        velocities[i] = Vec3(vel[I], vel[I+1], vel[I+2]) * parcas2si;
    }

    calc_statistics();

    if (calc_rms_distance()) {
        cluster = vector<int>(n_atoms, 0);
        coordination = vector<int>(n_atoms, 0);
        return true;
    }
    return false;
}

bool AtomReader::import_lammps(const int n_atoms, const double* const* xyz,
        const double* const* vel, const int* mask, const int groupbit)
{
    require(n_atoms > 0, "Zero input atoms detected!");

    // reserve memory for stored atoms
    int n_store = 0;
    for (int i = 0; i < n_atoms; ++i)
        if (mask[i] & groupbit)
            n_store++;
    atoms.resize(n_store);
    velocities.resize(n_store);

    // TODO In case of groupbit != all, ID value might exceed n_store.
    // Figure out how to store the atoms properly or consider ignoring groupbit.
    n_store = 0;
    for (int i = 0; i < n_atoms; ++i) {
        if (mask[i] & groupbit) {
            atoms[n_store] = Atom(i, Point3(xyz[i][0], xyz[i][1], xyz[i][2]), TYPES.BULK);
            velocities[n_store] = Vec3(vel[i][0], vel[i][1], vel[i][2]) * 1e-3;  // Lammps velocities are in Angstrom / ps
            n_store++;
        }
    }

    calc_statistics();

    if (calc_rms_distance()) {
        cluster = vector<int>(n_store, 0);
        coordination = vector<int>(n_store, 0);
        return true;
    }
    return false;
}

bool AtomReader::import_file(const string &file_name, const bool add_noise) {
    string file_type = get_file_type(file_name);

    if (file_type == "xyz")
        import_xyz(file_name);
    else if (file_type == "ckx")
        import_ckx(file_name);
    else
        require(false, "Unimplemented file type: " + file_type);

    if (add_noise) {
        // initialize random seed
        srand (time(NULL));
        // add random number that is close to the distance between nn atoms to all the coordinates
        const double eps = 0.1 * get_point(0).distance(get_point(1));
        for (int i = 0; i < size(); ++i)
            atoms[i].point += eps * static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
    }

    calc_statistics();
    return calc_rms_distance();
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

// =================================
// *** WRITERS: ***************

void AtomReader::write_ckx(ofstream &out) const {
    // write the start of xyz header
    FileWriter::write_xyz(out);

    // write Ovito header
    out << "Medium properties=type:I:1:pos:R:3" << endl;

    // write data
    const int n_atoms = size();
    for (int i = 0; i < n_atoms; ++i)
        out << atoms[i].marker << " " << atoms[i].point << endl;
}

void AtomReader::write_restart(ofstream &out) const {
    const int n_atoms = size();
    if (velocities.size() != n_atoms) return;

    // start xyz header
    FileWriter::write_xyz(out);

    // write Ovito header
    out << "properties=species:S:1:pos:R:3:type:I:1:id:I:1:velocity:R:3" << endl;

    // TODO obtain labels from conf

    // write data
    for (int i = 0; i < n_atoms; ++i) {
        int marker = get_marker(i);
        if (marker == TYPES.EVAPORATED || marker == TYPES.CLUSTER)
            out << "Ev " << atoms[i].point << " 0 " << i+1 << " " << velocities[i] << "\n";
        else if (marker == TYPES.FIXED)
            out << "Cu " << atoms[i].point << " -1 " << i+1 << " " << velocities[i] << "\n";
        else
            out << "Cu " << atoms[i].point << " 1 " << i+1 << " " << velocities[i] << "\n";
    }
}

} /* namespace femocs */
