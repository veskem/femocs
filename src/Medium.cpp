/*
 * Medium.cpp
 *
 *  Created on: 21.1.2016
 *      Author: veske
 */

#include "Medium.h"
#include <float.h>
#include <fstream>
#include <numeric>

#if USE_CGAL
#include <CGAL/hilbert_sort.h>
#endif

using namespace std;
namespace femocs {

Medium::Medium() {
    init_statistics();
    reserve(0);
}

Medium::Medium(const int n_atoms) {
    init_statistics();
    reserve(n_atoms);
}

void Medium::calc_linked_list(const double r_cut, const bool lat_periodic) {
    const int n_atoms = size();

    Point3 simubox_size(sizes.xbox, sizes.ybox, sizes.zbox);
    for (int j = 0; j < 3; ++j) {
        nborbox_size[j] = ceil(1.0 * simubox_size[j] / r_cut);
        require(nborbox_size[j] > 0,
                "Invalid " + to_string(j) + "th nborbox size: " + to_string(nborbox_size[j]));
    }

    head = vector<int>(nborbox_size[0]*nborbox_size[1]*nborbox_size[2], -1);
    list = vector<int>(n_atoms, -1);
    nborbox_indices.clear();
    nborbox_indices.reserve(n_atoms);

    // calculate linked list for the atoms
    for (int i = 0; i < n_atoms; ++i) {
        Point3 dx = atoms[i].point + simubox_size / 2;

        // Check that we are inside lateral boundaries
        if (lat_periodic) {
            if (dx.x < 0) dx.x += simubox_size.x;
            if (dx.x > simubox_size.x) dx.x -= simubox_size.x;
            if (dx.y < 0) dx.y += simubox_size.y;
            if (dx.y > simubox_size.y) dx.y -= simubox_size.y;
            if (dx.z < 0) dx.z += simubox_size.z;
            if (dx.z > simubox_size.z) dx.z -= simubox_size.z;
        }

        array<int,3> point_index;
        for (int j = 0; j < 3; ++j)
            point_index[j] = int( (dx[j] / simubox_size[j]) * nborbox_size[j] );

        // If not periodic, let border cells continue to infinity
        if (!lat_periodic)
            for (int j = 0; j < 3; ++j) {
                point_index[j] = max(0, point_index[j]);
                point_index[j] = min(nborbox_size[j]-1, point_index[j]);
            }

        int i_cell = (point_index[2] * nborbox_size[1] + point_index[1]) * nborbox_size[0] + point_index[0];
        for (int j = 0; j < 3; ++j) {
            require(point_index[j] >= 0 && point_index[j] < nborbox_size[j],
                    "Invalid " + to_string(j) + "th point nbor index: " + to_string(point_index[j]));
        }
        require(i_cell >= 0 && i_cell < head.size(), "Invalid neighbouring cell index: " + to_string(i_cell));

        nborbox_indices.push_back(point_index);
        list[i] = head[i_cell];
        head[i_cell] = i;
        set_marker(i, i_cell);
    }
}

void Medium::calc_linked_list2(const double r_cut) {
    const int n_atoms = size();

    array<double,3> rc = {sizes.xbox, sizes.ybox, sizes.zbox};
    array<int,3> lc {ceil(rc[0]/r_cut), ceil(rc[1]/r_cut), ceil(rc[2]/r_cut)};
    array<int,3> mc;

    int lcyz = lc[1]*lc[2];
    int lcxyz = lcyz * lc[0];

    vector<int> lscl(n_atoms);

    /* Reset the headers, head */
    vector<int> head(lcxyz, -1);


    /* Scan atoms to construct headers, head, & linked lists, lscl */
    for (int i = 0; i < n_atoms; i++) {
        /* Vector cell index to which this atom belongs */
        Point3 r = get_point(i);
        for (int a=0; a<3; a++)
            mc[a] = r[a] / rc[a];
        /* Translate the vector cell index, mc, to a scalar cell index */
        int c = mc[0]*lcyz+mc[1]*lc[2]+mc[2];
        /* Link to the previous occupant (or EMPTY if you're the 1st) */
        lscl[i] = head[c];
        /* The last one goes to the header */
        head[c] = i;
    }
}

// algorithm from
// http://cacs.usc.edu/education/cs596/01-1LinkedListCell.pdf
// useful Google:
// linked list to verlet neighborlist c++
void Medium::calc_nborlist2(double r_cut) {
    const int n_atoms = size;

    vector<int> head;
    vector<int> lscl;
    vector<vector<int>> nborlist(n_atoms);

    array<double,3> rc = {sizes.xbox, sizes.ybox, sizes.zbox};
    array<int,3> lc {ceil(rc[0]/r_cut), ceil(rc[1]/r_cut), ceil(rc[2]/r_cut)};
    array<int,3> mc, mc1;
    Point3 rshift;

    int lcyz = lc[1]*lc[2];
    int lcxyz = lcyz * lc[0];

    /* Scan inner cells */
     for (mc[0]=0; mc[0]<lc[0]; (mc[0])++)
     for (mc[1]=0; mc[1]<lc[1]; (mc[1])++)
     for (mc[2]=0; mc[2]<lc[2]; (mc[2])++) {
         /* Calculate a scalar cell index */
         int c = mc[0] * lcyz+mc[1] * lc[2] + mc[2];
         /* Scan the neighbor cells (including itself) of cell c */
         for (mc1[0]=mc[0]-1; mc1[0]<=mc[0]+1; (mc1[0])++)
         for (mc1[1]=mc[1]-1; mc1[1]<=mc[1]+1; (mc1[1])++)
         for (mc1[2]=mc[2]-1; mc1[2]<=mc[2]+1; (mc1[2])++) {

             /* Periodic boundary condition by shifting coordinates */
             for (int a = 0; a < 3; a++) {
                 if (mc1[a] < 0)
                     rshift[a] = -rc[a];
                 else if (mc1[a]>=lc[a])
                     rshift[a] = rc[a];
                 else
                     rshift[a] = 0.0;
             }
             /* Calculate the scalar cell index of the neighbor cell */
             int c1 = ((mc1[0]+lc[0])%lc[0])*lcyz
                     +((mc1[1]+lc[1])%lc[1])*lc[2]
                     +((mc1[2]+lc[2])%lc[2]);
             /* Scan atom i in cell c */
             int i = head[c];
             while (i != -1) {
                 Point3 ri = get_point(i);
                 /* Scan atom j in cell c1 */
                 int j = head[c1];
                 while (j != -1) {
                     if (i < j) { /* Avoid double counting of pair (i, j) */
                         Point3 rj = get_point(j);
                         /* Image-corrected relative pair position */
                         if (ri.distance2(rj+rshift) < r_cut) {
                             // Compute forces on pair (i, j)...
                             nborlist[i].push_back(j);
                             nborlist[j].push_back(i);
                         }

                     }
                     j = lscl[j];
                 }
                 i = lscl[i];
             }
         }
     }
}

void Medium::calc_nborlist(vector<vector<int>>& nborlist, const double r_cut, const bool lat_periodic) {
    const int n_atoms = size();
    const double r_cut2 = r_cut * r_cut;

    Point3 simubox_size(sizes.xbox, sizes.ybox, sizes.zbox);
    for (int j = 0; j < 3; ++j) {
        nborbox_size[j] = ceil(1.0 * simubox_size[j] / r_cut);
        require(nborbox_size[j] > 0,
                "Invalid " + to_string(j) + "th nborbox size: " + to_string(nborbox_size[j]));
    }

    head = vector<int>(nborbox_size[0]*nborbox_size[1]*nborbox_size[2], -1);
    list = vector<int>(n_atoms, -1);
    nborbox_indices.clear();
    nborbox_indices.reserve(n_atoms);

    // calculate linked list for the atoms
    for (int i = 0; i < n_atoms; ++i) {
        Point3 dx = atoms[i].point + simubox_size / 2;

        // Check that we are inside lateral boundaries
        if (lat_periodic) {
            if (dx.x < 0) dx.x += simubox_size.x;
            if (dx.x > simubox_size.x) dx.x -= simubox_size.x;
            if (dx.y < 0) dx.y += simubox_size.y;
            if (dx.y > simubox_size.y) dx.y -= simubox_size.y;
            if (dx.z < 0) dx.z += simubox_size.z;
            if (dx.z > simubox_size.z) dx.z -= simubox_size.z;
        }

        array<int,3> point_index;
        for (int j = 0; j < 3; ++j)
            point_index[j] = int( (dx[j] / simubox_size[j]) * nborbox_size[j] );

        // If not periodic, let border cells continue to infinity
        if (!lat_periodic)
            for (int j = 0; j < 3; ++j) {
                point_index[j] = max(0, point_index[j]);
                point_index[j] = min(nborbox_size[j]-1, point_index[j]);
            }

        int i_cell = (point_index[2] * nborbox_size[1] + point_index[1]) * nborbox_size[0] + point_index[0];
        for (int j = 0; j < 3; ++j) {
            require(point_index[j] >= 0 && point_index[j] < nborbox_size[j],
                    "Invalid " + to_string(j) + "th point nbor index: " + to_string(point_index[j]));
        }
        require(i_cell >= 0 && i_cell < head.size(), "Invalid neighbouring cell index: " + to_string(i_cell));

        nborbox_indices.push_back(point_index);
        list[i] = head[i_cell];
        head[i_cell] = i;
        set_marker(i, i_cell);
    }



    nborlist = vector<vector<int>>(n_atoms);

    for (int i = 0; i < n_atoms; ++i) {
        array<int,3> i_atom = nborbox_indices[i];
        Point3 point1 = atoms[i].point;

        // loop through the boxes where the neighbours are located; there are up to 3^3=27 boxes
        for (int iz = i_atom[2]-1; iz <= i_atom[2]+1; ++iz) {
            // some of the iterations can be skipped if the box is on the simu box boundary
            if (iz < 0 || iz >= nborbox_size[2]) continue;
            for (int iy = i_atom[1]-1; iy <= i_atom[1]+1; ++iy) {
                if (iy < 0 || iy >= nborbox_size[1]) continue;
                for (int ix = i_atom[0]-1; ix <= i_atom[0]+1; ++ix) {
                    if (ix < 0 || ix >= nborbox_size[0]) continue;

                    // transform volumetric neighbour box index to linear one
                    int i_cell = (iz * nborbox_size[1] + iy) * nborbox_size[0] + ix;
                    require(i_cell >= 0 && i_cell < head.size(), "Invalid neighbouring cell index: " + to_string(i_cell));

                    // get the index of first atom in given neighbouring cell
                    int j = head[i_cell];

                    // loop through all atoms in a given neighbouring cell
                    while(j >= 0) {
                        if (i != j && point1.distance2(get_point(j)) <= r_cut2) {
                            nborlist[i].push_back(j);
                        }
                        j = list[j];
                    }
                }
            }
        }
    }
}





// Sort the atoms by their cartesian or radial coordinate
void Medium::sort_atoms(const int coord, const string& direction) {
    require(coord >= 0 && coord <= 3, "Invalid coordinate: " + to_string(coord));

    if (size() < 2) return;

    if (coord == 3) {
        Point2 origin(sizes.xmid, sizes.ymid);
        for (int i = 0; i < size(); ++i)
            set_marker(i, 10000 * origin.distance2(get_point2(i)));
        if (direction == "up" || direction == "asc")
            sort(atoms.begin(), atoms.end(), Atom::sort_marker_up());
        else
            sort(atoms.begin(), atoms.end(), Atom::sort_marker_down());
    }

    else if (direction == "up" || direction == "asc")
        sort(atoms.begin(), atoms.end(), Atom::sort_up(coord));
    else if (direction == "down" || direction == "desc")
        sort(atoms.begin(), atoms.end(), Atom::sort_down(coord));
}

// Sort the atoms first by 1st and then by 2nd cartesian coordinate
void Medium::sort_atoms(const int x1, const int x2, const string& direction) {
    require(x1 >= 0 && x1 <= 2 && x2 >= 0 && x2 <= 2, "Invalid coordinates: " + to_string(x1) + ", " + to_string(x2));

    if (direction == "up" || direction == "asc")
        sort( atoms.begin(), atoms.end(), Atom::sort_up2(x1, x2) );
    else if (direction == "down" || direction == "desc")
        sort( atoms.begin(), atoms.end(), Atom::sort_down2(x1, x2) );
}

// Perform spatial sorting by ordering atoms along Hilbert curve
void Medium::sort_spatial() {
#if USE_CGAL
    CGAL::hilbert_sort( atoms.begin(), atoms.end(), Atom::sort_spatial(), CGAL::Hilbert_sort_middle_policy() );
#endif
}

// Reserve memory for data vectors
void Medium::reserve(const int n_atoms) {
    require(n_atoms >= 0, "Invalid number of atoms: " + to_string(n_atoms));
    atoms.clear();
    atoms.reserve(n_atoms);
}

// Reserve memory for data vectors
void Medium::resize(const int n_atoms) {
    require(n_atoms >= 0, "Invalid number of atoms: " + to_string(n_atoms));
    atoms.reserve(n_atoms);
}

// Define the addition of two Mediums
Medium& Medium::operator +=(const Medium &m) {
    atoms.insert(atoms.end(), m.atoms.begin(), m.atoms.end());
    calc_statistics();
    return *this;
}

// Add atom to atoms vector
void Medium::append(const Atom& atom) {
    expect((unsigned)size() < atoms.capacity(), "Allocated vector size exceeded!");
    atoms.push_back(atom);
}

// Add atom with defalt id and marker to atoms vector
void Medium::append(const Point3& point) {
    expect((unsigned)size() < atoms.capacity(), "Allocated vector sizes exceeded!");
    atoms.push_back(Atom(-1, point, 0));
}

// Initialize statistics about Medium
void Medium::init_statistics() {
    sizes.xmin = sizes.ymin = sizes.zmin = DBL_MAX;
    sizes.xmax = sizes.ymax = sizes.zmax =-DBL_MAX;
    sizes.xmean = sizes.ymean = sizes.zmean = 0.0;
    sizes.xmid = sizes.ymid = sizes.zmid = 0.0;

    sizes.xbox = sizes.ybox = sizes.zbox = 0;
    sizes.zminbox = DBL_MAX;
    sizes.zmaxbox =-DBL_MAX;
}

// Calculate the statistics about Medium
void Medium::calc_statistics() {
    int n_atoms = size();
    init_statistics();
    if (n_atoms <= 0) {
        expect(false, "Can't calculate statistics for empty set of atoms!");
        return;
    }

    Point3 average(0,0,0);

    // Find min and max coordinates
    for (int i = 0; i < n_atoms; ++i) {
        Point3 point = get_point(i);
        average += point;
        sizes.xmin = min(sizes.xmin, point.x);
        sizes.xmax = max(sizes.xmax, point.x);
        sizes.ymin = min(sizes.ymin, point.y);
        sizes.ymax = max(sizes.ymax, point.y);
        sizes.zmin = min(sizes.zmin, point.z);
        sizes.zmax = max(sizes.zmax, point.z);
    }

    // Define average coordinates
    average *= (1.0 / n_atoms);
    sizes.xmean = average.x;
    sizes.ymean = average.y;
    sizes.zmean = average.z;

    // Define size of simubox
    sizes.xbox = sizes.xmax - sizes.xmin;
    sizes.ybox = sizes.ymax - sizes.ymin;
    sizes.zbox = sizes.zmax - sizes.zmin;
    sizes.zminbox = sizes.zmin;
    sizes.zmaxbox = sizes.zmax;

    // Define the centre of simubox
    sizes.xmid = (sizes.xmax + sizes.xmin) / 2;
    sizes.ymid = (sizes.ymax + sizes.ymin) / 2;
    sizes.zmid = (sizes.zmax + sizes.zmin) / 2;
}

// Get number of atoms in Medium
int Medium::size() const {
    return atoms.size();
}

// Get 2-dimensional coordinates of i-th atom
Point2 Medium::get_point2(const int i) const {
    require(i >= 0 && i < size(), "Index out of bounds: " + to_string(i));
    return Point2(atoms[i].point.x, atoms[i].point.y);
}

// Get 3-dimensional coordinates of i-th atom
Point3 Medium::get_point(const int i) const {
    require(i >= 0 && i < size(), "Index out of bounds: " + to_string(i) + "/" + to_string(size()));
    return atoms[i].point;
}

// Get atom ID
int Medium::get_id(const int i) const {
    require(i >= 0 && i < size(), "Index out of bounds: " + to_string(i));
    return atoms[i].id;
}

// Get atom marker
int Medium::get_marker(const int i) const {
    require(i >= 0 && i < size(), "Index out of bounds: " + to_string(i));
    return atoms[i].marker;
}

// Get i-th Atom
Atom Medium::get_atom(const int i) const {
    require(i >= 0 && i < size(), "Index out of bounds: " + to_string(i));
    return atoms[i];
}

// Set entry to id-s vector
void Medium::set_id(const int i, const int id) {
    require(i >= 0 && i < size(), "Index out of bounds: " + to_string(i));
    atoms[i].id = id;
}

// Set entry to point-s vector
void Medium::set_point(const int i, const Point3& p) {
    require(i >= 0 && i < size(), "Index out of bounds: " + to_string(i));
    atoms[i].point = p;
}

// Set entry to x coordinate vector
void Medium::set_x(const int i, const double x) {
    require(i >= 0 && i < size(), "Index out of bounds: " + to_string(i));
    atoms[i].point.x = x;
}

// Set entry to y coordinate vector
void Medium::set_y(const int i, const double y) {
    require(i >= 0 && i < size(), "Index out of bounds: " + to_string(i));
    atoms[i].point.y = y;
}

// Set entry to z coordinate vector
void Medium::set_z(const int i, const double z) {
    require(i >= 0 && i < size(), "Index out of bounds: " + to_string(i));
    atoms[i].point.z = z;
}

// Set atom marker
void Medium::set_marker(const int i, const int m) {
    require(i >= 0 && i < size(), "Index out of bounds: " + to_string(i));
    atoms[i].marker = m;
}

// Pick the suitable write function based on the file type
// Function is active only when file write is enabled
void Medium::write(const string &file_name) const {
    if (!MODES.WRITEFILE) return;
    
    int n_atoms = size();
    expect(n_atoms > 0, "Zero atoms detected!");
    string ftype = get_file_type(file_name);
    
    ofstream outfile;
//    outfile << fixed;

    outfile.setf(std::ios::scientific);
    outfile.precision(6);

    if (ftype == "movie") outfile.open(file_name, ios_base::app);
    else outfile.open(file_name);
    require(outfile.is_open(), "Can't open a file " + file_name);
    
    if (ftype == "xyz" || ftype == "movie")
        write_xyz(outfile, n_atoms);
    else if (ftype == "vtk")
        write_vtk(outfile, n_atoms);
    else if (ftype == "ckx")
        write_ckx(outfile, n_atoms);
    else    
        require(false, "Unsupported file type: " + ftype);

    outfile.close();
}

// Compile data string from the data vectors
string Medium::get_data_string(const int i) const {
    if(i < 0) return "Medium properties=id:I:1:pos:R:3:marker:I:1";

    ostringstream strs; strs << fixed;
    strs << atoms[i];
    return strs.str();
}

// Output atom data in .xyz format
void Medium::write_xyz(ofstream& out, const int n_atoms) const {
    out << n_atoms << "\n";
    out << get_data_string(-1) << endl;

    for (int i = 0; i < n_atoms; ++i)
        out << get_data_string(i) << endl;
}

// Output atom data in .vtk format
void Medium::write_vtk(ofstream& out, const int n_atoms) const {
    out << "# vtk DataFile Version 3.0\n";
    out << "# Medium data\n";
    out << "ASCII\n";
    out << "DATASET UNSTRUCTURED_GRID\n\n";

    // Output the point coordinates
    out << "POINTS " << n_atoms << " double\n";
    for (int i = 0; i < n_atoms; ++i)
        out << get_point(i) << "\n";

    get_cell_data(out);
}

// Output atom data in .ckx format
void Medium::write_ckx(ofstream &out, const int n_atoms) const {
    out << n_atoms << "\n";
    out << "Medium properties=type:I:1:pos:R:3" << endl;
    for (int i = 0; i < n_atoms; ++i)
        out << atoms[i].marker << " " << atoms[i].point << endl;
}

// Get scalar and vector data associated with vtk cells
void Medium::get_cell_data(ofstream& out) const {
    const int celltype = 1;     // type of the cell in vtk format; 1-vertex, 10-tetrahedron
    const int dim = 1;          // number of vertices in the cell
    const int n_cells = size();
    const int n_atoms = size();

    // Output the vertices
    out << "\nCELLS " << n_cells << " " << (1+dim) * n_cells << "\n";
    for (int i = 0; i < n_cells; ++i)
        out << dim << " " << i << "\n";

    // Output cell types
    out << "\nCELL_TYPES " << n_cells << "\n";
    for (int i = 0; i < n_cells; ++i)
        out << celltype << "\n";

    out << "\nPOINT_DATA " << n_atoms << "\n";

    // write IDs of atoms
    out << "SCALARS ID int\nLOOKUP_TABLE default\n";
    for (int i = 0; i < n_atoms; ++i)
        out << atoms[i].id << "\n";

    // write atom markers
    out << "SCALARS marker int\nLOOKUP_TABLE default\n";
    for (int i = 0; i < n_atoms; ++i)
        out << atoms[i].marker << "\n";
}

// Copy statistics from another Medium
void Medium::copy_statistics(const Medium& m) {
    require(this->sizes.size() == m.sizes.size() , "Incompatible statistics!");
    for (int i = 0; i < m.sizes.size(); ++i)
        (&this->sizes.xmin)[i] = (&m.sizes.xmin)[i];
}

} /* namespace femocs */
