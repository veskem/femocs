/*
 * Mesh.cpp
 *
 *  Created on: 1.2.2016
 *      Author: veske
 */

#include "Mesh.h"

using namespace std;
namespace femocs {

Centre::Centre(double x, double y, double z) {
    this->x = x;
    this->y = y;
    this->z = z;
    //this->r2 = x*x + y*y + z*z;
}

bool Centre::is_equal(Centre centre) {
    const double eps = 0.1;
    bool dif1 = fabs(this->x - centre.x) < eps;
    bool dif2 = fabs(this->y - centre.y) < eps;
    bool dif3 = fabs(this->z - centre.z) < eps;
    return dif1 && dif2 && dif3;
}

// Mesh constructor
Mesh::Mesh() {
    inodes = 0;
    ielems = 0;
    ifaces = 0;
    i_nodemarker = 0;

    stat.Vmin = stat.Vmax = stat.Vmedian = stat.Vaverage = 0;
}

// Mesh destructor
Mesh::~Mesh() {
//    delete tetIO;
//    delete stat;
//    delete tetgenbeh;
//
//    delete nodemarkers;
//    delete facemarkers;
//    delete elemmarkers;
//    delete volumes;
//    delete centres;
}

// =================================
// *** GETTERS: ***************

const double Mesh::get_x(int i) {
    return tetIO.pointlist[3*i+0];
}

const double Mesh::get_y(int i) {
    return tetIO.pointlist[3*i+1];
}

const double Mesh::get_z(int i) {
    return tetIO.pointlist[3*i+2];
}

const double Mesh::get_node(int i, int xyz) {
#if DEBUGMODE
    if(  (i < 0) || (i >= get_n_nodes()) || (xyz < 0) || (xyz > 2) )
        cout << "Index out of number of nodes or coordinates bounds!" << endl;
#endif
    return tetIO.pointlist[3*i+xyz];
}

const int Mesh::get_face(int i, int node) {
#if DEBUGMODE
    if(  (i < 0) || (i >= get_n_faces()) || (node < 0) || (node > 2) ) {
        cout << "Index out of number of faces or nodes bounds!" << endl;
        return -1;
    }
#endif
    return tetIO.trifacelist[3*i+node];
}

const int Mesh::get_elem(int i, int node) {
#if DEBUGMODE
    if( (i < 0) || (i >= get_n_elems()) || (node < 0) || (node > 3) ) {
        cout << "Index out of number of elements or nodes bounds!" << endl;
        return -1;
    }
#endif
    return tetIO.tetrahedronlist[4*i+node];
}

const double Mesh::get_face_centre(int i, int xyz) {
#if DEBUGMODE
    if( (i < 0) || (i >= get_n_faces()) || (xyz < 0) || (xyz > 2) ) {
        cout << "Index out of number of faces or coordinates bounds!" << endl;
        return -1.0;
    }
#endif
    double node1 = get_node(get_face(i, 0), xyz);
    double node2 = get_node(get_face(i, 1), xyz);
    double node3 = get_node(get_face(i, 2), xyz);
    return (node1 + node2 + node3) / 3.0;
}

const double Mesh::get_elem_centre(int i, int xyz) {
#if DEBUGMODE
    if( (i < 0) || (i >= get_n_elems()) || (xyz < 0) || (xyz > 2) ) {
        cout << "Index out of number of elements or coordinates bounds!" << endl;
        return -1.0;
    }
#endif
    double node1 = get_node(get_elem(i, 0), xyz);
    double node2 = get_node(get_elem(i, 1), xyz);
    double node3 = get_node(get_elem(i, 2), xyz);
    double node4 = get_node(get_elem(i, 3), xyz);
    return (node1 + node2 + node3 + node4) / 4.0;
}

double* Mesh::get_nodes() {
    return tetIO.pointlist;
}

int* Mesh::get_faces() {
    return tetIO.trifacelist;
}

int* Mesh::get_elems() {
    return tetIO.tetrahedronlist;
}

double Mesh::get_volume(const int i) {
    return volumes[i];
}

Centre Mesh::get_centre(const int i) {
    return centres[i];
}

const int Mesh::get_nodemarker(const int i) {
#if DEBUGMODE
    if (i >= tetIO.numberofpoints) {
        cout << "Index exceeds node attributes size!" << endl;
        return -1;
    }
#endif
    return tetIO.pointmarkerlist[i];
    //    return nodemarkers[i];
}

const int Mesh::get_facemarker(const int i) {
    return facemarkers[i];
}

const int Mesh::get_elemmarker(const int i) {
    return elemmarkers[i];
}

int* Mesh::get_nodemarkers() {
    return tetIO.pointmarkerlist;
}

vector<int>* Mesh::get_facemarkers() {
    return &facemarkers;
}

vector<int>* Mesh::get_elemmarkers() {
    return &elemmarkers;
}

void Mesh::set_facemarker(const int i, const int m) {
    facemarkers[i] = m;
}

const int Mesh::get_n_nodes() {
    return tetIO.numberofpoints;
}

const int Mesh::get_n_elems() {
    return tetIO.numberoftetrahedra;
}

const int Mesh::get_n_faces() {
    return tetIO.numberoftrifaces;
}

const int Mesh::get_n_nodemarkers() {
    return tetIO.numberofpoints;
    //    return nodemarkers.size();
}

const int Mesh::get_n_facemarkers() {
    return facemarkers.size();
}

const int Mesh::get_n_elemmarkers() {
    return elemmarkers.size();
}

const int Mesh::get_n_volumes() {
    return volumes.size();
}

// =================================
// *** INITIALIZERS: ***************

void Mesh::init_nodemarkers(const int N) {
//    nodemarkers.reserve(N);
//    tetIO.numberofpointattributes = N;
    tetIO.pointmarkerlist = new int[N];
    i_nodemarker = 0;
}

void Mesh::init_facemarkers(const int N) {
    facemarkers.reserve(N);
}

void Mesh::init_elemmarkers(const int N) {
    elemmarkers.reserve(N);
}

void Mesh::init_nodes(const int N) {
    inodes = 0;
    tetIO.numberofpoints = N;
    tetIO.pointlist = new REAL[3 * N];
}

void Mesh::init_faces(const int N) {
    ifaces = 0;
    tetIO.numberoftrifaces = N;
    tetIO.trifacelist = new int[3 * N];
}

void Mesh::init_elems(const int N) {
    ielems = 0;
    tetIO.numberoftetrahedra = N;
    tetIO.tetrahedronlist = new int[4 * N];
}

void Mesh::init_volumes(const int N) {
    volumes.reserve(N);
}

void Mesh::init_centres(const int N) {
    centres.reserve(N);
}

// =================================
// *** ADDERS: ***************

void Mesh::add_centre(const double x, const double y, const double z) {
    centres.push_back(Centre(x,y,z));
}

void Mesh::add_volume(const double V) {
    volumes.push_back(V);
}

void Mesh::add_nodemarker(const int m) {
#if DEBUGMODE
    if (i_nodemarker >= tetIO.numberofpoints) {
        cout << "Node marker list is full!" << endl;
        return;
    }
#endif
    tetIO.pointmarkerlist[i_nodemarker] = m;
    i_nodemarker++;
//    nodemarkers.push_back(m);
}

void Mesh::add_facemarker(const int m) {
    facemarkers.push_back(m);
}

void Mesh::add_elemmarker(const int m) {
    elemmarkers.push_back(m);
}

void Mesh::add_elem(const int e1, const int e2, const int e3, const int e4) {
    int i = 4 * ielems;
    tetIO.tetrahedronlist[i + 0] = e1;
    tetIO.tetrahedronlist[i + 1] = e2;
    tetIO.tetrahedronlist[i + 2] = e3;
    tetIO.tetrahedronlist[i + 3] = e4;
    ielems++;
}

void Mesh::add_face(const int f1, const int f2, const int f3) {
    int i = 3 * ifaces;
    tetIO.trifacelist[i + 0] = f1;
    tetIO.trifacelist[i + 1] = f2;
    tetIO.trifacelist[i + 2] = f3;
    ifaces++;
}

void Mesh::add_node(const double x, const double y, const double z) {
    int i = 3 * inodes;
    tetIO.pointlist[i + 0] = (REAL) x;
    tetIO.pointlist[i + 1] = (REAL) y;
    tetIO.pointlist[i + 2] = (REAL) z;
    inodes++;
}

// =================================
// *** REPLICATORS: ***************

void Mesh::copy_statistics(Mesh* mesh) {
    stat.Vmin = mesh->stat.Vmin;
    stat.Vmax = mesh->stat.Vmax;
    stat.Vaverage = mesh->stat.Vaverage;
    stat.Vmedian = mesh->stat.Vmedian;
}

void Mesh::copy_nodes(Mesh* mesh) {
    int N = mesh->get_n_nodes();
    for (int i = 0; i < 3 * N; ++i)
        tetIO.pointlist[i] = mesh->get_nodes()[i];
    inodes = N;
}

void Mesh::copy_faces(Mesh* mesh, const int offset) {
    int N = mesh->get_n_faces();
    if (offset == 0)
        for (int i = 0; i < 3 * N; ++i)
            tetIO.trifacelist[i] = mesh->get_faces()[i];
    else
        for (int i = 0; i < 3 * N; ++i)
            tetIO.trifacelist[i] = offset + mesh->get_faces()[i];
    ifaces = N;
}

void Mesh::copy_elems(Mesh* mesh, const int offset) {
    int N = mesh->get_n_elems();
    if (offset == 0)
        for (int i = 0; i < 4 * N; ++i)
            tetIO.tetrahedronlist[i] = mesh->get_elems()[i];
    else
        for (int i = 0; i < 4 * N; ++i)
            tetIO.tetrahedronlist[i] = offset + mesh->get_elems()[i];
    ielems = N;
}

void Mesh::copy_nodemarkers(Mesh* mesh) {
    int N = mesh->get_n_nodemarkers();
    for (int i = 0; i < N; ++i)
        tetIO.pointmarkerlist[i] = mesh->get_nodemarker(i);
    i_nodemarker = N;
}

void Mesh::copy_facemarkers(Mesh* mesh) {
    int N = mesh->get_n_facemarkers();
    for (int i = 0; i < N; ++i)
        facemarkers.push_back(mesh->get_facemarker(i));
}

void Mesh::copy_elemmarkers(Mesh* mesh) {
    int N = mesh->get_n_elemmarkers();
    for (int i = 0; i < N; ++i)
        elemmarkers.push_back(mesh->get_elemmarker(i));
}

// =================================
// *** VARIA: ***************

void Mesh::calc_centres() {
    int i, j, k, n1, n2, n3, n4;
    double xyz[3];

    const int ncoords = 3; // nr of coordinates
    const int nnodes = 4;  // nr of nodes per element
    int N = get_n_elems();
    double* nodes = get_nodes();
    int* elems = get_elems();

    init_centres(N);

    // Loop through the elements
    for (i = 0; i < N; ++i) {
        j = nnodes * i;
        // Loop through x, y and z coordinates
        for (k = 0; k < ncoords; ++k) {
            n1 = ncoords * elems[j + 0] + k; // index of x, y or z coordinate of 1st node
            n2 = ncoords * elems[j + 1] + k; // ..2nd node
            n3 = ncoords * elems[j + 2] + k; // ..3rd node
            n4 = ncoords * elems[j + 3] + k; // ..4th node
            // Calculate the centre coordinate of element
            xyz[k] = (nodes[n1] + nodes[n2] + nodes[n3] + nodes[n4]) / nnodes;
        }
        add_centre(xyz[0], xyz[1], xyz[2]);
    }
}

void Mesh::calc_volumes() {
    int i, j, k, n1, n2, n3, n4;
    double V;
    double u[3], v[3], w[3];

    const int ncoords = 3; // nr of coordinates
    const int nnodes = 4;  // nr of nodes per element
    int N = get_n_elems();
    double* nodes = get_nodes();
    int* elems = get_elems();

    // Loop through the elements
    for (i = 0; i < N; ++i) {
        j = nnodes * i;
        // Loop through x, y and z coordinates
        for (k = 0; k < ncoords; ++k) {
            n1 = ncoords * elems[j + 0] + k; // index of x, y or z coordinate of 1st node
            n2 = ncoords * elems[j + 1] + k; // ..2nd node
            n3 = ncoords * elems[j + 2] + k; // ..3rd node
            n4 = ncoords * elems[j + 3] + k; // ..4th node

            u[k] = nodes[n1] - nodes[n2];
            v[k] = nodes[n1] - nodes[n3];
            w[k] = nodes[n1] - nodes[n4];
        }
        V = u[0] * (v[1] * w[2] - v[2] * w[1]) - u[1] * (v[0] * w[2] - v[2] * w[0])
                + u[2] * (v[0] * w[1] - v[1] * w[0]);
        add_volume(fabs(V) / 6.0);
    }
}

void Mesh::calc_volume_statistics() {
    size_t size = volumes.size();
    // Make a copy of the volumes vector
    vector<double> tempvec;
    tempvec.reserve(size);
    copy(volumes.begin(), volumes.end(), back_inserter(tempvec));

    sort(tempvec.begin(), tempvec.end());

    stat.Vmin = tempvec[0];
    stat.Vmax = tempvec[size - 1];
    stat.Vaverage = accumulate(tempvec.begin(), tempvec.end(), 0) / size;

    if (size % 2 == 0)
        stat.Vmedian = (tempvec[size / 2 - 1] + tempvec[size / 2]) / 2;
    else
        stat.Vmedian = tempvec[size / 2];
}

void Mesh::transform_elemmarkers() {
    int N = tetIO.numberoftetrahedronattributes;

    cout << "nelemmarkers=" << N << endl;

    init_elemmarkers(N);
    for (int i = 0; i < N; ++i)
        add_elemmarker((int) 10.0 * tetIO.tetrahedronattributelist[i]);
}

// Function to perform tetgen calculation on input and output data
void Mesh::recalc(const string cmd) {
    tetgenbeh.parse_commandline(const_cast<char*>(cmd.c_str()));
    tetrahedralize(&tetgenbeh, &tetIO, &tetIO);
}

void Mesh::output(const string cmd) {
    tetgenbeh.parse_commandline(const_cast<char*>(cmd.c_str()));
    tetrahedralize(&tetgenbeh, &tetIO, NULL);
}

// Function to output mesh in .vtk format
void Mesh::write_vtk(const string file_name, const int nnodes, const int ncells,
        const int nnodemarkers, const int nmarkers, const REAL* nodes, const int* cells,
        const int* nodemarkers, const vector<int>* markers, const int celltype,
        const int nnodes_in_cell) {
    int i, j;
    char file_name_char[1024];
    strcpy(file_name_char, file_name.c_str());

    FILE *outfile;
    outfile = fopen(file_name_char, "w");
    if (outfile == (FILE *) NULL) {
        printf("File I/O Error:  Cannot create file %s.\n", file_name_char);
        return;
    }

//    fprintf(outfile, "# vtk DataFile Version 2.0\n");
    fprintf(outfile, "# vtk DataFile Version 3.0\n");
    fprintf(outfile, "# Unstructured grid\n");
    fprintf(outfile, "ASCII\n"); // another option is BINARY
    fprintf(outfile, "DATASET UNSTRUCTURED_GRID\n\n");

    // Output the nodes
    if (nnodes > 0) {
        fprintf(outfile, "POINTS %d double\n", nnodes);
        for (i = 0; i < 3 * nnodes; i += 3)
            fprintf(outfile, "%.8g %.8g %.8g\n", nodes[i + 0], nodes[i + 1], nodes[i + 2]);
        fprintf(outfile, "\n");
    }

    // Output the cells (tetrahedra or triangles)
    if (ncells > 0) {
        fprintf(outfile, "CELLS %d %d\n", ncells, ncells * (nnodes_in_cell + 1));
        for (i = 0; i < nnodes_in_cell * ncells; i += nnodes_in_cell) {
            fprintf(outfile, "%d ", nnodes_in_cell);
            for (j = 0; j < nnodes_in_cell; ++j)
                fprintf(outfile, "%d ", cells[i + j]);
            fprintf(outfile, "\n");
        }
        fprintf(outfile, "\n");
    }

    // Output the types of cells, 10=tetrahedron, 5=triangle
    if (ncells > 0) {
        fprintf(outfile, "CELL_TYPES %d\n", ncells);
        for (i = 0; i < ncells; ++i)
            fprintf(outfile, "%d ", celltype);
//            fprintf(outfile, "%d\n", celltype);
        fprintf(outfile, "\n\n");
    }

    // Output point attributes
//    if (nnodemarkers > 0) {
//        fprintf(outfile, "POINT_DATA %d\n", nnodemarkers);
//        fprintf(outfile, "SCALARS Point_markers int\n");
//        fprintf(outfile, "LOOKUP_TABLE default\n");
//        for (i = 0; i < nnodemarkers; ++i)
//            fprintf(outfile, "%d\n", nodemarkers[i]);
//        fprintf(outfile, "\n");
//    }

    // Output cell attributes
    if (nmarkers > 0) {
        fprintf(outfile, "CELL_DATA %d\n", nmarkers);
        fprintf(outfile, "SCALARS Cell_markers int\n");
        fprintf(outfile, "LOOKUP_TABLE default\n");
        for (i = 0; i < nmarkers; ++i)
            fprintf(outfile, "%d\n", (*markers)[i]);
        fprintf(outfile, "\n");
    }

    fclose(outfile);
}

// Function to output faces in .vtk format
void Mesh::write_faces(const string file_name) {
#if not DEBUGMODE
    return;
#endif

    const int celltype = 5; // 5-triangle, 10-tetrahedron
    const int nnodes_in_cell = 3;

    int nnodes = get_n_nodes();
    int nfaces = get_n_faces();
    int nmarkers = get_n_facemarkers();
    int nnodemarkers = get_n_nodemarkers();
    REAL* nodes = get_nodes();          // pointer to nodes
    int* faces = get_faces();         // pointer to face nodes
    int* nodemakers = get_nodemarkers();
    vector<int>* cellmarkers = &facemarkers;    // pointer to face markers

    write_vtk(file_name, nnodes, nfaces, nnodemarkers, nmarkers, nodes, faces, nodemakers,
            cellmarkers, celltype, nnodes_in_cell);
}

// Function to output faces in .vtk format
void Mesh::write_elems(const string file_name) {
#if not DEBUGMODE
    return;
#endif
    const int celltype = 10; // 5-triangle, 10-tetrahedron
    const int nnodes_in_cell = 4;

    int nnodes = get_n_nodes();
    int nelems = get_n_elems();
    int nmarkers = get_n_elemmarkers();
    int nnodemarkers = get_n_nodemarkers();
    REAL* nodes = get_nodes();          // pointer to nodes
    int* elems = get_elems();     // pointer to faces nodes
    int* nodemakers = get_nodemarkers();
    vector<int>* cellmarkers = &elemmarkers;    // pointer to element markers

    write_vtk(file_name, nnodes, nelems, nnodemarkers, nmarkers, nodes, elems, nodemakers,
            cellmarkers, celltype, nnodes_in_cell);
}

// Function to output nodes in .xyz format
void Mesh::write_nodes(const string file_name) {
#if not DEBUGMODE
    return;
#endif
    int n_nodes = get_n_nodes();
    int n_markers = get_n_elemmarkers();
    REAL* nodes = get_nodes();          // pointer to nodes

    ofstream myfile;
    myfile.open(file_name);
    myfile << n_nodes << "\n";
    myfile << "Nodes of a mesh\n";

    for (int i = 0; i < n_nodes; ++i) {
        myfile << i << " ";
        myfile << get_x(i) << " ";
        myfile << get_y(i) << " ";
        myfile << get_z(i) << " ";
        myfile << get_nodemarker(i) << endl;
    }
    myfile.close();
}

// =================================
} /* namespace femocs */
