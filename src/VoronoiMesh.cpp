/*
 * VoronoiMesh.cpp
 *
 *  Created on: 6.3.2017
 *      Author: veske
 */

#include "VoronoiMesh.h"
#include <float.h>

using namespace std;
namespace femocs {

/* =====================================================================
 *  ============================= Voronoi =============================
 * ===================================================================== */

// Get the area of the face
double VoronoiFace::area() {
    const int n_nodes = size();
    calc_verts();

    Vec3 total(0.0);
    for (int i = 0; i < n_nodes; ++i) {
        int j = i+1;
        if (j == n_nodes) j = 0;
        total += verts[i].crossProduct(verts[j]);
    }

    return total.norm() * 0.5;
}

// Get the centroid coordinates of the face
Point3 VoronoiFace::centroid() {
    calc_verts();
    Vec3 centroid(0.0);

    for (Vec3& v : verts)
        centroid += v;
    centroid *= (1.0 / size());

    return Point3(centroid.x, centroid.y, centroid.z);
}

// Return the neighbouring cell for the caller cell
int VoronoiFace::nborcell(const int caller_id) {
    const int c1 = data->vfacetlist[id].c1;
    if (c1 == caller_id)
        return data->vfacetlist[id].c2;
    return c1;
}

// Calculate the unique node that is associated with the edge
int VoronoiFace::get_node(const int edge) {
    static int node = -1;

    const int v1 = data->vedgelist[edge].v1;
    if (node != v1) node = v1;
    else node = data->vedgelist[edge].v2;

    return node;
}

// Transform the node data from tetgenio into easily accessible form
void VoronoiFace::calc_verts() {
    if (verts.size() > 0) return;
    verts.reserve(size());

    for (int edge : *this) {
        int n = 3 * get_node(edge);
        verts.push_back( Vec3(data->vpointlist[n+0], data->vpointlist[n+1], data->vpointlist[n+2]) );
    }
}

// Get the norm vector of the face
Vec3 VoronoiFace::norm() {
    Vec3 edge1 = verts[1] - verts[0];
    Vec3 edge2 = verts[2] - verts[0];
    return edge2.crossProduct(edge1).normalize();
}

vector<int> VoronoiCell::get_neighbours() const {
    vector<int> nbors; nbors.reserve(size());
    for (VoronoiFace face : *this)
        nbors.push_back(face.nborcell(id));
    return nbors;
}

/* =====================================================================
 *  ============================ Voronois ============================
 * ===================================================================== */

void VoronoiCells::write_cells(ofstream& out) const {
    // Get total number of faces and their vertices
    size_t n_faces = 0, n_verts = 0;
    for (VoronoiCell cell : *this)
        for (VoronoiFace face : cell) {
            n_verts += face.size();
            n_faces++;
        }

    // Write the Voronoi cells as polygons
    out << "\nCELLS " << n_faces << " " << n_faces + n_verts << "\n";
    for (VoronoiCell cell : *this)
        for (VoronoiFace face : cell)
            out << face.size() << " " << face << endl;

    // Output cell types
    out << "\nCELL_TYPES " << n_faces << "\n";
    for (size_t i = 0; i < n_faces; ++i)
        out << celltype << "\n";

    // Output scalar data associated with Voronoi cells
    out << "\nCELL_DATA " << n_faces << "\n";

    // Write the markers of cells
    out << "SCALARS marker int\nLOOKUP_TABLE default\n";
    for (VoronoiCell cell : *this)
        for (VoronoiFace face : cell)
            out << get_marker(cell.id) << "\n";

    // Write the IDs of cells
    out << "SCALARS ID int\nLOOKUP_TABLE default\n";
    for (VoronoiCell cell : *this)
        for (VoronoiFace face : cell)
            out << cell.id << "\n";
}

// Calculate the face neighbours - faces, that share the same edge
void VoronoiFaces::calc_neighbours() {
    const int n_faces = size();
    const int n_edges = tetio->numberofvedges;

    // Produce edge to face mapping
    vector<vector<int>> edge2face(n_edges);
    for (int i = 0; i < n_edges; ++i)
        edge2face[i].reserve(10);

    for (VoronoiFace face : *this)
        for (int i : face)
            edge2face[i].push_back(face.id);

    // Initialize face neighbour list
    neighbours.clear();
    neighbours.resize(n_faces);
    for (int i = 0; i < n_faces; ++i)
        neighbours[i].reserve(10);

    // Convert edge to face list to face neighbour list
    for (vector<int> faces : edge2face)
        for (int i : faces)
            for (int j : faces) {
                if (i == j) continue;
                neighbours[i].push_back(j);
            }
}

void VoronoiFaces::calc_centroids() {
    const int n_faces = size();
    centroids.clear();
    centroids.reserve(n_faces);

    for (VoronoiFace face : *this) {
        if (face.size() > 0)
            centroids.push_back(face.centroid());
        else
            centroids.push_back(Vec3());
    }
}

void VoronoiFaces::write_cells(ofstream& out) const {
    // Get total number of faces and their vertices
    size_t n_faces = 0, n_verts = 0;
    for (VoronoiFace face : *this) {
        n_verts += face.size();
        n_faces += face.size() > 0;
    }

    // Write the Voronoi faces as polygons
    out << "\nCELLS " << n_faces << " " << n_faces + n_verts << "\n";
    for (VoronoiFace face : *this)
        if (face.size() > 0)
            out << face.size() << " " << face << endl;

    // Output cell types
    out << "\nCELL_TYPES " << n_faces << "\n";
    for (size_t i = 0; i < n_faces; ++i)
        out << celltype << "\n";

    // Output scalar data associated with Voronoi faces
    out << "\nCELL_DATA " << n_faces << "\n";

    // write markers of faces
    out << "SCALARS marker int\nLOOKUP_TABLE default\n";
    for (VoronoiFace face : *this)
        if (face.size() > 0)
            out << get_marker(face.id) << "\n";

    // write IDs of faces
    out << "SCALARS ID int\nLOOKUP_TABLE default\n";
    for (VoronoiFace face : *this)
        if (face.size() > 0)
            out << face.id << "\n";

    // write face areas
    out << "SCALARS area double\nLOOKUP_TABLE default\n";
    for (VoronoiFace face : *this)
        if (face.size() > 0)
            out << face.area() << "\n";
}

/* ==================================================================
 *  ======================= VoronoiMesh ===========================
 * ================================================================== */

// Initialize Tetgen data
VoronoiMesh::VoronoiMesh() {
    tetIOin.initialize();
    tetIOout.initialize();
}

// Find the Voronoi cell that for sure belongs to the surface
int VoronoiMesh::get_seedcell() {
    double zmax = -1e100;
    int seed = -1;

    for (int i = nodes.indxs.surf_start; i <= nodes.indxs.surf_end; ++i) {
        double z = nodes[i].z;
        if (z > zmax) {
            zmax = z;
            seed = i;
        }
    }

    return seed;
}

// Mark the cell and faces that are certainly on the surface
int VoronoiMesh::mark_seed() {
    const int cell_max = nodes.indxs.surf_end;
    int seedface = -1, seedcell = get_seedcell();

    // mark the cell that surrounds the up-most atom
    voros.set_marker(seedcell, TYPES.ZMAX);

    Vec3 znorm(0, 0, 1);
    VoronoiCell cell = voros[seedcell];
    Vec3 centre = nodes.get_vec(seedcell);

    // mark the faces that are on the upper half of the cell

    vector<int> cell_nbors = cell.get_neighbours();
    for (size_t i = 0; i < cell_nbors.size(); ++i)
        if (cell_nbors[i] > cell_max) {
            VoronoiFace face = cell[i];
            // face is on the upper half of the cell if its norm is upwards
            Vec3 normal = nodes.get_vec(face.nborcell(cell.id)) - centre;
            if (znorm.dotProduct(normal) >= 0) {
                seedface = face.id;
                vfaces.set_marker(face.id, TYPES.SURFACE);
            }
        }
    require(seedface >= 0, "Finding seed Voronoi facet failed!");
    return seedface;
}

// Mark Voronoi faces that are on the vacuum-material boundary
void VoronoiMesh::mark_faces(const double zmin, const int seed) {
    const int cell_max = nodes.indxs.surf_end;
    vfaces.calc_neighbours();

    // Mark the faces that are associated only with potential surface cells
    // and therefore cannot be on the surface
    for (VoronoiCell cell : voros)
        for (VoronoiFace face : cell) {
            if ( (cell.id <= cell_max && face.nborcell(cell.id) <= cell_max)
                    || (cell.id > cell_max && face.nborcell(cell.id) > cell_max) )
                vfaces.set_marker(face.id, TYPES.PERIMETER);
        }

    // Mark the faces that have the node outside the boundaries of surface atoms
    for (VoronoiFace face : vfaces) {
        face.calc_verts();
        for (Vec3 vert : face.verts) {
            bool b1 = false; //vert.x < sizes.xmin || vert.x > sizes.xmax;
            bool b2 = false; //vert.y < sizes.ymin || vert.y > sizes.ymax;
            bool b3 = vert.z < zmin;
            if (b1 || b2 || b3) {
                vfaces.set_marker(face.id, TYPES.ZMIN);
                break;
            }
        }
    }

    vfaces.set_marker(seed, TYPES.ZMAX);

    // Mark the faces around the seed face
    vector<int> neighbours = vfaces.get_neighbours(seed);
    for (size_t i = 0; i < neighbours.size(); ++i) {
        int face = neighbours[i];
        if (vfaces.get_marker(face) == TYPES.NONE) {

            // Mark the face as surface
            vfaces.set_marker(face, TYPES.SURFACE);

            // Expand the list of possible surface faces
            vector<int> nbors = vfaces.get_neighbours(face);
            neighbours.insert(neighbours.end(), nbors.begin(), nbors.end());
        }
    }
}

void VoronoiMesh::calc_ranks(vector<int>& ranks, const int seedface) {
    const int n_nbor_layers = 4;  // number of nearest faces that act as a seed
    const int n_faces = vfaces.size();
    const double max_rank = 100.0;

    // initialise all the ranks to 0
    ranks = vector<int>(n_faces);

    // calculate the ranks from vacuum side
    vector<int> neighbours = vfaces.get_neighbours(seedface);
    for (size_t i = 0; i < neighbours.size(); ++i) {
        int vface = neighbours[i];
        if (vfaces.get_marker(vface) == TYPES.NONE && ranks[vface]++ == 0) {
            vector<int> nbors = vfaces.get_neighbours(vface);
            neighbours.insert(neighbours.end(), nbors.begin(),nbors.end());
        }
    }

    // normalise all the ranks with respect to the maximum rank
    double norm_factor = max_rank / *max_element(ranks.begin(), ranks.end());
    for (int& r : ranks)
        if (r > 0)
            r *= norm_factor;

    // force the ranks around the seed region to the maximum value
    vector<vector<int>> nbors(n_nbor_layers);
    ranks[seedface] = max_rank;
    for (int layer = 0; layer < n_nbor_layers; ++layer) {
        // build next layer of neighbour list
        if (layer == 0)
            nbors[0] = vfaces.get_neighbours(seedface);
        else {
            for (int nbor : nbors[layer-1])
                if (ranks[nbor] > 0) {
                    vector<int> tmp_nbors = vfaces.get_neighbours(nbor);
                    nbors[layer].insert(nbors[layer].end(), tmp_nbors.begin(), tmp_nbors.end());
                }
        }
        for (int nbor : nbors[layer])
            if (ranks[nbor] > 0)
                ranks[nbor] = max_rank;
    }
}

// Mark Voronoi faces that are on the vacuum-material boundary
void VoronoiMesh::mark_faces_vol2(const double zmin, const int seed) {
    const int cell_max = nodes.indxs.surf_end;
    const int min_rank = 40;

    vfaces.calc_neighbours();
    vfaces.set_marker(seed, TYPES.ZMAX);

    // Mark the faces that are associated only with potential surface cells
    // and therefore cannot be on the surface
    for (VoronoiCell cell : voros)
        for (VoronoiFace face : cell) {
            const bool b1 = cell.id <= cell_max;
            const bool b2 = face.nborcell(cell.id) <= cell_max;
            if ((b1 && b2) || (!b1 && !b2))
                vfaces.set_marker(face.id, TYPES.PERIMETER);
        }

    // Mark the faces that are in the bottom region
    for (int cell = 0; cell < voros.size(); ++cell)
        if (cell <= cell_max && nodes[cell].z < zmin) {
            for (VoronoiFace face : voros[cell])
                vfaces.set_marker(face.id, TYPES.ZMIN);
        }

//    vector<int> ranks;
//    calc_ranks(ranks, seed);

//    for (int face = 0; face < vfaces.size(); ++face)
//        if (vfaces.get_marker(face) == TYPES.NONE)
//            vfaces.set_marker(face, ranks[face]);
//
//    return;

    // Mark the faces around the seed face
    vector<int> neighbours = vfaces.get_neighbours(seed);
    for (size_t i = 0; i < neighbours.size(); ++i) {
        int face = neighbours[i];
        if (vfaces.get_marker(face) == TYPES.NONE) {
            vfaces.set_marker(face, TYPES.SURFACE);

            // Expand the list of possible surface faces
//            if (ranks[face] >= min_rank) {
//                vector<int> nbors = vfaces.get_neighbours(face);
//                neighbours.insert(neighbours.end(), nbors.begin(), nbors.end());
//            }
            vector<int> nbors = vfaces.get_neighbours(face);
            neighbours.insert(neighbours.end(), nbors.begin(), nbors.end());
        }
    }
}

// Mark Voronoi cells and nodes that are on the surface of material
void VoronoiMesh::mark_cells_and_nodes() {
    const int cell_max = nodes.indxs.surf_end;
    const int n_nodes = nodes.size();

    // Mark the cells
    for (VoronoiCell cell : voros) {
        if (cell.id <= cell_max)
            // Mark the cells that have at least one face on the surface
            for (VoronoiFace face : cell)
                if (vfaces.get_marker(face.id) == TYPES.SURFACE || vfaces.get_marker(face.id) == TYPES.ZMAX) {
                    voros.set_marker(cell.id, TYPES.SURFACE);
                    break;
                }
        // Mark the cells that for sure are not on the surface
        else voros.set_marker(cell.id, TYPES.PERIMETER);
    }

    // Mark the surface nodes
    nodes.init_markers(n_nodes);
    for (int i = 0; i < n_nodes; ++i)
        nodes.append_marker(voros.get_marker(i));
}

// Calculate minimum z-coordinate of medium and mark mesh
void VoronoiMesh::mark_mesh(const Medium& medium, const double latconst) {
    const int n_atoms = medium.size();
    require(n_atoms > 0, "Can't mark Voronoi mesh if no surface atoms are present!");
    const int n_average = sqrt(n_atoms);

    // Get the sorted array of the atom z-coordinates
    vector<double> zcoord; zcoord.reserve(n_atoms);
    for (int i = 0; i < n_atoms; ++i)
        zcoord.push_back(medium.get_point(i).z);

    sort(zcoord.begin(), zcoord.end());

    // Calculate the average z-coordinate of the lowest n_average atoms
    double zmin(0);
    for (int i = 0; i < n_average; ++i)
        zmin += zcoord[i];
    zmin /= n_average;
    zmin += 0.5 * latconst;

    mark_mesh(zmin);
}

// Mark the Voronoi cells, faces and nodes by their relative location against top and bottom surfaces
void VoronoiMesh::mark_mesh(const double zmin) {
    // Mark the cell and faces that are certainly on the surface
    int seedface = mark_seed();

    // Mark the rest of surface faces
//    mark_faces(zmin, seedface);
    mark_faces_vol2(zmin, seedface);

    // Using the marked faces, marks also the cells and nodes
    mark_cells_and_nodes();
}

// Extract the atoms and their areas whose Voronoi cells are exposed to vacuum
void VoronoiMesh::extract_surface(Medium& surface, vector<Vec3>& areas, const Medium& nanotip) {
    vector<bool> on_surface = vector_equal(voros.get_markers(), TYPES.SURFACE);
    const int n_substrate = surface.size();
    const int n_atoms = n_substrate + vector_sum(on_surface);
    
    areas.clear();
    areas.reserve(n_atoms);
    surface.resize(n_atoms);
    
    for (int i = 0; i < n_substrate; ++i)
        areas.push_back(Vec3(0));

    for (int i = 0; i <= nodes.indxs.surf_end; ++i)
        if (on_surface[i]) {           
            Vec3 centre = nodes.get_vec(i);
            Vec3 cell_area(0);

            // calculate the total area of the exposed part of cell
            for (VoronoiFace face : voros[i])
                if (vfaces.get_marker(face.id) == TYPES.SURFACE) {
                    Vec3 face_area = nodes.get_vec(face.nborcell(i)) - centre;
                    face_area.normalize();
                    face_area *= face.area();
                    cell_area += face_area;
                }

            areas.push_back(cell_area);
            surface.append(nanotip.get_atom(i));
        }
}

// Function to perform tetgen calculation on input and write_tetgen data
int VoronoiMesh::recalc(const string& cmd1, const string& cmd2, const string& cmd3) {
    try {
        tetgenio tetIOtemp1, tetIOtemp2;
        tetrahedralize(const_cast<char*>(cmd1.c_str()), &tetIOin, &tetIOtemp1);
        tetrahedralize(const_cast<char*>(cmd2.c_str()), &tetIOtemp1, &tetIOtemp2);
        tetrahedralize(const_cast<char*>(cmd3.c_str()), &tetIOtemp2, &tetIOout);
        nodes.set_counter(tetIOout.numberofpoints);
        elems.set_counter(tetIOout.numberoftetrahedra);
    }
    catch (int e) { return e; }
    return 0;
}

// Function to perform tetgen calculation on input and write_tetgen data
int VoronoiMesh::recalc(const string& cmd1, const string& cmd2) {
    try {
        tetgenio tetIOtemp1, tetIOtemp2;
        tetrahedralize(const_cast<char*>(cmd1.c_str()), &tetIOin, &tetIOtemp1);
        tetrahedralize(const_cast<char*>(cmd2.c_str()), &tetIOtemp1, &tetIOout);
        nodes.set_counter(tetIOout.numberofpoints);
        elems.set_counter(tetIOout.numberoftetrahedra);
    }
    catch (int e) { return e; }
    return 0;
}

int VoronoiMesh::recalc(const string& cmd1) {
    try {
        tetgenio tetIOtemp1, tetIOtemp2;
        tetrahedralize(const_cast<char*>(cmd1.c_str()), &tetIOin, &tetIOout);
        nodes.set_counter(tetIOout.numberofpoints);
        elems.set_counter(tetIOout.numberoftetrahedra);
    }
    catch (int e) { return e; }
    return 0;
}

// Generate Voronoi cells around surface atoms
int VoronoiMesh::generate_modi(const Medium& surf, const double latconst, const string& cmd1, const string& cmd2) {
    const double l = 1.0*latconst;

    Medium bulk(4), vacuum(4);
    bulk.append( Point3(surf.sizes.xmin-l, surf.sizes.ymin-l, surf.sizes.zmin-l) );
    bulk.append( Point3(surf.sizes.xmax+l, surf.sizes.ymin-l, surf.sizes.zmin-l) );
    bulk.append( Point3(surf.sizes.xmin-l, surf.sizes.ymax+l, surf.sizes.zmin-l) );
    bulk.append( Point3(surf.sizes.xmax+l, surf.sizes.ymax+l, surf.sizes.zmin-l) );

    vacuum.append( Point3(surf.sizes.xmin-l, surf.sizes.ymin-l, surf.sizes.zmax+l) );
    vacuum.append( Point3(surf.sizes.xmax+l, surf.sizes.ymin-l, surf.sizes.zmax+l) );
    vacuum.append( Point3(surf.sizes.xmin-l, surf.sizes.ymax+l, surf.sizes.zmax+l) );
    vacuum.append( Point3(surf.sizes.xmax+l, surf.sizes.ymax+l, surf.sizes.zmax+l) );

    const int n_bulk = bulk.size();
    const int n_surf = surf.size();
    const int n_vacuum = vacuum.size();
    nodes.init(n_bulk + n_surf + n_vacuum);

    // Add surface atoms first,...
    for (int i = 0; i < n_surf; ++i)
        nodes.append(surf.get_point(i));

    // ... bulk atoms second,...
    for (int i = 0; i < n_bulk; ++i)
        nodes.append(bulk.get_point(i));

    // ... and vacuum atoms last
    for (int i = 0; i < n_vacuum; ++i)
        nodes.append(vacuum.get_point(i));

    // Memorize the positions of different types of nodes to speed up later calculations
    nodes.save_indices(n_surf, n_bulk, n_vacuum);

    const int err_code = recalc("Q", cmd1);
    if (err_code) return err_code;

    elems.write("out/voro_tets1.vtk");

    vector<bool> is_surf; is_surf.reserve(elems.size());

    const int node_max = nodes.indxs.surf_end;
    for (SimpleElement elem : elems) {
        int ns = 0;
        int nv = 0;

        for (int e : elem) {
            if (e > node_max)
                nv++;
            else
                ns++;
        }

        is_surf.push_back(nv<=2 && ns>=2);
    }


    nodes.init(nodes.size() + vector_sum(is_surf));
    elems.init(elems.size());

    for (SimpleElement elem : elems)
        elems.append(elem);

    for (Point3 node : nodes)
        nodes.append(node);

    for (int i = 0; i < elems.size(); ++i)
        if (is_surf[i])
            nodes.append(elems.get_centroid(i));

    return recalc("rQ", cmd2);
//    return recalc("Q", cmd1, cmd2);
}


int VoronoiMesh::generate(const Medium& surf, const double latconst, const string& cmd1, const string& cmd2) {
    require(surf.size() > 0, "Invalid # generator points: " + to_string(surf.size()));
    const double l = 1.0*latconst;

    Medium bulk(4), vacuum(4);
    bulk.append( Point3(surf.sizes.xmin-l, surf.sizes.ymin-l, surf.sizes.zmin-l) );
    bulk.append( Point3(surf.sizes.xmax+l, surf.sizes.ymin-l, surf.sizes.zmin-l) );
    bulk.append( Point3(surf.sizes.xmin-l, surf.sizes.ymax+l, surf.sizes.zmin-l) );
    bulk.append( Point3(surf.sizes.xmax+l, surf.sizes.ymax+l, surf.sizes.zmin-l) );

    vacuum.append( Point3(surf.sizes.xmin-l, surf.sizes.ymin-l, surf.sizes.zmax+l) );
    vacuum.append( Point3(surf.sizes.xmax+l, surf.sizes.ymin-l, surf.sizes.zmax+l) );
    vacuum.append( Point3(surf.sizes.xmin-l, surf.sizes.ymax+l, surf.sizes.zmax+l) );
    vacuum.append( Point3(surf.sizes.xmax+l, surf.sizes.ymax+l, surf.sizes.zmax+l) );

    const int n_bulk = bulk.size();
    const int n_surf = surf.size();
    const int n_vacuum = vacuum.size();
    nodes.init(n_bulk + n_surf + n_vacuum);

    // Add surface atoms first,...
    for (int i = 0; i < n_surf; ++i)
        nodes.append(surf.get_point(i));

    // ... bulk atoms second,...
    for (int i = 0; i < n_bulk; ++i)
        nodes.append(bulk.get_point(i));

    // ... and vacuum atoms last
    for (int i = 0; i < n_vacuum; ++i)
        nodes.append(vacuum.get_point(i));

    // Memorize the positions of different types of nodes to speed up later calculations
    nodes.save_indices(n_surf, n_bulk, n_vacuum);

    return recalc("Q", cmd1, cmd2);
}

// Mark the cells and faces with nodes in the infinity
void VoronoiMesh::clean() {
    // Check for the infinite nodes on voro faces
    for (int f = 0; f < tetIOout.numberofvfacets; ++f) {
        tetgenio::vorofacet facet = tetIOout.vfacetlist[f];
        int& n_edges = tetIOout.vfacetlist[f].elist[0];
        for (int j = 1; j <= n_edges; ++j) {
            tetgenio::voroedge edge = tetIOout.vedgelist[facet.elist[j]];
            if (edge.v1 < 0 || edge.v2 < 0) { n_edges = 0; break; }
        }
    }

    // Check for the incomplete faces on voro cells
    for (int cell = 0; cell < tetIOout.numberofvcells; ++cell ) {
        int& n_faces = tetIOout.vcelllist[cell][0];
        for (int f = 1; f <= n_faces; ++f) {
            tetgenio::vorofacet facet = tetIOout.vfacetlist[tetIOout.vcelllist[cell][f]];
            int n_edges = facet.elist[0];
            if (n_edges <= 0) { n_faces = 0; break; }
        }
    }

    // Check for the unconnected faces
    for (int f = 0; f < tetIOout.numberofvfacets; ++f) {
        tetgenio::vorofacet facet = tetIOout.vfacetlist[f];
        int c1 = tetIOout.vcelllist[facet.c1][0];
        int c2 = tetIOout.vcelllist[facet.c2][0];
        if (c1 <= 0 && c2 <= 0) { tetIOout.vfacetlist[f].elist[0] = 0; }
    }

    voros.init_markers();
    vfaces.init_markers();
}

} /* namespace femocs */
