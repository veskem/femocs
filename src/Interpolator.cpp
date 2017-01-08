/*
 * Interpolator.cpp
 *
 *  Created on: 19.10.2016
 *      Author: veske
 */

#include "Interpolator.h"
#include "Macros.h"

#include <float.h>

using namespace std;
namespace femocs {

// Initialize data vectors
Interpolator::Interpolator() {
    reserve(0);
    reserve_precompute(0);
};

/* Return the mapping between tetrahedral & hexahedral mesh nodes,
   nodes & hexahedral elements and nodes & element's vertices  */
const void Interpolator::get_maps(DealII& fem, vector<int>& tet2hex, vector<int>& node2hex, vector<int>& node2vert) {
    const int n_tet_nodes = get_n_atoms();
    const int n_hex_nodes = fem.triangulation.n_used_vertices();

    tet2hex.resize(n_tet_nodes, -1);
    node2hex.resize(n_hex_nodes);
    node2vert.resize(n_hex_nodes);

    // Loop through the hexahedral mesh vertices
    typename Triangulation<DIM>::active_vertex_iterator vertex = fem.triangulation.begin_active_vertex();
    for (int j = 0; j < n_tet_nodes; ++j, ++vertex)
        // Loop through tetrahedral mesh vertices
        for (int i = 0; i < n_tet_nodes; ++i)
            if ( (tet2hex[i] < 0) && (get_point(i) == vertex->vertex()) ) {
                tet2hex[i] = vertex->vertex_index();
                break;
            }

    // Loop through the hexahedral mesh elements
    typename DoFHandler<DIM>::active_cell_iterator cell;
    for (cell = fem.dof_handler.begin_active(); cell != fem.dof_handler.end(); ++cell)
        // Loop through all the vertices in the element
        for (int i = 0; i < fem.n_verts_per_elem; ++i) {
            node2hex[cell->vertex_index(i)] = cell->active_cell_index();
            node2vert[cell->vertex_index(i)] = i;
        }
}

// Extract the electric potential and electric field values on tetrahedral mesh nodes from FEM solution
const void Interpolator::extract_solution(DealII &fem, const TetgenMesh &mesh) {
    const int n_nodes = mesh.nodes.size();

    // Precompute tetrahedra to make interpolation faster
    precompute_tetrahedra(mesh);

    // Copy the mesh nodes
    reserve(n_nodes);
    for (int node = 0; node < n_nodes; ++node)
        add_atom( Atom(node, mesh.nodes[node], -1) );

    // To make solution extraction faster, generate mapping between desired and available data sequences
    vector<int> tet2hex, node2hex, node2vert;
    get_maps(fem, tet2hex, node2hex, node2vert);

    // Generate lists of hexahedra and hexahedra nodes where the tetrahedra nodes are located
    vector<int> cell_indxs; cell_indxs.reserve(n_nodes);
    vector<int> vert_indxs; vert_indxs.reserve(n_nodes);
    for (int n : tet2hex)
        if (n >= 0) {
            cell_indxs.push_back(node2hex[n]);
            vert_indxs.push_back(node2vert[n]);
        }

    vector<Vec3> ef = fem.get_elfield(cell_indxs, vert_indxs);      // get list of electric fields
    vector<double> pot = fem.get_potential(cell_indxs, vert_indxs); // get list of electric potentials

    require( ef.size() == pot.size(), "Mismatch of vector sizes: "
            + to_string(ef.size())  + ", " + to_string(pot.size()) );

    int i = 0;
    for (int node = 0; node < n_nodes; ++node) {
        // If there is a common node between tet and hex meshes, store actual solution
        if (tet2hex[node] >= 0) {
            require(i < ef.size(), "Invalid index: " + to_string(i));
            solution.push_back( Solution(ef[i], ef[i].norm(), pot[i]) );
            i++;
        }

        // In case of non-common node, store solution with error value
        else
            solution.push_back( Solution(error_field) );
    }
}

/* Return the mapping between tetrahedral & hexahedral mesh nodes,
   nodes & hexahedral elements and nodes & element's vertices  */
const void Interpolator::get_maps(dealii::Triangulation<3>* tria, dealii::DoFHandler<3>* dofh,
        vector<int>& tet2hex, vector<int>& cell_indxs, vector<int>& vert_indxs) {
    
    require(tria->n_vertices() > 0, "Invalid triangulation size!");
    
    const int n_tet_nodes = get_n_atoms();
    const int n_hex_nodes = tria->n_used_vertices();
    const int n_verts_per_elem = dealii::GeometryInfo<3>::vertices_per_cell;
    
    vector<int> node2hex, node2vert;
    tet2hex.resize(n_tet_nodes, -1);
    node2hex.resize(n_hex_nodes);
    node2vert.resize(n_hex_nodes);

    // Loop through the hexahedral mesh vertices
    typename Triangulation<DIM>::active_vertex_iterator vertex = tria->begin_active_vertex();
    for (int j = 0; j < n_tet_nodes; ++j, ++vertex)
        // Loop through tetrahedral mesh vertices
        for (int i = 0; i < n_tet_nodes; ++i)
            if ( (tet2hex[i] < 0) && (get_point(i) == vertex->vertex()) ) {
                tet2hex[i] = vertex->vertex_index();
                break;
            }
                
    // Loop through the hexahedral mesh elements
    typename DoFHandler<DIM>::active_cell_iterator cell;
    for (cell = dofh->begin_active(); cell != dofh->end(); ++cell)
        // Loop through all the vertices in the element
        for (int i = 0; i < n_verts_per_elem; ++i) {
            node2hex[cell->vertex_index(i)] = cell->active_cell_index();
            node2vert[cell->vertex_index(i)] = i;
        }
       
    // Generate lists of hexahedra and hexahedra nodes where the tetrahedra nodes are located
    cell_indxs.reserve(n_tet_nodes);
    vert_indxs.reserve(n_tet_nodes);
    for (int n : tet2hex)
        if (n >= 0) {
            cell_indxs.push_back(node2hex[n]);
            vert_indxs.push_back(node2vert[n]);
        }     
}

// Extract the electric potential and electric field values on tetrahedral mesh nodes from FEM solution
const void Interpolator::extract_solution(fch::Laplace<3>* fem, const TetgenMesh &mesh) {
    require(fem, "NULL pointer can't be handled!");
    
    // Copy the mesh nodes
    int i = 0;
    reserve(mesh.nodes.size());
    for (Point3 node : mesh.nodes)
        add_atom( Atom(i++, node, -1) );

    // Precompute tetrahedra to make interpolation faster
    precompute_tetrahedra(mesh);
    
    // To make solution extraction faster, generate mapping between desired and available data sequences
    vector<int> tet_node2hex_node, cell_indxs, vert_indxs;
    get_maps(fem->get_triangulation(), fem->get_dof_handler(), tet_node2hex_node, cell_indxs, vert_indxs);

    vector<dealii::Tensor<1,3>> ef = fem->get_efield(cell_indxs, vert_indxs); // get list of electric fields
    vector<double> pot = fem->get_potential(cell_indxs, vert_indxs);          // get list of electric potentials

    require( ef.size() == pot.size(), "Mismatch of vector sizes: "
            + to_string(ef.size())  + ", " + to_string(pot.size()) );

    i = 0;
    for (int n : tet_node2hex_node) {
        // If there is a common node between tet and hex meshes, store actual solution
        if (n >= 0) {
            require(i < ef.size(), "Invalid index: " + to_string(i));
            Vec3 elfield(ef[i][0], ef[i][1], ef[i][2]);
            solution.push_back( Solution(elfield, elfield.norm(), pot[i++]) );
        }

        // In case of non-common node, store solution with error value
        else solution.push_back( Solution(error_field) );
    }
}

// Extract the electric potential and electric field values on tetrahedral mesh nodes from FEM solution
const void Interpolator::extract_solution(fch::CurrentsAndHeating<3>* fem, const TetgenMesh &mesh) {  
    require(fem, "NULL pointer can't be handled!");
    
    // Copy the mesh nodes
    int i = 0;
    reserve(mesh.nodes.size());
    for (Point3 node : mesh.nodes)
        add_atom( Atom(i++, node, -1) );

    // Precompute tetrahedra to make interpolation faster
    precompute_tetrahedra(mesh);
    
    // To make solution extraction faster, generate mapping between desired and available data sequences
    vector<int> tet_node2hex_node, cell_indxs, vert_indxs;
    get_maps(fem->get_triangulation(), fem->get_dof_handler(), tet_node2hex_node, cell_indxs, vert_indxs);
    
    vector<dealii::Tensor<1,3>> rho = fem->get_current(cell_indxs, vert_indxs); // get list of current densities
    vector<double> temperature = fem->get_temperature(cell_indxs, vert_indxs); // get list of temperatures

    require( rho.size() == temperature.size(), "Mismatch of vector sizes: "
            + to_string(rho.size())  + ", " + to_string(temperature.size()) );
    
    i = 0;
    for (int n : tet_node2hex_node) {
        // If there is a common node between tet and hex meshes, store actual solution
        if (n >= 0) {
            require(i < rho.size(), "Invalid index: " + to_string(i));
            Vec3 current(rho[i][0], rho[i][1], rho[i][2]);
            solution.push_back( Solution(current, current.norm(), temperature[i++]) );
        }

        // In case of non-common node, store solution with error value
        else solution.push_back( Solution(error_field) );
    }
}

// Reserve memory for interpolation data
const void Interpolator::reserve(const int N) {
    require(N >= 0, "Invalid number of points: " + to_string(N));
    atoms.clear();
    solution.clear();

    atoms.reserve(N);
    solution.reserve(N);
}

// Reserve memory for pre-compute data
const void Interpolator::reserve_precompute(const int N) {
    tetrahedra.clear();
    tetneighbours.clear();
    centroid.clear();
    det0.clear();
    det1.clear();
    det2.clear();
    det3.clear();
    det4.clear();
    tet_not_valid.clear();

    tetrahedra.reserve(N);
    tetneighbours.reserve(N);
    centroid.reserve(N);
    det0.reserve(N);
    det1.reserve(N);
    det2.reserve(N);
    det3.reserve(N);
    det4.reserve(N);
    tet_not_valid.reserve(N);
}

// Precompute the data to tetrahedra to make later bcc calculations faster
const void Interpolator::precompute_tetrahedra(const TetgenMesh &mesh) {
    const int n_elems = mesh.elems.size();
    double d0, d1, d2, d3, d4;

    reserve_precompute(n_elems);

    // Copy tetrahedra
    for (int i = 0; i < n_elems; ++i)
        tetrahedra.push_back(mesh.elems[i]);

    // Copy tetrahedra neighbours
    for (int i = 0; i < n_elems; ++i)
        tetneighbours.push_back(mesh.elems.get_neighbours(i));

    // Calculate centroids of elements
    for (int i = 0; i < n_elems; ++i)
        centroid.push_back(mesh.elems.get_centroid(i));

    /* Calculate main and minor determinants for 1st, 2nd, 3rd and 4th
     * barycentric coordinate of tetrahedra using the relations below */
    for (SimpleElement se : tetrahedra) {
        Vec3 v1 = mesh.nodes.get_vec(se[0]);
        Vec3 v2 = mesh.nodes.get_vec(se[1]);
        Vec3 v3 = mesh.nodes.get_vec(se[2]);
        Vec3 v4 = mesh.nodes.get_vec(se[3]);

        /* =====================================================================================
         * det0 = |x1 y1 z1 1|
                  |x2 y2 z2 1|
                  |x3 y3 z3 1|
                  |x4 y4 z4 1|  */
        d0 = determinant(v1, v2, v3, v4);
        tet_not_valid.push_back(fabs(d0) < epsilon);
        det0.push_back(1.0 / d0);

        /* =====================================================================================
         * det1 = |x  y  z  1| = + x * |y2 z2 1| - y * |x2 z2 1| + z * |x2 y2 1| - |x2 y2 z2|
                  |x2 y2 z2 1|         |y3 z3 1|       |x3 z3 1|       |x3 y3 1|   |x3 y3 z3|
                  |x3 y3 z3 1|         |y4 z4 1|       |x4 z4 1|       |x4 y4 1|   |x4 y4 z4|
                  |x4 y4 z4 1|  */
        d1 = determinant(Vec3(v2.y, v3.y, v4.y), Vec3(v2.z, v3.z, v4.z));
        d2 = determinant(Vec3(v2.x, v3.x, v4.x), Vec3(v2.z, v3.z, v4.z));
        d3 = determinant(Vec3(v2.x, v3.x, v4.x), Vec3(v2.y, v3.y, v4.y));
        d4 = determinant(Vec3(v2.x, v3.x, v4.x), Vec3(v2.y, v3.y, v4.y), Vec3(v2.z, v3.z, v4.z));

        det1.push_back(Vec4(d1, -d2, d3, -d4));

        /* =====================================================================================
         * det2 = |x1 y1 z1 1| = - x * |y1 z1 1| + y * |x1 z1 1| - z * |x1 y1 1| + |x1 y1 z1|
                  |x  y  z  1|         |y3 z3 1|       |x3 z3 1|       |x3 y3 1|   |x3 y3 z3|
                  |x3 y3 z3 1|         |y4 z4 1|       |x4 z4 1|       |x4 y4 1|   |x4 y4 z4|
                  |x4 y4 z4 1|  */
        d1 = determinant(Vec3(v1.y, v3.y, v4.y), Vec3(v1.z, v3.z, v4.z));
        d2 = determinant(Vec3(v1.x, v3.x, v4.x), Vec3(v1.z, v3.z, v4.z));
        d3 = determinant(Vec3(v1.x, v3.x, v4.x), Vec3(v1.y, v3.y, v4.y));
        d4 = determinant(Vec3(v1.x, v3.x, v4.x), Vec3(v1.y, v3.y, v4.y), Vec3(v1.z, v3.z, v4.z));

        det2.push_back(Vec4(-d1, d2, -d3, d4));

        /* =====================================================================================
         * det3 = |x1 y1 z1 1| = + x * |y1 z1 1| - y * |x1 z1 1| + z * |x1 y1 1| - |x1 y1 z1|
                  |x2 y2 z2 1|         |y2 z2 1|       |x2 z2 1|       |x2 y2 1|   |x2 y2 z2|
                  |x  y  z  1|         |y4 z4 1|       |x4 z4 1|       |x4 y4 1|   |x4 y4 z4|
                  |x4 y4 z4 1|  */
        d1 = determinant(Vec3(v1.y, v2.y, v4.y), Vec3(v1.z, v2.z, v4.z));
        d2 = determinant(Vec3(v1.x, v2.x, v4.x), Vec3(v1.z, v2.z, v4.z));
        d3 = determinant(Vec3(v1.x, v2.x, v4.x), Vec3(v1.y, v2.y, v4.y));
        d4 = determinant(Vec3(v1.x, v2.x, v4.x), Vec3(v1.y, v2.y, v4.y), Vec3(v1.z, v2.z, v4.z));

        det3.push_back(Vec4(d1, -d2, d3, -d4));

        /* =====================================================================================
         * det4 = |x1 y1 z1 1| = - x * |y1 z1 1| + y * |x1 z1 1| - z * |x1 y1 1| + |x1 y1 z1|
                  |x2 y2 z2 1|         |y2 z2 1|       |x2 z2 1|       |x2 y2 1|   |x2 y2 z2|
                  |x3 y3 z3 1|         |y3 z3 1|       |x3 z3 1|       |x3 y3 1|   |x3 y3 z3|
                  |x  y  z  1|  */

        d1 = determinant(Vec3(v1.y, v2.y, v3.y), Vec3(v1.z, v2.z, v3.z));
        d2 = determinant(Vec3(v1.x, v2.x, v3.x), Vec3(v1.z, v2.z, v3.z));
        d3 = determinant(Vec3(v1.x, v2.x, v3.x), Vec3(v1.y, v2.y, v3.y));
        d4 = determinant(v1, v2, v3);

        det4.push_back(Vec4(-d1, d2, -d3, d4));
    }
}

// Check with barycentric coordinates whether the point is inside the i-th tetrahedron
const bool Interpolator::point_in_tetrahedron(const Point3 &point, const int i) {
    require(i >= 0 && i < det0.size(), "Index out of bounds: " + to_string(i));

    // Ignore co-planar tetrahedra
    // no need to check because Tetgen guarantees non-co-planar tetrahedra
//    if (tet_not_valid[i]) return false;

    Vec4 pt(point, 1);

    // If one of the barycentric coordinates is < zero, the point is outside the tetrahedron
    // Source: http://steve.hollasch.net/cgindex/geometry/ptintet.html
    if (det0[i] * pt.dotProduct(det1[i]) < zero) return false;
    if (det0[i] * pt.dotProduct(det2[i]) < zero) return false;
    if (det0[i] * pt.dotProduct(det3[i]) < zero) return false;
    if (det0[i] * pt.dotProduct(det4[i]) < zero) return false;

    // All bcc-s are >= 0, so point is inside the tetrahedron
    return true;
}

// Calculate barycentric coordinates for point
const Vec4 Interpolator::get_bcc(const Point3 &point, const int elem) {
    require(elem >= 0 && elem < det0.size(), "Index out of bounds: " + to_string(elem));

    Vec4 pt(point, 1);
    double bcc1 = det0[elem] * pt.dotProduct(det1[elem]);
    double bcc2 = det0[elem] * pt.dotProduct(det2[elem]);
    double bcc3 = det0[elem] * pt.dotProduct(det3[elem]);
    double bcc4 = det0[elem] * pt.dotProduct(det4[elem]);

    return Vec4(bcc1, bcc2, bcc3, bcc4);
}

// Find the element which contains the point or is the closest to it
const int Interpolator::locate_element(const Point3 &point, const int elem_guess) {
    const int n_elems = det0.size();
            
    // Check the guessed element
    if (point_in_tetrahedron(point, elem_guess)) return elem_guess;

    // Check the 1st nearest neighbours of guessed element
    for (int nbor : tetneighbours[elem_guess])
        if (nbor >= 0 && point_in_tetrahedron(point, nbor))
            return nbor;
   
    // Check the 2nd, 3rd & 4th nearest neighbours of guessed element
    // going further than 4th neighbour starts slowing things down
    
    // Mark the neighbour ranks of elements wrt to guessed element
    vector<int> nbor_rank(n_elems);
    for (int nbor1 : tetneighbours[elem_guess]) {
        if (nbor1 < 0) continue;
        for (int nbor2 : tetneighbours[nbor1]) {
            if (nbor2 < 0) continue;
            for (int nbor3 : tetneighbours[nbor2]) {
                if (nbor3 < 0) continue;
                for (int nbor4 : tetneighbours[nbor3]) {
                    if (nbor4 < 0) continue;
                    nbor_rank[nbor4] = 4;
                }
                nbor_rank[nbor3] = 3;
            }
            nbor_rank[nbor2] = 2;            
        }
        nbor_rank[nbor1] = 1;        
    }
    
    // Perform the check on the neighbours
    // checking ranks separately doesn't give any benefit in speed
    for (int elem = 0; elem < n_elems; ++elem)
        if ( nbor_rank[elem] > 1 && point_in_tetrahedron(point, elem) )
            return elem;
        
    // If no success, loop through all the elements
    double min_distance2 = DBL_MAX;
    int min_index = 0;
    
    for (int elem = 0; elem < n_elems; ++elem) {
        // If correct element is found, we're done
        if (point_in_tetrahedron(point, elem)) return elem;

        // Otherwise look for the element whose centroid is closest to the point
        else {
            double distance2 = point.distance2(centroid[elem]);
            if (distance2 < min_distance2) {
                min_distance2 = distance2;
                min_index = elem;
            }
        }
    }

    // If no perfect element found, return the best
    // indicate the imperfectness with the minus sign
    return -min_index;
}

// Calculate interpolation for point inside or near the elem-th tetrahedron
const Solution Interpolator::get_interpolation(const Point3 &point, const int elem) {
    require(elem >= 0 && elem < tetrahedra.size(), "Index out of bounds: " + to_string(elem));

    // Get barycentric coordinates of point in tetrahedron
    Vec4 bcc = get_bcc(point, elem);

    SimpleElement selem = tetrahedra[elem];

    // Interpolate electric field
    Vec3 elfield_i(0.0);
    for (int i = 0; i < selem.size(); ++i)
        elfield_i += solution[selem[i]].elfield * bcc[i];

    // Interpolate potential
    double potential_i(0.0);
    for (int i = 0; i < selem.size(); ++i)
        potential_i += solution[selem[i]].potential * bcc[i];

    return Solution(elfield_i, elfield_i.norm(), potential_i);
}

// Compile data string from the data vectors for file output
const string Interpolator::get_data_string(const int i) const {
    if (i < 0) return "Interpolator data: id x y z dummy Ex Ey Ez Enorm potential";

    ostringstream strs;
    strs << atoms[i] << " " << solution[i];
    return strs.str();
}

// Determinant of 3x3 matrix which's last column consists of ones
const double Interpolator::determinant(const Vec3 &v1, const Vec3 &v2) {
    return v1.x * (v2.y - v2.z)
         - v1.y * (v2.x - v2.z)
         + v1.z * (v2.x - v2.y);
}

// Determinant of 3x3 matrix which's columns consist of Vec3-s
const double Interpolator::determinant(const Vec3 &v1, const Vec3 &v2, const Vec3 &v3) {
    return v1.x * (v2.y * v3.z - v3.y * v2.z)
         - v2.x * (v1.y * v3.z - v3.y * v1.z)
         + v3.x * (v1.y * v2.z - v2.y * v1.z);
}

// Determinant of 4x4 matrix which's last column consists of ones
const double Interpolator::determinant(const Vec3 &v1, const Vec3 &v2, const Vec3 &v3, const Vec3 &v4) {
    const double det1 = determinant(v2, v3, v4);
    const double det2 = determinant(v1, v3, v4);
    const double det3 = determinant(v1, v2, v4);
    const double det4 = determinant(v1, v2, v3);

    return det4 - det3 + det2 - det1;
}

// Determinant of 4x4 matrix which's columns consist of Vec4-s
const double Interpolator::determinant(const Vec4 &v1, const Vec4 &v2, const Vec4 &v3, const Vec4 &v4) {
    double det1 = determinant(Vec3(v1.y,v1.z,v1.w), Vec3(v2.y,v2.z,v2.w), Vec3(v3.y,v3.z,v3.w));
    double det2 = determinant(Vec3(v1.x,v1.z,v1.w), Vec3(v2.x,v2.z,v2.w), Vec3(v3.x,v3.z,v3.w));
    double det3 = determinant(Vec3(v1.x,v1.y,v1.w), Vec3(v2.x,v2.y,v2.w), Vec3(v3.x,v3.y,v3.w));
    double det4 = determinant(Vec3(v1.x,v1.y,v1.z), Vec3(v2.x,v2.y,v2.z), Vec3(v3.x,v3.y,v3.z));

    return v4.w * det4 - v4.z * det3 + v4.y * det2 - v4.x * det1;
}

} // namespace femocs
