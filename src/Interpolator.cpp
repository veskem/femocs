/*
 * Interpolator.cpp
 *
 *  Created on: 19.10.2016
 *      Author: veske
 */

#include "Interpolator.h"

#include "Macros.h"
#include "float.h"

using namespace std;
namespace femocs {

/* ==================================================================
 *  ======================== Interpolator ==========================
 * ================================================================== */

Interpolator::Interpolator(const string& nl, const string& sl) : 
    nodes(nl, sl),
    lintri(&nodes), lintet(&nodes),
    quadtri(&nodes, &lintri), quadtet(&nodes, &lintet),
    linquad(&nodes, &lintri), linhex(&nodes, &lintet),
    mesh(NULL), empty_value(0) 
    {}

// Force the solution on tetrahedral nodes to be the weighed average of the solutions on its
// surrounding hexahedral nodes
bool Interpolator::average_sharp_nodes(const bool vacuum) {

    return false;
    vector<vector<unsigned int>> nborlist;
    mesh->calc_pseudo_3D_vorocells(nborlist, vacuum);

    // loop through the tetrahedral nodes
    for (int i = 0; i < nborlist.size(); ++i) {
        if (nborlist[i].size() == 0) continue;

        Point3 tetnode = nodes.get_vertex(i);
        Vec3 vec(0);
        double w_sum = 0;

        // tetnode new solution will be the weighed average of the solutions on neighbouring nodes
        for (unsigned nbor : nborlist[i]) {
            double w = exp(lintet.decay_factor * tetnode.distance(nodes.get_vertex(nbor)));
            w_sum += w;
            vec += nodes.get_vector(nbor) * w;
        }

        if (w_sum > 0) {
            vec *= (1.0 / w_sum);
            nodes.set_solution(i, Solution(vec, nodes.get_scalar(i)));
        }
    }

    return false;
}

/* Calculate mapping between Femocs & deal.II mesh nodes,
 nodes & hexahedral elements and nodes & element's vertices */
void Interpolator::get_maps(vector<int>& femocs2deal, vector<int>& cell_indxs, vector<int>& vert_indxs,
        dealii::Triangulation<3>* tria, dealii::DoFHandler<3>* dofh) const {

    const int n_verts_per_elem = dealii::GeometryInfo<3>::vertices_per_cell;
    const double vertex_epsilon = 1e-5 * mesh->tets.stat.edgemin;

    require(tria->n_vertices() > 0, "Invalid triangulation size!");
    const int n_femocs_nodes = mesh->nodes.size();
    expect(n_femocs_nodes > 0, "Interpolator expects non-empty mesh!");
    const int n_dealii_nodes = tria->n_used_vertices();

    vector<int> node2hex(n_dealii_nodes), node2vert(n_dealii_nodes);
    femocs2deal = vector<int>(n_femocs_nodes, -1);

    typename dealii::Triangulation<3>::active_vertex_iterator vertex = tria->begin_active_vertex();
    // Loop through tetrahedral mesh vertices
    for (int i = 0; i < n_femocs_nodes && vertex != tria->end_vertex(); ++i)
        if ( mesh->nodes[i].distance2(vertex->vertex()) < vertex_epsilon ) {
            femocs2deal[i] = vertex->vertex_index();
            vertex++;
        }

    // Loop through the hexahedral mesh elements
    typename dealii::DoFHandler<3>::active_cell_iterator cell;
    for (cell = dofh->begin_active(); cell != dofh->end(); ++cell)
        // Loop through all the vertices in the element
        for (int i = 0; i < n_verts_per_elem; ++i) {
            node2hex[cell->vertex_index(i)] = cell->active_cell_index();
            node2vert[cell->vertex_index(i)] = i;
        }

    // Generate lists of hexahedra and hexahedra nodes where the tetrahedra nodes are located
    cell_indxs.reserve(n_femocs_nodes);
    vert_indxs.reserve(n_femocs_nodes);
    for (int n : femocs2deal)
        if (n >= 0) {
            cell_indxs.push_back(node2hex[n]);
            vert_indxs.push_back(node2vert[n]);
        }
}

void Interpolator::store_solution(const vector<int>& femocs2deal,
        const vector<dealii::Tensor<1, 3>> vec_data, const vector<double> scal_data) {

    require(vec_data.size() == scal_data.size(), "Mismatch of vector sizes: "
            + to_string(vec_data.size()) + ", " + to_string(scal_data.size()));

    int j = 0;
    for (int i = 0; i < femocs2deal.size(); ++i) {
        // If there is a common node between Femocs and deal.II meshes, store actual solution
        if (femocs2deal[i] >= 0) {
            require(j < vec_data.size(), "Invalid index: " + to_string(j));
            dealii::Tensor<1, 3> vec = vec_data[j]; // that step needed to avoid complaints from Valgrind
            nodes.set_solution(i, Solution(Vec3(vec[0], vec[1], vec[2]), scal_data[j++]) );
        }
        else
            nodes.set_solution(i, Solution(empty_value));
    }
}

int Interpolator::update_point_cell(const SuperParticle& particle) const {
    int femocs_cell = linhex.deal2femocs(particle.cell);
    femocs_cell = linhex.locate_cell(particle.pos, femocs_cell);
    if (femocs_cell < 0) return -1;
    return linhex.femocs2deal(femocs_cell);
}

void Interpolator::initialize(const TetgenMesh* m, const double empty_val) {
    // update mesh
    mesh = m;
    nodes.set_mesh(mesh);
    lintri.set_mesh(mesh);
    lintet.set_mesh(mesh);
    quadtri.set_mesh(mesh);
    quadtet.set_mesh(mesh);
    linquad.set_mesh(mesh);
    linhex.set_mesh(mesh);

    // Precompute cells to make interpolation faster
    nodes.precompute();
    lintri.precompute();
    quadtri.precompute();
    lintet.precompute();
    quadtet.precompute();
    linquad.precompute();
    linhex.precompute();

    empty_value = empty_val;

    const int n_atoms = nodes.size();
    for (int i = 0; i < n_atoms; ++i)
        nodes.append_solution(Solution(empty_val));
}

void Interpolator::extract_charge_density(fch::PoissonSolver<3>& fem) {

    // To make solution extraction faster, generate mapping between desired and available data sequences
    vector<int> femocs2deal, cell_indxs, vert_indxs;
    get_maps(femocs2deal, cell_indxs, vert_indxs, fem.get_triangulation(), fem.get_dof_handler());

    // Read and store the electric field and potential from FEM solver
    store_solution(femocs2deal, fem.get_efield(cell_indxs, vert_indxs),
            fem.get_charge_dens(cell_indxs, vert_indxs));
}

void Interpolator::extract_solution(fch::CurrentHeatSolver<3>& fem) {

    // To make solution extraction faster, generate mapping between desired and available data sequences
    vector<int> femocs2deal, cell_indxs, vert_indxs;
    get_maps(femocs2deal, cell_indxs, vert_indxs, fem.get_triangulation(), fem.current.get_dof_handler());

    // Read and store current densities and temperatures from FEM solver
    store_solution(femocs2deal, fem.get_current(cell_indxs, vert_indxs),
            fem.get_temperature(cell_indxs, vert_indxs));
}

void Interpolator::extract_solution(fch::PoissonSolver<3>& fem) {

    // To make solution extraction faster, generate mapping between desired and available data sequences
    vector<int> femocs2deal, cell_indxs, vert_indxs;
    get_maps(femocs2deal, cell_indxs, vert_indxs, fem.get_triangulation(), fem.get_dof_handler());

    // Read and store current densities and temperatures from FEM solver
    store_solution(femocs2deal, fem.get_efield(cell_indxs, vert_indxs),
            fem.get_potential(cell_indxs, vert_indxs));

    // Remove the spikes from the solution
    average_sharp_nodes(true);
}

} // namespace femocs
