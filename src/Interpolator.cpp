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
    lintet(&nodes), lintri(&nodes, &lintet),
    quadtet(&nodes, &lintet), quadtri(&nodes, &lintri, &quadtet),
    linhex(&nodes, &lintet), linquad(&nodes, &lintri, &linhex),
    mesh(NULL), empty_value(0) 
    {}

bool Interpolator::average_nodal_fields(const bool vacuum) {
    vector<vector<unsigned int>> nborlist;
    mesh->calc_pseudo_3D_vorocells(nborlist, vacuum);

    // loop through the tetrahedral nodes
    for (unsigned int i = 0; i < nborlist.size(); ++i) {
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

void Interpolator::get_maps(vector<int>& cell_indxs, vector<int>& vert_indxs,
        dealii::Triangulation<3>* tria, dealii::DoFHandler<3>* dofh)
{
    require(tria->n_vertices() > 0, "Invalid triangulation size!");
    const int n_verts_per_elem = dealii::GeometryInfo<3>::vertices_per_cell;
    const int n_femocs_nodes = mesh->nodes.size();
    expect(n_femocs_nodes > 0, "Interpolator expects non-empty mesh!");
    const int n_dealii_nodes = tria->n_used_vertices();

    vector<int> node2hex(n_dealii_nodes), node2vert(n_dealii_nodes);

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

void Interpolator::store_solution(const vector<dealii::Tensor<1, 3>> vec_data, const vector<double> scal_data) {

    require(vec_data.size() == scal_data.size(), "Mismatch of vector sizes: "
            + d2s(vec_data.size()) + ", " + d2s(scal_data.size()));

    unsigned int j = 0;
    for (unsigned int i = 0; i < femocs2deal.size(); ++i) {
        // If there is a common node between Femocs and deal.II meshes, store actual solution
        if (femocs2deal[i] >= 0) {
            require(j < vec_data.size(), "Invalid index: " + d2s(j));
            dealii::Tensor<1, 3> vec = vec_data[j]; // that step needed to avoid complaints from Valgrind
            nodes.set_solution(i, Solution(Vec3(vec[0], vec[1], vec[2]), scal_data[j++]) );
        }
        else
            nodes.set_solution(i, Solution(empty_value));
    }
}

void Interpolator::store_solution(DealSolver<3>& solver) {
    const unsigned int n_nodes = nodes.size();
    require(femocs2deal.size() == n_nodes, "Invalid femocs2deal size: " + d2s(femocs2deal.size()));

    vector<double> scalars;
    solver.get_nodal_solution(scalars);

    unsigned int j = 0;
    for (unsigned int i = 0; i < n_nodes; ++i)
        if (femocs2deal[i] >= 0) {
            require(j < scalars.size(), "Invalid Femocs nodes to Deal.II nodes mapping!");
            nodes.set_scalar(i, scalars[j++]);
        }
}

void Interpolator::store_elfield(const int node) {
    // due to linear elements, field on a node must be averaged
    // over all the hexahedra that are connected to the node,
    // as the field is discontinuous between hexahedra

    Vec3 mean_field(0);

    int n_fields = node2cells[node].size();
    if (n_fields > 0) {
        for (pair<int,int> p : node2cells[node])
            mean_field += linhex.interp_gradient(p.first, p.second);
        mean_field *= (1.0 / n_fields);
    }

    nodes.set_vector(node, mean_field);
}

int Interpolator::update_point_cell(const SuperParticle& particle) const {
    int femocs_cell = linhex.deal2femocs(particle.cell);
    femocs_cell = linhex.locate_cell(particle.pos, femocs_cell);
    if (femocs_cell < 0) return -1;
    return linhex.femocs2deal(femocs_cell);
}

void Interpolator::initialize(const TetgenMesh* m, DealSolver<3>& solver, const double empty_val) {
    // === update mesh
    mesh = m;
    nodes.set_mesh(mesh);
    lintri.set_mesh(mesh);
    lintet.set_mesh(mesh);
    quadtri.set_mesh(mesh);
    quadtet.set_mesh(mesh);
    linquad.set_mesh(mesh);
    linhex.set_mesh(mesh);

    // === precompute cells to make interpolation faster
    nodes.precompute();
    lintri.precompute();
    quadtri.precompute();
    lintet.precompute();
    quadtet.precompute();
    linquad.precompute();
    linhex.precompute();

    empty_value = empty_val;

    const int n_nodes = nodes.size();
    const int n_hexs = mesh->hexs.size();
    static constexpr int n_nodes_per_hex = 8;

    expect(n_nodes > 0, "Interpolator expects non-empty mesh!");

    // === initialize solution values
    for (int i = 0; i < n_nodes; ++i)
        nodes.append_solution(Solution(empty_val));

    // === store mapping between mesh nodes and hexahedra
    node2cells = vector<vector<pair<int,int>>>(n_nodes);
    for (int hex = 0; hex < n_hexs; ++hex)
        if (mesh->hexs.get_marker(hex) > 0) {
            SimpleHex shex = mesh->hexs[hex];
            for (int node = 0; node < n_nodes_per_hex; ++node)
                node2cells[shex[node]].push_back( make_pair(hex, node) );
        }

    // === store mapping between Femocs and Deal.II mesh nodes

    const double vertex_epsilon = 1e-5 * mesh->tets.stat.edgemin;
    dealii::Triangulation<3>* tria = solver.get_triangulation();
    typename dealii::Triangulation<3>::active_vertex_iterator vertex = tria->begin_active_vertex();

    femocs2deal = vector<int>(n_nodes, -1);
    for (int i = 0; i < n_nodes && vertex != tria->end_vertex(); ++i)
        if ( mesh->nodes[i].distance2(vertex->vertex()) < vertex_epsilon ) {
            femocs2deal[i] = vertex->vertex_index();
            vertex++;
        }
}

void Interpolator::extract_charge_density(PoissonSolver<3>& fem) {

    // To make solution extraction faster, generate mapping between desired and available data sequences
    vector<int> cell_indxs, vert_indxs;
    get_maps(cell_indxs, vert_indxs, fem.get_triangulation(), fem.get_dof_handler());

    // Read and store the electric field and potential from FEM solver
    store_solution(fem.get_efield(cell_indxs, vert_indxs), fem.get_charge_dens(cell_indxs, vert_indxs));
}

void Interpolator::extract_solution(CurrentHeatSolver<3>& fem) {

    // To make solution extraction faster, generate mapping between desired and available data sequences
    vector<int> cell_indxs, vert_indxs;
    get_maps(cell_indxs, vert_indxs, fem.get_triangulation(), fem.current.get_dof_handler());

    // Read and store current densities and temperatures from FEM solver
    store_solution(fem.get_current(cell_indxs, vert_indxs), fem.get_temperature(cell_indxs, vert_indxs));
}

void Interpolator::extract_solution_v2(PoissonSolver<3>& fem, const bool smoothen) {
    // To make solution extraction faster, generate mapping between desired and available data sequences
    vector<int> cell_indxs, vert_indxs;
    get_maps(cell_indxs, vert_indxs, fem.get_triangulation(), fem.get_dof_handler());

    // Read and store current densities and temperatures from FEM solver
    store_solution(fem.get_efield(cell_indxs, vert_indxs), fem.get_potential(cell_indxs, vert_indxs));

    // Remove the spikes from the solution
    if (smoothen) average_nodal_fields(true);
}

void Interpolator::extract_solution(PoissonSolver<3>& fem, const bool smoothen) {
    store_solution(fem);

    // calculate field in the location of mesh nodes by calculating minus gradient of potential
#pragma omp parallel for
    for (int node = 0; node < nodes.size(); ++node) {
        if (femocs2deal[node] >= 0)
            store_elfield(node);
    }

    if (smoothen) average_nodal_fields(true);
}

} // namespace femocs
