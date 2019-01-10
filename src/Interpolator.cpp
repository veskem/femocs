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

Interpolator::Interpolator(const string& vl, const string& nl, const string& sl) :
    nodes(vl, nl, sl),
    lintet(&nodes), lintri(&nodes, &lintet),
    quadtet(&nodes, &lintet), quadtri(&nodes, &lintri, &quadtet),
    linhex(&nodes, &lintet), linquad(&nodes, &lintri, &linhex),
    mesh(NULL), empty_value(0) 
    {}

void Interpolator::initialize(const TetgenMesh* m, double empty_val, int search_region) {
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
    nodes.precompute(search_region);
    lintri.precompute();
    quadtri.precompute();
    lintet.precompute();
    quadtet.precompute();
    linquad.precompute();
    linhex.precompute(search_region);

    lintet.narrow_search_to(search_region);
    empty_value = empty_val;

    const int n_nodes = nodes.size();
    const int n_hexs = mesh->hexs.size();

    expect(n_nodes > 0, "Interpolator expects non-empty mesh!");

    // === initialize solution values
    for (int i = 0; i < n_nodes; ++i)
        nodes.append_solution(Solution(empty_val));

    // === store mapping from mesh nodes to hexahedra & its node
    node2cells = vector<vector<pair<int,int>>>(n_nodes);
    if (search_region == TYPES.VACUUM) {
        for (int hex = 0; hex < n_hexs; ++hex)
            if (mesh->hexs.get_marker(hex) > 0) {
                SimpleHex shex = mesh->hexs[hex];
                for (int node = 0; node < n_nodes_per_hex; ++node)
                    node2cells[shex[node]].push_back( make_pair(hex, node) );
            }
    } else {
        for (int hex = 0; hex < n_hexs; ++hex)
            if (mesh->hexs.get_marker(hex) < 0) {
                SimpleHex shex = mesh->hexs[hex];
                for (int node = 0; node < n_nodes_per_hex; ++node)
                    node2cells[shex[node]].push_back( make_pair(hex, node) );
            }
    }
}

void Interpolator::store_solution(const vector<dealii::Tensor<1, 3>> vec_data, const vector<double> scal_data) {
    require(vec_data.size() == scal_data.size(), "Mismatch of vector sizes: "
            + d2s(vec_data.size()) + " vs " + d2s(scal_data.size()));

    const int n_femocs_nodes = nodes.size();
    const int n_dealii_nodes = vec_data.size();

    unsigned int j = 0;
    for (unsigned int i = 0; i < n_femocs_nodes; ++i) {
        // If there is a common node between Femocs and deal.II meshes, store actual solution
        if (nodes.femocs2deal(i) >= 0) {
            require(j < n_dealii_nodes, "Invalid index: " + d2s(j) + ">=" + d2s(n_dealii_nodes));
            dealii::Tensor<1, 3> vec = vec_data[j]; // that step needed to avoid complaints from Valgrind
            nodes.set_solution(i, Solution(Vec3(vec[0], vec[1], vec[2]), scal_data[j++]) );
        }
        else
            nodes.set_solution(i, Solution(empty_value));
    }
}

void Interpolator::store_solution(DealSolver<3>& solver) {
    vector<double> scalars;
    solver.export_solution(scalars);

    const int n_femocs_nodes = nodes.size();
    const int n_dealii_nodes = scalars.size();

    unsigned int j = 0;
    for (int i = 0; i < n_femocs_nodes; ++i)
        if (nodes.femocs2deal(i) >= 0) {
            require(j < n_dealii_nodes, "Invalid index: " + d2s(j) + ">=" + d2s(n_dealii_nodes));
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

bool Interpolator::average_nodal_fields(const bool vacuum) {
    require(mesh, "NULL mesh can't be used!");

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

void Interpolator::extract_solution(PoissonSolver<3>& fem, const bool smoothen) {
    // Read and store electric potentials from FEM solver
    store_solution(fem);

    // calculate field in the location of mesh nodes by calculating minus gradient of potential
#pragma omp parallel for
    for (int node = 0; node < nodes.size(); ++node) {
        if (nodes.femocs2deal(node) >= 0)
            store_elfield(node);
    }

    // Remove the spikes from the solution
    if (smoothen) average_nodal_fields(true);
}

void Interpolator::extract_solution_old(PoissonSolver<3>& fem, const bool smoothen) {
    // Read electric potentials and fields from FEM solver
    vector<double> potentials;
    vector<dealii::Tensor<1,3>> fields;
    fem.export_solution(potentials);
    fem.export_solution_grad(fields);

    // Store them
    store_solution(fields, potentials);

    // Remove the spikes from the solution
    if (smoothen) average_nodal_fields(true);
}

void Interpolator::extract_charge_density(PoissonSolver<3>& fem) {
    vector<double> charge_dens;
    vector<dealii::Tensor<1,3>> fields;
    fem.export_solution_grad(fields);
    fem.export_charge_dens(charge_dens);

    store_solution(fields, charge_dens);
}

void Interpolator::extract_solution(CurrentHeatSolver<3>& fem) {
    // Read current densities and temperatures from FEM solver
    vector<double> potentials, temperatures;
    vector<dealii::Tensor<1,3>> rhos;
    fem.current.export_solution(potentials);
    fem.export_temp_rho(temperatures, rhos);

    // Store them
    const int n_femocs_nodes = nodes.size();
    const int n_dealii_nodes = potentials.size();

    unsigned int j = 0;
    for (unsigned int i = 0; i < n_femocs_nodes; ++i) {
        if (nodes.femocs2deal(i) >= 0) {
            require(j < n_dealii_nodes, "Invalid index: " + d2s(j) + ">=" + d2s(n_dealii_nodes));
            dealii::Tensor<1, 3> rho = rhos[j]; // that step needed to avoid complaints from Valgrind
            nodes.set_solution(i, Solution( Vec3(rho[0], rho[1], rho[2]), potentials[j], temperatures[j] ));
            j++;
        }
        else
            nodes.set_solution(i, Solution(Vec3(0), 0, empty_value));
    }
}

} // namespace femocs
