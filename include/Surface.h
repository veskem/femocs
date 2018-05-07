/*
 * Media.h
 *
 *  Created on: 5.4.2016
 *      Author: veske
 */

#ifndef SURFACE_H_
#define SURFACE_H_

#include "Config.h"
#include "Macros.h"
#include "Coarseners.h"
#include "Medium.h"
#include "Interpolator.h"
#include "VoronoiMesh.h"

using namespace std;
namespace femocs {

/** Routines and data related to making a Surface */
class Surface: public Medium {
public:
    /** Initialise empty surface */
    Surface();

    /** Initialise empty surface and reserve memory for atoms */
    Surface(const int n_atoms);

    /** Generate surface with 4 atoms in simubox lateral corners */
    Surface(const Medium::Sizes& sizes, const double z);

    /** Generate surface with regular atom distribution along surface edges */
    Surface(const Medium::Sizes& sizes, const double z, const double dist);

    /** Pick suitable method for extending Surface */
    void extend(Surface& extension, const Config& conf);

    /** Generate nodal data that can be used as mesh generators */
    int generate_boundary_nodes(Surface& bulk, Surface& coarse_surf, Surface& vacuum,
            const Surface& extended_surf, const Config& conf, const bool first_time);

    /** Remove the atoms that are too far from surface faces */
    void clean_by_triangles(vector<int>& surf2face, Interpolator& interpolator, const TetgenMesh* mesh, const double r_cut);

private:
    Coarseners coarseners;  ///< atomistic coarsening data & routines

    /** Function used to smoothen the atoms */
    inline double smooth_function(const double distance, const double smooth_factor) const;

    /** Separate cylindrical region from substrate region */
    void get_nanotip(Surface& nanotip, const double radius);

    /** Extend the flat area by generating additional atoms */
    void extend(Surface &extension, const double latconst, const double box_width);

    /** Extend the flat area by reading additional atoms */
    void extend(Surface& extension, const string &file_name);

    /** Coarsen the atoms by generating additional boundary nodes and then running cleaner */
    void coarsen(Surface& surface);

    /** Coarsen the surface by using the linked list.
     * Cut-off radius is taken from the size of provided system. */
    void fast_coarsen(Surface &surface, const Medium::Sizes &s);

    /** Clean the surface from atoms that are too close to each other */
    void clean(Surface& surface);

    /** Clean atoms inside the region of interest and add the result to the provided surface */
    void add_cleaned_roi_to(Surface& surface);

    /** Smoothen the atoms inside the cylinder */
    void smoothen(const double radius, const double smooth_factor, const double r_cut);

    /** Smoothen all the atoms in the system */
    void smoothen(const double smooth_factor, const double r_cut);
};

} /* namespace femocs */

#endif /* SURFACE_H_ */
