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
#include "AtomReader.h"
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

    /** Generate surface with 4 atoms in simubox corners */
    Surface(const Medium::Sizes& sizes, const double z);

    /** Generate surface with regular atom distribution along surface edges */
    Surface(const Medium::Sizes& sizes, const double z, const double dist);

    /** Extract atom with desired types and clean it from lonely atoms;
     * atom is considered lonely if its coordination is lower than threshold.
     * @param reader  AtomReader holding the atoms and their types
     * @param type    type of the atoms that will be read
     * @param invert  if true all the atoms except the 'type'-ones will be stored
     */
    void extract(const AtomReader& reader, const int type, const bool invert=false);

    /** Extend the flat area by generating additional atoms */
    void extend(Surface &extension, Coarseners &cr, const double latconst, const double box_width);

    /** Extend the flat area by reading additional atoms */
    Surface extend(const string &file_name, Coarseners &coarseners);

    /** Coarsen the atoms by generating additional boundary nodes and then running cleaner */
    Surface coarsen(Coarseners &coarseners);

    /** Clean the surface from atoms that are too close to each other */
    Surface clean(Coarseners &coarseners);

    Surface clean_roi(Coarseners &coarseners);

    /** Remove the atoms that are too far from surface faces */
    void clean_by_triangles(vector<int>& surf2face, Interpolator& interpolator, const TetgenMesh* mesh, const double r_cut);

    /** Coarsen the surface by using the linked list.
     * Cut-off radius is taken from the size of current system. */
    void coarsen(Surface &surf, Coarseners &coarseners);

    /** Coarsen the surface by using the linked list.
     * Cut-off radius is taken from the size of provided system. */
    void coarsen(Surface &surf, Coarseners &coarseners, const Medium::Sizes &s);

    /** Smoothen the atoms inside the cylinder */
    void smoothen(const double radius, const double smooth_factor, const double r_cut);

    /** Smoothen all the atoms in the system */
    void smoothen(const double smooth_factor, const double r_cut);

private:
    /** Function used to smoothen the atoms */
    inline double smooth_function(const double distance, const double smooth_factor) const;

    /** Separate cylindrical region from substrate region */
    void get_nanotip(Surface& nanotip, const double radius);
};

} /* namespace femocs */

#endif /* SURFACE_H_ */
