/*
 * Media.h
 *
 *  Created on: 5.4.2016
 *      Author: veske
 */

#ifndef MEDIA_H_
#define MEDIA_H_

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

    /** Increase or decrease the total volume of system without altering the centre of mass */
    void transform(const double latconst);

    /** Remove the atoms that are too far from surface faces */
    void clean_by_triangles(vector<int>& surf2face, SurfaceInterpolator* interpolator, const double r_cut);

    int clean_by_voronois(const double radius, const double latconst, const string& mesh_quality);

    /** Smoothen the atoms inside the cylinder */
    void smoothen(const double radius, const double smooth_factor, const double r_cut);

    /** Smoothen all the atoms in the system */
    void smoothen(const double smooth_factor, const double r_cut);

    /** Smoothen the atoms using Taubin lambda|mu smoothing algorithm */
    void smoothen(const Config& conf, const double r_cut);

private:
    /** Function used to smoothen the atoms */
    inline double smooth_function(const double distance, const double smooth_factor) const;

    /** Calculate neighbour list for atoms.
     * Atoms are considered neighbours if the distance between them is no more than r_cut. */
    void calc_nborlist(vector<vector<unsigned>>& nborlist, const int nnn, const double r_cut) ;

    /** Smoothen the atoms using Taubin lambda|mu algorithm with inverse neighbour count weighting */
    void laplace_smooth(const double scale, const vector<vector<unsigned>>& nborlist);

    /** Separate cylindrical region from substrate region */
    void get_nanotip(Surface& nanotip, const double radius);

    int get_nanotip(Surface& nanotip, vector<bool>& atom_in_nanotip, const double radius);

    int calc_voronois(VoronoiMesh& voromesh, vector<bool>& node_in_nanotip,
            const double radius, const double latconst, const string& mesh_quality);
};

} /* namespace femocs */

#endif /* MEDIA_H_ */
