/*
 * Media.h
 *
 *  Created on: 5.4.2016
 *      Author: veske
 */

#ifndef MEDIA_H_
#define MEDIA_H_

#include "Coarseners.h"
#include "Macros.h"
#include "AtomReader.h"
#include "Medium.h"
#include "TetgenMesh.h"

using namespace std;
namespace femocs {

/** Routines and data related to making an Edge */
class Edge: public Medium {
public:
    /** Edge constructor */
    Edge();

    /** Extract the atoms near the simulation box sides */
    void extract(const Medium* atoms, const AtomReader::Sizes& ar_sizes, const double eps);
};

/** Routines and data related to making a Surface */
class Media: public Medium {
public:
    /** Empty Media constructor */
    Media();

    /** Initialise empty Media and reserve memory for atoms */
    Media(const int n_atoms);

    /** Generate simple Media with 4 atoms on simubox edges */
    Media(const Medium::Sizes& ar_sizes, const double z);

    void generate_simple(const Medium::Sizes& ar_sizes, const double z);

    void generate_middle(const Medium::Sizes& ar_sizes, const double z, const double dist);

    /** Extract atom with desired types and clean it from lonely atoms;
     * atom is considered lonely if its coordination is lower than threshold.
     * @param reader  AtomReader holding the atoms and their types
     * @param type    type of the atoms that will be read
     * @param invert  if true all the atoms except the 'type'-ones will be stored
     */
    void extract(const AtomReader& reader, const int type, const bool invert=false);

    Media get_nanotip(const double radius);

    Media coarsen(Coarseners &coarseners);

    Media extend(const double latconst, const double box_width, Coarseners &coarseners);

    Media extend(const string &file_name, Coarseners &coarseners);

    Media clean(Coarseners &coarseners);

    /** Remove the atoms that are too far from surface faces */
    void clean(const TetgenMesh& mesh, const double r_cut);

    void smoothen(const double radius, const double smooth_factor, const double r_cut);

    void smoothen(const Point2 &origin, const double radius, const double smooth_factor, const double r_cut);

private:
    inline double smooth_function(const double distance, const double smooth_factor) const;

    void smoothen(const double smooth_factor, const double r_cut);
};

} /* namespace femocs */

#endif /* MEDIA_H_ */
