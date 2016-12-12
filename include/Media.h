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

using namespace std;
namespace femocs {

/** Routines and data related to making an Edge */
class Edge: public Medium {
public:
    /** Edge constructor */
    Edge();

    /** Extract the atoms near the simulation box sides */
    const void extract(const Medium* atoms, const AtomReader::Sizes& ar_sizes, const double eps);
};

/** Routines and data related to making a Surface */
class Media: public Medium {
public:
    /** Empty Media constructor */
    Media();

    /** Initialise empty Media and reserve memory for atoms */
    Media(const int n_atoms);

    /** Generate simple Media with 4 atoms on simubox edges */
    Media(const AtomReader::Sizes& ar_sizes, const double z);

    const void generate_simple(const AtomReader::Sizes& ar_sizes, const double z);

    const void generate_middle(const AtomReader::Sizes& ar_sizes, const double z, const double r_cut);

    /** Extract atom with desired types
     * @param reader  AtomReader holding the atoms and their types
     * @param type    type of the atoms that will be read
     * @param invert  if true all the atoms except the 'type'-ones will be stored
     */
    const void extract(const AtomReader& reader, const int type, const bool invert=false);

    const Media coarsen(Coarseners &coarseners, const AtomReader::Sizes& ar_sizes);

    /** Function to flatten the atoms on the sides of simulation box */
    const Media rectangularize(const AtomReader::Sizes& ar_sizes, const double eps, const double latconst);

    const Media clean(Coarseners &coarseners);

    const Media clean_lonely_atoms(const double r_cut);

    const void smoothen(double radius, double smooth_factor, double r_cut);

    const void smoothen(const Point2 &origin, double radius, double smooth_factor, double r_cut);

private:
    inline double smooth_function(double distance, double smooth_factor) const;

    const void smoothen(double smooth_factor, double r_cut);
};

} /* namespace femocs */

#endif /* MEDIA_H_ */
