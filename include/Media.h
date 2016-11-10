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
    const void extract(const Medium* atoms, const AtomReader::Sizes* sizes, const double eps);

    const void generate_uniform(const AtomReader::Sizes* sizes, const double z, const double eps);
};


/** Routines and data related to making a Vacuum */
class Vacuum: public Medium {
public:
    /** Vacuum constructor */
    Vacuum();

    /** Generates Vacuum by adding four points to the top of simulation cell */
    const void generate_simple(const AtomReader::Sizes* sizes);
};

/** Routines and data related to making a Surface */
class Surface: public Medium {
public:
    /** Surface constructors */
    Surface(const int n_atoms);
    Surface();

    const void generate_simple(const AtomReader::Sizes* sizes, const double z);

    /** Extract surface by the atom types */
    const void extract(AtomReader* reader);

    const Surface coarsen(Coarseners &coarseners, const AtomReader::Sizes* reader);

    /** Function to flatten the atoms on the sides of simulation box */
    const Surface rectangularize(const AtomReader::Sizes* sizes, const double eps, const double latconst);

    const Surface clean(Coarseners &coarseners);

    const Surface clean_lonely_atoms(const double r_cut);

    const void smoothen(double radius, double smooth_factor, double r_cut);

    const void smoothen(const Point2 &origin, double radius, double smooth_factor, double r_cut);

private:
    inline double smooth_function(double distance, double smooth_factor) const;

    const void smoothen(double smooth_factor, double r_cut);
};

/** Routines and data related to making a Bulk */
class Bulk: public Medium {
public:
    Bulk();

    const void generate_simple(const AtomReader::Sizes* sizes);

    /** Function to extract bulk material from input atomistic data */
    const void extract(AtomReader* reader);

    const void rectangularize(const AtomReader::Sizes* sizes, const double eps, const double latconst);
};

} /* namespace femocs */

#endif /* MEDIA_H_ */
