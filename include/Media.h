/*
 * Media.h
 *
 *  Created on: 5.4.2016
 *      Author: veske
 */

#ifndef MEDIA_H_
#define MEDIA_H_

#include "Macros.h"
#include "AtomReader.h"
#include "Medium.h"

using namespace std;
namespace femocs {

/** Routines and data related to making an Edge */
class Edge: public Medium {
public:
    /** Constructor for Edge */
    Edge(const double latconst, const int nnn);
    Edge();

    /** Extract the atoms near the simulation box sides */
    const void extract_edge(Medium* atoms, const AtomReader::Sizes* sizes, const double r_cut);

    /** Make sure all the atoms in Edge are not closer than r_cut */
    const Edge clean(const double r_cut);

    const Edge clean(const double r_cut, const double coarse_factor);
};


/** Routines and data related to making a Vacuum */
class Vacuum: public Medium {
public:
    /** Vacuum constructor */
    Vacuum();

    /** Generates Vacuum by adding four points to the top of simulation cell */
    const void generate_simple(const AtomReader::Sizes* sizes);

private:
};

/** Routines and data related to making a Surface */
class Surface: public Medium {
public:
    /**
     * Surface constructor
     * @param latconst - lattice constant
     * @param nnn - number of nearest neighbors
     */
    Surface(const double latconst, const int nnn);
    Surface();

    /**
     * Function to choose suitable method to extract surface.
     * @param data - x, y and z coordinates and types of atoms
     * @param cell - simulation cell parameters
     */
    const void extract_surface(AtomReader* reader);

    const Surface coarsen(double r_cut, double const_r_cut, double coarse_factor, const AtomReader::Sizes* ar_sizes);

    /** Function to flatten the atoms on the sides of simulation box */
    const Surface rectangularize(const AtomReader::Sizes* sizes, const double r_cut);

    const Surface clean();
    const Surface clean(const double r_cut);
    const Surface clean(const Point3 &origin, double r_in, double r_out, double multiplier);

private:
    /** Extract surface by the atom types given by kMC simulation */
    const void extract_by_type(AtomReader* reader);

    /**
     * Extract surface by the coordination analysis - atoms having coordination less than the number
     * of nearest neighbours in given crystal are considered to belong to surface.
     */
    const void extract_by_coordination(AtomReader* reader);
};

/** Routines and data related to making a Bulk */
class Bulk: public Medium {
public:
    Bulk(const double latconst, const int nnn);
    Bulk();

    // Function to make bulk material with nodes on surface and on by-hand made bottom coordinates
    const void extract_reduced_bulk(Surface* surf, const AtomReader::Sizes* sizes);
    const void generate_simple(const AtomReader::Sizes* sizes, const AtomReader::Types* types);

    // Function to extract bulk material from input atomistic data
    const void extract_truncated_bulk(AtomReader* reader);

    const void extract_truncated_bulk_old(AtomReader* reader);

    // Function to extract bulk material from input atomistic data
    const void extract_bulk(AtomReader* reader);

    // Function to extract bulk material from input atomistic data
    const void rectangularize(const AtomReader::Sizes* sizes, const double r_cut);

private:
};

} /* namespace femocs */

#endif /* MEDIA_H_ */
