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
#include "Mesh.h"

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
    const Edge coarsen(const double r_cut);

    const Edge clean();

private:
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

    /**
     * Function to choose suitable method to extract surface.
     * @param data - x, y and z coordinates and types of atoms
     * @param cell - simulation cell parameters
     */
    const void extract_surface(AtomReader* reader);

    const Surface coarsen(const double r_cut, const double zmax, const AtomReader::Sizes* sizes);

    /** Function to flatten the atoms on the sides of simulation box */
    const void rectangularize(const AtomReader::Sizes* sizes, const double r_cut);

    const Surface rectangularize_vol2(const AtomReader::Sizes* sizes, const double r_cut);

    const bool contains(Point2d &p);
    const bool contains(Point3d &p);

    const Surface clean();

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
    // Determine whether an atom is near the edge of simulation box
    const bool on_edge(const double x, const double x_boundary);
};

} /* namespace femocs */

#endif /* MEDIA_H_ */
