/*
 * Media.h
 *
 *  Created on: 5.4.2016
 *      Author: veske
 */

#ifndef MEDIA_H_
#define MEDIA_H_

#include "Medium.h"
#include "Femocs.h"
#include "AtomReader.h"

using namespace std;
namespace femocs {


/** Routines and data related to making a Vacuum */
class Vacuum: public Medium {
public:
    /** Vacuum constructor */
    Vacuum();

    /** Generates Vacuum by adding four points to the top of simulation cell */
    const void generate_simple(const Femocs::SimuCell* cell);

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
    const void extract_surface(AtomReader* reader, const Femocs::SimuCell* cell);

private:
    /** Extract surface by the atom types given by kMC simulation */
    const void kmc_extract(AtomReader* reader, const Femocs::SimuCell* cell);

    /**
     * Extract surface by the coordination analysis - atoms having coordination less than the number
     * of nearest neighbours in given crystal are considered to belong to surface.
     */
    const void coordination_extract(AtomReader* reader, const Femocs::SimuCell* cell);
};


/** Routines and data related to making a Bulk */
class Bulk: public Medium {
public:
    Bulk(const double latconst, const int nnn);

    // Function to make bulk material with nodes on surface and on by-hand made bottom coordinates
    const void extract_reduced_bulk(Surface* surf, const Femocs::SimuCell* cell);

    // Function to extract bulk material from input atomistic data
    const void extract_truncated_bulk(AtomReader* reader, const Femocs::SimuCell* cell);

    // Function to extract bulk material from input atomistic data
    const void extract_bulk(AtomReader* reader, const Femocs::SimuCell* cell);

    // Function to extract bulk material from input atomistic data
    const void rectangularize(const Femocs::SimuCell* cell);

private:
    // Determine whether an atom is near the edge of simulation box
    const bool on_edge(const double x, const double x_boundary);
};


/** Routines and data related to making an Edge */
class Edge: public Medium {
public:
    /** Constructor for Edge */
    Edge(const double latconst, const int nnn);

    /** Extract the atoms near the simulation box sides */
    const void extract_edge(Surface* atoms, const Femocs::SimuCell* cell);

private:
    /** Determine whether an atom is near the edge of simulation box */
    const bool on_edge(const double x, const double x_boundary);
};

} /* namespace femocs */

#endif /* MEDIA_H_ */
