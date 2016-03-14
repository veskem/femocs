/*
 * SurfaceExtracter.h
 *
 *  Created on: 14.12.2015
 *      Author: veske
 */

#ifndef SURFACEEXTRACTOR_H_
#define SURFACEEXTRACTOR_H_

#include <memory>
#include <string>

#include "AtomReader.h"
#include "Femocs.h"
#include "Surface.h"

using namespace std;
namespace femocs {

/**
 * Class to extract surface atoms from bulk material.
 */
class SurfaceExtractor {
public:
    /**
     * Constructor for surface atom extractor from bulk material.
     * @param adapter - method used to extract the surface
     * @param cutoff - cutoff distance used in coordination calculation
     * @param latconst - lattice constant of extracted material
     * @param nnn - expected number of nearest neighbours
     * @param nt - number of OpenMP threads
     */
    SurfaceExtractor(const Femocs::Config* conf);
    virtual ~SurfaceExtractor() {
    }
    ;

    /**
     * Function to choose suitable method to extract surface.
     * @param data - x, y and z coordinates and types of atoms
     * @param cell - simulation cell parameters
     * @return extracted surface
     */
    const virtual shared_ptr<Surface> extract_surface(const AtomReader::Data* data,
            const Femocs::SimuCell* cell);

    const shared_ptr<Surface> extract_reduced_bulk(const shared_ptr<Surface> surf, const Femocs::SimuCell* cell);
    const shared_ptr<Surface> extract_truncated_bulk(const AtomReader::Data* data, const Femocs::SimuCell* cell);
    const shared_ptr<Surface> extract_bulk(const AtomReader::Data* data, const Femocs::SimuCell* cell);
    const shared_ptr<Surface> extract_edge(const shared_ptr<Surface> surf, const Femocs::SimuCell* cell);
    const void rectangularize_bulk(shared_ptr<Surface> bulk, const Femocs::SimuCell* cell);


private:
    string adapter;  //!< specify the method to extract the surface (coordination or centrosymmetry)
    double cutoff2;	    //!< squared cutoff radius for coordination and common neighbour analysis
    int nnn;            //!< number of nearest neighbours for given crystallographic structure
    double surfEps; //!< distance from the simulation cell edges that contain non-interesting surface
    int nrOfOmpThreads; //!< number of OpenMP threads
    double latconst;    //!< lattice constant of extractable material

    /**
     * Function to check whether the atom is on the boundary of simulation cell.
     * @param i - index of atom under investigation
     * @param dat - atom coordinates data
     * @param cell - simulation cell data
     * @return Boolean whether the atom is on the edge of horizontal plane or on the bottom of simulation cell
     */
    const virtual bool is_boundary(const int i, const AtomReader::Data* dat,
            const Femocs::SimuCell* cell);

    /**
     * Extract surface by the atom types given by kMC simulation.
     */
    const shared_ptr<Surface> kmc_extract(const AtomReader::Data* dat,
            const Femocs::SimuCell* cell);

    /**
     * Extract surface by the coordination analysis - atoms having coordination less than the number
     * of nearest neighbours in given crystal are considered to belong to surface.
     */
    const virtual shared_ptr<Surface> coordination_extract(const AtomReader::Data* dat,
            const Femocs::SimuCell* cell);

    /**
     * Extract surface by the centrosymmetry analysis - atoms having centrosymmetry less than
     * threshold are considered to belong to surface.
     * @todo Implement it!
     */
    const virtual shared_ptr<Surface> centrosymmetry_extract(const AtomReader::Data* dat,
            const Femocs::SimuCell* cell);

    const bool on_edge(const double x, const double y, const Femocs::SimuCell* cell);
    const bool on_edge_vol2(const double x, const double x_boundary);
};

} /* namespace femocs */

#endif /* SURFACEEXTRACTOR_H_ */
