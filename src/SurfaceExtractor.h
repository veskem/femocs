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
    const void rectangularize(shared_ptr<Surface> bulk, const Femocs::SimuCell* cell);

private:
    string adapter;  //!< specify the method to extract the surface (coordination or centrosymmetry)
    int nnn;            //!< number of nearest neighbours for given crystallographic structure
    int nrOfOmpThreads; //!< number of OpenMP threads
    double latconst;    //!< lattice constant of extractable material

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

    const bool on_edge(const double x, const double x_boundary);
};

} /* namespace femocs */

#endif /* SURFACEEXTRACTOR_H_ */
