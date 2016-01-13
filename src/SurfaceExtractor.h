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
#include "Medium.h"

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
     * @param nnn - expected number of nearest neighbours
     * @param nt - number of OpenMP threads
     */
    SurfaceExtractor(const string& adapter, const double cutoff, const int nnn, const int nt);
    virtual ~SurfaceExtractor() {};

    /**
     * Function to choose suitable method to extract surface.
     * @param data - x, y and z coordinates of atoms from AtomReader
     * @return extracted surface
     */
    const virtual shared_ptr<Surface> extract_surface(const AtomReader::Data data);

private:
    string adapter;     //!< specify the method to extract the surface (coordination or centrosymmetry)
    double cutoff2;	    //!< squared cutoff radius for coordination and common neighbour analysis
    int nnn;            //!< number of nearest neighbours for given crystallographic structure
    double surfEps;	    //!< distance from the simulation cell edges that contain non-interesting surface
    int nrOfOmpThreads; //!< number of OpenMP threads

    /**
     * Function to check whether the atom is on the boundary of simulation cell.
     * @param i - index of atom under investigation
     * @return Boolean whether the atom is on the edge of horizontal plane or on the bottom of simulation cell
     */
    const virtual bool isBoundary(const int i, const AtomReader::Data dat);

    /**
     * Extract surface by the coordination analysis - atoms having coordination less than the number
     * of nearest neighbours in given crystal are considered to belong to surface.
     */
    const virtual shared_ptr<Surface> coordination_extract(const AtomReader::Data dat);

    /**
     * Extract surface by the centrosymmetry analysis - atoms having centrosymmetry less than
     * threshold are considered to belong to surface.
     * @todo Implement it!
     */
    const virtual shared_ptr<Surface> centrosymmetry_extract(const AtomReader::Data dat);
};

} /* namespace femocs */

#endif /* SURFACEEXTRACTOR_H_ */
