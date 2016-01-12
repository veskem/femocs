/*
 * SurfaceExtracter.h
 *
 *  Created on: 14.12.2015
 *      Author: veske
 */

#ifndef SURFACEEXTRACTOR_H_
#define SURFACEEXTRACTOR_H_

#include <iostream>
#include <string>
#include <memory>
#include "Medium.h"
#include "AtomReader.h"
using namespace std;
namespace femocs {

/**
 * Classes to extract surface atoms from bulk material
 */
class SurfaceExtractor {
public:

	SurfaceExtractor(const string&, const double, const int, const int);
	virtual ~SurfaceExtractor(){};
	/**
	 * choose suitable function to extract surface
	 * @param data from AtomReader
	 * @return extracted surface
	 */
	virtual shared_ptr<Surface> extract_surface(const AtomReader::Data);

private:
	string adapter; /** specify the method to extract the surface (coordination or centrosymmetry) */
	double cutoff2;	/** squared cutoff radius for coordination and cna */
	int nnn;		/** number of nearest neighbours for given crystallographic structure */
	double surfEps;	/** distance from the simulation cell edges that contain non-interesting surface */
	int nrOfOmpThreads; /** number of OpenMP threads */

	/**
	 * @param coordinates of atom
	 * @return if the atom is on the edge of horizontal plane or on the bottom of simulation cell
	 */
	virtual bool isBoundary(const int, const AtomReader::Data dat);
	virtual shared_ptr<Surface> coordination_extract(const AtomReader::Data);
//	virtual shared_ptr<Surface> centrosymmetry_extract(const AtomReader::Data);
};

} /* namespace femocs */

#endif /* SURFACEEXTRACTOR_H_ */
