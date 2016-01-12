/*
 * Medium.h
 *
 *  Created on: 14.12.2015
 *      Author: veske
 */

#ifndef MEDIUM_H_
#define MEDIUM_H_

//#include "Smoother.h"
#include "AtomReader.h"
#include <vector>
using namespace std;

namespace femocs {

/**
 * Routines and data related to making a surface
 */
class Surface {
	public:
		/** Surface constructor */
		Surface(int);
		/** Surface destructor */
		~Surface(){};
		/** Function to reduce surface roughness */
		//void smooth(Smoother s);
		/** Function to add atom to surface */
		const void add(const double, const double, const double, const int, const int);
		const void output(const string&);
		const double getX(const int);
		const double getY(const int);
		const double getZ(const int);
		const int getType(const int);
		const int getN();

	private:
		vector <double> x, y, z;	//< Real coordinates
		vector <double> Sx, Sy, Sz; //< Surface area components
		vector <int> coordination;	//< Number of nearest neighbours in cutoff radius
		vector <int> type;			//< Atom type (-1-fixed, 1-bulk, 2-surface, 3-vacancy)
		vector <bool> isEvaporated; //< List of evaporated atoms
};

} /* namespace femocs */

#endif /* MEDIUM_H_ */
