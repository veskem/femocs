/*
 * SurfaceExtracter.cpp
 *
 *  Created on: 14.12.2015
 *      Author: veske
 */

#include "SurfaceExtractor.h"

#include <omp.h>
#include <cmath>
#include <vector>


using namespace std;
namespace femocs {

// SurfaceExtractor constructor
SurfaceExtractor::SurfaceExtractor(const string& adapter, const double cutoff, const int nnn,
        const int nt) {
    this->adapter = adapter;
    this->cutoff2 = cutoff * cutoff;
    this->nnn = nnn;
    this->surfEps = cutoff / 1.0;
    this->nrOfOmpThreads = nt;
}

// Function to pick suitable extraction function
// Shared pointer prevents coping of data and makes memory management automatic
const shared_ptr<Surface> SurfaceExtractor::extract_surface(const AtomReader::Data data) {
    if (adapter == "coordination")
        return SurfaceExtractor::coordination_extract(data);
    else
        return shared_ptr<Surface>();
}

// Detect the atoms on the edges of horizontal plane or at the bottom of simulation cell
const bool SurfaceExtractor::isBoundary(const int i, const AtomReader::Data dat) {
    return (dat.xmax - abs(dat.x[i]) < surfEps) | (dat.ymax - abs(dat.y[i]) < surfEps)
            | (dat.z[i] - dat.zmin < surfEps);
}

// Extract surface by coordination analysis
const shared_ptr<Surface> SurfaceExtractor::coordination_extract(const AtomReader::Data dat) {
    int nrOfAtoms = dat.x.size();
    int i, j, coord;
    double r2, dx, dy, dz;
    vector<bool> isSurface(nrOfAtoms);
    vector<int> coords(nrOfAtoms);

    omp_set_num_threads(this->nrOfOmpThreads);

    int nrOfSurfaceAtoms = 0;
    coord = 0;

    // cutoff == 0 turns off the surface extraction
    if (cutoff2 == 0) {
        for (i = 0; i < nrOfAtoms; ++i)
            if (dat.type[i] != 3) {
                isSurface[i] = 1;
                coords[i] = 0;
                nrOfSurfaceAtoms++;
            }

    } else {
#pragma omp parallel for shared(coords,isSurface) private(i,j,dx,dy,dz,r2) reduction(+:coord,nrOfSurfaceAtoms)
        for (i = 0; i < nrOfAtoms; ++i) {
            coord = 0;
            for (j = 0; j < nrOfAtoms; ++j) {
                if (i == j) continue;
                dx = dat.x[i] - dat.x[j];
                dy = dat.y[i] - dat.y[j];
                dz = dat.z[i] - dat.z[j];
                r2 = dx * dx + dy * dy + dz * dz;
                if (r2 <= cutoff2) coord++;
                if (coord >= this->nnn) break; // Coordination can't be bigger than the biggest expected
            }
            coords[i] = coord;
            if ((coord < this->nnn) && (!isBoundary(i, dat) & (dat.type[i] != 3))) {
                nrOfSurfaceAtoms++;
                isSurface[i] = 1;
            }
        }
    }

    // Create empty surface and add atoms by their coord
    shared_ptr<Surface> s(new Surface(nrOfSurfaceAtoms));
    for (i = 0; i < nrOfAtoms; ++i)
        if (isSurface[i]) s->add(dat.x[i], dat.y[i], dat.z[i], coords[i], dat.type[i]);

    return s;
}

// Extract surface by centrosymmetry analysis
const shared_ptr<Surface> SurfaceExtractor::centrosymmetry_extract(const AtomReader::Data dat) {
    return NULL;
}

} /* namespace femocs */
