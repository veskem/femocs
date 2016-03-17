/*
 * Medium.h
 *
 *  Created on: 21.1.2016
 *      Author: veske
 */

#ifndef MEDIUM_H_
#define MEDIUM_H_

#include <string>
#include <vector>

using namespace std;
namespace femocs {

class Medium {
public:
    /**
     * Function to export the data of Medium.
     * @param file_name - path for file to save the data
     */
    const void output(const string& file_name);

    /** Return x-coordinate of i-th atom */
    const double get_x(const int i);
    /** Return y-coordinate of i-th atom */
    const double get_y(const int i);
    /** Return z-coordinate of i-th atom */
    const double get_z(const int i);
    /** Return type of i-th atom */
    const int get_type(const int i);
    /** Return number of atoms in a Medium */
    const int get_n_atoms();
    /** Sort Medium atoms by their radial coordinate */
    void sort_atoms(vector<int>* permutation_indxs);

protected:
    vector<double> x, y, z; //!< Real coordinates
    vector<int> type;       //!< Atom type
    int N;                  //!< Nr of atoms in Medium
};

} /* namespace femocs */

#endif /* MEDIUM_H_ */
