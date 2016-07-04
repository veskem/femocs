/*
 * Macros.cpp
 *
 *  Created on: 4.4.2016
 *      Author: veske
 */

#include "Macros.h"

#include <omp.h>
#include <stdexcept>
#include <sstream>

using namespace std;

/* Template to convert data to string */
template<typename T>
inline string d2s(T data) {
    ostringstream o;
    if (!(o << data)) throw runtime_error("Bad conversion of data to string!");
    return o.str();
}

// Function to handle failed requirement
void __requirement_fails(const char *file, int line, string message) {
    string exc = "Exception:\nfile = " + string(file) + "\nline = " + d2s(line) + "\nmessage = "
            + message + "\n";
    throw runtime_error(exc);
}

// Function to handle failed expectation
void __expectation_fails(const char *file, int line, string message) {
    string exc = "Warning:\nfile = " + string(file) + "\nline = " + d2s(line) + "\nmessage = "
            + message + "\n";
    cout << exc << endl;
}

const double __start_msg(const char* message) {
    cout << endl << string(message);
    cout.flush();
    return omp_get_wtime();
}

const void __end_msg(const double t0) {
    cout << ", time: " << omp_get_wtime() - t0 << endl;
}

// Return mask of indices that doesn't contain the entry
const vector<bool> vector_not(const vector<int> *v, const int entry) {
    vector<bool> mask(v->size());
    for (int i = 0; i < v->size(); ++i)
        mask[i] = (*v)[i] != entry;
    return mask;
}

// Return mask of indices that do contain the entry
const vector<bool> vector_equal(vector<int> *v, const int entry) {
    vector<bool> mask(v->size());
    for (int i = 0; i < v->size(); ++i)
        mask[i] = (*v)[i] == entry;
    return mask;
}
