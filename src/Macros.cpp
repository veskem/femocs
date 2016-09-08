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
#include <string.h>
#include <algorithm>

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
    const int row_len = 45;
    const size_t msglen = strlen(message);

    int whitespace_len = 0;
    if (row_len > msglen && message[msglen-1] != '\n')
        whitespace_len = row_len - msglen;

    cout << endl << string(message) << string(whitespace_len, ' ');
    cout.flush();
    return omp_get_wtime();
}

const void __end_msg(const double t0) {
    cout << "time: " << omp_get_wtime() - t0 << endl;
}

// Return mask of indices that are not equal to the scalar
const vector<bool> vector_not(const vector<int> *v, const int s) {
    return __vector_compare<int, std::not_equal_to<int>>(v, s);
}

// Return mask of indices that are equal to the scalar
const vector<bool> vector_equal(const vector<int> *v, const int s) {
    return __vector_compare<int, std::equal_to<int>>(v, s);
}

// Return mask of indices that are greater than the scalar
const vector<bool> vector_greater(const vector<double> *v, const double s) {
    return __vector_compare<double, std::greater<double>>(v, s);
}
// Return mask of indices that are greater or equal than the entry
const vector<bool> vector_ge(const vector<double> *v, const double s) {
    return __vector_compare<double, std::greater_equal<double>>(v, s);
}

// Return mask of indices that are less than the scalar
const vector<bool> vector_less(const vector<double> *v, const double s) {
    return __vector_compare<double, std::less<double>>(v, s);
}

// Return mask of indices that are less or equal than the scalar
const vector<bool> vector_less_equal(const vector<double> *v, const double s) {
    return __vector_compare<double, std::less_equal<double>>(v, s);
}

vector<size_t> get_shuffle_indices(const int nmax, const int nmin) {
    // Generate indices from [nmin nmax]
    require(nmax >= nmin, "Invalid arguments: " + to_string(nmax) + ", " + to_string(nmin));

    vector<size_t> indxs(nmax - nmin + 1);
    size_t n(nmin);
    generate(indxs.begin(), indxs.end(), [&]{ return n++; });

    // Shuffle the generated indices
    random_shuffle(indxs.begin(), indxs.end());

    return indxs;
}
// Return sorting indexes for vector
vector<size_t> get_sort_indices(const vector<double> &v, const string& direction) {
    // initialize original index locations
    vector<size_t> idx(v.size());
    iota(idx.begin(), idx.end(), 0);

    // sort indexes based on comparing values in v
    if (direction == "asc" || direction == "up")
      sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});
    else if (direction == "desc" || direction == "down")
      sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) {return v[i1] > v[i2];});

    return idx;
}

vector<size_t> get_sort_indices(const vector<double> &v1, const vector<double> &v2, const string& direction) {
    // initialize original index locations
    require(v1.size() == v2.size(), "Input vectors must be with equal size!");

    vector<size_t> idx(v1.size());
    iota(idx.begin(), idx.end(), 0);

    // Get sort indexes by sorting first v1 and then v2
    if (direction == "asc" || direction == "up")
        sort( idx.begin(), idx.end(), [&v1, &v2](size_t a, size_t b)
                { return (v1[a] < v1[b]) || ((v1[a] == v1[b]) && (v2[a] < v2[b])); } );

    else if (direction == "desc" || direction == "down")
        sort( idx.begin(), idx.end(), [&v1, &v2](size_t a, size_t b)
                { return (v1[a] > v1[b]) || ((v1[a] == v1[b]) && (v2[a] > v2[b])); } );

    return idx;
}
