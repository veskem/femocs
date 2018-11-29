/*
 * Macros.cpp
 *
 *  Created on: 4.4.2016
 *      Author: veske
 */

#include <omp.h>
#include <algorithm>
#include <fstream>
#include <numeric>
#include "Macros.h"
#include "Globals.h"

using namespace std;
namespace femocs {

const string FEMOCSLOGPATH = "out/femocs.log";

/* Template to return mask of indices that satisfy the comparison condition with an entry */
template<typename T, typename Op>
vector<bool> __vector_compare(const vector<T> *v, const T entry) {
    vector<bool> mask(v->size());
    Op op;
    for (size_t i = 0; i < v->size(); ++i)
        mask[i] = op((*v)[i], entry);

    return mask;
}

void requirement_fails(const char *file, int line, string message) {
    string exc = "\nFEMOCS ERROR:\nfile = " + string(file) + "\nline = " + d2s(line) + "\n" + message;
    throw runtime_error(exc + "\n");
}

void expectation_fails(const char *file, int line, string message) {
    string exc = "\nFEMOCS WARNING:\nfile = " + string(file) + "\nline = " + d2s(line) + "\n" + message;
    cout << exc << endl;
}

double periodic_image(double p, double max, double min){
    require(max > min, "maxbox  must be  > minbox. max = " + to_string(max) + "min = " + to_string(min));
    double from_max = p - max;
    if (from_max > 0) return min + from_max;
    double from_min = p - min;
    if (from_min < 0) return max + from_min;
    return p;
}

void start_msg(double& t0, const string& message) {
    write_log(message);
    if (MODES.VERBOSE) {
        const int row_len = 42;
        const size_t msg_len = message.size();

        int whitespace_len = 0;
        if (row_len > msg_len && message[msg_len-1] != '\n')
            whitespace_len = row_len - msg_len;

        if (message[0] == '=')
            cout << "\n" << message;
        else
            cout << "\n=== " << message << "..." << string(whitespace_len, ' ');
        cout.flush();

        t0 = omp_get_wtime();
    }
}

void end_msg(const double t0) {
    if (MODES.VERBOSE) printf("time: %.3f\n", omp_get_wtime() - t0);
}

void write_silent_msg(const string& message) {
    write_log("\n" + message + "\n");
    if (!MODES.MUTE) cout << "\nFEMOCS: " << message << endl;
}

void write_verbose_msg(const string& message) {
    write_log("  " + message);
    if (MODES.VERBOSE) cout << "  " << message << endl;
}

void write_log(const string& message) {
    if (MODES.WRITELOG) {
        ofstream logfile(FEMOCSLOGPATH, ios_base::app);
        if (logfile) logfile << message << endl;
    }
}

void clear_log() {
    if (MODES.WRITELOG && MODES.SHORTLOG) {
        const string cmd = "rm -f " + FEMOCSLOGPATH;
        if ( system(cmd.c_str()) ) return;
    }
}

vector<bool> vector_not(const vector<int> *v, const int s) {
    return __vector_compare<int, std::not_equal_to<int>>(v, s);
}

vector<bool> vector_equal(const vector<int> *v, const int s) {
    return __vector_compare<int, std::equal_to<int>>(v, s);
}

vector<bool> vector_greater(const vector<double> *v, const double s) {
    return __vector_compare<double, std::greater<double>>(v, s);
}

vector<bool> vector_greater(const vector<int> *v, const int s) {
    return __vector_compare<int, std::greater<int>>(v, s);
}

vector<bool> vector_greater_equal(const vector<double> *v, const double s) {
    return __vector_compare<double, std::greater_equal<double>>(v, s);
}

vector<bool> vector_greater_equal(const vector<int> *v, const int s) {
    return __vector_compare<int, std::greater_equal<int>>(v, s);
}

vector<bool> vector_less(const vector<double> *v, const double s) {
    return __vector_compare<double, std::less<double>>(v, s);
}

vector<bool> vector_less(const vector<int> *v, const int s) {
    return __vector_compare<int, std::less<int>>(v, s);
}

vector<bool> vector_less_equal(const vector<double> *v, const double s) {
    return __vector_compare<double, std::less_equal<double>>(v, s);
}

vector<bool> vector_less_equal(const vector<int> *v, const int s) {
    return __vector_compare<int, std::less_equal<int>>(v, s);
}

vector<int> get_sort_indices(const vector<int> &v, const string& direction) {
    // initialize original index locations
    vector<int> idx(v.size());
    iota(idx.begin(), idx.end(), 0);

    // sort indexes based on comparing values in v
    if (direction == "asc" || direction == "up")
        sort( idx.begin(), idx.end(), [&v](size_t i1, size_t i2) {return v[i1] < v[i2];} );
    else if (direction == "desc" || direction == "down")
        sort( idx.begin(), idx.end(), [&v](size_t i1, size_t i2) {return v[i1] > v[i2];} );

    return idx;
}

int vector_sum(const vector<bool> &v) {
    return accumulate(v.begin(), v.end(), 0);
}

int vector_sum(const vector<int> &v) {
    return accumulate(v.begin(), v.end(), 0);
}

double vector_sum(const vector<double> &v) {
    return accumulate(v.begin(), v.end(), 0.0);
}

bool on_boundary(const double val, const double boundary1, const double boundary2, const double eps) {
    return (fabs(val - boundary1) <= eps) || (fabs(val - boundary2) <= eps);
}

bool on_boundary(const double val, const double boundary, const double eps) {
    return fabs(val - boundary) <= eps;
}

string get_file_type(const string& file_name) {
    return file_name.substr(file_name.find_last_of('.') + 1);
}

string return_class_name(const string& typeid_name) {
    size_t start = typeid_name.find_last_of("1234567890") + 1;
    size_t end = typeid_name.size() - start - 1;
    return typeid_name.substr(start, end);
}

} /* namespace femocs */
