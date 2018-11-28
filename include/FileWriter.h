/*
 * FileWriter.h
 *
 *  Created on: 26.11.2018
 *      Author: veske
 */

#ifndef FILEWRITER_H_
#define FILEWRITER_H_

#include <string>

using namespace std;
namespace femocs {

/** General class for holding common routines for writing files */
class FileWriter {
public:
    FileWriter();
    FileWriter(const double delta_time);
    virtual ~FileWriter() {}

    /** By reading the file extension, pick the routine to write the data into file.
     * @param file   path to output file
     * @param force  omit the control for last write time */
    void write(const string &file, bool force=false);

    /** Size of the data vector */
    virtual int size() const { return 0; }

    /** Specify the min time in fs between consequent file writes.
     * 0 enables writing without time restriction, < 0 turns writing off. */
    void set_write_period(const double delta_time);

protected:
    bool not_write_time() const;

    /** Check if given file type is implemented for given class */
    virtual bool valid_extension(const string &extension) const = 0;

    // General routines for writing file with given extension
    virtual void write_xyz(ofstream &out) const;
    virtual void write_vtk(ofstream &out) const;
    virtual void write_msh(ofstream &out) const;
    virtual void write_restart(ofstream &out) const;
    virtual void write_bin(ofstream &out) const;
    virtual void write_dat(ofstream &out) const;
    virtual void write_ckx(ofstream &out) const;

    /** Label for restart data */
    virtual string get_restart_label() const;

    /* Routines to write data in VTK format */
    virtual void write_vtk_points_and_cells(ofstream &out) const {}
    virtual void write_vtk_point_data(ofstream &out) const {}
    virtual void write_vtk_cell_data(ofstream &out) const {}

private:
    double delta_time;       ///< Min time between consequent file writes [fs]
    double last_write_time;  ///< Last time file was written [fs]
};


} // namespace femocs

#endif /* FILEWRITER_H_ */
